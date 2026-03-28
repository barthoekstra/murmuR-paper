# =============================================================================
# Comparison of murmuR metrics with Roy et al. (2025)
# =============================================================================
#
# Roy et al. (2025) "Enhanced forecasting of bird nocturnal migration intensity
# in relation to previous days and synoptic weather patterns"
# International Journal of Biometeorology 69:1617-1630
# https://doi.org/10.1007/s00484-025-02917-4
#
# Key methodological differences:
#   1. Response variable: Roy predicts log-transformed density; we predict raw
#      density. Roy transforms predictions back to log(ind.km-3) profiles
#      before computing R2, so both evaluate R2 on log density profiles.
#   2. Temporal resolution: Roy uses 3 profiles/night (20:00, 23:00, 02:00 UTC,
#      each averaged from 3 quarter-hourly scans); we use ~12 hourly profiles.
#   3. Height bins: Roy uses 250m (12 bins); we use 50m (60 bins, 50-3000m).
#   4. Peak definition: Roy uses top 5% of density values; we use top-10 nights.
#   5. F1 level: Roy's F1 is on individual density values (or possibly VID);
#      ours is night-level.
#
# To make a fair comparison, we compute metrics at three resolution levels:
#   - 50m, all hours (our native resolution)
#   - 250m, all hours (Roy's height bins, our temporal resolution)
#   - 250m, 3 profiles/night at 20:00/23:00/02:00 UTC ("Roy-resolution")
#
# R2 is computed on log(density) excluding zeros (1-4% of observations,
# depending on radar/season). We use log() rather than log1p because log1p
# inflates R2 (compresses the near-zero range; R2 varies from 0.29 to 0.86
# depending on the offset). Roy likely has no zeros after filtering
# incomplete profiles, interpolating to 250m, and excluding insect events.
#
# F1 is computed at the 5% threshold (top 5% = peak) at three aggregation
# levels: per-observation (each height-bin x time-step), per-profile
# (vertically integrated density, VID, per time-step), and per-night
# (nightly VID summed from per-profile VIDs, then ranked).
#
# This script recomputes our metrics to be as comparable as possible to Roy's.
# =============================================================================

library(data.table)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

PRED_PATH <- "/data/birdcloudstorage-tvm/ibm-ml/data/models/predictions"
SAMPLING  <- "dynamic"  # primary sampling strategy from our paper
RADARS    <- c("nlhrw", "nldhl")
SEASONS   <- c("autumn", "spring")
YEARS     <- 2017:2024

# Roy et al. Table 1 results (their full model, averaged across 9 radars)
roy_results <- data.table(
  season = c("spring", "autumn"),
  roy_r2_mean       = c(0.47, 0.61),
  roy_r2_sd         = c(0.051, 0.029),
  roy_f1_mean       = c(0.47, 0.48),
  roy_f1_sd         = c(0.045, 0.088),
  # Their model without synoptic & previous days (local + instantaneous only)
  roy_r2_local_mean = c(0.47, 0.57),
  roy_f1_local_mean = c(0.432, 0.37),
  # Their phenology-only baseline
  roy_r2_pheno_mean = c(0.35, 0.51),
  roy_f1_pheno_mean = c(0.28, 0.36),
  # Van Doren & Horton (2018) baseline
  roy_r2_vdh_mean   = c(0.44, 0.55),
  roy_f1_vdh_mean   = c(0.40, 0.37)
)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

#' Load prediction file for a given radar/season/year
load_preds <- function(radar, season, year, sampling = SAMPLING) {
  fname <- glue::glue("{radar}_{season}_{year}_{sampling}_preds.RDS")
  fpath <- file.path(PRED_PATH, glue::glue("regr.rmse_{sampling}"), radar, fname)
  if (!file.exists(fpath)) return(NULL)
  readRDS(fpath)
}

#' Compute R-squared
rsq <- function(obs, pred) {
  ss_res <- sum((obs - pred)^2, na.rm = TRUE)
  ss_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  if (ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

#' Compute F1 score for binary peak classification
#' @param obs_peak Logical vector: is observation a true peak?
#' @param pred_peak Logical vector: is observation a predicted peak?
f1_score <- function(obs_peak, pred_peak) {
  tp <- sum(obs_peak & pred_peak, na.rm = TRUE)
  fp <- sum(!obs_peak & pred_peak, na.rm = TRUE)
  fn <- sum(obs_peak & !pred_peak, na.rm = TRUE)
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
  recall    <- if ((tp + fn) > 0) tp / (tp + fn) else 0
  if ((precision + recall) == 0) return(0)
  2 * precision * recall / (precision + recall)
}

#' Coarsen 50m bins to 250m bins by averaging density within each group.
#' Our bins: 50,100,...,3000 (60 bins). ceiling(height / 250) * 250 maps
#' 50-250 -> 250, 300-500 -> 500, ..., 2750-3000 -> 3000. Yields 12 bins.
coarsen_250m <- function(dt) {
  dt[, height_bin250 := ceiling(height / 250) * 250]
  dt250 <- dt[, .(
    dens     = mean(dens),
    response = mean(response)
  ), by = .(datetime, night_num, height_bin250)]
  setnames(dt250, "height_bin250", "height")
  dt250
}

#' Subsample to Roy's 3 profiles per night: 20:00, 23:00, 02:00 UTC.
#' Roy averages 15-min radar scans at h+0, h+15, h+30 for each of these hours.
#' Our data is hourly, so we select the nearest matching hours directly.
subsample_3profiles <- function(dt) {
  dt[, hour := as.integer(format(datetime, "%H"))]
  dt_sub <- dt[hour %in% c(20, 23, 2)]
  dt_sub[, hour := NULL]
  dt_sub
}

# ---------------------------------------------------------------------------
# 1. Compute metrics across all radar/season/year combinations
# ---------------------------------------------------------------------------

results <- rbindlist(lapply(RADARS, function(radar) {
  rbindlist(lapply(SEASONS, function(season) {
    rbindlist(lapply(YEARS, function(year) {
      dt <- load_preds(radar, season, year)
      if (is.null(dt)) return(NULL)

      dt <- dt[!is.na(dens) & !is.na(response)]

      n_obs <- nrow(dt)
      n_nights <- uniqueN(dt$night_num)

      # ===================================================================
      # A. R-squared on different scales
      # ===================================================================

      # A1. R2 on raw density (per-height-hour) -- our standard metric
      r2_raw <- rsq(dt$dens, dt$response)

      # A2. R2 on log(density) excluding zeros (per-height-hour)
      # Roy uses log(density) directly. Their data has no zeros because they
      # filter incomplete profiles and insect-dominated events. We exclude
      # zeros (typically ~4% of observations) to match.
      # NOTE: log1p vs log(x+epsilon) choice has massive impact on R2
      # (R2 ranges from 0.29 to 0.86 depending on offset). Excluding zeros
      # and using log() is most comparable to Roy's approach.
      dt_nz <- dt[dens > 0 & response > 0]
      r2_log_nozero <- rsq(log(dt_nz$dens), log(dt_nz$response))
      pct_zero <- 1 - nrow(dt_nz) / nrow(dt)

      # A2b. R2 with log1p (for reference -- inflated by offset compressing
      # the range of near-zero values)
      r2_log1p <- rsq(log1p(dt$dens), log1p(dt$response))

      # ===================================================================
      # B. F1 scores at different thresholds and aggregation levels
      # ===================================================================

      # B1. Per-observation F1 at 5% threshold (Roy's likely method)
      # "peaks = top 5% of the highest bird density values"
      obs_q95 <- quantile(dt$dens, 0.95)
      pred_q95 <- quantile(dt$response, 0.95)
      f1_obs_5pct <- f1_score(dt$dens >= obs_q95, dt$response >= pred_q95)

      # B1b. Per-profile F1 at 5% threshold on vertically integrated density
      # VID = sum(dens) * 0.05 km (50m bin width) per datetime profile.
      # Include night_num for subsequent nightly aggregation in B2.
      profiles <- dt[, .(
        vid_obs  = sum(dens) * 0.05,
        vid_pred = sum(response) * 0.05
      ), by = .(datetime, night_num)]
      vid_obs_q95  <- quantile(profiles$vid_obs, 0.95)
      vid_pred_q95 <- quantile(profiles$vid_pred, 0.95)
      f1_vid_5pct <- f1_score(profiles$vid_obs >= vid_obs_q95,
                               profiles$vid_pred >= vid_pred_q95)

      # B2. Per-night F1 at 5% threshold
      # Nightly VID: first integrate per profile (VID per datetime in B1b),
      # then sum across profiles per night. Matches add_night_ranks() in
      # R/measures.R. (Numerically equivalent to summing all height-hours at
      # once, but conceptually clearer: VID is a per-profile quantity.)
      nightly <- profiles[, .(
        vid_obs  = sum(vid_obs),
        vid_pred = sum(vid_pred)
      ), by = .(night_num)]
      n_peak_nights_5pct <- max(1, round(n_nights * 0.05))
      nightly[, truth_rank := rank(-vid_obs, ties.method = "min")]
      nightly[, pred_rank  := rank(-vid_pred, ties.method = "min")]
      f1_night_5pct <- f1_score(
        nightly$truth_rank <= n_peak_nights_5pct,
        nightly$pred_rank  <= n_peak_nights_5pct
      )

      # B3. Per-night F1 at top-10 nights (~9.4% for ~106 nights)
      # This is our paper's metric (peak nights identified / 10)
      top_n <- min(10, n_nights)
      f1_night_top10 <- f1_score(
        nightly$truth_rank <= top_n,
        nightly$pred_rank  <= top_n
      )
      peak_nights_identified <- sum(
        nightly[truth_rank <= top_n]$pred_rank <= top_n
      )

      # B4. Per-observation F1 at 10% threshold (for sensitivity check)
      obs_q90 <- quantile(dt$dens, 0.90)
      pred_q90 <- quantile(dt$response, 0.90)
      f1_obs_10pct <- f1_score(dt$dens >= obs_q90, dt$response >= pred_q90)

      # ===================================================================
      # C. Same metrics at 250m resolution (Roy's bin size)
      # ===================================================================
      # Average 50m bins into 250m bins, then recompute R2 and F1.
      # This removes the spatial autocorrelation advantage of finer bins.

      dt250 <- coarsen_250m(dt)

      # C1. R2 at 250m
      r2_raw_250m <- rsq(dt250$dens, dt250$response)
      dt250_nz <- dt250[dens > 0 & response > 0]
      r2_log_nozero_250m <- rsq(log(dt250_nz$dens), log(dt250_nz$response))

      # C2. Per-observation F1 at 5% on 250m bins
      obs_q95_250m  <- quantile(dt250$dens, 0.95)
      pred_q95_250m <- quantile(dt250$response, 0.95)
      f1_obs_5pct_250m <- f1_score(dt250$dens >= obs_q95_250m,
                                    dt250$response >= pred_q95_250m)

      # C3. VID-based F1 at 5% on 250m profiles
      # VID with 250m bins: sum(dens) * 0.25 km
      profiles250 <- dt250[, .(
        vid_obs  = sum(dens) * 0.25,
        vid_pred = sum(response) * 0.25
      ), by = .(datetime)]
      vid_obs_q95_250m  <- quantile(profiles250$vid_obs, 0.95)
      vid_pred_q95_250m <- quantile(profiles250$vid_pred, 0.95)
      f1_vid_5pct_250m <- f1_score(profiles250$vid_obs >= vid_obs_q95_250m,
                                    profiles250$vid_pred >= vid_pred_q95_250m)

      # ===================================================================
      # D. Roy-resolution: 250m bins + 3 profiles/night (20:00, 23:00, 02:00)
      # ===================================================================
      # This controls for both spatial and temporal autocorrelation differences.

      dt_roy <- subsample_3profiles(dt)
      dt_roy <- coarsen_250m(dt_roy)

      # D1. R2 at Roy resolution
      r2_raw_roy <- rsq(dt_roy$dens, dt_roy$response)
      dt_roy_nz <- dt_roy[dens > 0 & response > 0]
      r2_log_nozero_roy <- rsq(log(dt_roy_nz$dens), log(dt_roy_nz$response))

      # D2. Per-observation F1 at 5% at Roy resolution
      obs_q95_roy  <- quantile(dt_roy$dens, 0.95)
      pred_q95_roy <- quantile(dt_roy$response, 0.95)
      f1_obs_5pct_roy <- f1_score(dt_roy$dens >= obs_q95_roy,
                                   dt_roy$response >= pred_q95_roy)

      # D3. VID-based F1 at 5% at Roy resolution
      profiles_roy <- dt_roy[, .(
        vid_obs  = sum(dens) * 0.25,
        vid_pred = sum(response) * 0.25
      ), by = .(datetime, night_num)]
      vid_obs_q95_roy  <- quantile(profiles_roy$vid_obs, 0.95)
      vid_pred_q95_roy <- quantile(profiles_roy$vid_pred, 0.95)
      f1_vid_5pct_roy <- f1_score(profiles_roy$vid_obs >= vid_obs_q95_roy,
                                   profiles_roy$vid_pred >= vid_pred_q95_roy)

      # D4. Per-night F1 at Roy resolution
      # Nightly VID from only 3 profiles (vs ~12 in our native resolution).
      # This changes the VID distribution and thus the 5% threshold.
      nightly_roy <- profiles_roy[, .(
        vid_obs  = sum(vid_obs),
        vid_pred = sum(vid_pred)
      ), by = .(night_num)]
      n_nights_roy <- nrow(nightly_roy)
      n_peak_nights_5pct_roy <- max(1, round(n_nights_roy * 0.05))
      nightly_roy[, truth_rank := rank(-vid_obs, ties.method = "min")]
      nightly_roy[, pred_rank  := rank(-vid_pred, ties.method = "min")]
      f1_night_5pct_roy <- f1_score(
        nightly_roy$truth_rank <= n_peak_nights_5pct_roy,
        nightly_roy$pred_rank  <= n_peak_nights_5pct_roy
      )
      top_n_roy <- min(10, n_nights_roy)
      f1_night_top10_roy <- f1_score(
        nightly_roy$truth_rank <= top_n_roy,
        nightly_roy$pred_rank  <= top_n_roy
      )
      peak_nights_id_roy <- sum(
        nightly_roy[truth_rank <= top_n_roy]$pred_rank <= top_n_roy
      )

      data.table(
        radar   = radar,
        season  = season,
        year    = year,
        n_obs   = n_obs,
        n_nights = n_nights,
        pct_zero = round(pct_zero, 3),
        # R-squared variants (50m)
        r2_raw          = r2_raw,
        r2_log_nozero   = r2_log_nozero,
        r2_log1p        = r2_log1p,
        # R-squared variants (250m)
        r2_raw_250m        = r2_raw_250m,
        r2_log_nozero_250m = r2_log_nozero_250m,
        # R-squared variants (Roy resolution: 250m + 3 profiles/night)
        r2_raw_roy          = r2_raw_roy,
        r2_log_nozero_roy   = r2_log_nozero_roy,
        # F1 variants (50m)
        f1_obs_5pct     = f1_obs_5pct,
        f1_vid_5pct     = f1_vid_5pct,
        f1_obs_10pct    = f1_obs_10pct,
        f1_night_5pct   = f1_night_5pct,
        f1_night_top10  = f1_night_top10,
        peak_nights_id  = peak_nights_identified,
        # F1 variants (250m)
        f1_obs_5pct_250m = f1_obs_5pct_250m,
        f1_vid_5pct_250m = f1_vid_5pct_250m,
        # F1 variants (Roy resolution: 250m + 3 profiles/night)
        f1_obs_5pct_roy      = f1_obs_5pct_roy,
        f1_vid_5pct_roy      = f1_vid_5pct_roy,
        f1_night_5pct_roy    = f1_night_5pct_roy,
        f1_night_top10_roy   = f1_night_top10_roy,
        peak_nights_id_roy   = peak_nights_id_roy,
        n_peak_5pct          = n_peak_nights_5pct
      )
    }))
  }))
}))

# ---------------------------------------------------------------------------
# 2. R-squared summary
# ---------------------------------------------------------------------------
# Most comparable: r2_log_nozero_roy (250m + 3 profiles/night).
#
# R2 increases when coarsening to 250m (averaging smooths noise) and barely
# changes further when subsampling to 3 profiles/night:
#   nlhrw autumn:  0.75 -> 0.80 -> 0.80 (50m -> 250m -> Roy). Roy: 0.61.
#   nlhrw spring:  0.74 -> 0.79 -> 0.78.                      Roy: 0.47.
#   nldhl autumn:  0.50 -> 0.57 -> 0.58.                      Roy: 0.61.
#   nldhl spring:  0.43 -> 0.52 -> 0.52.                      Roy: 0.47.
# nlhrw clearly outperforms Roy. nldhl is below Roy in autumn but above in
# spring. Our finer resolution does not inflate R2.

# Mean R2 across years, by radar and season
r2_summary <- results[, .(
  r2_log_nozero      = round(mean(r2_log_nozero), 3),
  r2_log_nozero_250m = round(mean(r2_log_nozero_250m), 3),
  r2_log_nozero_roy  = round(mean(r2_log_nozero_roy), 3),
  pct_zero           = round(mean(pct_zero), 3),
  n_years            = .N
), by = .(radar, season)]
print(r2_summary)

# R2 by year
r2_by_year <- results[, .(radar, season, year, r2_log_nozero,
                           r2_log_nozero_250m, r2_log_nozero_roy, pct_zero)]
r2_by_year[, (4:6) := lapply(.SD, round, 3), .SDcols = 4:6]
print(r2_by_year)

# ---------------------------------------------------------------------------
# 3. F1 score summary
# ---------------------------------------------------------------------------
# Roy et al. F1: top 5% of "highest bird density values" at 5% threshold.
# Our paper: top-10 nights out of ~106 per season (~9.4% threshold).
#
# The aggregation level matters far more than spatial/temporal resolution:
#   Per-observation F1 (each height-bin x timestep classified independently):
#     50m: 0.69-0.78.  250m: 0.65-0.71.  Roy-res: 0.64-0.71.  Roy: 0.47-0.48.
#     We outperform Roy at all resolutions, even at Roy's exact setup.
#   VID F1 (one VID value per profile, classified as peak or not):
#     50m: 0.28-0.51.  250m: 0.26-0.50.  Roy-res: 0.18-0.50.  Roy: 0.47-0.48.
#     nlhrw autumn (0.50 all-years, 0.47 in 2017-2022 overlap) matches Roy.
#     Other radar/season combos fall below, especially nldhl spring (0.18).
#
# VID and nightly F1 use a two-step aggregation matching add_night_ranks()
# in R/measures.R: first integrate density per profile (VID per datetime),
# then for nightly F1, sum VIDs per night before ranking.
#
# Nightly top-10 F1 at Roy-resolution (3 profiles/night):
#     nlhrw: 0.50 (autumn), 0.50 (spring). nldhl: 0.39, 0.36.
#   (vs 50m all-profiles: 0.63, 0.59, 0.49, 0.31.)
#   Fewer profiles per night generally reduces nightly F1, but nldhl spring
#   actually improves (0.31 -> 0.36) with fewer profiles.
#
# Subsampling to 3 profiles/night has a larger effect on VID F1 than
# coarsening to 250m, consistent with temporal autocorrelation in our
# hourly data boosting the 50m VID F1. Per-observation F1 is less affected
# because correlated height bins already dominate that metric.

# Mean F1 across years, by radar and season
f1_summary <- results[, .(
  f1_obs_5pct          = round(mean(f1_obs_5pct), 3),
  f1_vid_5pct          = round(mean(f1_vid_5pct), 3),
  f1_obs_5pct_250m     = round(mean(f1_obs_5pct_250m), 3),
  f1_vid_5pct_250m     = round(mean(f1_vid_5pct_250m), 3),
  f1_obs_5pct_roy      = round(mean(f1_obs_5pct_roy), 3),
  f1_vid_5pct_roy      = round(mean(f1_vid_5pct_roy), 3),
  f1_night_top10       = round(mean(f1_night_top10), 3),
  f1_night_top10_roy   = round(mean(f1_night_top10_roy), 3),
  peak_nights_id       = round(mean(peak_nights_id), 1),
  peak_nights_id_roy   = round(mean(peak_nights_id_roy), 1),
  n_years              = .N
), by = .(radar, season)]
print(f1_summary)

# F1 by year
f1_by_year <- results[, .(radar, season, year, f1_obs_5pct, f1_vid_5pct,
                           f1_obs_5pct_250m, f1_vid_5pct_250m,
                           f1_obs_5pct_roy, f1_vid_5pct_roy,
                           f1_night_top10, f1_night_top10_roy,
                           peak_nights_id, peak_nights_id_roy)]
f1_by_year[, (4:12) := lapply(.SD, round, 3), .SDcols = 4:12]
print(f1_by_year)

# ---------------------------------------------------------------------------
# 4. Year-by-year variability (2017-2022 overlap with Roy)
# ---------------------------------------------------------------------------
# Roy covers 2017-2022 (6 years), we cover 2017-2024 (8 years).
# Restrict to overlapping years for fairer comparison.

overlap_results <- results[year %in% 2017:2022]

overlap_summary <- overlap_results[, .(
  r2_log_nozero        = round(mean(r2_log_nozero), 3),
  r2_log_nozero_250m   = round(mean(r2_log_nozero_250m), 3),
  r2_log_nozero_roy    = round(mean(r2_log_nozero_roy), 3),
  f1_obs_5pct          = round(mean(f1_obs_5pct), 3),
  f1_vid_5pct          = round(mean(f1_vid_5pct), 3),
  f1_obs_5pct_250m     = round(mean(f1_obs_5pct_250m), 3),
  f1_vid_5pct_250m     = round(mean(f1_vid_5pct_250m), 3),
  f1_obs_5pct_roy      = round(mean(f1_obs_5pct_roy), 3),
  f1_vid_5pct_roy      = round(mean(f1_vid_5pct_roy), 3),
  f1_night_top10       = round(mean(f1_night_top10), 3),
  f1_night_top10_roy   = round(mean(f1_night_top10_roy), 3),
  peak_nights_id       = round(mean(peak_nights_id), 1),
  peak_nights_id_roy   = round(mean(peak_nights_id_roy), 1)
), by = .(radar, season)]
print(overlap_summary)

# ---------------------------------------------------------------------------
# 5. Sensitivity to threshold choice
# ---------------------------------------------------------------------------
# How does F1 change with different peak thresholds?
# This helps understand whether our higher top-10 F1 is an artifact
# of the more lenient threshold.

# Per-observation threshold sensitivity
threshold_results <- rbindlist(lapply(RADARS, function(radar) {
  rbindlist(lapply(SEASONS, function(season) {
    rbindlist(lapply(YEARS, function(year) {
      dt <- load_preds(radar, season, year)
      if (is.null(dt)) return(NULL)
      dt <- dt[!is.na(dens) & !is.na(response)]

      rbindlist(lapply(c(0.90, 0.925, 0.95, 0.975, 0.99), function(q) {
        obs_thr  <- quantile(dt$dens, q)
        pred_thr <- quantile(dt$response, q)
        f1 <- f1_score(dt$dens >= obs_thr, dt$response >= pred_thr)
        data.table(
          radar = radar, season = season, year = year,
          quantile = q, pct_label = paste0("top_", (1 - q) * 100, "pct"),
          f1 = f1
        )
      }))
    }))
  }))
}))

threshold_summary <- threshold_results[, .(
  f1_mean = round(mean(f1), 3),
  f1_sd   = round(sd(f1), 3)
), by = .(radar, season, pct_label, quantile)]
setorder(threshold_summary, radar, season, quantile)
print(threshold_summary)

# VID per-profile threshold sensitivity
vid_threshold_results <- rbindlist(lapply(RADARS, function(radar) {
  rbindlist(lapply(SEASONS, function(season) {
    rbindlist(lapply(YEARS, function(year) {
      dt <- load_preds(radar, season, year)
      if (is.null(dt)) return(NULL)
      dt <- dt[!is.na(dens) & !is.na(response)]

      profiles <- dt[, .(
        vid_obs  = sum(dens) * 0.05,
        vid_pred = sum(response) * 0.05
      ), by = .(datetime)]

      rbindlist(lapply(c(0.90, 0.925, 0.95, 0.975, 0.99), function(q) {
        obs_thr  <- quantile(profiles$vid_obs, q)
        pred_thr <- quantile(profiles$vid_pred, q)
        f1 <- f1_score(profiles$vid_obs >= obs_thr, profiles$vid_pred >= pred_thr)
        data.table(
          radar = radar, season = season, year = year,
          quantile = q, pct_label = paste0("top_", (1 - q) * 100, "pct"),
          f1 = f1
        )
      }))
    }))
  }))
}))

vid_threshold_summary <- vid_threshold_results[, .(
  f1_mean = round(mean(f1), 3),
  f1_sd   = round(sd(f1), 3)
), by = .(radar, season, pct_label, quantile)]
setorder(vid_threshold_summary, radar, season, quantile)
print(vid_threshold_summary)

# Night-level threshold sensitivity
night_threshold_results <- rbindlist(lapply(RADARS, function(radar) {
  rbindlist(lapply(SEASONS, function(season) {
    rbindlist(lapply(YEARS, function(year) {
      dt <- load_preds(radar, season, year)
      if (is.null(dt)) return(NULL)
      dt <- dt[!is.na(dens) & !is.na(response)]

      # Two-step nightly VID: per-profile first, then per-night
      profiles <- dt[, .(
        vid_obs  = sum(dens) * 0.05,
        vid_pred = sum(response) * 0.05
      ), by = .(datetime, night_num)]
      nightly <- profiles[, .(
        vid_obs  = sum(vid_obs),
        vid_pred = sum(vid_pred)
      ), by = .(night_num)]
      n_nights <- nrow(nightly)
      nightly[, truth_rank := rank(-vid_obs, ties.method = "min")]
      nightly[, pred_rank  := rank(-vid_pred, ties.method = "min")]

      rbindlist(lapply(c(3, 5, 7, 10, 15), function(top) {
        f1 <- f1_score(
          nightly$truth_rank <= top,
          nightly$pred_rank  <= top
        )
        overlap <- sum(nightly[truth_rank <= top]$pred_rank <= top)
        data.table(
          radar = radar, season = season, year = year,
          top_n = top, pct = round(top / n_nights * 100, 1),
          f1 = f1, overlap = overlap
        )
      }))
    }))
  }))
}))

night_threshold_summary <- night_threshold_results[, .(
  f1_mean  = round(mean(f1), 3),
  f1_sd    = round(sd(f1), 3),
  overlap  = round(mean(overlap), 1),
  mean_pct = round(mean(pct), 1)
), by = .(radar, season, top_n)]
setorder(night_threshold_summary, radar, season, top_n)
print(night_threshold_summary)

# ---------------------------------------------------------------------------
# 6. Combined comparison table (2017-2022 overlap period)
# ---------------------------------------------------------------------------
# At Roy-resolution (250m bins, 3 profiles/night at 20:00/23:00/02:00 UTC),
# restricted to 2017-2022 (Roy's study period):
#
#                R2 (log)    F1 obs-level   F1 VID-level   F1 nightly top-10
#                Ours  Roy   Ours  Roy      Ours  Roy      50m   Roy-res
# nlhrw autumn   0.79  0.61  0.70  0.48     0.47  0.48     0.62  0.50
# nldhl autumn   0.54  0.61  0.66  0.48     0.31  0.48     0.48  0.40
# nlhrw spring   0.78  0.47  0.67  0.47     0.31  0.47     0.57  0.47
# nldhl spring   0.51  0.47  0.67  0.47     0.16  0.47     0.27  0.35
#
# VID and nightly F1 use two-step aggregation: VID per profile (sum of
# density * bin width), then sum VIDs per night for the nightly metric.
# At Roy-resolution, nightly VID is summed from only 3 profiles (vs ~12).
#
# If Roy uses per-observation F1: we clearly outperform (0.66-0.70 vs 0.47-0.48).
# If Roy uses VID-based F1: only nlhrw autumn matches; other combos are below.
#
# Nightly top-10 F1 drops from 0.62 to 0.50 (nlhrw autumn) when using
# Roy's 3 profiles, because fewer profiles per night change the ranking.
# Interestingly, nldhl spring improves (0.27 -> 0.35), suggesting our model
# is better at ranking nights when temporal noise is removed.

combined <- overlap_summary[, .(radar, season,
  r2_log_nozero, r2_log_nozero_250m, r2_log_nozero_roy,
  f1_obs_5pct, f1_obs_5pct_250m, f1_obs_5pct_roy,
  f1_vid_5pct, f1_vid_5pct_250m, f1_vid_5pct_roy,
  f1_night_top10, f1_night_top10_roy)]
combined <- merge(combined, roy_results[, .(season, roy_r2_mean, roy_f1_mean)],
                  by = "season")
setcolorder(combined, c("radar", "season",
  "r2_log_nozero", "r2_log_nozero_250m", "r2_log_nozero_roy", "roy_r2_mean",
  "f1_obs_5pct", "f1_obs_5pct_250m", "f1_obs_5pct_roy",
  "f1_vid_5pct", "f1_vid_5pct_250m", "f1_vid_5pct_roy", "roy_f1_mean",
  "f1_night_top10", "f1_night_top10_roy"))
setnames(combined, c("radar", "season",
  "R2_50m", "R2_250m", "R2_250m3p", "roy_R2",
  "F1obs_50m", "F1obs_250m", "F1obs_250m3p",
  "F1vid_50m", "F1vid_250m", "F1vid_250m3p", "roy_F1",
  "F1ngt_50m", "F1ngt_250m3p"))
print(combined)

# Save results
outdir <- file.path(dirname(getwd()), "murmuR", "output")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
saveRDS(results, file.path(outdir, "comparison_roy_et_al_results.RDS"))

# CSV with human-readable headers
combined_csv <- copy(combined)
setnames(combined_csv, c(
  "Radar", "Season",
  "R2 (50m, all hours)", "R2 (250m, all hours)", "R2 (250m, 3 profiles)", "R2 Roy et al.",
  "F1 obs (50m, all hours)", "F1 obs (250m, all hours)", "F1 obs (250m, 3 profiles)",
  "F1 VID (50m, all hours)", "F1 VID (250m, all hours)", "F1 VID (250m, 3 profiles)",
  "F1 Roy et al.",
  "F1 nightly (50m, all hours)", "F1 nightly (250m, 3 profiles)"))
write_csv(combined_csv, file.path(outdir, "comparison_roy_et_al_combined.csv"))
