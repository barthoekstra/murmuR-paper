# ============================================================
# Sensitivity analysis: effect of number of IBM migrants on
# predictor stability + cohort proportions
# ============================================================

# --- Libraries -----------------------------------------------
library(murmuR)
library(data.table)
library(ggplot2)
library(glue)

# --- Options -------------------------------------------------
options(murmuR.data_path = "/data/birdcloudstorage-tvm/ibm-ml/data/")

# --- Constants -----------------------------------------------
radars   <- names(radars_named)
years    <- 2017:2024
seasons  <- c("autumn", "spring")

fractions  <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75)
n_repeats  <- 10L

# Representative variables spanning different feature types:
#   fl_inst_u10        — en-route weather (mean of few birds at radar)
#   stp_t2m_24hrs_slope — stopover weather trend (mean of many birds' histories)
#   stp_lat_mean       — spatial distribution of stopover origins
target_vars <- c("fl_inst_u10", "stp_t2m_24hrs_slope", "stp_lat_mean")

out_dir <- "output/sensitivity"
dir.create(file.path(out_dir, "predictor_convergence"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "cohort_proportions"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Helper functions
# ============================================================

#' Subsample migrating birds, keeping all non-migrating birds intact
#'
#' @param birds data.table from simulated_birds() with column `migration`
#' @param fraction Numeric in (0, 1]. Fraction of unique migrating bird IDs to keep.
#' @param seed Integer random seed.
#' @return Filtered birds data.table.
subsample_birds <- function(birds, fraction, seed) {
  mig_ids <- unique(birds[migration == TRUE, bird])
  n_keep  <- max(1L, round(length(mig_ids) * fraction))

  set.seed(seed)
  keep_ids <- sample(mig_ids, n_keep, replace = FALSE)

  rbindlist(list(
    birds[migration == TRUE & bird %in% keep_ids],
    birds[migration == FALSE]
  ))
}


#' Prepare pre-indexed lookup tables for fast per-timestep computation
#'
#' Extracts only the columns needed for the 3 target variables and sets keys
#' for fast indexed lookup. Call once per combo, then pass to
#' compute_predictor_values_fast().
#'
#' @return A named list of pre-indexed data.tables.
prepare_lookup_tables <- function(birds, cond_mig, cond_stp_mig) {
  # Which birds are at the radar per timestep (only migrating, state==1)
  birds_radar <- birds[state == 1 & radar == TRUE & migration == TRUE,
                       .(bird_ids = list(bird)), by = tidx]
  setkey(birds_radar, tidx)

  # fl_inst_u10: only need inst_u10 from instantaneous conditions
  # Strip units class upfront to avoid rbindlist class-mismatch errors
  cond_u10 <- cond_mig[, .(tidx, bird, inst_u10 = as.numeric(inst_u10))]
  setkey(cond_u10, bird, tidx)

  # stp_t2m_24hrs_slope: one value per bird from stopover conditions
  stp_slope <- cond_stp_mig[, .(bird, stp_t2m_24hrs_slope)]
  setkey(stp_slope, bird)

  # stp_lat_mean: stopover locations (state==0) of migrating birds
  stp_locs <- birds[migration == TRUE & state == 0, .(bird, y)]
  setkey(stp_locs, bird)

  list(birds_radar = birds_radar, cond_u10 = cond_u10,
       stp_slope = stp_slope, stp_locs = stp_locs)
}


#' Compute only the 3 target predictor values for timesteps with birds
#'
#' Lean replacement for process_tidx() that computes only fl_inst_u10,
#' stp_t2m_24hrs_slope, and stp_lat_mean — ~13x faster because it avoids
#' computing the other ~200 condition columns.
#'
#' @param lookups Named list from prepare_lookup_tables().
#' @param tidx_with_birds Integer vector of timesteps to process (only those
#'   where the baseline has birds over radar).
#' @return data.table with columns: tidx, nr_birds, fl_inst_u10,
#'   stp_t2m_24hrs_slope, stp_lat_mean.
compute_predictor_values_fast <- function(lookups, tidx_with_birds) {
  birds_radar <- lookups$birds_radar
  cond_u10    <- lookups$cond_u10
  stp_slope   <- lookups$stp_slope
  stp_locs    <- lookups$stp_locs

  results <- lapply(tidx_with_birds, function(i) {
    br <- birds_radar[.(i)]
    bird_ids <- br$bird_ids[[1L]]

    if (is.null(bird_ids) || length(bird_ids) == 0L) {
      return(data.table(tidx = i, nr_birds = 0L,
                        fl_inst_u10 = NA_real_,
                        stp_t2m_24hrs_slope = NA_real_,
                        stp_lat_mean = NA_real_))
    }

    nr <- length(bird_ids)

    # fl_inst_u10: mean u10 across flight history (tidx <= i, skip 1st per bird)
    fl <- cond_u10[.(bird_ids)][tidx <= i]
    fl <- fl[, .SD[2:.N], by = bird]
    fl_u10 <- mean(fl$inst_u10, na.rm = TRUE)

    # stp_t2m_24hrs_slope: mean across birds' stopover conditions
    stp_vals <- stp_slope[.(bird_ids), stp_t2m_24hrs_slope, nomatch = NULL]
    stp_val  <- mean(stp_vals[is.finite(stp_vals)], na.rm = TRUE)

    # stp_lat_mean: mean latitude of all stopover locations of radar birds
    locs     <- stp_locs[.(bird_ids), y, nomatch = NULL]
    lat_mean <- mean(locs, na.rm = TRUE)

    data.table(tidx = i, nr_birds = nr,
               fl_inst_u10 = fl_u10,
               stp_t2m_24hrs_slope = stp_val,
               stp_lat_mean = lat_mean)
  })

  rbindlist(results)
}


# ============================================================
# Part 1: Predictor convergence analysis
# ============================================================

all_convergence <- list()
combo_idx <- 0L

for (radar in radars) {
  for (year in years) {
    for (season in seasons) {
      combo_idx <- combo_idx + 1L
      combo_label <- glue("{radar}_{year}_{season}")

      rds_out <- file.path(out_dir, "predictor_convergence",
                           glue("convergence_{combo_label}.RDS"))
      if (file.exists(rds_out)) {
        message(glue("[{combo_idx}/32] Skipping {combo_label} (already exists)"))
        all_convergence[[combo_label]] <- readRDS(rds_out)
        next
      }

      message(glue("[{combo_idx}/32] Loading data: {combo_label}"))

      # Load data once per combo (only migrating conditions needed —
      # we only compare timesteps where both have birds, so the non-migrating
      # fallback path is never used)
      birds        <- simulated_birds(radar, year, season)
      cond_mig     <- load_instantaneous_conditions(radar, year, season)
      cond_stp_mig <- load_stopover_conditions(radar, year, season)

      n_mig_birds <- length(unique(birds[migration == TRUE, bird]))

      # Only process timesteps where the baseline has migrating birds at radar.
      # These are the only timesteps that enter the correlation comparison.
      tidx_with_birds <- sort(birds[state == 1 & radar == TRUE & migration == TRUE,
                                    unique(tidx)])

      message(glue("  {length(tidx_with_birds)} timesteps with birds, ",
                    "{n_mig_birds} migrating birds"))

      # --- Baseline (100%) ---
      message("  Computing baseline (100%)...")
      lookups_base <- prepare_lookup_tables(birds, cond_mig, cond_stp_mig)
      baseline     <- compute_predictor_values_fast(lookups_base, tidx_with_birds)

      # --- Subsampled fractions ---
      combo_results <- list()

      for (frac in fractions) {
        for (rep in seq_len(n_repeats)) {
          seed <- as.integer(frac * 1000 + rep)
          message(glue("  fraction={frac}, repeat={rep}/{n_repeats}"))

          birds_sub    <- subsample_birds(birds, frac, seed)
          lookups_sub  <- prepare_lookup_tables(birds_sub, cond_mig, cond_stp_mig)
          sub_vals     <- compute_predictor_values_fast(lookups_sub, tidx_with_birds)

          # Compare against baseline on timesteps where both have birds > 0
          merged <- merge(baseline, sub_vals, by = "tidx", suffixes = c("_base", "_sub"))

          # Only compare timesteps where both have migrating birds over radar.
          # When nr_birds drops to 0, process_tidx falls back to non-migrating
          # conditions — a fundamentally different quantity, not a degraded estimate.
          both_have_birds <- merged$nr_birds_base > 0 & merged$nr_birds_sub > 0

          stats_list <- lapply(target_vars, function(v) {
            base_col <- paste0(v, "_base")
            sub_col  <- paste0(v, "_sub")

            valid <- both_have_birds &
                     !is.na(merged[[base_col]]) & !is.na(merged[[sub_col]])
            if (sum(valid) < 3L) {
              return(data.table(
                variable = v, n_valid = sum(valid),
                pearson_r = NA_real_, rmse = NA_real_, mae = NA_real_
              ))
            }

            bv <- merged[[base_col]][valid]
            sv <- merged[[sub_col]][valid]

            data.table(
              variable  = v,
              n_valid   = sum(valid),
              pearson_r = cor(bv, sv, use = "complete.obs"),
              rmse      = sqrt(mean((bv - sv)^2)),
              mae       = mean(abs(bv - sv))
            )
          })

          stats <- rbindlist(stats_list)
          # Count timesteps where baseline has birds but subsample lost them all
          n_base_has_birds <- sum(merged$nr_birds_base > 0)
          n_both_have      <- sum(both_have_birds)
          n_lost           <- n_base_has_birds - n_both_have

          stats[, `:=`(
            radar             = radar,
            year              = year,
            season            = season,
            fraction          = frac,
            repeat_id         = rep,
            seed              = seed,
            n_mig_birds       = n_mig_birds,
            n_mig_sampled     = max(1L, round(n_mig_birds * frac)),
            nr_birds_base_mean = mean(merged$nr_birds_base, na.rm = TRUE),
            nr_birds_sub_mean  = mean(merged$nr_birds_sub, na.rm = TRUE),
            n_timesteps_with_birds = n_base_has_birds,
            n_timesteps_lost       = n_lost
          )]

          combo_results[[length(combo_results) + 1L]] <- stats
        }
      }

      combo_dt <- rbindlist(combo_results)
      saveRDS(combo_dt, rds_out)
      all_convergence[[combo_label]] <- combo_dt
      message(glue("  Saved: {rds_out}"))

      # Free memory
      rm(birds, cond_mig, cond_stp_mig, lookups_base,
         baseline, combo_results, combo_dt)
      gc()
    }
  }
}

# Combine all results
convergence_summary <- rbindlist(all_convergence)
saveRDS(convergence_summary,
        file.path(out_dir, "predictor_convergence", "convergence_summary.RDS"))
message("Saved convergence_summary.RDS")


# ============================================================
# Part 1: Convergence plots
# ============================================================

convergence_summary <- readRDS(
  file.path(out_dir, "predictor_convergence", "convergence_summary.RDS")
)

# Compute starting population equivalent for x-axis annotation
# 150,000 starting migrants → n_mig_birds processed migrants per combo
convergence_summary[, starting_pop := round(150000 * fraction)]

# Variable labels
var_labels <- c(
  fl_inst_u10        = "En-route 10m u-wind\n(fl_inst_u10)",
  stp_t2m_24hrs_slope = "Stopover 24h temp. slope\n(stp_t2m_24hrs_slope)",
  stp_lat_mean       = "Mean stopover latitude\n(stp_lat_mean)"
)
convergence_summary[, variable_label := var_labels[variable]]

# Summarise across all combos per fraction × variable
plot_data <- convergence_summary[
  , .(
    r_mean   = mean(pearson_r, na.rm = TRUE),
    r_median = median(pearson_r, na.rm = TRUE),
    r_q05    = quantile(pearson_r, 0.05, na.rm = TRUE),
    r_q25    = quantile(pearson_r, 0.25, na.rm = TRUE),
    r_q75    = quantile(pearson_r, 0.75, na.rm = TRUE),
    r_q95    = quantile(pearson_r, 0.95, na.rm = TRUE)
  ),
  by = .(fraction, variable, variable_label)
]

# Add the 100% baseline point (r = 1.0 by definition)
baseline_rows <- unique(plot_data[, .(variable, variable_label)])
baseline_rows[, `:=`(fraction = 1.0, r_mean = 1, r_median = 1,
                      r_q05 = 1, r_q25 = 1, r_q75 = 1, r_q95 = 1)]
plot_data <- rbindlist(list(plot_data, baseline_rows), use.names = TRUE)

# Main convergence plot: all variables, aggregated across radar × year × season
p_convergence <- ggplot(plot_data, aes(x = fraction, color = variable_label,
                                       fill = variable_label)) +
  # geom_ribbon(aes(ymin = r_q05, ymax = r_q95), alpha = 0.15, color = NA) +
  geom_ribbon(aes(ymin = r_q25, ymax = r_q75), alpha = 0.3, color = NA) +
  geom_line(aes(y = r_median), linewidth = 1) +
  geom_point(aes(y = r_median), size = 2) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
  annotate("text", x = 0.02, y = 0.98, label = "r = 0.95",
           color = "grey40", hjust = 0, vjust = 1, size = 3) +
  scale_x_continuous(
    breaks = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 1.00),
    labels = scales::percent,
    sec.axis = sec_axis(
      ~ . * 150000,
      name = "Equivalent starting population",
      breaks = c(1500, 7500, 15000, 37500, 75000, 112500, 150000)
    )
  ) +
  scale_color_brewer(palette = "Set2", name = "Predictor variable") +
  scale_fill_brewer(palette = "Set2", name = "Predictor variable") +
  labs(
    x = "Proportion of original migrant population",
    y = "Pearson correlation with baseline (100%)",
    title = "Predictor stability vs. number of simulated migrants",
    subtitle = "Median with IQR band across all radar-year-season combinations (32 combos, 10 repeats each)"
  ) +
  coord_cartesian(ylim = c(0, 1.01)) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "plots", "predictor_convergence.pdf"),
       plot = p_convergence, width = 10, height = 7)
message("Saved predictor_convergence.pdf")

# Faceted version by radar × season
plot_data_faceted <- convergence_summary[
  , .(
    r_mean   = mean(pearson_r, na.rm = TRUE),
    r_median = median(pearson_r, na.rm = TRUE),
    r_q05    = quantile(pearson_r, 0.05, na.rm = TRUE),
    r_q25    = quantile(pearson_r, 0.25, na.rm = TRUE),
    r_q75    = quantile(pearson_r, 0.75, na.rm = TRUE),
    r_q95    = quantile(pearson_r, 0.95, na.rm = TRUE)
  ),
  by = .(fraction, variable, variable_label, radar, season)
]

radar_labels <- c(nlhrw = "Herwijnen", nldhl = "Den Helder")
plot_data_faceted[, radar_label  := radar_labels[radar]]
plot_data_faceted[, season_label := tools::toTitleCase(season)]

p_convergence_faceted <- ggplot(
  plot_data_faceted,
  aes(x = fraction, color = variable_label, fill = variable_label)
) +
  geom_ribbon(aes(ymin = r_q05, ymax = r_q95), alpha = 0.15, color = NA) +
  geom_ribbon(aes(ymin = r_q25, ymax = r_q75), alpha = 0.3, color = NA) +
  geom_line(aes(y = r_median), linewidth = 0.8) +
  geom_point(aes(y = r_median), size = 1.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey40") +
  facet_grid(radar_label ~ season_label) +
  scale_x_continuous(breaks = fractions, labels = scales::percent) +
  scale_color_brewer(palette = "Set2", name = "Predictor") +
  scale_fill_brewer(palette = "Set2", name = "Predictor") +
  labs(
    x = "Fraction of original migrant population",
    y = "Pearson correlation with baseline",
    title = "Predictor convergence by radar and season"
  ) +
  coord_cartesian(ylim = c(NA, 1.01)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "plots", "predictor_convergence_faceted.pdf"),
       plot = p_convergence_faceted, width = 12, height = 8)
message("Saved predictor_convergence_faceted.pdf")

# Identify stability threshold per variable
threshold_summary <- convergence_summary[
  , .(
    r_median = median(pearson_r, na.rm = TRUE),
    r_q05    = quantile(pearson_r, 0.05, na.rm = TRUE)
  ),
  by = .(fraction, variable)
][order(variable, fraction)]

message("\n=== Stability threshold (median r across all combos) ===")
for (v in target_vars) {
  sub <- threshold_summary[variable == v]
  above_95 <- sub[r_median >= 0.95, min(fraction)]
  above_95_q05 <- sub[r_q05 >= 0.95, min(fraction)]
  message(glue(
    "  {v}: median r >= 0.95 at {scales::percent(above_95)} ",
    "(= {scales::comma(above_95 * 150000)} starting migrants); ",
    "95th pct safe at {scales::percent(above_95_q05, accuracy = 1)}"
  ))
}


# ============================================================
# Part 2: Cohort proportion analysis
# ============================================================

message("\n=== Part 2: Cohort proportions ===")

cohort_results <- list()

for (radar in radars) {
  for (year in years) {
    for (season in seasons) {
      message(glue("  Cohort proportions: {radar} {year} {season}"))

      birds <- simulated_birds(radar, year, season)

      # Filter to migrating birds flying over radar
      over_radar <- birds[state == 1 & radar == TRUE]

      n_total <- length(unique(over_radar$bird))
      if (n_total == 0L) next

      # Count per cohort (pref_dir is the cohort identifier)
      counts <- over_radar[, .(n_birds = uniqueN(bird)), by = pref_dir]
      counts[, `:=`(
        radar      = radar,
        year       = year,
        season     = season,
        total      = n_total,
        proportion = n_birds / n_total
      )]

      # Map pref_dir to cohort names
      cohort_map <- rbindlist(lapply(names(bird_routes[[season]]), function(cn) {
        data.table(
          pref_dir   = bird_routes[[season]][[cn]]$endogenous_heading,
          cohort     = cn
        )
      }))
      counts <- merge(counts, cohort_map, by = "pref_dir", all.x = TRUE)

      cohort_results[[length(cohort_results) + 1L]] <- counts
      rm(birds, over_radar, counts)
    }
  }
}

cohort_dt <- rbindlist(cohort_results)
saveRDS(cohort_dt, file.path(out_dir, "cohort_proportions", "cohort_proportions.RDS"))
message("Saved cohort_proportions.RDS")

# Summary table
cohort_summary <- cohort_dt[
  , .(
    mean_n     = mean(n_birds),
    mean_total = mean(total),
    mean_prop  = mean(proportion),
    sd_prop    = sd(proportion)
  ),
  by = .(radar, season, cohort, pref_dir)
]
message("\nCohort proportion summary:")
print(cohort_summary)

# Verify proportions sum to 1
prop_check <- cohort_dt[, .(sum_prop = sum(proportion)), by = .(radar, year, season)]
stopifnot(all(abs(prop_check$sum_prop - 1.0) < 1e-10))
message("Proportion check passed (all sum to 1.0)")

# Cohort proportion plot
radar_labels <- c(nlhrw = "Herwijnen", nldhl = "Den Helder")
cohort_dt[, radar_label  := radar_labels[radar]]
cohort_dt[, season_label := tools::toTitleCase(season)]
cohort_dt[, cohort_label := paste0(cohort, " (", pref_dir, " deg)")]

p_cohort <- ggplot(cohort_dt, aes(x = factor(year), y = proportion,
                                   fill = cohort_label)) +
  geom_col(position = "stack", width = 0.7) +
  facet_grid(radar_label ~ season_label) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set1", name = "Cohort") +
  labs(
    x = "Year",
    y = "Proportion of birds over radar",
    title = "Cohort contributions to radar passages"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(out_dir, "plots", "cohort_proportions.pdf"),
       plot = p_cohort, width = 10, height = 7)
message("Saved cohort_proportions.pdf")

# ============================================================
# Part 2b: Per-timestep cohort counts & nr_birds relationship
# ============================================================

message("\n=== Part 2b: Cohort nr_birds distributions ===")

tidx_cohort_results <- list()

for (radar in radars) {
  for (year in years) {
    for (season in seasons) {
      message(glue("  Timestep cohort counts: {radar} {year} {season}"))

      birds <- simulated_birds(radar, year, season)
      over_radar <- birds[state == 1 & radar == TRUE]
      if (nrow(over_radar) == 0L) next

      # Per-timestep counts by cohort
      per_tidx <- over_radar[, .(n_cohort = uniqueN(bird)), by = .(tidx, pref_dir)]

      # Total per timestep
      totals <- over_radar[, .(nr_birds = uniqueN(bird)), by = tidx]
      per_tidx <- merge(per_tidx, totals, by = "tidx")

      # Map pref_dir to cohort names
      cohort_map <- rbindlist(lapply(names(bird_routes[[season]]), function(cn) {
        data.table(pref_dir = bird_routes[[season]][[cn]]$endogenous_heading,
                   cohort   = cn)
      }))
      per_tidx <- merge(per_tidx, cohort_map, by = "pref_dir", all.x = TRUE)

      per_tidx[, `:=`(radar = radar, year = year, season = season)]
      tidx_cohort_results[[length(tidx_cohort_results) + 1L]] <- per_tidx
      rm(birds, over_radar, per_tidx, totals)
    }
  }
}

tidx_cohort_dt <- rbindlist(tidx_cohort_results)
tidx_cohort_dt[, prop_cohort := n_cohort / nr_birds]
tidx_cohort_dt[, radar_label  := radar_labels[radar]]
tidx_cohort_dt[, season_label := tools::toTitleCase(season)]
tidx_cohort_dt[, cohort_label := paste0(cohort, " (", pref_dir, " deg)")]

saveRDS(tidx_cohort_dt, file.path(out_dir, "cohort_proportions", "cohort_nr_birds.RDS"))
message("Saved cohort_nr_birds.RDS")

# --- Plot A: Distribution of per-timestep bird counts by cohort ---
p_cohort_density <- ggplot(tidx_cohort_dt, aes(x = n_cohort, fill = cohort_label)) +
  geom_histogram(aes(y = after_stat(density)), position = "identity",
                 alpha = 0.5, binwidth = 1) +
  facet_grid(radar_label ~ season_label) +
  scale_x_continuous(limits = c(0, 50)) +
  scale_fill_brewer(palette = "Set1", name = "Cohort") +
  labs(
    x = "Birds per timestep (per cohort)",
    y = "Density",
    title = "Distribution of per-timestep bird counts by cohort"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "plots", "cohort_nr_birds_density.pdf"),
       plot = p_cohort_density, width = 10, height = 7)
message("Saved cohort_nr_birds_density.pdf")

# --- Plot B: Cohort proportion as a function of total nr_birds ---
# Bin total nr_birds and compute mean cohort proportions per bin
tidx_cohort_dt[, nr_birds_bin := cut(nr_birds,
  breaks = c(0, 2, 5, 10, 20, 50, Inf),
  labels = c("1-2", "3-5", "6-10", "11-20", "21-50", "50+")
)]

# Aggregate total cohort birds at each integer nr_birds value
nbirds_cohort <- tidx_cohort_dt[, .(
  total_cohort = sum(n_cohort)
), by = .(radar_label, season_label, cohort_label, nr_birds)]

# Clip x-axis at 99th percentile per facet; fill gaps; rolling sum to smooth
clip <- tidx_cohort_dt[, .(x_max = ceiling(quantile(nr_birds, 0.99))),
                        by = .(radar_label, season_label)]

filled <- nbirds_cohort[, {
  xmax <- clip[radar_label == .BY$radar_label & season_label == .BY$season_label, x_max]
  cohorts <- unique(cohort_label)
  grid <- CJ(nr_birds = 1:xmax, cohort_label = cohorts)
  vals <- .SD[grid, on = .(nr_birds, cohort_label)]
  vals[is.na(total_cohort), total_cohort := 0]
  vals
}, by = .(radar_label, season_label)]

setorder(filled, radar_label, season_label, cohort_label, nr_birds)
filled[, total_rollsum := zoo::rollsum(total_cohort, k = 21, fill = NA, align = "center"),
       by = .(radar_label, season_label, cohort_label)]
filled <- filled[!is.na(total_rollsum)]

p_cohort_by_nbirds <- ggplot(filled,
    aes(x = nr_birds, y = total_rollsum, fill = cohort_label)) +
  geom_area(position = "fill") +
  facet_grid(radar_label ~ season_label, scales = "free_x") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set1", name = "Cohort") +
  labs(
    x = "Total birds per timestep",
    y = "Cohort proportion",
    title = "Cohort composition by migration intensity",
    subtitle = "Proportion of birds from each cohort as a function of nr_birds"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "plots", "cohort_by_nr_birds.pdf"),
       plot = p_cohort_by_nbirds, width = 10, height = 7)
message("Saved cohort_by_nr_birds.pdf")

message("\n=== Sensitivity analysis complete ===")
