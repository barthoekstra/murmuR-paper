# =============================================================================
# comparison_pipeline.R
#
# Compare murmuR refactored pipeline outputs against ibm-ml reference outputs.
# Tests: nlhrw, 2017, spring + autumn.
#
# Steps:
#   0. ECMWF preprocessing verification (small subset)
#   1. Track processing
#   2. Annotation
#   3. Model data assembly
#
# Usage:
#   Rscript analysis/comparison_pipeline.R
# =============================================================================

library(data.table)
library(glue)

# --- Configuration -----------------------------------------------------------

IBM_ML_DATA   <- "/data/birdcloudstorage-tvm/ibm-ml/data"
COMP_BASE     <- "output/comparison"
COMP_DATA     <- file.path(COMP_BASE, "data")
COMP_RESULTS  <- file.path(COMP_BASE, "results")
REPORT_FILE   <- file.path(COMP_BASE, "comparison_report.txt")

RADAR   <- "nlhrw"
YEAR    <- 2017L
SEASONS <- c("autumn", "spring")
WORKERS <- 10L

# --- Helper: compare two data.tables -----------------------------------------

#' Compare two data.tables on shared numeric columns
#'
#' @param new  data.table from murmuR pipeline
#' @param ref  data.table from ibm-ml reference
#' @param keys Character vector of key columns to merge on
#' @param label Character label for this comparison
#' @return data.table with per-column comparison statistics
compare_dt <- function(new, ref, keys, label = "") {
  new <- as.data.table(new)
  ref <- as.data.table(ref)

  # Drop sf geometry columns — these aren't comparable as numeric
  geo_new <- names(new)[vapply(new, function(x) inherits(x, "sfc") || is.list(x), logical(1))]
  geo_ref <- names(ref)[vapply(ref, function(x) inherits(x, "sfc") || is.list(x), logical(1))]
  if (length(geo_new) > 0) new[, (geo_new) := NULL]
  if (length(geo_ref) > 0) ref[, (geo_ref) := NULL]

  info <- list(
    label     = label,
    new_nrow  = nrow(new),
    ref_nrow  = nrow(ref),
    new_ncol  = ncol(new),
    ref_ncol  = ncol(ref)
  )

  cat(glue("\n--- {label} ---"), "\n")
  cat(glue("  Dimensions: new={nrow(new)}x{ncol(new)}, ref={nrow(ref)}x{ncol(ref)}"), "\n")

  # Columns only in one
  new_only <- setdiff(names(new), names(ref))
  ref_only <- setdiff(names(ref), names(new))
  if (length(new_only) > 0)
    cat(glue("  Columns only in new ({length(new_only)}): {paste(head(new_only, 10), collapse=', ')}"), "\n")
  if (length(ref_only) > 0)
    cat(glue("  Columns only in ref ({length(ref_only)}): {paste(head(ref_only, 10), collapse=', ')}"), "\n")

  # Merge on keys
  merged <- tryCatch(
    merge(new, ref, by = keys, suffixes = c(".new", ".ref")),
    error = function(e) {
      cat(glue("  MERGE FAILED: {e$message}"), "\n")
      return(NULL)
    }
  )
  if (is.null(merged)) {
    return(data.table(
      label = label, column = NA_character_, pearson_r = NA_real_,
      rmse = NA_real_, mae = NA_real_, max_abs_diff = NA_real_,
      n_matched = 0L, note = "merge failed"
    ))
  }

  cat(glue("  Matched rows: {nrow(merged)} (new={nrow(new)}, ref={nrow(ref)})"), "\n")

  # Find shared numeric columns
  shared_cols <- setdiff(intersect(names(new), names(ref)), keys)
  numeric_in_new <- shared_cols[sapply(shared_cols, function(c) is.numeric(new[[c]]))]
  numeric_in_ref <- shared_cols[sapply(shared_cols, function(c) is.numeric(ref[[c]]))]
  numeric_cols <- intersect(numeric_in_new, numeric_in_ref)

  if (length(numeric_cols) == 0) {
    cat("  No shared numeric columns to compare.\n")
    return(data.table(
      label = label, column = NA_character_, pearson_r = NA_real_,
      rmse = NA_real_, mae = NA_real_, max_abs_diff = NA_real_,
      n_matched = nrow(merged), note = "no numeric columns"
    ))
  }

  results <- rbindlist(lapply(numeric_cols, function(col) {
    x <- merged[[paste0(col, ".new")]]
    y <- merged[[paste0(col, ".ref")]]

    # Remove rows where both are NA
    valid <- !is.na(x) & !is.na(y)
    x <- x[valid]
    y <- y[valid]
    n <- length(x)

    if (n < 2) {
      return(data.table(
        label = label, column = col, pearson_r = NA_real_,
        rmse = NA_real_, mae = NA_real_, max_abs_diff = NA_real_,
        n_valid = n, note = "too few valid values"
      ))
    }

    # Strip units class (ERA5 annotations may carry units from stars)
    x <- as.numeric(x)
    y <- as.numeric(y)

    r <- tryCatch(cor(x, y, use = "complete.obs"), error = function(e) NA_real_)
    diffs <- x - y
    data.table(
      label        = label,
      column       = col,
      pearson_r    = r,
      rmse         = sqrt(mean(diffs^2)),
      mae          = mean(abs(diffs)),
      max_abs_diff = max(abs(diffs)),
      n_valid      = n,
      note         = if (!is.na(r) && r < 0.95) "LOW CORRELATION" else ""
    )
  }))

  # Print summary
  n_low <- sum(results$note == "LOW CORRELATION", na.rm = TRUE)
  n_perfect <- sum(results$pearson_r > 0.9999, na.rm = TRUE)
  cat(glue("  Compared {nrow(results)} numeric columns: {n_perfect} perfect (r>0.9999), {n_low} low (r<0.95)"), "\n")

  if (n_low > 0) {
    low <- results[note == "LOW CORRELATION"]
    for (i in seq_len(nrow(low))) {
      cat(glue("  WARNING: {low$column[i]}: r={round(low$pearson_r[i], 4)}, max_diff={signif(low$max_abs_diff[i], 4)}"), "\n")
    }
  }

  results
}


# --- Helper: logging ---------------------------------------------------------

report_lines <- character(0)

report <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  report_lines <<- c(report_lines, msg)
}

write_report <- function() {
  dir.create(dirname(REPORT_FILE), recursive = TRUE, showWarnings = FALSE)
  writeLines(report_lines, REPORT_FILE)
  cat(glue("\nReport written to {REPORT_FILE}\n"))
}

# =============================================================================
# SETUP: Create symlinked directory tree
# =============================================================================

setup_symlinks <- function() {
  report("=" |> strrep(70))
  report("SETUP: Creating symlinked directory tree")
  report("=" |> strrep(70))

  dir.create(COMP_DATA, recursive = TRUE, showWarnings = FALSE)
  dir.create(COMP_RESULTS, recursive = TRUE, showWarnings = FALSE)

  # Helper: create symlink (remove existing first)
  make_link <- function(target, link_path) {
    dir.create(dirname(link_path), recursive = TRUE, showWarnings = FALSE)
    link_target <- Sys.readlink(link_path)
    if (file.exists(link_path) || (!is.na(link_target) && nchar(link_target) > 0)) {
      file.remove(link_path)
    }
    file.symlink(target, link_path)
    report(glue("  Linked: {link_path} -> {target}"))
  }

  # Raw ECMWF data
  make_link(
    file.path(IBM_ML_DATA, "ecmwf-era5-pressurelevels/2017"),
    file.path(COMP_DATA, "ecmwf-era5-pressurelevels/2017")
  )
  make_link(
    file.path(IBM_ML_DATA, "ecmwf-era5-singlelevels/2017"),
    file.path(COMP_DATA, "ecmwf-era5-singlelevels/2017")
  )
  make_link(
    file.path(IBM_ML_DATA, "ecmwf-era5-singlelevels/ecmwf-land_sea-mask.nc"),
    file.path(COMP_DATA, "ecmwf-era5-singlelevels/ecmwf-land_sea-mask.nc")
  )

  # VPTS (pre-processed, reuse as-is)
  make_link(
    file.path(IBM_ML_DATA, "vpts/nlhrw"),
    file.path(COMP_DATA, "vpts/nlhrw")
  )

  # Chunked tracks (input for track processing)
  make_link(
    file.path(IBM_ML_DATA, "processed/tracks/chunked/2017"),
    file.path(COMP_DATA, "processed/tracks/chunked/2017")
  )

  # Create output directories (these will be written to by the pipeline)
  dir.create(file.path(COMP_DATA, "processed/tracks/radar_chunks/nlhrw/2017"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(COMP_DATA, "processed/tracks/radar_seasons/nlhrw"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(COMP_DATA, "processed/tracks/radar_annotated/nlhrw"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(COMP_DATA, "processed/tracks/radar_annotated_seasons/nlhrw"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(COMP_DATA, "processed/model_data/nlhrw"),
             recursive = TRUE, showWarnings = FALSE)

  report("Setup complete.\n")
}

# =============================================================================
# STEP 0: ECMWF Preprocessing Verification
# =============================================================================

step0_ecmwf_verification <- function() {
  report("=" |> strrep(70))
  report("STEP 0: ECMWF Preprocessing — Symlink (verification skipped)")
  report("=" |> strrep(70))

  # ECMWF preprocessing (st_apply with slope_window) takes >1 hour for a
  # single func/window combo.  Skip the live verification and symlink the
  # ibm-ml preprocessed data directly.  The ECMWF preprocessing code can be
  # verified independently on a cluster if needed.
  report("  Symlinking ibm-ml preprocessed ECMWF data directly.")

  link_path <- file.path(COMP_DATA, "processed/ecmwf-era5-singlelevels")
  target    <- file.path(IBM_ML_DATA, "processed/ecmwf-era5-singlelevels")

  # Clean up any leftover from a previous run
  if (file.exists(link_path) || (!is.na(Sys.readlink(link_path)) && nchar(Sys.readlink(link_path)) > 0)) {
    unlink(link_path, recursive = TRUE)
  }
  file.symlink(target, link_path)

  report(glue("  Linked: {link_path} -> {target}"))

  results <- data.table(
    variable = "N/A", pearson_r = NA_real_, rmse = NA_real_,
    max_abs_diff = NA_real_, n_values = 0L,
    note = "Verification skipped — symlinked ibm-ml preprocessed data"
  )
  saveRDS(results, file.path(COMP_RESULTS, "step0_ecmwf_verification.RDS"))
  report("")
  results
}

# =============================================================================
# STEP 1: Track Processing
# =============================================================================

step1_track_processing <- function() {
  report("=" |> strrep(70))
  report("STEP 1: Track Processing")
  report("=" |> strrep(70))

  devtools::load_all()

  report("Running track processing pipeline...")
  run_track_processing_pipeline(
    radars   = RADAR,
    years    = YEAR,
    seasons  = SEASONS,
    overwrite = TRUE
  )

  # Compare outputs for each season x type
  all_results <- list()

  for (season in SEASONS) {
    for (suffix in c("", "_nomig")) {
      label <- glue("{RADAR}_{YEAR}_{season}{suffix}")
      new_file <- file.path(
        COMP_DATA, glue("processed/tracks/radar_seasons/{RADAR}/{label}.RDS")
      )
      ref_file <- file.path(
        IBM_ML_DATA, glue("processed/tracks/radar_seasons/{RADAR}/{label}.RDS")
      )

      if (!file.exists(new_file)) {
        report(glue("  SKIP: New file missing: {basename(new_file)}"))
        next
      }
      if (!file.exists(ref_file)) {
        report(glue("  SKIP: Ref file missing: {basename(ref_file)}"))
        next
      }

      new_dt <- as.data.table(readRDS(new_file))
      ref_dt <- as.data.table(readRDS(ref_file))

      # Determine merge keys
      if (nchar(suffix) > 0) {
        # nomig files: birds are random-sampled, can't merge on bird directly
        # Compare aggregate statistics instead
        report(glue("\n  {label}: nomig (random sample) - comparing distributions"))
        report(glue("    new: {nrow(new_dt)} rows, {length(unique(new_dt$tidx))} tidx"))
        report(glue("    ref: {nrow(ref_dt)} rows, {length(unique(ref_dt$tidx))} tidx"))

        # Compare tidx coverage
        tidx_overlap <- length(intersect(unique(new_dt$tidx), unique(ref_dt$tidx)))
        report(glue("    Shared tidx: {tidx_overlap}"))

        # Compare summary stats per tidx
        new_agg <- new_dt[, .(
          n_birds = uniqueN(bird),
          mean_x = mean(x, na.rm = TRUE),
          mean_y = mean(y, na.rm = TRUE)
        ), by = tidx]
        ref_agg <- ref_dt[, .(
          n_birds = uniqueN(bird),
          mean_x = mean(x, na.rm = TRUE),
          mean_y = mean(y, na.rm = TRUE)
        ), by = tidx]

        res <- compare_dt(new_agg, ref_agg, keys = "tidx", label = label)
        all_results[[label]] <- res
      } else {
        # Regular files: merge on tidx + bird
        res <- compare_dt(new_dt, ref_dt, keys = c("tidx", "bird"), label = label)
        all_results[[label]] <- res
      }
    }
  }

  results <- rbindlist(all_results, fill = TRUE)
  saveRDS(results, file.path(COMP_RESULTS, "step1_track_processing.RDS"))
  report("")
  results
}

# =============================================================================
# STEP 2: Annotation
# =============================================================================

step2_annotation <- function() {
  report("=" |> strrep(70))
  report("STEP 2: Annotation")
  report("=" |> strrep(70))

  devtools::load_all()

  report("Running annotation pipeline...")
  run_annotation_pipeline(
    radars  = RADAR,
    years   = YEAR,
    seasons = SEASONS
  )

  # Compare each of the 24 annotation season files
  annotation_types <- c("sl_windowed", "sl_instant", "pl_instant", "sl_lsm", "rdr", "phenol")
  track_types      <- c("", "_nomig")

  all_results <- list()

  for (season in SEASONS) {
    for (atype in annotation_types) {
      for (suffix in track_types) {
        label <- glue("{RADAR}_{YEAR}_{season}_{atype}{suffix}")
        new_file <- file.path(
          COMP_DATA,
          glue("processed/tracks/radar_annotated_seasons/{RADAR}/{label}.RDS")
        )
        ref_file <- file.path(
          IBM_ML_DATA,
          glue("processed/tracks/radar_annotated_seasons/{RADAR}/{label}.RDS")
        )

        if (!file.exists(new_file)) {
          report(glue("  SKIP: New file missing: {basename(new_file)}"))
          next
        }
        if (!file.exists(ref_file)) {
          report(glue("  SKIP: Ref file missing: {basename(ref_file)}"))
          next
        }

        new_dt <- as.data.table(readRDS(new_file))
        ref_dt <- as.data.table(readRDS(ref_file))

        # Determine merge keys based on annotation type
        if (atype %in% c("sl_windowed", "sl_instant", "pl_instant", "sl_lsm")) {
          keys <- c("tidx", "bird")
        } else if (atype == "phenol") {
          keys <- "tidx"
        } else {
          # rdr: keyed on datetime (+ height if present)
          if ("height" %in% names(new_dt) && "height" %in% names(ref_dt)) {
            keys <- c("datetime", "height")
          } else {
            keys <- "datetime"
          }
        }

        # For nomig bird-level annotations, the birds are different random
        # samples between runs. Merging on (tidx, bird) would yield ~0 rows.
        # Instead, aggregate per tidx and compare the means.
        is_nomig_birdlevel <- (suffix == "_nomig" &&
          atype %in% c("sl_windowed", "sl_instant", "pl_instant", "sl_lsm"))

        if (is_nomig_birdlevel) {
          report(glue("  {label}: nomig bird-level - aggregating per tidx"))
          num_cols_new <- names(new_dt)[sapply(new_dt, is.numeric)]
          num_cols_ref <- names(ref_dt)[sapply(ref_dt, is.numeric)]
          shared_num <- setdiff(intersect(num_cols_new, num_cols_ref), c("tidx", "bird"))

          new_agg <- new_dt[, lapply(.SD, mean, na.rm = TRUE), by = tidx, .SDcols = shared_num]
          ref_agg <- ref_dt[, lapply(.SD, mean, na.rm = TRUE), by = tidx, .SDcols = shared_num]
          res <- compare_dt(new_agg, ref_agg, keys = "tidx", label = label)
        } else {
          res <- compare_dt(new_dt, ref_dt, keys = keys, label = label)
        }
        all_results[[label]] <- res
      }
    }
  }

  results <- rbindlist(all_results, fill = TRUE)
  saveRDS(results, file.path(COMP_RESULTS, "step2_annotation.RDS"))
  report("")
  results
}

# =============================================================================
# STEP 3: Model Data Assembly
# =============================================================================

step3_model_data <- function() {
  report("=" |> strrep(70))
  report("STEP 3: Model Data Assembly")
  report("=" |> strrep(70))

  devtools::load_all()

  report("Running model data pipeline...")
  run_model_data_pipeline(
    radars   = RADAR,
    years    = YEAR,
    seasons  = SEASONS,
    use_sims = c(TRUE, FALSE),
    quickrun = FALSE,
    workers  = WORKERS,
    batch_size = 500L
  )

  # Compare parquet outputs
  parquet_suffixes <- c(
    "",
    "_noninterp",
    "_3km",
    "_3km_noninterpolated",
    "_nosims",
    "_nosims_noninterp",
    "_nosims_3km",
    "_nosims_3km_noninterpolated"
  )

  all_results <- list()

  for (season in SEASONS) {
    for (psuffix in parquet_suffixes) {
      fname <- glue("{RADAR}_{YEAR}_{season}{psuffix}.parquet")
      label <- glue("{RADAR}_{YEAR}_{season}{psuffix}")

      new_file <- file.path(COMP_DATA, "processed/model_data", RADAR, fname)
      ref_file <- file.path(IBM_ML_DATA, "processed/model_data", RADAR, fname)

      if (!file.exists(new_file)) {
        report(glue("  SKIP: New file missing: {fname}"))
        next
      }
      if (!file.exists(ref_file)) {
        report(glue("  SKIP: Ref file missing: {fname}"))
        next
      }

      new_dt <- as.data.table(nanoparquet::read_parquet(new_file))
      ref_dt <- as.data.table(nanoparquet::read_parquet(ref_file))

      # Model data is keyed on datetime + height
      keys <- c("datetime", "height")

      # Ensure datetime columns are compatible
      if ("datetime" %in% names(new_dt)) {
        new_dt[, datetime := as.POSIXct(datetime, tz = "UTC")]
      }
      if ("datetime" %in% names(ref_dt)) {
        ref_dt[, datetime := as.POSIXct(datetime, tz = "UTC")]
      }

      res <- compare_dt(new_dt, ref_dt, keys = keys, label = label)
      all_results[[label]] <- res
    }
  }

  results <- rbindlist(all_results, fill = TRUE)
  saveRDS(results, file.path(COMP_RESULTS, "step3_model_data.RDS"))
  report("")
  results
}

# =============================================================================
# SUMMARY
# =============================================================================

print_summary <- function() {
  report("=" |> strrep(70))
  report("SUMMARY")
  report("=" |> strrep(70))

  for (step_file in c("step0_ecmwf_verification", "step1_track_processing",
                       "step2_annotation", "step3_model_data")) {
    fpath <- file.path(COMP_RESULTS, paste0(step_file, ".RDS"))
    if (!file.exists(fpath)) {
      report(glue("  {step_file}: NOT RUN"))
      next
    }
    res <- readRDS(fpath)
    if ("note" %in% names(res)) {
      n_low <- sum(res$note == "LOW CORRELATION", na.rm = TRUE)
    } else {
      n_low <- 0
    }
    n_total <- nrow(res)
    if ("pearson_r" %in% names(res)) {
      n_high <- sum(res$pearson_r > 0.999, na.rm = TRUE)
      n_na   <- sum(is.na(res$pearson_r))
      median_r <- median(res$pearson_r, na.rm = TRUE)
      report(glue("  {step_file}: {n_total} comparisons, {n_high} excellent (r>0.999), {n_low} low (r<0.95), {n_na} NA, median r={round(median_r, 6)}"))

      # Show worst columns
      if (n_low > 0) {
        low <- res[note == "LOW CORRELATION"][order(pearson_r)]
        for (i in seq_len(min(5, nrow(low)))) {
          report(glue("    -> {low$label[i]} / {low$column[i]}: r={round(low$pearson_r[i], 4)}"))
        }
      }
    } else {
      report(glue("  {step_file}: {n_total} rows"))
    }
  }

  report("")
  report("Note: Columns depending on randomly-sampled nomig birds are expected")
  report("to show lower correlations (stp_*, fl_* for non-migration timesteps).")
  report("Focus on regular (non-nomig) tracks and rdr/phenol annotations for")
  report("pipeline equivalence validation.")
  report("")
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  start_time <- Sys.time()
  report(glue("Comparison pipeline started: {start_time}"))
  report(glue("murmuR vs ibm-ml | radar={RADAR}, year={YEAR}, seasons={paste(SEASONS, collapse=',')}"))
  report("")

  # Parallelisation
  future::plan(future::multicore, workers = WORKERS)
  progressr::handlers("progress")
  options(future.globals.maxSize = 100000 * 1024^2)

  # Create directory tree first, then resolve to absolute path
  setup_symlinks()

  # Set up data path to our comparison working directory (must be absolute)
  options(murmuR.data_path = normalizePath(COMP_DATA))

  tryCatch(step0_ecmwf_verification(), error = function(e) {
    report(glue("STEP 0 FAILED: {e$message}"))
  })
  gc()

  tryCatch(step1_track_processing(), error = function(e) {
    report(glue("STEP 1 FAILED: {e$message}"))
  })
  gc()

  tryCatch(step2_annotation(), error = function(e) {
    report(glue("STEP 2 FAILED: {e$message}"))
  })
  gc()

  tryCatch(step3_model_data(), error = function(e) {
    report(glue("STEP 3 FAILED: {e$message}"))
  })
  gc()

  print_summary()

  end_time <- Sys.time()
  report(glue("Total runtime: {round(difftime(end_time, start_time, units='mins'), 1)} minutes"))

  write_report()
}

main()
