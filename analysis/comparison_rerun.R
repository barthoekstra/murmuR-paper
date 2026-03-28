# Re-run Steps 2 (comparison only) and Step 3 (pipeline + comparison)
# after fixing:
#   - compare_dt: strip units class before numeric comparison
#   - write_model_data: DT[interpolated == FALSE] instead of DT[!interpolated]

library(data.table)
library(glue)

IBM_ML_DATA   <- "/data/birdcloudstorage-tvm/ibm-ml/data"
COMP_BASE     <- "output/comparison"
COMP_DATA     <- file.path(COMP_BASE, "data")
COMP_RESULTS  <- file.path(COMP_BASE, "results")
REPORT_FILE   <- file.path(COMP_BASE, "comparison_report_rerun.txt")

RADAR   <- "nlhrw"
YEAR    <- 2017L
SEASONS <- c("autumn", "spring")
WORKERS <- 10L

# --- Helper: compare two data.tables -----------------------------------------

compare_dt <- function(new, ref, keys, label = "") {
  new <- as.data.table(new)
  ref <- as.data.table(ref)

  # Drop sf geometry columns
  geo_new <- names(new)[vapply(new, function(x) inherits(x, "sfc") || is.list(x), logical(1))]
  geo_ref <- names(ref)[vapply(ref, function(x) inherits(x, "sfc") || is.list(x), logical(1))]
  if (length(geo_new) > 0) new[, (geo_new) := NULL]
  if (length(geo_ref) > 0) ref[, (geo_ref) := NULL]

  cat(glue("\n--- {label} ---"), "\n")
  cat(glue("  Dimensions: new={nrow(new)}x{ncol(new)}, ref={nrow(ref)}x{ncol(ref)}"), "\n")

  new_only <- setdiff(names(new), names(ref))
  ref_only <- setdiff(names(ref), names(new))
  if (length(new_only) > 0)
    cat(glue("  Columns only in new ({length(new_only)}): {paste(head(new_only, 10), collapse=', ')}"), "\n")
  if (length(ref_only) > 0)
    cat(glue("  Columns only in ref ({length(ref_only)}): {paste(head(ref_only, 10), collapse=', ')}"), "\n")

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

report_lines <- character(0)
report <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  report_lines <<- c(report_lines, msg)
}

# =============================================================================
# STEP 2: Annotation comparison only (files already generated)
# =============================================================================

step2_comparison <- function() {
  report("=" |> strrep(70))
  report("STEP 2: Annotation Comparison (files already exist)")
  report("=" |> strrep(70))

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

        if (!file.exists(new_file)) { report(glue("  SKIP: {basename(new_file)}")); next }
        if (!file.exists(ref_file)) { report(glue("  SKIP ref: {basename(ref_file)}")); next }

        new_dt <- as.data.table(readRDS(new_file))
        ref_dt <- as.data.table(readRDS(ref_file))

        if (atype %in% c("sl_windowed", "sl_instant", "pl_instant", "sl_lsm")) {
          keys <- c("tidx", "bird")
        } else if (atype == "phenol") {
          keys <- "tidx"
        } else {
          if ("height" %in% names(new_dt) && "height" %in% names(ref_dt)) {
            keys <- c("datetime", "height")
          } else {
            keys <- "datetime"
          }
        }

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

        # Free memory
        rm(new_dt, ref_dt)
        gc(verbose = FALSE)
      }
    }
  }

  results <- rbindlist(all_results, fill = TRUE)
  saveRDS(results, file.path(COMP_RESULTS, "step2_annotation.RDS"))
  report("")
  results
}

# =============================================================================
# STEP 3: Model Data Assembly (re-run pipeline + comparison)
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

  gc()

  parquet_suffixes <- c(
    "", "_noninterp", "_3km", "_3km_noninterpolated",
    "_nosims", "_nosims_noninterp", "_nosims_3km", "_nosims_3km_noninterpolated"
  )
  all_results <- list()

  for (season in SEASONS) {
    for (psuffix in parquet_suffixes) {
      fname <- glue("{RADAR}_{YEAR}_{season}{psuffix}.parquet")
      label <- glue("{RADAR}_{YEAR}_{season}{psuffix}")

      new_file <- file.path(COMP_DATA, "processed/model_data", RADAR, fname)
      ref_file <- file.path(IBM_ML_DATA, "processed/model_data", RADAR, fname)

      if (!file.exists(new_file)) { report(glue("  SKIP: {fname}")); next }
      if (!file.exists(ref_file)) { report(glue("  SKIP ref: {fname}")); next }

      new_dt <- as.data.table(nanoparquet::read_parquet(new_file))
      ref_dt <- as.data.table(nanoparquet::read_parquet(ref_file))

      keys <- c("datetime", "height")
      if ("datetime" %in% names(new_dt))
        new_dt[, datetime := as.POSIXct(datetime, tz = "UTC")]
      if ("datetime" %in% names(ref_dt))
        ref_dt[, datetime := as.POSIXct(datetime, tz = "UTC")]

      res <- compare_dt(new_dt, ref_dt, keys = keys, label = label)
      all_results[[label]] <- res

      rm(new_dt, ref_dt)
      gc(verbose = FALSE)
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

  for (step_file in c("step1_track_processing", "step2_annotation", "step3_model_data")) {
    fpath <- file.path(COMP_RESULTS, paste0(step_file, ".RDS"))
    if (!file.exists(fpath)) { report(glue("  {step_file}: NOT RUN")); next }
    res <- readRDS(fpath)
    n_low <- if ("note" %in% names(res)) sum(res$note == "LOW CORRELATION", na.rm = TRUE) else 0
    n_total <- nrow(res)
    if ("pearson_r" %in% names(res)) {
      n_high <- sum(res$pearson_r > 0.999, na.rm = TRUE)
      n_na   <- sum(is.na(res$pearson_r))
      median_r <- median(res$pearson_r, na.rm = TRUE)
      report(glue("  {step_file}: {n_total} comparisons, {n_high} excellent (r>0.999), {n_low} low (r<0.95), {n_na} NA, median r={round(median_r, 6)}"))
      if (n_low > 0) {
        low <- res[note == "LOW CORRELATION"][order(pearson_r)]
        for (i in seq_len(min(10, nrow(low)))) {
          report(glue("    -> {low$label[i]} / {low$column[i]}: r={round(low$pearson_r[i], 4)}"))
        }
      }
    }
  }

  report("")
  report("Note: nomig columns expected to show lower correlations (random sampling).")
  report("")
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  start_time <- Sys.time()
  report(glue("Re-run started: {start_time}"))

  future::plan(future::multicore, workers = WORKERS)
  progressr::handlers("progress")
  options(future.globals.maxSize = 100000 * 1024^2)
  options(murmuR.data_path = normalizePath(COMP_DATA))

  tryCatch(step2_comparison(), error = function(e) {
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

  dir.create(dirname(REPORT_FILE), recursive = TRUE, showWarnings = FALSE)
  writeLines(report_lines, REPORT_FILE)
  cat(glue("\nReport written to {REPORT_FILE}\n"))
}

main()
