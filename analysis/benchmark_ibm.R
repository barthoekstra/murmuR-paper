# ============================================================
# Benchmark: IBM simulation runtime
# ============================================================
#
# Simulates 1000 migrants (100 per replicate × 10 replicates)
# across 10 cores for 7 days (168 hourly timesteps).

library(murmuR)
library(future)

era5_file <- "/data/birdcloudstorage-tvm/ibm-ml/data/ecmwf-era5-pressurelevels/2020/202008.nc"

br <- bird_routes[["autumn"]][["cohort1"]]
cohort <- list(
  name = "cohort1",
  season = "autumn",
  endogenous_heading = br$endogenous_heading,
  target_lat = br$target_pos$target_lat,
  target_lon = br$target_pos$target_lon
)

n_runs   <- 5
nbirds   <- 100
nsims    <- 10
workers  <- 10
timestep_range <- 1:168  # 7 days at hourly resolution

plan("multisession", workers = workers)

timings <- numeric(n_runs)
for (i in seq_len(n_runs)) {
  cat(sprintf("\n=== Run %d/%d ===\n", i, n_runs))
  t0 <- Sys.time()
  out <- run_ibm(
    era5_files = era5_file,
    cohort = cohort,
    nbirds = nbirds,
    nsims = nsims,
    timesteps = timestep_range,
    randomseed = 42 + i
  )
  t1 <- Sys.time()
  timings[i] <- as.numeric(difftime(t1, t0, units = "secs"))
  cat(sprintf("Run %d: %.1f seconds (%d rows)\n", i, timings[i], nrow(out)))
  rm(out); gc(verbose = FALSE)
}

cat("\n=== Benchmark summary ===\n")
cat(sprintf("Runs:   %d\n", length(timings)))
cat(sprintf("Mean:   %.1f s\n", mean(timings)))
cat(sprintf("Median: %.1f s\n", median(timings)))
cat(sprintf("Min:    %.1f s\n", min(timings)))
cat(sprintf("Max:    %.1f s\n", max(timings)))
cat(sprintf("SD:     %.1f s\n", sd(timings)))
cat(sprintf("\nConfig: %d birds (%d x %d replicates), %d workers, %d timesteps (%d days)\n",
            nbirds * nsims, nbirds, nsims, workers,
            length(timestep_range), length(timestep_range) / 24))

# Save results
results <- list(
  timings_s = timings,
  config = list(
    nbirds = nbirds,
    nsims = nsims,
    total_birds = nbirds * nsims,
    workers = workers,
    timesteps = length(timestep_range),
    days = length(timestep_range) / 24,
    era5_file = era5_file,
    cohort = cohort,
    date = Sys.time()
  ),
  summary = list(
    mean_s   = mean(timings),
    median_s = median(timings),
    min_s    = min(timings),
    max_s    = max(timings),
    sd_s     = sd(timings)
  )
)

saveRDS(results, "output/benchmark_ibm.RDS")
cat("\nSaved output/benchmark_ibm.RDS\n")
