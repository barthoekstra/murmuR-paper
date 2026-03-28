# ============================================================
# Groundspeed summary of simulated migrants
# ============================================================

library(murmuR)
library(data.table)

options(murmuR.data_path = "/data/birdcloudstorage-tvm/ibm-ml/data/")

# Sample across radars, years, and seasons
combos <- CJ(
  radar  = names(radars_named),
  year   = 2017:2024,
  season = c("autumn", "spring")
)

results <- lapply(seq_len(nrow(combos)), function(i) {
  r <- combos[i]
  birds <- simulated_birds(r$radar, r$year, r$season)
  flying <- birds[state == 1]

  data.table(
    radar     = r$radar,
    year      = r$year,
    season    = r$season,
    n_birds   = uniqueN(flying$bird),
    n_rows    = nrow(flying),
    gs_mean   = mean(flying$groundspeed, na.rm = TRUE),
    gs_median = median(flying$groundspeed, na.rm = TRUE),
    gs_sd     = sd(flying$groundspeed, na.rm = TRUE),
    gs_q05    = quantile(flying$groundspeed, 0.05, na.rm = TRUE),
    gs_q95    = quantile(flying$groundspeed, 0.95, na.rm = TRUE)
  )
})
dt <- rbindlist(results)

cat("=== Per-combo groundspeed summary (m/s) ===\n\n")
print(dt, digits = 2)

cat("\n=== By season ===\n")
season_summary <- dt[, .(
  gs_mean   = mean(gs_mean),
  gs_median = mean(gs_median),
  gs_q05    = mean(gs_q05),
  gs_q95    = mean(gs_q95)
), by = season]
print(season_summary, digits = 2)

cat("\n=== Overall (m/s and km/h) ===\n")
cat(sprintf("Mean:   %.1f m/s (%.1f km/h)\n", mean(dt$gs_mean), mean(dt$gs_mean) * 3.6))
cat(sprintf("Median: %.1f m/s (%.1f km/h)\n", mean(dt$gs_median), mean(dt$gs_median) * 3.6))
cat(sprintf("5th-95th pctile: %.1f-%.1f m/s (%.0f-%.0f km/h)\n",
            mean(dt$gs_q05), mean(dt$gs_q95),
            mean(dt$gs_q05) * 3.6, mean(dt$gs_q95) * 3.6))
