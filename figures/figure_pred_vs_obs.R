library(tidyverse)
library(data.table)
library(ggpointdensity)
library(ggdensity)
library(scattermore)
library(ggh4x)

library(murmuR)
options(murmuR.data_path = "/data/birdcloudstorage-tvm/ibm-ml/data/")

out_path <- function(...) file.path("/data/birdcloudstorage-tvm/murmuR/output", ...)

# --- Load predictions for dynamic sampling ---
load_preds <- function(sampling) {
  preds_list <- list()
  for (radar in names(radars_named)) {
    for (season in names(months_seasons)) {
      pattern <- sprintf(
        murmuR_path("models/predictions/regr.rmse_%s/%s/%s_%s_*_%s_preds.RDS"),
        sampling, radar, radar, season, sampling
      )
      files <- Sys.glob(pattern)
      for (f in files) {
        dt <- readRDS(f)[, .(year, height, truth, response, datetime)]
        dt[, `:=`(radar = radar, season = season, sampling = sampling)]
        preds_list[[length(preds_list) + 1]] <- dt
      }
    }
  }
  rbindlist(preds_list)
}

preds <- rbindlist(lapply(c("dynamic", "radar"), load_preds))

# Human-readable labels
preds[, radar_human := fifelse(radar == "nlhrw", "Herwijnen", "Den Helder")]
preds[, season_human := str_to_sentence(season)]
preds[, sampling_human := fifelse(sampling == "dynamic", "Dynamic", "Radar only")]

preds[, radar_human := factor(radar_human, levels = c("Den Helder", "Herwijnen"))]
preds[, season_human := factor(season_human, levels = c("Spring", "Autumn"))]
preds[, sampling_human := factor(sampling_human, levels = c("Dynamic", "Radar only"))]

# --- Identify top 10 peak nights (by observed VID) per radar × season × year ---
# Peak nights are defined on truth, so consistent across samplings
night_ranks <- preds[sampling == "dynamic",
                     .(nightly_vid = sum(truth) * 0.05),
                     by = .(radar, season, year, datetime)]
night_ranks[, date := as.Date(datetime)]
night_ranks <- night_ranks[, .(nightly_vid = sum(nightly_vid)),
                           by = .(radar, season, year, date)]
night_ranks[, vid_rank := rank(-nightly_vid, ties.method = "min"),
            by = .(radar, season, year)]
night_ranks[, is_peak := vid_rank <= 10]

# Merge peak flag back
preds[, date := as.Date(datetime)]
preds <- merge(preds, night_ranks[, .(radar, season, year, date, is_peak)],
               by = c("radar", "season", "year", "date"), all.x = TRUE)

cbrt_trans <- scales::new_transform(
  "cbrt",
  transform = function(x) sign(x) * abs(x)^(1/3),
  inverse   = function(x) sign(x) * abs(x)^3
)

set.seed(42)
preds_sub <- preds[truth <= 500 & response <= 500]
preds_dynamic <- preds_sub[sampling == "dynamic"][
  , .SD[sample(.N, min(.N, 100000))], by = .(radar_human, season_human)]

x_breaks <- c(0, 50, 100, 200, 300)
y_breaks <- c(0, 50, 100, 200, 300, 400, 500)

ggplot(preds_dynamic, aes(x = response, y = truth, color = truth)) +
  geom_scattermore(pointsize = 8, alpha = 0.5, pixels = c(2048, 2048), interpolate = TRUE) +
  scale_color_viridis_c(
    option = "inferno", trans = cbrt_trans,
    breaks = c(0, 50, 100, 200, 400),
    name = expression("Observed density [" * birds ~ km^{-3} * "]")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = 3, color = "white", linewidth = 1) +
  geom_smooth(method = "lm", color = "red", linewidth = 0.8) +
  facet_grid(rows = vars(radar_human), cols = vars(season_human)) +
  scale_x_continuous(transform = "sqrt", breaks = x_breaks, labels = x_breaks) +
  scale_y_continuous(transform = "sqrt", breaks = y_breaks, labels = y_breaks) +
  coord_equal(expand = TRUE, clip = "off") +
  labs(
    x = expression("Predicted density [" * birds ~ km^{-3} * "]"),
    y = expression("Observed density [" * birds ~ km^{-3} * "]")
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.98),
    legend.justification = c(1, 1),
    legend.direction = "horizontal"
  ) +
  guides(color = guide_colorbar(barwidth = 15, barheight = 1.2,
                                title.position = "top", title.hjust = 0))

ggsave(
  out_path("models/plots/pred_vs_obs/pred_vs_obs_dynamic.pdf"),
  width = 10, height = 8, create.dir = TRUE
)

