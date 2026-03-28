library(tidyverse)
library(mlr3verse)
library(mlr3db)
library(ggdist)
library(ggsci)
library(patchwork)
library(ggstats)
library(ggh4x)
library(ggnewscale)

library(murmuR)
options(murmuR.data_path = "/data/birdcloudstorage-tvm/ibm-ml-refactor/data/")

task <- create_task(radar = "nlhrw", season = "autumn", sampling = "dynamic")
orig_features <- task$feature_names

importances_nlhrw <- Sys.glob(murmuR_path("models/predictions/regr.rmse_dynamic/nlhrw/*_importances.RDS"))
importances_nldhl <- Sys.glob(murmuR_path("models/predictions/regr.rmse_dynamic/nldhl/*_importances.RDS"))
importances <- c(importances_nlhrw, importances_nldhl)

vi <- importances %>%
  lapply(readRDS) %>%
  bind_rows() %>%
  pivot_longer(-c(radar, season, sampling, year), names_to = "variable", values_to = "importance") %>%
  drop_na() %>%
  mutate(
    radar = case_when(
      radar == "nlhrw" ~ "Herwijnen",
      radar == "nldhl" ~ "Den Helder",
    ),
    radar = fct_relevel(as.factor(radar), c("Herwijnen", "Den Helder")),
    season = case_when(
      season == "autumn" ~ "Autumn",
      season == "spring" ~ "Spring"
    ),
    season = factor(season, levels = c("Spring", "Autumn")),
    class = case_when(
      str_detect(variable, "rdr_") ~ "Radar-level",
      str_detect(variable, "fl_") ~ "En-route",
      str_detect(variable, "stp_") ~ "Stopover",
      str_detect(variable, "yday|solar_elev|solar_azim|solar_elevation_change") ~ "Phenology",
      variable == "height" ~ "Radar-level",
      variable == "nr_birds" ~ "Radar-level"
    ),
    class = fct_rev(factor(class, levels = c("Radar-level", "En-route", "Stopover", "Phenology"))),
    envir_condition = case_when(
      str_detect(variable, "_u|_v") ~ "Wind",
      str_detect(variable, "_t|_t2m") ~ "Temperature",
      str_detect(variable, "_lsrr|_crr") ~ "Precipitation",
      str_detect(variable, "_msl") ~ "Pressure",
      str_detect(variable, "_blh") ~ "Boundary layer height",
      str_detect(variable, "_cc|_lcc|_mcc|_hcc|_cbh") ~ "Cloud cover",
      str_detect(variable, "_r|_q") ~ "Humidity",
      str_detect(variable, "_lsm") ~ "Land-sea mask",
      variable == "height" ~ "Height",
      variable == "yday" ~ "Day of year",
      str_detect(variable, "solar_") ~ "Time of night",
      str_detect(variable, "_lat|_lon") ~ "Take-off location",
      variable == "nr_birds" ~ "Nr. of simulated migrants",
      .default = NA
    ),
    period = case_when(
      str_detect(variable, "_3hrs") ~ "3 hours",
      str_detect(variable, "_24hrs") ~ "24 hours",
      str_detect(variable, "_72hrs") ~ "3 days",
      str_detect(variable, "_168hrs") ~ "1 week",
      .default = "Instantaneous"
    ),
    importance = importance * 100
  )

# Aggregated importance of certain metric classes
observable_colors <- pal_observable(palette = c("observable10"), alpha = 1)(8)
observable_colors <- c(observable_colors, "#000000")
observable_four <- rev(pal_observable(palette = c("observable10"), alpha = 1)(8)[c(1, 3, 4, 2)])

vi %>%
  filter(variable != "height") %>%
  group_by(radar, class, season, year) %>%
  summarise(summed_importance = sum(importance)) %>%
  ggplot(aes(x = summed_importance, y = class, color = class, group = radar, shape = radar)) +
  stat_pointinterval(point_size = 4, position = "dodge") +
  geom_stripped_rows() +
  scale_color_manual(name = "Class", guide = "none", values = observable_four) +
  # scale_color_observable(name = "Class", guide = "none") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_shape(name = "Radar") +
  facet_grid2(cols = vars(season)) +
  theme_bw(base_size = 14) +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  labs(x = "Summed variable importance [%]", y = "") -> plot_vi_overall

vi %>%
  filter(variable != "height") %>%
  mutate(class = fct_rev(class)) %>%
  group_by(radar, class, season, year, envir_condition) %>%
  mutate(avg_importance = mean(importance, na.rm = TRUE)) %>%
  ggplot(aes(x = avg_importance, y = envir_condition, color = class, group = radar, shape = radar)) +
  stat_pointinterval(point_size = 4, position = "dodge") +
  geom_stripped_rows() +
  # scale_color_observable(name = "Class") +
  scale_color_manual(name = "Class", values = rev(observable_four)) +
  scale_alpha_manual(values = c(1, 0)) +
  scale_shape(name = "Radar") +
  facet_grid2(rows = vars(class), cols = vars(season), scales = "free", space = "free_y") +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0)) +
  labs(x = "Average variable importance [%]", y = "") -> plot_vi_variables

vi %>%
  filter(variable != "height", class == "Stopover") %>%
  mutate(period = fct_rev(factor(period, levels = c("Instantaneous", "3 hours", "24 hours", "3 days", "1 week")))) %>%
  group_by(radar, season, year, period) %>%
  mutate(avg_importance = mean(importance, na.rm = TRUE)) %>%
  ggplot(aes(x = avg_importance, y = period, color = class, group = radar, shape = radar)) +
  stat_pointinterval(point_size = 4, position = "dodge") +
  geom_stripped_rows() +
  scale_color_manual(values = observable_colors[4], guide = "none", name = "Class") +
  scale_shape(name = "Radar") +
  facet_grid2(rows = vars(class), cols = vars(season), scales = "free", space = "free_y") +
  theme_bw(base_size = 14) +
  theme(strip.text.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0)) +
  labs(x = "Average variable importance [%]", y = "") -> plot_vi_temporal

plot_vi_overall /
  plot_vi_variables /
  plot_vi_temporal + plot_layout(guides = "collect", axes = "collect", axis_titles = "collect", heights = c(0.3,2,0.3)) +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position = "top")

ggsave(murmuR_path("models/plots/varimp/varimp.pdf"), width = 12, height = 15)
