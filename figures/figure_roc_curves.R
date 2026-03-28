library(tidyverse)
library(data.table)
library(ggpp)
library(ggsci)
library(suncalc)
library(plotROC)
library(patchwork)

library(murmuR)
options(murmuR.data_path = "/data/birdcloudstorage-tvm/ibm-ml-refactor/data/")

preds_autumn <- Sys.glob(murmuR_path("models/predictions/regr.rmse_dynamic/nlhrw/nlhrw_autumn_*.RDS")) %>%
  lapply(readRDS) %>%
  bind_rows() %>%
  mutate(season = "autumn")

preds_spring <- Sys.glob(murmuR_path("models/predictions/regr.rmse_dynamic/nlhrw/nlhrw_spring_*.RDS")) %>%
  lapply(readRDS) %>%
  bind_rows() %>%
  mutate(season = "spring")

# Threshold identification
plot_roc_curves <- function(preds, height_cutoff = 3000) {
  dt <- preds
  dt <- dt[height <= height_cutoff]
  dt[, date := as.Date(datetime)]

  dt_nights <- dt %>%
    select(date) %>%
    distinct() %>%
    mutate(
      lat = radars_named$nlhrw$lat,
      lon = radars_named$nlhrw$lon,
      nextday = date + 1
    ) %>% as.data.frame()

  sunset_times <- getSunlightTimes(data = dt_nights, keep = "sunset", tz = "UTC")
  sunrise_times <- getSunlightTimes(data = dt_nights %>% mutate(date = date + 1), keep = "sunrise", tz = "UTC") %>%
    mutate(date = date - 1)
  night_times <- left_join(sunset_times, sunrise_times) %>%
    mutate(hrs = difftime(sunrise, sunset))

  dt_agg_night <- dt[, .(dens250 = sum(dens, na.rm = TRUE) * 0.05),
                     by = .(year, night_num, date)]
  dt_agg_night <- merge(dt_agg_night, night_times, by = "date")

  dt_agg_hourly <- dt[, .(dens250 = sum(dens, na.rm = TRUE) * 0.05),
                      by = .(datetime, year)]

  get_total_hours_full_nights <- function(threshold, data) {
    data[, is_high := dens250 > threshold]
    bursts <- data[is_high == TRUE]
    total_hours <- bursts[, .(sum_hours = sum(hrs), threshold = threshold), by = .(year)]
    total_hours[, avg_hours := mean(sum_hours)]
    return(total_hours)
  }

  get_total_hours <- function(threshold, data) {
    data[, is_high := dens250 > threshold]
    highs <- data[is_high == TRUE]
    total_hours <- highs[, .N, by = .(year)][, .(year, sum_hours = N, threshold = threshold)]
    total_hours[, avg_hours := mean(sum_hours)]
    return(total_hours)
  }

  thresholds_night <- quantile(dt_agg_night$dens250, probs = seq(from = 0.5, to = 1, by = 0.0001), na.rm = TRUE)
  results_night <- lapply(thresholds_night, get_total_hours_full_nights, data = dt_agg_night) %>%
    bind_rows()

  thresholds_hourly <- quantile(dt_agg_hourly$dens250, probs = seq(from = 0.5, to = 1, by = 0.0001), na.rm = TRUE)
  results_hourly <- lapply(thresholds_hourly, get_total_hours, data = dt_agg_hourly) %>%
    bind_rows()

  # ggplot(results_night) +
  #   geom_point(aes(x = threshold, y = sum_hours, color = year)) +
  #   geom_line(aes(x = threshold, y = avg_hours), color = "black") +
  #   geom_smooth(aes(x = threshold, y = avg_hours)) +
  #   geom_hline(aes(yintercept = 100)) -> plot_threshold_hours_night
  #
  # ggplot(results_hourly) +
  #   geom_point(aes(x = threshold, y = sum_hours, color = year)) +
  #   geom_line(aes(x = threshold, y = avg_hours), color = "black") +
  #   geom_smooth(aes(x = threshold, y = avg_hours)) +
  #   geom_hline(aes(yintercept = 100)) -> plot_threshold_hours_hourly

  best_threshold_nightly <- results_night[which.min(abs(results_night$avg_hours - 100)), "threshold"]
  print(best_threshold_nightly)

  best_threshold_hourly <- results_hourly[which.min(abs(results_hourly$avg_hours - 100)), "threshold"]
  print(best_threshold_hourly)

  preds %>%
    select(year, height, truth, response, datetime, night_num) %>%
    {if (!is.null(height_cutoff)) filter(., height <= height_cutoff) else .} -> altitude_data

  altitude_data %>%
    group_by(year) %>%
    arrange(desc(truth)) %>%
    slice(100) %>%
    ungroup() %>%
    summarize(threshold = mean(truth)) -> altitude_threshold

  print(altitude_threshold)

  preds %>%
    select(year, height, truth, response, datetime, night_num) %>%
    filter(height <= !!height_cutoff) %>%
    mutate(year = as.factor(year)) %>%
    group_by(year, night_num, datetime) %>%
    summarise(
      vid_pred = sum(response) * 0.05,  # Thickness of VP layer (50m)
      vid_truth = sum(truth) * 0.05
    ) %>%
    group_by(year, night_num) %>%
    summarise(
      sum_vid_pred = sum(vid_pred),
      sum_vid_truth = sum(vid_truth)
    ) %>%
    identity() -> nightly_data

  preds %>%
    select(year, height, truth, response, datetime, night_num) %>%
    filter(height <= !!height_cutoff) %>%
    mutate(year = as.factor(year)) %>%
    group_by(year, night_num, datetime) %>%
    summarise(
      vid_pred = sum(response) * 0.05,  # Thickness of VP layer (50m)
      vid_truth = sum(truth) * 0.05
    ) %>%
    group_by(year, datetime) %>%
    summarise(
      sum_vid_pred = sum(vid_pred),
      sum_vid_truth = sum(vid_truth)
    ) %>%
    identity() -> hourly_data

  ggplot() +
    geom_roc(aes(d = truth > as.vector(altitude_threshold), m = response, color = as.factor(year)), n.cuts = 0, data = altitude_data) +
    ggsci::scale_color_observable(name = "Year") +
    ggtitle(
      paste0("Altitude-specific forecast (50-", height_cutoff, "m)"),
      bquote("100-hr threshold value: " ~ .(round(as.numeric(altitude_threshold), 0)) ~ "[birds" ~ km^{-3}*"]")
    ) +
    coord_equal() -> plot_altitude

  altitude_roc <- layer_data(plot_altitude) %>%
    select(true_positive_fraction, false_positive_fraction, label, group) %>%
    mutate(type = "altitude")

  ggplot() +
    geom_roc(aes(d = sum_vid_truth > !!as.numeric(best_threshold_nightly), m = sum_vid_pred, color = year), n.cuts = 0, data = nightly_data) +
    ggsci::scale_color_observable(name = "Year") +
    # labs(
    #   title = paste0("Optimal 100-hr threshold (vertical integration 50-", height_cutoff, "m)"),
    #   subtitle = paste0("Threshold value of nightly summed VID: ", round(best_threshold_nightly, 0), " [birds/km2]")
    # ) +
    ggtitle(
      paste0("Vertically integrated forecast (50-", height_cutoff, "m)"),
      bquote("100-hr threshold value: " ~ .(round(as.numeric(best_threshold_nightly), 0)) ~ "[birds" ~ km^{-2} ~ night^{-1}*"]")
    ) +
    coord_equal() -> plot_nightly

  nightly_roc <- layer_data(plot_nightly) %>%
    select(true_positive_fraction, false_positive_fraction, label, group) %>%
    mutate(type = "nightly")

  ggplot() +
    geom_roc(aes(d = sum_vid_truth > !!as.numeric(best_threshold_hourly), m = sum_vid_pred, color = year), n.cuts = 0, data = hourly_data) +
    ggsci::scale_color_observable(name = "Year") +
    ggtitle(
      paste0("Vertically integrated forecast (50-", height_cutoff, "m)"),
      bquote("100-hr threshold value: " ~ .(round(as.numeric(best_threshold_hourly), 0)) ~ "[birds" ~ km^{-2} ~ hour^{-1}*"]")
    ) +
    coord_equal() -> plot_hourly

  hourly_roc <- layer_data(plot_hourly) %>%
    select(true_positive_fraction, false_positive_fraction, label, group) %>%
    mutate(type = "hourly")

  bind_rows(altitude_roc, nightly_roc, hourly_roc) %>%
    mutate(year = as.factor(group + min(years) - 1))
}

interpolate_roc <- function(roc_curves) {
  common_fpf <- seq(0, 1, by = 1e-3)
  roc_curves %>%
    group_by(type, year) %>%
    group_split() %>%
    map_dfr(function(df) {
      df <- df[order(df$false_positive_fraction), ] # Order just in case
      tibble(
        false_positive_fraction = common_fpf,
        true_positive_fraction = approx(
          x = df$false_positive_fraction,
          y = df$true_positive_fraction,
          xout = common_fpf,
          method = "constant",   # nearest-neighbor interpolation
          f = 0,                 # "left" constant
          rule = 2               # extend to ends
        )$y,
        type = df$type[1],
        year = df$year[1]
      )
    }) %>%
    mutate(
      true_positive_fraction = if_else(false_positive_fraction == 0, 0, true_positive_fraction),
      true_positive_fraction = if_else(false_positive_fraction == 1, 1, true_positive_fraction)
    )
}

mean_roc <- function(interpolated_roc_curves) {
  interpolated_roc_curves %>%
    group_by(type, false_positive_fraction) %>%
    summarize(mean_true_positive_fraction = mean(true_positive_fraction, na.rm = TRUE), .groups = "drop")
}

comp_auc <- function(df, yearly = FALSE) {
  # Modified from plotROC function
  auc <- 0
  for (i in 2:nrow(df)) {
    if (yearly) {
      auc <- auc + 0.5 * (df$false_positive_fraction[i] - df$false_positive_fraction[i - 1]) *
        (df$true_positive_fraction[i] + df$true_positive_fraction[i - 1])
    } else {
      auc <- auc + 0.5 * (df$false_positive_fraction[i] - df$false_positive_fraction[i - 1]) *
        (df$mean_true_positive_fraction[i] + df$mean_true_positive_fraction[i - 1])
    }
  }
  if (yearly) {
    return(tibble(type = df$type[1], year = df$year[1], AUC = auc, season = df$season[1], altitudes = df$altitudes[1]))
  } else {
    return(tibble(type = df$type[1], AUC = auc, season = df$season[1], altitudes = df$altitudes[1]))
  }
}

roc_curves_fullprofile_autumn <- plot_roc_curves(preds_autumn, height_cutoff = 3000)
interpolated_roc_autumn_fullprofile <- interpolate_roc(roc_curves_fullprofile_autumn) %>% mutate(season = "Autumn")
mean_roc_autumn_fullprofile <- mean_roc(interpolated_roc_autumn_fullprofile) %>% mutate(season = "Autumn")

roc_curves_250m_autumn <- plot_roc_curves(preds_autumn, height_cutoff = 250)
interpolated_roc_autumn_250m <- interpolate_roc(roc_curves_250m_autumn) %>% mutate(season = "Autumn")
mean_roc_autumn_250m <- mean_roc(interpolated_roc_autumn_250m) %>% mutate(season = "Autumn")

roc_curves_fullprofile_spring <- plot_roc_curves(preds_spring, height_cutoff = 3000)
interpolated_roc_spring_fullprofile <- interpolate_roc(roc_curves_fullprofile_spring) %>% mutate(season = "Spring")
mean_roc_spring_fullprofile <- mean_roc(interpolated_roc_spring_fullprofile) %>% mutate(season = "Spring")

roc_curves_250m_spring <- plot_roc_curves(preds_spring, height_cutoff = 250)
interpolated_roc_spring_250m <- interpolate_roc(roc_curves_250m_spring) %>% mutate(season = "Spring")
mean_roc_spring_250m <- mean_roc(interpolated_roc_spring_250m) %>% mutate(season = "Spring")

observable_colors <- pal_observable(palette = c("observable10"), alpha = 1)(10)
scales::show_col(observable_colors[c(1,2,5)])
observable_colors <- observable_colors[c(1, 2, 5)]

interpolated_roc_fullprofile <- bind_rows(interpolated_roc_autumn_fullprofile, interpolated_roc_spring_fullprofile) %>% mutate(altitudes = "50-3000 m")
mean_roc_fullprofile <- bind_rows(mean_roc_autumn_fullprofile, mean_roc_spring_fullprofile) %>% mutate(altitudes = "50-3000 m")

interpolated_roc_250m <- bind_rows(interpolated_roc_autumn_250m, interpolated_roc_spring_250m) %>% mutate(altitudes = "<300 m")
mean_roc_250m <- bind_rows(mean_roc_autumn_250m, mean_roc_spring_250m) %>% mutate(altitudes = "<300 m")

interpolated_roc <- bind_rows(interpolated_roc_fullprofile, interpolated_roc_250m) %>%
  mutate(altitudes = factor(altitudes, levels = c("50-3000 m", "<300 m")))
mean_rocs <- bind_rows(mean_roc_fullprofile, mean_roc_250m) %>%
  mutate(altitudes = factor(altitudes, levels = c("50-3000 m", "<300 m")))

mean_rocs %>%
  group_by(type, season, altitudes) %>%
  group_split() %>%
  map_df(comp_auc) -> auc_type

interpolated_roc %>%
  group_by(type, year, season, altitudes) %>%
  group_split() %>%
  map_df(comp_auc, yearly = TRUE) %>%
  group_by(type, season, altitudes) %>%
  mutate(min_AUC = min(AUC),
         max_AUC = max(AUC)) %>%
  select(type, min_AUC, max_AUC, season, altitudes) %>%
  distinct() -> auc_type_minmax

left_join(auc_type, auc_type_minmax) %>%
  mutate(label_long = sprintf("%s %.3f (%.3f-%.3f)", season, round(AUC, digits = 3), min_AUC, max_AUC),
         label_short = sprintf("%s: %.3f", season, round(AUC, digits = 3)),
         label_long = str_replace_all(label_long, c("\\(0\\." = "\\(\\.", "-0\\." = "-\\.", "0\\." = "\\.")),
         label_short = str_replace(label_short, "\\(0\\.", "\\.")) %>%
  select(type, season, label_long, label_short, altitudes) %>%
  mutate(x = 1, y = c(0.34, 0.34, 0.28, 0.28, 0.20, 0.20, 0.14, 0.14, 0.06, 0.06, 0, 0)) -> auc_scores

roc_major_breaks <- c(0, .1, .25, .5, .75, .9, 1)
roc_minor_breaks <- c(seq(0, .1, by = .01), seq(.9, 1, by = .01))
roc_x_minor_breaks <- c(seq(0, .1, by = .01))
roc_y_minor_breaks <- c(seq(.9, 1, by = .01))

ggplot() +
  geom_path(data = interpolated_roc,
            aes(x = false_positive_fraction,
                y = true_positive_fraction,
                group = interaction(type, year, season),
                color = type,
                linetype = season),
            alpha = 0.4) +
  geom_path(data = mean_rocs,
            aes(x = false_positive_fraction,
                y = mean_true_positive_fraction,
                group = interaction(type, season),
                color = type,
                linetype = season),
            linewidth = 1.2) +
  theme_bw(base_size = 14) +
  facet_wrap(vars(altitudes)) +
  scale_x_continuous(breaks = roc_major_breaks, minor_breaks = roc_x_minor_breaks) +
  scale_y_continuous(breaks = roc_major_breaks, minor_breaks = roc_y_minor_breaks) +
  scale_linetype(name = "Season") +
  labs(x = "False positive fraction", y = "True positive fraction") +
  coord_fixed() +
  scale_color_manual(values = c(observable_colors),
                     labels = c("Altitude-specific", "Hourly VID", "Nightly VID"), name = "Threshold") -> roc_baseplot

roc_baseplot +
  geom_text(aes(x = x, y = y, label = label_long, color = type), data = auc_scores, hjust = 1, vjust = 0, size = 6, show.legend = FALSE) +
  annotate(geom = "text", x = 1, y = 0.41, hjust = 1, vjust = 0, label = "Average AUC (min-max)", size = 6, color = "gray10")

ggsave(murmuR_path("models/plots/roc/roc_curves_250_3000_longlabels.pdf"), width = 12, height = 6)

roc_baseplot +
  geom_text(aes(x = x, y = y, label = label_short, color = type), data = auc_scores, hjust = 1, vjust = 0, size = 6, show.legend = FALSE) +
  annotate(geom = "text", x = 1, y = 0.41, hjust = 1, vjust = 0, label = "Average AUC", size = 6, color = "gray10")

ggsave(murmuR_path("models/plots/roc/roc_curves_250_3000_shortlabels.pdf"), width = 12, height = 6)
