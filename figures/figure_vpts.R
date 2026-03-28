library(tidyverse)
library(ggpp)
library(patchwork)
library(ggrastr)
library(ggh4x)
library(ggdark)

library(murmuR)
options(murmuR.data_path = "/data/birdcloudstorage-tvm/ibm-ml-refactor/data/")

annotate_peaknights <- function(preds, height_cutoff = FALSE, top_n = 10) {
  preds %>%
    { if (height_cutoff) filter(., height <= 250) else . } %>%
    group_by(datetime) %>%
    mutate(
      vid_pred = sum(response) * 0.05,
      vid_truth = sum(truth) * 0.05
    ) %>%
    ungroup() %>%
    group_by(night_num) %>%
    summarise(
      sum_vid_pred = sum(vid_pred),
      sum_vid_truth = sum(vid_truth),
      x_min = min(datetime),
      x_max = max(datetime)
    ) %>%
    mutate(
      vid_pred_rank = rank(-sum_vid_pred, ties.method = "min"),
      vid_truth_rank = rank(-sum_vid_truth, ties.method = "min"),
      peak_pred = as.numeric(vid_pred_rank <= top_n),
      peak_truth = as.numeric(vid_truth_rank <= top_n),
      peak_both = as.integer(peak_pred & peak_truth),
      peak_missed = as.integer(!peak_pred & peak_truth)
    ) -> peaknights
  preds %>%
    left_join(peaknights)
}

for (radar in names(radars_named)) {
for (season in names(months_seasons)) {
for (year in years) {

  preds_file <- sprintf(murmuR_path("models/predictions/regr.rmse_dynamic/%s/%s_%s_%s_dynamic_preds.RDS"), radar, radar, season, year)

  preds_dynamic <- readRDS(preds_file) %>%
    mutate(pred_type = "dynamic") %>%
    # rename(pred_dynamic = response) %>%
    annotate_peaknights() %>%
    identity()

  # preds_dynamic %>%
  #   mutate(datetime = as.factor(datetime),
  #          nightchange = c(0, diff(night_num)),
  #          peak_label = peak_pred + peak_both,
  #          peak_label = factor(peak_label, levels = c(0, 1, 2), labels = c("None", "Peak", "Both"))) %>%
  #   identity() -> plot_data_dynamic

  preds_dynamic %>%
    mutate(datetime = as.factor(datetime),
           nightchange = c(0, diff(night_num)),
           peak_label = case_when(
             peak_pred == 0 & peak_truth == 0 ~ "TN",
             peak_pred == 1 & peak_truth == 0 ~ "FP",
             peak_pred == 0 & peak_truth == 1 ~ "FN",
             peak_pred == 1 & peak_truth == 1 ~ "TP"
           ),
           peak_label = factor(peak_label, levels = c("TN", "FP", "FN", "TP"), labels = c("True negative", "False positive", "False negative", "True positive"))) %>%
    identity() -> plot_data_dynamic

  preds_dynamic %>%
    mutate(datetime = as.factor(datetime),
           nightchange = c(0, diff(night_num)),
           # peak_label = case_when(
           #   peak_pred == 0 & peak_truth == 0 ~ "TN",
           #   peak_pred == 1 & peak_truth == 0 ~ "FP",
           #   peak_pred == 0 & peak_truth == 1 ~ "FN",
           #   peak_pred == 1 & peak_truth == 1 ~ "TP"
           # ),
           peak_label = case_when(
             peak_truth == 0 ~ "TN",
             peak_truth == 1 ~ "TP"
           ),
           peak_label = factor(peak_label, levels = c("TN", "FP", "FN", "TP"), labels = c("True negative", "False positive", "False negative", "True positive"))) %>%
    identity() -> plot_data_truth

  nightly_breaks_dynamic <- plot_data_dynamic %>%
    filter(nightchange > 0) %>%
    select(datetime) %>%
    pull(datetime)

  nightly_breaks_diff <- sample(nightly_breaks_dynamic, 10)

  nightly_breaks_truth <- plot_data_truth %>%
    filter(nightchange > 0) %>%
    select(datetime) %>%
    pull(datetime)

  metrics_description_dynamic <- sprintf("Peak nights: %d/10\nPeak hours: %d/100\nProp. of migration: %.0f%%",
                                         unique(plot_data_dynamic$regr.peaknights),
                                         unique(plot_data_dynamic$regr.peakhours),
                                         unique(plot_data_dynamic$regr.propmigrtt) * 100)

  # Color scale from colorbrewer: https://observablehq.com/@kanitw/categorical-color-plots
  ggplot(plot_data_dynamic) +
    geom_tile(aes(x = datetime, y = height, fill = response)) +
    geom_point(aes(x = datetime, color = peak_label), y = -100, shape = "|", size = 6, position = position_nudge(x = 2)) + # Nudge position to align with ticks
    # geom_point(aes(x = datetime, y = rdr_blh), color = "white", size = 0.25) + # BLH
    scale_fill_viridis_c(trans = "sqrt", option = "magma", limits = c(0, 300), na.value = 0, expression("Bird density [" * birds ~ km^{-3} * "]"),
                         oob = scales::squish, guide = "none") +
    scale_x_discrete(breaks = nightly_breaks_dynamic,
                     labels = function(x) {
                       dates <- as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S")
                       condition <- day(dates) %in% c(1, 15) | (day(dates) %in% c(30) & month(dates) %in% c(11))
                       newdates <- format(dates, "%b %-d")
                       newdates[!condition] <- ""
                       newdates[newdates == "NA"] <- ""
                       newdates
                     }) +
    annotate(geom = "text_npc", npcx = 0.01, npcy = 0.95, label = "Predicted", color = "white", size = 8, hjust = 0, vjust = 1, fontface = "bold") +
    annotate(geom = "text_npc", npcx = 0.01, npcy = 0.80, label = metrics_description_dynamic, color = "white", size = 5, hjust = 0, vjust = 1) +
    scale_color_manual(values = c("True negative" = NA, "False negative" = "#ff7f00", "True positive" = "#33a02c", "False positive" = "#e31a1c"),
                       labels = c("True negative" = "", "False negative" = "False negative", "True positive" = "True positive", "False positive" = "False positive"),
                       breaks = c(NA, "False negative", "True positive", "False positive"),
                       name = "Peak night", na.value = "#00000000") +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "black"),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white"),
          axis.ticks = element_line(colour = "white", linewidth = 0.5),
      panel.background = element_rect(fill = "black"),
      panel.grid = element_blank(),
      legend.position.inside = c(.98, .98),
      legend.justification = c("right", "top"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(color = "white", vjust = 0.5),
      legend.title = element_text(color = "white", vjust = 0.5),
      legend.key.width = unit(2, "lines"),
      axis.title.x=element_blank(),
      #axis.text.x=element_blank(),
      # axis.ticks = element_line(colour = "black", linewidth = 0.5)
    ) +
    guides(#fill = guide_colorbar(position = "inside", direction = "horizontal"),
           color = guide_legend(position = "inside", override.aes = list(size = 6, shape = 15), direction = "horizontal")) +
    coord_cartesian(ylim = c(-75, 3000), expand = FALSE) +
    labs(x = "",
         y = "Height [m]") -> plot_dynamic

  ggplot(plot_data_truth) +
    geom_tile(aes(x = datetime, y = height, fill = truth)) +
    geom_point(aes(x = datetime, color = peak_label), y = -100, shape = "|", size = 6, position = position_nudge(x = 2)) + # Nudge position to align with ticks
    # geom_point(aes(x = datetime, y = rdr_blh), color = "white", size = 0.25) + # BLH
    scale_fill_viridis_c(trans = "sqrt", option = "magma", limits = c(0, 300), name = expression("Bird density [" * birds ~ km^{-3} * "]"), guide = "none",
                         oob = scales::squish) +
    scale_x_discrete(breaks = nightly_breaks_dynamic,
                     labels = function(x) {
                       dates <- as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S")
                       condition <- day(dates) %in% c(1, 15) | (day(dates) %in% c(30) & month(dates) %in% c(11))
                       newdates <- format(dates, "%b %-d")
                       newdates[!condition] <- ""
                       newdates[newdates == "NA"] <- ""
                       newdates
                     }) +
    annotate(geom = "text_npc", npcx = 0.01, npcy = 0.95, label = "Measured", color = "white", size = 8, hjust = 0, vjust = 1, fontface = "bold") +
    scale_color_manual(values = c("True negative" = NA, "False negative" = "#ff7f00", "True positive" = "#33a02c", "False positive" = "#e31a1c"),
                       labels = c("True negative" = "", "False negative" = "False negative", "True positive" = "True positive", "False positive" = "False positive"),
                       breaks = c(NA, "False negative", "True positive", "False positive"),
                       name = "Peak night", guide = "none", na.value = "#00000000") +
    theme_classic(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "black"),
          axis.title = element_text(color = "white"),
          axis.text = element_text(color = "white"),
          axis.ticks = element_line(colour = "white", linewidth = 0.5),
          panel.background = element_rect(fill = "black"),
          panel.grid = element_blank(),
          legend.position.inside = c(.98, .98),
          legend.justification = c("right", "top"),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white", vjust = 1),
          legend.key.width = unit(2, "lines"),
          axis.title.x=element_blank(),
          # axis.text.x=element_blank(),
          #axis.ticks = element_line(colour = "black", linewidth = 0.5)
    ) +
    guides(fill = guide_colorbar(position = "inside", direction = "horizontal")) +
    coord_cartesian(ylim = c(-75, 3000), expand = FALSE) +
    labs(x = "",
         y = "Height [m]") -> plot_truth

  plot_vpts <-
    plot_truth / plot_dynamic  +
    plot_layout(axis_titles = "collect") & theme(plot.background = element_rect(fill = "black", colour = "black"))

  file_out <- sprintf(murmuR_path("models/plots/vpts_dark/vpts_%s_%s_%s.pdf"), radar, season, year)
  file_out_png <- str_replace(file_out, ".pdf", ".png")

  ggsave(filename = file_out, plot = plot_vpts, width = 12, height = 8)
  ggsave(filename = file_out_png, plot = plot_vpts, width = 12, height = 8)

  print(paste0("Finished writing: ", file_out))
}
}
}


preds_dynamic %>%
  select(vid_pred_rank, vid_truth_rank, night_num) %>%
  distinct(night_num, vid_pred_rank, vid_truth_rank) %>%
  # left_join(preds_dynamic %>% select(night_num, datetime) %>% mutate(date = lubridate::date(datetime))) %>%
  # distinct(night_num, vid_pred_rank, vid_truth_rank, date) %>%
  arrange(night_num) %>%
  mutate(date = seq(from = as.Date("2018-08-01"), to = as.Date("2018-12-01"), by = "1 day")) %>%
  identity() -> night_ranks

selected_nights <- as.Date(c("2018-09-25", "2018-09-28", "2018-10-06", "2018-10-10", "2018-10-11", "2018-10-13", "2018-10-14", "2018-10-15",
                     "2018-10-16", "2018-10-17", "2018-10-18", "2018-10-19", "2018-10-27", "2018-10-28", "2018-11-03", "2018-11-04",
                     "2018-11-05", "2018-11-06", "2018-11-08", "2018-11-11", "2018-11-15", "2018-11-16", "2018-11-17"))

data.frame(date = selected_nights, status = "peak") -> selected_peaks

night_ranks %>%
  pivot_longer(cols = c(vid_truth_rank, vid_pred_rank), values_to = "rank", names_to = "ranktype") %>%
  left_join(selected_peaks) %>%
  mutate(status = replace_na(status, "no peak"),
         status = case_when(
           status == "no peak" ~ "No peak",
           status == "peak" ~ "Peak"
         ),
         ranktype = case_when(
           ranktype == "vid_pred_rank" ~ "Predicted rank",
           ranktype == "vid_truth_rank" ~ "Measured rank",
         )) %>%
  # identity() -> a
  ggplot(aes(x = date, y = rank, color = ranktype)) +
  geom_point(aes(size = status)) +
  geom_line() +
  scale_y_reverse() +
  scale_color_discrete(name = "Herwijnen night rank") +
  scale_size_discrete(name = "Eemshaven peak night", range = c(1, 3)) +
  labs(x = "Date", y = "Ranked night (1 = most migration; 123 = least migration)") +
  theme_bw(base_size = 10)

ggsave("ranking.pdf", width = 14, height = 5)

cor(night_ranks$vid_pred_rank, night_ranks$vid_truth_rank)

peak_nights_top_20 <- night_ranks %>%
  filter(vid_truth_rank <= 20) %>%
  mutate(diff = vid_truth_rank - vid_pred_rank) %>%
  summarise(m_diff = mean(diff),
            min_diff = min(diff),
            max_diff = max(diff))

cor(peak_nights_top_20$vid_truth_rank, peak_nights_top_20$vid_pred_rank)