library(tidyverse)
library(data.table)
library(ggpp)
library(ggsci)
library(ggrepel)

library(murmuR)
options(murmuR.data_path = "/data/birdcloudstorage-tvm/ibm-ml-refactor/data/")

for (sampling in c("dynamic", "radar")) {
for (radar in names(radars_named)) {
  for (season in names(months_seasons)) {
    print(sprintf("Parsing %s %s %s", radar, season, sampling))
    preds_files <- sprintf(murmuR_path("models/predictions/regr.rmse_%s/%s/*_%s_*_preds.RDS"), sampling, radar, season)
    preds <- Sys.glob(preds_files) %>%
      mclapply(readRDS) %>%
      bind_rows()

    preds %>%
      group_by(year, datetime, mtr) %>%
      summarise(vid_truth = sum(truth) * 0.05) %>%
      ungroup() %>%
      arrange(year, desc(vid_truth)) %>%
      group_by(year) %>%
      mutate(cumsum_mtr_truth = cumsum(mtr),
             cumsum_mtr_truth = (cumsum_mtr_truth / max(cumsum_mtr_truth)) * 100,
             ranked_hour = row_number(),
             year = as.factor(year)) -> mtr_truth

    preds %>%
      group_by(year, datetime, mtr) %>%
      summarise(vid_pred = sum(response) * 0.05) %>%
      ungroup() %>%
      arrange(year, desc(vid_pred)) %>%
      group_by(year) %>%
      mutate(cumsum_mtr_pred = cumsum(mtr),
             cumsum_mtr_pred = (cumsum_mtr_pred / max(cumsum_mtr_pred)) * 100,
             ranked_hour = row_number(),
             year = as.factor(year)) -> mtr_pred

    mtr_pred %>%
      mutate(cumsum50 = abs(50 - cumsum_mtr_pred)) %>%
      filter(cumsum50 == min(cumsum50)) %>%
      group_by(year) %>%
      slice_head(n = 1) %>%
      mutate(label = paste0(ranked_hour, "hr")) %>%
      select(year, cumsum_mtr_pred, ranked_hour, label) -> cumsum50_hrs_pred

    mtr_pred %>%
      filter(ranked_hour == 100) %>%
      mutate(label = paste0(round(cumsum_mtr_pred, digits = 0), "%")) %>%
      select(year, cumsum_mtr_pred, ranked_hour, label) -> cumsum_100hrs_pred

    mtr_truth %>%
      mutate(cumsum50 = abs(50 - cumsum_mtr_truth)) %>%
      filter(cumsum50 == min(cumsum50)) %>%
      group_by(year) %>%
      slice_head(n = 1) %>%
      mutate(label = paste0(ranked_hour, "hr")) %>%
      select(year, cumsum_mtr_truth, ranked_hour, label) -> cumsum50_hrs_truth

    mtr_truth %>%
      filter(ranked_hour == 100) %>%
      mutate(label = paste0(round(cumsum_mtr_truth, digits = 0), "%")) %>%
      select(year, cumsum_mtr_truth, ranked_hour, label) -> cumsum_100hrs_truth

    mtr_truth %>%
      left_join(mtr_pred, by = join_by(year, ranked_hour)) %>%
      mutate(diff = cumsum_mtr_truth - cumsum_mtr_pred) %>%
      identity() -> mtr_differences

    bind_rows(
      cumsum_100hrs_pred %>% rename(cumsum_mtr = cumsum_mtr_pred) %>% mutate(type = "pred", color = year, class = "100hr"),
      cumsum_100hrs_truth %>% rename(cumsum_mtr = cumsum_mtr_truth) %>% mutate(type = "truth", color = "black", class = "100hr"),
      cumsum50_hrs_pred %>% rename(cumsum_mtr = cumsum_mtr_pred) %>% mutate(type = "pred", color = year, class = "50%"),
      cumsum50_hrs_truth %>% rename(cumsum_mtr = cumsum_mtr_truth) %>% mutate(type = "truth", color = "black", class = "50%")
    ) -> label_points

    label_points <- label_points %>%
      mutate(radar = !!radar,
             season = !!season,
             sampling = !!sampling)

    label_points_out <- sprintf(murmuR_path("models/cum_mtr_ranked_hours/label_points_%s_%s_%s.RDS"), radar, season, sampling)
    saveRDS(label_points, file = label_points_out)

    label_points %>%
      filter(class == "100hr") %>%
      group_by(type) %>%
      summarise(mean_100hrs = mean(cumsum_mtr),
                min_100hrs = min(cumsum_mtr),
                max_100hrs = max(cumsum_mtr))

    label_points %>%
      filter(class == "50%") %>%
      group_by(type) %>%
      summarise(mean_50perc = mean(ranked_hour),
                min_50perc = min(ranked_hour),
                max_50perc = max(ranked_hour))

    observable_colors <- pal_observable(palette = c("observable10"), alpha = 1)(8)
    observable_colors <- c(observable_colors, "#000000")

    ggplot() +
      geom_col(aes(x = ranked_hour, y = diff, group = year), data = mtr_differences, color = "grey") +
      # Model
      geom_path(aes(x = ranked_hour, y = cumsum_mtr_pred, group = year, color = year), data = mtr_pred, linewidth = 1) +
      # geom_segment(aes(y = cumsum_mtr_pred, yend = -Inf, x = ranked_hour, xend = ranked_hour, color = year), data = cumsum50_hrs_pred, linetype = 3) +
      geom_segment(aes(y = cumsum_mtr_pred, yend = cumsum_mtr_pred, x = ranked_hour, xend = -Inf, color = year), data = cumsum50_hrs_pred, linetype = 4) +
      geom_segment(aes(y = cumsum_mtr_pred, yend = -Inf, x = ranked_hour, xend = ranked_hour, color = year), data = cumsum_100hrs_pred, linetype = 2) +
      # geom_segment(aes(y = cumsum_mtr_pred, yend = cumsum_mtr_pred, x = ranked_hour, xend = -Inf, color = year), data = cumsum_100hrs_pred, linetype = 2) +
      geom_point(aes(x = ranked_hour, y = cumsum_mtr_pred, color = year), data = cumsum50_hrs_pred) +
      # Measurements
      geom_path(aes(x = ranked_hour, y = cumsum_mtr_truth, group = year), data = mtr_truth, color = "black", linewidth = 1) +
      # geom_segment(aes(y = cumsum_mtr_truth, yend = -Inf, x = ranked_hour, xend = ranked_hour), data = cumsum50_hrs_truth, linetype = 3) +
      geom_segment(aes(y = cumsum_mtr_truth, yend = cumsum_mtr_truth, x = ranked_hour, xend = -Inf), data = cumsum50_hrs_truth, linetype = 4) +
      geom_segment(aes(y = cumsum_mtr_truth, yend = -Inf, x = ranked_hour, xend = ranked_hour), data = cumsum_100hrs_truth, linetype = 2) +
      # geom_segment(aes(y = cumsum_mtr_truth, yend = cumsum_mtr_truth, x = ranked_hour, xend = -Inf), data = cumsum_100hrs_truth, linetype = 2) +
      geom_point(aes(x = ranked_hour, y = cumsum_mtr_truth), data = cumsum50_hrs_truth) +
      # Labels
      ggrepel::geom_label_repel(aes(x = ranked_hour, y = cumsum_mtr, label = label, color = color), data = label_points, fontface = "bold",
                                min.segment.length = 0, force = 2, force_pull = 0.5,
                                nudge_x = if_else(label_points$type == "truth", -50, 50),
                                nudge_y = if_else(label_points$type == "truth", 5, -5)) +
      # Plots
      scale_color_manual(values = observable_colors, guide = "none") +
      facet_wrap(vars(year), nrow = 2) +
      # coord_cartesian(expand = FALSE) +
      coord_cartesian(expand = FALSE, xlim = c(0, 500), ylim = c(0, 87.5)) +
      # scale_x_continuous(breaks = c(0, 100, 500, 1000, 1500), labels = c(0, 100, 500, 1000, "")) +
      scale_x_continuous(breaks = c(0, 100, 200, 300, 400, 500), labels = c(0, 100, 200, 300, 400, "")) +
      scale_y_continuous(breaks = c(0, 25, 50, 75), labels = c(0, 25, 50, 75)) +
      labs(x = "Ranked hours",
           y = "Cumulative migration traffic [%]") +
      theme_bw(base_size = 14) -> main


    # Insets, following https://r4photobiology.info/galleries/plot-insets.html
    make_inset <- function(year) {
      color <- observable_colors[which(year == years)]
      ggplot() +
        geom_path(aes(x = ranked_hour, y = cumsum_mtr_pred, group = year, color = year), data = mtr_pred %>% filter(year == !!year), linewidth = 1) +
        geom_path(aes(x = ranked_hour, y = cumsum_mtr_truth, group = year), data = mtr_truth %>% filter(year == !!year), color = "black", linewidth = 1) +
        scale_color_manual(values = color, guide = "none") +
        scale_x_continuous(breaks = c(0, 500, 1000, 1500), labels = c(0, "", 1000, "")) +
        scale_y_continuous(breaks = c(0, 50, 100)) +
        theme_classic() +
        theme(panel.background = element_blank(),
              plot.background = element_blank(),
              plot.margin = unit(c(0,0,0,0), "cm")) +
        labs(x = "",
             y = "")
    }

    # list of box plots as insets
    p.inset.list <- lapply(unique(mtr_pred$year), make_inset)

    # data frame to be used as data for plotting insets
    inset.df <- data.frame(year = unique(mtr_pred$year),
                           plots =  I(p.inset.list),
                           npcx = 0.95,
                           npcy = 0.01)

    main +
      geom_plot_npc(data = inset.df, mapping = aes(npcx = npcx, npcy = npcy, label = plots), vp.width = 0.45, vp.height = 0.4, vjust = 0)

    ggsave(filename = sprintf(murmuR_path("models/plots/cumulative_migration_ranked_hours_%s_%s_%s.pdf"), radar, season, sampling), width = 12, height = 8)

  }
}
}

Sys.glob(murmuR_path("models/cum_mtr_ranked_hours/*")) %>%
  lapply(readRDS) %>%
  bind_rows() -> hours_cum_mtr

hours_cum_mtr %>%
  filter(class == "100hr") %>%
  group_by(radar, season, sampling, type) %>%
  summarize(mean_100hr = mean(cumsum_mtr),
            min_100hr = min(cumsum_mtr),
            max_100hr = max(cumsum_mtr)) %>%
  mutate(mt_in_100hr = sprintf("%.f (%.f-%.f)", mean_100hr, min_100hr, max_100hr)) %>%
  select(radar, season, type, sampling, mt_in_100hr) %>%
  pivot_wider(names_from = type, values_from = mt_in_100hr) %>%
  select(radar, season, sampling, truth_100hr = truth, pred_100hr = pred) %>%
  arrange(sampling) %>%
  identity() -> mt_in_100hrs

# hours_cum_mtr %>%
#   filter(class == "100hr") %>%
#   group_by(radar, season, type, sampling) %>%
#   summarize(mean_100hr = mean(cumsum_mtr),
#             min_100hr = min(cumsum_mtr),
#             max_100hr = max(cumsum_mtr)) %>%
#   # mutate(mt_in_100_hr = sprintf("%.f (%.f-%.f)", mean_100hr, min_100hr, max_100hr)) %>%
#   select(radar, season, type, sampling, mean_100hr, min_100hr, max_100hr) %>%
#   pivot_wider(
#     id_cols = c(radar, season, type),
#     names_from = sampling,
#     values_from = c(mean_100hr, min_100hr, max_100hr),
#     names_glue = "{.value}_{sampling}"
#   ) %>%
#   mutate(
#     mean_improvement = mean_100hr_dynamic - mean_100hr_radar,
#   ) %>%
#   group_by(radar, season) %>%
#   mutate(
#     truth_mean = if_else(type == "truth", mean_100hr_dynamic, NA_real_),
#   ) %>%
#   fill(truth_mean, .direction = "downup") %>% # truth_min, truth_max, .direction = "downup") %>%
#   filter(type == "pred") %>%
#   # ungroup() -> a
#   identity() %>%
#   mutate(mt_100hrs_pred = sprintf("%.f (%.f-%.f)", mean_100hr_dynamic, min_100hr_dynamic, max_100hr_dynamic),
#          mt_100hrs_truth = sprintf("%.f (%.f-%.f)", mean_100hr, min_100hr_dynamic, max_100hr_dynamic)) %>%
#   identity() -> a
#   # select(radar, season, )
#   # select(radar, season, type, truth_100hr = truth, pred_100hr = pred) %>%
#   # arrange(sampling) %>%
#   # identity() -> mt_in_100hrs

hours_cum_mtr %>%
  filter(class == "50%") %>%
  group_by(radar, season, sampling, type) %>%
  summarize(mean_50perc = mean(ranked_hour),
            min_50perc = min(ranked_hour),
            max_50perc = max(ranked_hour)) %>%
  mutate(hrs_for_50perc = sprintf("%.f (%.f-%.f)", mean_50perc, min_50perc, max_50perc)) %>%
  select(radar, season, type, sampling, hrs_for_50perc) %>%
  pivot_wider(names_from = type, values_from = hrs_for_50perc) %>%
  select(radar, season, sampling, truth_50perc = truth, pred_50perc = pred) %>%
  arrange(sampling) %>%
  identity() -> hrs_for_50perc

left_join(mt_in_100hrs, hrs_for_50perc) %>%
  mutate(radar_human = case_when(
    radar == "nldhl" ~ "Den Helder",
    radar == "nlhrw" ~ "Herwijnen",
  ),
  season_human = case_when(
    season == "spring" ~ "Spring",
    season == "autumn" ~ "Autumn"
  ),
  sampling_human = case_when(
    sampling == "radar" ~ "Radar-level",
    sampling == "dynamic" ~ "Dynamic"
  )) %>%
  ungroup() %>%
  select(Radar = radar_human, Season = season_human, sampling_human,
         `100hr measurements` = truth_100hr, `100hr model` = pred_100hr,
         `50% measurements` = truth_50perc, `50% model` = pred_50perc) %>%
  # identity()
  tt() # %>%
  # print("latex")
