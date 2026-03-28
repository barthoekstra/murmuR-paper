library(murmuR)
options(murmuR.data_path = "/data/birdcloudstorage-tvm/ibm-ml-refactor/data/")

library(data.table)
library(tidyverse)
library(ggdist)
library(ggh4x)

df <- Sys.glob(murmuR_path("models/sampling_benchmarks/*.RDS")) %>%
  lapply(readRDS) %>%
  bind_rows() %>%
  mutate(
    radar_human = case_when(
      radar == "nlhrw" ~ "Herwijnen",
      radar == "nldhl" ~ "Den Helder"
    ),
    season_human = str_to_sentence(season),
    sampling_human = case_when(
      sampling == "radar" ~ "Radar only",
      sampling == "static_stopover" ~ "+ Fix. Stopover",
      sampling == "static_stopover_enroute" ~ "+ Fix. Stopover + Fix. En-route",
      sampling == "dynamic" ~ "+ Dyn. Stopover + Dyn. En-route"
    ),
    sampling_human_short = case_when(
      sampling == "radar" ~ "R",
      sampling == "static_stopover" ~ "+Fix. S",
      sampling == "static_stopover_enroute" ~ "+Fix. S+E",
      sampling == "dynamic" ~ "+Dyn. S+E"
    ),
    radar_human = factor(radar_human, levels = c("Den Helder", "Herwijnen")),
    season_human = factor(season_human, levels = c("Spring", "Autumn")),
    sampling_human = factor(sampling_human, levels = c(
      "+ Dyn. Stopover + Dyn. En-route",
      "+ Fix. Stopover + Fix. En-route",
      "+ Fix. Stopover",
      "Radar only"
    )),
    sampling_human_short = factor(sampling_human_short, levels = c(
      "+Dyn. S+E", "+Fix. S+E", "+Fix. S", "R"
    )),
  ) %>%
  pivot_longer(cols = c("regr.rmse", "regr.rsq",
                        "regr.peaknights",
                        "regr.propmig",
                        "regr.propmigrtt",
                        "regr.peakhours"),
               names_to = "metric",
               values_to = "score") %>%
  filter(metric %in% c(
    "regr.rmse", "regr.rsq",
    "regr.peaknights", "regr.peaknightslowalt",
    "regr.propmigrtt", "regr.propmiglowaltrtt",
    "regr.peakhours", "regr.peakhourslowalt")) %>%
  mutate(
    metric = case_when(
      metric == "regr.rmse" ~ "RMSE",
      metric == "regr.rsq" ~ "R2",
      metric == "regr.peaknights" ~ "Peak nights (of 10)",
      metric == "regr.peakhours" ~ "Peak hours (of 100)",
      metric == "regr.propmigrtt" ~ "% of migration",
    ),
    metric = fct_rev(metric)) %>%
  drop_na(metric)

# Define which metrics are higher-is-better
higher_is_better <- c("R2", "% of migration", "Peak nights (of 10)", "Peak hours (of 100)")

# Compute best-performing model per radar × season × metric
best_df <- df %>%
  group_by(radar_human, season_human, metric, sampling_human, sampling_human_short) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  group_by(radar_human, season_human, metric) %>%
  filter(
    (metric %in% higher_is_better & mean_score == max(mean_score, na.rm = TRUE)) |
      (!metric %in% higher_is_better & mean_score == min(mean_score, na.rm = TRUE))
  )

# Create combined x-axis label
df <- df %>%
  # mutate(x_label = interaction(sampling_human, season_human, sep = " - ")) %>%
  mutate(x_label = sampling_human_short)
best_df <- best_df %>%
  # mutate(x_label = interaction(sampling_human, season_human, sep = " - "))
  mutate(x_label = sampling_human_short)

# Plot
ggplot(df, aes(y = x_label, x = score, color = sampling_human)) +
  stat_pointinterval(point_interval = "mean_qi") +
  geom_point(data = best_df,
             aes(y = x_label, x = mean_score),
             color = "white", size = 1, inherit.aes = FALSE) +
  ggsci::scale_color_observable(name = "Environ. sampling") +
  facet_grid2(
    vars(radar_human, season_human),
    vars(metric),
    strip = strip_nested(bleed = TRUE),
    scales = "free",
    labeller = labeller(metric = as_labeller(function(x) {
      sapply(x, function(label) {
        switch(label,
               "R2" = "R^2",
               "RMSE" = "RMSE",
               "Peak nights (of 10)" = "\"Peak nights (of 10)\"",
               "Peak hours (of 100)" = "\"Peak hours (of 100)\"",
               "% of migration" = "\"% of migration\"",
               label)
      })
    }, default = label_parsed))
  ) +
  labs(x = "Score",
       y = "Environmental sampling") +
  theme_bw(base_size = 14) +
  theme(plot.title.position = "plot",
        legend.key.size = unit(1, "cm")) +
  guides(color = guide_legend(position = "top", direction = "horizontal", reverse = TRUE))

ggsave(murmuR_path("models/plots/sampling_benchmark/benchmark.pdf"), width = 12, height = 8, create.dir = TRUE)

df %>% filter(metric == "R2") %>% group_by(radar_human, sampling_human) %>% summarise(meanr2 = mean(score))
