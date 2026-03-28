# ============================================================
# Training analysis: init benchmarks → feature selection →
# sampling benchmarks → HPO results →
# predictions & importances → VPTS visualisations → ROC
# ============================================================

# --- Libraries -----------------------------------------------
library(murmuR)
library(tidyverse)
library(data.table)
library(patchwork)
library(tinytable)
library(plotROC)
library(ggh4x)
library(ggpp)
library(ggdist)
library(tidytext)
library(tictoc)
library(future.apply)
# library(mlr3)
# library(mlr3extralearners)

# --- Options -------------------------------------------------
# options(murmuR.data_path = "/data/birdcloudstorage-tvm/ibm-ml-refactor/data/")
options(murmuR.data_path = "/data/birdcloudstorage-tvm/ibm-ml/data/")

# Separate output path so writes never touch the read-only source data.
# Uses the same subdirectory structure as the source store.
out_path <- function(...) file.path("/data/birdcloudstorage-tvm/ibm-ml-refactor/data", ...)
future::plan(future::multisession)
progressr::handlers("progress")
options(future.globals.maxSize = 8000 * 1024^2)

# --- Shared constants ----------------------------------------
dataset      <- "hourly"
env_sampling <- c("dynamic", "static_stopover_enroute", "static_stopover", "radar")

metrics <- c(
  "regr.rmse", "regr.rsq",
  "regr.peaknights", "regr.peaknightslowalt",
  "regr.propmig", "regr.propmiglowalt",
  "regr.propmigrtt", "regr.propmiglowaltrtt",
  "regr.peakhours", "regr.peakhourslowalt",
  "regr.nightsneeded", "regr.nightsneededlowalt",
  "regr.nightsneededrtt", "regr.nightsneededlowaltrtt"
)

tuned_metrics <- c(
  "regr.rmse", "regr.rsq",
  "regr.peaknights", "regr.peaknightslowalt",
  "regr.propmigrtt", "regr.propmiglowaltrtt",
  "regr.peakhours", "regr.peakhourslowalt"
)

# ============================================================
# Helper functions
# ============================================================

annotate_peaknights <- function(preds, height_cutoff = FALSE, top_n = 10) {
  preds %>%
    { if (height_cutoff) filter(., height <= 250) else . } %>%
    group_by(datetime) %>%
    mutate(
      vid_pred  = sum(response) * 0.05,
      vid_truth = sum(truth) * 0.05
    ) %>%
    ungroup() %>%
    group_by(night_num) %>%
    summarise(
      sum_vid_pred  = sum(vid_pred),
      sum_vid_truth = sum(vid_truth),
      x_min = min(datetime),
      x_max = max(datetime)
    ) %>%
    mutate(
      vid_pred_rank  = rank(-sum_vid_pred,  ties.method = "min"),
      vid_truth_rank = rank(-sum_vid_truth, ties.method = "min"),
      peak_pred  = as.numeric(vid_pred_rank  <= top_n),
      peak_truth = as.numeric(vid_truth_rank <= top_n),
      peak_both  = as.integer(peak_pred & peak_truth)
    ) -> peaknights

  preds %>% left_join(peaknights)
}

plot_init_model_performance <- function(radar_metric_combo) {
  radar  <- radar_metric_combo$radar
  metric <- radar_metric_combo$metric

  plot_data <- bmr_scores %>%
    group_by(metric, type_human, learner_id, radar_human, season, sampling) %>%
    mutate(
      rank_score = case_when(
        metric == "regr.rmse"              ~ rank(score,  ties.method = "min"),
        metric == "regr.nightsneeded"      ~ rank(score,  ties.method = "min"),
        metric == "regr.nightsneededlowalt"~ rank(score,  ties.method = "min"),
        TRUE                               ~ rank(-score, ties.method = "min")
      )
    ) %>%
    filter(radar == !!radar, metric == !!metric) %>%
    ungroup() %>%
    mutate(learner_human = fct_reorder(learner_human, score, .fun = median, .desc = FALSE))

  if (nrow(plot_data) == 0) return(NULL)

  tweedie_medians <- plot_data %>%
    filter(learner_human == "Tweedie") %>%
    group_by(metric, season_human) %>%
    summarise(median_score = median(score), .groups = "drop")

  xlab <- if (unique(plot_data$metric) %in% c("regr.rmse", "regr.nightsneeded", "regr.nightsneededlowalt")) {
    "RMSE (lower is better)"
  } else {
    sprintf("%s (higher is better)", unique(plot_data$metric_human))
  }

  plot <- ggplot(plot_data, aes(x = score, y = learner_human, color = type_human)) +
    geom_boxplot() +
    scale_color_discrete(name = "Learner class") +
    geom_jitter(width = 0.1, alpha = 0.5) +
    geom_vline(data = tweedie_medians, aes(xintercept = median_score),
               linetype = "dashed", color = "black") +
    facet_grid(cols = vars(season_human), scales = "free", space = "free") +
    labs(x = xlab, y = "Learner")

  file_out <- sprintf(
    out_path("models/plots/init_models/%s_%s.pdf"),
    radar, str_replace(metric, "regr\\.", "")
  )
  ggsave(plot = plot, filename = file_out, width = 10, height = 10)
  return(file_out)
}

# ============================================================
# 1. Initial benchmark analysis
# ============================================================

bmr_outputs <- Sys.glob(murmuR_path("models/init_benchmarks/*_score.RDS"))
bmr_scores <- lapply(bmr_outputs, readRDS) %>%
  bind_rows() %>%
  rename(regr.rsq = rsq) %>%
  mutate(
    type_human = case_when(
      str_detect(learner_id, "quant")  ~ "Quantile",
      str_detect(learner_id, "poiss")  ~ "Poisson",
      str_detect(learner_id, "regrl1") ~ "L1 Regr.",
      str_detect(learner_id, "regr")   ~ "L2 Regr.",
      str_detect(learner_id, "twdie")  ~ "Tweedie",
      str_detect(learner_id, "gamma")  ~ "Gamma"
    ),
    trafo_human = case_when(
      str_detect(learner_id, "sqrt")  ~ "square root",
      str_detect(learner_id, "cbrt")  ~ "cube root",
      str_detect(learner_id, "log")   ~ "log",
      str_detect(learner_id, "log10") ~ "log10",
      .default = ""
    ),
    season_human  = str_to_title(season),
    trafo_human   = if_else(type_human == "Gamma" & trafo_human != "",
                            paste0(trafo_human, "+1"), trafo_human),
    learner_human = if_else(trafo_human != "",
                            sprintf("%s (%s)", type_human, trafo_human), type_human),
    radar_human   = case_when(
      radar == "nlhrw" ~ "Herwijnen",
      radar == "nldhl" ~ "Den Helder"
    ),
    sampling = case_when(sampling == "static" ~ "static_stopover_enroute", .default = sampling)
  ) %>%
  pivot_longer(cols = starts_with("regr"), names_to = "metric", values_to = "score") %>%
  mutate(metric_human = case_when(
    metric == "regr.rmse"                  ~ "RMSE",
    metric == "regr.rsq"                   ~ "R-squared",
    metric == "regr.peaknights"            ~ "# of peak nights",
    metric == "regr.peaknightslowalt"      ~ "# of peak nights (<250m)",
    metric == "regr.propmig"               ~ "% of migration",
    metric == "regr.propmiglowalt"         ~ "% of migration (<250m)",
    metric == "regr.propmigrtt"            ~ "% of peak night migration",
    metric == "regr.propmiglowaltrtt"      ~ "% of peak night migration (<250m)",
    metric == "regr.peakhours"             ~ "# of peak hours",
    metric == "regr.peakhourslowalt"       ~ "# of peak hours (<250m)",
    metric == "regr.nightsneeded"          ~ "Extra nights needed for 50% vs. truth",
    metric == "regr.nightsneededlowalt"    ~ "Extra nights needed for 50% vs. truth (<250m)",
    metric == "regr.nightsneededrtt"       ~ "Proportion real nights vs. prediction nights needed",
    metric == "regr.nightsneededrttlowalt" ~ "Proportion real nights vs. prediction nights needed (<250m)"
  ))

# Generate per-metric performance plots
jobs <- expand.grid(radar = names(radars_named), metric = metrics, stringsAsFactors = FALSE)
results <- future_lapply(seq_len(nrow(jobs)), function(i) {
  plot_init_model_performance(jobs[i, ])
})


# ============================================================
# 2. Feature selection analysis
# ============================================================

fs_archive_files <- Sys.glob(murmuR_path("models/feature_selection/summarised/*.RDS"))
fs_archive <- lapply(fs_archive_files, readRDS) %>%
  bind_rows() %>%
  rename(regr.rsq = rsq) %>%
  mutate(
    n_features  = as.numeric(n_features),
    sampling    = if_else(is.na(sampling), "radar", sampling),
    radar_human = case_when(
      radar == "nlhrw" ~ "Herwijnen",
      radar == "nldhl" ~ "Den Helder"
    ),
    season_human   = str_to_sentence(season),
    sampling_human = case_when(
      sampling == "radar"                  ~ "Radar only",
      sampling == "static_stopover"        ~ "Stat. Stopover",
      sampling == "static_stopover_enroute"~ "Stat. Stopover + Enroute",
      sampling == "dynamic"                ~ "Dyn. Stopover + Enroute"
    ),
    learner_id_human = case_when(
      learner_id == "lgbm_twdie"                                       ~ "No transformation",
      learner_id == "targetmutate.lgbm_twdie.cbrttargetinvert"         ~ "Cube-root transformed",
      learner_id == "targetmutate.lgbm_twdie.log10targetinvert"        ~ "Log10-transformed",
      learner_id == "targetmutate.lgbm_twdie.logtargetinvert"          ~ "Log-transformed",
      learner_id == "targetmutate.lgbm_twdie.sqrttargetinvert"         ~ "Square-root transformed"
    ),
    radar_human      = factor(radar_human,      levels = c("Den Helder", "Herwijnen")),
    season_human     = factor(season_human,     levels = c("Spring", "Autumn")),
    sampling_human   = factor(sampling_human,   levels = c(
      "Radar only", "Stat. Stopover", "Stat. Stopover + Enroute", "Dyn. Stopover + Enroute"
    )),
    learner_id_human = factor(learner_id_human, levels = c(
      "No transformation", "Square-root transformed", "Cube-root transformed",
      "Log-transformed", "Log10-transformed"
    ))
  )

# Performance vs. number of features
fs_archive %>%
  ggplot(aes(x = n_features, y = regr.nightsneededrtt, color = sampling_human)) +
  geom_line() + geom_point() +
  facet_wrap(vars(radar_human, season_human), nrow = 2, scales = "free") +
  labs(x = "Number of features in model", y = expression("RMSE [birds/km"^3*"]"))
ggsave(out_path("models/feature_selection/summarised/plots/fs_nr_features_rmse.pdf"),
       width = 20, height = 14)

ggplot(fs_archive, aes(x = n_features, y = regr.peaknights, color = learner_id_human)) +
  geom_smooth() + geom_point() +
  scale_color_discrete(name = "Target\nTransformation") +
  facet_wrap(vars(radar_human, season_human, sampling_human), scales = "free") +
  labs(x = "Number of features in model", y = "Identified peak nights [#]")
ggsave(out_path("models/feature_selection/summarised/plots/fs_nr_features_peaknights.pdf"),
       width = 20, height = 14)

ggplot(fs_archive, aes(x = n_features, y = regr.nightsneededlowaltrtt, color = learner_id_human)) +
  geom_smooth() + geom_point() +
  scale_color_discrete(name = "Target\nTransformation") +
  facet_wrap(vars(radar_human, season_human, sampling_human), scales = "free") +
  labs(x = "Number of features in model", y = "Identified peak hours")
ggsave(out_path("models/feature_selection/summarised/plots/fs_nr_features_peakhours.pdf"),
       width = 20, height = 14)

ggplot(fs_archive, aes(x = n_features, y = regr.propmigrtt, color = learner_id_human)) +
  geom_smooth() + geom_point() +
  scale_color_discrete(name = "Target\nTransformation") +
  facet_wrap(vars(radar_human, season_human, sampling_human), scales = "free") +
  labs(x = "Number of features in model",
       y = "Proportion of migration in predicted peak nights / actual peak nights [%]")
ggsave(out_path("models/feature_selection/summarised/plots/fs_nr_features_propmigrtt.pdf"),
       width = 20, height = 14)

ggplot(fs_archive, aes(x = n_features, y = regr.rsq, color = learner_id_human)) +
  geom_smooth() + geom_point() +
  scale_color_discrete(name = "Target\nTransformation") +
  facet_wrap(vars(radar_human, season_human, sampling_human), scales = "free") +
  labs(x = "Number of features in model", y = expression("R"^2*""))
ggsave(out_path("models/feature_selection/summarised/plots/fs_nr_features_rsquared.pdf"),
       width = 20, height = 14)

# Identify best-performing n_features ranges per group
fs_archive_betterperforming <- fs_archive %>%
  pivot_longer(cols = all_of(metrics), names_to = "metric", values_to = "score") %>%
  group_by(learner_id, metric, radar, season, sampling) %>%
  mutate(score_rank = case_when(
    metric == "regr.rmse"              ~ rank(score,  ties.method = "min"),
    metric == "regr.nightsneeded"      ~ rank(score,  ties.method = "min"),
    metric == "regr.nightsneededlowalt"~ rank(score,  ties.method = "min"),
    metric != "regr.rmse"              ~ rank(-score, ties.method = "min")
  )) %>%
  filter(score_rank <= 10)

fs_archive_betterperforming %>%
  filter(
    learner_id == "lgbm_twdie",
    metric %in% c("regr.rmse", "regr.rsq", "regr.peaknights",
                  "regr.propmigrtt", "regr.peakhours", "regr.nightsneededrtt")
  ) %>%
  ggplot(aes(x = interaction(season_human, sampling_human, sep = " - "),
             y = score, color = sampling_human)) +
  geom_boxplot(outlier.shape = NA, width = 0.4) +
  geom_jitter(height = 0.15, alpha = 0.5, size = 0.8) +
  scale_color_discrete(name = "Env. sampling") +
  coord_flip() +
  facet_wrap(vars(radar_human, metric), scales = "free_x", nrow = 2) +
  labs(y = "Score", x = "Environmental sampling method")
ggsave(out_path("models/feature_selection/summarised/plots/fs_sampling_method.pdf"),
       width = 20, height = 14)

# Top feature importances
fs_top_features <- fs_archive_betterperforming %>%
  filter(learner_id == "lgbm_twdie") %>%
  mutate(imp = purrr::map(importance, ~ enframe(.x, name = "feature", value = "importance"))) %>%
  select(-importance) %>%
  unnest(imp) %>%
  group_by(season_human, season, radar_human, radar, sampling_human, sampling, feature) %>%
  summarise(count = n(), feature_importance = sum(importance), .groups = "drop_last") %>%
  arrange(desc(feature_importance)) %>%
  mutate(rank = row_number()) %>%
  slice_max(order_by = feature_importance, n = 50) %>%
  group_by(season_human, season, radar_human, radar, sampling_human, sampling) %>%
  mutate(
    total_importance = feature_importance / sum(feature_importance),
    feature_type = case_when(
      str_starts(feature, "stp") ~ "Stop-over",
      str_starts(feature, "fl")  ~ "En-route",
      str_starts(feature, "rdr") ~ "Radar-level",
      feature == "height"        ~ "Height",
      .default                   = "Phenology"
    )
  )

saveRDS(fs_top_features, file = out_path("models/feature_selection/fs_top_features.RDS"))

fs_top_features %>%
  filter(rank <= 20) %>%
  ggplot() +
  geom_point(aes(
    x = total_importance * 100,
    y = tidytext::reorder_within(feature, total_importance,
                                 list(sampling_human, radar_human, season_human)),
    color = feature_type
  )) +
  facet_wrap(vars(sampling_human, radar_human, season_human), scales = "free_y") +
  tidytext::scale_y_reordered() +
  scale_color_discrete(name = "Feature type") +
  labs(x = "Relative importance [%]", y = "Top 20 features")
ggsave(out_path("models/feature_selection/summarised/plots/fs_top_features.pdf"),
       width = 20, height = 14)

# Save selected feature sets per group
fs_top_features %>%
  group_by(season, radar, sampling) %>%
  group_walk(~ {
    features <- .x$feature
    filename <- sprintf(
      out_path("models/hpo/fs_selected_top_features_%s_%s_%s.RDS"),
      .y$radar, .y$season, .y$sampling
    )
    saveRDS(features, file = filename)
  })


# ============================================================
# 3. Sampling benchmark
# ============================================================

for (radar in names(radars_named)) {
  for (season in names(months_seasons)) {
    for (sampling in env_sampling) {
      file_out <- sprintf(
        out_path("models/sampling_benchmarks/sampling_benchmark_tweedie_%s_%s_%s.RDS"),
        radar, season, sampling
      )

      task <- create_task(radar, season, sampling = sampling, dataset = "hourly")
      selected_features <- readRDS(sprintf(
        out_path("models/hpo/fs_selected_top_features_%s_%s_%s.RDS"),
        radar, season, sampling
      ))
      task$select(selected_features)

      learner <- lrn("regr.lightgbm",
                     objective  = "tweedie",
                     verbose    = 1,
                     num_threads = future::availableCores() - 1,
                     id         = paste0("lgbm_twdie_", sampling))

      years <- unique(task$data(cols = "year")$year)
      set.seed(42)
      cv <- rsmp("cv", folds = length(years))
      cv$instantiate(task)

      design <- benchmark_grid(task, learner, cv)

      set.seed(42)
      tic("Benchmark run")
      bmr <- benchmark(design, store_models = FALSE)
      toc()

      bmr_score <- bmr$score(measures = msrs(list(
        "regr.rmse", "regr.rsq",
        "regr.peaknights", "regr.peaknightslowalt",
        "regr.propmig", "regr.propmiglowalt",
        "regr.propmigrtt", "regr.propmiglowaltrtt",
        "regr.peakhours", "regr.peakhourslowalt"
      )))
      bmr_score$radar    <- radar
      bmr_score$season   <- season
      bmr_score$sampling <- sampling
      bmr_score$learner_id <- learner$id
      bmr_score$year     <- cv$instance$row_id

      bmr_score <- select(bmr_score, -c(uhash, task, learner, resampling, prediction_test))
      saveRDS(bmr_score, file_out)
    }
  }
}


# ============================================================
# 4. HPO results
# ============================================================

hpo_params <- tidyr::expand_grid(
  radar  = names(radars_named),
  season = names(months_seasons),
  metric = tuned_metrics
)

# Models tuned for their own metric
hpo_results <- pmap_dfr(hpo_params, function(radar, season, metric) {
  file_scores <- paste0(
    murmuR_path("models/hpo_mbo/fs_twdie_"),
    radar, "_", season, "_dynamic_", metric, "_score.RDS"
  )
  if (!file.exists(file_scores)) { warning(paste("Missing:", file_scores)); return(NULL) }

  alt_metric <- if_else(metric == "regr.rsq", "rsq", metric)
  scores <- readRDS(file_scores) %>% pull(!!alt_metric)
  if (str_detect(metric, "nightsneeded")) scores <- scores^(-1)

  score_mean <- mean(scores); score_min <- min(scores); score_max <- max(scores)
  display <- if_else(
    str_detect(metric, "nightsneeded"),
    paste0(round(score_mean, 2), "x |(", round(score_min, 2), "-", round(score_max, 2), ")"),
    paste0(round(score_mean, 2), "|(",   round(score_min, 2), "-", round(score_max, 2), ")")
  )
  tibble(radar = radar, season = season, metric = metric,
         score_mean = score_mean, score_min = score_min, score_max = score_max,
         display = display)
})

hpo_results %>%
  select(radar, season, metric, display) %>%
  pivot_wider(names_from = c(radar, season), values_from = display) %>%
  tt() %>%
  group_tt(j = list("Herwijnen" = 2:3, "Den Helder" = 4:5)) %>%
  print("latex")

# Models tuned for RMSE only (cross-metric comparison)
hpo_results_rmse <- pmap_dfr(hpo_params, function(radar, season, metric) {
  file_scores <- paste0(
    murmuR_path("models/hpo_mbo/fs_twdie_"),
    radar, "_", season, "_dynamic_regr.rmse_score.RDS"
  )
  if (!file.exists(file_scores)) { warning(paste("Missing:", file_scores)); return(NULL) }

  alt_metric <- if_else(metric == "regr.rsq", "rsq", metric)
  scores <- readRDS(file_scores) %>% pull(!!alt_metric)
  if (str_detect(metric, "nightsneeded")) scores <- scores^(-1)

  score_mean <- mean(scores); score_min <- min(scores); score_max <- max(scores)
  display <- if_else(
    str_detect(metric, "nightsneeded"),
    paste0(round(score_mean, 2), "x |(", round(score_min, 2), "-", round(score_max, 2), ")"),
    paste0(round(score_mean, 2), "|(",   round(score_min, 2), "-", round(score_max, 2), ")")
  )
  tibble(radar = radar, season = season, metric = metric,
         score_mean = score_mean, score_min = score_min, score_max = score_max,
         display = display)
})

hpo_results_rmse %>%
  select(radar, season, metric, display) %>%
  pivot_wider(names_from = c(radar, season), values_from = display) %>%
  tt() %>%
  group_tt(j = list("Herwijnen" = 2:3, "Den Helder" = 4:5)) %>%
  print("latex")


# ============================================================
# 5. Predictions and feature importances
# ============================================================

radars    <- c("nlhrw", "nldhl")
seasons   <- c("autumn", "spring")
samplings <- c("dynamic", "radar")

for (season in seasons) {
  for (radar in radars) {
    for (sampling in samplings) {
      importances_file_out <- sprintf(
        out_path("models/predictions/regr.rmse_%s/%s/%s_%s_%s_importances.RDS"),
        sampling, radar, radar, season, sampling
      )
      if (file.exists(importances_file_out)) next

      task <- create_task(radar, season, sampling = sampling, dataset = "hourly")
      train_features <- readRDS(sprintf(
        out_path("models/hpo/fs_selected_top_features_%s_%s_%s.RDS"),
        radar, season, sampling
      ))
      task$select(train_features)

      years <- unique(task$data(cols = "year")$year)
      set.seed(42)
      cv <- rsmp("cv", folds = length(years))
      cv$instantiate(task)

      for (cvid in seq_along(years)) {
        year    <- cv$instance$row_id[cvid]
        file_out <- sprintf(
          out_path("models/predictions/regr.rmse_%s/%s/%s_%s_%s_%s_preds.RDS"),
          sampling, radar, radar, season, year, sampling
        )
        if (file.exists(file_out)) next
        dir.create(dirname(file_out), recursive = TRUE, showWarnings = FALSE)

        original_data <- task$data(rows = cv$test_set(cvid), cols = task$col_info$id)
        message(sprintf("Model %s %s %s %s", radar, year, sampling, season))

        learner <- lrn("regr.lightgbm", objective = "tweedie", verbose = 2, num_threads = 10)
        best_params <- readRDS(sprintf(
          murmuR_path("models/hpo_mbo/fs_twdie_%s_%s_%s_regr.rmse_best_params_regr.rmse.RDS"),
          radar, season, "dynamic"
        ))
        learner$param_set$values <- best_params

        tic(paste0("Training model for ", radar, " ", season, " ", year))
        learner$train(task, row_ids = cv$train_set(cvid))
        toc()

        tic("Predicting & scoring")
        preds  <- learner$predict(task, row_ids = cv$test_set(cvid))
        scores <- preds$score(measures = msrs(metrics), task = task)
        toc()

        tic("Importance calculation")
        importance <- learner$importance()
        saveRDS(importance, file = str_replace(file_out, "_preds.RDS", "_importance.RDS"))
        toc()

        preds <- cbind(
          original_data,
          preds %>% as.data.table() %>% select(-c(row_ids)),
          as.data.frame(t(scores))
        )
        saveRDS(preds, file = file_out)
      }

      importances_files <- sprintf(
        out_path("models/predictions/regr.rmse_%s/%s/%s_%s_*_%s_importance.RDS"),
        sampling, radar, radar, season, sampling
      )
      importances          <- Sys.glob(importances_files) %>% lapply(readRDS) %>% bind_rows()
      importances$year     <- years
      importances$radar    <- radar
      importances$season   <- season
      importances$sampling <- sampling
      saveRDS(importances, file = importances_file_out)
    }
  }
}


# ============================================================
# 6. VPTS visualisations
# ============================================================

task  <- create_task("nlhrw", "autumn")
years <- unique(task$data(cols = "year")$year)

for (year in years) {
  preds_dynamic <- readRDS(sprintf(
    out_path("models/predictions/regr.rmse_dynamic/nlhrw/nlhrw_autumn_%s_dynamic_preds.RDS"), year
  )) %>%
    mutate(pred_type = "dynamic") %>%
    annotate_peaknights()

  preds_radar <- readRDS(sprintf(
    out_path("models/predictions/regr.rmse_radar/nlhrw/nlhrw_autumn_%s_radar_preds.RDS"), year
  )) %>%
    mutate(pred_type = "radar") %>%
    annotate_peaknights()

  plot_data_dynamic <- preds_dynamic %>%
    mutate(
      datetime    = as.factor(datetime),
      nightchange = c(0, diff(night_num)),
      peak_label  = factor(peak_pred + peak_both, levels = c(0, 1, 2),
                           labels = c("None", "Peak", "Both"))
    )

  plot_data_truth <- preds_dynamic %>%
    mutate(
      datetime    = as.factor(datetime),
      nightchange = c(0, diff(night_num)),
      peak_label  = factor(peak_truth, levels = c(0, 1, 2),
                           labels = c("None", "Peak", "Both"))
    )

  nightly_breaks_dynamic <- plot_data_dynamic %>%
    filter(nightchange > 0) %>%
    pull(datetime)

  date_labeller <- function(x) {
    dates     <- as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S")
    condition <- day(dates) %in% c(1, 15) | (day(dates) %in% 30 & month(dates) %in% 11)
    newdates  <- format(dates, "%b %d")
    newdates[!condition] <- ""
    newdates[newdates == "NA"] <- ""
    newdates
  }

  metrics_description_dynamic <- sprintf(
    "Peak nights: %d/10\nPeak hours: %d/100\nProp. of migration: %.0f%%",
    unique(plot_data_dynamic$regr.peaknights),
    unique(plot_data_dynamic$regr.peakhours),
    unique(plot_data_dynamic$regr.propmigrtt) * 100
  )

  vpts_theme <- list(
    theme_classic(base_size = 14),
    theme(
      panel.background        = element_rect(fill = "black"),
      panel.grid              = element_blank(),
      legend.position.inside  = c(.98, .98),
      legend.justification    = c("right", "top"),
      legend.background       = element_blank(),
      legend.box.background   = element_blank(),
      legend.key              = element_blank(),
      legend.text             = element_text(color = "white"),
      legend.title            = element_text(color = "white", vjust = 0.75),
      legend.key.width        = unit(2, "lines"),
      axis.ticks              = element_line(colour = "black", linewidth = 0.5)
    )
  )

  plot_dynamic <- ggplot(plot_data_dynamic) +
    geom_tile(aes(x = datetime, y = height, fill = response)) +
    geom_point(aes(x = datetime, color = peak_label), y = 0, shape = "|", size = 3) +
    scale_fill_viridis_c(trans = "sqrt", option = "magma", limits = c(0, 300), na.value = 0,
                         name = expression("Bird density [" * birds ~ km^{-3} * "]"),
                         oob = scales::squish) +
    scale_x_discrete(breaks = nightly_breaks_dynamic, labels = date_labeller) +
    annotate("text_npc", npcx = 0.01, npcy = 0.95, label = "Predicted",
             color = "white", size = 8, hjust = 0, vjust = 1, fontface = "bold") +
    annotate("text_npc", npcx = 0.01, npcy = 0.80, label = metrics_description_dynamic,
             color = "white", size = 5, hjust = 0, vjust = 1) +
    scale_color_manual(values = c("None" = "black", "Peak" = "red", "Both" = "yellow"),
                       name = "Migration peak", guide = "none") +
    vpts_theme +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
    guides(fill = guide_colorbar(position = "inside", direction = "horizontal")) +
    coord_cartesian(ylim = c(0, 3000), expand = FALSE) +
    labs(x = "", y = "Height [m]")

  plot_truth <- ggplot(plot_data_truth) +
    geom_tile(aes(x = datetime, y = height, fill = truth)) +
    geom_point(aes(x = datetime, color = peak_label), y = 0, shape = "|", size = 3) +
    scale_fill_viridis_c(trans = "sqrt", option = "magma", limits = c(0, 300),
                         name = expression("Bird density [" * birds ~ km^{-3} * "]"),
                         guide = "none", oob = scales::squish) +
    scale_x_discrete(breaks = nightly_breaks_dynamic, labels = date_labeller) +
    annotate("text_npc", npcx = 0.01, npcy = 0.95, label = "Measured",
             color = "white", size = 8, hjust = 0, vjust = 1, fontface = "bold") +
    scale_color_manual(values = c("None" = "black", "Peak" = "yellow", "Both" = "red"),
                       name = "Migration peak", guide = "none") +
    vpts_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, 3000), expand = FALSE) +
    labs(x = "", y = "Height [m]")

  plot_vpts <- plot_dynamic / plot_truth + plot_layout(axis_titles = "collect")
  ggsave(
    filename = sprintf(out_path("models/plots/vpts/vpts_%s_dynamic_vs_radar.pdf"), year),
    plot = plot_vpts, width = 12, height = 8
  )
}


# ============================================================
# 7. ROC curves
# ============================================================

preds_all <- Sys.glob(out_path("models/predictions/regr.rmse_dynamic/nlhrw/*.RDS")) %>%
  future_lapply(readRDS) %>%
  bind_rows()

preds_all %>%
  select(year, height, truth, response, datetime, night_num) %>%
  mutate(year = as.factor(year)) %>%
  group_by(year, night_num, datetime) %>%
  summarise(
    vid_pred  = sum(response) * 0.05,
    vid_truth = sum(truth) * 0.05
  ) %>%
  ggplot() +
  geom_roc(aes(d = vid_truth > 25, m = vid_pred, color = year), n.cuts = 0) +
  style_roc()
