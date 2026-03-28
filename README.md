# murmuR-paper

Analysis and figure scripts accompanying the paper:

> **Simulation-informed forecasts of peak nocturnal bird migration events to reduce human-wildlife conflicts**
>
> Bart Hoekstra, Emiel van Loon, Judy Shamoun-Baranes

This repository contains the code to reproduce all analyses and figures presented in the paper. The core modeling framework is implemented in the [murmuR](https://github.com/barthoekstra/murmuR) R package, which this repository depends on.

## Setup

### Install murmuR

```r
# install.packages("remotes")
remotes::install_github("barthoekstra/murmuR@v1.0.0")
```

### Configure data path

All scripts expect the data path to be set before running:

```r
library(murmuR)
options(murmuR.data_path = "/path/to/your/data/")
```

## Repository structure

### `analysis/`

| Script | Description |
|--------|-------------|
| `training.R` | Full training pipeline: initial benchmarks, feature selection, sampling benchmarks, HPO, predictions, variable importance, VPTS visualisations, ROC analysis |
| `benchmark_ibm.R` | IBM simulation runtime benchmark (Appendix Table) |
| `comparison_roy_et_al.R` | Harmonized comparison with Roy et al. (2025) (Appendix) |
| `sensitivity_nbirds.R` | Predictor stability as a function of IBM population size (Appendix) |
| `groundspeed_summary.R` | Groundspeed statistics of simulated migrants |
| `comparison_pipeline.R` | Verification of refactored murmuR package against original outputs |
| `comparison_rerun.R` | Re-run of comparison pipeline after bug fixes |

### `figures/`

| Script | Description |
|--------|-------------|
| `figure_vpts.R` | Vertical profile time series of bird densities (Figure 3) |
| `figure_sampling_comparison.R` | Comparison of environmental sampling strategies (Figure 1) |
| `figure_variable_importance.R` | Variable importance across predictor groups (Figure 4) |
| `figure_cum_mtr_ranked_hours.R` | Cumulative migration traffic by ranked hours |
| `figure_roc_curves.R` | ROC curves for peak night classification |
| `figure_pred_vs_obs.R` | Predicted vs. observed bird density scatter plots (Figure 2) |
| `figure_catchment_areas.R` | Radar catchment area map (Appendix Figure) |
| `custom_geom_stripes.R` | Custom ggplot2 geom for alternating background stripes |
