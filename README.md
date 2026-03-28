# murmuR-paper

Analysis and figure scripts accompanying the paper:

> **Simulation-informed forecasts of peak nocturnal bird migration events to reduce human-wildlife conflicts**
>
> Bart Hoekstra, Emiel van Loon, Judy Shamoun-Baranes

This repository contains the code to reproduce all analyses and figures presented in the paper. Most scripts depend on the [murmuR](https://github.com/barthoekstra/murmuR) R package (v1.0.0) for data paths, constants, and model infrastructure.

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
| `benchmark_ibm.R` | IBM simulation runtime benchmark (Appendix) |
| `comparison_roy_et_al.R` | Harmonized comparison with Roy et al. (2025) (Appendix) |
| `sensitivity_nbirds.R` | Predictor stability as a function of IBM population size (Appendix) |

### `figures/`

| Script | Description |
|--------|-------------|
| `figure_sampling_comparison.R` | Comparison of environmental sampling strategies (Figure 1) |
| `figure_pred_vs_obs.R` | Predicted vs. observed bird density scatter plots (Figure 2) |
| `figure_vpts.R` | Vertical profile time series of bird densities (Figure 3) |
| `figure_variable_importance.R` | Variable importance across predictor groups (Figure 4) |
| `figure_cum_mtr_ranked_hours.R` | Cumulative migration traffic by ranked hours |
| `figure_roc_curves.R` | ROC curves for peak night classification |
| `figure_catchment_areas.R` | Radar catchment area map (Appendix) |
| `custom_geom_stripes.R` | Custom ggplot2 geom for alternating background stripes |
