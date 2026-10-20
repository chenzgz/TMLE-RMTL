# TMLE-RMTL

This repository contains R code for targeted maximum likelihood estimation (TMLE) of treatment effects on restricted mean time lost (RMTL) under competing risks.

The main function is `contmle.rmtl()`. It implements continuous-time one-step TMLE using cause-specific hazard models, inverse probability of censoring weighting, treatment weighting, and efficient influence curve-based standard errors.

## Installation

This code is provided as analysis code rather than a formal R package. Clone or download the repository, then source the main R file:

```r
source("R/TMLE-RMTL.R")
```

Required R packages:

```r
install.packages(c(
  "data.table",
  "survival",
  "Matrix",
  "zoo",
  "glmnet",
  "stringr"
))
```

## Data Format

The input dataset should contain one row per subject and include:

| Variable | Description |
| --- | --- |
| `time` | Observed follow-up time |
| `delta` | Event indicator: `0` for censoring, `1` and `2` for competing event types |
| `A` | Binary treatment or exposure indicator |
| Covariates | Baseline covariates used in the treatment, event, and censoring models |

The default example assumes covariates named `L1`, `L2`, and `L3`, but these can be changed through the model formulas.

## Basic Usage

```r
library(data.table)
source("R/TMLE-RMTL.R")
source("examples/simulation_dgm.R")

example_data <- sim.data(
  n = 500,
  seed = 20261013,
  setting = 1,
  censoring = 0.15,
  censoring.informative = FALSE
)

fit <- contmle.rmtl(
  dt = example_data,
  estimation = list(
    event1 = list(
      fit = "cox",
      model = Surv(time, delta == 1) ~ A + L1 + L2 + L3
    ),
    event2 = list(
      fit = "cox",
      model = Surv(time, delta == 2) ~ A + L1 + L2 + L3
    ),
    cens = list(
      fit = "cox",
      model = Surv(time, delta == 0) ~ 1
    )
  ),
  treat.model = A ~ L1 + L2 + L3,
  target = 1:2,
  tau = 2,
  treat.effect = "ate"
)
```

See `examples/simulation_dgm.R` for the simulation data-generating mechanisms used in the manuscript settings, and `examples/example_usage.R` for a minimal reproducible analysis example.

## Main Arguments

| Argument | Description |
| --- | --- |
| `dt` | Input dataset as a `data.table` or data frame |
| `estimation` | List specifying initial models for event 1, event 2, and censoring |
| `target` | Event type for the RMTL estimand |
| `treat.model` | Logistic regression model for the treatment mechanism |
| `treat.effect` | Treatment contrast: `"ate"`, `"1"`, or `"0"` |
| `tau` | Restriction time point(s) for RMTL |
| `cut.off.cens` | Lower truncation value for censoring survival probabilities |
| `hal.tuning` | Whether to use cross-validated HAL tuning |
| `verbose` | Whether to print model-fitting details |

## Notes

- The public interface keeps the final analysis options used in the manuscript and removes development-only switches.
- Initial models can be specified as Cox, HAL, Cox-HAL, or Kaplan-Meier style fits, depending on the value of `fit`.
- The repository is intended for reproducible research use. For a formal R package release, the code can later be reorganized with roxygen2 documentation and unit tests.
