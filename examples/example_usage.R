library(data.table)
library(survival)

source("R/TMLE-RMTL.R")
source("examples/simulation_dgm.R")

set.seed(20261013)

setting <- 1
censoring.rate <- 0.15
censoring.informative <- FALSE
tau <- 2

example_data <- sim.data(
  n = 500,
  seed = 20261013,
  setting = setting,
  censoring = censoring.rate,
  censoring.informative = censoring.informative
)

cens.model <- if (censoring.informative) {
  Surv(time, delta == 0) ~ A + L1 + L2 + L3
} else {
  Surv(time, delta == 0) ~ 1
}

estimation <- list(
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
    model = cens.model
  )
)

fit <- contmle.rmtl(
  dt = example_data,
  estimation = estimation,
  target = 1:2,
  treat.model = A ~ L1 + L2 + L3,
  treat.effect = "ate",
  tau = tau,
  cut.off.cens = 1e-2,
  verbose = FALSE
)

print(fit$tmle)
