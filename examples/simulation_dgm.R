library(data.table)

sim.data1 <- function(n,
                      betaA1 = log(0.7),
                      betaA2 = log(1.2),
                      nu = 1.5,
                      eta = 0.8,
                      censoring = 0,
                      seed = sample(4034244, 1),
                      verbose = FALSE,
                      censoring.informative = FALSE,
                      intervention.A = NULL) {
  set.seed(seed)

  if (censoring == 0.15) {
    censoring.alpha <- -1.1
  } else if (censoring == 0.3) {
    censoring.alpha <- 0.02
  }

  if (length(intervention.A) > 0) {
    censoring <- 0
  }

  L1 <- runif(n, 0, 1)
  L2 <- rnorm(n, 0, 1)
  L3 <- rnorm(n, 0, 1)

  A <- rbinom(n, 1, plogis(
    -0.11 + log(0.8) * L1 + log(2) * L2 + log(1.5) * L3
  ))

  if (verbose) {
    message("Treatment model: logit(P(A=1)) = -0.11 + log(0.8)*L1 + log(2)*L2 + log(1.5)*L3")
  }

  if (length(intervention.A) > 0 && is.numeric(intervention.A)) {
    A <- intervention.A
  }

  phiT1 <- function(t, A, L1, L2, L3, betaA1) {
    exp(betaA1 * A + log(2.5) * L1 + log(2) * L2 + log(1.5) * L3)
  }

  lambdaT1 <- function(t, A, L1, L2, L3, betaA1, eta, nu) {
    phiT1(t, A, L1, L2, L3, betaA1) * eta * nu * t^(nu - 1)
  }

  if (censoring > 0) {
    if (!censoring.informative) {
      phiC <- function(t, A, L1, L2, L3) {
        exp(censoring.alpha)
      }
    } else {
      phiC <- function(t, A, L1, L2, L3) {
        exp(log(2) * L1 + log(0.5) * L3 - 0.68 + censoring.alpha)
      }
    }
  }

  lambdaC <- function(t, A, L1, L2, L3, eta, nu) {
    phiC(t, A, L1, L2, L3) * eta * nu * t^(nu - 1)
  }

  phiT2 <- function(t, A, L1, L2, L3, betaA2) {
    exp(betaA2 * A + log(0.6) * L1 + log(1.6) * L2 + log(2) * L3)
  }

  lambdaT2 <- function(t, A, L1, L2, L3, betaA2, eta, nu) {
    phiT2(t, A, L1, L2, L3, betaA2) * eta * nu * t^(nu - 1)
  }

  if (censoring > 0) {
    phi <- function(t, A, L1, L2, L3, betaA1, betaA2) {
      phiT1(t, A, L1, L2, L3, betaA1) +
        phiC(t, A, L1, L2, L3) +
        phiT2(t, A, L1, L2, L3, betaA2)
    }
  } else {
    phi <- function(t, A, L1, L2, L3, betaA1, betaA2) {
      phiT1(t, A, L1, L2, L3, betaA1) +
        phiT2(t, A, L1, L2, L3, betaA2)
    }
  }

  Lambda.inv <- function(u, t, A, L1, L2, L3, nu, eta) {
    ((u + eta * phi(t, A, L1, L2, L3, betaA1, betaA2) * t^nu) /
       (eta * phi(t, A, L1, L2, L3, betaA1, betaA2)))^(1 / nu) - t
  }

  T.lat <- Lambda.inv(rexp(n), 0, A, L1, L2, L3, nu, eta)
  C.admin <- 5
  T.obs <- pmin(T.lat, C.admin)

  h1 <- lambdaT1(T.lat, A, L1, L2, L3, betaA1, eta, nu)
  h2 <- lambdaT2(T.lat, A, L1, L2, L3, betaA2, eta, nu)
  hc <- if (censoring > 0) lambdaC(T.lat, A, L1, L2, L3, eta, nu) else rep(0, n)

  h.total <- h1 + h2 + hc
  probs <- runif(n)
  delta <- rep(0, n)

  p1 <- h1 / h.total
  p2 <- h2 / h.total

  delta[probs <= p1] <- 1
  delta[probs > p1 & probs <= p1 + p2] <- 2
  delta[T.lat > C.admin] <- 0

  data.table(time = T.obs, delta = delta, A = A, L1 = L1, L2 = L2, L3 = L3)
}

sim.data2 <- function(n,
                      betaA1.pre = log(1.20),
                      betaA1.post = log(0.80),
                      betaA2 = log(1.50),
                      nu = 1.5,
                      eta = 0.8,
                      censoring = 0,
                      seed = sample(4034244, 1),
                      t0 = 0.40,
                      verbose = FALSE,
                      censoring.informative = FALSE,
                      intervention.A = NULL) {
  set.seed(seed)

  if (censoring == 0.15 && !censoring.informative) {
    censoring.alpha <- -0.911457840
  } else if (censoring == 0.3 && !censoring.informative) {
    censoring.alpha <- 0.085848638
  } else if (censoring == 0.15 && censoring.informative) {
    censoring.alpha <- -1.335163174
  } else if (censoring == 0.3 && censoring.informative) {
    censoring.alpha <- -0.231629312
  }

  if (length(intervention.A) > 0) {
    censoring <- 0
  }

  L1 <- runif(n, 0, 1)
  L2 <- rnorm(n, 0, 1)
  L3 <- rnorm(n, 0, 1)
  z2sq <- pmin(L2^2, 4) - 1
  z3hi <- as.numeric(L3 > 0.55)
  z1hi <- as.numeric(L1 > 0.65)
  int_l2_l3 <- L2 * L3

  A <- rbinom(n, 1, plogis(
    0.289 +
      log(0.80) * L1 +
      log(2.00) * L2 +
      log(1.50) * L3 +
      log(1.50) * z2sq +
      log(1.00) * z3hi +
      log(0.80) * z1hi +
      log(1.20) * int_l2_l3
  ))

  if (verbose) {
    message("Treatment model uses L1, L2, L3, a truncated L2 quadratic term, threshold terms, and L2-L3 interaction.")
  }

  if (length(intervention.A) > 0 && is.numeric(intervention.A)) {
    A <- intervention.A
  }

  phiT1.fixed <- function(A, L1, L2, L3, betaA1.current) {
    z2sq <- pmin(L2^2, 4) - 1
    z3hi <- as.numeric(L3 > 0.55)
    z1hi <- as.numeric(L1 > 0.65)
    int_l2_l3 <- L2 * L3

    exp(betaA1.current * A +
          log(2.50) * L1 +
          log(2.00) * L2 +
          log(1.00) * z1hi +
          log(0.80) * z2sq +
          log(1.00) * z3hi +
          log(1.50) * L3 +
          log(1.00) * int_l2_l3)
  }

  phiT1 <- function(t, A, L1, L2, L3, betaA1.pre, betaA1.post) {
    betaA1.current <- ifelse(t <= t0, betaA1.pre, betaA1.post)
    phiT1.fixed(A, L1, L2, L3, betaA1.current)
  }

  lambdaT1 <- function(t, A, L1, L2, L3, betaA1.pre, betaA1.post, eta, nu) {
    phiT1(t, A, L1, L2, L3, betaA1.pre, betaA1.post) * eta * nu * t^(nu - 1)
  }

  if (censoring > 0) {
    if (!censoring.informative) {
      phiC <- function(t, A, L1, L2, L3) {
        exp(censoring.alpha)
      }
    } else {
      phiC <- function(t, A, L1, L2, L3) {
        exp(log(1.5) * L1 + log(0.6) * L3 + censoring.alpha)
      }
    }
  }

  lambdaC <- function(t, A, L1, L2, L3, eta, nu) {
    phiC(t, A, L1, L2, L3) * eta * nu * t^(nu - 1)
  }

  phiT2 <- function(t, A, L1, L2, L3, betaA2) {
    z2sq <- pmin(L2^2, 4) - 1
    z3hi <- as.numeric(L3 > 0.55)
    z1hi <- as.numeric(L1 > 0.65)
    int_l2_l3 <- L2 * L3

    exp(betaA2 * A +
          log(0.60) * L1 +
          log(1.50) * z1hi +
          log(1.20) * z2sq +
          log(1.00) * z3hi +
          log(1.60) * L2 +
          log(2.00) * L3 +
          log(1.00) * int_l2_l3)
  }

  lambdaT2 <- function(t, A, L1, L2, L3, betaA2, eta, nu) {
    phiT2(t, A, L1, L2, L3, betaA2) * eta * nu * t^(nu - 1)
  }

  if (censoring > 0) {
    phi <- function(A, L1, L2, L3, betaA1.current, betaA2) {
      phiT1.fixed(A, L1, L2, L3, betaA1.current) +
        phiC(0, A, L1, L2, L3) +
        phiT2(0, A, L1, L2, L3, betaA2)
    }
  } else {
    phi <- function(A, L1, L2, L3, betaA1.current, betaA2) {
      phiT1.fixed(A, L1, L2, L3, betaA1.current) +
        phiT2(0, A, L1, L2, L3, betaA2)
    }
  }

  Lambda.inv <- function(u, t, A, L1, L2, L3, nu, eta) {
    if (any(t != 0)) {
      stop("Lambda.inv is implemented for inversion from time 0 in this simulation.")
    }

    phi.pre <- phi(A, L1, L2, L3, betaA1.pre, betaA2)
    phi.post <- phi(A, L1, L2, L3, betaA1.post, betaA2)
    cutoff <- eta * phi.pre * t0^nu

    T.lat <- numeric(length(u))
    before.t0 <- u <= cutoff
    T.lat[before.t0] <- (u[before.t0] / (eta * phi.pre[before.t0]))^(1 / nu)
    T.lat[!before.t0] <- ((u[!before.t0] -
                             eta * phi.pre[!before.t0] * t0^nu +
                             eta * phi.post[!before.t0] * t0^nu) /
                            (eta * phi.post[!before.t0]))^(1 / nu)

    T.lat
  }

  T.lat <- Lambda.inv(rexp(n), 0, A, L1, L2, L3, nu, eta)
  C.admin <- 5
  T.obs <- pmin(T.lat, C.admin)

  h1 <- lambdaT1(T.lat, A, L1, L2, L3, betaA1.pre, betaA1.post, eta, nu)
  h2 <- lambdaT2(T.lat, A, L1, L2, L3, betaA2, eta, nu)
  hc <- if (censoring > 0) lambdaC(T.lat, A, L1, L2, L3, eta, nu) else rep(0, n)

  h.total <- h1 + h2 + hc
  probs <- runif(n)
  delta <- rep(0, n)

  p1 <- h1 / h.total
  p2 <- h2 / h.total

  delta[probs <= p1] <- 1
  delta[probs > p1 & probs <= p1 + p2] <- 2
  delta[T.lat > C.admin] <- 0

  data.table(time = T.obs, delta = delta, A = A, L1 = L1, L2 = L2, L3 = L3)
}

sim.data <- function(n,
                     censoring = 0,
                     seed = sample(4034244, 1),
                     setting = 1,
                     verbose = FALSE,
                     censoring.informative = FALSE,
                     intervention.A = NULL) {
  if (setting == 1) {
    return(sim.data1(
      n = n,
      censoring = censoring,
      seed = seed,
      verbose = verbose,
      censoring.informative = censoring.informative,
      intervention.A = intervention.A
    ))
  }

  if (setting == 2) {
    return(sim.data2(
      n = n,
      censoring = censoring,
      seed = seed,
      verbose = verbose,
      censoring.informative = censoring.informative,
      intervention.A = intervention.A
    ))
  }

  stop("setting must be 1 or 2.")
}
