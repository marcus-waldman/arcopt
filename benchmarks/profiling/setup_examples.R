# benchmarks/profiling/setup_examples.R
# ======================================
# Defines three problem builders that each return a list with elements
#   x0   : starting vector
#   fn   : objective (numeric)
#   gr   : gradient (numeric vector)
#   hess : Hessian (matrix)
#   n    : length(x0)
#   label: short identifier used in result filenames
#
# These wrap the same problems used in the JSS manuscript (mixture
# saddle, GMM, DPM Heckman MAP) so the profiling sweep matches what
# users actually run. Mixture and GMM are lifted directly from the
# manuscript chunks; DPM Heckman wraps benchmarks/dpm_heckman/setup.R.

# --------------------------------------------------------------------
# Mixture saddle (n = 2)
# Source: manuscript/arcopt-jss.qmd chunk `mixture-setup`
# --------------------------------------------------------------------
make_mixture_problem <- function(seed = 42L) {
  set.seed(seed)
  n <- 200
  z <- rbinom(n, 1, 0.5)
  y <- ifelse(z == 0, rnorm(n, -2, 1), rnorm(n, 2, 1))

  mixture_nll <- function(theta) {
    mu1 <- theta[1]
    mu2 <- theta[2]
    log_lik <- sapply(y, function(yi) {
      log(0.5 * dnorm(yi, mu1, 1) + 0.5 * dnorm(yi, mu2, 1))
    })
    -sum(log_lik)
  }

  mixture_grad <- function(theta) {
    mu1 <- theta[1]
    mu2 <- theta[2]
    g <- c(0, 0)
    for (yi in y) {
      p1 <- 0.5 * dnorm(yi, mu1, 1)
      p2 <- 0.5 * dnorm(yi, mu2, 1)
      w1 <- p1 / (p1 + p2)
      w2 <- p2 / (p1 + p2)
      g[1] <- g[1] + w1 * (yi - mu1)
      g[2] <- g[2] + w2 * (yi - mu2)
    }
    -g
  }

  mixture_hess <- function(theta, h = 1e-5) {
    d <- length(theta)
    H <- matrix(0, d, d)
    for (i in seq_len(d)) {
      ei <- rep(0, d)
      ei[i] <- h
      gp <- mixture_grad(theta + ei)
      gm <- mixture_grad(theta - ei)
      H[, i] <- (gp - gm) / (2 * h)
    }
    0.5 * (H + t(H))
  }

  list(
    label = "mixture",
    n = 2L,
    x0 = c(0.01, 0.01),
    fn = mixture_nll,
    gr = mixture_grad,
    hess = mixture_hess
  )
}

# --------------------------------------------------------------------
# Growth Mixture Model saddle (n = 6)
# Source: manuscript/arcopt-jss.qmd chunk `gmm-setup`
# --------------------------------------------------------------------
make_gmm_problem <- function(seed = 42L) {
  set.seed(seed)
  N <- 40
  T_obs <- 6
  t_vec <- seq(0, 5, length.out = T_obs)

  b01_true <- 5
  b11_true <- -0.5
  b02_true <- 0
  b12_true <- 2.0
  pi_true <- 0.4
  sigma_true <- 0.5

  z <- rbinom(N, 1, 1 - pi_true)
  Y <- matrix(NA, N, T_obs)
  for (i in seq_len(N)) {
    mu_i <- if (z[i] == 0) b01_true + b11_true * t_vec
            else b02_true + b12_true * t_vec
    Y[i, ] <- mu_i + rnorm(T_obs, 0, sigma_true)
  }

  gmm_nll <- function(theta) {
    if (any(!is.finite(theta))) return(1e12)
    b01 <- theta[1]; b11 <- theta[2]
    b02 <- theta[3]; b12 <- theta[4]
    logit_pi <- theta[5]; log_sigma <- theta[6]
    pi_val <- 1 / (1 + exp(-logit_pi)); sigma <- exp(log_sigma)
    mu1 <- b01 + b11 * t_vec; mu2 <- b02 + b12 * t_vec
    log_phi1 <- rowSums(dnorm(Y, matrix(mu1, N, T_obs, byrow = TRUE),
                              sigma, log = TRUE))
    log_phi2 <- rowSums(dnorm(Y, matrix(mu2, N, T_obs, byrow = TRUE),
                              sigma, log = TRUE))
    a <- log(pi_val) + log_phi1
    b <- log(1 - pi_val) + log_phi2
    mx <- pmax(a, b)
    ll_i <- mx + log(exp(a - mx) + exp(b - mx))
    -sum(ll_i)
  }

  gmm_grad <- function(theta, h = 1e-5) {
    d <- length(theta); g <- numeric(d)
    for (i in seq_len(d)) {
      ei <- rep(0, d); ei[i] <- h
      g[i] <- (gmm_nll(theta + ei) - gmm_nll(theta - ei)) / (2 * h)
    }
    g
  }

  gmm_hess <- function(theta, h = 1e-4) {
    d <- length(theta); H <- matrix(0, d, d)
    for (i in seq_len(d)) {
      ei <- rep(0, d); ei[i] <- h
      H[, i] <- (gmm_grad(theta + ei) - gmm_grad(theta - ei)) / (2 * h)
    }
    0.5 * (H + t(H))
  }

  grand_int <- mean(Y[, 1])
  grand_slope <- (mean(Y[, T_obs]) - mean(Y[, 1])) /
                 (t_vec[T_obs] - t_vec[1])
  x0 <- c(grand_int, grand_slope, grand_int, grand_slope, 0, log(1))

  list(
    label = "gmm",
    n = 6L,
    x0 = x0,
    fn = gmm_nll,
    gr = gmm_grad,
    hess = gmm_hess
  )
}

# --------------------------------------------------------------------
# DPM Heckman MAP (n = 299)
# Source: benchmarks/dpm_heckman/{setup.R, cache_for_manuscript.R}
# Uses BridgeStan-compiled Stan model and AD Hessian, so the profile
# reflects arcopt-internal time -- not user R callbacks.
# --------------------------------------------------------------------
make_dpm_heckman_problem <- function(seed = 101L,
                                     n_obs = 500L,
                                     k_true = 2L,
                                     p_out = 5L,
                                     p_sel = 7L,
                                     n_excl = 2L,
                                     k_trunc = 20L,
                                     setup_path = file.path(
                                       "benchmarks", "dpm_heckman",
                                       "setup.R")) {
  if (!requireNamespace("bridgestan", quietly = TRUE)) {
    stop("DPM Heckman problem requires the 'bridgestan' package")
  }
  if (!file.exists(setup_path)) {
    stop("could not find ", setup_path,
         " (run from repo root, or supply setup_path=)")
  }
  source(setup_path, local = FALSE)
  model_so <- compile_dpm_heckman_model()
  data <- simulate_dpm_heckman(n_obs = n_obs, k_true = k_true,
                               p_out = p_out, p_sel = p_sel,
                               n_excl = n_excl, seed = seed)
  fns <- make_dpm_heckman_fns(data, model_so, k_trunc = k_trunc)

  list(
    label = "dpm_heckman",
    n = fns$n_params,
    x0 = fns$x0,
    fn = fns$fn,
    gr = fns$gr,
    hess = fns$hess_ad
  )
}
