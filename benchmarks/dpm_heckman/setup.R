# DPM Mixture Heckman — simulation + BridgeStan wrapper
# =====================================================
#
# Simulates data from a K_true-component mixture Heckman model with
# exclusion restriction, then returns fn / gr / hess_ad / hess_fd
# callables suitable for passing to arcopt / optim / nlminb / trust.

suppressPackageStartupMessages({
  library(bridgestan)
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("package 'jsonlite' required")
  }
})

# --------------------------------------------------------------------
# Simulation
#
# Generates:
#   X: N x P_out covariate matrix for the outcome equation, with
#      intercept in column 1
#   W: N x P_sel covariate matrix for the selection equation, with
#      intercept in column 1; the last two columns are the exclusion
#      restriction (do not appear in X)
#   Z: selection indicator
#   Y: outcome, valid for Z=1 rows (set to 0 otherwise)
#   class: true component label
# True parameters:
#   beta_k: P_out vector per class, with substantial between-class
#           separation in the non-intercept coefficients
#   gamma_k: P_sel vector per class; last two entries are the
#            exclusion-restriction coefficients (nonzero)
#   sigma_k = 1 for all k
#   rho_k in {-0.3, 0.3, ...} alternating
# --------------------------------------------------------------------
simulate_dpm_heckman <- function(n_obs = 500L,
                                 k_true = 2L,
                                 p_out = 5L,
                                 p_sel = 7L,
                                 n_excl = 2L,
                                 seed = 1L) {
  stopifnot(p_sel >= p_out + n_excl - 0L)   # at least n_excl extra cols in W
  set.seed(seed)

  # Class sizes roughly equal
  class_probs <- rep(1 / k_true, k_true)
  class_labels <- sample.int(k_true, n_obs, replace = TRUE,
                             prob = class_probs)

  # Shared covariates (intercept + noise)
  x_noise <- matrix(stats::rnorm(n_obs * (p_out - 1L)), n_obs, p_out - 1L)
  X <- cbind(1, x_noise)
  # W copies the first (p_sel - n_excl) columns of X (so those are
  # common) and appends n_excl exclusion-only columns
  shared_cols <- p_sel - n_excl
  if (shared_cols > p_out) {
    # Not enough columns in X to cover shared cols; pad with noise
    extra_shared <- matrix(stats::rnorm(n_obs * (shared_cols - p_out)),
                           n_obs, shared_cols - p_out)
    w_shared <- cbind(X, extra_shared)
  } else {
    w_shared <- X[, seq_len(shared_cols), drop = FALSE]
  }
  w_excl <- matrix(stats::rnorm(n_obs * n_excl), n_obs, n_excl)
  W <- cbind(w_shared, w_excl)

  # True beta, gamma: well-separated across classes
  beta_true <- matrix(0, k_true, p_out)
  gamma_true <- matrix(0, k_true, p_sel)
  for (k in seq_len(k_true)) {
    beta_true[k, ] <- stats::rnorm(p_out, sd = 0.8) *
      ((-1)^(k - 1L)) * (1 + (k - 1L) * 0.7)
    gamma_true[k, ] <- stats::rnorm(p_sel, sd = 0.6) *
      ((-1)^(k - 1L)) * (1 + (k - 1L) * 0.5)
    # Ensure exclusion-restriction coefficients are substantively nonzero
    gamma_true[k, (p_sel - n_excl + 1L):p_sel] <-
      stats::rnorm(n_excl, mean = 0.6 * ((-1)^(k - 1L)), sd = 0.2)
  }
  sigma_true <- rep(1.0, k_true)
  rho_true <- c(0.3, -0.3, 0.4, -0.4, 0.2, -0.2, 0.1,
                -0.1, 0.0, 0.0)[seq_len(k_true)]

  # Generate correlated (eps, u) per observation given class
  eps <- numeric(n_obs)
  u <- numeric(n_obs)
  for (i in seq_len(n_obs)) {
    k <- class_labels[i]
    rho_k <- rho_true[k]
    z_pair <- stats::rnorm(2)
    # cov(eps/sigma, u) = rho_k; simple factorization
    u[i] <- z_pair[1]
    eps[i] <- sigma_true[k] * (rho_k * z_pair[1] +
                                 sqrt(1 - rho_k^2) * z_pair[2])
  }

  mu_sel <- numeric(n_obs)
  mu_out <- numeric(n_obs)
  for (i in seq_len(n_obs)) {
    k <- class_labels[i]
    mu_sel[i] <- sum(W[i, ] * gamma_true[k, ])
    mu_out[i] <- sum(X[i, ] * beta_true[k, ])
  }
  Z <- as.integer((mu_sel + u) > 0)
  Y_full <- mu_out + eps
  Y <- ifelse(Z == 1L, Y_full, 0)   # 0 is a placeholder for unselected

  list(
    n_obs = n_obs, p_out = p_out, p_sel = p_sel, n_excl = n_excl,
    k_true = k_true,
    X = X, W = W, Z = Z, Y = Y,
    beta_true = beta_true, gamma_true = gamma_true,
    sigma_true = sigma_true, rho_true = rho_true,
    class_labels = class_labels,
    selection_rate = mean(Z)
  )
}

# --------------------------------------------------------------------
# Compile the Stan model (once). Subsequent calls reuse the compiled .so.
# --------------------------------------------------------------------
compile_dpm_heckman_model <- function(
    stan_file = file.path("benchmarks", "dpm_heckman",
                          "dpm_heckman.stan")) {
  stan_file <- normalizePath(stan_file)
  bridgestan::compile_model(stan_file)  # returns path to .so
}

# --------------------------------------------------------------------
# Build fn / gr / hess callables from a compiled Stan model and data.
#
# Stan gives us the LOG posterior density. arcopt / optim / nlminb all
# MINIMIZE, so we negate.
#
# hess_ad uses Stan's AD Hessian; hess_fd uses numDeriv::jacobian on gr.
# --------------------------------------------------------------------
make_dpm_heckman_fns <- function(data, model_so,
                                 k_trunc = 20L,
                                 sigma_beta = 5.0,
                                 sigma_gamma = 5.0,
                                 alpha = 1.0,
                                 tmpdir = tempdir()) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("package 'numDeriv' required")
  }
  data_list <- list(
    N = data$n_obs,
    P_out = data$p_out,
    P_sel = data$p_sel,
    K = as.integer(k_trunc),
    X = data$X,
    W = data$W,
    Z = as.integer(data$Z),
    Y = as.numeric(data$Y),
    alpha = as.numeric(alpha),
    sigma_beta = as.numeric(sigma_beta),
    sigma_gamma = as.numeric(sigma_gamma)
  )
  data_path <- file.path(tmpdir,
                         sprintf("dpm_heckman_data_%s.json",
                                 format(Sys.time(), "%Y%m%d%H%M%OS3")))
  writeLines(
    jsonlite::toJSON(data_list, auto_unbox = TRUE, matrix = "rowmajor"),
    data_path
  )

  model <- bridgestan::StanModel$new(lib = model_so, data = data_path,
                                     seed = 1L)

  fn <- function(x) {
    -model$log_density(x)
  }
  gr <- function(x) {
    -model$log_density_gradient(x)$gradient
  }
  hess_ad <- function(x) {
    # Stan returns the Hessian of log density; flip sign to match fn.
    -model$log_density_hessian(x)$hessian
  }
  hess_fd <- function(x) {
    numDeriv::jacobian(gr, x)
  }

  # Initialize at a constrained-zero point and unconstrain it so we have
  # a concrete starting vector that respects parameter bounds.
  k_int <- as.integer(k_trunc)
  p_out <- data$p_out
  p_sel <- data$p_sel
  # Constrained init: beta=0, gamma=0, sigma=1, rho=0, V=0.5
  init_constrained <- list(
    beta = matrix(0, k_int, p_out),
    gamma = matrix(0, k_int, p_sel),
    sigma = rep(1.0, k_int),
    rho = rep(0, k_int),
    V = rep(0.5, k_int - 1L)
  )
  init_path <- file.path(tmpdir, "dpm_heckman_init.json")
  writeLines(jsonlite::toJSON(init_constrained, auto_unbox = TRUE,
                              matrix = "rowmajor"), init_path)
  x0 <- model$param_unconstrain_json(readLines(init_path))
  # A light random perturbation away from the exact constrained zero
  set.seed(101L)
  x0 <- x0 + stats::rnorm(length(x0), sd = 0.05)

  list(
    fn = fn, gr = gr, hess_ad = hess_ad, hess_fd = hess_fd,
    n_params = model$param_unc_num(),
    model = model,
    x0 = x0,
    data_path = data_path,
    init_path = init_path
  )
}
