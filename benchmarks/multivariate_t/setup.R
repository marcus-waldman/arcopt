# Multivariate-t MLE
# ==================
#
# Density (Liu & Rubin 1995):
#   f(y; mu, Sigma, nu) = Gamma((nu+p)/2) /
#                         (Gamma(nu/2) (nu pi)^{p/2} |Sigma|^{1/2}) *
#                         (1 + d(y)/nu)^{-(nu+p)/2}
# where d(y) = (y - mu)^T Sigma^{-1} (y - mu).
#
# Negative log-likelihood:
#   -ell = n [ log Gamma(nu/2) - log Gamma((nu+p)/2)
#              + (p/2) log(nu pi) + (1/2) log|Sigma| ]
#          + ((nu+p)/2) sum_i log(1 + d_i/nu)
#
# Parameterization (unconstrained vector x):
#   mu_k          location          (p params)
#   l_{kk}        log diagonal of L (p params)
#   L_{kl}, k>l   off-diagonal of L (p(p-1)/2 params)
#   eta = log nu  log degrees of freedom (1 param)
#
# where Sigma = L L^T (lower-triangular Cholesky). Total dimension:
#   n_par(p) = p + p(p+1)/2 + 1
#   p=5  -> n_par = 21      (manuscript headline)
#   p=10 -> n_par = 66      (stress)
#
# This is the §4.4 "healthy-basin" example: at the MLE the Hessian is
# strongly positive-definite, so qn_polish should be able to switch out
# of cubic mode and replace per-iteration Hessian calls with Wolfe
# line-search BFGS steps.

# ---- pack / unpack ----------------------------------------------------

mvt_npar <- function(p) {
  as.integer(p + p * (p + 1L) / 2L + 1L)
}

# Off-diagonal entries are stored in column-major order, matching how R
# fills `L[lower.tri(L)] <- off`.
mvt_pack <- function(mu, sigma, nu) {
  l_chol <- t(chol(sigma))
  diag_log <- log(diag(l_chol))
  off_diag <- l_chol[lower.tri(l_chol)]
  c(mu, diag_log, off_diag, log(nu))
}

mvt_unpack <- function(x, p) {
  mu <- x[seq_len(p)]
  diag_log <- x[p + seq_len(p)]
  n_off <- p * (p - 1L) / 2L
  off_diag <- if (n_off > 0L) x[2L * p + seq_len(n_off)] else numeric(0L)
  log_nu <- x[length(x)]

  l_chol <- matrix(0, p, p)
  diag(l_chol) <- exp(diag_log)
  if (n_off > 0L) l_chol[lower.tri(l_chol)] <- off_diag
  list(mu = mu, l_chol = l_chol, nu = exp(log_nu))
}

# ---- simulation -------------------------------------------------------

simulate_mvt <- function(p = 5L, n = 200L, nu = 5, seed = 42L) {
  set.seed(seed)
  mu_true <- rep(0, p)
  # Sigma_true: AR(1)-style correlation structure, unit marginal scale
  rho <- 0.5
  sigma_true <- rho^abs(outer(seq_len(p), seq_len(p), `-`))
  l_chol_true <- t(chol(sigma_true))

  # Sample y_i = mu + L z_i / sqrt(g_i), z_i ~ N(0, I), g_i ~ chi^2_nu / nu
  z_mat <- matrix(stats::rnorm(n * p), n, p)
  g_vec <- stats::rchisq(n, df = nu) / nu
  y_mat <- t(l_chol_true %*% t(z_mat)) / sqrt(g_vec)
  y_mat <- sweep(y_mat, 2, mu_true, FUN = "+")

  list(
    p = p,
    n = n,
    y_mat = y_mat,
    mu_true = mu_true,
    sigma_true = sigma_true,
    nu_true = nu,
    x_true = mvt_pack(mu_true, sigma_true, nu)
  )
}

# ---- objective, gradient, Hessian ------------------------------------

make_mvt_fns <- function(data) {
  y_mat <- data$y_mat
  n <- data$n
  p <- data$p

  fn <- function(x) {
    pars <- mvt_unpack(x, p)
    mu <- pars$mu
    l_chol <- pars$l_chol
    nu <- pars$nu

    yc <- sweep(y_mat, 2, mu)                          # n x p
    u_mat <- forwardsolve(l_chol, t(yc))                # p x n
    d_vec <- colSums(u_mat^2)                           # length n

    log_det_sigma <- 2 * sum(log(diag(l_chol)))

    n * (lgamma(nu / 2) - lgamma((nu + p) / 2) +
           (p / 2) * log(nu * pi) +
           0.5 * log_det_sigma) +
      0.5 * (nu + p) * sum(log1p(d_vec / nu))
  }

  gr <- function(x) {
    pars <- mvt_unpack(x, p)
    mu <- pars$mu
    l_chol <- pars$l_chol
    nu <- pars$nu

    yc <- sweep(y_mat, 2, mu)                          # n x p
    u_mat <- forwardsolve(l_chol, t(yc))                # p x n; u_i = L^{-1} (y_i - mu)
    d_vec <- colSums(u_mat^2)                           # length n
    v_mat <- backsolve(t(l_chol), u_mat)                # p x n; v_i = Sigma^{-1} (y_i - mu)

    w_vec <- (nu + p) / (nu + d_vec)                    # length n

    # d(-ell)/d mu = -sum_i w_i v_i
    g_mu <- -as.numeric(v_mat %*% w_vec)

    # d(-ell)/d L_{kl} (k >= l) = n (L^{-1})_{lk} - sum_i w_i (v_i)_k (u_i)_l
    l_inv <- forwardsolve(l_chol, diag(p))
    w_mat <- matrix(w_vec, nrow = p, ncol = n, byrow = TRUE)
    m_full <- -tcrossprod(v_mat * w_mat, u_mat) + n * t(l_inv)

    g_diag <- diag(l_chol) * diag(m_full)               # chain rule for log-diag
    g_off <- m_full[lower.tri(m_full)]

    # d(-ell)/d nu (then chain to log-nu by * nu)
    g_nu <- (n / 2) * (digamma(nu / 2) - digamma((nu + p) / 2)) +
      (n * p) / (2 * nu) +
      0.5 * sum(log1p(d_vec / nu)) -
      sum(w_vec * d_vec) / (2 * nu)
    g_log_nu <- nu * g_nu

    c(g_mu, g_diag, g_off, g_log_nu)
  }

  hess <- function(x) {
    if (!requireNamespace("numDeriv", quietly = TRUE)) {
      stop("package 'numDeriv' is required for the Hessian")
    }
    h_mat <- numDeriv::jacobian(gr, x, method = "Richardson")
    (h_mat + t(h_mat)) / 2
  }

  list(fn = fn, gr = gr, hess = hess)
}

# ---- starting point ---------------------------------------------------
#
# Method-of-moments style start that does not peek at the truth: sample
# mean for mu, sample covariance for Sigma (sample covariance over-states
# scale for heavy-tailed t but is in the right ballpark), nu = 10 (a
# generic moderate-tail default; the MLE will pull it toward truth).

mvt_start <- function(data) {
  y_mat <- data$y_mat
  mu0 <- colMeans(y_mat)
  sigma0 <- stats::cov(y_mat)
  mvt_pack(mu0, sigma0, nu = 10)
}
