# 3PL IRT Model, Marginal Maximum Likelihood
# ==========================================
#
# Item response function:
#   P_i(theta) = c_i + (1 - c_i) / (1 + exp(-a_i (theta - b_i)))
#
# Marginal likelihood (ability theta_j ~ N(0, 1) integrated out via
# Gauss-Hermite quadrature):
#   L(a, b, c) = prod_{j=1..N} [ sum_q w_q prod_i p_iq^{y_ij} (1-p_iq)^{1-y_ij} ]
#
# We minimize the negative log-likelihood. Parameters are stored
# unconstrained via:
#   alpha_i = log(a_i)        (a > 0)
#   beta_i  = b_i             (unconstrained)
#   gamma_i = logit(c_i)      (c in (0, 1))
# Ordering of the parameter vector x is (alpha_1..J, beta_1..J,
# gamma_1..J) -- grouped by parameter type for cleaner vectorization.
#
# 3PL without a prior on c is a classical pathology: the c parameters
# are weakly identified when items are easy (most theta_j above b_i so
# the lower asymptote is unsampled). We test "hard items" (b_i large)
# first as a should-work baseline, then will switch to "easy items" to
# probe the pathology.

sigmoid <- function(x) 1 / (1 + exp(-x))

irt3pl_pack <- function(a, b, c) {
  c(log(a), b, stats::qlogis(c))
}

irt3pl_unpack <- function(x, j_items) {
  list(
    a = exp(x[seq_len(j_items)]),
    b = x[j_items + seq_len(j_items)],
    c = stats::plogis(x[2L * j_items + seq_len(j_items)])
  )
}

# Simulate a 3PL data set.
#   items: "hard" (b in [0.5, 2.5]) or "easy" (b in [-2.5, -0.5]) or
#          "mixed" (b in [-1.5, 1.5])
simulate_irt3pl <- function(j_items = 10L, n_examinees = 1000L,
                            items = c("hard", "easy", "mixed"),
                            seed = 42L) {
  items <- match.arg(items)
  set.seed(seed)
  a_true <- stats::runif(j_items, 0.8, 2.0)
  b_true <- switch(items,
                   hard  = stats::runif(j_items,  0.5,  2.5),
                   easy  = stats::runif(j_items, -2.5, -0.5),
                   mixed = stats::runif(j_items, -1.5,  1.5))
  c_true <- stats::runif(j_items, 0.1, 0.25)
  theta <- stats::rnorm(n_examinees)
  # P_ij has shape (n_examinees, j_items)
  eta <- outer(theta, b_true, `-`) * matrix(a_true,
                                            nrow = n_examinees,
                                            ncol = j_items,
                                            byrow = TRUE)
  p_mat <- matrix(c_true, nrow = n_examinees, ncol = j_items,
                  byrow = TRUE) +
    (1 - matrix(c_true, nrow = n_examinees, ncol = j_items,
                byrow = TRUE)) * sigmoid(eta)
  y <- matrix(stats::rbinom(n_examinees * j_items, 1, as.vector(p_mat)),
              n_examinees, j_items)

  list(
    j_items = j_items,
    n_examinees = n_examinees,
    items = items,
    a_true = a_true,
    b_true = b_true,
    c_true = c_true,
    theta_true = theta,
    y = y,
    x_true = irt3pl_pack(a_true, b_true, c_true)
  )
}

# Build negative log-likelihood fn and gradient gr on the unconstrained
# parameter vector x = (alpha, beta, gamma) via Gauss-Hermite.
#
# Numerical stability: sigmoid(eta) exactly equals 1 for eta >~ 37 in
# double precision, which makes log(1 - p) = -Inf and dp/(1-p) = NaN.
# We therefore compute log p_iq and log(1 - p_iq) directly in log-space:
#   log p         = log(c + (1-c) sigmoid(eta))
#                 = LSE(log c, log(1-c) + log_sig(eta))
#   log (1 - p)   = log(1-c) + log(1 - sigmoid(eta))
#                 = log(1-c) + log_sig(-eta)
# where log_sig(eta) = plogis(eta, log.p = TRUE) = -softplus(-eta).
#
# For the gradient, we use analytical simplifications that remove the
# removable singularities:
#   dp/da / (1 - p) = s (theta - b)       (the (1-c)(1-s) factors cancel)
#   dp/db / (1 - p) = -s a
#   dp/dc / (1 - p) = 1 / (1 - c)
#   dp/da / p       = (1-c) s (1-s) (theta - b) / p
#   dp/db / p       = -(1-c) s (1-s) a / p
#   dp/dc / p       = (1 - s) / p
# All of these are finite even when sigmoid saturates.
make_irt3pl_fns <- function(data, q_nodes = 40L) {
  if (!requireNamespace("statmod", quietly = TRUE)) {
    stop("package 'statmod' required for Gauss-Hermite quadrature")
  }
  gh <- statmod::gauss.quad.prob(q_nodes, dist = "normal")
  nodes <- gh$nodes
  log_w <- log(gh$weights)
  j_items <- data$j_items
  n_examinees <- data$n_examinees
  y <- data$y
  y_t <- t(y)                                # J x N

  logspace_add <- function(a, b) {
    # element-wise log(exp(a) + exp(b)), vectorized on matrices of same
    # shape
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
  }

  compute_all <- function(a, b, cc) {
    eta <- outer(nodes, b, `-`) *
      matrix(a, nrow = q_nodes, ncol = j_items, byrow = TRUE)
    log_sig <- stats::plogis(eta, log.p = TRUE)            # Q x J
    log_1msig <- stats::plogis(-eta, log.p = TRUE)         # Q x J
    sig <- exp(log_sig)                                    # Q x J
    c_mat <- matrix(cc, nrow = q_nodes, ncol = j_items, byrow = TRUE)
    log_c <- matrix(log(cc), nrow = q_nodes, ncol = j_items,
                    byrow = TRUE)
    log_1mc <- matrix(log(1 - cc), nrow = q_nodes, ncol = j_items,
                      byrow = TRUE)
    # log p = LSE(log c, log(1-c) + log_sig)
    log_p <- logspace_add(log_c, log_1mc + log_sig)
    log_1mp <- log_1mc + log_1msig
    list(eta = eta, sig = sig, c_mat = c_mat,
         log_p = log_p, log_1mp = log_1mp)
  }

  compute_ll_q_n <- function(log_p, log_1mp) {
    log_p %*% y_t + log_1mp %*% (1 - y_t)                  # Q x N
  }

  fn <- function(x) {
    pars <- irt3pl_unpack(x, j_items)
    out <- compute_all(pars$a, pars$b, pars$c)
    ll_q_n <- compute_ll_q_n(out$log_p, out$log_1mp)
    m <- ll_q_n + log_w
    mx <- apply(m, 2, max)
    log_marg <- mx + log(colSums(
      exp(m - matrix(mx, q_nodes, n_examinees, byrow = TRUE))
    ))
    -sum(log_marg)
  }

  gr <- function(x) {
    pars <- irt3pl_unpack(x, j_items)
    a <- pars$a; b <- pars$b; cc <- pars$c
    out <- compute_all(a, b, cc)
    sig <- out$sig
    log_p <- out$log_p
    c_mat <- out$c_mat

    ll_q_n <- compute_ll_q_n(log_p, out$log_1mp)

    # Posterior weights pi_qj over Q nodes
    m <- ll_q_n + log_w
    mx <- apply(m, 2, max)
    log_denom <- mx + log(colSums(
      exp(m - matrix(mx, q_nodes, n_examinees, byrow = TRUE))
    ))
    pi_q_n <- exp(m - matrix(log_denom, q_nodes, n_examinees,
                             byrow = TRUE))                # Q x N

    node_minus_b <- outer(nodes, b, `-`)                    # Q x J
    a_mat <- matrix(a, nrow = q_nodes, ncol = j_items, byrow = TRUE)
    s1ms <- sig * (1 - sig)                                 # Q x J
    one_minus_c <- 1 - c_mat                                 # Q x J
    p_mat <- exp(log_p)                                      # Q x J

    # Stable score matrices (all Q x J)
    sc_a_pos <- one_minus_c * s1ms * node_minus_b / p_mat
    sc_a_neg <- sig * node_minus_b                  # analytically simplified
    sc_b_pos <- -one_minus_c * s1ms * a_mat / p_mat
    sc_b_neg <- -sig * a_mat
    sc_c_pos <- (1 - sig) / p_mat
    sc_c_neg <- 1 / one_minus_c

    a_pos_n <- crossprod(sc_a_pos, pi_q_n)
    a_neg_n <- crossprod(sc_a_neg, pi_q_n)
    b_pos_n <- crossprod(sc_b_pos, pi_q_n)
    b_neg_n <- crossprod(sc_b_neg, pi_q_n)
    c_pos_n <- crossprod(sc_c_pos, pi_q_n)
    c_neg_n <- crossprod(sc_c_neg, pi_q_n)

    g_a <- rowSums(y_t * a_pos_n - (1 - y_t) * a_neg_n)
    g_b <- rowSums(y_t * b_pos_n - (1 - y_t) * b_neg_n)
    g_c <- rowSums(y_t * c_pos_n - (1 - y_t) * c_neg_n)

    # Chain rule to unconstrained params; negate for neg log-lik.
    -c(g_a * a, g_b, g_c * cc * (1 - cc))
  }

  list(fn = fn, gr = gr, nodes = nodes, log_w = log_w)
}
