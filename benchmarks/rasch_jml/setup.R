# Rasch joint maximum likelihood (JML)
# ====================================
#
# Item-response model:
#   P(Y_ij = 1) = sigma(theta_i - beta_j)
# with i = 1..N participants, j = 1..J items, sigma(z) = 1/(1+exp(-z)).
#
# Joint log-likelihood:
#   ell(theta, beta) = sum_ij [ Y_ij (theta_i - beta_j)
#                              - log(1 + exp(theta_i - beta_j)) ]
#
# Identification: ell is invariant under (theta, beta) -> (theta + c,
# beta + c). We anchor beta_1 = 0, so the free parameter vector is
#   x = (theta_1, ..., theta_N, beta_2, ..., beta_J),  n_par = N + J - 1.
# At N = 400, J = 100, n_par = 499.
#
# This is a STRONGLY CONVEX MLE (sum of negative-log-sigmoids of linear
# combinations) with a unique optimum once anchored. It is also a
# canonical large-parameter statistical problem -- N + J grows linearly
# with study size. Note: JML for Rasch is known to be inconsistent for
# beta as N -> infty with J fixed (Andersen 1973); CML or MML are the
# recommended inferential strategies. JML is used here purely as an
# *optimization* benchmark -- a real-world example of a healthy convex
# basin with many parameters and an expensive Hessian.
#
# Hessian structure (negative log-likelihood = NLL):
#   d_i = sum_j p_ij (1 - p_ij)        > 0
#   e_j = sum_i p_ij (1 - p_ij)        > 0
#   W_ij = p_ij (1 - p_ij)
#   H_theta_theta = diag(d)            (N x N)
#   H_beta_beta   = diag(e[-1])        ((J-1) x (J-1) after anchoring)
#   H_theta_beta  = -W[, -1]           (N x (J-1))
# H is positive semi-definite; after anchoring beta_1 it is positive
# definite (one null direction removed).

# ---- pack / unpack ---------------------------------------------------

rasch_npar <- function(n_persons, n_items) {
  as.integer(n_persons + n_items - 1L)
}

rasch_pack <- function(theta, beta) {
  # beta_1 is anchored to 0; pack only beta[-1].
  c(theta, beta[-1L])
}

rasch_unpack <- function(x, n_persons, n_items) {
  theta <- x[seq_len(n_persons)]
  beta_free <- x[n_persons + seq_len(n_items - 1L)]
  beta <- c(0, beta_free)
  list(theta = theta, beta = beta)
}

# ---- simulation ------------------------------------------------------

simulate_rasch <- function(n_persons = 400L, n_items = 100L, seed = 42L) {
  set.seed(seed)
  theta_true <- stats::rnorm(n_persons, mean = 0, sd = 1)
  # Center beta to match the (later) beta_1 = 0 anchor without biasing
  # the overall difficulty distribution.
  beta_raw <- stats::rnorm(n_items, mean = 0, sd = 1)
  beta_true <- beta_raw - beta_raw[1L]

  z_mat <- outer(theta_true, beta_true, `-`)
  p_mat <- stats::plogis(z_mat)
  y_mat <- matrix(stats::rbinom(n_persons * n_items, 1L, as.vector(p_mat)),
                  n_persons, n_items)

  list(
    n_persons = n_persons,
    n_items = n_items,
    y_mat = y_mat,
    theta_true = theta_true,
    beta_true = beta_true,
    x_true = rasch_pack(theta_true, beta_true)
  )
}

# ---- objective, gradient, Hessian ------------------------------------

make_rasch_fns <- function(data) {
  y_mat <- data$y_mat
  n_persons <- data$n_persons
  n_items <- data$n_items
  row_sum_y <- rowSums(y_mat)              # length N (sufficient stat for theta_i)
  col_sum_y <- colSums(y_mat)              # length J (sufficient stat for beta_j)

  fn <- function(x) {
    pars <- rasch_unpack(x, n_persons, n_items)
    theta <- pars$theta
    beta <- pars$beta

    # NLL = -sum_ij [ Y_ij z_ij - softplus(z_ij) ]
    #     = -sum_i theta_i row_sum_y_i + sum_j beta_j col_sum_y_j
    #         + sum_ij softplus(theta_i - beta_j)
    z_mat <- outer(theta, beta, `-`)            # N x J
    sum(-theta * row_sum_y) + sum(beta * col_sum_y) +
      sum(log1p(exp(-abs(z_mat))) + pmax(z_mat, 0))   # numerically stable softplus
  }

  gr <- function(x) {
    pars <- rasch_unpack(x, n_persons, n_items)
    theta <- pars$theta
    beta <- pars$beta

    z_mat <- outer(theta, beta, `-`)
    p_mat <- stats::plogis(z_mat)               # N x J

    # d(NLL)/d theta_i = -row_sum_y_i + sum_j p_ij = sum_j (p_ij - Y_ij)
    g_theta <- rowSums(p_mat) - row_sum_y       # length N
    # d(NLL)/d beta_j = col_sum_y_j - sum_i p_ij = sum_i (Y_ij - p_ij)
    g_beta_full <- col_sum_y - colSums(p_mat)   # length J
    # Drop g_beta[1] (anchored)
    c(g_theta, g_beta_full[-1L])
  }

  hess <- function(x) {
    pars <- rasch_unpack(x, n_persons, n_items)
    theta <- pars$theta
    beta <- pars$beta

    z_mat <- outer(theta, beta, `-`)
    p_mat <- stats::plogis(z_mat)
    w_mat <- p_mat * (1 - p_mat)                # N x J

    d_vec <- rowSums(w_mat)                     # length N
    e_vec <- colSums(w_mat)[-1L]                # length J-1

    n_par <- n_persons + n_items - 1L
    h_mat <- matrix(0, n_par, n_par)
    # H_theta_theta = diag(d)
    diag(h_mat)[seq_len(n_persons)] <- d_vec
    # H_beta_beta   = diag(e_free)
    beta_idx <- n_persons + seq_len(n_items - 1L)
    diag(h_mat)[beta_idx] <- e_vec
    # H_theta_beta  = -W[, -1]; symmetrize
    cross <- -w_mat[, -1L, drop = FALSE]        # N x (J-1)
    h_mat[seq_len(n_persons), beta_idx] <- cross
    h_mat[beta_idx, seq_len(n_persons)] <- t(cross)

    h_mat
  }

  list(fn = fn, gr = gr, hess = hess,
       row_sum_y = row_sum_y, col_sum_y = col_sum_y)
}

# ---- starting point --------------------------------------------------
#
# Cold start: theta_i = 0, beta_j = 0 for all. This is the natural
# uninformative initialization for a Rasch JML solver.

rasch_cold_start <- function(n_persons, n_items) {
  numeric(n_persons + n_items - 1L)
}
