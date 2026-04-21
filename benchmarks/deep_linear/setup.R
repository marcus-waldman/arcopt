# Deep Linear Network
# ===================
#
# Problem: minimize  f(W_1, W_2, W_3) = 0.5 * || Y - W_3 W_2 W_1 X ||_F^2
#
# with W_1 in R^{h1 x d_in}, W_2 in R^{h2 x h1}, W_3 in R^{d_out x h2}.
#
# Properties (Kawaguchi 2016, "Deep Learning without Poor Local Minima"):
#   - All critical points are either global minima or saddles
#   - Global minima form a continuous orbit under left/right scalings of
#     intermediate factors (multiplicative ambiguity)
#   - Saddles appear at rank-deficient products and at origin
#
# This makes it the cleanest theoretical proxy for deep neural networks
# on saddle-escape questions.
#
# Default sizing: d_in = d_out = 5, h1 = h2 = 18 -> n = 18*5 + 18*18 +
# 5*18 = 504 parameters.
#
# Parameter packing (column-major for each layer):
#   x[1..(h1*d_in)]                       -> W_1 (h1 x d_in)
#   x[(h1*d_in + 1)..(h1*d_in + h2*h1)]   -> W_2 (h2 x h1)
#   x[...end]                              -> W_3 (d_out x h2)

make_deep_linear_problem <- function(d_in = 5L, h1 = 18L, h2 = 18L,
                                     d_out = 5L,
                                     n_samples = 20L,
                                     w_star_seed = 1L,
                                     data_seed = 2L) {
  dims <- list(
    d_in = d_in, h1 = h1, h2 = h2, d_out = d_out,
    n_w1 = h1 * d_in,
    n_w2 = h2 * h1,
    n_w3 = d_out * h2,
    n_samples = n_samples
  )
  dims$n <- dims$n_w1 + dims$n_w2 + dims$n_w3

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }

  set.seed(w_star_seed)
  w_star <- matrix(stats::rnorm(d_out * d_in), d_out, d_in)

  set.seed(data_seed)
  x_data <- matrix(stats::rnorm(d_in * n_samples), d_in, n_samples)
  y_data <- w_star %*% x_data

  if (!is.null(old_seed)) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }

  unpack <- function(x) {
    w1 <- matrix(x[seq_len(dims$n_w1)], h1, d_in)
    w2 <- matrix(x[dims$n_w1 + seq_len(dims$n_w2)], h2, h1)
    w3 <- matrix(x[dims$n_w1 + dims$n_w2 + seq_len(dims$n_w3)], d_out, h2)
    list(w1 = w1, w2 = w2, w3 = w3)
  }

  fn <- function(x) {
    w <- unpack(x)
    pred <- w$w3 %*% w$w2 %*% w$w1 %*% x_data
    resid_mat <- y_data - pred
    0.5 * sum(resid_mat * resid_mat)
  }

  gr <- function(x) {
    w <- unpack(x)
    b2 <- w$w1 %*% x_data             # h1 x N
    b3 <- w$w2 %*% b2                 # h2 x N
    pred <- w$w3 %*% b3               # d_out x N
    resid_mat <- y_data - pred        # d_out x N

    # dL / dW_3 = -R B_3'
    g_w3 <- -resid_mat %*% t(b3)

    # dL / dW_2 = -W_3' R B_2'
    g_w2 <- -t(w$w3) %*% resid_mat %*% t(b2)

    # dL / dW_1 = -W_2' W_3' R X'
    g_w1 <- -t(w$w2) %*% t(w$w3) %*% resid_mat %*% t(x_data)

    c(as.vector(g_w1), as.vector(g_w2), as.vector(g_w3))
  }

  list(
    dims = dims,
    w_star = w_star,
    x_data = x_data,
    y_data = y_data,
    unpack = unpack,
    fn = fn,
    gr = gr
  )
}

make_deep_linear_start <- function(dims, seed, sd = 0.3) {
  set.seed(seed)
  stats::rnorm(dims$n, sd = sd)
}
