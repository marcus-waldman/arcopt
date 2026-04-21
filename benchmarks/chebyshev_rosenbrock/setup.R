# Nesterov's Chebyshev-Rosenbrock Function
# ========================================
#
# f(x) = (1/4) (x_1 - 1)^2 + sum_{i=1}^{n-1} (x_{i+1} - 2 x_i^2 + 1)^2
#
# Properties:
#   - Global minimum at x* = (1, 1, ..., 1), f* = 0
#   - Has 2^(n-1) critical points, the vast majority of which are
#     genuine local minima, not saddles
#   - Chain structure: each residual couples only x_i, x_{i+1}, so the
#     Hessian is tridiagonal (diagonal + sub/super-diagonal), very
#     sparse despite n = 500
#
# Standard benchmark in the ARC literature (Cartis, Gould, Toint 2010,
# "Adaptive cubic regularisation methods ... Part I ... numerical
# results"; Nesterov 1983 original).
#
# Starting point x_0 = -1 is adversarial: every residual starts at
# -2 (f_0 = 4n - 3 = 1997 at n = 500), and the chain must propagate
# sign corrections from left to right to reach the global minimum.

make_chebrosen_problem <- function(n) {
  stopifnot(n >= 2L)
  n <- as.integer(n)

  fn <- function(x) {
    n <- length(x)
    r_vec <- x[-1] - 2 * x[-n]^2 + 1     # length n - 1
    0.25 * (x[1] - 1)^2 + sum(r_vec * r_vec)
  }

  gr <- function(x) {
    n <- length(x)
    r_vec <- x[-1] - 2 * x[-n]^2 + 1      # length n - 1
    g <- numeric(n)
    # Contribution from (1/4)(x_1 - 1)^2
    g[1] <- 0.5 * (x[1] - 1)
    # Contribution from r_i^2 via dr_i/dx_i = -4 x_i
    g[1:(n - 1)] <- g[1:(n - 1)] + (-8) * x[1:(n - 1)] * r_vec
    # Contribution from r_i^2 via dr_i/dx_{i+1} = 1
    g[2:n] <- g[2:n] + 2 * r_vec
    g
  }

  hess <- function(x) {
    n <- length(x)
    r_vec <- x[-1] - 2 * x[-n]^2 + 1
    diag_vals <- numeric(n)
    # H[1,1] = 1/2 + (-8 r_1 + 32 x_1^2)
    diag_vals[1] <- 0.5 + (-8 * r_vec[1] + 32 * x[1]^2)
    # H[j,j] = 2 + (-8 r_j + 32 x_j^2)  for j = 2..n-1
    if (n >= 3L) {
      idx <- 2:(n - 1)
      diag_vals[idx] <- 2 + (-8 * r_vec[idx] + 32 * x[idx]^2)
    }
    # H[n,n] = 2 (only from r_{n-1}^2 term)
    diag_vals[n] <- 2
    # Off-diagonals: H[i, i+1] = -8 x_i for i = 1..n-1
    offdiag <- -8 * x[1:(n - 1)]

    h <- diag(diag_vals)
    for (i in seq_len(n - 1L)) {
      h[i, i + 1L] <- offdiag[i]
      h[i + 1L, i] <- offdiag[i]
    }
    h
  }

  list(n = n, fn = fn, gr = gr, hess = hess)
}

make_chebrosen_start <- function(n, kind = c("adversarial", "origin",
                                             "random"),
                                 seed = 1L, sd = 1.0) {
  kind <- match.arg(kind)
  if (kind == "adversarial") {
    return(rep(-1, n))
  }
  if (kind == "origin") {
    return(rep(0, n))
  }
  set.seed(seed)
  stats::rnorm(n, sd = sd)
}
