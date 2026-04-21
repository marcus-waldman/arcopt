# Brown Almost-Linear Function (MGH #27)
# ======================================
#
# Residuals:
#   f_i(x) = x_i + sum(x) - (n+1),  i = 1, ..., n-1
#   f_n(x) = prod(x) - 1
#
# Objective:
#   F(x) = sum_{i=1}^{n} f_i(x)^2
#
# Global minimum at x* = (1, 1, ..., 1) with F* = 0. Additional real
# critical points exist as roots of the degree-n polynomial
#   n alpha^n - (n+1) alpha^{n-1} + 1 = 0,
# which has alpha = 1 as a root plus up to n-1 other real roots.
#
# Standard start x_0 = (0.5, ..., 0.5). Classic MGH dimension is n = 10;
# the function scales formally to any n.
#
# Notes on large n: prod(x_0) = 0.5^n underflows for n > ~1000 in double
# precision, but remains meaningful for n up to a few hundred. Near the
# solution the multiplicative term sum_{j != k} 1/(x_j x_k) creates
# dense coupling that produces near-rank-deficient curvature, the
# classical failure mode this problem was designed to probe.

make_brown_problem <- function(n) {
  stopifnot(n >= 2L)
  n <- as.integer(n)

  fn <- function(x) {
    s <- sum(x)
    p <- prod(x)
    r_first <- x[-n] + s - (n + 1)  # wait: the first n-1 residuals are over all indices
    # f_i = x_i + s - (n+1)  for i = 1..n-1 means we take the first n-1 x's
    r_first <- x[seq_len(n - 1L)] + s - (n + 1)
    r_last <- p - 1
    sum(r_first * r_first) + r_last * r_last
  }

  gr <- function(x) {
    n_len <- length(x)
    s <- sum(x)
    p <- prod(x)
    f_first <- x[seq_len(n_len - 1L)] + s - (n_len + 1)
    a_sum <- sum(f_first)
    f_last <- p - 1

    # Sum over first n-1 residuals: a_sum + f_first[j] for j<=n-1, else a_sum.
    # Product residual gives 2 f_n (p/x_j).
    f_extended <- c(f_first, 0)  # zero in position n
    2 * (a_sum + f_extended) + 2 * f_last * (p / x)
  }

  hess <- function(x) {
    n_len <- length(x)
    p <- prod(x)
    f_last <- p - 1

    # Contribution from the n-1 linear residuals (all summed via
    # Sum_i ([i==j]+1)([i==k]+1)):
    # base_{jk} = (n-1) + [j<=n-1] + [k<=n-1] + [j==k<=n-1]
    indic <- c(rep(1, n_len - 1L), 0)  # length n
    # outer part: (n-1) + [j<=n-1] + [k<=n-1] as a rank-1 + constant
    base <- (n_len - 1L) + outer(indic, rep(1, n_len)) +
      outer(rep(1, n_len), indic)
    diag(base)[seq_len(n_len - 1L)] <-
      diag(base)[seq_len(n_len - 1L)] + 1

    # Contribution from f_n^2 and f_n * d2f_n/dx_j dx_k:
    # 2 (p/x_j)(p/x_k) + 2 f_n (p/(x_j x_k)) * [j != k]
    inv_x <- 1 / x
    inv_outer <- outer(inv_x, inv_x)  # 1/(x_j x_k)
    prod_contrib <- 2 * p * p * inv_outer  # 2 p^2 / (x_j x_k)
    cross <- 2 * f_last * p * inv_outer   # 2 f_n p / (x_j x_k)
    diag(cross) <- 0                       # only off-diagonal

    2 * base + prod_contrib + cross
  }

  list(n = n, fn = fn, gr = gr, hess = hess)
}

make_brown_start <- function(n, kind = c("standard", "random"),
                             seed = 1L, sd = 0.3) {
  kind <- match.arg(kind)
  if (kind == "standard") {
    return(rep(0.5, n))
  }
  set.seed(seed)
  # Random but positive to avoid sign-flip degeneracies in the product term
  abs(stats::rnorm(n, mean = 1, sd = sd))
}
