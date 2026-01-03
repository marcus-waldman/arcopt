# Finite Difference Utilities for Derivative Validation
# ======================================================
#
# Helper functions to verify analytic derivatives against numerical
# approximations using central finite differences.

#' Check Gradient via Finite Differences
#'
#' Compares an analytic gradient against a numerical gradient computed
#' via central finite differences.
#'
#' @param fn Function that computes f(x)
#' @param gr Function that computes gradient ∇f(x)
#' @param x Point at which to check derivatives (numeric vector)
#' @param tol Tolerance for maximum absolute difference (default: 1e-5)
#' @param h Step size for finite differences (default: 1e-8)
#'
#' @return List with:
#'   - passed: Logical, TRUE if max difference < tol
#'   - max_diff: Maximum absolute difference between analytic and numeric
#'   - g_analytic: Analytic gradient at x
#'   - g_numeric: Numerical gradient at x
check_gradient <- function(fn, gr, x, tol = 1e-5, h = 1e-8) {
  # Compute analytic gradient
  g_analytic <- gr(x)

  # Compute numerical gradient via central differences
  # g[i] ≈ (f(x + h*e_i) - f(x - h*e_i)) / (2h)
  n <- length(x)
  g_numeric <- numeric(n)

  for (i in 1:n) {
    x_plus <- x
    x_minus <- x
    x_plus[i] <- x[i] + h
    x_minus[i] <- x[i] - h
    g_numeric[i] <- (fn(x_plus) - fn(x_minus)) / (2 * h)
  }

  # Check maximum absolute difference
  max_diff <- max(abs(g_analytic - g_numeric))

  list(
    passed = max_diff < tol,
    max_diff = max_diff,
    g_analytic = g_analytic,
    g_numeric = g_numeric
  )
}

#' Check Hessian via Finite Differences
#'
#' Compares an analytic Hessian against a numerical Hessian computed
#' via central finite differences of the gradient.
#'
#' @param gr Function that computes gradient ∇f(x)
#' @param hess Function that computes Hessian H(x)
#' @param x Point at which to check derivatives (numeric vector)
#' @param tol Tolerance for maximum absolute difference (default: 1e-4)
#' @param h Step size for finite differences (default: 1e-6)
#'
#' @return List with:
#'   - passed: Logical, TRUE if max difference < tol
#'   - max_diff: Maximum absolute difference between analytic and numeric
#'   - H_analytic: Analytic Hessian at x
#'   - H_numeric: Numerical Hessian at x
check_hessian <- function(gr, hess, x, tol = 1e-4, h = 1e-6) {
  # Compute analytic Hessian
  H_analytic <- hess(x)

  # Compute numerical Hessian via central differences of gradient
  # H[i,j] ≈ (∇f_j(x + h*e_i) - ∇f_j(x - h*e_i)) / (2h)
  n <- length(x)
  H_numeric <- matrix(0, n, n)

  for (i in 1:n) {
    x_plus <- x
    x_minus <- x
    x_plus[i] <- x[i] + h
    x_minus[i] <- x[i] - h
    H_numeric[, i] <- (gr(x_plus) - gr(x_minus)) / (2 * h)
  }

  # Symmetrize (should already be symmetric in theory, but enforce numerically)
  H_numeric <- (H_numeric + t(H_numeric)) / 2

  # Check maximum absolute difference
  max_diff <- max(abs(H_analytic - H_numeric))

  list(
    passed = max_diff < tol,
    max_diff = max_diff,
    H_analytic = H_analytic,
    H_numeric = H_numeric
  )
}
