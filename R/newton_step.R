#' Try Newton Step with Gershgorin Check
#'
#' Attempts to compute and validate a Newton step when the Hessian is positive
#' definite. Uses Gershgorin circle theorem for a cheap O(n^2) check before
#' attempting Cholesky factorization.
#'
#' @param g Gradient vector (length n)
#' @param H Hessian matrix (n x n symmetric matrix)
#'
#' @return List with components:
#'   \item{success}{Logical; TRUE if Newton step computed successfully}
#'   \item{s}{Newton step vector if successful, NULL otherwise}
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Computes Gershgorin lower bound on minimum eigenvalue
#'   \item If bound <= 0, returns failure (H may be indefinite)
#'   \item Attempts Cholesky factorization of H
#'   \item If successful, solves H * s = -g for Newton step
#' }
#'
#' The Gershgorin theorem states that for symmetric matrix H, the minimum
#' eigenvalue satisfies: lambda_min >= min over i of (H_ii - sum of abs(H_ij) for j != i).
#' If this lower bound is positive, H is guaranteed positive definite.
#'
#' @keywords internal
try_newton_step <- function(g, H) {
  n <- length(g)

  # STEP 1: Gershgorin pre-check (cheap O(n^2) test)
  # Compute lower bound on minimum eigenvalue
  row_sums <- sapply(1:n, function(i) sum(abs(H[i, -i])))
  gershgorin_bounds <- diag(H) - row_sums
  lambda_min_bound <- min(gershgorin_bounds)

  if (lambda_min_bound <= 0) {
    # H may be indefinite; skip Newton
    return(list(success = FALSE, s = NULL))
  }

  # STEP 2: Attempt Cholesky factorization
  chol_result <- tryCatch(
    {
      chol_factor <- chol(H)
      list(success = TRUE, chol_factor = chol_factor)
    },
    error = function(e) {
      list(success = FALSE, chol_factor = NULL)
    }
  )

  if (!chol_result$success) {
    # Cholesky failed; H is indefinite
    return(list(success = FALSE, s = NULL))
  }

  # STEP 3: Solve Newton system: H * s = -g
  # Since chol gives upper triangular U where H = U^T U:
  # U^T U s = -g
  chol_factor <- chol_result$chol_factor
  s_newton <- backsolve(chol_factor, backsolve(chol_factor, -g, transpose = TRUE))

  list(success = TRUE, s = s_newton)
}
