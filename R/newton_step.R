#' Try Newton Step
#'
#' Attempts to compute and validate a Newton step when the Hessian is positive
#' definite. Uses Cholesky factorization to check positive definiteness and
#' solve the Newton system.
#'
#' @param g Gradient vector (length n)
#' @param H Hessian matrix (n x n symmetric matrix)
#'
#' @return List with components:
#'   \item{success}{Logical; TRUE if Newton step computed successfully}
#'   \item{s}{Newton step vector if successful, NULL otherwise}
#'   \item{pred_reduction}{Predicted reduction from quadratic model if successful, NULL otherwise}
#'
#' @details
#' The algorithm:
#' \enumerate{
#'   \item Attempts Cholesky factorization of H (fails if H is indefinite)
#'   \item If successful, solves H * s = -g for Newton step
#'   \item Computes predicted reduction from quadratic model
#' }
#'
#' No pre-filtering based on sigma is required - the ratio test in the main
#' ARC loop naturally determines whether Newton steps are appropriate. See
#' literature/consensus_reviews/No explicit σ-threshold for Newton Step.txt
#'
#' @keywords internal
try_newton_step <- function(g, H) {
  # STEP 1: Attempt Cholesky factorization
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
    return(list(success = FALSE, s = NULL, pred_reduction = NULL))
  }

  # STEP 2: Solve Newton system: H * s = -g
  # Since chol gives upper triangular U where H = U^T U:
  # U^T U s = -g
  chol_factor <- chol_result$chol_factor
  s_newton <- backsolve(chol_factor, backsolve(chol_factor, -g, transpose = TRUE))

  # STEP 3: Compute predicted reduction from quadratic model
  # Newton minimizes m_quad(s) = f + g^T s + (1/2) s^T H s
  # Predicted reduction = m(0) - m(s) = -(g^T s + (1/2) s^T H s)
  Hs <- as.vector(H %*% s_newton)
  pred_reduction <- -(sum(g * s_newton) + 0.5 * sum(s_newton * Hs))

  # STEP 4: Safeguard - ensure non-negative predicted reduction
  # Zero is allowed (happens when g = 0, already at stationary point)
  if (pred_reduction < 0) {
    # This shouldn't happen with PD H and s = -H^{-1}g, but check anyway
    return(list(success = FALSE, s = NULL, pred_reduction = NULL))
  }

  list(success = TRUE, s = s_newton, pred_reduction = pred_reduction)
}
