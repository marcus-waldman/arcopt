#' Check Convergence Criteria (Stan-style)
#'
#' Checks multiple convergence criteria following Stan's L-BFGS implementation.
#' Termination occurs when **any** criterion is satisfied, preventing premature
#' termination on ill-scaled problems while ensuring robust detection.
#'
#' @param g_current Current gradient vector (length Q)
#' @param f_current Current objective function value (scalar)
#' @param f_previous Previous objective function value (scalar), NA for first iteration
#' @param x_current Current parameter vector (length Q)
#' @param x_previous Previous parameter vector (length Q), NULL for first iteration
#' @param iter Current iteration number (0-indexed)
#' @param tol_grad Absolute gradient tolerance
#' @param tol_rel_grad Relative gradient tolerance
#' @param tol_obj Absolute objective change tolerance
#' @param tol_rel_obj Relative objective change tolerance
#' @param tol_param Parameter change tolerance (infinity norm)
#' @param max_iter Maximum iterations allowed
#'
#' @return List with two components:
#'   \item{converged}{Logical; TRUE if any criterion satisfied}
#'   \item{reason}{Character; description of which criterion triggered, or "" if not converged}
#'
#' @details
#' The function checks 6 convergence criteria in order:
#' \enumerate{
#'   \item \strong{max_iter}: Iteration limit reached
#'   \item \strong{gradient_abs}: max(abs(g)) < tol_grad
#'   \item \strong{gradient_rel}: max(abs(g)) < tol_rel_grad * max(1, abs(f))
#'   \item \strong{objective_abs}: abs(f_current - f_previous) < tol_obj
#'   \item \strong{objective_rel}: abs(f_current - f_previous) < tol_rel_obj * max(1, abs(f_current))
#'   \item \strong{parameter}: max(abs(x_current - x_previous)) < tol_param
#' }
#'
#' @keywords internal
check_convergence <- function(g_current, f_current, f_previous,
                              x_current, x_previous, iter,
                              tol_grad = 1e-8,
                              tol_rel_grad = 1e-6,
                              tol_obj = 1e-12,
                              tol_rel_obj = 1e-8,
                              tol_param = 1e-8,
                              max_iter = 1000) {
  # 1. CHECK ITERATION LIMIT
  if (iter >= max_iter) {
    return(list(converged = TRUE, reason = "max_iter"))
  }

  # Compute gradient infinity norm (max absolute value)
  g_norm_inf <- max(abs(g_current))

  # 2. CHECK GRADIENT NORM (absolute)
  if (g_norm_inf < tol_grad) {
    return(list(converged = TRUE, reason = "gradient_abs"))
  }

  # 3. CHECK GRADIENT NORM (relative)
  f_scale <- max(1, abs(f_current))
  if (g_norm_inf < tol_rel_grad * f_scale) {
    return(list(converged = TRUE, reason = "gradient_rel"))
  }

  # Remaining checks require previous iteration (skip on first iteration)
  if (iter == 0 || is.na(f_previous) || is.null(x_previous)) {
    return(list(converged = FALSE, reason = ""))
  }

  # Compute objective change
  obj_change <- abs(f_current - f_previous)

  # For objective/parameter convergence criteria, require gradient to be
  # reasonably small (10x the strict tolerance) to avoid false positives
  # on non-stationary points
  gradient_reasonably_small <- g_norm_inf < 10 * tol_grad

  # 4. CHECK OBJECTIVE CHANGE (absolute)
  if (obj_change < tol_obj && gradient_reasonably_small) {
    return(list(converged = TRUE, reason = "objective_abs"))
  }

  # 5. CHECK OBJECTIVE CHANGE (relative)
  if (obj_change < tol_rel_obj * max(1, abs(f_current)) &&
      gradient_reasonably_small) {
    return(list(converged = TRUE, reason = "objective_rel"))
  }

  # 6. CHECK PARAMETER CHANGE (infinity norm)
  param_change <- max(abs(x_current - x_previous))
  if (param_change < tol_param && gradient_reasonably_small) {
    return(list(converged = TRUE, reason = "parameter"))
  }

  # No criteria satisfied
  return(list(converged = FALSE, reason = ""))
}
