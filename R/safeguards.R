#' Detect Stagnation and Recommend Action
#'
#' Monitors optimization progress and detects stagnation based on small consecutive
#' steps or function changes. Returns recommended action: continue, refresh Hessian,
#' or stop.
#'
#' @param step_norms Vector of recent step norms (most recent last)
#' @param f_values Vector of recent function values (most recent last)
#' @param tol_step Step size tolerance for stagnation (default: 1e-12)
#' @param tol_f Function change tolerance for stagnation (default: 1e-14)
#' @param max_stagnant Maximum consecutive stagnant iterations before action (default: 5)
#' @param is_qn Logical; TRUE if using quasi-Newton (enables Hessian refresh option)
#' @param already_refreshed Logical; TRUE if Hessian was already refreshed once
#'
#' @return Character; one of "continue", "refresh_hessian", or "stop"
#'
#' @details
#' Stagnation is detected when either:
#' \itemize{
#'   \item Step norms are consecutively below tol_step for max_stagnant iterations
#'   \item Function changes are consecutively below tol_f for max_stagnant iterations
#' }
#'
#' When stagnation is detected:
#' \itemize{
#'   \item If using quasi-Newton and not already refreshed: return "refresh_hessian"
#'   \item Otherwise: return "stop" (declare convergence failure)
#' }
#'
#' @keywords internal
detect_stagnation <- function(step_norms,
                              f_values,
                              tol_step = 1e-12,
                              tol_f = 1e-14,
                              max_stagnant = 5,
                              is_qn = FALSE,
                              already_refreshed = FALSE) {
  # Need at least max_stagnant values to detect stagnation
  n_steps <- length(step_norms)
  n_f <- length(f_values)

  # STEP 1: Check step size stagnation
  stagnant_steps <- 0
  if (n_steps >= max_stagnant) {
    # Count consecutive small steps from the end
    for (i in n_steps:max(1, n_steps - max_stagnant + 1)) {
      if (step_norms[i] < tol_step) {
        stagnant_steps <- stagnant_steps + 1
      } else {
        break
      }
    }
  }

  # STEP 2: Check function value stagnation
  stagnant_f <- 0
  if (n_f >= max_stagnant + 1) {  # Need n+1 values for n changes
    # Count consecutive small function changes from the end
    for (i in n_f:max(2, n_f - max_stagnant + 1)) {
      f_change <- abs(f_values[i] - f_values[i - 1])
      if (f_change < tol_f) {
        stagnant_f <- stagnant_f + 1
      } else {
        break
      }
    }
  }

  # STEP 3: Determine action
  if (stagnant_steps >= max_stagnant || stagnant_f >= max_stagnant) {
    # Stagnation detected
    if (is_qn && !already_refreshed) {
      return("refresh_hessian")
    } else {
      return("stop")
    }
  }

  "continue"
}


#' Protect Against NaN/Inf in Function Evaluations
#'
#' Detects NaN or Inf values in function and gradient evaluations and attempts
#' recovery by increasing regularization to produce smaller steps.
#'
#' @param f_value Function value to check
#' @param g_value Gradient vector to check
#'
#' @return Logical; TRUE if values are finite and valid, FALSE if NaN/Inf detected
#'
#' @details
#' This is a simple check function. The recovery logic (re-solving with larger sigma)
#' should be implemented in the main optimization loop.
#'
#' Checks that:
#' \itemize{
#'   \item Function value is finite (not NaN or Inf)
#'   \item All gradient components are finite
#' }
#'
#' @keywords internal
check_finite <- function(f_value, g_value) {
  if (!is.finite(f_value)) {
    return(FALSE)
  }

  if (any(!is.finite(g_value))) {
    return(FALSE)
  }

  TRUE
}
