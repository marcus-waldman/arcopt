#' Apply Box Constraints to Step
#'
#' Truncates a proposed step to ensure the new iterate stays within box
#' constraints, then projects for numerical safety.
#'
#' @param x_current Current iterate
#' @param s Proposed unconstrained step
#' @param lower Lower bounds (use -Inf for unbounded)
#' @param upper Upper bounds (use Inf for unbounded)
#'
#' @return List with components:
#'   \item{s_bounded}{Truncated step that respects bounds}
#'   \item{x_new}{New iterate projected to feasible region}
#'
#' @details
#' This implements Algorithm 7a from the design specification. The approach:
#' \enumerate{
#'   \item Compute maximum feasible step length alpha_max
#'   \item Truncate step: s_bounded = alpha_max * s
#'   \item Project x_k + s_bounded to box for numerical safety
#' }
#'
#' @keywords internal
apply_box_constraints <- function(x_current, s, lower, upper) {
  n <- length(x_current)

  # STEP 1: Compute maximum feasible step length
  alpha_max <- 1.0

  for (i in seq_len(n)) {
    if (s[i] > 0 && is.finite(upper[i])) {
      # Moving toward upper bound
      alpha_max <- min(alpha_max, (upper[i] - x_current[i]) / s[i])
    } else if (s[i] < 0 && is.finite(lower[i])) {
      # Moving toward lower bound
      alpha_max <- min(alpha_max, (lower[i] - x_current[i]) / s[i])
    }
  }

  # STEP 2: Truncate step
  s_bounded <- alpha_max * s

  # STEP 3: Project for numerical safety
  x_new <- project_to_box(x_current + s_bounded, lower, upper)

  list(s_bounded = s_bounded, x_new = x_new)
}


#' Project Point to Box Constraints
#'
#' Projects a point to the feasible region defined by box constraints.
#'
#' @param x Point to project
#' @param lower Lower bounds
#' @param upper Upper bounds
#'
#' @return Projected point within bounds
#'
#' @details
#' Component-wise projection: each component is clamped to its bounds
#'
#' @keywords internal
project_to_box <- function(x, lower, upper) {
  pmax(lower, pmin(upper, x))
}


#' Validate and Project Initial Point
#'
#' Validates bounds and projects initial point to feasible region if needed.
#'
#' @param x0 User-provided initial point
#' @param lower Lower bounds
#' @param upper Upper bounds
#'
#' @return Feasible initial point
#'
#' @details
#' This implements Algorithm 7b from the design specification:
#' \enumerate{
#'   \item Validate bounds (lower <= upper, correct dimensions)
#'   \item Check if x0 is feasible
#'   \item Project to box if infeasible (with warning)
#' }
#'
#' @keywords internal
validate_and_project_initial <- function(x0, lower, upper) {
  n <- length(x0)

  # STEP 1: Validate bounds
  if (length(lower) != n || length(upper) != n) {
    stop("Bounds dimension mismatch: lower and upper must have same length as x0")
  }

  if (any(lower >= upper, na.rm = TRUE)) {
    stop("Invalid bounds: lower must be strictly less than upper for all parameters")
  }

  # STEP 2: Check feasibility
  is_feasible <- all(lower <= x0 & x0 <= upper)

  # STEP 3: Project if needed
  if (!is_feasible) {
    warning("Initial point infeasible, projecting to bounds")
    x0 <- project_to_box(x0, lower, upper)
  }

  x0
}
