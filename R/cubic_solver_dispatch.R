#' Solve Cubic Subproblem with Automatic Solver Selection
#'
#' Unified dispatcher for cubic subproblem solvers with automatic solver
#' selection based on problem size and Hessian representation.
#'
#' @param g Gradient vector (length n)
#' @param H Hessian matrix (n x n symmetric, may be indefinite)
#' @param sigma Regularization parameter (positive scalar)
#' @param solver Solver to use: "auto" (default) or "eigen".
#'   Auto-selection uses eigendecomposition (Algorithm 5a).
#' @param ... Additional arguments passed to specific solvers
#'
#' @return List with:
#'   \item{s}{Solution vector}
#'   \item{lambda}{Lagrange multiplier}
#'   \item{pred_reduction}{Predicted reduction}
#'   \item{converged}{TRUE if solver converged}
#'   \item{solver_used}{Name of solver that was used}
#'
#' @details
#' Uses eigendecomposition-based solver (Algorithm 5a) which provides robust
#' handling of indefinite Hessians and hard cases.
#'
#' For large-scale matrix-free optimization (n > 500), see design/scalable-arcs.qmd
#' for deferred large-scale solver implementations.
#'
#' @keywords internal
solve_cubic_subproblem_dispatch <- function(g, H = NULL, sigma, solver = "auto", ...) {
  # Validate inputs
  if (is.null(H)) {
    stop("Hessian matrix 'H' must be provided")
  }

  # Check for removed LDL solver
  if (!missing(solver) && solver == "ldl") {
    stop(
      "LDL cubic solver has been removed from arcopt.\n",
      "Use solver='eigen' (or 'auto') for the eigendecomposition solver.\n",
      "The LDL approach was removed in favor of the more robust Algorithm 5a.",
      call. = FALSE
    )
  }

  if (!solver %in% c("auto", "eigen")) {
    stop("solver must be one of: 'auto', 'eigen'")
  }

  # AUTO-SELECTION LOGIC
  if (solver == "auto") {
    solver <- "eigen"
  }

  # DISPATCH TO SELECTED SOLVER
  # Currently only eigendecomposition, but architecture preserved for future solvers
  if (solver == "eigen") {
    result <- solve_cubic_eigen(g, H, sigma, ...)
  } else {
    stop("Unknown solver: ", solver)
  }

  # Add solver_used field to result
  result$solver_used <- solver

  return(result)
}
