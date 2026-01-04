#' Solve Cubic Subproblem with Automatic Solver Selection
#'
#' Unified dispatcher for cubic subproblem solvers with automatic solver
#' selection based on problem size and Hessian representation.
#'
#' @param g Gradient vector (length n)
#' @param H Hessian matrix (n x n symmetric, may be indefinite). Optional if
#'   hess_vec is provided.
#' @param hess_vec Function computing Hessian-vector products: hess_vec(v) -> Hv.
#'   Optional if H is provided.
#' @param sigma Regularization parameter (positive scalar)
#' @param solver Solver to use: "auto" (default), "ldl", "eigen", or "cg"
#' @param solver_threshold Threshold for auto-selection (default: 500).
#'   Problems with n <= threshold use eigendecomposition, larger problems use CG.
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
#' Auto-selection logic:
#' - If hess_vec is provided: use "cg" (matrix-free optimization)
#' - Else if n <= solver_threshold: use "eigen" (full eigendecomposition)
#' - Else: construct hess_vec from H and use "cg" (large-scale)
#'
#' Manual solver options:
#' - "ldl": LDL factorization with modified Cholesky (legacy solver)
#' - "eigen": Eigendecomposition with explicit hard-case handling (Algorithm 5a)
#'   **Recommended for most applications**
#' - "cg": ARCqK multi-shift CG-Lanczos (Algorithm 5b)
#'   **EXPERIMENTAL - may perform poorly on small problems (n < 500)**
#'
#' @section Performance Notes:
#' The CG solver is designed for large-scale problems (n > 500) and matrix-free
#' optimization. Benchmarks show it performs poorly on small problems due to
#' discrete shift selection and Lanczos breakdown issues. For best results:
#' - Use "auto" (default) to let the dispatcher choose
#' - Use "eigen" explicitly for problems with n <= 500
#' - Use "cg" only for large-scale or matrix-free problems
#'
#' @keywords internal
solve_cubic_subproblem_dispatch <- function(g, H = NULL, hess_vec = NULL,
                                            sigma,
                                            solver = "auto",
                                            solver_threshold = 500,
                                            ...) {
  n <- length(g)

  # Validate inputs
  if (is.null(H) && is.null(hess_vec)) {
    stop("At least one of 'H' or 'hess_vec' must be provided")
  }

  if (!solver %in% c("auto", "ldl", "eigen", "cg")) {
    stop("solver must be one of: 'auto', 'ldl', 'eigen', 'cg'")
  }

  # AUTO-SELECTION LOGIC
  if (solver == "auto") {
    if (!is.null(hess_vec)) {
      # User provided hess_vec → use CG (matrix-free)
      solver_used <- "cg"
    } else if (n <= solver_threshold) {
      # Small problem → use eigendecomposition
      solver_used <- "eigen"
    } else {
      # Large problem → use CG with finite-difference hess_vec
      solver_used <- "cg"
    }
  } else {
    solver_used <- solver
  }

  # DISPATCH TO SELECTED SOLVER
  if (solver_used == "ldl") {
    # Legacy LDL solver
    if (is.null(H)) {
      stop("LDL solver requires full Hessian matrix H")
    }
    result <- solve_cubic_subproblem(g, H, sigma, ...)
  } else if (solver_used == "eigen") {
    # Eigendecomposition solver (Algorithm 5a)
    if (is.null(H)) {
      stop("Eigendecomposition solver requires full Hessian matrix H")
    }
    result <- solve_cubic_eigen(g, H, sigma, ...)
  } else if (solver_used == "cg") {
    # ARCqK CG-Lanczos solver (Algorithm 5b)
    if (is.null(hess_vec)) {
      # Construct hess_vec from H
      if (is.null(H)) {
        stop("CG solver requires either H or hess_vec")
      }
      hess_vec <- hess_vec_fd(H)
    }
    result <- solve_cubic_cg(g, hess_vec, sigma, ...)
  } else {
    stop("Unknown solver: ", solver_used)
  }

  # Add solver_used field to result
  result$solver_used <- solver_used

  return(result)
}
