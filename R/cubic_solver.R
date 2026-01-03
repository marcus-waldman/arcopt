#' Solve Cubic Regularization Subproblem
#'
#' Solves the cubic subproblem: min_s m(s) = g^T s + 1/2 s^T H s + sigma/3 ||s||^3
#' using a modified Cholesky factorization approach with secular equation iteration.
#'
#' @param g Gradient vector (length n)
#' @param H Hessian matrix (n x n symmetric matrix, may be indefinite)
#' @param sigma Regularization parameter (positive scalar)
#' @param max_lambda_iter Maximum iterations for lambda secular equation solver (default: 50)
#' @param lambda_tol Convergence tolerance for secular equation (default: 1e-10)
#'
#' @return List with components:
#'   \item{s}{Solution vector (step)}
#'   \item{lambda}{Final Lagrange multiplier}
#'   \item{pred_reduction}{Predicted reduction in objective: -g^T s - 1/2 s^T H s - sigma/3 ||s||^3}
#'   \item{converged}{Logical; TRUE if secular equation converged}
#'
#' @details
#' The algorithm solves the secular equation norm(s(lambda)) = lambda/sigma where
#' s(lambda) = -(H + lambda*I)^(-1) * g. This is done using:
#' \enumerate{
#'   \item Cholesky factorization with iterative diagonal shift for indefinite H
#'   \item Newton-Raphson iteration on the secular equation phi(lambda) = ||s|| - lambda/sigma
#' }
#'
#' When H is positive definite, the algorithm starts with lambda = 0. When H is
#' indefinite, it automatically adds a diagonal shift until factorization succeeds.
#'
#' @keywords internal
solve_cubic_subproblem <- function(g, H, sigma,
                                   max_lambda_iter = 50,
                                   lambda_tol = 1e-10) {
  n <- length(g)

  # STEP 1: Try initial factorization with lambda = 0
  lambda <- 0
  eps <- 1e-12  # Small shift for indefinite case

  # Try Cholesky of H
  chol_result <- tryCatch(
    {
      chol_factor <- chol(H)  # Returns upper triangular; R's chol gives U where H = U^T U
      list(success = TRUE, chol_factor = chol_factor, lambda_lower = 0)
    },
    error = function(e) {
      # H is not positive definite, need to find minimum shift
      # Use Gershgorin theorem: smallest eigenvalue >= min_i(H[i,i] - sum_j|H[i,j]|)
      row_sums <- sapply(1:n, function(i) sum(abs(H[i, -i])))
      lambda_lower <- max(0, -min(diag(H) - row_sums))
      lambda <- lambda_lower + eps

      # Try factorization with shift
      chol_factor <- tryCatch(
        chol(H + lambda * diag(n)),
        error = function(e2) NULL
      )

      if (is.null(chol_factor)) {
        # Even with shift failed; use larger shift
        lambda <- lambda_lower + 1.0
        chol_factor <- chol(H + lambda * diag(n))
      }

      list(success = TRUE, chol_factor = chol_factor, lambda_lower = lambda_lower)
    }
  )

  chol_factor <- chol_result$chol_factor
  lambda_lower <- chol_result$lambda_lower

  # STEP 2: Solve initial system (H + lambda*I) s = -g
  # Since R's chol gives upper triangular U where H + lambda*I = U^T U:
  # U^T U s = -g
  # First solve U^T w = -g for w, then solve U s = w for s
  if (lambda == 0) {
    s <- backsolve(chol_factor, backsolve(chol_factor, -g, transpose = TRUE))
  } else {
    s <- backsolve(chol_factor, backsolve(chol_factor, -g, transpose = TRUE))
  }

  # STEP 3: Iterate secular equation ||s(lambda)|| = lambda/sigma
  converged <- FALSE
  for (iter in 1:max_lambda_iter) {
    s_norm <- sqrt(sum(s^2))

    # Secular equation residual
    phi <- s_norm - lambda / sigma

    # Check convergence
    if (abs(phi) < lambda_tol * max(1, lambda / sigma)) {
      converged <- TRUE
      break
    }

    # Newton-Raphson step on phi(lambda)
    # phi'(lambda) = -s^T (H + lambda*I)^{-1} s / ||s|| - 1/sigma
    # Let w = chol_factor^{-1} s, then s^T (H + lambda*I)^{-1} s = ||w||^2
    w <- backsolve(chol_factor, s, transpose = TRUE)  # Solve chol_factor^T w = s
    w_norm_sq <- sum(w^2)

    phi_prime <- -w_norm_sq / s_norm - 1 / sigma

    # Update lambda
    lambda_new <- lambda - phi / phi_prime
    lambda_new <- max(lambda_new, lambda_lower + eps)  # Stay in valid region

    # Recompute factorization and solution for new lambda
    chol_factor <- chol(H + lambda_new * diag(n))
    s <- backsolve(chol_factor, backsolve(chol_factor, -g, transpose = TRUE))
    lambda <- lambda_new
  }

  # STEP 4: Compute predicted reduction
  # pred_red = -g^T s - 1/2 s^T H s - sigma/3 ||s||^3
  s_norm <- sqrt(sum(s^2))
  h_s <- as.vector(H %*% s)
  pred_reduction <- -sum(g * s) - 0.5 * sum(s * h_s) - (sigma / 3) * s_norm^3

  return(list(
    s = s,
    lambda = lambda,
    pred_reduction = pred_reduction,
    converged = converged
  ))
}
