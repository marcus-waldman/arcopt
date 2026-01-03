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
#'   \item Modified Cholesky (fastmatrix::mchol) to determine minimum shift lambda_lower
#'   \item Newton-Raphson iteration on the secular equation phi(lambda) = ||s|| - lambda/sigma
#' }
#'
#' Uses fastmatrix::ldl() to compute the LDL decomposition. The diagonal D contains
#' eigenvalue information - negative diagonal elements indicate indefiniteness. The
#' minimum shift lambda_lower = max(0, -min(D)) ensures positive definiteness.
#'
#' @keywords internal
solve_cubic_subproblem <- function(g, H, sigma,
                                   max_lambda_iter = 50,
                                   lambda_tol = 1e-10) {
  n <- length(g)

  # Handle zero gradient case: optimal step is s = 0
  g_norm <- sqrt(sum(g^2))
  if (g_norm < .Machine$double.eps) {
    return(list(
      s = rep(0, n),
      lambda = 0,
      pred_reduction = 0,
      converged = TRUE
    ))
  }

  # STEP 1: Modified Cholesky via LDL to get perturbation
  ldl_result <- fastmatrix::ldl(H)
  d_diag <- ldl_result$d

  # Minimum eigenvalue approximation from LDL diagonal
  # If any diagonal element is negative, that's the minimum eigenvalue
  lambda_lower <- max(0, -min(d_diag))

  # Use sqrt(machine epsilon) as robust buffer

  # sqrt(.Machine$double.eps) ≈ 1.49e-8 - standard numerical tolerance
  # This ensures H + lambda*I is strictly positive definite (not just semi-definite)
  sqrt_eps <- sqrt(.Machine$double.eps)
  eps <- max(sqrt_eps, sqrt_eps * lambda_lower)
  lambda <- lambda_lower + eps

  # Get Cholesky factor with retry on failure
  chol_factor <- NULL
  for (attempt in 1:5) {
    chol_factor <- tryCatch(
      chol(H + lambda * diag(n)),
      error = function(e) NULL
    )
    if (!is.null(chol_factor)) break
    lambda <- lambda * 2  # Double lambda and retry
  }
  if (is.null(chol_factor)) {
    stop("Failed to compute Cholesky factorization of regularized Hessian")
  }

  # STEP 2: Solve initial system (H + lambda*I) s = -g
  # Since R's chol gives upper triangular U where H + lambda*I = U^T U:
  # U^T U s = -g
  # First solve U^T w = -g for w, then solve U s = w for s
  s <- backsolve(chol_factor, backsolve(chol_factor, -g, transpose = TRUE))

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
    lambda_new <- max(lambda_new, lambda_lower)  # Stay in valid region

    # Recompute factorization and solution for new lambda with retry on failure
    chol_success <- FALSE
    for (chol_attempt in 1:5) {
      chol_result <- tryCatch(
        chol(H + lambda_new * diag(n)),
        error = function(e) NULL
      )
      if (!is.null(chol_result)) {
        chol_factor <- chol_result
        chol_success <- TRUE
        break
      }
      lambda_new <- lambda_new * 2  # Double lambda and retry
    }
    if (!chol_success) {
      # If all retries fail, stop iterating and use current solution
      break
    }
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
