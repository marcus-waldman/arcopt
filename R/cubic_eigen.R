#' Solve Cubic Subproblem via Eigendecomposition (Algorithm 5a)
#'
#' Solves min_s m(s) = g^T s + 1/2 s^T H s + sigma/3 ||s||^3 using full
#' eigendecomposition with explicit hard-case detection.
#'
#' @param g Gradient vector (length n)
#' @param H Hessian matrix (n x n symmetric, may be indefinite)
#' @param sigma Regularization parameter (positive scalar)
#' @param max_lambda_iter Maximum secular equation iterations (default: 50)
#' @param lambda_tol Secular equation convergence tolerance (default: 1e-10)
#' @param hard_case_tol Hard case detection threshold (default: 1e-8)
#'
#' @return List with:
#'   \item{s}{Solution vector}
#'   \item{lambda}{Lagrange multiplier}
#'   \item{pred_reduction}{Predicted reduction: -g^T s - 1/2 s^T H s - sigma/3 ||s||^3}
#'   \item{converged}{TRUE if secular equation converged}
#'   \item{case_type}{"easy", "hard", or "zero_gradient"}
#'   \item{eigenvalues}{Eigenvalues of `H` sorted increasing, or `NULL` for
#'     the zero-gradient early-return path. Exposed so callers can reuse
#'     the spectrum (e.g. for `lambda_min` in end-of-iteration detectors)
#'     instead of re-decomposing `H`.}
#'
#' @details
#' Uses eigen() for O(n^3) eigendecomposition. Hard case occurs when gradient
#' is nearly orthogonal to smallest eigenvector (|g_1| < tol_hard * ||g||).
#' Recommended for n <= 500 due to cubic computational cost.
#'
#' @keywords internal
solve_cubic_eigen <- function(g, H, sigma,
                              max_lambda_iter = 50,
                              lambda_tol = 1e-10,
                              hard_case_tol = 1e-8) {
  n <- length(g)

  # STEP 0: HANDLE ZERO GRADIENT
  g_norm <- sqrt(sum(g^2))
  if (g_norm < .Machine$double.eps) {
    return(list(
      s = rep(0, n),
      lambda = 0,
      pred_reduction = 0,
      converged = TRUE,
      case_type = "zero_gradient",
      eigenvalues = NULL
    ))
  }

  # STEP 1: EIGENDECOMPOSITION
  # eigen() returns eigenvalues in decreasing order, so reverse to get increasing
  eig <- eigen(H, symmetric = TRUE)
  eigenvalues <- rev(eig$values)  # Now lambda_1 <= ... <= lambda_n
  eigenvectors <- eig$vectors[, rev(seq_len(n))]  # Reverse columns to match

  # Transform gradient to spectral coordinates
  g_tilde <- as.vector(t(eigenvectors) %*% g)

  # STEP 2: DETERMINE LOWER BOUND
  lambda_lower <- max(0, -eigenvalues[1])

  # STEP 3: CHECK FOR HARD CASE
  is_hard_case <- (eigenvalues[1] < 0) &&
                  (abs(g_tilde[1]) < hard_case_tol * g_norm)

  if (is_hard_case) {
    # HARD CASE: g nearly orthogonal to smallest eigenvector
    lambda_star <- -eigenvalues[1]
    s_tilde <- numeric(n)

    # Compute s_tilde[i] = -g_tilde[i]/(lambda_i + lambda_star) for i >= 2
    for (i in 2:n) {
      denom <- eigenvalues[i] + lambda_star
      if (abs(denom) > .Machine$double.eps) {
        s_tilde[i] <- -g_tilde[i] / denom
      } else {
        s_tilde[i] <- 0
      }
    }

    # Add component along v_1 to satisfy secular equation
    s_partial_norm <- sqrt(sum(s_tilde[2:n]^2))
    target_norm <- lambda_star / sigma

    if (target_norm^2 > s_partial_norm^2) {
      tau <- sqrt(target_norm^2 - s_partial_norm^2)
      s_tilde[1] <- tau  # Direction doesn't matter
    } else {
      s_tilde[1] <- 0
    }

    # Transform back to original coordinates
    s <- as.vector(eigenvectors %*% s_tilde)

    # Compute predicted reduction
    h_s <- as.vector(H %*% s)
    s_norm <- sqrt(sum(s^2))
    pred_reduction <- -sum(g * s) - 0.5 * sum(s * h_s) - (sigma / 3) * s_norm^3

    return(list(
      s = s,
      lambda = lambda_star,
      pred_reduction = pred_reduction,
      converged = TRUE,
      case_type = "hard",
      eigenvalues = eigenvalues
    ))
  }

  # STEP 4: EASY CASE - SECULAR EQUATION ITERATION
  # Initialize lambda > lambda_lower to ensure positive definiteness
  lambda <- lambda_lower + 0.1 * max(1, abs(lambda_lower))
  converged <- FALSE
  s_tilde <- numeric(n)

  for (iter in seq_len(max_lambda_iter)) {
    # Compute step in spectral coordinates: s_tilde[i] = -g_tilde[i]/(lambda_i + lambda)
    for (i in seq_len(n)) {
      denom <- eigenvalues[i] + lambda
      if (abs(denom) > .Machine$double.eps) {
        s_tilde[i] <- -g_tilde[i] / denom
      } else {
        s_tilde[i] <- 0
      }
    }

    s_norm <- sqrt(sum(s_tilde^2))

    # Check for numerical issues early
    if (!is.finite(s_norm) || s_norm < .Machine$double.eps) {
      break
    }

    # Secular equation residual
    phi <- s_norm - lambda / sigma

    # Check convergence
    if (!is.finite(phi)) {
      break
    }

    if (abs(phi) < lambda_tol * max(1, lambda / sigma)) {
      converged <- TRUE
      break
    }

    # Newton step on phi(lambda)
    # phi'(lambda) = -sum_i s_tilde[i]^2 / (lambda_i + lambda) / s_norm - 1/sigma
    sum_term <- 0
    for (i in seq_len(n)) {
      sum_term <- sum_term + s_tilde[i]^2 / (eigenvalues[i] + lambda)
    }
    phi_prime <- -sum_term / s_norm - 1 / sigma

    # Check for numerical issues
    if (!is.finite(phi_prime)) {
      break
    }

    lambda_new <- lambda - phi / phi_prime
    lambda <- max(lambda_new, lambda_lower + .Machine$double.eps)
  }

  # Transform back to original coordinates
  s <- as.vector(eigenvectors %*% s_tilde)

  # STEP 5: COMPUTE PREDICTED REDUCTION
  h_s <- as.vector(H %*% s)
  s_norm <- sqrt(sum(s^2))
  pred_reduction <- -sum(g * s) - 0.5 * sum(s * h_s) - (sigma / 3) * s_norm^3

  return(list(
    s = s,
    lambda = lambda,
    pred_reduction = pred_reduction,
    converged = converged,
    case_type = "easy",
    eigenvalues = eigenvalues
  ))
}
