# Woodbury-Based Cubic Solver (Low-Rank QN)
# ==========================================
#
# Implements Algorithm 5b: Efficient cubic subproblem solver exploiting
# low-rank QN structure. Reduces complexity from O(n^3) to O(m^2 n + m^3).


#' Woodbury-Based Cubic Subproblem Solver
#'
#' Solves the cubic regularization subproblem when the Hessian approximation
#' has low-rank structure from quasi-Newton updates: B = gamma * I + w' c_mat w.
#'
#' @param g Gradient vector (length n).
#' @param gamma Positive scalar for initial Hessian scaling (B0 = gamma * I).
#' @param w Low-rank update matrix (m x n, where m << n). If NULL, uses
#'   eigendecomposition fallback.
#' @param c_mat Small coefficient matrix (m x m, symmetric). May be indefinite
#'   for SR1 updates.
#' @param sigma Positive cubic regularization parameter.
#' @param tol Convergence tolerance for secular equation (default: 1e-10).
#' @param max_iter Maximum Newton/bisection iterations (default: 50).
#'
#' @return A list with components:
#'   \item{s}{Step vector (length n).}
#'   \item{lambda}{Regularization parameter satisfying secular equation.}
#'   \item{iterations}{Number of iterations used.}
#'   \item{converged}{Logical indicating convergence.}
#'
#' @details
#' For limited-memory quasi-Newton methods (L-BFGS, L-SR1), the Hessian
#' approximation has the form B = gamma * I + w' c_mat w where w contains the
#' secant vectors and c_mat is a small m x m matrix.
#'
#' This solver uses the Woodbury identity to compute:
#' (B + lambda * I)^-1 g = (1/(gamma + lambda)) g
#'                        - (1/(gamma + lambda)^2) w' M^-1 w g
#' where M = c_mat^-1 + (1/(gamma + lambda)) w w'.
#'
#' Complexity is O(m^2 n + m^3 log(1/tol)) vs O(n^3) for eigendecomposition.
#'
#' @references
#' Algorithm 5b in design/pseudocode.qmd
#'
#' @keywords internal
solve_cubic_woodbury <- function(g, gamma, w, c_mat, sigma,
                                 tol = 1e-10, max_iter = 50L) {
  n <- length(g)
  g_norm <- sqrt(sum(g^2))

  # Step 0: Handle zero gradient
  if (g_norm < tol) {
    return(list(
      s = rep(0, n),
      lambda = 0,
      iterations = 0L,
      converged = TRUE
    ))
  }

  # Validate inputs
  if (gamma <= 0) {
    stop("gamma must be positive")
  }
  if (sigma <= 0) {
    stop("sigma must be positive")
  }

  # Handle case when w is NULL or empty (no QN updates yet)
  if (is.null(w) || nrow(w) == 0) {
    # B = gamma * I, use simple closed-form
    return(solve_cubic_simple(g, gamma, sigma, tol, max_iter))
  }

  m <- nrow(w)

  # Step 1: Precompute Woodbury components
  w_g <- as.vector(w %*% g)  # m x 1
  w_wt <- w %*% t(w)         # m x m (cache for reuse)

  # Precompute c_mat^{-1} with regularization for numerical stability
  c_inv <- tryCatch(
    solve(c_mat),
    error = function(e) {
      # c_mat may be singular or ill-conditioned (especially for SR1)
      # Use pseudoinverse with regularization
      eig <- eigen(c_mat, symmetric = TRUE)
      tol_eig <- max(abs(eig$values)) * .Machine$double.eps * m
      valid <- abs(eig$values) > tol_eig
      if (sum(valid) == 0) {
        # c_mat is essentially zero, treat as B = gamma * I
        return(NULL)
      }
      eig$vectors[, valid, drop = FALSE] %*%
        diag(1 / eig$values[valid], nrow = sum(valid)) %*%
        t(eig$vectors[, valid, drop = FALSE])
    }
  )

  # If c_mat^{-1} computation failed completely, fall back to simple solver
  if (is.null(c_inv)) {
    return(solve_cubic_simple(g, gamma, sigma, tol, max_iter))
  }

  # Step 4: Determine search interval
  lambda_init <- sigma * sqrt(g_norm)  # Initial guess
  lambda_lo <- 0
  lambda_hi <- max(lambda_init, 10 * sigma * g_norm)  # Conservative upper bound

  # Step 5: Newton/bisection iteration for secular equation
  for (iter in seq_len(max_iter)) {
    # Compute step at current lambda
    step_result <- woodbury_step(lambda_init, g, gamma, w_g, w_wt, c_inv, w)
    s <- step_result$s
    norm_s <- sqrt(sum(s^2))

    # Secular equation value: psi = ||s|| - lambda/sigma
    psi <- norm_s - lambda_init / sigma

    # Check convergence
    if (abs(psi) < tol * (1 + norm_s)) {
      return(list(
        s = s,
        lambda = lambda_init,
        iterations = iter,
        converged = TRUE
      ))
    }

    # Update bounds based on psi sign
    if (psi > 0) {
      # ||s|| too large, need more lambda
      lambda_lo <- lambda_init
    } else {
      # ||s|| too small, need less lambda
      lambda_hi <- lambda_init
    }

    # Newton step with safeguards
    # Derivative: psi'(lambda) = -s' (B + lambda I)^{-2} g / ||s|| - 1/sigma
    s_deriv <- woodbury_step_deriv(lambda_init, s, gamma, w_wt, c_inv, w)
    psi_prime <- -sum(s * s_deriv) / norm_s - 1 / sigma

    lambda_newton <- lambda_init - psi / psi_prime
    lambda_bisect <- (lambda_lo + lambda_hi) / 2

    # Choose Newton if in valid interval, else bisection
    if (lambda_lo < lambda_newton && lambda_newton < lambda_hi) {
      lambda_init <- lambda_newton
    } else {
      lambda_init <- lambda_bisect
    }
  }

  # Max iterations reached - return best approximation
  s_final <- woodbury_step(lambda_init, g, gamma, w_g, w_wt, c_inv, w)$s

  warning("Woodbury solver: max iterations reached")
  list(
    s = s_final,
    lambda = lambda_init,
    iterations = max_iter,
    converged = FALSE
  )
}


#' Compute Step Using Woodbury Identity
#'
#' Computes s(lambda) = -(B + lambda * I) inverse times g using the Woodbury formula.
#'
#' @param lambda Current regularization parameter.
#' @param g Gradient vector.
#' @param gamma Scaling parameter (B0 = gamma * I).
#' @param w_g Precomputed w times g.
#' @param w_wt Precomputed w times t(w).
#' @param c_inv Precomputed inverse of c_mat.
#' @param w Low-rank update matrix.
#'
#' @return List with s (step vector) and m_mat (for potential reuse).
#'
#' @keywords internal
woodbury_step <- function(lambda, g, gamma, w_g, w_wt, c_inv, w) {
  alpha <- gamma + lambda

  # M = c_mat^{-1} + (1/alpha) w w'
  m_mat <- c_inv + w_wt / alpha

  # Solve M z = w g
  z <- tryCatch(
    solve(m_mat, w_g),
    error = function(e) {
      # Fallback: use pseudoinverse
      eig <- eigen(m_mat, symmetric = TRUE)
      tol_eig <- max(abs(eig$values)) * .Machine$double.eps * length(w_g)
      valid <- abs(eig$values) > tol_eig
      if (sum(valid) == 0) {
        return(rep(0, length(w_g)))
      }
      eig$vectors[, valid, drop = FALSE] %*%
        (t(eig$vectors[, valid, drop = FALSE]) %*% w_g /
           eig$values[valid])
    }
  )

  # s = (-1/alpha) g + (1/alpha^2) w' z
  s <- (-1 / alpha) * g + (1 / alpha^2) * as.vector(t(w) %*% z)

  list(s = s, m_mat = m_mat)
}


#' Compute Step Derivative for Newton Update
#'
#' Computes the second-order inverse (B + lambda*I)^-2 times g, which equals
#' (B + lambda*I)^-1 times s, for the Newton derivative of the secular equation.
#'
#' @param lambda Current regularization parameter.
#' @param s Current step vector.
#' @param gamma Scaling parameter.
#' @param w_wt Precomputed w times t(w).
#' @param c_inv Precomputed inverse of c_mat.
#' @param w Low-rank update matrix.
#'
#' @return Vector equal to (B + lambda*I) inverse times s.
#'
#' @keywords internal
woodbury_step_deriv <- function(lambda, s, gamma, w_wt, c_inv, w) {
  alpha <- gamma + lambda

  # Compute w s
  w_s <- as.vector(w %*% s)

  # M = c_mat^{-1} + (1/alpha) w w' (same as in woodbury_step)
  m_mat <- c_inv + w_wt / alpha

  # Solve M z = w s
  z <- tryCatch(
    solve(m_mat, w_s),
    error = function(e) {
      rep(0, length(w_s))
    }
  )

  # (B + lambda I)^{-1} s = (-1/alpha) s + (1/alpha^2) w' z
  (-1 / alpha) * s + (1 / alpha^2) * as.vector(t(w) %*% z)
}


#' Simple Cubic Solver for Scaled Identity Hessian
#'
#' Solves the cubic subproblem when B = gamma * I (no QN updates).
#' Uses closed-form solution: s = -g / (gamma + lambda) with
#' ||s|| = lambda / sigma.
#'
#' @param g Gradient vector.
#' @param gamma Positive scalar (B = gamma * I).
#' @param sigma Cubic regularization parameter.
#' @param tol Convergence tolerance.
#' @param max_iter Maximum iterations.
#'
#' @return List with s, lambda, iterations, converged.
#'
#' @keywords internal
solve_cubic_simple <- function(g, gamma, sigma, tol = 1e-10, max_iter = 50L) {
  n <- length(g)
  g_norm <- sqrt(sum(g^2))

  if (g_norm < tol) {
    return(list(
      s = rep(0, n),
      lambda = 0,
      iterations = 0L,
      converged = TRUE
    ))
  }

  # For B = gamma * I: s = -g / (gamma + lambda)
  # Secular equation: ||s|| = ||g|| / (gamma + lambda) = lambda / sigma
  # Solving: sigma * ||g|| = lambda * (gamma + lambda)
  #          lambda^2 + gamma * lambda - sigma * ||g|| = 0
  # Using quadratic formula:
  lambda <- (-gamma + sqrt(gamma^2 + 4 * sigma * g_norm)) / 2

  # Ensure lambda >= 0
  lambda <- max(lambda, 0)

  s <- -g / (gamma + lambda)

  list(
    s = s,
    lambda = lambda,
    iterations = 1L,
    converged = TRUE
  )
}


#' Build L-BFGS Compact Form Matrices
#'
#' Constructs the w and c_mat matrices for the compact representation
#' B = gamma * I + w' c_mat w from L-BFGS history.
#'
#' @param history L-BFGS history from update_lbfgs.
#'
#' @return List with w (m x n matrix) and c_mat (m x m matrix), or NULL
#'   if history is empty.
#'
#' @details
#' The compact L-BFGS form uses a simplified SR1-like formulation:
#' B = gamma * I + Psi (Psi' S)^-1 Psi' where Psi = Y - gamma * S.
#'
#' @references
#' Nocedal, J., & Wright, S. J. (2006). Numerical Optimization (2nd ed.).
#' Springer. Section 7.2.
#'
#' @keywords internal
lbfgs_compact_form <- function(history) {
  if (is.null(history) || length(history$s) == 0) {
    return(NULL)
  }

  m <- length(history$s)
  gamma <- history$gamma

  # Build s_mat and y_mat matrices (n x m)
  s_mat <- do.call(cbind, history$s)
  y_mat <- do.call(cbind, history$y)

  # Psi = Y - gamma * S (n x m)
  psi <- y_mat - gamma * s_mat

  # w = Psi' (m x n)
  w <- t(psi)

  # Compute denominators for SR1-like form
  denoms <- sapply(seq_len(m), function(i) {
    sum(history$s[[i]] * psi[, i])
  })

  # Filter out near-zero denominators
  valid <- abs(denoms) > .Machine$double.eps * sqrt(sum(psi^2))

  if (sum(valid) == 0) {
    return(NULL)
  }

  w <- w[valid, , drop = FALSE]
  c_mat <- diag(1 / denoms[valid], nrow = sum(valid))

  list(w = w, c_mat = c_mat)
}


#' Build L-SR1 Compact Form Matrices
#'
#' Constructs the w and c_mat matrices for the compact representation
#' B = gamma * I + w' c_mat w from L-SR1 history.
#'
#' @param history L-SR1 history from update_lsr1.
#'
#' @return List with w (m x n matrix) and c_mat (m x m matrix), or NULL
#'   if history is empty.
#'
#' @details
#' For L-SR1, the compact form is:
#' B = gamma * I + Psi * (Psi' S)^-1 * Psi'
#'
#' where Psi = Y - gamma * S. This can be written as B = gamma*I + w' c_mat w
#' with w = Psi' (so w is m x n) and c_mat = (Psi' S)^-1.
#'
#' @keywords internal
lsr1_compact_form <- function(history) {
  if (is.null(history) || length(history$s) == 0) {
    return(NULL)
  }

  m <- length(history$s)
  gamma <- history$gamma

  # Build s_mat and y_mat matrices (n x m)
  s_mat <- do.call(cbind, history$s)
  y_mat <- do.call(cbind, history$y)

  # Psi = Y - gamma * S (n x m)
  psi <- y_mat - gamma * s_mat

  # w = Psi' (m x n)
  w <- t(psi)

  # c_mat = (Psi' S)^{-1} (m x m)
  psi_s <- t(psi) %*% s_mat  # m x m

  c_mat <- tryCatch(
    solve(psi_s),
    error = function(e) {
      # Use SVD-based pseudoinverse for robustness
      svd_result <- svd(psi_s)
      tol_svd <- max(svd_result$d) * .Machine$double.eps * m
      valid <- svd_result$d > tol_svd
      if (sum(valid) == 0) {
        return(NULL)
      }
      # Pseudoinverse: V * D^{-1} * U'
      svd_result$v[, valid, drop = FALSE] %*%
        diag(1 / svd_result$d[valid], nrow = sum(valid)) %*%
        t(svd_result$u[, valid, drop = FALSE])
    }
  )

  if (is.null(c_mat)) {
    return(NULL)
  }

  list(w = w, c_mat = c_mat)
}
