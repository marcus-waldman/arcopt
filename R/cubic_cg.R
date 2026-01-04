#' Solve Cubic Subproblem via ARCqK Multi-Shift CG-Lanczos (Algorithm 5b)
#'
#' Solves min_s m(s) = g^T s + 1/2 s^T H s + sigma/3 ||s||^3 using multi-shift
#' CG-Lanczos with automatic negative curvature detection.
#'
#' @param g Gradient vector (length n)
#' @param hess_vec Function computing Hessian-vector products: hess_vec(v) -> Hv
#' @param sigma Regularization parameter (positive scalar)
#' @param shift_psi Shift spacing ratio (default: sqrt(10) ≈ 3.16)
#' @param shift_e_lower Min shift exponent (default: -10 → 10^-10)
#' @param shift_e_upper Max shift exponent (default: 20 → 10^20)
#' @param max_cg_iter Maximum CG iterations (default: min(100, 2*n))
#' @param residual_zeta Residual tolerance exponent (default: 0.5)
#' @param residual_xi Residual tolerance constant (default: 1e-6)
#'
#' @return List with:
#'   \item{s}{Solution vector}
#'   \item{lambda}{Lagrange multiplier (selected shift)}
#'   \item{pred_reduction}{Predicted reduction: -g^T s - 1/2 s^T H s - sigma/3 ||s||^3}
#'   \item{converged}{TRUE if CG converged for selected shift}
#'   \item{cg_iters}{Number of CG iterations performed}
#'
#' @details
#' **EXPERIMENTAL SOLVER** - This solver is designed for large-scale problems
#' (n > 500) and may perform poorly on small problems. Benchmarks show it can
#' produce large errors and numerical instability on problems with n < 100.
#' Use solve_cubic_eigen() for most applications.
#'
#' Uses CG-Lanczos with multiple shifts (typically 61) geometrically spaced
#' from 10^-10 to 10^20. Single Hessian-vector product per iteration shared
#' across all shifts. Negative curvature detected via gamma_j < 0. Recommended
#' only for large-scale problems (n > 500) or matrix-free optimization.
#'
#' @keywords internal
solve_cubic_cg <- function(g, hess_vec, sigma,
                           shift_psi = sqrt(10),
                           shift_e_lower = -10,
                           shift_e_upper = 20,
                           max_cg_iter = NULL,
                           residual_zeta = 0.5,
                           residual_xi = 1e-6) {
  n <- length(g)

  # STEP 0: HANDLE ZERO GRADIENT
  g_norm <- sqrt(sum(g^2))
  if (g_norm < .Machine$double.eps) {
    return(list(
      s = rep(0, n),
      lambda = 0,
      pred_reduction = 0,
      converged = TRUE,
      cg_iters = 0
    ))
  }

  # Set default max_cg_iter
  if (is.null(max_cg_iter)) {
    max_cg_iter <- min(100, 2 * n)
  }

  # STEP 1: DISCRETIZE SHIFTS
  # Number of shifts to span from 10^e_L to 10^e_U with spacing psi
  m <- ceiling((shift_e_upper - shift_e_lower) / log10(shift_psi))
  lambda_0 <- 10^shift_e_lower  # Minimum shift (e.g., 10^-10)

  shifts <- numeric(m + 1)
  active <- logical(m + 1)
  converged <- logical(m + 1)

  for (i in seq_len(m + 1)) {
    shifts[i] <- lambda_0 * shift_psi^(i - 1)
    active[i] <- TRUE
    converged[i] <- FALSE
  }

  # STEP 2: INITIALIZE CG-LANCZOS
  # Solution vectors: x^(i) for each shift
  x <- matrix(0, n, m + 1)

  # Lanczos vectors
  mu_0 <- g_norm
  v_j <- g / mu_0
  v_prev <- rep(0, n)

  # Search directions: p^(i) for each shift
  # Initialize with negative gradient for minimization
  p <- matrix(0, n, m + 1)
  for (i in seq_len(m + 1)) {
    p[, i] <- -g
  }

  # Residual norms and ||p||^2 values
  sigma_res <- rep(mu_0, m + 1)
  pi_vals <- rep(mu_0^2, m + 1)

  # Recurrence variables (per shift)
  omega_prev <- rep(0, m + 1)
  gamma_prev <- rep(1, m + 1)
  omega <- numeric(m + 1)
  gamma <- numeric(m + 1)

  # STEP 3: CG-LANCZOS ITERATION
  cg_iters <- 0
  for (j in seq_len(max_cg_iter)) {
    cg_iters <- j

    # === LANCZOS PART (shift-independent) ===
    Hv_j <- hess_vec(v_j)
    delta_j <- sum(v_j * Hv_j)
    temp <- Hv_j - delta_j * v_j - mu_0 * v_prev

    # Update for next iteration
    mu_j <- mu_0
    mu_0 <- sqrt(sum(temp^2))

    if (mu_0 < .Machine$double.eps) {
      break  # Lanczos breakdown
    }

    v_prev <- v_j
    v_j <- temp / mu_0

    # === CG PART (for each active shift) ===
    any_active <- FALSE

    for (i in seq_len(m + 1)) {
      if (!active[i]) next

      # Shifted diagonal element
      delta_j_i <- delta_j + shifts[i]

      # CG recurrence coefficient
      if (abs(gamma_prev[i]) > .Machine$double.eps) {
        gamma[i] <- 1 / (delta_j_i - omega_prev[i] / gamma_prev[i])
      } else {
        gamma[i] <- 0
      }

      # NEGATIVE CURVATURE CHECK
      if (gamma[i] < 0) {
        active[i] <- FALSE
        next
      }

      # Update recurrence variables
      omega[i] <- (mu_0 * gamma[i])^2

      # Update residual norm
      sigma_res_new <- -mu_0 * gamma[i] * sigma_res[i]

      # Update solution
      x[, i] <- x[, i] + gamma[i] * p[, i]

      # Update search direction
      p[, i] <- sigma_res_new * v_j + omega[i] * p[, i]

      # Update ||p||^2
      pi_vals[i] <- sigma_res_new^2 + omega[i]^2 * pi_vals[i]

      # CONVERGENCE CHECK
      x_norm <- sqrt(sum(x[, i]^2))
      tol_i <- residual_xi * min(g_norm, x_norm)^(1 + residual_zeta)

      if (abs(sigma_res_new) <= tol_i) {
        converged[i] <- TRUE
      }

      # Update residual norm for next iteration
      sigma_res[i] <- sigma_res_new

      any_active <- TRUE
    }

    # Update previous iteration values for next iteration
    omega_prev <- omega
    gamma_prev <- gamma

    if (!any_active) {
      break  # All systems terminated
    }
  }

  # STEP 4: SELECT BEST SHIFT
  # Find smallest non-indefinite shift index
  i_plus <- which(active | converged)[1]
  if (is.na(i_plus)) {
    # Fallback: all shifts became indefinite, use first one
    i_plus <- 1
  }

  # Among valid shifts, find closest to satisfying secular equation
  best_i <- i_plus
  best_error <- Inf

  for (i in i_plus:(m + 1)) {
    if (active[i] || converged[i]) {
      x_norm <- sqrt(sum(x[, i]^2))
      error <- abs(sigma * x_norm - shifts[i])
      if (error < best_error) {
        best_error <- error
        best_i <- i
      }
    }
  }

  s <- x[, best_i]

  # STEP 5: COMPUTE PREDICTED REDUCTION
  Hs <- hess_vec(s)
  s_norm <- sqrt(sum(s^2))
  pred_reduction <- -sum(g * s) - 0.5 * sum(s * Hs) - (sigma / 3) * s_norm^3

  return(list(
    s = s,
    lambda = shifts[best_i],
    pred_reduction = pred_reduction,
    converged = converged[best_i],
    cg_iters = cg_iters
  ))
}

#' Create Hessian-Vector Product Function from Full Matrix
#'
#' Wraps a full Hessian matrix in a function interface for use with
#' solve_cubic_cg.
#'
#' @param H Hessian matrix (n x n)
#'
#' @return Function that computes H %*% v for any vector v
#'
#' @keywords internal
hess_vec_fd <- function(H) {
  function(v) as.vector(H %*% v)
}
