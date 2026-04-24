#' Solve Trust-Region Subproblem via Eigendecomposition (Algorithm 5c)
#'
#' Solves min_s m(s) = g^T s + 1/2 s^T H s  subject to  ||s|| <= radius
#' using full eigendecomposition with explicit hard-case detection.
#'
#' @param g Gradient vector (length n)
#' @param H Hessian matrix (n x n symmetric, may be indefinite)
#' @param radius Trust-region radius (positive scalar)
#' @param max_lambda_iter Maximum secular equation iterations (default: 50)
#' @param lambda_tol Secular equation convergence tolerance (default: 1e-10)
#' @param hard_case_tol Hard case detection threshold (default: 1e-8)
#'
#' @return List with:
#'   \item{s}{Solution vector}
#'   \item{lambda}{Lagrange multiplier (0 if interior Newton step)}
#'   \item{pred_reduction}{Predicted reduction: -g^T s - 1/2 s^T H s}
#'   \item{converged}{TRUE if solver returned a valid step}
#'   \item{case_type}{"interior", "easy", "hard", or "zero_gradient"}
#'   \item{on_boundary}{TRUE if ||s|| == radius (active constraint)}
#'
#' @details
#' This is the trust-region counterpart to \code{solve_cubic_eigen}. The only
#' structural difference is the constraint: the cubic solver enforces
#' \code{lambda = sigma * ||s||}, while the trust-region solver enforces
#' \code{||s|| <= radius}. Both reduce to \code{(H + lambda*I) s = -g} with
#' \code{lambda >= max(0, -lambda_min(H))}.
#'
#' Three cases are handled:
#' \itemize{
#'   \item \strong{Interior}: H is positive definite and the unconstrained
#'     Newton step \code{-H^{-1} g} lies inside the trust region. Returned
#'     with \code{lambda = 0} and \code{on_boundary = FALSE}.
#'   \item \strong{Easy boundary}: A unique \code{lambda > max(0, -lambda_min)}
#'     places \code{||s(lambda)|| = radius}. Found by Newton on the secular
#'     equation.
#'   \item \strong{Hard case}: Hessian is indefinite and the gradient is
#'     nearly orthogonal to the smallest-eigenvalue eigenvector. Resolved
#'     per Moré-Sorensen by setting \code{lambda = -lambda_min(H)} and adding
#'     a component along the smallest eigenvector to reach the boundary.
#' }
#'
#' Recommended for \code{n <= 500} due to O(n^3) eigendecomposition cost.
#'
#' @keywords internal
solve_tr_eigen <- function(g, H, radius,
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
      on_boundary = FALSE
    ))
  }

  # STEP 1: EIGENDECOMPOSITION
  # eigen() returns eigenvalues in decreasing order, reverse to ascending
  eig <- eigen(H, symmetric = TRUE)
  eigenvalues <- rev(eig$values) # Now lambda_1 <= ... <= lambda_n
  eigenvectors <- eig$vectors[, rev(seq_len(n))]

  # Transform gradient to spectral coordinates
  g_tilde <- as.vector(t(eigenvectors) %*% g)

  lambda_1 <- eigenvalues[1]
  lambda_lower <- max(0, -lambda_1)

  # STEP 2: INTERIOR CASE
  # If H is positive definite, try unconstrained Newton step first.
  # If it fits inside the trust region, that is the global model minimizer.
  if (lambda_1 > .Machine$double.eps) {
    s_newton_tilde <- -g_tilde / eigenvalues
    s_newton_norm <- sqrt(sum(s_newton_tilde^2))

    if (is.finite(s_newton_norm) && s_newton_norm <= radius) {
      s <- as.vector(eigenvectors %*% s_newton_tilde)
      h_s <- as.vector(H %*% s)
      pred_reduction <- -sum(g * s) - 0.5 * sum(s * h_s)
      return(list(
        s = s,
        lambda = 0,
        pred_reduction = pred_reduction,
        converged = TRUE,
        case_type = "interior",
        on_boundary = FALSE
      ))
    }
  }

  # STEP 3: HARD CASE DETECTION
  # Hard case: H indefinite and g nearly orthogonal to v_1 (smallest-eig vector).
  # At lambda = -lambda_1 the "standard" step s(lambda) skips the v_1 component
  # and has bounded norm. If that norm is <= radius, we add a v_1 component
  # of magnitude tau = sqrt(radius^2 - ||s_{i>=2}||^2) to reach the boundary.
  is_hard_case <- (lambda_1 <= 0) &&
    (abs(g_tilde[1]) < hard_case_tol * g_norm)

  if (is_hard_case) {
    lambda_star <- -lambda_1
    s_tilde <- numeric(n)

    # Build s along i >= 2 (skip the degenerate v_1 direction)
    if (n >= 2) {
      for (i in 2:n) {
        denom <- eigenvalues[i] + lambda_star
        if (abs(denom) > .Machine$double.eps) {
          s_tilde[i] <- -g_tilde[i] / denom
        } else {
          s_tilde[i] <- 0
        }
      }
    }

    s_partial_norm <- sqrt(sum(s_tilde^2))

    if (s_partial_norm <= radius) {
      # Add v_1 component to reach boundary
      tau <- sqrt(max(0, radius^2 - s_partial_norm^2))
      s_tilde[1] <- tau # Sign is arbitrary; v_1 is a null direction for Hs
    }
    # If s_partial_norm > radius, the "hard-hard" boundary is violated; fall
    # through to the secular equation below to find a larger lambda.

    if (s_partial_norm <= radius) {
      s <- as.vector(eigenvectors %*% s_tilde)
      h_s <- as.vector(H %*% s)
      pred_reduction <- -sum(g * s) - 0.5 * sum(s * h_s)
      return(list(
        s = s,
        lambda = lambda_star,
        pred_reduction = pred_reduction,
        converged = TRUE,
        case_type = "hard",
        on_boundary = TRUE
      ))
    }
  }

  # STEP 4: EASY BOUNDARY CASE - SECULAR EQUATION via uniroot
  # Find lambda > lambda_lower such that ||s(lambda)|| = radius, where
  # s(lambda) = -(H + lambda I)^{-1} g. Uses `uniroot` (bisection-safeguarded
  # Brent) on phi(lambda) = ||s(lambda)|| - r. phi is monotonically
  # decreasing in lambda (from +infinity at lambda = lambda_lower^+ to
  # -radius as lambda -> infinity), so uniroot is unconditionally stable
  # given a valid bracket. The previous plain-Newton iterator was
  # unstable on ill-conditioned indefinite Hessians because Newton steps
  # can overshoot into the lambda_lower neighborhood, where s explodes.

  s_tilde_at <- function(lam) {
    denom <- eigenvalues + lam
    ifelse(abs(denom) > .Machine$double.eps, -g_tilde / denom, 0)
  }

  # Lower bracket: small eps above lambda_lower; s_tilde blows up, phi > 0
  lam_lo <- lambda_lower + max(1e-8, 1e-8 * (abs(lambda_lower) + 1))
  s_lo <- s_tilde_at(lam_lo)
  norm_lo <- sqrt(sum(s_lo^2))

  # If even at tiny shift above the lower bound the step norm is already
  # below the radius, the constraint is not binding. That should have
  # been caught by the interior case, but with indefinite H we arrived
  # here from the hard-case fall-through. Use the lam_lo step as-is.
  if (!is.finite(norm_lo) || norm_lo <= radius) {
    s <- as.vector(eigenvectors %*% s_lo)
    h_s <- as.vector(H %*% s)
    pred_reduction <- -sum(g * s) - 0.5 * sum(s * h_s)
    return(list(
      s = s, lambda = lam_lo,
      pred_reduction = pred_reduction,
      converged = TRUE,
      case_type = "easy",
      on_boundary = is.finite(norm_lo) && abs(norm_lo - radius) <=
        lambda_tol * max(1, radius)
    ))
  }

  # Upper bracket: grow lambda until phi < 0 (i.e., ||s|| < radius).
  lam_hi <- max(1, abs(lambda_lower)) * 10
  for (expand in seq_len(60)) {
    s_hi <- s_tilde_at(lam_hi)
    norm_hi <- sqrt(sum(s_hi^2))
    if (is.finite(norm_hi) && norm_hi < radius) break
    lam_hi <- lam_hi * 10
  }

  phi <- function(lam) {
    sqrt(sum(s_tilde_at(lam)^2)) - radius
  }

  lambda <- tryCatch(
    stats::uniroot(phi, lower = lam_lo, upper = lam_hi,
                   tol = lambda_tol * max(1, radius),
                   maxiter = max_lambda_iter * 4)$root,
    error = function(e) NA_real_
  )

  if (is.na(lambda)) {
    # Fall back: use a best-effort shifted solve at lam_hi
    lambda <- lam_hi
    converged <- FALSE
  } else {
    converged <- TRUE
  }

  s_tilde_final <- s_tilde_at(lambda)
  s <- as.vector(eigenvectors %*% s_tilde_final)
  h_s <- as.vector(H %*% s)
  pred_reduction <- -sum(g * s) - 0.5 * sum(s * h_s)
  s_norm_final <- sqrt(sum(s_tilde_final^2))

  list(
    s = s,
    lambda = lambda,
    pred_reduction = pred_reduction,
    converged = converged,
    case_type = "easy",
    on_boundary = is.finite(s_norm_final) &&
      abs(s_norm_final - radius) <= lambda_tol * max(1, radius) + 1e-6
  )
}
