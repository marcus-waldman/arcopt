#' Wolfe Line Search
#'
#' Performs a line search along a descent direction `d` from point `x`,
#' finding a step length `alpha` that satisfies the strong Wolfe conditions.
#' Implements the bracketing + zoom scheme from Nocedal & Wright (2006),
#' Algorithms 3.5 (bracketing) and 3.6 (zoom).
#'
#' @param fn Objective function. Called as `fn(x_new, ...)`.
#' @param gr Gradient function. Called as `gr(x_new, ...)`.
#' @param x Numeric vector; current iterate.
#' @param d Numeric vector; search direction (must be a descent direction,
#'   i.e., `sum(g_x * d) < 0`).
#' @param f_x Current objective value `fn(x)`.
#' @param g_x Current gradient `gr(x)`.
#' @param c1 Armijo (sufficient decrease) constant; default `1e-4`.
#' @param c2 Curvature constant (strong Wolfe); default `0.9` (standard for
#'   quasi-Newton). Smaller values give stricter curvature conditions; for
#'   linear-CG use `0.1`.
#' @param alpha_max Maximum step length; default `1.0`. Line search will
#'   try `alpha_curr = alpha_max` first (appropriate for Newton-like
#'   directions).
#' @param max_iter Maximum total evaluations (bracketing + zoom); default
#'   `20`.
#' @param ... Extra arguments forwarded to `fn` and `gr`.
#'
#' @return A list with components:
#'   \item{success}{TRUE if strong Wolfe was satisfied; FALSE otherwise.}
#'   \item{alpha}{Final step length (best available even on failure).}
#'   \item{x_new}{`x + alpha * d`.}
#'   \item{f_new}{`fn(x_new)`.}
#'   \item{g_new}{`gr(x_new)`.}
#'   \item{evals_f}{Number of objective evaluations performed.}
#'   \item{evals_g}{Number of gradient evaluations performed.}
#'   \item{reason}{String describing success or failure mode.}
#'
#' @details
#' Strong Wolfe conditions:
#' \itemize{
#'   \item Armijo (sufficient decrease):
#'     \eqn{f(x + \alpha d) \le f(x) + c_1 \alpha g_x^\top d}
#'   \item Curvature:
#'     \eqn{|g(x + \alpha d)^\top d| \le c_2 |g_x^\top d|}
#' }
#' The algorithm first brackets an interval containing a Wolfe point, then
#' bisects within the bracket to locate one. Bisection is used for
#' simplicity; quadratic or cubic interpolation as in Nocedal & Wright §3.5
#' can be added as a future refinement.
#'
#' Returns `success = FALSE` when no strong Wolfe alpha is found within
#' `max_iter` evaluations. The caller can still consult `alpha` (and the
#' associated trial point) to decide whether to take the best-known step
#' or reject the iteration.
#'
#' @keywords internal
wolfe_line_search <- function(fn, gr, x, d, f_x, g_x,
                              c1 = 1e-4, c2 = 0.9,
                              alpha_max = 1.0, max_iter = 20, ...) {
  # Directional derivative at alpha = 0
  phi_0 <- f_x
  phi_prime_0 <- sum(g_x * d)

  if (!is.finite(phi_prime_0) || phi_prime_0 >= 0) {
    # Not a descent direction (or numerical issue); caller should handle
    return(list(
      success = FALSE,
      alpha = 0,
      x_new = x,
      f_new = f_x,
      g_new = g_x,
      evals_f = 0L,
      evals_g = 0L,
      reason = "not_descent_direction"
    ))
  }

  # Evaluation bookkeeping
  evals_f <- 0L
  evals_g <- 0L

  # Evaluate phi(alpha) = f(x + alpha*d) and phi'(alpha) = g(x + alpha*d)^T d
  eval_alpha <- function(alpha) {
    x_trial <- x + alpha * d
    f_trial <- fn(x_trial, ...)
    evals_f <<- evals_f + 1L
    if (!is.finite(f_trial)) {
      return(list(
        phi = f_trial, phi_prime = NA_real_,
        x_new = x_trial, f_new = f_trial, g_new = g_x,
        finite = FALSE
      ))
    }
    g_trial <- gr(x_trial, ...)
    evals_g <<- evals_g + 1L
    list(
      phi = f_trial,
      phi_prime = sum(g_trial * d),
      x_new = x_trial,
      f_new = f_trial,
      g_new = g_trial,
      finite = is.finite(sum(g_trial))
    )
  }

  # Zoom phase (Nocedal-Wright Algorithm 3.6) -- bisection within bracket.
  # info_lo corresponds to alpha_lo (must satisfy Armijo); info_hi to alpha_hi.
  zoom <- function(alpha_lo, alpha_hi, info_lo, iter_remaining) {
    info_last <- info_lo
    alpha_last <- alpha_lo
    for (zi in seq_len(iter_remaining)) {
      if (evals_f >= max_iter) {
        break
      }
      # Bisection
      alpha <- (alpha_lo + alpha_hi) / 2
      if (!is.finite(alpha) || abs(alpha_hi - alpha_lo) < .Machine$double.eps) {
        break
      }
      info <- eval_alpha(alpha)
      info_last <- info
      alpha_last <- alpha

      if (!info$finite) {
        # Non-finite in the middle of bracket: shrink toward lo
        alpha_hi <- alpha
        next
      }

      # Armijo violated OR non-monotone relative to lo?
      if (info$phi > phi_0 + c1 * alpha * phi_prime_0 ||
            info$phi >= info_lo$phi) {
        alpha_hi <- alpha
      } else {
        # Armijo OK -- check curvature
        if (abs(info$phi_prime) <= -c2 * phi_prime_0) {
          return(list(success = TRUE, alpha = alpha, info = info,
                      reason = "wolfe_satisfied"))
        }
        # Curvature not met -- shift bracket toward the side of the minimum
        if (info$phi_prime * (alpha_hi - alpha_lo) >= 0) {
          alpha_hi <- alpha_lo
        }
        alpha_lo <- alpha
        info_lo <- info
      }
    }
    list(success = FALSE, alpha = alpha_last, info = info_last,
         reason = "zoom_max_iter")
  }

  # Bracketing phase (Nocedal-Wright Algorithm 3.5)
  alpha_prev <- 0
  info_prev <- list(phi = phi_0, phi_prime = phi_prime_0,
                    x_new = x, f_new = f_x, g_new = g_x,
                    finite = TRUE)
  alpha_curr <- alpha_max

  result <- NULL
  for (iter in seq_len(max_iter)) {
    if (evals_f >= max_iter) {
      break
    }

    info_curr <- eval_alpha(alpha_curr)

    if (!info_curr$finite) {
      # Non-finite trial: bracket [0, alpha_curr] and zoom
      z <- zoom(alpha_prev, alpha_curr, info_prev, max_iter - evals_f)
      result <- z
      break
    }

    # Armijo violated, or non-monotone -> zoom in [alpha_prev, alpha_curr]
    if (info_curr$phi > phi_0 + c1 * alpha_curr * phi_prime_0 ||
          (iter > 1 && info_curr$phi >= info_prev$phi)) {
      z <- zoom(alpha_prev, alpha_curr, info_prev, max_iter - evals_f)
      result <- z
      break
    }

    # Strong Wolfe check
    if (abs(info_curr$phi_prime) <= -c2 * phi_prime_0) {
      result <- list(success = TRUE, alpha = alpha_curr, info = info_curr,
                     reason = "wolfe_satisfied")
      break
    }

    # Gradient flipped sign -> minimum is between alpha_prev and alpha_curr
    if (info_curr$phi_prime >= 0) {
      z <- zoom(alpha_curr, alpha_prev, info_curr, max_iter - evals_f)
      result <- z
      break
    }

    # Extend the bracket; standard choice is 2 * alpha_prev, modestly
    # capped so we don't walk off the problem
    alpha_prev <- alpha_curr
    info_prev <- info_curr
    alpha_curr <- min(2 * alpha_curr, alpha_max * 10)
    if (alpha_curr == alpha_prev) {
      # Hit the extension cap without finding a Wolfe point: return
      # best info so far (this is an Armijo-but-not-curvature case).
      result <- list(
        success = FALSE, alpha = alpha_prev, info = info_prev,
        reason = "alpha_extension_capped"
      )
      break
    }
  }

  if (is.null(result)) {
    result <- list(
      success = FALSE, alpha = alpha_prev, info = info_prev,
      reason = "bracket_max_iter"
    )
  }

  info <- result$info
  list(
    success = result$success,
    alpha = result$alpha,
    x_new = info$x_new,
    f_new = info$f_new,
    g_new = info$g_new,
    evals_f = evals_f,
    evals_g = evals_g,
    reason = result$reason
  )
}
