#' Flat-Ridge Detector for Cubic-to-Trust-Region Fallback
#'
#' Detects when cubic regularization has entered a flat-ridge stagnation
#' regime: gradient stalls in a near-singular positive-definite region
#' where the cubic penalty has no bite. Triggers the trust-region fallback
#' in the main arcopt iteration when fired.
#'
#' @details
#' The detector maintains a sliding window of runtime diagnostics and fires
#' when \strong{all four} of the following signals hold throughout the window:
#' \enumerate{
#'   \item \strong{Sigma pinned at floor}: \code{sigma_k <= 10 * sigma_min}
#'     for every iteration in the window.
#'   \item \strong{Near-perfect model}: \code{|rho_k - 1| < rho_tol} for
#'     every iteration in the window (model reduction matches actual
#'     reduction).
#'   \item \strong{Gradient stagnant}: The ratio of the most recent
#'     \code{||g||_inf} to the oldest in the window exceeds
#'     \code{grad_decrease_max} (less than 10 percent decrease over the
#'     window by default).
#'   \item \strong{Hessian ill-conditioned}: The smallest Hessian
#'     eigenvalue at the most recent iteration satisfies
#'     \code{lambda_min < tol_ridge}. This includes both the classical
#'     "flat-ridge" regime (small positive \code{lambda_min}) and the
#'     "stuck-at-indefinite-saddle" regime (any negative
#'     \code{lambda_min}), because cubic regularization loses its grip
#'     in both.
#' }
#'
#' An absolute floor \code{g_inf_floor} is also enforced: the trigger will
#' not fire if \code{||g||_inf} is already below the typical convergence
#' tolerance, preventing spurious triggers at true local minima.
#'
#' The detector is an empirical proxy for local-error-bound (EB) violation
#' (Yue, Zhou, So 2018). Under EB, cubic regularization attains Q-quadratic
#' convergence even at degenerate minima; the detector fires precisely in
#' the regime where EB is practically vacuous and no theoretical convergence
#' guarantee applies.
#'
#' @name flat_ridge
NULL

#' Initialize Flat-Ridge Detector State
#'
#' @param window Integer sliding-window length (default 10)
#'
#' @return A list with empty diagnostic vectors and the window size.
#'
#' @seealso \code{\link{update_flat_ridge_state}},
#'   \code{\link{check_flat_ridge_trigger}}
#' @keywords internal
init_flat_ridge_state <- function(window = 10) {
  stopifnot(is.numeric(window), length(window) == 1, window >= 1)
  list(
    window = as.integer(window),
    sigma = numeric(0),
    rho = numeric(0),
    g_inf = numeric(0),
    lambda_min = numeric(0)
  )
}

#' Update Flat-Ridge Detector State with One Iteration's Diagnostics
#'
#' Pushes a new record onto the sliding window and trims the oldest entry
#' once the window is full.
#'
#' @param state List returned by \code{\link{init_flat_ridge_state}}.
#' @param sigma Current regularization parameter.
#' @param rho Current step-acceptance ratio.
#' @param g_inf Current \code{||g||_inf}.
#' @param lambda_min Current smallest Hessian eigenvalue (may be negative
#'   for indefinite Hessians; pass \code{NA_real_} if unavailable).
#'
#' @return The updated state list.
#'
#' @keywords internal
update_flat_ridge_state <- function(state, sigma, rho, g_inf, lambda_min) {
  state$sigma <- c(state$sigma, sigma)
  state$rho <- c(state$rho, rho)
  state$g_inf <- c(state$g_inf, g_inf)
  state$lambda_min <- c(state$lambda_min, lambda_min)

  w <- state$window
  if (length(state$sigma) > w) {
    state$sigma <- utils::tail(state$sigma, w)
    state$rho <- utils::tail(state$rho, w)
    state$g_inf <- utils::tail(state$g_inf, w)
    state$lambda_min <- utils::tail(state$lambda_min, w)
  }
  state
}

#' Check Flat-Ridge Trigger
#'
#' Evaluates the four-signal rule against the current sliding window.
#'
#' @param state List returned by \code{\link{init_flat_ridge_state}} and
#'   maintained by \code{\link{update_flat_ridge_state}}.
#' @param sigma_min Minimum regularization floor used by the orchestrator.
#' @param tol_ridge Upper bound on \code{lambda_min} that counts as
#'   "poorly conditioned" (default \code{1e-3}).
#' @param rho_tol Maximum \code{|rho - 1|} that counts as "near-perfect
#'   model" (default \code{0.1}).
#' @param grad_decrease_max Ratio of latest to oldest \code{||g||_inf}
#'   above which the gradient is deemed stagnant (default \code{0.9}).
#' @param g_inf_floor Absolute lower bound on \code{||g||_inf}; the trigger
#'   will not fire below this value to avoid firing at true local minima
#'   (default \code{1e-6}).
#'
#' @return \code{TRUE} if the trust-region fallback should be activated,
#'   \code{FALSE} otherwise. Returns \code{FALSE} whenever the window is
#'   not yet full or any diagnostic is non-finite.
#'
#' @keywords internal
check_flat_ridge_trigger <- function(state,
                                     sigma_min,
                                     tol_ridge = 1e-3,
                                     rho_tol = 0.1,
                                     grad_decrease_max = 0.9,
                                     g_inf_floor = 1e-6) {
  w <- state$window
  if (length(state$sigma) < w) {
    return(FALSE)
  }

  # Absolute gradient floor: bail out if we are already near convergence
  latest_g <- state$g_inf[w]
  if (!is.finite(latest_g) || latest_g < g_inf_floor) {
    return(FALSE)
  }

  # Signal 1: sigma pinned at floor throughout the window
  sigma_threshold <- 10 * sigma_min
  if (!all(is.finite(state$sigma)) ||
        !all(state$sigma <= sigma_threshold)) {
    return(FALSE)
  }

  # Signal 2: near-perfect model (|rho - 1| < rho_tol) throughout
  if (!all(is.finite(state$rho)) ||
        !all(abs(state$rho - 1) < rho_tol)) {
    return(FALSE)
  }

  # Signal 3: gradient stagnant (decrease < 1 - grad_decrease_max over window)
  oldest_g <- state$g_inf[1]
  if (!is.finite(oldest_g) || oldest_g < .Machine$double.eps) {
    return(FALSE)
  }
  if (latest_g / oldest_g <= grad_decrease_max) {
    return(FALSE)
  }

  # Signal 4: Hessian ill-conditioned at the current iteration.
  # Covers both "flat-ridge" (0 < lambda_min < tol_ridge) and
  # "stuck-indefinite-saddle" (lambda_min < 0) regimes. Cubic
  # regularization loses grip in both; well-conditioned PD
  # (lambda_min >= tol_ridge) does not fire.
  latest_lmin <- state$lambda_min[w]
  if (!is.finite(latest_lmin) || latest_lmin >= tol_ridge) {
    return(FALSE)
  }

  TRUE
}
