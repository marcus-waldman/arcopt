#' Healthy-Basin Detector for Cubic-to-QN-Polish Transition
#'
#' Detects when the iterate has entered the quadratic attraction basin of a
#' strict local minimum, at which point arcopt can switch from cubic
#' regularization to line-search BFGS (qn_polish mode) and skip further
#' Hessian evaluations for the remainder of the run.
#'
#' @details
#' The detector maintains a sliding window of runtime diagnostics and fires
#' when \strong{all five} of the following signals hold throughout the
#' window:
#' \enumerate{
#'   \item \strong{Newton step accepted}: \code{step_type == "newton"} at
#'     every iteration in the window. Indicates the cubic subproblem has
#'     not been invoked (Cholesky succeeded and the Newton step passed the
#'     ratio test).
#'   \item \strong{High ratio}: \code{rho_k >= rho_polish} (default
#'     \code{0.9}) throughout the window. The quadratic model is accurately
#'     predicting actual reduction.
#'   \item \strong{Hessian well-conditioned PD}: \code{lambda_min(H_k) >=
#'     lambda_min_polish} (default \code{1e-3}) throughout the window. Not
#'     just "Cholesky succeeded" but strictly bounded away from zero.
#'   \item \strong{Superlinear gradient decay}: \code{||g_k||_inf /
#'     ||g_{k-1}||_inf <= g_decay_polish} (default \code{0.5}) for every
#'     consecutive pair in the window. This is the strongest signal of the
#'     quadratic attraction basin.
#'   \item \strong{Above convergence floor}: \code{||g_0||_inf >
#'     g_inf_floor_polish} at window start (default \code{1e-8}). Prevents
#'     firing at convergence where Newton-accept and gradient-decay are
#'     trivially satisfied.
#' }
#'
#' Contrast with the flat-ridge detector (\code{\link{flat_ridge}}): that
#' detector fires on \emph{stagnation} (gradient not decreasing); this one
#' fires on \emph{healthy convergence} (gradient decreasing super-
#' linearly). Both are expected to be mutually exclusive in practice.
#'
#' @name healthy_basin
#' @keywords internal
NULL

#' Initialize Healthy-Basin Detector State
#'
#' @param window Integer sliding-window length (default 5; shorter than
#'   the flat-ridge detector's 10 because the basin signal is stronger
#'   and revert-on-false-positive is cheap).
#'
#' @return A list with empty diagnostic vectors and the window size.
#'
#' @seealso \code{\link{update_healthy_basin_state}},
#'   \code{\link{check_healthy_basin_trigger}}
#' @keywords internal
init_healthy_basin_state <- function(window = 5) {
  stopifnot(is.numeric(window), length(window) == 1, window >= 1)
  list(
    window = as.integer(window),
    used_newton = logical(0),
    rho = numeric(0),
    lambda_min = numeric(0),
    g_inf = numeric(0)
  )
}

#' Update Healthy-Basin Detector State with One Iteration's Diagnostics
#'
#' Pushes a new record onto the sliding window and trims the oldest entry
#' once the window is full.
#'
#' @param state List returned by \code{\link{init_healthy_basin_state}}.
#' @param used_newton Logical; \code{TRUE} if the current iteration
#'   accepted a Newton step (i.e., \code{step_type == "newton"}).
#' @param rho Current step-acceptance ratio (NA if step was rejected).
#' @param lambda_min Current smallest Hessian eigenvalue.
#' @param g_inf Current \code{||g||_inf}.
#'
#' @return The updated state list.
#'
#' @keywords internal
update_healthy_basin_state <- function(state, used_newton, rho,
                                       lambda_min, g_inf) {
  state$used_newton <- c(state$used_newton, isTRUE(used_newton))
  state$rho <- c(state$rho, rho)
  state$lambda_min <- c(state$lambda_min, lambda_min)
  state$g_inf <- c(state$g_inf, g_inf)

  w <- state$window
  if (length(state$used_newton) > w) {
    state$used_newton <- utils::tail(state$used_newton, w)
    state$rho <- utils::tail(state$rho, w)
    state$lambda_min <- utils::tail(state$lambda_min, w)
    state$g_inf <- utils::tail(state$g_inf, w)
  }
  state
}

#' Check Healthy-Basin Trigger
#'
#' Evaluates the five-signal rule against the current sliding window.
#'
#' @param state List returned by \code{\link{init_healthy_basin_state}} and
#'   maintained by \code{\link{update_healthy_basin_state}}.
#' @param rho_polish Minimum acceptable \code{rho} throughout the window
#'   (default \code{0.9}).
#' @param lambda_min_polish Minimum acceptable \code{lambda_min(H)}
#'   throughout the window (default \code{1e-3}).
#' @param g_decay_polish Maximum acceptable ratio of consecutive
#'   \code{||g||_inf} values (default \code{0.5} = 2x-per-iter contraction).
#' @param g_inf_floor_polish Absolute lower bound on \code{||g||_inf} at
#'   window start (default \code{1e-8}); prevents firing at convergence.
#'
#' @return \code{TRUE} if the qn_polish transition should be activated,
#'   \code{FALSE} otherwise. Returns \code{FALSE} whenever the window is
#'   not yet full or any diagnostic is non-finite.
#'
#' @keywords internal
check_healthy_basin_trigger <- function(state,
                                        rho_polish = 0.9,
                                        lambda_min_polish = 1e-3,
                                        g_decay_polish = 0.5,
                                        g_inf_floor_polish = 1e-8) {
  w <- state$window
  if (length(state$used_newton) < w) {
    return(FALSE)
  }

  # Absolute gradient floor at window start: must be non-trivial
  oldest_g <- state$g_inf[1]
  if (!is.finite(oldest_g) || oldest_g <= g_inf_floor_polish) {
    return(FALSE)
  }

  # Signal 1: Newton step accepted throughout the window
  if (!all(state$used_newton)) {
    return(FALSE)
  }

  # Signal 2: high rho throughout
  if (!all(is.finite(state$rho)) || !all(state$rho >= rho_polish)) {
    return(FALSE)
  }

  # Signal 3: lambda_min bounded above lambda_min_polish throughout
  if (!all(is.finite(state$lambda_min)) ||
        !all(state$lambda_min >= lambda_min_polish)) {
    return(FALSE)
  }

  # Signal 4: superlinear gradient decay between consecutive iterations.
  # Requires ratio r_i = g_inf[i] / g_inf[i-1] <= g_decay_polish for every
  # i = 2..w. Gradient must be strictly positive.
  if (!all(is.finite(state$g_inf)) || any(state$g_inf <= 0)) {
    return(FALSE)
  }
  ratios <- state$g_inf[-1] / state$g_inf[-w]
  if (!all(is.finite(ratios)) || !all(ratios <= g_decay_polish)) {
    return(FALSE)
  }

  TRUE
}
