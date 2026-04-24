# Advanced control parameters topic page
# =======================================
#
# Documentation-only file. Provides ?arcopt_advanced_controls covering
# Tier 2 controls — cubic regularization tuning, the trust-region
# fallback, the quasi-Newton polish mode, and the QN-variant routing
# parameters. ?arcopt itself documents only the Tier 1 controls.

#' Advanced Control Parameters for arcopt
#'
#' The main \code{\link{arcopt}} help page documents only the user-facing
#' tolerances and switches (`maxit`, `gtol_abs`, `ftol_abs`, `xtol_abs`,
#' `trace`, `verbose`, `use_qn`). This page documents every other entry
#' that the `control` list of `arcopt()` (and the routed-to
#' `arcopt:::arcopt_qn()` variant) recognizes, organized by the
#' subsystem each parameter governs.
#'
#' @section Cubic regularization:
#' These parameters control the adaptive cubic model
#' \eqn{m_k(s) = f_k + g_k^\top s + \tfrac{1}{2} s^\top H_k s +
#' \tfrac{\sigma_k}{3} \|s\|^3} and the sigma adaptation rule (Algorithm
#' 2a/2b of `design/pseudocode.qmd`).
#' \describe{
#'   \item{`sigma0`}{Initial regularization parameter (default `1.0`).}
#'   \item{`sigma_min`}{Floor on `sigma_k` (default `1e-6`). Prevents the
#'     cubic term from vanishing entirely on flat regions.}
#'   \item{`sigma_max`}{Ceiling on `sigma_k` (default `1e12`). Triggers
#'     emergency-stop behavior when reached.}
#'   \item{`eta1`}{Step-acceptance threshold; steps with \eqn{\rho_k \ge
#'     \texttt{eta1}} are accepted (default `0.1`).}
#'   \item{`eta2`}{Very-successful threshold; \eqn{\rho_k \ge
#'     \texttt{eta2}} triggers `sigma` shrinkage (default `0.9`).}
#'   \item{`gamma1`}{Multiplicative shrink factor on a very-successful
#'     step (default `0.5`).}
#'   \item{`gamma2`}{Multiplicative grow factor on an unsuccessful step
#'     (default `2.0`).}
#' }
#'
#' @section Trust-region fallback (cubic to TR on flat-ridge detection):
#' Cubic regularization can stagnate in "flat-ridge" regimes — iterations
#' with the regularization floor pinned, model predictions matching the
#' objective (\eqn{\rho \approx 1}), gradient stalling, and a Hessian
#' that is positive-definite but nearly singular. This is outside the
#' local-error-bound condition of Yue, Zhou & So (2018) under which
#' cubic regularization is guaranteed to converge quadratically at
#' degenerate minimizers. arcopt detects the regime and switches once
#' from the cubic subproblem to a trust-region subproblem.
#' \describe{
#'   \item{`tr_fallback_enabled`}{One-way cubic to TR switch
#'     (default `TRUE`).}
#'   \item{`tr_fallback_window`}{Sliding-window length for the detector
#'     (default `10`).}
#'   \item{`tr_fallback_tol_ridge`}{`lambda_min(H)` threshold defining a
#'     "near-singular PD" Hessian (default `1e-3`).}
#'   \item{`tr_fallback_rho_tol`}{Tolerance on `|rho - 1|` for the
#'     "near-perfect model" signal (default `0.1`).}
#'   \item{`tr_fallback_grad_decrease_max`}{Ratio of latest to oldest
#'     `||g||_inf` above which the gradient counts as stagnant
#'     (default `0.9`).}
#'   \item{`tr_fallback_g_inf_floor`}{Absolute lower bound on `||g||_inf`
#'     below which the switch will not fire — keeps the hybrid from
#'     triggering at true local minima (default `1e-6`).}
#'   \item{`tr_r0`}{Initial trust-region radius at the switch
#'     (default `1.0`).}
#'   \item{`tr_rmax`}{Maximum trust-region radius (default `100`).}
#'   \item{`tr_eta1`}{TR step-acceptance threshold (default `0.25`).}
#'   \item{`tr_eta2`}{TR expansion threshold (default `0.75`).}
#'   \item{`tr_gamma_shrink`}{Radius shrink factor on a poor step
#'     (default `0.25`).}
#'   \item{`tr_gamma_grow`}{Radius grow factor on a very-good
#'     boundary step (default `2.0`).}
#' }
#'
#' @section Quasi-Newton polish (cubic to BFGS line-search):
#' Once the iterate enters the quadratic attraction basin of a strict
#' local minimum, the cubic regularization penalty has decayed to its
#' floor and contributes negligible damping, but arcopt still evaluates
#' `hess()` every iteration. For expensive Hessians (analytic AD via
#' Stan, finite differences) this dominates wall-clock time. Polish
#' mode replaces the cubic subproblem with a Wolfe line search along the
#' BFGS-approximated Newton direction, skipping further `hess()` calls
#' until convergence or until the BFGS approximation drifts.
#'
#' Off by default in v0.2.0 because the existing manuscript and benchmark
#' problems converge before the five-signal healthy-basin detector
#' accumulates a full window. Enable opt-in for long-running smooth
#' problems with expensive Hessians.
#' \describe{
#'   \item{`qn_polish_enabled`}{Enable the cubic to qn_polish bidirectional
#'     switch (default `FALSE`).}
#'   \item{`qn_polish_window`}{Sliding-window length for the healthy-
#'     basin detector (default `5`).}
#'   \item{`qn_polish_rho`}{Minimum `rho_k` required throughout the
#'     window (default `0.9`).}
#'   \item{`qn_polish_lambda_min`}{Minimum `lambda_min(H_k)` required
#'     throughout the window (default `1e-3`).}
#'   \item{`qn_polish_g_decay`}{Maximum ratio of consecutive `||g||_inf`
#'     values; e.g. `0.5` requires 2x-per-iteration contraction
#'     (default `0.5`).}
#'   \item{`qn_polish_g_inf_floor`}{Absolute lower bound on `||g||_inf`
#'     at window start; prevents firing at convergence (default `1e-8`).}
#'   \item{`qn_polish_c1`, `qn_polish_c2`}{Wolfe line-search constants
#'     (defaults `1e-4` and `0.9`).}
#'   \item{`qn_polish_alpha_max`}{Initial step length tried by the line
#'     search (default `1.0`).}
#'   \item{`qn_polish_max_ls_iter`}{Maximum line-search evaluations per
#'     iteration (default `20`).}
#'   \item{`qn_polish_max_fail`}{Consecutive line-search failures that
#'     trigger a revert to cubic mode (default `3`).}
#'   \item{`qn_polish_reenter_delay`}{Cubic iterations required after a
#'     revert before qn_polish may re-fire (default `5`).}
#'   \item{`qn_polish_curv_eps`}{Curvature threshold for skipping BFGS
#'     updates to preserve PD (default `1e-10`).}
#' }
#'
#' @section Quasi-Newton variant (`use_qn = TRUE`):
#' When `control$use_qn = TRUE`, `arcopt()` routes to an internal
#' quasi-Newton variant that approximates `H_k` via SR1/BFGS updates and
#' does not require a `hess` function. The following parameters apply
#' only when `use_qn = TRUE`.
#' \describe{
#'   \item{`qn_method`}{One of `"hybrid"` (default), `"sr1"`, `"bfgs"`.
#'     `"hybrid"` uses state-aware routing between SR1-first and
#'     BFGS-first orderings based on the current B's eigenstructure and
#'     recent rho values.}
#'   \item{`bfgs_tol`}{Curvature tolerance for the BFGS update
#'     (default `1e-10`).}
#'   \item{`sr1_skip_tol`}{SR1 skip-test tolerance (default `1e-8`).}
#'   \item{`sr1_restart_threshold`}{Consecutive SR1 skips before restart
#'     (default `5`).}
#'   \item{`qn_route_demote_rho`}{`rho` below this counts as a "bad"
#'     step in the routing FSM (default `0.25`).}
#'   \item{`qn_route_promote_rho`}{`rho` above this counts as a "good"
#'     step (default `0.5`).}
#'   \item{`qn_route_demote_k`}{Consecutive bad steps in `"pd"` routing
#'     mode that demote back to `"indefinite"` (default `2`).}
#'   \item{`qn_route_promote_k`}{Consecutive good PD steps in
#'     `"indefinite"` mode that promote to `"pd"` (default `3`).}
#'   \item{`qn_fd_refresh_k`}{Consecutive bad-rho iterations in
#'     `"indefinite"` mode that trigger an FD-Hessian refresh of B
#'     (default `3`).}
#'   \item{`qn_stuck_refresh_k`}{Iterations stuck in `"indefinite"` mode
#'     without promotion that force a refresh (default `100`).}
#'   \item{`use_accel_qn`}{**EXPERIMENTAL.** Enable Nesterov acceleration
#'     in the QN path. May improve convergence on strongly convex
#'     problems but can hurt nonconvex (default `FALSE`).}
#' }
#'
#' @section Diagnostics in `result$diagnostics`:
#' Mode-dispatch diagnostics are nested under `result$diagnostics` so the
#' primary return list stays compact.
#' \describe{
#'   \item{`solver_mode_final`}{`"cubic"`, `"tr"`, or `"qn_polish"` —
#'     which subproblem solver was active at termination.}
#'   \item{`ridge_switches`}{Integer count of cubic to TR transitions
#'     (`0` or `1` in v1; the switch is one-way).}
#'   \item{`radius_final`}{Final trust-region radius (`NA` if the solver
#'     never switched to TR mode).}
#'   \item{`qn_polish_switches`}{Integer count of cubic to qn_polish
#'     transitions (bidirectional; may be `> 1`).}
#'   \item{`qn_polish_reverts`}{Integer count of qn_polish to cubic
#'     reversions.}
#'   \item{`hess_evals_at_polish_switch`}{`evaluations$hess` at the first
#'     polish switch; compare against final `evaluations$hess` to
#'     quantify Hessian-evaluation savings.}
#' }
#' QN-variant runs add `qn_updates`, `qn_skips`, `qn_restarts`, and
#' `qn_fd_refreshes` to the same sublist.
#'
#' @seealso \code{\link{arcopt}} for the user-facing entry point.
#'
#' @name arcopt_advanced_controls
NULL
