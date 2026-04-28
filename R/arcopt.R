#' Adaptive Regularization using Cubics Optimizer
#'
#' Minimizes a nonlinear objective function using Adaptive Regularization with
#' Cubics (ARC). Designed for robust optimization of ill-conditioned, nonconvex,
#' and indefinite Hessian problems common in statistical applications.
#'
#' @param x0 Numeric vector of initial parameter values (length Q).
#' @param fn Function that computes the objective function value. Should take
#'   a numeric vector of length Q and return a scalar.
#' @param gr Function that computes the gradient. Should take a numeric
#'   vector of length Q and return a numeric vector of length Q. Required.
#' @param hess Function that computes the Hessian matrix. Should take a
#'   numeric vector of length Q and return a Q-by-Q symmetric matrix.
#'   Required (unless `control$use_qn = TRUE`; see Details).
#' @param ... Additional arguments passed to `fn`, `gr`, and `hess`.
#' @param lower Numeric vector of lower bounds (length Q). Use `-Inf` for
#'   unbounded parameters. Default: all `-Inf`.
#' @param upper Numeric vector of upper bounds (length Q). Use `Inf` for
#'   unbounded parameters. Default: all `Inf`.
#' @param control A named list of control parameters. The user-facing
#'   tolerances and switches are documented below; advanced regularization
#'   tuning, the trust-region fallback, and the quasi-Newton polish mode
#'   live on a separate help page (see `\link{arcopt_advanced_controls}`).
#'   Recognized entries:
#'   \describe{
#'     \item{`maxit`}{Maximum number of iterations (default `1000`).}
#'     \item{`gtol_abs`}{Absolute gradient-norm tolerance for convergence
#'       (default `1e-5`).}
#'     \item{`ftol_abs`}{Absolute objective-value tolerance (default
#'       `1e-8`).}
#'     \item{`xtol_abs`}{Absolute step-size tolerance (default `1e-8`).}
#'     \item{`trace`}{Integer in `0:3`. Depth of per-iteration data captured
#'       into `result$trace`: `0` collects nothing, `1` (the default)
#'       collects function value and gradient norm, `2` adds sigma, rho,
#'       step type and reciprocal Hessian condition number, `3` adds the
#'       full iterate, step, Hessian and convergence-criterion record.
#'       This flag controls *saved* data only -- for live console output
#'       see `verbose`.}
#'     \item{`verbose`}{Logical. If `TRUE`, prints one line per iteration
#'       to the console showing iteration number, objective value,
#'       `||g||_inf`, ratio rho, regularization scale, and active solver
#'       mode (default `FALSE`). Orthogonal to `trace`.}
#'     \item{`use_qn`}{Logical. If `TRUE`, route to the quasi-Newton ARC
#'       variant, which approximates the Hessian via SR1/BFGS updates and
#'       does not require an analytic `hess` function. See the
#'       advanced-controls page for QN-specific parameters (default
#'       `FALSE`).}
#'   }
#'   See `\link{arcopt_advanced_controls}` for the full set of advanced
#'   tuning parameters governing the cubic regularization, the trust-region
#'   fallback, and the quasi-Newton polish mode.
#'
#' @details
#' The ARC algorithm iteratively minimizes a cubic regularization model:
#' \deqn{m_k(s) = f_k + g_k^\top s + \frac{1}{2} s^\top H_k s +
#'   \frac{\sigma_k}{3} \|s\|^3}
#' where \eqn{\sigma_k} is adapted from observed model accuracy. arcopt
#' may transparently fall back to a trust-region subproblem in flat-ridge
#' regimes and (optionally, opt-in) to a line-search BFGS polish in the
#' quadratic attraction basin. The transitions are observable via
#' `result$diagnostics`; the algorithmic details and tunable thresholds
#' are documented under `\link{arcopt_advanced_controls}`.
#'
#' ## Hessian Requirement
#' arcopt is Hessian-centric: an analytic `hess` function is strongly
#' recommended. If the analytic form is unavailable, set
#' `control$use_qn = TRUE` to obtain Hessian-free quasi-Newton updates
#' (see the advanced-controls page).
#'
#' @return A list with components:
#' * `par`: Optimal parameter vector.
#' * `value`: Objective value at `par`.
#' * `gradient`: Gradient at `par`.
#' * `hessian`: Hessian at `par` (or the final BFGS approximation if the
#'   run ended in qn_polish mode).
#' * `sigma`: Final cubic regularization parameter.
#' * `converged`: Logical; whether convergence criteria were met.
#' * `iterations`: Number of iterations performed.
#' * `evaluations`: Named list of `fn`, `gr`, and `hess` evaluation counts.
#' * `message`: Convergence reason.
#' * `trace`: Per-iteration trace data (depth controlled by
#'   `control$trace`); `NULL` when `trace = 0`.
#' * `diagnostics`: Sublist of internal mode-dispatch diagnostics --
#'   `solver_mode_final`, `ridge_switches`, `radius_final`,
#'   `qn_polish_switches`, `qn_polish_reverts`, and
#'   `hess_evals_at_polish_switch`. See `\link{arcopt_advanced_controls}`
#'   for the meaning of each field. Most users do not need to inspect
#'   this; it is preserved for diagnostic and benchmarking use.
#'
#' @seealso \code{\link{arcopt_advanced_controls}} for advanced tuning of
#'   the cubic regularization, trust-region fallback, and quasi-Newton
#'   polish mode.
#'
#' @importFrom utils modifyList
#' @export
#' @examples
#' \donttest{
#' # Rosenbrock function
#' rosenbrock <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
#'
#' rosenbrock_gr <- function(x) {
#'   c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
#'     200 * (x[2] - x[1]^2))
#' }
#'
#' rosenbrock_hess <- function(x) {
#'   matrix(c(
#'     1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
#'     -400 * x[1], 200
#'   ), 2, 2)
#' }
#'
#' result <- arcopt(
#'   x0 = c(-1.2, 1),
#'   fn = rosenbrock,
#'   gr = rosenbrock_gr,
#'   hess = rosenbrock_hess
#' )
#'
#' print(result$par)      # Should be near c(1, 1)
#' print(result$value)    # Should be near 0
#' }
arcopt <- function(x0, fn, gr, hess = NULL, ...,
                   lower = rep(-Inf, length(x0)),
                   upper = rep(Inf, length(x0)),
                   control = list()) {
  # Input validation
  if (!is.numeric(x0) || any(!is.finite(x0))) {
    stop("x0 must be a finite numeric vector")
  }

  q <- length(x0)

  if (!is.function(fn)) {
    stop("fn must be a function")
  }

  if (!is.function(gr)) {
    stop("gr (gradient function) is required")
  }

  if (length(lower) != q || length(upper) != q) {
    stop("lower and upper must have same length as x0")
  }

  if (any(lower >= upper)) {
    stop("lower must be strictly less than upper for all parameters")
  }

  # Route to QN optimizer if requested
  use_qn <- if (!is.null(control$use_qn)) control$use_qn else FALSE
  if (use_qn) {
    return(arcopt_qn(x0, fn, gr, hess, ...,
                     lower = lower, upper = upper, control = control))
  }

  # Default control parameters. The user-facing "tier 1" entries
  # (tolerances, trace, verbose, use_qn) are documented in ?arcopt; the
  # advanced regularization, trust-region, and qn-polish entries are
  # documented under ?arcopt_advanced_controls.
  control_defaults <- list(
    # Tier 1: tolerances and switches
    maxit = 1000,
    gtol_abs = 1e-5,
    ftol_abs = 1e-8,
    xtol_abs = 1e-8,
    trace = 1,
    verbose = FALSE,
    # Tier 2: cubic-regularization tuning
    sigma0 = 1.0,
    sigma_min = 1e-6,
    sigma_max = 1e12,
    eta1 = 0.1,
    eta2 = 0.9,
    gamma1 = 0.5,
    gamma2 = 2.0,
    # Tier 2: trust-region fallback (cubic -> TR on flat-ridge detection)
    tr_fallback_enabled = TRUE,
    tr_fallback_window = 10,
    tr_fallback_tol_ridge = 1e-3,
    tr_fallback_rho_tol = 0.1,
    tr_fallback_grad_decrease_max = 0.9,
    tr_fallback_g_inf_floor = 1e-6,
    tr_r0 = 1.0,
    tr_rmax = 100,
    tr_eta1 = 0.25,
    tr_eta2 = 0.75,
    tr_gamma_shrink = 0.25,
    tr_gamma_grow = 2.0,
    # Tier 2: qn-polish fallback (cubic <-> line-search BFGS once inside
    # the quadratic attraction basin). Off by default in v0.2.0; enable
    # via control$qn_polish_enabled = TRUE for expensive-Hessian problems.
    qn_polish_enabled = FALSE,
    qn_polish_window = 5,
    qn_polish_rho = 0.9,
    qn_polish_lambda_min = 1e-3,
    qn_polish_g_decay = 0.5,
    qn_polish_g_inf_floor = 1e-8,
    qn_polish_c1 = 1e-4,
    qn_polish_c2 = 0.9,
    qn_polish_alpha_max = 1.0,
    qn_polish_max_ls_iter = 20,
    qn_polish_max_fail = 3,
    qn_polish_reenter_delay = 5,
    qn_polish_curv_eps = 1e-10
  )

  control <- modifyList(control_defaults, control)

  # Validate and project initial point to bounds
  x_current <- validate_and_project_initial(x0, lower, upper)

  # Initialize tracking variables
  sigma_current <- control$sigma0
  iter <- 0
  prev_rejected <- FALSE
  prev_used_newton <- FALSE  # Track if previous iteration used Newton
  fn_evals <- 0
  gr_evals <- 0
  hess_evals <- 0

  # Safeguards: Track history for stagnation detection
  step_norms <- numeric(0)
  f_values <- numeric(0)
  already_refreshed <- FALSE

  # Trust-region fallback state (one-way switch cubic -> tr)
  solver_mode <- "cubic"
  ridge_state <- init_flat_ridge_state(window = control$tr_fallback_window)
  radius_current <- NA_real_
  ridge_switches <- 0L
  s_current <- NULL # last step; used to seed TR radius at switch

  # QN-polish state (bidirectional cubic <-> qn_polish with cooldown)
  basin_state <- init_healthy_basin_state(window = control$qn_polish_window)
  b_current_polish <- NULL # BFGS approximation; non-NULL iff in qn_polish mode
  qn_polish_switches <- 0L
  qn_polish_reverts <- 0L
  qn_polish_fail_count <- 0L
  qn_polish_cooldown <- 0L
  hess_evals_at_polish_switch <- NA_integer_

  # Initial evaluation
  f_current <- fn(x_current, ...)
  g_current <- gr(x_current, ...)
  fn_evals <- fn_evals + 1
  gr_evals <- gr_evals + 1

  # Check for NaN/Inf at initial point
  if (!check_finite(f_current, g_current)) {
    stop("Initial point yields NaN or Inf in function or gradient")
  }

  # Validate Hessian specification
  if (is.null(hess)) {
    stop("Hessian function 'hess' must be provided")
  }

  # Initialize Hessian
  h_current <- hess(x_current, ...)
  hess_evals <- hess_evals + 1

  # Initialize history with initial values
  f_values <- c(f_values, f_current)

  # Store previous values for convergence checking
  f_previous <- NA
  x_previous <- NULL

  converged <- FALSE
  conv_reason <- ""

  # Initialize trace storage based on trace level
  trace_data <- NULL
  if (control$trace >= 1) {
    trace_data <- list(
      f = numeric(0),
      g_norm = numeric(0),
      converge_reason = character(0)
    )

    if (control$trace >= 2) {
      trace_data$sigma <- numeric(0)
      trace_data$rho <- numeric(0)
      trace_data$step_type <- character(0)
      trace_data$solver_used <- character(0)
      trace_data$rcond <- numeric(0)
      trace_data$lambda <- numeric(0)
      trace_data$time <- numeric(0)
    }

    if (control$trace >= 3) {
      n <- length(x0)
      trace_data$x <- matrix(NA_real_, nrow = 0, ncol = n)
      trace_data$s <- matrix(NA_real_, nrow = 0, ncol = n)
      trace_data$H <- array(NA_real_, dim = c(0, n, n))
      trace_data$converge_criteria <- list()
    }
  }

  # Initialize iteration timer
  iter_start_time <- NULL

  # Main optimization loop
  while (iter < control$maxit) {
    # STEP 1: Check convergence
    conv_result <- check_convergence(
      g_current = g_current,
      f_current = f_current,
      f_previous = f_previous,
      x_current = x_current,
      x_previous = x_previous,
      iter = iter,
      tol_grad = control$gtol_abs,
      tol_rel_grad = 1e-6,
      tol_obj = control$ftol_abs,
      tol_rel_obj = 1e-8,
      tol_param = control$xtol_abs,
      max_iter = control$maxit
    )

    if (conv_result$converged) {
      converged <- TRUE
      conv_reason <- conv_result$reason
      break
    }

    # STEP 1b: Check for stagnation
    stagnation_result <- detect_stagnation(
      step_norms = step_norms,
      f_values = f_values,
      max_stagnant = 5,
      is_qn = FALSE,  # Currently no QN support
      already_refreshed = already_refreshed
    )

    if (stagnation_result == "stop") {
      # Stagnation detected (no progress), but only declare convergence if gradient is small
      g_norm <- sqrt(sum(g_current^2))
      if (g_norm < control$gtol_abs) {
        converged <- TRUE
        conv_reason <- "stagnation"
        break
      } else {
        # Stagnation but gradient not small: unusual situation
        # Prevent infinite loop by not re-triggering stagnation detection
        already_refreshed <- TRUE
      }
    } else if (stagnation_result == "refresh_hessian") {
      # arcopt uses analytic or FD Hessians (not quasi-Newton)
      # Stagnation refresh mechanism reserved for future quasi-Newton support
      already_refreshed <- TRUE
    }

    # STEP 2 / 3 in cubic mode; TR branch in tr mode
    rho <- NA_real_ # will be set by whichever branch runs

    if (solver_mode == "tr") {
      # --- Trust-region fallback branch ---
      step_type <- "tr"
      tr_result <- solve_tr_eigen(g_current, h_current, radius_current)
      s_current <- tr_result$s
      pred_reduction <- tr_result$pred_reduction
      on_boundary <- isTRUE(tr_result$on_boundary)

      box_result <- apply_box_constraints(x_current, s_current, lower, upper)
      s_current <- box_result$s_bounded
      x_trial <- box_result$x_new

      f_trial <- fn(x_trial, ...)
      fn_evals <- fn_evals + 1

      if (!check_finite(f_trial, rep(0, length(x_trial)))) {
        # Non-finite trial: shrink radius and try again next iter
        radius_current <- max(radius_current * control$tr_gamma_shrink,
                              .Machine$double.eps)
        rho <- -Inf
      } else {
        actual_reduction <- f_current - f_trial
        rho <- actual_reduction / pred_reduction

        if (is.finite(rho) && rho >= control$tr_eta1) {
          x_previous <- x_current
          f_previous <- f_current
          x_current <- x_trial
          f_current <- f_trial
          g_current <- gr(x_current, ...)
          gr_evals <- gr_evals + 1
          if (!check_finite(f_current, g_current)) {
            stop("NaN or Inf detected in gradient at accepted TR point")
          }
          h_current <- hess(x_current, ...)
          hess_evals <- hess_evals + 1
          step_norms <- c(step_norms, sqrt(sum(s_current^2)))
          f_values <- c(f_values, f_current)
        }

        # Radius update
        if (!is.finite(rho) || rho < control$tr_eta1) {
          radius_current <- max(radius_current * control$tr_gamma_shrink,
                                .Machine$double.eps)
        } else if (rho >= control$tr_eta2 && on_boundary) {
          radius_current <- min(radius_current * control$tr_gamma_grow,
                                control$tr_rmax)
        }
      }
    }

    if (solver_mode == "qn_polish") {
      # --- QN-polish branch: line-search BFGS with warmstarted B_0 = H
      # (at switch). No further hess() evaluations unless we revert. ---
      step_type <- "qn_polish"
      used_newton <- FALSE

      # 1. Compute direction d = -B^{-1} g via Cholesky.
      chol_B <- tryCatch(chol(b_current_polish), error = function(e) NULL)
      if (is.null(chol_B)) {
        # B drifted non-PD; revert to cubic
        solver_mode <- "cubic"
        h_current <- hess(x_current, ...)
        hess_evals <- hess_evals + 1
        qn_polish_reverts <- qn_polish_reverts + 1L
        b_current_polish <- NULL
        qn_polish_fail_count <- 0L
        qn_polish_cooldown <- control$qn_polish_reenter_delay
      } else {
        d_polish <- as.vector(backsolve(
          chol_B, backsolve(chol_B, -g_current, transpose = TRUE)
        ))

        # Descent-direction sanity check (should hold with PD B)
        g_dot_d <- sum(g_current * d_polish)
        if (!is.finite(g_dot_d) || g_dot_d >= 0) {
          solver_mode <- "cubic"
          h_current <- hess(x_current, ...)
          hess_evals <- hess_evals + 1
          qn_polish_reverts <- qn_polish_reverts + 1L
          b_current_polish <- NULL
          qn_polish_fail_count <- 0L
          qn_polish_cooldown <- control$qn_polish_reenter_delay
        } else {
          # 2. Wolfe line search along d
          ls <- wolfe_line_search(
            fn = fn, gr = gr,
            x = x_current, d = d_polish,
            f_x = f_current, g_x = g_current,
            c1 = control$qn_polish_c1,
            c2 = control$qn_polish_c2,
            alpha_max = control$qn_polish_alpha_max,
            max_iter = control$qn_polish_max_ls_iter,
            ...
          )
          fn_evals <- fn_evals + ls$evals_f
          gr_evals <- gr_evals + ls$evals_g

          if (!ls$success) {
            # Line search failed; count, maybe revert
            qn_polish_fail_count <- qn_polish_fail_count + 1L
            if (qn_polish_fail_count >= control$qn_polish_max_fail) {
              solver_mode <- "cubic"
              h_current <- hess(x_current, ...)
              hess_evals <- hess_evals + 1
              qn_polish_reverts <- qn_polish_reverts + 1L
              b_current_polish <- NULL
              qn_polish_fail_count <- 0L
              qn_polish_cooldown <- control$qn_polish_reenter_delay
            }
            # rho meaningless on LS failure; use NaN for trace
            rho <- NaN
            pred_reduction <- 0
          } else {
            # Line search succeeded -- accept step
            qn_polish_fail_count <- 0L
            s_polish <- ls$alpha * d_polish

            # Apply box constraints (same convention as cubic path)
            box_result <- apply_box_constraints(x_current, s_polish,
                                                lower, upper)
            s_polish <- box_result$s_bounded
            x_trial <- box_result$x_new

            # If box truncation changed the step substantially, re-eval
            if (max(abs(x_trial - ls$x_new)) > .Machine$double.eps) {
              f_trial <- fn(x_trial, ...)
              fn_evals <- fn_evals + 1
              if (check_finite(f_trial, rep(0, length(x_trial)))) {
                g_trial <- gr(x_trial, ...)
                gr_evals <- gr_evals + 1
              } else {
                g_trial <- ls$g_new
                f_trial <- ls$f_new
              }
            } else {
              f_trial <- ls$f_new
              g_trial <- ls$g_new
            }

            # Track rho for trace and detector (quadratic model vs actual)
            actual_reduction <- f_current - f_trial
            h_s <- as.vector(b_current_polish %*% s_polish)
            pred_reduction <- -sum(g_current * s_polish) -
              0.5 * sum(s_polish * h_s)
            rho <- if (abs(pred_reduction) < .Machine$double.eps) {
              if (actual_reduction >= 0) 1.0 else 0.0
            } else {
              actual_reduction / pred_reduction
            }

            # Promote trial to current
            g_previous <- g_current
            x_previous <- x_current
            f_previous <- f_current
            x_current <- x_trial
            f_current <- f_trial
            g_current <- g_trial
            s_current <- s_polish
            if (!check_finite(f_current, g_current)) {
              stop("NaN/Inf in gradient at accepted qn_polish step")
            }

            # BFGS update: skip if curvature condition fails
            y_vec <- g_current - g_previous
            s_dot_y <- sum(s_polish * y_vec)
            curv_bound <- control$qn_polish_curv_eps *
              sqrt(sum(s_polish^2) * sum(y_vec^2))
            if (is.finite(s_dot_y) && s_dot_y > curv_bound) {
              update_result <- update_bfgs(b_current_polish,
                                           s_polish, y_vec)
              b_current_polish <- update_result$b
            }
            # else: skip update, preserves PD of B

            step_norms <- c(step_norms, sqrt(sum(s_polish^2)))
            f_values <- c(f_values, f_current)
          }
        }
      }
    }

    if (solver_mode == "cubic") {
    # STEP 2: Try Newton step first (if not rejected and hess provided)
    step_type <- "cubic"
    used_newton <- FALSE
    cubic_result <- NULL

    if (!prev_rejected && !is.null(hess)) {
      newton_result <- try_newton_step(g_current, h_current)

      # If previous iteration used Newton but current Cholesky fails (indefinite H),
      # reset sigma to initial value to provide adequate regularization
      if (!newton_result$success && prev_used_newton) {
        sigma_current <- control$sigma0
      }

      if (newton_result$success) {
        # Newton step succeeded (H is positive definite)
        s_current <- newton_result$s
        pred_reduction <- newton_result$pred_reduction

        # Apply box constraints to Newton step
        box_result <- apply_box_constraints(x_current, s_current, lower, upper)
        s_current <- box_result$s_bounded
        x_trial <- box_result$x_new

        # Evaluate trial point
        f_trial <- fn(x_trial, ...)
        fn_evals <- fn_evals + 1

        # Check for NaN/Inf in Newton trial evaluation
        if (!check_finite(f_trial, rep(0, length(x_trial)))) {
          # NaN/Inf detected: treat as unsuccessful, will use cubic
          prev_rejected <- TRUE
          used_newton <- FALSE
        } else {
          # Compute actual reduction and ratio
          actual_reduction <- f_current - f_trial
          rho <- actual_reduction / pred_reduction

          # Accept or reject based on ratio test (same as cubic)
          if (is.finite(rho) && rho >= control$eta1) {
            # Accept Newton step
            x_previous <- x_current
            f_previous <- f_current

            # Trial point before momentum
            y_new <- x_trial
            f_new <- f_trial
            g_new <- gr(y_new, ...)
            gr_evals <- gr_evals + 1

            # Check for NaN/Inf in new gradient
            if (!check_finite(f_new, g_new)) {
              stop("NaN or Inf detected in gradient at accepted Newton point")
            }

            x_current <- y_new
            f_current <- f_new
            g_current <- g_new

            h_current <- hess(x_current, ...)
            hess_evals <- hess_evals + 1

            # Track step norm and function value
            step_norms <- c(step_norms, sqrt(sum(s_current^2)))
            f_values <- c(f_values, f_current)

            # Update sigma using standard algorithm (SAME AS CUBIC)
            sigma_current <- update_sigma_cgt(
              sigma_current = sigma_current,
              rho = rho,
              eta1 = control$eta1,
              eta2 = control$eta2,
              gamma1 = control$gamma1,
              gamma2 = control$gamma2,
              sigma_min = control$sigma_min,
              sigma_max = control$sigma_max
            )

            prev_rejected <- FALSE
            used_newton <- TRUE
            step_type <- "newton"
          } else {
            # Reject Newton step - will use cubic solver
            prev_rejected <- TRUE
            used_newton <- FALSE
          }
        }
      }
    }

    # STEP 3: If Newton failed/skipped, use cubic solver
    if (!used_newton) {
      # Solve cubic subproblem via the (internal) dispatcher; "auto"
      # currently always selects the eigendecomposition solver.
      cubic_result <- solve_cubic_subproblem_dispatch(
        g = g_current,
        H = h_current,
        sigma = sigma_current
      )

      s_current <- cubic_result$s
      pred_reduction <- cubic_result$pred_reduction

      # Apply box constraints to cubic step
      box_result <- apply_box_constraints(x_current, s_current, lower, upper)
      s_current <- box_result$s_bounded
      x_trial <- box_result$x_new

      # Evaluate trial point
      f_trial <- fn(x_trial, ...)
      fn_evals <- fn_evals + 1

      # Check for NaN/Inf in trial evaluation
      if (!check_finite(f_trial, rep(0, length(x_trial)))) {
        # NaN/Inf detected: increase sigma and reject step
        sigma_current <- min(10 * sigma_current, 1e16)
        prev_rejected <- TRUE
        iter <- iter + 1
        next  # Skip to next iteration with larger sigma
      }

      # Compute actual reduction and ratio
      actual_reduction <- f_current - f_trial
      rho <- actual_reduction / pred_reduction

      # Accept or reject step (check for finite rho)
      if (is.finite(rho) && rho >= control$eta1) {
        # Accept step
        x_previous <- x_current
        f_previous <- f_current

        x_current <- x_trial
        f_current <- f_trial
        g_current <- gr(x_current, ...)
        gr_evals <- gr_evals + 1

        if (!check_finite(f_current, g_current)) {
          stop("NaN or Inf detected in gradient at accepted point")
        }

        h_current <- hess(x_current, ...)
        hess_evals <- hess_evals + 1

        step_norms <- c(step_norms, sqrt(sum(s_current^2)))
        f_values <- c(f_values, f_current)

        prev_rejected <- FALSE
      } else {
        prev_rejected <- TRUE
      }

      # Update sigma based on step quality
      sigma_current <- update_sigma_cgt(
        sigma_current = sigma_current,
        rho = rho,
        eta1 = control$eta1,
        eta2 = control$eta2,
        gamma1 = control$gamma1,
        gamma2 = control$gamma2,
        sigma_min = control$sigma_min,
        sigma_max = control$sigma_max
      )
    }

    # Update both end-of-cubic-iteration detectors. Reuse the cubic
    # solver's eigendecomposition (sorted increasing) if it ran this
    # iteration; otherwise fall back to a fresh `eigen()` call.
    if (is.finite(rho) &&
          (control$tr_fallback_enabled || control$qn_polish_enabled)) {
      lambda_min_current <- if (!is.null(cubic_result$eigenvalues)) {
        cubic_result$eigenvalues[1]
      } else {
        tryCatch(
          min(eigen(h_current, symmetric = TRUE,
                    only.values = TRUE)$values),
          error = function(e) NA_real_
        )
      }
      g_inf_current <- max(abs(g_current))
    } else {
      lambda_min_current <- NA_real_
      g_inf_current <- NA_real_
    }

    # Flat-ridge detector: switch to TR mode if all four signals fire
    if (control$tr_fallback_enabled && is.finite(rho)) {
      ridge_state <- update_flat_ridge_state(
        ridge_state,
        sigma = sigma_current,
        rho = rho,
        g_inf = g_inf_current,
        lambda_min = lambda_min_current
      )
      if (check_flat_ridge_trigger(
        ridge_state,
        sigma_min = control$sigma_min,
        tol_ridge = control$tr_fallback_tol_ridge,
        rho_tol = control$tr_fallback_rho_tol,
        grad_decrease_max = control$tr_fallback_grad_decrease_max,
        g_inf_floor = control$tr_fallback_g_inf_floor
      )) {
        solver_mode <- "tr"
        # Start TR with the control's `tr_r0` (default 1.0, matching
        # `trust::trust`'s `rinit`). The last cubic step norm is NOT a
        # reliable initializer because cubic often takes huge huge steps
        # (scale 1/sigma) right before the detector fires. TR then grows
        # or shrinks the radius adaptively based on the subsequent rhos.
        radius_current <- min(control$tr_r0, control$tr_rmax)
        ridge_switches <- ridge_switches + 1L
      }
    }

    # Healthy-basin detector: switch to qn_polish mode if all five signals
    # fire and we're not in cooldown after a recent revert.
    if (control$qn_polish_enabled && is.finite(rho) &&
          solver_mode == "cubic") {
      basin_state <- update_healthy_basin_state(
        basin_state,
        used_newton = used_newton,
        rho = rho,
        lambda_min = lambda_min_current,
        g_inf = g_inf_current
      )
      if (qn_polish_cooldown > 0L) {
        qn_polish_cooldown <- qn_polish_cooldown - 1L
      } else if (check_healthy_basin_trigger(
        basin_state,
        rho_polish = control$qn_polish_rho,
        lambda_min_polish = control$qn_polish_lambda_min,
        g_decay_polish = control$qn_polish_g_decay,
        g_inf_floor_polish = control$qn_polish_g_inf_floor
      )) {
        solver_mode <- "qn_polish"
        b_current_polish <- h_current # warmstart B_0 = current Hessian
        qn_polish_switches <- qn_polish_switches + 1L
        hess_evals_at_polish_switch <- hess_evals
      }
    }
    } # end if (solver_mode == "cubic")

    # Capture trace data for this iteration
    if (control$trace >= 1) {
      iter_end_time <- Sys.time()

      # Level 1
      trace_data$f <- c(trace_data$f, f_current)
      trace_data$g_norm <- c(trace_data$g_norm, max(abs(g_current)))
      trace_data$converge_reason <- c(trace_data$converge_reason,
                                       if (iter == 0) "init" else "continuing")

      if (control$trace >= 2) {
        # Level 2
        trace_data$sigma <- c(trace_data$sigma, sigma_current)
        trace_data$rho <- c(trace_data$rho, if (iter == 0) NA_real_ else rho)
        trace_data$step_type <- c(trace_data$step_type,
                                   if (iter == 0) "init" else step_type)
        trace_data$solver_used <- c(trace_data$solver_used,
                                     if (iter == 0 || step_type == "newton")
                                       "cholesky" else "eigen")

        # Compute reciprocal condition number
        if (!is.null(h_current)) {
          eig_vals <- eigen(h_current, symmetric = TRUE, only.values = TRUE)$values
          rcond_val <- if (all(eig_vals > 0)) {
            min(eig_vals) / max(eig_vals)
          } else {
            0  # Indefinite
          }
        } else {
          rcond_val <- NA_real_
        }
        trace_data$rcond <- c(trace_data$rcond, rcond_val)

        # Lambda = sigma * ||s||
        lambda_val <- if (iter == 0) NA_real_ else {
          sigma_current * sqrt(sum(s_current^2))
        }
        trace_data$lambda <- c(trace_data$lambda, lambda_val)

        # Iteration time
        iter_time <- if (is.null(iter_start_time)) 0 else {
          as.numeric(iter_end_time - iter_start_time, units = "secs")
        }
        trace_data$time <- c(trace_data$time, iter_time)
      }

      if (control$trace >= 3) {
        # Level 3
        trace_data$x <- rbind(trace_data$x, x_current)
        trace_data$s <- rbind(trace_data$s,
                              if (iter == 0) rep(NA_real_, length(x_current)) else s_current)

        # Hessian (3D array binding)
        if (is.null(h_current)) {
          h_to_store <- matrix(NA_real_, length(x_current), length(x_current))
        } else {
          h_to_store <- h_current
        }
        trace_data$H <- array(c(trace_data$H, h_to_store),
                              dim = c(dim(trace_data$H)[1] + 1, length(x_current), length(x_current)))

        # Convergence criteria (all 6 checks)
        trace_data$converge_criteria[[iter + 1]] <- list(
          gtol_abs = max(abs(g_current)),
          gtol_rel = max(abs(g_current * x_current / max(1, abs(f_current)))),
          gtol_scaled = max(abs(g_current / pmax(abs(f_current), 1))),
          xtol_abs = if (iter == 0) NA_real_ else max(abs(s_current)),
          xtol_rel = if (iter == 0) NA_real_ else {
            max(abs(s_current / pmax(abs(x_current), 1)))
          },
          ftol_rel = if (iter == 0 || is.na(f_previous)) NA_real_ else {
            abs(f_current - f_previous) / pmax(abs(f_current), 1)
          }
        )
      }

      # Reset timer for next iteration
      iter_start_time <- Sys.time()
    }

    # Track if current iteration used Newton for next iteration's sigma reset logic
    prev_used_newton <- used_newton

    if (isTRUE(control$verbose)) {
      scale_label <- if (solver_mode == "tr") "radius" else "sigma"
      scale_value <- if (solver_mode == "tr") radius_current else sigma_current
      message(sprintf(
        "iter %4d  f = %12.6e  |g|_inf = %9.3e  rho = %6.3f  %s = %9.3e  mode = %s",
        iter + 1L, f_current, max(abs(g_current)),
        if (is.finite(rho)) rho else NA_real_,
        scale_label,
        if (is.null(scale_value) || !is.finite(scale_value)) NA_real_ else scale_value,
        step_type
      ))
    }

    iter <- iter + 1
  }

  # Check if we hit max iterations without converging
  if (!converged && iter >= control$maxit) {
    converged <- TRUE
    conv_reason <- "max_iter"
  }

  result <- list(
    par = x_current,
    value = f_current,
    gradient = g_current,
    hessian = if (solver_mode == "qn_polish") b_current_polish else h_current,
    sigma = sigma_current,
    converged = converged,
    iterations = iter,
    evaluations = list(fn = fn_evals, gr = gr_evals, hess = hess_evals),
    message = conv_reason,
    diagnostics = list(
      solver_mode_final = solver_mode,
      ridge_switches = ridge_switches,
      radius_final = radius_current,
      qn_polish_switches = qn_polish_switches,
      qn_polish_reverts = qn_polish_reverts,
      hess_evals_at_polish_switch = hess_evals_at_polish_switch
    )
  )

  # Add trace data if collected
  if (control$trace >= 1) {
    result$trace <- trace_data
  }

  return(result)
}
