# Quasi-Newton ARC Optimizer (Algorithm 4)
# =========================================
#
# QN-ARC: Uses quasi-Newton Hessian approximations instead of exact Hessians.
# Supports SR1, BFGS, L-BFGS, and L-SR1 update methods.


#' Quasi-Newton ARC Optimizer
#'
#' Minimizes a nonlinear objective function using Adaptive Regularization with
#' Cubics and quasi-Newton Hessian approximations. Does not require an exact
#' Hessian function.
#'
#' @param x0 Numeric vector of initial parameter values (length n).
#' @param fn Function that computes the objective function value.
#' @param gr Function that computes the gradient.
#' @param hess Optional Hessian function for hybrid mode. If provided,
#'   initializes B_0 = H(x_0) and can refresh when approximation degrades.
#' @param ... Additional arguments passed to `fn`, `gr`, and `hess`.
#' @param lower Numeric vector of lower bounds (length n). Default: all -Inf.
#' @param upper Numeric vector of upper bounds (length n). Default: all Inf.
#' @param control List of control parameters (see Details).
#'
#' @details
#' The QN-ARC algorithm maintains a quasi-Newton Hessian approximation B_k
#' and solves the cubic regularization subproblem:
#' \deqn{m_k(s) = f_k + g_k^T s + 1/2 s^T B_k s + sigma_k/3 ||s||^3}
#'
#' ## Initialization of B_0
#' When `hess` is not supplied, `arcopt_qn()` seeds the initial Hessian
#' approximation with a one-time finite-difference Hessian computed from
#' the supplied gradient at `x0` (cost: `2 * length(x0)` gradient
#' evaluations once at startup). This provides negative-curvature
#' information that the pure-identity initialization used by classical
#' quasi-Newton methods would miss on saddle-prone problems. To opt out
#' and recover the identity-initialization convention, pass
#' `hess = function(x) diag(length(x))` explicitly.
#'
#' ## Control Parameters
#' In addition to standard arcopt controls:
#' * `qn_method`: Update method (default: "hybrid"):
#'   - "hybrid": State-aware routing between SR1-first and BFGS-first
#'     orderings based on the current B_k's eigenstructure and recent
#'     rho_k values (see below). Includes automatic FD refresh of B_k
#'     when SR1 is stuck. Recommended for saddle-prone problems.
#'   - "sr1": Symmetric Rank-1 (allows indefinite Hessians)
#'   - "bfgs": Standard BFGS (maintains positive definiteness)
#' * `bfgs_tol`: Curvature tolerance for BFGS in hybrid mode (default: 1e-10)
#' * `sr1_skip_tol`: SR1 skip test tolerance (default: 1e-8)
#' * `sr1_restart_threshold`: Consecutive skips before restart (default: 5)
#' * `use_accel_qn`: **EXPERIMENTAL** Enable Nesterov acceleration (default:
#'   FALSE). May improve convergence on strongly convex problems but can hurt
#'   performance on nonconvex problems. Use with caution.
#'
#' ## Hybrid Routing Parameters
#' The `"hybrid"` method maintains an internal mode flag that starts in
#' `"indefinite"` (SR1-first priority) and promotes to `"pd"` (BFGS-first)
#' once the B_k approximation has been reliably positive definite and
#' cubic-model predictions have tracked the true objective:
#' * `qn_route_demote_rho` (default 0.25): rho_k below this counts as a
#'   "bad" step
#' * `qn_route_promote_rho` (default 0.5): rho_k above this counts as a
#'   "good" step
#' * `qn_route_demote_k` (default 2): consecutive bad steps in "pd" mode
#'   demote back to "indefinite"
#' * `qn_route_promote_k` (default 3): consecutive good PD steps in
#'   "indefinite" mode promote to "pd"
#' * `qn_fd_refresh_k` (default 3): while in "indefinite" mode, this many
#'   consecutive bad rho steps rebuild B_k from a fresh FD Hessian
#' * `qn_stuck_refresh_k` (default 100): while in "indefinite" mode,
#'   this many iterations without promoting also triggers an FD refresh
#'   (safety net for secondary-saddle stalls)
#'
#' @return Same structure as arcopt, plus:
#' * `qn_updates`: Number of successful QN updates
#' * `qn_skips`: Number of skipped updates
#' * `qn_restarts`: Number of approximation restarts
#' * `qn_fd_refreshes`: Number of FD Hessian refreshes performed (hybrid
#'   mode only; 0 for other qn_method values)
#'
#' @keywords internal
arcopt_qn <- function(x0, fn, gr, hess = NULL, ...,
                      lower = rep(-Inf, length(x0)),
                      upper = rep(Inf, length(x0)),
                      control = list()) {
  # Input validation
 if (!is.numeric(x0) || any(!is.finite(x0))) {
    stop("x0 must be a finite numeric vector")
  }

  n <- length(x0)

  if (!is.function(fn)) {
    stop("fn must be a function")
  }
  if (!is.function(gr)) {
    stop("gr (gradient function) is required")
  }
  if (length(lower) != n || length(upper) != n) {
    stop("lower and upper must have same length as x0")
  }
  if (any(lower >= upper)) {
    stop("lower must be strictly less than upper for all parameters")
  }

  # Default control parameters
  control_defaults <- list(
    maxit = 1000,
    ftol_abs = 1e-8,
    gtol_abs = 1e-5,
    xtol_abs = 1e-8,
    sigma0 = 1.0,
    sigma_min = 1e-6,
    sigma_max = 1e12,
    eta1 = 0.1,
    eta2 = 0.9,
    gamma1 = 0.5,
    gamma2 = 2.0,
    # QN-specific parameters
    qn_method = "hybrid",
    bfgs_tol = 1e-10,
    sr1_skip_tol = 1e-8,
    sr1_restart_threshold = 5L,
    # State-aware hybrid routing: switches SR1-first when B is indefinite
    # or when BFGS predictions drift; switches BFGS-first when B has been
    # positive definite and predictions have tracked well.
    qn_route_demote_rho = 0.25,   # rho below this counts as a bad step
    qn_route_promote_rho = 0.5,   # rho above this counts as a good step
    qn_route_demote_k = 2L,       # consecutive bad steps -> indefinite mode
    qn_route_promote_k = 3L,      # consecutive good PD steps -> pd mode
    # FD refresh: if bad rho persists while already in indefinite mode,
    # recompute B from an FD Hessian at the current iterate (2*n gradient
    # evals). Handles cases where the SR1 approximation has drifted too
    # far from the true Hessian to recover via secant updates alone.
    qn_fd_refresh_k = 3L,         # consecutive bad rho (in indefinite mode) -> refresh
    # Safety-net refresh: if the iteration stays in indefinite mode for
    # this many iterations without promoting to pd, force a refresh.
    # Catches "stuck at secondary saddle" cases where rho is OK per-step
    # but the iterate is not making global progress to a local minimum.
    # Deliberately large so that SR1 has time to accumulate secant
    # information before we discard it.
    qn_stuck_refresh_k = 100L,    # iterations stuck in indefinite mode -> refresh
    # Nesterov acceleration (Algorithm 4b) - EXPERIMENTAL, disabled by default
    # Can hurt convergence on nonconvex problems; needs adaptive restart logic
    use_accel_qn = FALSE,
    trace = 1
  )

  control <- modifyList(control_defaults, control)

  # Validate QN method
  valid_methods <- c("hybrid", "sr1", "bfgs")
  if (!control$qn_method %in% valid_methods) {
    stop("qn_method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  # Validate and project initial point to bounds
  x_current <- validate_and_project_initial(x0, lower, upper)

  # Initialize tracking variables
  sigma_current <- control$sigma0
  iter <- 0
  fn_evals <- 0
  gr_evals <- 0
  hess_evals <- 0
  qn_updates <- 0
  qn_skips <- 0
  qn_restarts <- 0
  skip_count <- 0L

  # State-aware hybrid routing. Initialized to "indefinite" because B_0
  # is the FD Hessian (or user-supplied hess), which can have negative
  # eigenvalues on saddle-prone problems; SR1-first is the safer default
  # until we confirm B is reliably positive definite and accurate.
  routing_mode <- "indefinite"
  bad_rho_count <- 0L
  pd_count <- 0L
  # Counter for FD refresh trigger: tracks consecutive bad rho while
  # we are already in indefinite mode. Separate from bad_rho_count so
  # that a pd -> indefinite demotion does not immediately arm the
  # refresh.
  indef_bad_rho_count <- 0L
  # Counter for safety-net refresh: tracks iterations spent in
  # indefinite mode. Reset on promotion to pd or on any refresh.
  indef_iter_count <- 0L
  qn_fd_refreshes <- 0L

  # Initial evaluation
  f_current <- fn(x_current, ...)
  g_current <- gr(x_current, ...)
  fn_evals <- fn_evals + 1
  gr_evals <- gr_evals + 1

  if (!check_finite(f_current, g_current)) {
    stop("Initial point yields NaN or Inf in function or gradient")
  }

  # Finite-difference Hessian from the gradient -- used to seed B_0 when
  # an exact Hessian is not supplied. Costs 2*n gradient evaluations once
  # at startup. Seeding from a real Hessian (rather than the identity) is
  # essential on saddle-prone problems: it injects negative-curvature
  # information into B_0 immediately, so the cubic subproblem can step
  # off the symmetric axis in iteration 1.
  fd_hess_from_grad <- function(theta, h = 1e-5) {
    d <- length(theta)
    fd_mat <- matrix(0, d, d)
    for (i in seq_len(d)) {
      ei <- rep(0, d)
      ei[i] <- h
      fd_mat[, i] <- (gr(theta + ei, ...) - gr(theta - ei, ...)) / (2 * h)
    }
    0.5 * (fd_mat + t(fd_mat))
  }

  # Initialize full-matrix Hessian approximation (B is n x n)
  if (!is.null(hess)) {
    # Exact Hessian: initialize from H(x_0)
    b_current <- hess(x_current, ...)
    hess_evals <- hess_evals + 1
  } else {
    # Default: FD Hessian from the gradient (one-time at x_0)
    b_current <- fd_hess_from_grad(x_current)
    gr_evals <- gr_evals + 2 * n
  }

  # Store previous values
  x_previous <- NULL
  f_previous <- NA
  g_previous <- NULL

  converged <- FALSE
  conv_reason <- ""

  # Stagnation tracking
  step_norms <- numeric(0)
  f_values <- c(f_current)

  # Nesterov acceleration variables (Algorithm 4b)
  use_accel <- control$use_accel_qn
  if (use_accel) {
    v_current <- x_current  # Auxiliary sequence
    a_current <- 1.0        # Acceleration parameter a_k
    a_sum <- 1.0            # Cumulative sum A_k = sum of a_j
  }

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

    # STEP 1b: Nesterov acceleration - compute interpolation point (Algorithm 4b)
    if (use_accel) {
      # Compute acceleration parameters: a_{k+1}^2 = A_k + a_{k+1}
      a_next <- (1 + sqrt(1 + 4 * a_sum)) / 2
      tau_k <- a_next / (a_sum + a_next)

      # Interpolation point: y_k = tau_k * v_k + (1 - tau_k) * x_k
      y_current <- tau_k * v_current + (1 - tau_k) * x_current

      # Use y_current for gradient evaluation
      g_eval <- gr(y_current, ...)
      gr_evals <- gr_evals + 1
    } else {
      y_current <- x_current
      g_eval <- g_current
    }

    # STEP 2: Solve cubic subproblem via eigendecomposition on the
    # full-matrix approximation. (Use g_eval under acceleration,
    # g_current otherwise.)
    g_subproblem <- g_eval
    cubic_result <- solve_cubic_subproblem_dispatch(
      g = g_subproblem,
      H = b_current,
      sigma = sigma_current,
      solver = "eigen"
    )
    s_current <- cubic_result$s

    # Compute predicted reduction
    b_s <- as.vector(b_current %*% s_current)

    s_norm <- sqrt(sum(s_current^2))
    pred_reduction <- -(
      sum(g_subproblem * s_current) +
        0.5 * sum(s_current * b_s) +
        (sigma_current / 3) * s_norm^3
    )

    # Apply box constraints (from y_current for acceleration, x_current otherwise)
    box_result <- apply_box_constraints(y_current, s_current, lower, upper)
    s_current <- box_result$s_bounded
    x_trial <- box_result$x_new

    # Evaluate trial point
    f_trial <- fn(x_trial, ...)
    fn_evals <- fn_evals + 1

    if (!check_finite(f_trial, rep(0, n))) {
      # NaN/Inf: increase sigma and reject
      sigma_current <- min(10 * sigma_current, control$sigma_max)
      iter <- iter + 1
      next
    }

    # Compute actual reduction and ratio
    actual_reduction <- f_current - f_trial

    # Guard against division issues
    if (abs(pred_reduction) < .Machine$double.eps) {
      rho <- if (actual_reduction >= 0) 1.0 else 0.0
    } else {
      rho <- actual_reduction / pred_reduction
    }

    # Accept or reject step
    if (is.finite(rho) && rho >= control$eta1) {
      # Accept step
      x_previous <- x_current
      f_previous <- f_current
      g_previous <- g_current

      x_current <- x_trial
      f_current <- f_trial
      g_current <- gr(x_current, ...)
      gr_evals <- gr_evals + 1

      if (!check_finite(f_current, g_current)) {
        stop("NaN or Inf detected in gradient at accepted point")
      }

      # Update Nesterov auxiliary sequence (Algorithm 4b)
      if (use_accel) {
        # v_{k+1} = x_{k+1} + (a_k - 1) / a_{k+1} * (x_{k+1} - x_{previous})
        momentum_coef <- (a_current - 1) / a_next
        v_current <- x_current + momentum_coef * (x_current - x_previous)
        a_current <- a_next
        a_sum <- a_sum + a_next
      }

      # STEP 3: Update QN approximation
      s_vec <- x_current - x_previous
      y_vec <- g_current - g_previous

      # Full-matrix QN update
      if (control$qn_method == "hybrid") {
        # State-aware routing: choose update priority based on whether
        # B_k is currently indefinite or whether BFGS predictions have
        # been drifting. See qn_route_* control parameters.
        hybrid_fn <- if (routing_mode == "indefinite") {
          update_hybrid_sr1_first
        } else {
          update_hybrid
        }
        update_result <- hybrid_fn(
          b_current, s_vec, y_vec,
          bfgs_tol = control$bfgs_tol,
          sr1_skip_tol = control$sr1_skip_tol,
          skip_count = skip_count,
          restart_threshold = control$sr1_restart_threshold
        )
        b_current <- update_result$b
        skip_count <- update_result$skip_count

        # Update routing state based on lambda_min(B) and recent rho.
        # lambda_min is cheap to obtain via Cholesky inertia: if chol()
        # succeeds, B is positive definite; otherwise at least one
        # eigenvalue is non-positive.
        b_is_pd <- !inherits(
          try(chol(b_current), silent = TRUE), "try-error"
        )

        # rho tracking: count consecutive "bad" predictions
        if (is.finite(rho) && rho < control$qn_route_demote_rho) {
          bad_rho_count <- bad_rho_count + 1L
        } else {
          bad_rho_count <- 0L
        }

        # Promote pd_count only when B is PD AND current step's rho
        # looks healthy
        if (b_is_pd && is.finite(rho) && rho > control$qn_route_promote_rho) {
          pd_count <- pd_count + 1L
        } else {
          pd_count <- 0L
        }

        # Mode transitions
        if (!b_is_pd) {
          # Always enter indefinite mode when B has non-PD directions
          routing_mode <- "indefinite"
        } else if (routing_mode == "pd" &&
                   bad_rho_count >= control$qn_route_demote_k) {
          # Drift in BFGS predictions: fall back to SR1-first
          routing_mode <- "indefinite"
          bad_rho_count <- 0L
        } else if (routing_mode == "indefinite" &&
                   pd_count >= control$qn_route_promote_k) {
          # B has been reliably PD and predictions good: promote
          routing_mode <- "pd"
          pd_count <- 0L
        }

        # FD refresh: two complementary triggers while routing is in
        # indefinite mode.
        #   (a) Drift: rho < demote_rho for K consecutive steps means
        #       the current B's predictions don't match the objective.
        #   (b) Stall: indefinite mode persists for qn_stuck_refresh_k
        #       iterations without promoting to pd, which indicates
        #       the iteration is tracking a persistent saddle (local
        #       rho may be acceptable but global progress has stalled).
        # Either trigger rebuilds B from a fresh FD Hessian at the
        # current iterate (2*n gradient evals).
        if (routing_mode == "indefinite") {
          indef_iter_count <- indef_iter_count + 1L
          if (is.finite(rho) && rho < control$qn_route_demote_rho) {
            indef_bad_rho_count <- indef_bad_rho_count + 1L
          } else {
            indef_bad_rho_count <- 0L
          }
          trigger_drift <- indef_bad_rho_count >= control$qn_fd_refresh_k
          trigger_stall <- indef_iter_count >= control$qn_stuck_refresh_k
          if (trigger_drift || trigger_stall) {
            b_current <- fd_hess_from_grad(x_current)
            gr_evals <- gr_evals + 2 * n
            qn_fd_refreshes <- qn_fd_refreshes + 1L
            indef_bad_rho_count <- 0L
            indef_iter_count <- 0L
            bad_rho_count <- 0L
            pd_count <- 0L
            skip_count <- 0L
          }
        } else {
          indef_bad_rho_count <- 0L
          indef_iter_count <- 0L
        }

        if (update_result$restarted) {
          qn_restarts <- qn_restarts + 1
        }
        if (update_result$skipped) {
          qn_skips <- qn_skips + 1
        } else {
          qn_updates <- qn_updates + 1
        }
      } else if (control$qn_method == "sr1") {
        update_result <- update_sr1(
          b_current, s_vec, y_vec,
          skip_tol = control$sr1_skip_tol,
          skip_count = skip_count,
          restart_threshold = control$sr1_restart_threshold
        )
        b_current <- update_result$b
        skip_count <- update_result$skip_count

        if (update_result$restarted) {
          qn_restarts <- qn_restarts + 1
        }
        if (update_result$skipped) {
          qn_skips <- qn_skips + 1
        } else {
          qn_updates <- qn_updates + 1
        }
      } else {
        # BFGS update
        update_result <- update_bfgs(b_current, s_vec, y_vec)
        b_current <- update_result$b

        if (update_result$skipped) {
          qn_skips <- qn_skips + 1
        } else {
          qn_updates <- qn_updates + 1
        }
      }

      # Track for stagnation detection
      step_norms <- c(step_norms, sqrt(sum(s_vec^2)))
      f_values <- c(f_values, f_current)
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

    iter <- iter + 1

    # Print progress if tracing
    if (control$trace >= 1 && iter %% 10 == 0) {
      g_norm <- sqrt(sum(g_current^2))
      message(sprintf(
        "Iter %4d: f = %.6e, |g| = %.6e, sigma = %.2e",
        iter, f_current, g_norm, sigma_current
      ))
    }
  }

  # Check if max iterations reached
  if (iter >= control$maxit && !converged) {
    conv_reason <- "maxit"
  }

  # Build final Hessian approximation for output
  h_final <- b_current

  list(
    par = x_current,
    value = f_current,
    gradient = g_current,
    hessian = h_final,
    converged = converged,
    iterations = iter,
    evaluations = list(
      fn = fn_evals,
      gr = gr_evals,
      hess = hess_evals
    ),
    message = conv_reason,
    qn_updates = qn_updates,
    qn_skips = qn_skips,
    qn_restarts = qn_restarts,
    qn_fd_refreshes = qn_fd_refreshes
  )
}
