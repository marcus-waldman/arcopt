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
#' @param lower Numeric vector of lower bounds (length n). Default: all -Inf.
#' @param upper Numeric vector of upper bounds (length n). Default: all Inf.
#' @param control List of control parameters (see Details).
#'
#' @details
#' The QN-ARC algorithm maintains a quasi-Newton Hessian approximation B_k
#' and solves the cubic regularization subproblem:
#' \deqn{m_k(s) = f_k + g_k^T s + 1/2 s^T B_k s + sigma_k/3 ||s||^3}
#'
#' ## Control Parameters
#' In addition to standard arcopt controls:
#' * `qn_method`: Update method (default: "hybrid"):
#'   - "hybrid": Tries BFGS, falls back to SR1, then Powell-damped BFGS
#'   - "sr1": Symmetric Rank-1 (allows indefinite Hessians)
#'   - "bfgs": Standard BFGS (maintains positive definiteness)
#'   - "lbfgs": Limited-memory BFGS
#'   - "lsr1": Limited-memory SR1
#' * `bfgs_tol`: Curvature tolerance for BFGS in hybrid mode (default: 1e-10)
#' * `qn_memory`: History size for limited-memory methods (default: 10)
#' * `sr1_skip_tol`: SR1 skip test tolerance (default: 1e-8)
#' * `sr1_restart_threshold`: Consecutive skips before restart (default: 5)
#' * `use_accel_qn`: **EXPERIMENTAL** Enable Nesterov acceleration (default:
#'   FALSE). May improve convergence on strongly convex problems but can hurt
#'   performance on nonconvex problems. Use with caution.
#'
#' @return Same structure as arcopt, plus:
#' * `qn_updates`: Number of successful QN updates
#' * `qn_skips`: Number of skipped updates
#' * `qn_restarts`: Number of approximation restarts
#'
#' @keywords internal
arcopt_qn <- function(x0, fn, gr, hess = NULL,
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
    qn_memory = 10L,
    bfgs_tol = 1e-10,
    sr1_skip_tol = 1e-8,
    sr1_restart_threshold = 5L,
    # Nesterov acceleration (Algorithm 4b) - EXPERIMENTAL, disabled by default
    # Can hurt convergence on nonconvex problems; needs adaptive restart logic
    use_accel_qn = FALSE,
    trace = 1
  )

  control <- modifyList(control_defaults, control)

  # Validate QN method
  valid_methods <- c("hybrid", "sr1", "bfgs", "lbfgs", "lsr1")
  if (!control$qn_method %in% valid_methods) {
    stop("qn_method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  # Determine if using limited-memory method
  use_limited_memory <- control$qn_method %in% c("lbfgs", "lsr1")

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

  # Initial evaluation
  f_current <- fn(x_current)
  g_current <- gr(x_current)
  fn_evals <- fn_evals + 1
  gr_evals <- gr_evals + 1

  if (!check_finite(f_current, g_current)) {
    stop("Initial point yields NaN or Inf in function or gradient")
  }

  # Initialize Hessian approximation
  if (use_limited_memory) {
    # Limited-memory: use history structure
    qn_history <- NULL
    gamma <- 1.0  # Initial scaling

    # If Hessian provided, use it to set initial gamma
    if (!is.null(hess)) {
      h_init <- hess(x_current)
      hess_evals <- hess_evals + 1
      # Use average diagonal as gamma
      gamma <- mean(diag(h_init))
      if (gamma <= 0) gamma <- 1.0
    }
  } else {
    # Full matrix: B is n x n
    if (!is.null(hess)) {
      # Hybrid mode: initialize from exact Hessian
      b_current <- hess(x_current)
      hess_evals <- hess_evals + 1
    } else {
      # Pure QN: initialize with scaled identity
      b_current <- diag(n)
    }
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
      g_eval <- gr(y_current)
      gr_evals <- gr_evals + 1
    } else {
      y_current <- x_current
      g_eval <- g_current
    }

    # STEP 2: Solve cubic subproblem (use g_eval for acceleration, g_current otherwise)
    g_subproblem <- g_eval
    if (use_limited_memory) {
      # Build compact form and use Woodbury solver
      if (control$qn_method == "lsr1") {
        compact <- lsr1_compact_form(qn_history)
      } else {
        compact <- lbfgs_compact_form(qn_history)
      }

      gamma_current <- if (!is.null(qn_history)) qn_history$gamma else gamma
      # Ensure gamma is positive (required by Woodbury solver)
      if (gamma_current <= 0) gamma_current <- 1.0

      if (!is.null(compact)) {
        cubic_result <- solve_cubic_woodbury(
          g = g_subproblem,
          gamma = gamma_current,
          w = compact$w,
          c_mat = compact$c_mat,
          sigma = sigma_current
        )
      } else {
        # No history yet, use simple solver
        cubic_result <- solve_cubic_simple(
          g = g_subproblem,
          gamma = gamma_current,
          sigma = sigma_current
        )
      }

      s_current <- cubic_result$s
    } else {
      # Full QN: use eigendecomposition solver
      cubic_result <- solve_cubic_subproblem_dispatch(
        g = g_subproblem,
        H = b_current,
        sigma = sigma_current,
        solver = "eigen"
      )

      s_current <- cubic_result$s
    }

    # Compute predicted reduction
    if (use_limited_memory) {
      # For limited-memory, B*s computed via compact form
      if (!is.null(compact)) {
        # B*s = gamma*s + w' * c_mat * w * s
        w_s <- as.vector(compact$w %*% s_current)
        b_s <- gamma_current * s_current +
          as.vector(t(compact$w) %*% (compact$c_mat %*% w_s))
      } else {
        b_s <- gamma_current * s_current
      }
    } else {
      b_s <- as.vector(b_current %*% s_current)
    }

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
    f_trial <- fn(x_trial)
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
      g_current <- gr(x_current)
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

      if (use_limited_memory) {
        if (control$qn_method == "lsr1") {
          # L-SR1 update
          b_times_s <- if (!is.null(qn_history)) {
            lsr1_multiply(qn_history, s_vec)
          } else {
            gamma * s_vec
          }
          qn_history <- update_lsr1(
            qn_history, s_vec, y_vec,
            m = control$qn_memory,
            skip_tol = control$sr1_skip_tol,
            b_times_s = b_times_s
          )
          if (qn_history$skipped) {
            qn_skips <- qn_skips + 1
          } else {
            qn_updates <- qn_updates + 1
          }
        } else {
          # L-BFGS update
          qn_history <- update_lbfgs(
            qn_history, s_vec, y_vec,
            m = control$qn_memory
          )
          # L-BFGS always updates if curvature condition satisfied
          ys <- sum(y_vec * s_vec)
          if (ys > .Machine$double.eps) {
            qn_updates <- qn_updates + 1
          } else {
            qn_skips <- qn_skips + 1
          }
        }
      } else {
        # Full matrix update
        if (control$qn_method == "hybrid") {
          # Hybrid: BFGS -> SR1 -> Powell-damped BFGS
          update_result <- update_hybrid(
            b_current, s_vec, y_vec,
            bfgs_tol = control$bfgs_tol,
            sr1_skip_tol = control$sr1_skip_tol,
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
  if (use_limited_memory) {
    h_final <- NULL  # Cannot easily return full matrix
  } else {
    h_final <- b_current
  }

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
    qn_restarts = qn_restarts
  )
}
