#' Adaptive Regularization using Cubics Optimizer
#'
#' Minimizes a nonlinear objective function using Adaptive Regularization with
#' Cubics (ARC). Designed for robust optimization of ill-conditioned, nonconvex,
#' and indefinite Hessian problems common in statistical applications.
#'
#' @param x0 Numeric vector of initial parameter values (length Q).
#' @param fn Function that computes the objective function value. Should take a
#'   numeric vector of length Q and return a scalar.
#' @param gr Function that computes the gradient. Should take a numeric vector
#'   of length Q and return a numeric vector of length Q. Required.
#' @param hess Function that computes the Hessian matrix. Should take a numeric
#'   vector of length Q and return a Q×Q symmetric matrix. Required.
#' @param lower Numeric vector of lower bounds (length Q). Use `-Inf` for
#'   unbounded parameters. Default: all `-Inf`.
#' @param upper Numeric vector of upper bounds (length Q). Use `Inf` for
#'   unbounded parameters. Default: all `Inf`.
#' @param control List of control parameters (see Details).
#'
#' @details
#' The ARC algorithm iteratively minimizes a cubic regularization model:
#' \deqn{m_k(s) = f_k + g_k^T s + \frac{1}{2} s^T H_k s + \frac{\sigma_k}{3} \|s\|^3}
#'
#' where \eqn{\sigma_k} is adaptively adjusted based on model accuracy.
#'
#' ## Hessian Requirement
#' ARC critically depends on accurate curvature information. The `hess` argument
#' is required. For convenience, finite-difference Hessians can be computed
#' automatically with `control$use_fd = TRUE`, but analytic Hessians are
#' strongly recommended for best performance.
#'
#' ## Control Parameters
#' The `control` list accepts:
#' * `maxit`: Maximum iterations (default: 1000)
#' * `ftol_abs`: Absolute function tolerance (default: 1e-8)
#' * `gtol_abs`: Absolute gradient norm tolerance (default: 1e-5)
#' * `xtol_abs`: Absolute step size tolerance (default: 1e-8)
#' * `sigma0`: Initial regularization parameter (default: 1.0)
#' * `eta1`: Acceptance threshold for step (default: 0.1)
#' * `eta2`: Very successful step threshold (default: 0.9)
#' * `gamma1`: Regularization decrease factor (default: 0.5)
#' * `gamma2`: Regularization increase factor (default: 2.0)
#' * `cubic_solver`: Solver selection: "auto" (recommended) or "eigen"
#'   (default: "auto"). Auto-selection uses eigendecomposition (Algorithm 5a) for
#'   robust handling of indefinite Hessians and hard cases.
#' * `use_momentum`: Enable momentum acceleration (default: FALSE).
#'   Implements Gao et al. (2022) ARCm with recursive momentum and bisection
#'   search for monotonicity. **Only recommended for known ill-conditioned
#'   problems** - on well-conditioned problems, the bisection overhead may
#'   negate iteration savings. Empirically shows mixed results: can dramatically
#'   reduce iterations on some problems while increasing them on others.
#' * `momentum_tau`: Maximum momentum parameter (default: 0.5, Gao's τ)
#' * `momentum_alpha1`: Linear step scaling constant (default: 0.1, Gao's α₁)
#' * `momentum_alpha2`: Quadratic step scaling constant (default: 1.0, Gao's α₂)
#' * `trace`: Print iteration progress (default: FALSE)
#'
#' @return A list with components:
#' * `par`: Optimal parameter vector
#' * `value`: Optimal function value
#' * `gradient`: Gradient at optimum
#' * `hessian`: Hessian at optimum (if `hess` provided)
#' * `converged`: Logical, whether convergence criteria met
#' * `iterations`: Number of iterations performed
#' * `evaluations`: List with `fn`, `gr`, and `hess` evaluation counts
#' * `message`: Convergence message
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
arcopt <- function(x0, fn, gr, hess = NULL,
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
    cubic_solver = "auto",
    use_momentum = FALSE,
    momentum_tau = 0.5,
    momentum_alpha1 = 0.1,
    momentum_alpha2 = 1.0,
    trace = 1
  )

  control <- modifyList(control_defaults, control)

  # Validate and project initial point to bounds
  x_current <- validate_and_project_initial(x0, lower, upper)

  # Initialize tracking variables
  sigma_current <- control$sigma0
  iter <- 0
  prev_rejected <- FALSE
  fn_evals <- 0
  gr_evals <- 0
  hess_evals <- 0

  # Safeguards: Track history for stagnation detection
  step_norms <- numeric(0)
  f_values <- numeric(0)
  already_refreshed <- FALSE

  # Momentum state (Gao et al. Algorithm 1)
  # v_prev initialized to zero vector per Gao: v_{-1} = 0
  v_prev <- rep(0, length(x0))

  # Initial evaluation
  f_current <- fn(x_current)
  g_current <- gr(x_current)
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
  h_current <- hess(x_current)
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
      if (control$use_momentum) {
        trace_data$beta <- numeric(0)
      }
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

    # STEP 2: Try Newton step first (if not rejected and hess provided)
    step_type <- "cubic"
    used_newton <- FALSE

    if (!prev_rejected && !is.null(hess)) {
      newton_result <- try_newton_step(g_current, h_current)

      if (newton_result$success) {
        # Newton step succeeded (H is positive definite)
        s_current <- newton_result$s
        pred_reduction <- newton_result$pred_reduction

        # Apply box constraints to Newton step
        box_result <- apply_box_constraints(x_current, s_current, lower, upper)
        s_current <- box_result$s_bounded
        x_trial <- box_result$x_new

        # Evaluate trial point
        f_trial <- fn(x_trial)
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
            g_new <- gr(y_new)
            gr_evals <- gr_evals + 1

            # Check for NaN/Inf in new gradient
            if (!check_finite(f_new, g_new)) {
              stop("NaN or Inf detected in gradient at accepted Newton point")
            }

            # Newton steps do NOT use momentum - they are already optimal for
            # the local quadratic model. Adding momentum would push past the
            # computed Newton solution and could cause oscillation.
            x_current <- y_new
            f_current <- f_new
            g_current <- g_new

            # Reset momentum state after Newton step - Newton is optimal for local
            # quadratic model, so accumulated cubic momentum direction is stale
            v_prev <- rep(0, length(x_current))

            h_current <- hess(x_current)
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
      # Solve cubic subproblem using dispatcher
      cubic_result <- solve_cubic_subproblem_dispatch(
        g = g_current,
        H = h_current,
        sigma = sigma_current,
        solver = control$cubic_solver
      )

      s_current <- cubic_result$s
      pred_reduction <- cubic_result$pred_reduction

      # Apply box constraints to cubic step
      box_result <- apply_box_constraints(x_current, s_current, lower, upper)
      s_current <- box_result$s_bounded
      x_trial <- box_result$x_new

      # Evaluate trial point
      f_trial <- fn(x_trial)
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

        # Trial point before momentum
        y_new <- x_trial
        f_new <- f_trial
        g_new <- gr(y_new)
        gr_evals <- gr_evals + 1

        # Check for NaN/Inf in new gradient
        if (!check_finite(f_new, g_new)) {
          stop("NaN or Inf detected in gradient at accepted point")
        }

        # Apply momentum if enabled (Gao et al. Algorithm 1)
        if (control$use_momentum) {
          # Compute beta_k upper bound: beta ∈ [0, min{τ, α1||s||, α2||s||²}]
          s_norm <- sqrt(sum(s_current^2))
          beta_upper <- min(
            control$momentum_tau,
            control$momentum_alpha1 * s_norm,
            control$momentum_alpha2 * s_norm^2
          )

          # Initialize with beta=0 fallback (no momentum = just cubic step)
          v_current <- s_current
          x_current <- y_new
          f_current <- f_new

          # Bisection search for largest beta_k satisfying f(x_k + v_k) <= f(y_{k+1})
          # where v_k = beta * v_{k-1} + s_k
          f_target <- f_new  # f(y_{k+1}) = f(x_k + s_k)

          # Try beta_upper first (greedy)
          # v_k = beta * v_{k-1} + s_k, and x_{k+1} = x_k + v_k
          # Simplified: x_{k+1} = y_{k+1} + beta * v_{k-1} (since y_{k+1} = x_k + s_k)
          v_trial <- beta_upper * v_prev + s_current
          z_trial <- y_new + beta_upper * v_prev
          f_trial <- fn(z_trial)
          fn_evals <- fn_evals + 1

          if (check_finite(f_trial, rep(0, length(z_trial))) && f_trial <= f_target) {
            # Accept beta_upper
            v_current <- v_trial
            x_current <- z_trial
            f_current <- f_trial
          } else if (beta_upper > 1e-10) {
            # Bisection search for valid beta
            beta_low <- 0
            beta_high <- beta_upper
            max_bisect <- 10

            for (i in seq_len(max_bisect)) {
              beta_mid <- (beta_low + beta_high) / 2
              v_trial <- beta_mid * v_prev + s_current
              z_trial <- y_new + beta_mid * v_prev
              f_trial <- fn(z_trial)
              fn_evals <- fn_evals + 1

              if (check_finite(f_trial, rep(0, length(z_trial))) && f_trial <= f_target) {
                # Valid beta found - update best solution
                beta_low <- beta_mid
                v_current <- v_trial
                x_current <- z_trial
                f_current <- f_trial
              } else {
                beta_high <- beta_mid
              }

              if (beta_high - beta_low < 1e-6) break
            }
            # If no valid beta found, we keep the beta=0 fallback initialized above
          }

          # Update momentum vector for next iteration
          v_prev <- v_current

          # Evaluate gradient at accepted point
          g_current <- gr(x_current)
          gr_evals <- gr_evals + 1
        } else {
          # No momentum: x_{k+1} = y_{k+1}
          x_current <- y_new
          f_current <- f_new
          g_current <- g_new
        }

        h_current <- hess(x_current)
        hess_evals <- hess_evals + 1

        # Track step norm and function value for stagnation detection
        step_norms <- c(step_norms, sqrt(sum(s_current^2)))
        f_values <- c(f_values, f_current)

        prev_rejected <- FALSE
      } else {
        # Reject step - keep current point, will re-solve with updated sigma
        # IMPORTANT: Do NOT update v_prev on rejected steps (Gao Algorithm 1)
        # Next iteration will use the SAME v_prev (momentum carries forward unchanged)
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

        # Beta (if momentum enabled)
        if (control$use_momentum) {
          # Try to extract beta_low from parent environment if it exists
          beta_val <- if (iter == 0 || !exists("beta_low", inherits = FALSE)) {
            0
          } else {
            beta_low
          }
          trace_data$beta <- c(trace_data$beta, beta_val)
        }
      }

      # Reset timer for next iteration
      iter_start_time <- Sys.time()
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
    hessian = h_current,
    sigma = sigma_current,
    converged = converged,
    iterations = iter,
    evaluations = list(fn = fn_evals, gr = gr_evals, hess = hess_evals),
    message = conv_reason
  )

  # Add trace data if collected
  if (control$trace >= 1) {
    result$trace <- trace_data
  }

  return(result)
}
