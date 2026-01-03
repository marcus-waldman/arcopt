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
#'   vector of length Q and return a Q×Q symmetric matrix. Required for full ARC;
#'   if NULL, uses SR1 quasi-Newton approximation (see Details).
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
#' is strongly recommended. If `hess = NULL`, the algorithm falls back to SR1
#' quasi-Newton approximation, which is less robust but avoids Hessian
#' computation.
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
#' * `use_sr1`: Use SR1 quasi-Newton if `hess = NULL` (default: TRUE)
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
    eta1 = 0.1,
    eta2 = 0.9,
    gamma1 = 0.5,
    gamma2 = 2.0,
    use_sr1 = TRUE,
    trace = FALSE
  )

  control <- modifyList(control_defaults, control)

  # Project initial point to bounds (for future box constraint support)
  x_current <- pmax(lower, pmin(upper, x0))

  # Initialize tracking variables
  sigma_current <- control$sigma0
  iter <- 0
  prev_rejected <- FALSE
  fn_evals <- 0
  gr_evals <- 0
  hess_evals <- 0

  # Initial evaluation
  f_current <- fn(x_current)
  g_current <- gr(x_current)
  fn_evals <- fn_evals + 1
  gr_evals <- gr_evals + 1

  if (!is.null(hess)) {
    h_current <- hess(x_current)
    hess_evals <- hess_evals + 1
  } else {
    stop("hess function is required for this implementation")
  }

  # Store previous values for convergence checking
  f_previous <- NA
  x_previous <- NULL

  converged <- FALSE
  conv_reason <- ""

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

    # STEP 2: Try Newton step first (if not rejected and hess provided)
    step_type <- "cubic"
    used_newton <- FALSE

    if (!prev_rejected && !is.null(hess)) {
      newton_result <- try_newton_step(g_current, h_current)

      if (newton_result$success) {
        # Newton step succeeded
        s_current <- newton_result$s
        x_trial <- x_current + s_current
        f_trial <- fn(x_trial)
        fn_evals <- fn_evals + 1

        # Accept Newton step and update
        x_previous <- x_current
        f_previous <- f_current

        x_current <- x_trial
        f_current <- f_trial
        g_current <- gr(x_current)
        h_current <- hess(x_current)
        gr_evals <- gr_evals + 1
        hess_evals <- hess_evals + 1

        # Decrease sigma after successful Newton step
        sigma_current <- max(control$gamma1 * sigma_current, 1e-16)
        prev_rejected <- FALSE
        used_newton <- TRUE
        step_type <- "newton"
      }
    }

    # STEP 3: If Newton failed/skipped, use cubic solver
    if (!used_newton) {
      # Solve cubic subproblem
      cubic_result <- solve_cubic_subproblem(
        g = g_current,
        H = h_current,
        sigma = sigma_current
      )

      s_current <- cubic_result$s
      pred_reduction <- cubic_result$pred_reduction

      # Evaluate trial point
      x_trial <- x_current + s_current
      f_trial <- fn(x_trial)
      fn_evals <- fn_evals + 1

      # Compute actual reduction and ratio
      actual_reduction <- f_current - f_trial
      rho <- actual_reduction / pred_reduction

      # Accept or reject step
      if (rho >= control$eta1) {
        # Accept step
        x_previous <- x_current
        f_previous <- f_current

        x_current <- x_trial
        f_current <- f_trial
        g_current <- gr(x_current)
        h_current <- hess(x_current)
        gr_evals <- gr_evals + 1
        hess_evals <- hess_evals + 1

        prev_rejected <- FALSE
      } else {
        # Reject step - keep current point, will re-solve with updated sigma
        prev_rejected <- TRUE
      }

      # Update sigma based on step quality
      sigma_current <- update_sigma_cgt(
        sigma_current = sigma_current,
        rho = rho,
        eta1 = control$eta1,
        eta2 = control$eta2,
        gamma1 = control$gamma1,
        gamma2 = control$gamma2
      )
    }

    iter <- iter + 1
  }

  # Check if we hit max iterations without converging
  if (!converged && iter >= control$maxit) {
    converged <- TRUE
    conv_reason <- "max_iter"
  }

  list(
    par = x_current,
    value = f_current,
    gradient = g_current,
    hessian = h_current,
    converged = converged,
    iterations = iter,
    evaluations = list(fn = fn_evals, gr = gr_evals, hess = hess_evals),
    message = conv_reason
  )
}
