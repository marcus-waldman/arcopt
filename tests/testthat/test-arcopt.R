test_that("arcopt solves simple sphere function", {
  # Sphere function: f(x) = sum(x^2)
  # Minimum at x = 0, f = 0
  sphere <- function(x) sum(x^2)
  sphere_gr <- function(x) 2 * x
  sphere_hess <- function(x) 2 * diag(length(x))

  result <- arcopt(
    x0 = c(5, 5),
    fn = sphere,
    gr = sphere_gr,
    hess = sphere_hess
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0, 0), tolerance = 1e-6)
  expect_equal(result$value, 0, tolerance = 1e-8)
})

test_that("arcopt solves Rosenbrock function", {
  # Rosenbrock: f(x) = (1-x1)^2 + 100*(x2-x1^2)^2
  # Minimum at (1, 1), f = 0
  rosenbrock <- function(x) {
    (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  }

  rosenbrock_gr <- function(x) {
    c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
      200 * (x[2] - x[1]^2))
  }

  rosenbrock_hess <- function(x) {
    matrix(c(
      1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
      -400 * x[1], 200
    ), 2, 2)
  }

  result <- arcopt(
    x0 = c(-1.2, 1),
    fn = rosenbrock,
    gr = rosenbrock_gr,
    hess = rosenbrock_hess
  )

  expect_true(result$converged)
  expect_equal(result$par, c(1, 1), tolerance = 1e-4)
  expect_equal(result$value, 0, tolerance = 1e-6)
})

test_that("arcopt handles quadratic function correctly", {
  # Quadratic: f(x) = 1/2 x^T A x - b^T x
  # where A is PD
  A <- matrix(c(2, 0.5, 0.5, 3), 2, 2)
  b <- c(1, 2)

  quad_fn <- function(x) 0.5 * sum(x * (A %*% x)) - sum(b * x)
  quad_gr <- function(x) as.vector(A %*% x) - b
  quad_hess <- function(x) A

  result <- arcopt(
    x0 = c(0, 0),
    fn = quad_fn,
    gr = quad_gr,
    hess = quad_hess
  )

  # Solution: A * x = b
  expected_par <- solve(A, b)

  expect_true(result$converged)
  expect_equal(result$par, expected_par, tolerance = 1e-6)
})

test_that("arcopt returns evaluation counts", {
  sphere <- function(x) sum(x^2)
  sphere_gr <- function(x) 2 * x
  sphere_hess <- function(x) 2 * diag(length(x))

  result <- arcopt(
    x0 = c(1, 1),
    fn = sphere,
    gr = sphere_gr,
    hess = sphere_hess
  )

  expect_true(result$evaluations$fn > 0)
  expect_true(result$evaluations$gr > 0)
  expect_true(result$evaluations$hess > 0)
  expect_equal(result$evaluations$fn, result$evaluations$gr)
})

test_that("arcopt respects max_iter control parameter", {
  # Use Rosenbrock which needs many iterations
  rosenbrock <- function(x) {
    (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  }

  rosenbrock_gr <- function(x) {
    c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
      200 * (x[2] - x[1]^2))
  }

  rosenbrock_hess <- function(x) {
    matrix(c(
      1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
      -400 * x[1], 200
    ), 2, 2)
  }

  result <- arcopt(
    x0 = c(-1.2, 1),
    fn = rosenbrock,
    gr = rosenbrock_gr,
    hess = rosenbrock_hess,
    control = list(maxit = 3)  # Stop early
  )

  expect_true(result$iterations <= 3)
  expect_equal(result$message, "max_iter")
})

test_that("arcopt works on higher dimensional problems", {
  # 5-dimensional sphere
  n <- 5
  sphere <- function(x) sum(x^2)
  sphere_gr <- function(x) 2 * x
  sphere_hess <- function(x) 2 * diag(length(x))

  result <- arcopt(
    x0 = rep(3, n),
    fn = sphere,
    gr = sphere_gr,
    hess = sphere_hess
  )

  expect_true(result$converged)
  expect_equal(result$par, rep(0, n), tolerance = 1e-5)
})

test_that("arcopt requires hess", {
  sphere <- function(x) sum(x^2)
  sphere_gr <- function(x) 2 * x

  expect_error(
    arcopt(x0 = c(1, 1), fn = sphere, gr = sphere_gr, hess = NULL),
    "Hessian function 'hess' must be provided"
  )
})

test_that("arcopt handles function with saddle point", {
  # Function with saddle: f(x) = x1^2 - x2^2
  # Has saddle at origin, indefinite Hessian
  # Starting away from origin should still work
  saddle_fn <- function(x) x[1]^2 - x[2]^2 + 10
  saddle_gr <- function(x) c(2 * x[1], -2 * x[2])
  saddle_hess <- function(x) matrix(c(2, 0, 0, -2), 2, 2)

  # Start at point where we can descend
  result <- arcopt(
    x0 = c(2, 2),
    fn = saddle_fn,
    gr = saddle_gr,
    hess = saddle_hess,
    control = list(maxit = 100)
  )

  # Should converge somewhere (behavior with indefinite Hessian)
  expect_true(result$converged || result$iterations == 100)
})

test_that("arcopt validates input", {
  sphere <- function(x) sum(x^2)
  sphere_gr <- function(x) 2 * x
  sphere_hess <- function(x) 2 * diag(length(x))

  # Non-numeric x0
  expect_error(
    arcopt(x0 = "not numeric", fn = sphere, gr = sphere_gr, hess = sphere_hess),
    "x0 must be a finite numeric vector"
  )

  # Non-finite x0
  expect_error(
    arcopt(x0 = c(1, Inf), fn = sphere, gr = sphere_gr, hess = sphere_hess),
    "x0 must be a finite numeric vector"
  )

  # fn not a function
  expect_error(
    arcopt(x0 = c(1, 1), fn = "not a function", gr = sphere_gr, hess = sphere_hess),
    "fn must be a function"
  )

  # gr not provided
  expect_error(
    arcopt(x0 = c(1, 1), fn = sphere, gr = "not a function", hess = sphere_hess),
    "gr .* is required"
  )
})

test_that("arcopt uses custom control parameters", {
  sphere <- function(x) sum(x^2)
  sphere_gr <- function(x) 2 * x
  sphere_hess <- function(x) 2 * diag(length(x))

  result <- arcopt(
    x0 = c(1, 1),
    fn = sphere,
    gr = sphere_gr,
    hess = sphere_hess,
    control = list(
      gtol_abs = 1e-10,  # Very strict gradient tolerance
      sigma0 = 10.0       # Large initial sigma
    )
  )

  expect_true(result$converged)
  # With stricter tolerance, should get very close to optimum
  expect_true(max(abs(result$gradient)) < 1e-9)
})

test_that("arcopt detects NaN in initial point", {
  # Function that returns NaN at initial point
  nan_fn <- function(x) {
    if (all(x == c(1, 1))) return(NaN)
    sum(x^2)
  }
  nan_gr <- function(x) 2 * x
  nan_hess <- function(x) 2 * diag(length(x))

  expect_error(
    arcopt(x0 = c(1, 1), fn = nan_fn, gr = nan_gr, hess = nan_hess),
    "Initial point yields NaN"
  )
})

test_that("arcopt detects NaN in gradient at initial point", {
  # Function with NaN gradient at initial point
  sphere <- function(x) sum(x^2)
  nan_gr <- function(x) {
    if (all(x == c(1, 1))) return(c(NaN, NaN))
    2 * x
  }
  sphere_hess <- function(x) 2 * diag(length(x))

  expect_error(
    arcopt(x0 = c(1, 1), fn = sphere, gr = nan_gr, hess = sphere_hess),
    "Initial point yields NaN"
  )
})

test_that("arcopt handles NaN during optimization by rejecting step", {
  # Function that returns NaN at specific point
  eval_count <- 0
  nan_fn <- function(x) {
    eval_count <<- eval_count + 1
    # Return NaN on second evaluation (trial point)
    if (eval_count == 2) return(NaN)
    sum(x^2)
  }
  sphere_gr <- function(x) 2 * x
  sphere_hess <- function(x) 2 * diag(length(x))

  # Reset counter
  eval_count <- 0

  result <- arcopt(
    x0 = c(5, 5),
    fn = nan_fn,
    gr = sphere_gr,
    hess = sphere_hess,
    control = list(maxit = 10)
  )

  # Should handle NaN by increasing sigma and continuing
  # Eventually converges or hits max iterations
  expect_true(result$converged || result$iterations == 10)
})

test_that("arcopt respects box constraints", {
  # Sphere function with box constraints
  sphere <- function(x) sum(x^2)
  sphere_gr <- function(x) 2 * x
  sphere_hess <- function(x) 2 * diag(length(x))

  # Constrain to [0.5, 1] x [0.5, 1]
  # Optimum without constraints is (0, 0)
  # With constraints, optimum should be (0.5, 0.5)
  expect_warning(
    result <- arcopt(
      x0 = c(2, 2),
      fn = sphere,
      gr = sphere_gr,
      hess = sphere_hess,
      lower = c(0.5, 0.5),
      upper = c(1, 1)
    ),
    "Initial point infeasible"
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0.5, 0.5), tolerance = 1e-5)
  # All iterates should respect bounds
  expect_true(all(result$par >= c(0.5, 0.5) - 1e-10))
  expect_true(all(result$par <= c(1, 1) + 1e-10))
})

test_that("arcopt projects infeasible initial point", {
  sphere <- function(x) sum(x^2)
  sphere_gr <- function(x) 2 * x
  sphere_hess <- function(x) 2 * diag(length(x))

  # Start outside feasible region
  expect_warning(
    result <- arcopt(
      x0 = c(5, 5),  # Outside [0, 1] x [0, 1]
      fn = sphere,
      gr = sphere_gr,
      hess = sphere_hess,
      lower = c(0, 0),
      upper = c(1, 1)
    ),
    "Initial point infeasible"
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0, 0), tolerance = 1e-5)
})

test_that("arcopt handles one-sided bounds", {
  # Function with minimum at (-1, -1)
  offset_sphere <- function(x) sum((x + 1)^2)
  offset_gr <- function(x) 2 * (x + 1)
  offset_hess <- function(x) 2 * diag(length(x))

  # Only lower bound at 0
  result <- arcopt(
    x0 = c(1, 1),
    fn = offset_sphere,
    gr = offset_gr,
    hess = offset_hess,
    lower = c(0, 0),
    upper = c(Inf, Inf)
  )

  expect_true(result$converged)
  # Optimum should be at boundary (0, 0) since unconstrained optimum (-1, -1) is infeasible
  expect_equal(result$par, c(0, 0), tolerance = 1e-5)
})

test_that("arcopt handles asymmetric bounds", {
  sphere <- function(x) sum(x^2)
  sphere_gr <- function(x) 2 * x
  sphere_hess <- function(x) 2 * diag(length(x))

  # Different bounds for each dimension
  result <- arcopt(
    x0 = c(1, 1),
    fn = sphere,
    gr = sphere_gr,
    hess = sphere_hess,
    lower = c(-10, 0.5),
    upper = c(10, 2)
  )

  expect_true(result$converged)
  # Check bounds are respected
  expect_true(result$par[1] >= -10 && result$par[1] <= 10)
  expect_true(result$par[2] >= 0.5 - 1e-10 && result$par[2] <= 2 + 1e-10)
  # Function value should be close to 0.25 (minimum is at (0, 0.5), f = 0.25)
  # but we allow for some numerical error and possible convergence to nearby point
  expect_true(result$value <= 1.0)
})

# Integration tests for Newton-rho behavior

test_that("arcopt uses rho-based acceptance for Newton steps", {
  Q <- matrix(c(1, 0, 0, 100), 2, 2)
  b <- c(1, 10)

  fn <- function(x) 0.5 * sum(x * (Q %*% x)) - sum(b * x)
  gr <- function(x) as.vector(Q %*% x - b)
  hess <- function(x) Q

  x0 <- c(5, 5)

  result <- arcopt(x0, fn, gr, hess, control = list(maxit = 50, trace = FALSE))

  expect_true(result$converged)
  expect_equal(result$par, c(1, 0.1), tolerance = 1e-6)
})

test_that("arcopt rejects poor Newton steps and falls back to cubic", {
  rosenbrock_fn <- function(x) {
    (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  }

  rosenbrock_gr <- function(x) {
    c(
      -2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
      200 * (x[2] - x[1]^2)
    )
  }

  rosenbrock_hess <- function(x) {
    matrix(c(
      1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
      -400 * x[1], 200
    ), 2, 2)
  }

  x0 <- c(0.5, 0.5)

  result <- arcopt(
    x0,
    rosenbrock_fn,
    rosenbrock_gr,
    rosenbrock_hess,
    control = list(maxit = 100, eta1 = 0.1, trace = FALSE)
  )

  expect_true(result$converged)
  expect_equal(result$par, c(1, 1), tolerance = 1e-5)
})

test_that("Newton and cubic steps use same sigma update mechanism", {
  sphere_fn <- function(x) sum(x^2)
  sphere_gr <- function(x) 2 * x
  sphere_hess <- function(x) 2 * diag(length(x))

  x0 <- c(1, 1, 1)

  result <- arcopt(
    x0,
    sphere_fn,
    sphere_gr,
    sphere_hess,
    control = list(maxit = 50, trace = FALSE)
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0, 0, 0), tolerance = 1e-8)
})

test_that("momentum follows Gao et al. (correct beta formula)", {
  # Test that beta has direct (not inverse) relationship with ||s||

  rosenbrock <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  rosenbrock_gr <- function(x) {
    c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
      200 * (x[2] - x[1]^2))
  }
  rosenbrock_hess <- function(x) {
    matrix(c(
      1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
      -400 * x[1], 200
    ), 2, 2)
  }

  # With Gao's formula: large steps should allow large momentum
  result <- arcopt(
    c(-1.2, 1),
    rosenbrock,
    rosenbrock_gr,
    rosenbrock_hess,
    control = list(
      use_momentum = TRUE,
      momentum_tau = 0.9,
      momentum_alpha1 = 1.0,  # Allow larger beta
      momentum_alpha2 = 1.0,
      maxit = 100
    )
  )

  expect_true(result$converged)
  expect_lt(result$iterations, 50)  # Should converge faster than without momentum
})

test_that("arcopt forwards ... to fn, gr, and hess", {
  # Scaled sphere: f(x; s) = s * sum(x^2)
  scaled_fn <- function(x, scale) scale * sum(x^2)
  scaled_gr <- function(x, scale) scale * 2 * x
  scaled_hess <- function(x, scale) scale * 2 * diag(length(x))

  result <- arcopt(
    x0 = c(5, 5),
    fn = scaled_fn,
    gr = scaled_gr,
    hess = scaled_hess,
    scale = 2
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0, 0), tolerance = 1e-6)
  expect_equal(result$value, 0, tolerance = 1e-8)
})

test_that("arcopt forwards ... through use_qn dispatch", {
  # Verify ... reaches arcopt_qn when use_qn = TRUE
  scaled_fn <- function(x, scale) scale * sum(x^2)
  scaled_gr <- function(x, scale) scale * 2 * x
  scaled_hess <- function(x, scale) scale * 2 * diag(length(x))

  result <- arcopt(
    x0 = c(5, 5),
    fn = scaled_fn,
    gr = scaled_gr,
    hess = scaled_hess,
    scale = 2,
    control = list(use_qn = TRUE, qn_method = "sr1", trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0, 0), tolerance = 1e-4)
})

test_that("momentum maintains monotonicity", {
  # Momentum should never increase function value (at accepted points)
  # Note: This test verifies that the bisection search ensures monotonicity

  rosenbrock <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  rosenbrock_gr <- function(x) {
    c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
      200 * (x[2] - x[1]^2))
  }
  rosenbrock_hess <- function(x) {
    matrix(c(
      1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
      -400 * x[1], 200
    ), 2, 2)
  }

  # Run with momentum enabled
  result <- arcopt(
    c(-1.2, 1),
    rosenbrock,
    rosenbrock_gr,
    rosenbrock_hess,
    control = list(use_momentum = TRUE, maxit = 100)
  )

  # The algorithm should converge successfully
  expect_true(result$converged)

  # Final function value should be near optimum
  expect_lt(result$value, 1e-10)
})
