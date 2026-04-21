# Tests for Quasi-Newton ARC Optimizer (Algorithm 4)
# ==================================================

# =============================================================================
# Test Functions
# =============================================================================

# Sphere function (convex quadratic)
sphere_fn <- function(x) sum(x^2)
sphere_gr <- function(x) 2 * x
sphere_hess <- function(x) 2 * diag(length(x))

# Rosenbrock function (n-dimensional)
rosenbrock_fn <- function(x) {
  n <- length(x)
  sum(100 * (x[2:n] - x[1:(n - 1)]^2)^2 + (1 - x[1:(n - 1)])^2)
}

rosenbrock_gr <- function(x) {
  n <- length(x)
  g <- numeric(n)
  for (i in 1:(n - 1)) {
    g[i] <- g[i] - 400 * x[i] * (x[i + 1] - x[i]^2) - 2 * (1 - x[i])
    g[i + 1] <- g[i + 1] + 200 * (x[i + 1] - x[i]^2)
  }
  g
}

rosenbrock_hess <- function(x) {
  n <- length(x)
  h <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    h[i, i] <- h[i, i] + 1200 * x[i]^2 - 400 * x[i + 1] + 2
    h[i, i + 1] <- h[i, i + 1] - 400 * x[i]
    h[i + 1, i] <- h[i + 1, i] - 400 * x[i]
    h[i + 1, i + 1] <- h[i + 1, i + 1] + 200
  }
  h
}

# Simple quadratic with known minimum
quadratic_fn <- function(x) {
  0.5 * (x[1]^2 + 10 * x[2]^2)  # Ill-conditioned
}
quadratic_gr <- function(x) c(x[1], 10 * x[2])
quadratic_hess <- function(x) diag(c(1, 10))


# =============================================================================
# Basic Functionality Tests
# =============================================================================

test_that("arcopt_qn forwards ... to fn and gr", {
  # Scaled sphere: f(x; s) = s * sum(x^2)
  scaled_fn <- function(x, scale) scale * sum(x^2)
  scaled_gr <- function(x, scale) scale * 2 * x

  result <- arcopt_qn(
    x0 = c(5, 5),
    fn = scaled_fn,
    gr = scaled_gr,
    scale = 2,
    control = list(qn_method = "sr1", trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0, 0), tolerance = 1e-4)
})

test_that("arcopt_qn forwards ... to hess in hybrid mode", {
  # Scaled sphere with Hessian
  scaled_fn <- function(x, scale) scale * sum(x^2)
  scaled_gr <- function(x, scale) scale * 2 * x
  scaled_hess <- function(x, scale) scale * 2 * diag(length(x))

  result <- arcopt_qn(
    x0 = c(5, 5),
    fn = scaled_fn,
    gr = scaled_gr,
    hess = scaled_hess,
    scale = 2,
    control = list(qn_method = "sr1", trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0, 0), tolerance = 1e-4)
  # Hessian should have been evaluated once for initialization
  expect_equal(result$evaluations$hess, 1)
})

test_that("arcopt_qn validates inputs correctly", {
  # Invalid x0

expect_error(
    arcopt_qn(c(1, NA), sphere_fn, sphere_gr),
    "x0 must be a finite numeric vector"
  )

  # Missing gradient
  expect_error(
    arcopt_qn(c(1, 1), sphere_fn, gr = NULL),
    "gr"
  )

  # Invalid QN method
  expect_error(
    arcopt_qn(c(1, 1), sphere_fn, sphere_gr,
              control = list(qn_method = "invalid")),
    "qn_method must be one of"
  )
})

test_that("arcopt_qn returns correct output structure", {
  result <- arcopt_qn(
    c(2, 2), sphere_fn, sphere_gr,
    control = list(maxit = 10, trace = 0)
  )

  expect_true(is.list(result))
  expect_true("par" %in% names(result))
  expect_true("value" %in% names(result))
  expect_true("gradient" %in% names(result))
  expect_true("converged" %in% names(result))
  expect_true("iterations" %in% names(result))
  expect_true("evaluations" %in% names(result))
  expect_true("qn_updates" %in% names(result))
  expect_true("qn_skips" %in% names(result))
  expect_true("qn_restarts" %in% names(result))

  expect_equal(length(result$par), 2)
  expect_true(is.numeric(result$value))
  expect_equal(length(result$gradient), 2)
})


# =============================================================================
# Pure QN Mode Tests (no Hessian provided)
# =============================================================================

test_that("arcopt_qn converges on sphere function with SR1", {
  x0 <- c(5, -3, 2)
  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "sr1", trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, rep(0, 3), tolerance = 1e-4)
  expect_equal(result$value, 0, tolerance = 1e-8)
})

test_that("arcopt_qn converges on sphere function with BFGS", {
  x0 <- c(5, -3, 2)
  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "bfgs", trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, rep(0, 3), tolerance = 1e-4)
})

test_that("arcopt_qn converges on sphere function with L-BFGS", {
  x0 <- c(5, -3, 2)
  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "lbfgs", qn_memory = 5, trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, rep(0, 3), tolerance = 1e-4)
})

test_that("arcopt_qn converges on sphere function with L-SR1", {
  x0 <- c(5, -3, 2)
  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "lsr1", qn_memory = 5, trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, rep(0, 3), tolerance = 1e-4)
})

test_that("arcopt_qn converges on ill-conditioned quadratic", {
  x0 <- c(10, 10)
  result <- arcopt_qn(
    x0, quadratic_fn, quadratic_gr,
    control = list(qn_method = "sr1", trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0, 0), tolerance = 1e-4)
})


# =============================================================================
# Hybrid Mode Tests (with Hessian)
# =============================================================================

test_that("arcopt_qn uses provided Hessian for initialization", {
  x0 <- c(2, 2)

  # With Hessian - should use exact initialization
  result_hybrid <- arcopt_qn(
    x0, sphere_fn, sphere_gr, sphere_hess,
    control = list(qn_method = "sr1", trace = 0)
  )

  # Without Hessian - uses scaled identity
  result_pure <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "sr1", trace = 0)
  )

  # Both should converge
  expect_true(result_hybrid$converged)
  expect_true(result_pure$converged)

  # Hybrid should use Hessian evaluation
  expect_equal(result_hybrid$evaluations$hess, 1)
  expect_equal(result_pure$evaluations$hess, 0)
})

test_that("arcopt_qn hybrid mode converges on Rosenbrock 2D", {
  x0 <- c(-1.2, 1)
  result <- arcopt_qn(
    x0, rosenbrock_fn, rosenbrock_gr, rosenbrock_hess,
    control = list(qn_method = "sr1", maxit = 500, trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, c(1, 1), tolerance = 1e-3)
})


# =============================================================================
# QN Update Tracking Tests
# =============================================================================

test_that("arcopt_qn tracks QN updates correctly", {
  x0 <- c(5, 5)
  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "sr1", trace = 0)
  )

  # Should have some updates
  expect_true(result$qn_updates >= 0)
  expect_true(result$qn_skips >= 0)
  expect_true(result$qn_restarts >= 0)

  # Total should be related to iterations
  # (each successful step attempts an update)
  expect_true(result$qn_updates + result$qn_skips <= result$iterations)
})

test_that("BFGS counts updates and skips", {
  # BFGS may skip when curvature condition violated
  x0 <- c(5, 5)
  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "bfgs", trace = 0)
  )

  expect_true(result$qn_updates >= 0)
  expect_true(result$qn_skips >= 0)
})


# =============================================================================
# Box Constraint Tests
# =============================================================================

test_that("arcopt_qn respects box constraints", {
  x0 <- c(0.5, 0.5)
  lower <- c(0.1, 0.1)
  upper <- c(0.9, 0.9)

  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    lower = lower, upper = upper,
    control = list(qn_method = "sr1", maxit = 500, trace = 0)
  )

  # Solution should respect bounds
  expect_true(all(result$par >= lower - 1e-8))
  expect_true(all(result$par <= upper + 1e-8))

  # Solution should be at lower bounds (minimum of sphere in [0.1, 0.9]^2)
  # Allow some tolerance since QN may not converge exactly at boundary
  expect_equal(result$par, lower, tolerance = 1e-3)
})

test_that("arcopt_qn validates bound dimensions", {
  expect_error(
    arcopt_qn(c(1, 1), sphere_fn, sphere_gr,
              lower = c(0), upper = c(2, 2)),
    "lower and upper must have same length"
  )

  expect_error(
    arcopt_qn(c(1, 1), sphere_fn, sphere_gr,
              lower = c(0, 0, 0), upper = c(2, 2, 2)),
    "lower and upper must have same length"
  )
})

test_that("arcopt_qn validates bound ordering", {
  expect_error(
    arcopt_qn(c(1, 1), sphere_fn, sphere_gr,
              lower = c(2, 0), upper = c(1, 2)),
    "lower must be strictly less than upper"
  )
})


# =============================================================================
# Rosenbrock Convergence Tests
# =============================================================================

test_that("arcopt_qn converges on Rosenbrock 2D (full-matrix methods)", {
  x0 <- c(-1.2, 1)

  # Only test full-matrix methods on Rosenbrock
  # Limited-memory methods (lbfgs, lsr1) struggle with ill-conditioning
  methods <- c("sr1", "bfgs")

  for (method in methods) {
    result <- arcopt_qn(
      x0, rosenbrock_fn, rosenbrock_gr,
      control = list(qn_method = method, maxit = 1000, trace = 0)
    )

    expect_true(
      result$converged || sqrt(sum(result$gradient^2)) < 1e-3,
      info = paste("Method:", method)
    )

    # Should be close to (1, 1)
    expect_equal(
      result$par, c(1, 1), tolerance = 0.1,
      info = paste("Method:", method)
    )
  }
})


# =============================================================================
# Limited Memory Tests
# =============================================================================

test_that("L-BFGS respects memory limit", {
  x0 <- c(5, -3, 2, 1, -1)

  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "lbfgs", qn_memory = 3, trace = 0)
  )

  expect_true(result$converged)
  # Hessian should be NULL for limited-memory
  expect_null(result$hessian)
})

test_that("L-SR1 respects memory limit", {
  x0 <- c(5, -3, 2, 1, -1)

  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "lsr1", qn_memory = 3, trace = 0)
  )

  expect_true(result$converged)
  expect_null(result$hessian)
})


# =============================================================================
# Comparison with Standard ARC
# =============================================================================

test_that("arcopt routes correctly with use_qn", {
  x0 <- c(2, 2)

  # Standard ARC (requires Hessian)
  result_std <- arcopt(
    x0, sphere_fn, sphere_gr, sphere_hess,
    control = list(trace = 0)
  )

  # QN-ARC via router
  result_qn <- arcopt(
    x0, sphere_fn, sphere_gr, sphere_hess,
    control = list(use_qn = TRUE, qn_method = "sr1", trace = 0)
  )

  # Both should converge
  expect_true(result_std$converged)
  expect_true(result_qn$converged)

  # Both should find the same minimum
  expect_equal(result_std$par, result_qn$par, tolerance = 1e-4)
})


# =============================================================================
# Edge Cases
# =============================================================================

test_that("arcopt_qn handles starting at optimum", {
  x0 <- c(0, 0, 0)

  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "sr1", trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0, 0, 0), tolerance = 1e-10)
  expect_equal(result$iterations, 0)  # Should converge immediately
})

test_that("arcopt_qn handles single dimension", {
  fn <- function(x) x^2
  gr <- function(x) 2 * x

  result <- arcopt_qn(
    5, fn, gr,
    control = list(qn_method = "sr1", trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, 0, tolerance = 1e-6)
})

test_that("arcopt_qn stops at maxit", {
  # Use very restrictive tolerance
  x0 <- c(100, 100)

  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "sr1", maxit = 5, gtol_abs = 1e-50, trace = 0)
  )

  expect_false(result$converged)
  expect_equal(result$iterations, 5)
  expect_equal(result$message, "maxit")
})


# =============================================================================
# Numerical Robustness
# =============================================================================

test_that("arcopt_qn handles NaN in objective", {
  bad_fn <- function(x) {
    if (x[1] < 0.1) return(NaN)
    sum(x^2)
  }
  bad_gr <- function(x) 2 * x

  x0 <- c(1, 1)

  # Should handle NaN by increasing sigma and rejecting step
  result <- arcopt_qn(
    x0, bad_fn, bad_gr,
    control = list(qn_method = "sr1", maxit = 100, trace = 0)
  )

  # May or may not converge, but shouldn't crash
  expect_true(is.list(result))
  expect_true(all(is.finite(result$par)))
})

test_that("arcopt_qn handles very small sigma_min", {
  x0 <- c(2, 2)

  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(
      qn_method = "sr1",
      sigma_min = 1e-12,
      trace = 0
    )
  )

  expect_true(result$converged)
})

test_that("arcopt_qn handles large sigma_max", {
  x0 <- c(2, 2)

  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(
      qn_method = "sr1",
      sigma_max = 1e20,
      trace = 0
    )
  )

  expect_true(result$converged)
})


# =============================================================================
# Higher Dimensional Tests
# =============================================================================

test_that("arcopt_qn handles moderate dimensions with L-BFGS", {
  n <- 20
  x0 <- rep(5, n)

  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "lbfgs", qn_memory = 10, trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, rep(0, n), tolerance = 1e-4)
})

test_that("arcopt_qn handles moderate dimensions with L-SR1", {
  n <- 20
  x0 <- rep(5, n)

  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "lsr1", qn_memory = 10, trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, rep(0, n), tolerance = 1e-4)
})


# =============================================================================
# Nesterov Acceleration Tests (Algorithm 4b)
# =============================================================================

test_that("arcopt_qn with acceleration finds near-optimal on sphere", {
  x0 <- c(5, -3, 2)
  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "bfgs", use_accel_qn = TRUE, maxit = 2000, trace = 0)
  )

  # Check that solution is near optimal (may not fully converge due to
  # acceleration extrapolation, but should be close)
  expect_equal(result$par, rep(0, 3), tolerance = 1e-3)
  expect_true(sqrt(sum(result$gradient^2)) < 0.01)
})

test_that("arcopt_qn with acceleration finds near-optimal on quadratic", {
  x0 <- c(10, 10)
  result <- arcopt_qn(
    x0, quadratic_fn, quadratic_gr,
    control = list(qn_method = "bfgs", use_accel_qn = TRUE, maxit = 2000, trace = 0)
  )

  # Check that solution is near optimal
  expect_equal(result$par, c(0, 0), tolerance = 1e-3)
  expect_true(sqrt(sum(result$gradient^2)) < 0.01)
})

test_that("arcopt_qn acceleration reaches similar solution as non-accelerated", {
  x0 <- c(10, 10)

  # Non-accelerated
  result_normal <- arcopt_qn(
    x0, quadratic_fn, quadratic_gr,
    control = list(qn_method = "bfgs", use_accel_qn = FALSE, trace = 0)
  )

  # Accelerated
  result_accel <- arcopt_qn(
    x0, quadratic_fn, quadratic_gr,
    control = list(qn_method = "bfgs", use_accel_qn = TRUE, maxit = 2000, trace = 0)
  )

  # Normal should converge
  expect_true(result_normal$converged)

  # Both should reach similar solution
  expect_equal(result_normal$par, result_accel$par, tolerance = 1e-2)
})

test_that("arcopt_qn acceleration works with SR1", {
  x0 <- c(5, 5, 5)
  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "sr1", use_accel_qn = TRUE, maxit = 2000, trace = 0)
  )

  expect_equal(result$par, rep(0, 3), tolerance = 1e-3)
})

test_that("arcopt_qn acceleration respects box constraints", {
  x0 <- c(0.5, 0.5)
  lower <- c(0.1, 0.1)
  upper <- c(0.9, 0.9)

  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    lower = lower, upper = upper,
    control = list(qn_method = "bfgs", use_accel_qn = TRUE, maxit = 1000, trace = 0)
  )

  # Solution should respect bounds
  expect_true(all(result$par >= lower - 1e-8))
  expect_true(all(result$par <= upper + 1e-8))
})

test_that("arcopt_qn acceleration handles single dimension", {
  fn <- function(x) x^2
  gr <- function(x) 2 * x

  result <- arcopt_qn(
    5, fn, gr,
    control = list(qn_method = "bfgs", use_accel_qn = TRUE, maxit = 2000, trace = 0)
  )

  expect_equal(result$par, 0, tolerance = 1e-3)
})


# =============================================================================
# Hybrid Method Tests
# =============================================================================

test_that("hybrid is the default qn_method", {
  # Verify default by checking result without specifying method
  x0 <- c(2, 2)
  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0, 0), tolerance = 1e-4)
})

test_that("arcopt_qn converges on sphere with hybrid method", {
  x0 <- c(5, -3, 2)
  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "hybrid", trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, rep(0, 3), tolerance = 1e-4)
})

test_that("arcopt_qn hybrid converges on Rosenbrock", {
  x0 <- c(-1.2, 1)
  result <- arcopt_qn(
    x0, rosenbrock_fn, rosenbrock_gr,
    control = list(qn_method = "hybrid", maxit = 1000, trace = 0)
  )

  # Rosenbrock minimum is at (1, 1)
  expect_true(result$converged || sqrt(sum(result$gradient^2)) < 1e-3)
  expect_equal(result$par, c(1, 1), tolerance = 0.05)
})

test_that("arcopt_qn hybrid handles ill-conditioned quadratic", {
  x0 <- c(10, 10)
  result <- arcopt_qn(
    x0, quadratic_fn, quadratic_gr,
    control = list(qn_method = "hybrid", trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0, 0), tolerance = 1e-4)
})

test_that("arcopt_qn hybrid respects bfgs_tol parameter", {
  x0 <- c(5, 5)

  # With very low bfgs_tol, BFGS should be used more often
  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "hybrid", bfgs_tol = 1e-15, trace = 0)
  )

  expect_true(result$converged)
})

test_that("arcopt_qn hybrid respects box constraints", {
  x0 <- c(0.5, 0.5)
  lower <- c(0.1, 0.1)
  upper <- c(0.9, 0.9)

  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    lower = lower, upper = upper,
    control = list(qn_method = "hybrid", trace = 0)
  )

  # Solution should respect bounds
  expect_true(all(result$par >= lower - 1e-8))
  expect_true(all(result$par <= upper + 1e-8))
  # Minimum is at (0,0) but constrained, so should be at (0.1, 0.1)
  expect_equal(result$par, c(0.1, 0.1), tolerance = 1e-4)
})

test_that("arcopt_qn hybrid tracks QN updates correctly", {
  x0 <- c(5, 5)
  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "hybrid", trace = 0)
  )

  # Should have some updates
  expect_true(result$qn_updates > 0)
  # Hybrid method should rarely skip since Powell is a fallback
  expect_true(result$qn_skips <= result$qn_updates)
})


# =============================================================================
# L-Hybrid Integration Tests
# =============================================================================

test_that("arcopt_qn converges with lhybrid on sphere", {
  x0 <- c(5, 5)
  result <- arcopt_qn(
    x0, sphere_fn, sphere_gr,
    control = list(qn_method = "lhybrid", trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0, 0), tolerance = 1e-4)
  expect_true(result$qn_updates > 0)
})

test_that("arcopt_qn converges with lhybrid on Rosenbrock", {
  x0 <- c(-1.2, 1)
  result <- arcopt_qn(
    x0, rosenbrock_fn, rosenbrock_gr,
    control = list(qn_method = "lhybrid", trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, c(1, 1), tolerance = 1e-4)
})

test_that("arcopt_qn converges with lhybrid on quadratic", {
  x0 <- c(10, 10)
  result <- arcopt_qn(
    x0, quadratic_fn, quadratic_gr,
    control = list(qn_method = "lhybrid", trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, c(0, 0), tolerance = 1e-4)
})

test_that("arcopt_qn lhybrid respects memory limit", {
  x0 <- rep(1, 20)  # Higher dimension
  result <- arcopt_qn(
    x0,
    fn = function(x) sum(x^2),
    gr = function(x) 2 * x,
    control = list(qn_method = "lhybrid", qn_memory = 5, trace = 0)
  )

  expect_true(result$converged)
  expect_equal(result$par, rep(0, 20), tolerance = 1e-3)
})

test_that("arcopt_qn lhybrid handles non-convex functions", {
  # Beale function (non-convex)
  beale_fn <- function(x) {
    (1.5 - x[1] * (1 - x[2]))^2 +
      (2.25 - x[1] * (1 - x[2]^2))^2 +
      (2.625 - x[1] * (1 - x[2]^3))^2
  }
  beale_gr <- function(x) {
    g <- numeric(2)
    y1 <- 1 - x[2]
    y2 <- 1 - x[2]^2
    y3 <- 1 - x[2]^3
    r1 <- 1.5 - x[1] * y1
    r2 <- 2.25 - x[1] * y2
    r3 <- 2.625 - x[1] * y3
    g[1] <- -2 * (r1 * y1 + r2 * y2 + r3 * y3)
    g[2] <- 2 * x[1] * (r1 + 2 * r2 * x[2] + 3 * r3 * x[2]^2)
    g
  }

  x0 <- c(0, 0)
  result <- arcopt_qn(
    x0, beale_fn, beale_gr,
    control = list(qn_method = "lhybrid", trace = 0, maxit = 200)
  )

  # Optimal at (3, 0.5)
  expect_true(result$converged)
  expect_equal(result$par, c(3, 0.5), tolerance = 1e-3)
})
