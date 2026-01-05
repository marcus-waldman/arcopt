# Tests for Woodbury-Based Cubic Solver (Algorithm 5b)
# =====================================================

# =============================================================================
# Basic Functionality Tests
# =============================================================================

test_that("solve_cubic_simple handles zero gradient", {
  g <- c(0, 0, 0)
  result <- solve_cubic_simple(g, gamma = 1.0, sigma = 1.0)

  expect_equal(result$s, c(0, 0, 0))
  expect_equal(result$lambda, 0)
  expect_true(result$converged)
})

test_that("solve_cubic_simple returns correct solution for identity Hessian", {
  g <- c(1, 2, 3)
  gamma <- 1.0
  sigma <- 1.0

  result <- solve_cubic_simple(g, gamma, sigma)

  # Check secular equation: ||s|| = lambda / sigma
  norm_s <- sqrt(sum(result$s^2))
  expect_equal(norm_s, result$lambda / sigma, tolerance = 1e-8)

  # Check step direction: s = -g / (gamma + lambda)
  expected_s <- -g / (gamma + result$lambda)
  expect_equal(result$s, expected_s, tolerance = 1e-10)
})

test_that("solve_cubic_woodbury handles NULL w matrix", {
  g <- c(1, 2, 3)
  gamma <- 1.0
  sigma <- 1.0

  result <- solve_cubic_woodbury(g, gamma, w = NULL, c_mat = NULL, sigma)

  # Should fall back to simple solver
  expect_true(result$converged)
  norm_s <- sqrt(sum(result$s^2))
  expect_equal(norm_s, result$lambda / sigma, tolerance = 1e-8)
})

test_that("solve_cubic_woodbury handles empty w matrix", {
  g <- c(1, 2, 3)
  gamma <- 1.0
  sigma <- 1.0
  w <- matrix(0, nrow = 0, ncol = 3)
  c_mat <- matrix(0, nrow = 0, ncol = 0)

  result <- solve_cubic_woodbury(g, gamma, w, c_mat, sigma)

  expect_true(result$converged)
})


# =============================================================================
# Accuracy Tests: Compare with Eigendecomposition
# =============================================================================

test_that("solve_cubic_woodbury matches eigen solver on small problem", {
  set.seed(42)
  n <- 10
  m <- 3
  gamma <- 1.5
  sigma <- 2.0

  # Generate random w and c_mat
  w <- matrix(rnorm(m * n), m, n)
  c_mat <- diag(runif(m, 0.5, 2))  # Positive definite C

  # Build explicit B = gamma * I + w' c_mat w
  b_explicit <- gamma * diag(n) + t(w) %*% c_mat %*% w

  # Random gradient
  g <- rnorm(n)

  # Solve with Woodbury
  result_wood <- solve_cubic_woodbury(g, gamma, w, c_mat, sigma)

  # Solve with eigendecomposition (for comparison)
  result_eigen <- solve_cubic_subproblem_dispatch(
    g, b_explicit, sigma,
    solver = "eigen"
  )

  # Compare solutions
  expect_equal(result_wood$s, result_eigen$s, tolerance = 1e-6)
  expect_equal(result_wood$lambda, result_eigen$lambda, tolerance = 1e-6)
})

test_that("solve_cubic_woodbury satisfies secular equation", {
  set.seed(123)
  n <- 20
  m <- 5
  gamma <- 2.0
  sigma <- 1.5

  w <- matrix(rnorm(m * n), m, n)
  c_mat <- diag(runif(m, 0.1, 3))
  g <- rnorm(n)

  result <- solve_cubic_woodbury(g, gamma, w, c_mat, sigma)

  # Check secular equation: ||s|| = lambda / sigma
  norm_s <- sqrt(sum(result$s^2))
  expect_equal(norm_s, result$lambda / sigma, tolerance = 1e-8)
})


# =============================================================================
# Indefinite C Matrix Tests (SR1)
# =============================================================================

test_that("solve_cubic_woodbury handles indefinite c_mat from SR1", {
  set.seed(456)
  n <- 10
  m <- 3
  gamma <- 1.0
  sigma <- 1.0

  # Create indefinite C (typical for SR1)
  w <- matrix(rnorm(m * n), m, n)
  c_mat <- diag(c(-1, 2, -0.5))  # Indefinite

  g <- rnorm(n)

  # Should not crash and should satisfy secular equation
  result <- solve_cubic_woodbury(g, gamma, w, c_mat, sigma)

  expect_true(result$converged || result$iterations <= 50)
  norm_s <- sqrt(sum(result$s^2))

  # For indefinite case, secular equation should still approximately hold
  expect_equal(norm_s, result$lambda / sigma, tolerance = 1e-4)
})


# =============================================================================
# Compact Form Construction Tests
# =============================================================================

test_that("lbfgs_compact_form returns NULL for empty history", {
  history <- list(s = list(), y = list(), gamma = 1.0)
  result <- lbfgs_compact_form(history)
  expect_null(result)
})

test_that("lbfgs_compact_form builds correct matrices", {
  # Create simple L-BFGS history
  history <- list(
    s = list(c(1, 0, 0), c(0, 1, 0)),
    y = list(c(2, 0.1, 0), c(0.1, 3, 0.2)),
    gamma = 1.5
  )

  result <- lbfgs_compact_form(history)

  expect_false(is.null(result))
  expect_true(is.matrix(result$w))
  expect_true(is.matrix(result$c_mat))

  # w should be m x n
  expect_equal(ncol(result$w), 3)
})

test_that("lsr1_compact_form returns NULL for empty history", {
  history <- list(s = list(), y = list(), gamma = 1.0)
  result <- lsr1_compact_form(history)
  expect_null(result)
})

test_that("lsr1_compact_form builds valid matrices", {
  # Create simple L-SR1 history
  history <- list(
    s = list(c(1, 0, 0), c(0, 1, 0)),
    y = list(c(2.5, 0, 0), c(0, 2.5, 0)),
    gamma = 1.0
  )

  result <- lsr1_compact_form(history)

  expect_false(is.null(result))
  expect_true(is.matrix(result$w))
  expect_true(is.matrix(result$c_mat))
})


# =============================================================================
# Integration: L-SR1 History to Woodbury Solver
# =============================================================================

test_that("Woodbury solver works with L-SR1 compact form", {
  set.seed(789)
  n <- 5

  # Build L-SR1 history from simulated optimization
  history <- NULL
  x <- rnorm(n)
  g <- rnorm(n)

  for (i in 1:3) {
    s <- rnorm(n) * 0.1  # Small steps
    y <- rnorm(n) + 1.5 * s  # Ensure non-trivial update

    b_times_s <- if (is.null(history)) s else lsr1_multiply(history, s)
    history <- update_lsr1(history, s, y, m = 10, b_times_s = b_times_s)
  }

  # Get compact form
  compact <- lsr1_compact_form(history)

  if (!is.null(compact)) {
    # Solve cubic subproblem
    g <- rnorm(n)
    sigma <- 1.0

    result <- solve_cubic_woodbury(
      g, history$gamma, compact$w, compact$c_mat, sigma
    )

    expect_true(result$converged)

    # Verify secular equation
    norm_s <- sqrt(sum(result$s^2))
    expect_equal(norm_s, result$lambda / sigma, tolerance = 1e-6)
  }
})


# =============================================================================
# Edge Cases
# =============================================================================

test_that("solve_cubic_woodbury handles very small gradient", {
  g <- c(1e-12, 1e-12, 1e-12)
  gamma <- 1.0
  sigma <- 1.0
  w <- matrix(c(1, 0, 0), 1, 3)
  c_mat <- matrix(1, 1, 1)

  result <- solve_cubic_woodbury(g, gamma, w, c_mat, sigma)

  expect_true(result$converged)
  expect_true(sqrt(sum(result$s^2)) < 1e-6)
})

test_that("solve_cubic_woodbury validates gamma", {
  g <- c(1, 2, 3)
  expect_error(
    solve_cubic_woodbury(g, gamma = 0, w = NULL, c_mat = NULL, sigma = 1),
    "gamma must be positive"
  )
  expect_error(
    solve_cubic_woodbury(g, gamma = -1, w = NULL, c_mat = NULL, sigma = 1),
    "gamma must be positive"
  )
})

test_that("solve_cubic_woodbury validates sigma", {
  g <- c(1, 2, 3)
  expect_error(
    solve_cubic_woodbury(g, gamma = 1, w = NULL, c_mat = NULL, sigma = 0),
    "sigma must be positive"
  )
})

test_that("solve_cubic_simple converges with various sigma values", {
  g <- c(1, 2, 3)
  gamma <- 1.0

  for (sigma in c(0.01, 0.1, 1, 10, 100)) {
    result <- solve_cubic_simple(g, gamma, sigma)
    expect_true(result$converged)

    norm_s <- sqrt(sum(result$s^2))
    expect_equal(norm_s, result$lambda / sigma, tolerance = 1e-8)
  }
})


# =============================================================================
# Convergence Tests
# =============================================================================

test_that("solve_cubic_woodbury converges within max_iter", {
  set.seed(101)
  n <- 15
  m <- 4

  for (trial in 1:5) {
    w <- matrix(rnorm(m * n), m, n)
    c_mat <- diag(runif(m, 0.5, 3))
    g <- rnorm(n)
    gamma <- runif(1, 0.5, 3)
    sigma <- runif(1, 0.5, 5)

    result <- solve_cubic_woodbury(g, gamma, w, c_mat, sigma)

    # Either converged or used exactly max_iter
    expect_true(result$converged || result$iterations == 50)

    # Solution should satisfy secular equation reasonably well
    if (result$converged) {
      norm_s <- sqrt(sum(result$s^2))
      expect_equal(norm_s, result$lambda / sigma, tolerance = 1e-6)
    }
  }
})
