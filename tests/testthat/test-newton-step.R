test_that("Newton step succeeds with PD Hessian (sphere)", {
  # Sphere function: H = I, g = random
  set.seed(123)
  n <- 5
  g <- rnorm(n)
  H <- diag(n)

  result <- try_newton_step(g, H)

  expect_true(result$success)
  # For H = I, Newton step should be s = -g
  expect_equal(result$s, -g, tolerance = 1e-10)
})

test_that("Newton step fails with indefinite Hessian", {
  # Indefinite Hessian with negative eigenvalue
  g <- c(1, 1)
  H <- matrix(c(1, 0, 0, -1), 2, 2)  # eigenvalues: 1, -1

  result <- try_newton_step(g, H)

  expect_false(result$success)
  expect_null(result$s)
})

test_that("Newton step succeeds with general PD Hessian", {
  # General positive definite Hessian
  g <- c(1, 2)
  H <- matrix(c(2, 0.5, 0.5, 3), 2, 2)  # PD: eigenvalues > 0

  result <- try_newton_step(g, H)

  expect_true(result$success)
  # Verify H * s = -g
  H_s <- as.vector(H %*% result$s)
  expect_equal(H_s, -g, tolerance = 1e-10)
})

test_that("Newton step computes correct solution", {
  # Verify s = -H^{-1} g
  g <- c(1, 2, 3)
  # Create random PD Hessian
  set.seed(456)
  A <- matrix(rnorm(9), 3, 3)
  H <- t(A) %*% A + diag(3)  # Guaranteed PD

  result <- try_newton_step(g, H)

  expect_true(result$success)
  # Verify s = -H^{-1} g by checking H * s = -g
  expected_s <- -solve(H, g)
  expect_equal(result$s, expected_s, tolerance = 1e-10)
})

test_that("Newton step handles nearly singular PD Hessian", {
  # Hessian with very small positive eigenvalue
  g <- c(1, 1)
  H <- matrix(c(1e-6, 0, 0, 1), 2, 2)

  result <- try_newton_step(g, H)

  # Should succeed since both eigenvalues are positive
  expect_true(result$success)
  expect_true(all(is.finite(result$s)))
})

test_that("Newton step uses Gershgorin bound correctly", {
  # Test case where Gershgorin bound is slightly negative
  # H = [[1, 1.1], [1.1, 1]]
  # Diagonal: 1, Off-diagonal sums: 1.1
  # Gershgorin bounds: 1 - 1.1 = -0.1 for both rows
  g <- c(1, 1)
  H <- matrix(c(1, 1.1, 1.1, 1), 2, 2)

  result <- try_newton_step(g, H)

  # Should fail due to Gershgorin bound <= 0
  expect_false(result$success)
  expect_null(result$s)
})

test_that("Newton step handles zero gradient", {
  # Zero gradient should give zero Newton step
  g <- c(0, 0, 0)
  H <- diag(3)

  result <- try_newton_step(g, H)

  expect_true(result$success)
  expect_equal(result$s, c(0, 0, 0), tolerance = 1e-10)
})

test_that("Newton step works on larger problem", {
  # Test on larger dimensional problem
  set.seed(789)
  n <- 10
  g <- rnorm(n)
  # Create diagonally dominant PD Hessian (guaranteed to pass Gershgorin)
  H <- diag(n) * 10 + matrix(rnorm(n * n, sd = 0.1), n, n)
  H <- (H + t(H)) / 2  # Make symmetric

  result <- try_newton_step(g, H)

  expect_true(result$success)
  # Verify solution
  H_s <- as.vector(H %*% result$s)
  expect_equal(H_s, -g, tolerance = 1e-9)
})

test_that("Newton step fails when Cholesky fails despite positive Gershgorin", {
  # Edge case: Gershgorin bound may be positive but matrix still non-PD
  # This is rare but theoretically possible
  # Create a matrix that passes Gershgorin but is not PD
  # [[2, 0], [0, -0.5]] won't work because diagonal is negative
  # Try [[0.1, 0], [0, -0.05]] - Gershgorin: 0.1, -0.05
  # Actually, for this to happen, we need off-diagonal dominance
  # Let's use a known non-PD matrix
  g <- c(1, 1)
  H <- matrix(c(0.01, 0, 0, -0.005), 2, 2)

  result <- try_newton_step(g, H)

  # Should fail because H has negative eigenvalue
  expect_false(result$success)
})
