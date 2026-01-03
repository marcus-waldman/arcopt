test_that("cubic solver handles PD Hessian (sphere function)", {
  # Sphere function: H = I, g = random
  # Should approximately get Newton direction s = -g when sigma is small
  set.seed(123)
  n <- 5
  g <- rnorm(n)
  H <- diag(n)
  sigma <- 0.1

  result <- solve_cubic_subproblem(g, H, sigma)

  expect_true(result$converged)
  # For small sigma and PD H, should be in same direction as Newton direction
  # and have reasonable magnitude
  newton_step <- -g
  # Check direction is similar (dot product close to ||s|| * ||newton||)
  cosine_similarity <- sum(result$s * newton_step) / (sqrt(sum(result$s^2)) * sqrt(sum(newton_step^2)))
  expect_true(cosine_similarity > 0.9)
})

test_that("cubic solver handles indefinite Hessian", {
  # Indefinite Hessian: one negative eigenvalue
  g <- c(1, 1)
  H <- matrix(c(1, 0, 0, -1), 2, 2)  # eigenvalues: 1, -1
  sigma <- 1.0

  result <- solve_cubic_subproblem(g, H, sigma)

  expect_true(result$converged)
  expect_equal(length(result$s), 2)
  # Lambda should be positive to make H + lambda*I positive definite
  expect_true(result$lambda > 0)
  # Verify secular equation is satisfied: ||s|| ≈ lambda/sigma
  s_norm <- sqrt(sum(result$s^2))
  expect_equal(s_norm, result$lambda / sigma, tolerance = 1e-6)
})

test_that("cubic solver produces short steps with large sigma", {
  # Large sigma should produce small steps (high regularization)
  g <- c(1, 1)
  H <- diag(2)
  sigma_large <- 100

  result <- solve_cubic_subproblem(g, H, sigma_large)

  expect_true(result$converged)
  s_norm <- sqrt(sum(result$s^2))
  # With large sigma, step should be small
  expect_true(s_norm < 0.5)
})

test_that("cubic solver produces longer steps with small sigma", {
  # Small sigma should produce larger steps (less regularization)
  g <- c(1, 1)
  H <- diag(2)
  sigma_small <- 0.01

  result <- solve_cubic_subproblem(g, H, sigma_small)

  expect_true(result$converged)
  s_norm <- sqrt(sum(result$s^2))
  # With small sigma, step should be closer to Newton direction
  # Newton step would be s = -g = c(-1, -1), norm = sqrt(2) ≈ 1.414
  expect_true(s_norm > 1.0)
})

test_that("cubic solver handles zero gradient", {
  # Zero gradient should give zero step
  g <- c(0, 0, 0)
  H <- diag(3)
  sigma <- 1.0

  result <- solve_cubic_subproblem(g, H, sigma)

  expect_true(result$converged)
  expect_equal(result$s, c(0, 0, 0), tolerance = 1e-10)
  expect_equal(result$lambda, 0, tolerance = 1e-10)
})

test_that("cubic solver computes correct predicted reduction", {
  # Verify pred_red = -g^T s - 1/2 s^T H s - sigma/3 ||s||^3
  g <- c(1, 2)
  H <- matrix(c(2, 0.5, 0.5, 3), 2, 2)
  sigma <- 1.0

  result <- solve_cubic_subproblem(g, H, sigma)

  # Manually compute predicted reduction
  s <- result$s
  h_s <- as.vector(H %*% s)
  s_norm <- sqrt(sum(s^2))
  expected_pred_red <- -sum(g * s) - 0.5 * sum(s * h_s) - (sigma / 3) * s_norm^3

  expect_equal(result$pred_reduction, expected_pred_red, tolerance = 1e-10)
  # Predicted reduction should be positive (we're minimizing)
  expect_true(result$pred_reduction > 0)
})

test_that("cubic solver secular equation is satisfied", {
  # For any valid solution, ||s(lambda)|| = lambda/sigma must hold
  set.seed(456)
  n <- 4
  g <- rnorm(n)
  # Create random PD Hessian
  A <- matrix(rnorm(n * n), n, n)
  H <- t(A) %*% A + diag(n)  # Guaranteed PD
  sigma <- 0.5

  result <- solve_cubic_subproblem(g, H, sigma)

  expect_true(result$converged)
  s_norm <- sqrt(sum(result$s^2))
  # Secular equation: ||s|| = lambda/sigma
  expect_equal(s_norm, result$lambda / sigma, tolerance = 1e-8)
})

test_that("cubic solver handles nearly singular PD Hessian", {
  # Hessian with very small positive eigenvalue
  g <- c(1, 1)
  H <- matrix(c(1e-6, 0, 0, 1), 2, 2)
  sigma <- 1.0

  result <- solve_cubic_subproblem(g, H, sigma)

  expect_true(result$converged)
  expect_equal(length(result$s), 2)
  # Should still produce valid step
  expect_true(all(is.finite(result$s)))
})

test_that("cubic solver handles larger problem", {
  # Test on larger dimensional problem
  set.seed(789)

  n <- 10
  g <- rnorm(n)
  # Create PD Hessian
  A <- matrix(rnorm(n * n), n, n)
  H <- t(A) %*% A + diag(n)
  sigma <- 1.0

  result <- solve_cubic_subproblem(g, H, sigma)

  expect_true(result$converged)
  expect_equal(length(result$s), n)
  # Verify secular equation
  s_norm <- sqrt(sum(result$s^2))
  expect_equal(s_norm, result$lambda / sigma, tolerance = 1e-6)
})

test_that("cubic solver handles diagonal indefinite Hessian requiring exact shift", {

  # This is the edge case that previously caused "leading minor not positive" error
  # H = diag(1, -1) requires lambda = 1 + eps to make H + lambda*I strictly PD
  # (not just semi-definite)
  g <- c(1, 1)
  H <- matrix(c(1, 0, 0, -1), 2, 2)  # eigenvalues: 1, -1
  sigma <- 1.0

  result <- solve_cubic_subproblem(g, H, sigma)

  expect_true(result$converged)
  # Lambda must be > 1 (the absolute value of the most negative eigenvalue)
  # to make H + lambda*I strictly positive definite
  expect_true(result$lambda > 1)
  # Verify secular equation
  s_norm <- sqrt(sum(result$s^2))
  expect_equal(s_norm, result$lambda / sigma, tolerance = 1e-6)
})

test_that("cubic solver handles strongly indefinite Hessian", {
  # More challenging indefinite Hessian with off-diagonal elements
  # This tests the LDL-based shift estimation
  g <- c(1, 1)
  H <- matrix(c(2, 1, 1, -2), 2, 2)  # eigenvalues: ~2.24, ~-2.24
  sigma <- 1.0

  result <- solve_cubic_subproblem(g, H, sigma)

  expect_true(result$converged)
  # Lambda must be > 2.24 (abs of most negative eigenvalue)
  expect_true(result$lambda > 2.2)
  # Verify secular equation
  s_norm <- sqrt(sum(result$s^2))
  expect_equal(s_norm, result$lambda / sigma, tolerance = 1e-6)
})
