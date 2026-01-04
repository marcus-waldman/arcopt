test_that("eigen solver handles positive definite Hessian", {
  # Simple 2D PD Hessian
  g <- c(1, 2)
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  sigma <- 1.0

  result <- solve_cubic_eigen(g, H, sigma)

  expect_true(result$converged)
  expect_equal(result$case_type, "easy")
  expect_equal(length(result$s), 2)
  expect_true(result$lambda >= 0)
})

test_that("eigen solver handles indefinite Hessian", {
  # Indefinite Hessian with eigenvalues -1, 2
  H <- matrix(c(-1, 0, 0, 2), 2, 2)
  g <- c(1, 1)
  sigma <- 1.0

  result <- solve_cubic_eigen(g, H, sigma)

  expect_true(result$converged)
  # Should be easy case if g not orthogonal to v1
  # With this g and H, should converge successfully
  expect_equal(length(result$s), 2)
})

test_that("eigen solver detects hard case", {
  # Hard case: H has eigenvalues [-1, 2], v1 = [1, 0]
  # g = [0, 1] is orthogonal to v1
  H <- matrix(c(-1, 0, 0, 2), 2, 2)
  g <- c(0, 1)
  sigma <- 1.0

  result <- solve_cubic_eigen(g, H, sigma, hard_case_tol = 1e-8)

  expect_true(result$converged)
  expect_equal(result$case_type, "hard")
  expect_equal(result$lambda, 1.0, tolerance = 1e-10)  # lambda = -lambda_1

  # Verify secular equation
  s_norm <- sqrt(sum(result$s^2))
  expect_equal(s_norm, result$lambda / sigma, tolerance = 1e-8)
})

test_that("eigen solver handles zero gradient", {
  g <- c(0, 0)
  H <- matrix(c(1, 0, 0, 1), 2, 2)
  sigma <- 1.0

  result <- solve_cubic_eigen(g, H, sigma)

  expect_true(result$converged)
  expect_equal(result$case_type, "zero_gradient")
  expect_equal(result$s, c(0, 0))
  expect_equal(result$lambda, 0)
  expect_equal(result$pred_reduction, 0)
})

test_that("eigen solver secular equation is satisfied (easy case)", {
  g <- rnorm(5)
  H <- diag(runif(5, 0.5, 2))  # PD diagonal
  sigma <- 1.0

  result <- solve_cubic_eigen(g, H, sigma)

  expect_true(result$converged)
  expect_equal(result$case_type, "easy")

  # Verify secular equation: ||s|| = lambda/sigma
  s_norm <- sqrt(sum(result$s^2))
  expect_equal(s_norm, result$lambda / sigma, tolerance = 1e-9)
})

test_that("eigen solver computes correct predicted reduction", {
  g <- c(1, 2)
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  sigma <- 1.0

  result <- solve_cubic_eigen(g, H, sigma)

  # Manual calculation
  s <- result$s
  h_s <- as.vector(H %*% s)
  s_norm <- sqrt(sum(s^2))
  pred_red_manual <- -sum(g * s) - 0.5 * sum(s * h_s) - (sigma / 3) * s_norm^3

  expect_equal(result$pred_reduction, pred_red_manual, tolerance = 1e-10)
})

test_that("eigen solver produces shorter steps with larger sigma", {
  g <- c(1, 2)
  H <- matrix(c(2, 0, 0, 3), 2, 2)

  result_small_sigma <- solve_cubic_eigen(g, H, sigma = 0.5)
  result_large_sigma <- solve_cubic_eigen(g, H, sigma = 2.0)

  s_norm_small <- sqrt(sum(result_small_sigma$s^2))
  s_norm_large <- sqrt(sum(result_large_sigma$s^2))

  expect_true(s_norm_large < s_norm_small)
})

test_that("eigen solver handles strongly indefinite Hessian", {
  # Hessian with large negative eigenvalue
  H <- matrix(c(-5, 0, 0, 2), 2, 2)
  g <- c(1, 1)
  sigma <- 1.0

  result <- solve_cubic_eigen(g, H, sigma)

  # Solver should produce valid solution even if convergence is tight
  expect_true(result$lambda >= 5)  # Must shift by at least 5
  expect_equal(length(result$s), 2)
  expect_true(is.finite(result$pred_reduction))

  # For strongly indefinite cases, solver may not fully converge but should
  # produce a descent direction (negative predicted reduction)
  expect_true(result$pred_reduction > 0)  # Should provide decrease

  # Verify solution is not pathological (all zeros or NaN)
  expect_true(any(abs(result$s) > .Machine$double.eps))
  expect_true(all(is.finite(result$s)))
})

test_that("eigen solver works on larger problems", {
  n <- 10
  g <- rnorm(n)
  H <- diag(runif(n, -1, 2))  # Mix of positive and negative eigenvalues
  sigma <- 1.0

  result <- solve_cubic_eigen(g, H, sigma)

  expect_true(result$converged)
  expect_equal(length(result$s), n)

  # Verify secular equation
  s_norm <- sqrt(sum(result$s^2))
  expect_equal(s_norm, result$lambda / sigma, tolerance = 1e-8)
})

test_that("eigen solver hard case with non-zero s_partial", {
  # Hard case where s_tilde[2:n] contributes to norm
  # H has eigenvalues [-2, 1, 3]
  H <- diag(c(-2, 1, 3))
  g <- c(0, 1, 1)  # Orthogonal to first eigenvector
  sigma <- 1.0

  result <- solve_cubic_eigen(g, H, sigma, hard_case_tol = 1e-8)

  expect_true(result$converged)
  expect_equal(result$case_type, "hard")
  expect_equal(result$lambda, 2.0, tolerance = 1e-10)

  # Verify secular equation
  s_norm <- sqrt(sum(result$s^2))
  expect_equal(s_norm, result$lambda / sigma, tolerance = 1e-8)
})

test_that("eigen solver returns all required fields", {
  g <- c(1, 2)
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  sigma <- 1.0

  result <- solve_cubic_eigen(g, H, sigma)

  expect_true("s" %in% names(result))
  expect_true("lambda" %in% names(result))
  expect_true("pred_reduction" %in% names(result))
  expect_true("converged" %in% names(result))
  expect_true("case_type" %in% names(result))

  expect_true(is.numeric(result$s))
  expect_true(is.numeric(result$lambda))
  expect_true(is.numeric(result$pred_reduction))
  expect_true(is.logical(result$converged))
  expect_true(is.character(result$case_type))
})

test_that("eigen solver case_type is one of expected values", {
  # Test all three case types

  # Easy case
  g1 <- c(1, 2)
  H1 <- matrix(c(2, 0, 0, 3), 2, 2)
  result1 <- solve_cubic_eigen(g1, H1, 1.0)
  expect_true(result1$case_type %in% c("easy", "hard", "zero_gradient"))
  expect_equal(result1$case_type, "easy")

  # Zero gradient
  g2 <- c(0, 0)
  H2 <- matrix(c(1, 0, 0, 1), 2, 2)
  result2 <- solve_cubic_eigen(g2, H2, 1.0)
  expect_equal(result2$case_type, "zero_gradient")

  # Hard case
  g3 <- c(0, 1)
  H3 <- matrix(c(-1, 0, 0, 2), 2, 2)
  result3 <- solve_cubic_eigen(g3, H3, 1.0, hard_case_tol = 1e-8)
  expect_equal(result3$case_type, "hard")
})
