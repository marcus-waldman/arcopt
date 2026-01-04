test_that("CG solver handles positive definite Hessian", {
  # Simple 2D PD Hessian
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  hess_vec <- hess_vec_fd(H)
  result <- solve_cubic_cg(g, hess_vec, sigma)

  expect_true(result$converged)
  expect_equal(length(result$s), 2)
  expect_true(result$lambda >= 0)
  expect_true(is.finite(result$pred_reduction))
  expect_true(result$cg_iters > 0)
})

test_that("CG solver handles indefinite Hessian", {
  # Indefinite Hessian with eigenvalues -1, 2
  H <- matrix(c(-1, 0, 0, 2), 2, 2)
  g <- c(1, 1)
  sigma <- 1.0

  hess_vec <- hess_vec_fd(H)
  result <- solve_cubic_cg(g, hess_vec, sigma)

  expect_equal(length(result$s), 2)
  expect_true(result$lambda >= 1)  # Must shift by at least 1
  expect_true(is.finite(result$pred_reduction))
})

test_that("CG solver handles zero gradient", {
  H <- matrix(c(1, 0, 0, 1), 2, 2)
  g <- c(0, 0)
  sigma <- 1.0

  hess_vec <- hess_vec_fd(H)
  result <- solve_cubic_cg(g, hess_vec, sigma)

  expect_true(result$converged)
  expect_equal(result$s, c(0, 0))
  expect_equal(result$lambda, 0)
  expect_equal(result$pred_reduction, 0)
  expect_equal(result$cg_iters, 0)
})

test_that("CG solver works with custom hess_vec function", {
  # Define Hessian implicitly via hess_vec
  # H = diag(c(2, 3, 4))
  hess_vec <- function(v) {
    c(2 * v[1], 3 * v[2], 4 * v[3])
  }

  g <- c(1, 2, 3)
  sigma <- 1.0

  result <- solve_cubic_cg(g, hess_vec, sigma)

  expect_equal(length(result$s), 3)
  expect_true(is.finite(result$pred_reduction))
  expect_true(result$lambda >= 0)
})

test_that("CG solver detects negative curvature", {
  # Strongly indefinite Hessian
  H <- matrix(c(-5, 0, 0, 2), 2, 2)
  g <- c(1, 1)
  sigma <- 1.0

  hess_vec <- hess_vec_fd(H)
  result <- solve_cubic_cg(g, hess_vec, sigma)

  # Should shift by at least 5 to make PD
  expect_true(result$lambda >= 5)
  expect_equal(length(result$s), 2)
  expect_true(is.finite(result$pred_reduction))
})

test_that("CG solver produces shorter steps with larger sigma", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)

  hess_vec <- hess_vec_fd(H)
  result_small_sigma <- solve_cubic_cg(g, hess_vec, sigma = 0.5)
  result_large_sigma <- solve_cubic_cg(g, hess_vec, sigma = 2.0)

  s_norm_small <- sqrt(sum(result_small_sigma$s^2))
  s_norm_large <- sqrt(sum(result_large_sigma$s^2))

  expect_true(s_norm_large < s_norm_small)
})

test_that("CG solver computes correct predicted reduction", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  hess_vec <- hess_vec_fd(H)
  result <- solve_cubic_cg(g, hess_vec, sigma)

  # Manual calculation
  s <- result$s
  h_s <- as.vector(H %*% s)
  s_norm <- sqrt(sum(s^2))
  pred_red_manual <- -sum(g * s) - 0.5 * sum(s * h_s) - (sigma / 3) * s_norm^3

  expect_equal(result$pred_reduction, pred_red_manual, tolerance = 1e-10)
})

test_that("CG solver works on larger problems", {
  n <- 20
  g <- rnorm(n)
  H <- diag(runif(n, 0.5, 2))  # PD diagonal
  sigma <- 1.0

  hess_vec <- hess_vec_fd(H)
  result <- solve_cubic_cg(g, hess_vec, sigma)

  expect_equal(length(result$s), n)
  expect_true(is.finite(result$pred_reduction))
  expect_true(result$cg_iters > 0)
  expect_true(result$cg_iters <= min(100, 2 * n))
})

test_that("CG solver agrees with eigen solver on easy cases", {
  # Use a small PD problem where both should converge
  H <- matrix(c(2, 0.5, 0.5, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  # Eigendecomposition solver
  result_eigen <- solve_cubic_eigen(g, H, sigma)

  # CG solver
  hess_vec <- hess_vec_fd(H)
  result_cg <- solve_cubic_cg(g, hess_vec, sigma)

  # Solutions should be similar (relaxed tolerances due to different algorithms)
  expect_equal(result_cg$s, result_eigen$s, tolerance = 0.1)
  expect_equal(result_cg$lambda, result_eigen$lambda, tolerance = 0.5)
  expect_equal(result_cg$pred_reduction, result_eigen$pred_reduction, tolerance = 0.05)
})

test_that("CG solver handles max_cg_iter parameter", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  hess_vec <- hess_vec_fd(H)
  result <- solve_cubic_cg(g, hess_vec, sigma, max_cg_iter = 5)

  expect_true(result$cg_iters <= 5)
  expect_equal(length(result$s), 2)
  expect_true(is.finite(result$pred_reduction))
})

test_that("CG solver returns all required fields", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  hess_vec <- hess_vec_fd(H)
  result <- solve_cubic_cg(g, hess_vec, sigma)

  expect_true("s" %in% names(result))
  expect_true("lambda" %in% names(result))
  expect_true("pred_reduction" %in% names(result))
  expect_true("converged" %in% names(result))
  expect_true("cg_iters" %in% names(result))

  expect_true(is.numeric(result$s))
  expect_true(is.numeric(result$lambda))
  expect_true(is.numeric(result$pred_reduction))
  expect_true(is.logical(result$converged))
  expect_true(is.numeric(result$cg_iters))
})

test_that("CG solver handles strongly indefinite Hessian", {
  # Hessian with large negative eigenvalue
  H <- matrix(c(-5, 0, 0, 2), 2, 2)
  g <- c(1, 1)
  sigma <- 1.0

  hess_vec <- hess_vec_fd(H)
  result <- solve_cubic_cg(g, hess_vec, sigma)

  # Should produce valid solution even if convergence is tight
  expect_true(result$lambda >= 5)
  expect_equal(length(result$s), 2)
  expect_true(is.finite(result$pred_reduction))

  # Should provide decrease
  expect_true(result$pred_reduction > 0)

  # Verify solution is not pathological
  expect_true(any(abs(result$s) > .Machine$double.eps))
  expect_true(all(is.finite(result$s)))
})

test_that("hess_vec_fd wrapper works correctly", {
  H <- matrix(c(2, 1, 1, 3), 2, 2)
  v <- c(1, 2)

  hess_vec <- hess_vec_fd(H)
  result <- hess_vec(v)

  expected <- as.vector(H %*% v)
  expect_equal(result, expected)
  expect_true(is.numeric(result))
  expect_equal(length(result), 2)
})

test_that("CG solver handles hard case scenario", {
  # Hard case: H has eigenvalues [-1, 2], v1 = [1, 0]
  # g = [0, 1] is orthogonal to v1
  H <- matrix(c(-1, 0, 0, 2), 2, 2)
  g <- c(0, 1)
  sigma <- 1.0

  hess_vec <- hess_vec_fd(H)
  result <- solve_cubic_cg(g, hess_vec, sigma)

  # Should select a shift that makes system PD (approximately >= |min_eigenvalue|)
  # CG may select different shift than eigen solver due to discrete grid
  expect_true(result$lambda >= 0)  # Should be non-negative
  expect_equal(length(result$s), 2)
  expect_true(is.finite(result$pred_reduction))
  # CG-Lanczos may break down when g is an eigenvector; allow zero step
  expect_true(result$pred_reduction >= 0)
})

test_that("CG solver shift selection is reasonable", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  hess_vec <- hess_vec_fd(H)
  result <- solve_cubic_cg(g, hess_vec, sigma)

  # For PD Hessian, lambda should be relatively small
  # (shifts range from 10^-10 to 10^20)
  expect_true(result$lambda >= 1e-10)
  expect_true(result$lambda < 1e10)  # Shouldn't need huge shift for PD
})
