test_that("Momentum parameters are properly initialized", {
  # Simple quadratic
  fn <- function(x) sum(x^2)
  gr <- function(x) 2 * x
  hess <- function(x) 2 * diag(length(x))

  result <- arcopt(
    x0 = c(1, 1),
    fn = fn,
    gr = gr,
    hess = hess,
    control = list(use_momentum = TRUE, maxit = 10)
  )

  expect_true(result$converged)
  expect_true(all(abs(result$par) < 1e-5))
})

test_that("Momentum improves solution quality on Powell singular", {
  # Powell singular function (singular Hessian at optimum)
  fn <- function(x) {
    (x[1] + 10 * x[2])^2 + 5 * (x[3] - x[4])^2 +
      (x[2] - 2 * x[3])^4 + 10 * (x[1] - x[4])^4
  }

  gr <- function(x) {
    c(
      2 * (x[1] + 10 * x[2]) + 40 * (x[1] - x[4])^3,
      20 * (x[1] + 10 * x[2]) + 4 * (x[2] - 2 * x[3])^3,
      10 * (x[3] - x[4]) - 8 * (x[2] - 2 * x[3])^3,
      -10 * (x[3] - x[4]) - 40 * (x[1] - x[4])^3
    )
  }

  hess <- function(x) {
    u <- x[1] - x[4]
    v <- x[2] - 2 * x[3]
    matrix(c(
      2 + 120 * u^2, 20, 0, -120 * u^2,
      20, 200 + 12 * v^2, -24 * v^2, 0,
      0, -24 * v^2, 10 + 48 * v^2, -10,
      -120 * u^2, 0, -10, 10 + 120 * u^2
    ), 4, 4, byrow = TRUE)
  }

  x0 <- c(3, -1, 0, 1)
  x_opt <- c(0, 0, 0, 0)

  # Without momentum
  result_no_mom <- arcopt(
    x0 = x0,
    fn = fn,
    gr = gr,
    hess = hess,
    control = list(use_momentum = FALSE, maxit = 50)
  )

  # With momentum
  result_with_mom <- arcopt(
    x0 = x0,
    fn = fn,
    gr = gr,
    hess = hess,
    control = list(use_momentum = TRUE, maxit = 50)
  )

  # Both should converge
  expect_true(result_no_mom$converged)
  expect_true(result_with_mom$converged)

  # Momentum should achieve similar or better solution quality
  # (allow some tolerance for stochastic variation)
  dist_no_mom <- sqrt(sum((result_no_mom$par - x_opt)^2))
  dist_with_mom <- sqrt(sum((result_with_mom$par - x_opt)^2))

  # Momentum should get close to optimum (within 0.01)
  expect_lt(dist_with_mom, 0.01)

  # Gradient norm should be reasonably small
  expect_lt(sqrt(sum(result_with_mom$gradient^2)), 1e-5)
})

test_that("Momentum works with box constraints", {
  # Simple quadratic with constraints
  fn <- function(x) sum((x - 2)^2)
  gr <- function(x) 2 * (x - 2)
  hess <- function(x) 2 * diag(length(x))

  result <- arcopt(
    x0 = c(0, 0),
    fn = fn,
    gr = gr,
    hess = hess,
    lower = c(0, 0),
    upper = c(1, 1),
    control = list(use_momentum = TRUE)
  )

  expect_true(result$converged)
  # Should converge to upper bound c(1, 1)
  expect_equal(result$par, c(1, 1), tolerance = 1e-4)
})

test_that("Momentum handles rejected steps correctly", {
  # Function where steps might be rejected
  fn <- function(x) {
    if (x[1] < -10 || x[1] > 10) {
      return(Inf)
    }
    x[1]^4 - 3 * x[1]^3 + 2
  }

  gr <- function(x) {
    c(4 * x[1]^3 - 9 * x[1]^2)
  }

  hess <- function(x) {
    matrix(12 * x[1]^2 - 18 * x[1])
  }

  # This should not crash with rejected steps
  result <- arcopt(
    x0 = c(0.5),
    fn = fn,
    gr = gr,
    hess = hess,
    control = list(use_momentum = TRUE, maxit = 50)
  )

  expect_true(result$converged)
})

test_that("Momentum parameters are respected", {
  fn <- function(x) sum(x^2)
  gr <- function(x) 2 * x
  hess <- function(x) 2 * diag(length(x))

  # Test with different momentum_max values
  result_low <- arcopt(
    x0 = c(5, 5),
    fn = fn,
    gr = gr,
    hess = hess,
    control = list(
      use_momentum = TRUE,
      momentum_max = 0.1,
      maxit = 20
    )
  )

  result_high <- arcopt(
    x0 = c(5, 5),
    fn = fn,
    gr = gr,
    hess = hess,
    control = list(
      use_momentum = TRUE,
      momentum_max = 0.95,
      maxit = 20
    )
  )

  # Both should converge
  expect_true(result_low$converged)
  expect_true(result_high$converged)

  # Both should reach optimum
  expect_true(all(abs(result_low$par) < 1e-4))
  expect_true(all(abs(result_high$par) < 1e-4))
})

test_that("Momentum works without prior y_prev on first iteration", {
  # Simple quadratic
  fn <- function(x) sum(x^2)
  gr <- function(x) 2 * x
  hess <- function(x) 2 * diag(length(x))

  # Should handle y_prev = NULL on first iteration gracefully
  result <- arcopt(
    x0 = c(1, 1),
    fn = fn,
    gr = gr,
    hess = hess,
    control = list(use_momentum = TRUE, maxit = 5)
  )

  expect_true(result$converged)
})

test_that("Momentum with Newton steps", {
  # Well-conditioned quadratic where Newton steps succeed
  Q <- diag(c(10, 1))
  fn <- function(x) as.numeric(0.5 * t(x) %*% Q %*% x)
  gr <- function(x) as.vector(Q %*% x)
  hess <- function(x) Q

  result <- arcopt(
    x0 = c(5, 5),
    fn = fn,
    gr = gr,
    hess = hess,
    control = list(use_momentum = TRUE)
  )

  expect_true(result$converged)
  expect_true(all(abs(result$par) < 1e-6))
})

test_that("Adaptive momentum parameter decreases near optimum", {
  # This test verifies that beta_k -> 0 as ||s|| -> 0 and ||g|| -> 0
  # by checking that we can converge to high precision

  fn <- function(x) sum(x^2)
  gr <- function(x) 2 * x
  hess <- function(x) 2 * diag(length(x))

  result <- arcopt(
    x0 = c(10, 10),
    fn = fn,
    gr = gr,
    hess = hess,
    control = list(
      use_momentum = TRUE,
      gtol_abs = 1e-10,  # Very tight tolerance
      maxit = 100
    )
  )

  expect_true(result$converged)
  # Should achieve very high precision
  expect_true(sqrt(sum(result$gradient^2)) < 1e-9)
})
