test_that("TR solver returns interior Newton step when it fits", {
  # PD Hessian, very large radius -> unconstrained Newton step is feasible
  g <- c(1, 2)
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  radius <- 100

  result <- solve_tr_eigen(g, H, radius)

  expect_true(result$converged)
  expect_equal(result$case_type, "interior")
  expect_false(result$on_boundary)
  expect_equal(result$lambda, 0)

  # Verify Newton step satisfies H s == -g
  expect_equal(as.vector(H %*% result$s), -g, tolerance = 1e-10)

  # And lies strictly inside the radius
  expect_true(sqrt(sum(result$s^2)) < radius)
})

test_that("TR solver returns boundary solution when radius is binding", {
  # PD Hessian with tight radius
  g <- c(1, 2)
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  radius <- 0.1

  result <- solve_tr_eigen(g, H, radius)

  expect_true(result$converged)
  expect_equal(result$case_type, "easy")
  expect_true(result$on_boundary)
  expect_true(result$lambda > 0)

  # Verify ||s|| == radius
  expect_equal(sqrt(sum(result$s^2)), radius, tolerance = 1e-8)

  # Verify (H + lambda I) s == -g
  residual <- (H + result$lambda * diag(2)) %*% result$s + g
  expect_equal(as.vector(residual), c(0, 0), tolerance = 1e-8)
})

test_that("TR solver handles indefinite Hessian on boundary", {
  # Indefinite Hessian: eigenvalues -1, 2
  H <- matrix(c(-1, 0, 0, 2), 2, 2)
  g <- c(1, 1)
  radius <- 1.0

  result <- solve_tr_eigen(g, H, radius)

  expect_true(result$converged)
  expect_true(result$on_boundary)
  expect_true(result$lambda > 1.0) # must shift past |lambda_min|

  # ||s|| == radius
  expect_equal(sqrt(sum(result$s^2)), radius, tolerance = 1e-8)
})

test_that("TR solver detects hard case", {
  # H = diag(-1, 2), v_1 = e_1. g = [0, 1] is orthogonal to v_1.
  # Hard-case lambda = 1. Step along v_1 should be added to reach boundary.
  H <- matrix(c(-1, 0, 0, 2), 2, 2)
  g <- c(0, 1)
  radius <- 1.0

  result <- solve_tr_eigen(g, H, radius, hard_case_tol = 1e-8)

  expect_true(result$converged)
  expect_equal(result$case_type, "hard")
  expect_true(result$on_boundary)
  expect_equal(result$lambda, 1.0, tolerance = 1e-10)

  # ||s|| == radius
  expect_equal(sqrt(sum(result$s^2)), radius, tolerance = 1e-8)

  # Step should include a non-trivial v_1 = e_1 component
  expect_true(abs(result$s[1]) > 0)
})

test_that("TR solver handles zero gradient", {
  g <- c(0, 0)
  H <- matrix(c(1, 0, 0, 1), 2, 2)
  radius <- 1.0

  result <- solve_tr_eigen(g, H, radius)

  expect_true(result$converged)
  expect_equal(result$case_type, "zero_gradient")
  expect_equal(result$s, c(0, 0))
  expect_equal(result$lambda, 0)
  expect_equal(result$pred_reduction, 0)
  expect_false(result$on_boundary)
})

test_that("TR solver secular equation is satisfied on boundary", {
  set.seed(42)
  g <- rnorm(5)
  H <- diag(runif(5, 0.5, 2)) # PD diagonal
  radius <- 0.5 # tight enough to bind

  result <- solve_tr_eigen(g, H, radius)

  expect_true(result$converged)
  expect_true(result$on_boundary)
  expect_equal(result$case_type, "easy")

  expect_equal(sqrt(sum(result$s^2)), radius, tolerance = 1e-9)
})

test_that("TR solver computes correct predicted reduction", {
  g <- c(1, 2)
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  radius <- 0.3

  result <- solve_tr_eigen(g, H, radius)

  # Manual quadratic-model decrease: pred = -g's - 0.5 s'Hs
  s <- result$s
  h_s <- as.vector(H %*% s)
  pred_red_manual <- -sum(g * s) - 0.5 * sum(s * h_s)

  expect_equal(result$pred_reduction, pred_red_manual, tolerance = 1e-10)

  # Must provide decrease
  expect_true(result$pred_reduction > 0)
})

test_that("TR solver produces longer steps with larger radius", {
  g <- c(1, 2)
  H <- matrix(c(2, 0, 0, 3), 2, 2)

  result_tight <- solve_tr_eigen(g, H, radius = 0.1)
  result_loose <- solve_tr_eigen(g, H, radius = 0.5)

  s_tight <- sqrt(sum(result_tight$s^2))
  s_loose <- sqrt(sum(result_loose$s^2))

  expect_true(s_loose > s_tight)
})

test_that("TR solver handles strongly indefinite Hessian", {
  H <- matrix(c(-5, 0, 0, 2), 2, 2)
  g <- c(1, 1)
  radius <- 1.0

  result <- solve_tr_eigen(g, H, radius)

  # Must shift by at least 5 to make (H + lambda I) PD
  expect_true(result$lambda >= 5 - 1e-6)
  expect_equal(length(result$s), 2)
  expect_true(is.finite(result$pred_reduction))
  expect_true(result$pred_reduction > 0) # provides decrease
  expect_true(all(is.finite(result$s)))
})

test_that("TR solver works on larger problems", {
  set.seed(99)
  n <- 10
  g <- rnorm(n)
  H <- diag(runif(n, -1, 2)) # mix of signs
  radius <- 0.8

  result <- solve_tr_eigen(g, H, radius)

  expect_true(result$converged)
  expect_equal(length(result$s), n)
  expect_true(is.finite(result$pred_reduction))
})

test_that("TR solver hard case with non-zero s_partial (3D)", {
  # H has eigenvalues [-2, 1, 3]; v_1 = e_1. g orthogonal to e_1.
  H <- diag(c(-2, 1, 3))
  g <- c(0, 1, 1)
  radius <- 1.5

  result <- solve_tr_eigen(g, H, radius, hard_case_tol = 1e-8)

  expect_true(result$converged)
  expect_equal(result$case_type, "hard")
  expect_equal(result$lambda, 2.0, tolerance = 1e-10)
  expect_equal(sqrt(sum(result$s^2)), radius, tolerance = 1e-8)
})

test_that("TR solver returns all required fields", {
  g <- c(1, 2)
  H <- matrix(c(2, 0, 0, 3), 2, 2)

  result <- solve_tr_eigen(g, H, radius = 0.5)

  for (field in c("s", "lambda", "pred_reduction", "converged",
                  "case_type", "on_boundary")) {
    expect_true(field %in% names(result), info = paste("missing:", field))
  }

  expect_true(is.numeric(result$s))
  expect_true(is.numeric(result$lambda))
  expect_true(is.numeric(result$pred_reduction))
  expect_true(is.logical(result$converged))
  expect_true(is.character(result$case_type))
  expect_true(is.logical(result$on_boundary))
})

test_that("TR solver case_type covers all expected values", {
  # interior
  g1 <- c(1, 2)
  H1 <- matrix(c(2, 0, 0, 3), 2, 2)
  res1 <- solve_tr_eigen(g1, H1, radius = 100)
  expect_equal(res1$case_type, "interior")

  # easy boundary
  res2 <- solve_tr_eigen(g1, H1, radius = 0.1)
  expect_equal(res2$case_type, "easy")

  # zero_gradient
  res3 <- solve_tr_eigen(c(0, 0), H1, radius = 1)
  expect_equal(res3$case_type, "zero_gradient")

  # hard
  H4 <- matrix(c(-1, 0, 0, 2), 2, 2)
  g4 <- c(0, 1)
  res4 <- solve_tr_eigen(g4, H4, radius = 1.0, hard_case_tol = 1e-8)
  expect_equal(res4$case_type, "hard")
})

test_that("TR solver agrees with trust::trust on easy PD problem", {
  skip_if_not_installed("trust")
  # Quadratic objective: f(x) = g'x + 0.5 x'Hx, starting from 0.
  # For a quadratic, trust::trust converges in one accepted step.
  g <- c(1.2, -0.8)
  H <- matrix(c(2.5, 0.4, 0.4, 3.1), 2, 2)
  radius <- 0.4

  tr_fn <- function(x) {
    list(
      value = sum(g * x) + 0.5 * sum(x * (H %*% x)),
      gradient = g + as.vector(H %*% x),
      hessian = H
    )
  }

  # Use small rmax so the first step stays at the requested radius
  t_result <- trust::trust(tr_fn, parinit = c(0, 0),
                           rinit = radius, rmax = radius,
                           iterlim = 1, fterm = 0, mterm = 0)
  step_trust <- t_result$argument

  ours <- solve_tr_eigen(g, H, radius)

  # Step norms should agree to a few digits (both should be on boundary)
  expect_equal(sqrt(sum(step_trust^2)), sqrt(sum(ours$s^2)),
               tolerance = 1e-4)

  # Step directions should align (cosine > 0.999)
  cos_sim <- sum(step_trust * ours$s) /
    (sqrt(sum(step_trust^2)) * sqrt(sum(ours$s^2)))
  expect_true(cos_sim > 0.999)
})
