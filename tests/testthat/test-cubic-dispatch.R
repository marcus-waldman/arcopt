test_that("Dispatcher auto-selects CG when hess_vec is provided", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  hess_vec <- hess_vec_fd(H)
  result <- solve_cubic_subproblem_dispatch(g, H = H, hess_vec = hess_vec,
                                            sigma = sigma, solver = "auto")

  expect_equal(result$solver_used, "cg")
  expect_equal(length(result$s), 2)
  expect_true(is.finite(result$pred_reduction))
})

test_that("Dispatcher auto-selects eigen for small problems", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  result <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                            solver = "auto",
                                            solver_threshold = 10)

  expect_equal(result$solver_used, "eigen")
  expect_equal(length(result$s), 2)
  expect_true(is.finite(result$pred_reduction))
})

test_that("Dispatcher auto-selects CG for large problems", {
  n <- 100
  H <- diag(runif(n, 0.5, 2))
  g <- rnorm(n)
  sigma <- 1.0

  result <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                            solver = "auto",
                                            solver_threshold = 50)

  expect_equal(result$solver_used, "cg")
  expect_equal(length(result$s), n)
  expect_true(is.finite(result$pred_reduction))
})

test_that("Dispatcher allows manual solver selection: ldl", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  result <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                            solver = "ldl")

  expect_equal(result$solver_used, "ldl")
  expect_equal(length(result$s), 2)
})

test_that("Dispatcher allows manual solver selection: eigen", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  result <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                            solver = "eigen")

  expect_equal(result$solver_used, "eigen")
  expect_equal(length(result$s), 2)
})

test_that("Dispatcher allows manual solver selection: cg", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  result <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                            solver = "cg")

  expect_equal(result$solver_used, "cg")
  expect_equal(length(result$s), 2)
})

test_that("Dispatcher validates inputs: requires H or hess_vec", {
  g <- c(1, 2)
  sigma <- 1.0

  expect_error(
    solve_cubic_subproblem_dispatch(g, sigma = sigma),
    "At least one of 'H' or 'hess_vec' must be provided"
  )
})

test_that("Dispatcher validates solver parameter", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  expect_error(
    solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                   solver = "invalid"),
    "solver must be one of"
  )
})

test_that("Dispatcher returns solver_used field", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  result <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                            solver = "eigen")

  expect_true("solver_used" %in% names(result))
  expect_true(is.character(result$solver_used))
  expect_equal(result$solver_used, "eigen")
})

test_that("Dispatcher passes additional arguments to solvers", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  # Pass max_lambda_iter to eigen solver
  result <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                            solver = "eigen",
                                            max_lambda_iter = 10)

  expect_equal(result$solver_used, "eigen")
  expect_true(is.finite(result$pred_reduction))
})

test_that("Dispatcher works with hess_vec only (no H)", {
  # Define Hessian implicitly via hess_vec
  hess_vec <- function(v) {
    c(2 * v[1], 3 * v[2])
  }

  g <- c(1, 2)
  sigma <- 1.0

  result <- solve_cubic_subproblem_dispatch(g, hess_vec = hess_vec,
                                            sigma = sigma, solver = "auto")

  expect_equal(result$solver_used, "cg")
  expect_equal(length(result$s), 2)
  expect_true(is.finite(result$pred_reduction))
})

test_that("Dispatcher constructs hess_vec from H for CG when needed", {
  n <- 100
  H <- diag(runif(n, 0.5, 2))
  g <- rnorm(n)
  sigma <- 1.0

  # Force CG solver without providing hess_vec
  result <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                            solver = "cg")

  expect_equal(result$solver_used, "cg")
  expect_equal(length(result$s), n)
})

test_that("Dispatcher produces consistent results across solvers", {
  H <- matrix(c(2, 0.5, 0.5, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  result_ldl <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                                solver = "ldl")
  result_eigen <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                                  solver = "eigen")
  result_cg <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                              solver = "cg")

  # All should produce similar solutions (within reasonable tolerance)
  expect_equal(result_eigen$s, result_ldl$s, tolerance = 1e-3)
  expect_equal(result_cg$s, result_eigen$s, tolerance = 0.1)
})

test_that("Dispatcher default threshold is 500", {
  # Create a problem with n = 400 (< 500)
  n <- 400
  H <- diag(runif(n, 0.5, 2))
  g <- rnorm(n)
  sigma <- 1.0

  result <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                            solver = "auto")

  # Should use eigen for n < 500
  expect_equal(result$solver_used, "eigen")
})

test_that("Dispatcher handles indefinite Hessians", {
  # Indefinite Hessian
  H <- matrix(c(-1, 0, 0, 2), 2, 2)
  g <- c(1, 1)
  sigma <- 1.0

  result_eigen <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                                  solver = "eigen")
  result_cg <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                              solver = "cg")

  # Both should handle indefiniteness
  expect_true(result_eigen$lambda >= 1)
  expect_true(result_cg$lambda >= 0)
  expect_true(is.finite(result_eigen$pred_reduction))
  expect_true(is.finite(result_cg$pred_reduction))
})
