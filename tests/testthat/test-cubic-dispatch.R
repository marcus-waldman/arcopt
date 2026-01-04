test_that("Dispatcher auto-selects eigen by default", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  result <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                            solver = "auto")

  expect_equal(result$solver_used, "eigen")
  expect_equal(length(result$s), 2)
  expect_true(is.finite(result$pred_reduction))
})

test_that("Dispatcher errors on removed LDL solver", {
  H <- matrix(c(2, 0, 0, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  expect_error(
    solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma, solver = "ldl"),
    "LDL cubic solver has been removed"
  )
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

test_that("Dispatcher validates inputs: requires H", {
  g <- c(1, 2)
  sigma <- 1.0

  expect_error(
    solve_cubic_subproblem_dispatch(g, sigma = sigma),
    "Hessian matrix 'H' must be provided"
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

test_that("Auto solver selects eigendecomposition", {
  H <- matrix(c(2, 0.5, 0.5, 3), 2, 2)
  g <- c(1, 2)
  sigma <- 1.0

  result_auto <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                                 solver = "auto")
  result_eigen <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                                  solver = "eigen")

  # Auto should produce identical results to explicit eigen
  expect_equal(result_auto$s, result_eigen$s)
  expect_equal(result_auto$lambda, result_eigen$lambda)
  expect_equal(result_auto$solver_used, "eigen")
})

test_that("Dispatcher handles indefinite Hessians", {
  # Indefinite Hessian
  H <- matrix(c(-1, 0, 0, 2), 2, 2)
  g <- c(1, 1)
  sigma <- 1.0

  result_eigen <- solve_cubic_subproblem_dispatch(g, H = H, sigma = sigma,
                                                  solver = "eigen")

  # Should handle indefiniteness
  expect_true(result_eigen$lambda >= 1)
  expect_true(is.finite(result_eigen$pred_reduction))
})
