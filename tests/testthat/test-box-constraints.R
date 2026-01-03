test_that("project_to_box projects correctly", {
  # Test basic projection
  x <- c(-2, 0, 5, 10)
  lower <- c(-1, -1, -1, -1)
  upper <- c(1, 1, 1, 1)

  result <- project_to_box(x, lower, upper)

  expect_equal(result, c(-1, 0, 1, 1))
})

test_that("project_to_box handles infinite bounds", {
  x <- c(-10, 5, 20)
  lower <- c(-Inf, 0, -Inf)
  upper <- c(Inf, 10, 15)

  result <- project_to_box(x, lower, upper)

  expect_equal(result, c(-10, 5, 15))
})

test_that("project_to_box returns unchanged point if already feasible", {
  x <- c(0.5, -0.5, 0)
  lower <- c(-1, -1, -1)
  upper <- c(1, 1, 1)

  result <- project_to_box(x, lower, upper)

  expect_equal(result, x)
})

test_that("apply_box_constraints truncates step correctly", {
  # Current point near upper bound
  x_current <- c(0.8, 0)
  s <- c(0.5, 0.5)  # Would exceed upper bound
  lower <- c(-1, -1)
  upper <- c(1, 1)

  result <- apply_box_constraints(x_current, s, lower, upper)

  # Alpha should be (1 - 0.8) / 0.5 = 0.4
  # s_bounded should be 0.4 * c(0.5, 0.5) = c(0.2, 0.2)
  expect_equal(result$s_bounded, c(0.2, 0.2))
  expect_equal(result$x_new, c(1, 0.2))
})

test_that("apply_box_constraints handles negative steps", {
  # Current point near lower bound
  x_current <- c(-0.8, 0)
  s <- c(-0.5, 0.5)  # Would exceed lower bound in first component
  lower <- c(-1, -1)
  upper <- c(1, 1)

  result <- apply_box_constraints(x_current, s, lower, upper)

  # Alpha should be (-1 - (-0.8)) / (-0.5) = 0.4
  expect_equal(result$s_bounded, c(-0.2, 0.2))
  expect_equal(result$x_new, c(-1, 0.2))
})

test_that("apply_box_constraints returns full step when unconstrained", {
  x_current <- c(0, 0)
  s <- c(0.1, -0.1)
  lower <- c(-Inf, -Inf)
  upper <- c(Inf, Inf)

  result <- apply_box_constraints(x_current, s, lower, upper)

  # No constraints, should return full step
  expect_equal(result$s_bounded, s)
  expect_equal(result$x_new, x_current + s)
})

test_that("apply_box_constraints handles multiple active constraints", {
  x_current <- c(0.9, -0.9)
  s <- c(0.5, -0.5)  # Both would hit bounds
  lower <- c(-1, -1)
  upper <- c(1, 1)

  result <- apply_box_constraints(x_current, s, lower, upper)

  # Alpha for first: (1 - 0.9) / 0.5 = 0.2
  # Alpha for second: (-1 - (-0.9)) / (-0.5) = 0.2
  # Both give alpha = 0.2
  expect_equal(result$s_bounded, c(0.1, -0.1))
  expect_equal(result$x_new, c(1, -1))
})

test_that("validate_and_project_initial accepts feasible point", {
  x0 <- c(0.5, -0.5)
  lower <- c(-1, -1)
  upper <- c(1, 1)

  result <- validate_and_project_initial(x0, lower, upper)

  expect_equal(result, x0)
})

test_that("validate_and_project_initial projects infeasible point", {
  x0 <- c(2, -2)
  lower <- c(-1, -1)
  upper <- c(1, 1)

  expect_warning(
    result <- validate_and_project_initial(x0, lower, upper),
    "Initial point infeasible"
  )

  expect_equal(result, c(1, -1))
})

test_that("validate_and_project_initial detects bounds dimension mismatch", {
  x0 <- c(0, 0)
  lower <- c(-1, -1, -1)  # Wrong length
  upper <- c(1, 1)

  expect_error(
    validate_and_project_initial(x0, lower, upper),
    "Bounds dimension mismatch"
  )
})

test_that("validate_and_project_initial detects invalid bounds", {
  x0 <- c(0, 0)
  lower <- c(1, -1)  # lower[1] > upper[1]
  upper <- c(-1, 1)

  expect_error(
    validate_and_project_initial(x0, lower, upper),
    "Invalid bounds"
  )
})

test_that("apply_box_constraints handles zero step components", {
  x_current <- c(0, 0)
  s <- c(0, 0.5)  # Zero in first component
  lower <- c(-1, -1)
  upper <- c(1, 1)

  result <- apply_box_constraints(x_current, s, lower, upper)

  # Should return full step when components are zero
  expect_equal(result$s_bounded, s)
  expect_equal(result$x_new, c(0, 0.5))
})

test_that("apply_box_constraints projects at bounds", {
  # Start exactly at bound
  x_current <- c(1, -1)
  s <- c(0.1, -0.1)  # Would move beyond bounds
  lower <- c(-1, -1)
  upper <- c(1, 1)

  result <- apply_box_constraints(x_current, s, lower, upper)

  # Alpha should be 0 since already at bounds
  expect_equal(result$s_bounded, c(0, 0))
  expect_equal(result$x_new, c(1, -1))
})
