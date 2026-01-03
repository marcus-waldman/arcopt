test_that("stagnation detection identifies small consecutive steps", {
  # Create history with 5 consecutive small steps
  step_norms <- c(0.1, 0.05, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13)
  f_values <- c(10, 9, 8, 7, 6, 5, 4)

  result <- detect_stagnation(
    step_norms = step_norms,
    f_values = f_values,
    max_stagnant = 5
  )

  expect_equal(result, "stop")
})

test_that("stagnation detection identifies small consecutive function changes", {
  # Create history with 5 consecutive small function changes
  step_norms <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  f_values <- c(10, 10 + 1e-15, 10 + 2e-15, 10 + 3e-15, 10 + 4e-15, 10 + 5e-15)

  result <- detect_stagnation(
    step_norms = step_norms,
    f_values = f_values,
    max_stagnant = 5
  )

  expect_equal(result, "stop")
})

test_that("stagnation detection continues when progress is made", {
  # Normal optimization progress
  step_norms <- c(0.5, 0.3, 0.2, 0.1, 0.05)
  f_values <- c(10, 8, 6, 4, 2)

  result <- detect_stagnation(
    step_norms = step_norms,
    f_values = f_values,
    max_stagnant = 5
  )

  expect_equal(result, "continue")
})

test_that("stagnation detection suggests Hessian refresh for quasi-Newton", {
  # Stagnation with quasi-Newton
  step_norms <- c(0.1, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13)
  f_values <- c(10, 9, 8, 7, 6, 5)

  result <- detect_stagnation(
    step_norms = step_norms,
    f_values = f_values,
    max_stagnant = 5,
    is_qn = TRUE,
    already_refreshed = FALSE
  )

  expect_equal(result, "refresh_hessian")
})

test_that("stagnation detection stops after Hessian refresh attempted", {
  # Stagnation even after refresh
  step_norms <- c(0.1, 1e-13, 1e-13, 1e-13, 1e-13, 1e-13)
  f_values <- c(10, 9, 8, 7, 6, 5)

  result <- detect_stagnation(
    step_norms = step_norms,
    f_values = f_values,
    max_stagnant = 5,
    is_qn = TRUE,
    already_refreshed = TRUE
  )

  expect_equal(result, "stop")
})

test_that("stagnation detection handles insufficient history", {
  # Not enough iterations to detect stagnation
  step_norms <- c(1e-13, 1e-13)
  f_values <- c(10, 10 + 1e-15)

  result <- detect_stagnation(
    step_norms = step_norms,
    f_values = f_values,
    max_stagnant = 5
  )

  expect_equal(result, "continue")
})

test_that("stagnation detection resets after non-stagnant iteration", {
  # Had stagnation, then progress, then stagnation again (but not consecutive)
  step_norms <- c(1e-13, 1e-13, 0.1, 1e-13, 1e-13)
  f_values <- c(10, 10, 10, 9, 9)

  result <- detect_stagnation(
    step_norms = step_norms,
    f_values = f_values,
    max_stagnant = 5
  )

  # Only 2 consecutive small steps at end, not 5
  expect_equal(result, "continue")
})

test_that("stagnation detection uses custom tolerances", {
  step_norms <- c(0.1, 0.01, 0.01, 0.01, 0.01, 0.01)
  f_values <- c(10, 9, 8, 7, 6, 5)

  # With loose tolerance, these are "small"
  result <- detect_stagnation(
    step_norms = step_norms,
    f_values = f_values,
    tol_step = 0.05,  # Larger tolerance
    max_stagnant = 5
  )

  expect_equal(result, "stop")
})

test_that("check_finite detects NaN in function value", {
  f_value <- NaN
  g_value <- c(1, 2, 3)

  result <- check_finite(f_value, g_value)

  expect_false(result)
})

test_that("check_finite detects Inf in function value", {
  f_value <- Inf
  g_value <- c(1, 2, 3)

  result <- check_finite(f_value, g_value)

  expect_false(result)
})

test_that("check_finite detects NaN in gradient", {
  f_value <- 10.5
  g_value <- c(1, NaN, 3)

  result <- check_finite(f_value, g_value)

  expect_false(result)
})

test_that("check_finite detects Inf in gradient", {
  f_value <- 10.5
  g_value <- c(1, 2, -Inf)

  result <- check_finite(f_value, g_value)

  expect_false(result)
})

test_that("check_finite accepts valid finite values", {
  f_value <- 10.5
  g_value <- c(1, 2, 3)

  result <- check_finite(f_value, g_value)

  expect_true(result)
})

test_that("check_finite handles zero and negative values", {
  f_value <- -100.5
  g_value <- c(0, -1, -2)

  result <- check_finite(f_value, g_value)

  expect_true(result)
})
