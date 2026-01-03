test_that("convergence check detects max_iter criterion", {
  result <- check_convergence(
    g_current = c(0.1, 0.1),
    f_current = 10,
    f_previous = 11,
    x_current = c(1, 2),
    x_previous = c(0.9, 1.9),
    iter = 1000,
    max_iter = 1000
  )

  expect_true(result$converged)
  expect_equal(result$reason, "max_iter")
})

test_that("convergence check detects gradient_abs criterion", {
  result <- check_convergence(
    g_current = c(1e-10, -5e-11),  # max |g| = 1e-10 < 1e-8
    f_current = 10,
    f_previous = 11,
    x_current = c(1, 2),
    x_previous = c(0.9, 1.9),
    iter = 5,
    tol_grad = 1e-8
  )

  expect_true(result$converged)
  expect_equal(result$reason, "gradient_abs")
})

test_that("convergence check detects gradient_rel criterion", {
  result <- check_convergence(
    g_current = c(1e-5, -5e-6),  # max |g| = 1e-5
    f_current = 1e6,  # Large f value
    f_previous = 1e6 + 1,
    x_current = c(1, 2),
    x_previous = c(0.9, 1.9),
    iter = 5,
    tol_grad = 1e-6,  # Won't trigger abs criterion (1e-5 > 1e-6)
    tol_rel_grad = 1e-10  # But 1e-5 < 1e-10 * 1e6 = 1e-4
  )

  expect_true(result$converged)
  expect_equal(result$reason, "gradient_rel")
})

test_that("convergence check detects objective_abs criterion", {
  result <- check_convergence(
    g_current = c(0.1, 0.1),  # Large enough to not trigger gradient criterion
    f_current = 10.0,
    f_previous = 10.0 + 1e-13,  # Change < 1e-12
    x_current = c(1, 2),
    x_previous = c(0.9, 1.9),  # Large enough to not trigger parameter criterion
    iter = 5,
    tol_obj = 1e-12
  )

  expect_true(result$converged)
  expect_equal(result$reason, "objective_abs")
})

test_that("convergence check detects objective_rel criterion", {
  result <- check_convergence(
    g_current = c(10.0, 10.0),  # Large enough to avoid gradient criteria
    f_current = 1e6,
    f_previous = 1e6 + 1e-2,  # Change = 1e-2, relative = 1e-2/1e6 = 1e-8
    x_current = c(1, 2),
    x_previous = c(0.9, 1.9),  # Large enough to avoid parameter criterion
    iter = 5,
    tol_obj = 1e-20,  # Won't trigger abs criterion
    tol_rel_obj = 1e-7  # But 1e-8 < 1e-7
  )

  expect_true(result$converged)
  expect_equal(result$reason, "objective_rel")
})

test_that("convergence check detects parameter criterion", {
  result <- check_convergence(
    g_current = c(0.1, 0.1),
    f_current = 10,
    f_previous = 10.1,
    x_current = c(1.0, 2.0),
    x_previous = c(1.0 + 1e-10, 2.0 - 5e-11),  # max |dx| = 1e-10 < 1e-8
    iter = 5,
    tol_param = 1e-8
  )

  expect_true(result$converged)
  expect_equal(result$reason, "parameter")
})

test_that("convergence check returns FALSE when no criteria met", {
  result <- check_convergence(
    g_current = c(0.1, 0.1),
    f_current = 10,
    f_previous = 9,
    x_current = c(1, 2),
    x_previous = c(0.9, 1.9),
    iter = 5,
    max_iter = 1000,
    tol_grad = 1e-8,
    tol_rel_grad = 1e-6,
    tol_obj = 1e-12,
    tol_rel_obj = 1e-8,
    tol_param = 1e-8
  )

  expect_false(result$converged)
  expect_equal(result$reason, "")
})

test_that("convergence check handles first iteration correctly", {
  # First iteration: no previous values, should only check gradient and max_iter
  result <- check_convergence(
    g_current = c(0.1, 0.1),
    f_current = 10,
    f_previous = NA,
    x_current = c(1, 2),
    x_previous = NULL,
    iter = 0
  )

  expect_false(result$converged)
  expect_equal(result$reason, "")
})

test_that("convergence check uses infinity norm for gradient", {
  # Test that max(abs(g)) is used, not Euclidean norm
  result <- check_convergence(
    g_current = c(1e-10, 1e-10, 1e-10, 1e-7),  # max = 1e-7
    f_current = 10,
    f_previous = 11,
    x_current = rep(1, 4),
    x_previous = rep(0.9, 4),
    iter = 5,
    tol_grad = 1e-6  # 1e-7 < 1e-6, should converge
  )

  expect_true(result$converged)
  expect_equal(result$reason, "gradient_abs")
})

test_that("convergence check uses infinity norm for parameters", {
  # Test that max(abs(x - x_prev)) is used
  result <- check_convergence(
    g_current = c(0.1, 0.1, 0.1, 0.1),
    f_current = 10,
    f_previous = 9,
    x_current = c(1, 2, 3, 4),
    x_previous = c(1 + 1e-10, 2, 3, 4 + 1e-7),  # max change = 1e-7
    iter = 5,
    tol_param = 1e-6  # 1e-7 < 1e-6, should converge
  )

  expect_true(result$converged)
  expect_equal(result$reason, "parameter")
})
