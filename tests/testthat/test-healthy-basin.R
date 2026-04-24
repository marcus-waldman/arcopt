test_that("init_healthy_basin_state has expected structure", {
  s <- init_healthy_basin_state(window = 4)
  expect_equal(s$window, 4L)
  expect_equal(length(s$used_newton), 0)
  expect_equal(length(s$rho), 0)
  expect_equal(length(s$lambda_min), 0)
  expect_equal(length(s$g_inf), 0)
})

test_that("init_healthy_basin_state defaults to window = 5", {
  expect_equal(init_healthy_basin_state()$window, 5L)
})

test_that("init_healthy_basin_state rejects invalid window", {
  expect_error(init_healthy_basin_state(window = 0))
  expect_error(init_healthy_basin_state(window = -2))
  expect_error(init_healthy_basin_state(window = "five"))
})

test_that("update_healthy_basin_state accumulates and trims", {
  s <- init_healthy_basin_state(window = 3)
  for (i in 1:5) {
    s <- update_healthy_basin_state(s,
                                    used_newton = TRUE,
                                    rho = 0.95,
                                    lambda_min = 1,
                                    g_inf = 1 * 0.5^i)
  }
  expect_equal(length(s$used_newton), 3)
  expect_equal(s$g_inf, 1 * 0.5^(3:5))
})

test_that("trigger FALSE before window is full", {
  s <- init_healthy_basin_state(window = 5)
  for (i in 1:4) {
    s <- update_healthy_basin_state(s, TRUE, 0.99, 1,
                                    g_inf = 1 * 0.3^i)
  }
  expect_false(check_healthy_basin_trigger(s))
})

test_that("trigger TRUE when all 5 signals align over full window", {
  s <- init_healthy_basin_state(window = 5)
  # Gradient decaying by factor 0.3 each iteration (superlinear)
  gs <- 1 * 0.3^(0:4)
  for (i in 1:5) {
    s <- update_healthy_basin_state(s,
                                    used_newton = TRUE,
                                    rho = 0.95,
                                    lambda_min = 1,
                                    g_inf = gs[i])
  }
  expect_true(check_healthy_basin_trigger(s))
})

test_that("trigger FALSE if ANY Newton step NOT used in window", {
  s <- init_healthy_basin_state(window = 5)
  gs <- 1 * 0.3^(0:4)
  used <- c(TRUE, TRUE, FALSE, TRUE, TRUE)
  for (i in 1:5) {
    s <- update_healthy_basin_state(s,
                                    used_newton = used[i],
                                    rho = 0.95,
                                    lambda_min = 1,
                                    g_inf = gs[i])
  }
  expect_false(check_healthy_basin_trigger(s))
})

test_that("trigger FALSE if rho dips below rho_polish anywhere in window", {
  s <- init_healthy_basin_state(window = 5)
  gs <- 1 * 0.3^(0:4)
  rhos <- c(0.95, 0.95, 0.5, 0.95, 0.95)
  for (i in 1:5) {
    s <- update_healthy_basin_state(s,
                                    used_newton = TRUE,
                                    rho = rhos[i],
                                    lambda_min = 1,
                                    g_inf = gs[i])
  }
  expect_false(check_healthy_basin_trigger(s))
})

test_that("trigger FALSE if lambda_min dips below threshold anywhere", {
  s <- init_healthy_basin_state(window = 5)
  gs <- 1 * 0.3^(0:4)
  lmins <- c(1, 1, 1e-4, 1, 1)
  for (i in 1:5) {
    s <- update_healthy_basin_state(s,
                                    used_newton = TRUE,
                                    rho = 0.95,
                                    lambda_min = lmins[i],
                                    g_inf = gs[i])
  }
  expect_false(check_healthy_basin_trigger(s))
})

test_that("trigger FALSE if gradient decays too slowly", {
  s <- init_healthy_basin_state(window = 5)
  # Decay ratio 0.7 per iter — slower than default 0.5 threshold
  gs <- 1 * 0.7^(0:4)
  for (i in 1:5) {
    s <- update_healthy_basin_state(s,
                                    used_newton = TRUE,
                                    rho = 0.95,
                                    lambda_min = 1,
                                    g_inf = gs[i])
  }
  expect_false(check_healthy_basin_trigger(s))
})

test_that("trigger FALSE when ||g|| already below convergence floor", {
  s <- init_healthy_basin_state(window = 5)
  # Start at 1e-10 — already converged, all other signals trivially pass
  gs <- 1e-10 * 0.3^(0:4)
  for (i in 1:5) {
    s <- update_healthy_basin_state(s,
                                    used_newton = TRUE,
                                    rho = 0.95,
                                    lambda_min = 1,
                                    g_inf = gs[i])
  }
  expect_false(check_healthy_basin_trigger(s))
})

test_that("trigger handles NA/NaN/Inf safely", {
  s <- init_healthy_basin_state(window = 5)
  gs <- 1 * 0.3^(0:4)
  rhos <- c(0.95, 0.95, NA_real_, 0.95, 0.95)
  for (i in 1:5) {
    s <- update_healthy_basin_state(s,
                                    used_newton = TRUE,
                                    rho = rhos[i],
                                    lambda_min = 1,
                                    g_inf = gs[i])
  }
  expect_false(check_healthy_basin_trigger(s))
})

test_that("trigger respects custom rho_polish threshold", {
  s <- init_healthy_basin_state(window = 5)
  gs <- 1 * 0.3^(0:4)
  for (i in 1:5) {
    s <- update_healthy_basin_state(s,
                                    used_newton = TRUE,
                                    rho = 0.75,
                                    lambda_min = 1,
                                    g_inf = gs[i])
  }
  # Default rho_polish = 0.9 -> FALSE (rho 0.75 < 0.9)
  expect_false(check_healthy_basin_trigger(s))
  # Looser threshold 0.7 -> TRUE
  expect_true(check_healthy_basin_trigger(s, rho_polish = 0.7))
})

test_that("trigger respects custom lambda_min_polish threshold", {
  s <- init_healthy_basin_state(window = 5)
  gs <- 1 * 0.3^(0:4)
  for (i in 1:5) {
    s <- update_healthy_basin_state(s,
                                    used_newton = TRUE,
                                    rho = 0.95,
                                    lambda_min = 5e-4,
                                    g_inf = gs[i])
  }
  # Default threshold = 1e-3 -> FALSE
  expect_false(check_healthy_basin_trigger(s))
  # Looser threshold 1e-4 -> TRUE
  expect_true(check_healthy_basin_trigger(s, lambda_min_polish = 1e-4))
})

test_that("trigger respects custom g_decay_polish threshold", {
  s <- init_healthy_basin_state(window = 5)
  # Decay 0.7 per iter
  gs <- 1 * 0.7^(0:4)
  for (i in 1:5) {
    s <- update_healthy_basin_state(s,
                                    used_newton = TRUE,
                                    rho = 0.95,
                                    lambda_min = 1,
                                    g_inf = gs[i])
  }
  expect_false(check_healthy_basin_trigger(s))
  expect_true(check_healthy_basin_trigger(s, g_decay_polish = 0.8))
})

test_that("trigger FALSE if any g_inf is zero or negative (numerical safety)", {
  s <- init_healthy_basin_state(window = 5)
  gs <- c(1, 0.5, 0, 0.125, 0.0625)
  for (i in 1:5) {
    s <- update_healthy_basin_state(s,
                                    used_newton = TRUE,
                                    rho = 0.95,
                                    lambda_min = 1,
                                    g_inf = gs[i])
  }
  expect_false(check_healthy_basin_trigger(s))
})
