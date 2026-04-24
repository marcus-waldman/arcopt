test_that("init_flat_ridge_state has expected structure", {
  state <- init_flat_ridge_state(window = 5)
  expect_equal(state$window, 5L)
  expect_equal(length(state$sigma), 0)
  expect_equal(length(state$rho), 0)
  expect_equal(length(state$g_inf), 0)
  expect_equal(length(state$lambda_min), 0)
})

test_that("init_flat_ridge_state defaults to window = 10", {
  state <- init_flat_ridge_state()
  expect_equal(state$window, 10L)
})

test_that("init_flat_ridge_state rejects invalid window", {
  expect_error(init_flat_ridge_state(window = 0))
  expect_error(init_flat_ridge_state(window = -1))
  expect_error(init_flat_ridge_state(window = "ten"))
})

test_that("update_flat_ridge_state accumulates diagnostics", {
  state <- init_flat_ridge_state(window = 3)
  state <- update_flat_ridge_state(state, sigma = 1e-6, rho = 0.99,
                                   g_inf = 1.5, lambda_min = 0.001)
  expect_equal(length(state$sigma), 1)
  expect_equal(state$sigma, 1e-6)
  expect_equal(state$rho, 0.99)
  expect_equal(state$g_inf, 1.5)
  expect_equal(state$lambda_min, 0.001)
})

test_that("update_flat_ridge_state trims window to fixed length", {
  state <- init_flat_ridge_state(window = 3)
  for (i in 1:5) {
    state <- update_flat_ridge_state(state, sigma = i, rho = i,
                                     g_inf = i, lambda_min = i)
  }
  expect_equal(length(state$sigma), 3)
  expect_equal(state$sigma, c(3, 4, 5)) # oldest dropped
  expect_equal(state$rho, c(3, 4, 5))
  expect_equal(state$g_inf, c(3, 4, 5))
  expect_equal(state$lambda_min, c(3, 4, 5))
})

test_that("check_flat_ridge_trigger returns FALSE before window is full", {
  state <- init_flat_ridge_state(window = 10)
  # Fill only 9 entries (window requires 10)
  for (i in 1:9) {
    state <- update_flat_ridge_state(state, sigma = 1e-6, rho = 0.99,
                                     g_inf = 1.0, lambda_min = 1e-4)
  }
  expect_false(check_flat_ridge_trigger(state, sigma_min = 1e-6))
})

test_that("check_flat_ridge_trigger fires when all four signals align", {
  state <- init_flat_ridge_state(window = 10)
  for (i in 1:10) {
    state <- update_flat_ridge_state(
      state,
      sigma = 1e-6,          # pinned
      rho = 1.0,             # perfect model
      g_inf = 1.5,           # stagnant (all equal -> ratio 1.0)
      lambda_min = 1e-4      # PD, < tol_ridge
    )
  }
  expect_true(check_flat_ridge_trigger(state, sigma_min = 1e-6))
})

test_that("trigger does NOT fire at a normal local minimum (||g|| tiny)", {
  state <- init_flat_ridge_state(window = 10)
  # Gradient below g_inf_floor default (1e-6)
  for (i in 1:10) {
    state <- update_flat_ridge_state(
      state, sigma = 1e-6, rho = 1.0, g_inf = 1e-9, lambda_min = 1e-4
    )
  }
  expect_false(check_flat_ridge_trigger(state, sigma_min = 1e-6))
})

test_that("trigger DOES fire at stuck-indefinite-saddle", {
  # Revised Signal 4 semantics (2026-04-23): detector now fires on
  # BOTH classical flat-ridge (small positive lambda_min) AND stuck-
  # indefinite-saddle (negative lambda_min), since cubic regularization
  # loses grip in either regime when the other signals are present.
  state <- init_flat_ridge_state(window = 10)
  for (i in 1:10) {
    state <- update_flat_ridge_state(
      state, sigma = 1e-6, rho = 1.0, g_inf = 1.5, lambda_min = -0.5
    )
  }
  expect_true(check_flat_ridge_trigger(state, sigma_min = 1e-6))
})

test_that("trigger does NOT fire when lambda_min above tol_ridge", {
  state <- init_flat_ridge_state(window = 10)
  for (i in 1:10) {
    state <- update_flat_ridge_state(
      state, sigma = 1e-6, rho = 1.0, g_inf = 1.5, lambda_min = 0.5
    )
  }
  expect_false(check_flat_ridge_trigger(state, sigma_min = 1e-6))
})

test_that("trigger does NOT fire when sigma is still oscillating", {
  state <- init_flat_ridge_state(window = 10)
  # One iteration has sigma well above floor
  for (i in 1:10) {
    sigma_i <- if (i == 5) 1.0 else 1e-6
    state <- update_flat_ridge_state(
      state, sigma = sigma_i, rho = 1.0, g_inf = 1.5, lambda_min = 1e-4
    )
  }
  expect_false(check_flat_ridge_trigger(state, sigma_min = 1e-6))
})

test_that("trigger does NOT fire when rho is far from 1", {
  state <- init_flat_ridge_state(window = 10)
  for (i in 1:10) {
    state <- update_flat_ridge_state(
      state, sigma = 1e-6, rho = 0.3, g_inf = 1.5, lambda_min = 1e-4
    )
  }
  expect_false(check_flat_ridge_trigger(state, sigma_min = 1e-6))
})

test_that("trigger does NOT fire when gradient is still decreasing", {
  state <- init_flat_ridge_state(window = 10)
  # g_inf decays geometrically by 0.5 per iter -> latest/oldest = 0.5^9 ~ 2e-3
  for (i in 1:10) {
    state <- update_flat_ridge_state(
      state, sigma = 1e-6, rho = 1.0,
      g_inf = 1.5 * 0.5^(i - 1), lambda_min = 1e-4
    )
  }
  expect_false(check_flat_ridge_trigger(state, sigma_min = 1e-6))
})

test_that("trigger fires with >= 0.9 ratio but not with strictly < 0.9", {
  # Ratio exactly at 0.9: do NOT fire (strict inequality guards against false
  # triggers on a single lucky window)
  make_state <- function(ratio) {
    state <- init_flat_ridge_state(window = 10)
    g_old <- 1.0
    g_new <- ratio * g_old
    gs <- seq(g_old, g_new, length.out = 10)
    for (i in seq_along(gs)) {
      state <- update_flat_ridge_state(
        state, sigma = 1e-6, rho = 1.0,
        g_inf = gs[i], lambda_min = 1e-4
      )
    }
    state
  }
  # Just over 0.9 ratio -> stagnant -> fires
  expect_true(check_flat_ridge_trigger(make_state(0.95), sigma_min = 1e-6))
  # Exactly 0.9 -> does NOT fire (uses strict > 0.9)
  expect_false(check_flat_ridge_trigger(make_state(0.9), sigma_min = 1e-6))
  # Below 0.9 -> does NOT fire
  expect_false(check_flat_ridge_trigger(make_state(0.5), sigma_min = 1e-6))
})

test_that("trigger handles NA / non-finite diagnostics safely", {
  state <- init_flat_ridge_state(window = 10)
  for (i in 1:10) {
    lm_i <- if (i == 5) NA_real_ else 1e-4
    state <- update_flat_ridge_state(
      state, sigma = 1e-6, rho = 1.0, g_inf = 1.5, lambda_min = lm_i
    )
  }
  # Latest lambda_min is 1e-4 (finite) but earlier NA would still be tolerated
  # because only latest lambda_min is checked; however the window is length 10
  # and latest is i=10 (finite), so trigger could fire. Let's ensure it does.
  expect_true(check_flat_ridge_trigger(state, sigma_min = 1e-6))

  # But if the *latest* lambda_min is non-finite, must NOT fire
  state2 <- update_flat_ridge_state(
    state, sigma = 1e-6, rho = 1.0, g_inf = 1.5, lambda_min = NaN
  )
  expect_false(check_flat_ridge_trigger(state2, sigma_min = 1e-6))
})

test_that("custom tol_ridge changes trigger behavior", {
  state <- init_flat_ridge_state(window = 10)
  # lambda_min = 5e-3: above default tol_ridge (1e-3), below custom (1e-2)
  for (i in 1:10) {
    state <- update_flat_ridge_state(
      state, sigma = 1e-6, rho = 1.0, g_inf = 1.5, lambda_min = 5e-3
    )
  }
  expect_false(check_flat_ridge_trigger(state, sigma_min = 1e-6))
  expect_true(check_flat_ridge_trigger(state, sigma_min = 1e-6,
                                       tol_ridge = 1e-2))
})

test_that("custom g_inf_floor prevents or permits trigger", {
  state <- init_flat_ridge_state(window = 10)
  # ||g||_inf = 1e-4: above default floor (1e-6), below custom (1e-3)
  for (i in 1:10) {
    state <- update_flat_ridge_state(
      state, sigma = 1e-6, rho = 1.0, g_inf = 1e-4, lambda_min = 1e-4
    )
  }
  expect_true(check_flat_ridge_trigger(state, sigma_min = 1e-6))
  expect_false(check_flat_ridge_trigger(state, sigma_min = 1e-6,
                                        g_inf_floor = 1e-3))
})
