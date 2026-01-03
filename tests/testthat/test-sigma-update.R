test_that("sigma update handles very successful step", {
  # rho >= eta2 should decrease sigma
  sigma_current <- 1.0
  rho <- 0.95  # > eta2 = 0.9
  eta1 <- 0.1
  eta2 <- 0.9
  gamma1 <- 0.5

  result <- update_sigma_cgt(sigma_current, rho, eta1, eta2, gamma1)

  expect_equal(result, gamma1 * sigma_current)
})

test_that("sigma update handles moderately successful step", {
  # eta1 <= rho < eta2 should keep sigma unchanged
  sigma_current <- 1.0
  rho <- 0.5  # Between eta1=0.1 and eta2=0.9
  eta1 <- 0.1
  eta2 <- 0.9

  result <- update_sigma_cgt(sigma_current, rho, eta1, eta2)

  expect_equal(result, sigma_current)
})

test_that("sigma update handles unsuccessful step", {
  # rho < eta1 should increase sigma
  sigma_current <- 1.0
  rho <- 0.05  # < eta1 = 0.1
  eta1 <- 0.1
  eta2 <- 0.9
  gamma2 <- 2.0

  result <- update_sigma_cgt(sigma_current, rho, eta1, eta2, gamma2 = gamma2)

  expect_equal(result, gamma2 * sigma_current)
})

test_that("sigma update respects sigma_min bound", {
  # Very successful step should not decrease below sigma_min
  sigma_current <- 1e-15
  rho <- 0.95
  sigma_min <- 1e-16

  result <- update_sigma_cgt(sigma_current, rho, sigma_min = sigma_min)

  expect_equal(result, sigma_min)
  expect_true(result >= sigma_min)
})

test_that("sigma update respects sigma_max bound", {
  # Unsuccessful step should not increase above sigma_max
  sigma_current <- 1e16  # Already at max
  rho <- 0.05  # Unsuccessful, would normally double
  sigma_max <- 1e16
  gamma2 <- 2.0

  result <- update_sigma_cgt(sigma_current, rho, gamma2 = gamma2, sigma_max = sigma_max)

  expect_equal(result, sigma_max)
  expect_true(result <= sigma_max)
})

test_that("sigma update handles rho exactly at eta1", {
  # rho = eta1 should be treated as successful (>= condition)
  sigma_current <- 1.0
  rho <- 0.1
  eta1 <- 0.1
  eta2 <- 0.9

  result <- update_sigma_cgt(sigma_current, rho, eta1, eta2)

  expect_equal(result, sigma_current)  # Should keep sigma
})

test_that("sigma update handles rho exactly at eta2", {
  # rho = eta2 should be treated as very successful (>= condition)
  sigma_current <- 1.0
  rho <- 0.9
  eta1 <- 0.1
  eta2 <- 0.9
  gamma1 <- 0.5

  result <- update_sigma_cgt(sigma_current, rho, eta1, eta2, gamma1)

  expect_equal(result, gamma1 * sigma_current)  # Should decrease
})

test_that("sigma update handles negative rho", {
  # Negative rho (worse than no progress) should increase sigma
  sigma_current <- 1.0
  rho <- -0.5
  gamma2 <- 2.0

  result <- update_sigma_cgt(sigma_current, rho, gamma2 = gamma2)

  expect_equal(result, gamma2 * sigma_current)
})

test_that("sigma update handles rho > 1", {
  # rho > 1 (better than predicted) should still decrease sigma
  sigma_current <- 1.0
  rho <- 1.2
  gamma1 <- 0.5

  result <- update_sigma_cgt(sigma_current, rho, gamma1 = gamma1)

  expect_equal(result, gamma1 * sigma_current)
})

test_that("sigma update uses default parameters correctly", {
  # Test with all defaults
  sigma_current <- 1.0
  rho <- 0.95

  result <- update_sigma_cgt(sigma_current, rho)

  # Default gamma1 = 0.5
  expect_equal(result, 0.5)
})

# Tests for Algorithm 2b (Interpolation-Enhanced)

test_that("interpolation sigma update handles very successful step", {
  # rho >= eta2 should use interpolation with CGT decrease
  sigma_current <- 1.0
  rho <- 0.95
  s <- c(0.1, 0.1)
  f_trial <- 1.0
  m_trial <- 1.01  # Small discrepancy
  eta1 <- 0.1
  eta2 <- 0.9
  gamma1 <- 0.5

  result <- update_sigma_interpolation(sigma_current, rho, s, f_trial, m_trial,
                                       eta1, eta2, gamma1)

  # Should be max(l_hat/2, gamma1*sigma, sigma_min)
  # Result should be finite and positive
  expect_true(is.finite(result))
  expect_true(result > 0)
  expect_true(result >= 1e-16)  # Above sigma_min
})

test_that("interpolation sigma update handles moderately successful step", {
  # eta1 <= rho < eta2 should keep sigma unchanged (same as CGT)
  sigma_current <- 1.0
  rho <- 0.5
  s <- c(0.1, 0.1)
  f_trial <- 1.0
  m_trial <- 1.0
  eta1 <- 0.1
  eta2 <- 0.9

  result <- update_sigma_interpolation(sigma_current, rho, s, f_trial, m_trial,
                                       eta1, eta2)

  expect_equal(result, sigma_current)
})

test_that("interpolation sigma update handles unsuccessful step", {
  # rho < eta1 should use interpolation with CGT increase
  sigma_current <- 1.0
  rho <- 0.05
  s <- c(0.1, 0.1)
  f_trial <- 1.0
  m_trial <- 0.5  # Large discrepancy
  eta1 <- 0.1
  eta2 <- 0.9
  gamma2 <- 2.0

  result <- update_sigma_interpolation(sigma_current, rho, s, f_trial, m_trial,
                                       eta1, eta2, gamma2 = gamma2)

  # Should be min(max(l_hat/2, gamma2*sigma), sigma_max)
  expect_true(result >= sigma_current)  # Should increase
})

test_that("interpolation sigma update handles zero step", {
  # Very small step should not cause division by zero
  sigma_current <- 1.0
  rho <- 0.95
  s <- c(0, 0)  # Zero step
  f_trial <- 1.0
  m_trial <- 1.0

  result <- update_sigma_interpolation(sigma_current, rho, s, f_trial, m_trial)

  # Should fall back to current sigma or default behavior
  expect_true(is.finite(result))
  expect_true(result > 0)
})

test_that("interpolation sigma update respects sigma_min", {
  sigma_current <- 1e-15
  rho <- 0.95
  s <- c(0.1, 0.1)
  f_trial <- 1.0
  m_trial <- 1.0
  sigma_min <- 1e-16

  result <- update_sigma_interpolation(sigma_current, rho, s, f_trial, m_trial,
                                       sigma_min = sigma_min)

  expect_true(result >= sigma_min)
})

test_that("interpolation sigma update respects sigma_max", {
  sigma_current <- 1e16
  rho <- 0.05  # Unsuccessful
  s <- c(0.1, 0.1)
  f_trial <- 1.0
  m_trial <- 0.0  # Large discrepancy
  sigma_max <- 1e16

  result <- update_sigma_interpolation(sigma_current, rho, s, f_trial, m_trial,
                                       sigma_max = sigma_max)

  expect_true(result <= sigma_max)
})

test_that("interpolation sigma estimates Lipschitz constant", {
  # Test that l_hat is computed correctly
  sigma_current <- 1.0
  rho <- 0.95
  s <- c(1, 0)  # Step of length 1
  f_trial <- 1.0
  m_trial <- 0.5  # Discrepancy of 0.5

  # L_hat = 2 * |f_trial - m_trial| / ||s||^3 = 2 * 0.5 / 1 = 1.0
  # For very successful: sigma_new = max(l_hat/2, gamma1*sigma, sigma_min)
  #                                 = max(0.5, 0.5, 1e-16) = 0.5
  result <- update_sigma_interpolation(sigma_current, rho, s, f_trial, m_trial)

  expect_equal(result, 0.5)
})
