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
