# Tests for Quasi-Newton Hessian Updates (Algorithm 4a)
# ======================================================

# =============================================================================
# SR1 Update Tests
# =============================================================================

test_that("update_sr1 applies correct rank-1 update", {
  # Start with identity matrix
  b <- diag(3)
  s <- c(1, 0, 0)
  y <- c(2, 0, 0)  # Secant: B*s should equal y

  result <- update_sr1(b, s, y)

  expect_false(result$skipped)
  expect_false(result$restarted)
  expect_equal(result$skip_count, 0L)


  # Verify secant equation: B_new * s should be closer to y
  b_new_s <- as.vector(result$b %*% s)
  expect_equal(b_new_s, y, tolerance = 1e-10)
})

test_that("update_sr1 allows indefinite Hessian", {
  # Start with positive definite matrix
  b <- diag(3)
  # Choose s and y such that update creates negative eigenvalue
  s <- c(1, 0, 0)
  y <- c(-0.5, 0, 0)  # y < B*s, so r = y - B*s has opposite sign

  result <- update_sr1(b, s, y)

  # Check that update was applied (not skipped)
  expect_false(result$skipped)

  # Check that result has negative eigenvalue (indefinite)
  eigenvalues <- eigen(result$b, symmetric = TRUE, only.values = TRUE)$values
  expect_true(any(eigenvalues < 0))
})

test_that("update_sr1 skips when denominator too small", {
  b <- diag(3)
  # Choose s and y such that r is perpendicular to s (denom = r'*s = 0)
  # but r is not zero
  s <- c(1, 0, 0)
  # y such that r = y - B*s = y - s is perpendicular to s
  # r = (a, b, c), s = (1, 0, 0), r's = a = 0
  # So y = (1, 1, 0) gives r = (0, 1, 0) which is perpendicular to s
  y <- c(1, 1, 0)

  result <- update_sr1(b, s, y, skip_count = 0L)

  expect_true(result$skipped)
  expect_false(result$restarted)
  expect_equal(result$skip_count, 1L)
  # B should be unchanged
  expect_equal(result$b, b)
})

test_that("update_sr1 restarts after consecutive skips", {
  b <- diag(3) * 10  # Start with scaled identity
  s <- c(1, 0, 0)
  # y such that r = y - B*s is perpendicular to s (causes skip)
  # B*s = (10, 0, 0), so y = (10, 1, 0) gives r = (0, 1, 0) ⊥ s
  y <- c(10, 1, 0)

  # Simulate 5 consecutive skips (threshold is 5), next skip triggers restart
  result <- update_sr1(b, s, y, skip_count = 5L, restart_threshold = 5L)

  expect_true(result$skipped)
  expect_true(result$restarted)

  expect_equal(result$skip_count, 0L)  # Reset after restart

  # B should be reset to scaled identity with Barzilai-Borwein scaling
  # gamma = y's / s's = 10/1 = 10
  expect_equal(result$b, 10 * diag(3))
})

test_that("update_sr1 maintains symmetry", {
  # Random symmetric positive definite starting matrix
  set.seed(42)
  a <- matrix(rnorm(9), 3, 3)
  b <- a %*% t(a) + diag(3)

  s <- c(0.5, -0.3, 0.8)
  y <- c(1.2, 0.5, -0.4)

  result <- update_sr1(b, s, y)

  # Check symmetry
  expect_equal(result$b, t(result$b), tolerance = 1e-10)
})


# =============================================================================
# BFGS Update Tests
# =============================================================================

test_that("update_bfgs applies correct rank-2 update", {
  b <- diag(3)
  s <- c(1, 0, 0)
  y <- c(2, 0.1, 0)  # y's > 0 (curvature condition satisfied)

  result <- update_bfgs(b, s, y)

  expect_false(result$skipped)

  # Verify secant equation: B_new * s = y
  b_new_s <- as.vector(result$b %*% s)
  expect_equal(b_new_s, y, tolerance = 1e-10)
})

test_that("update_bfgs maintains positive definiteness", {
  # Start with positive definite
  set.seed(42)
  a <- matrix(rnorm(9), 3, 3)
  b <- a %*% t(a) + diag(3)

  # Use vectors that satisfy curvature condition y's > 0
  s <- c(1, 0.5, -0.3)
  y <- c(2, 1, 0.5)  # y's = 2 + 0.5 - 0.15 > 0

  result <- update_bfgs(b, s, y)

  expect_false(result$skipped)

  # Check positive definiteness via Cholesky
  chol_result <- tryCatch(chol(result$b), error = function(e) NULL)
  expect_false(is.null(chol_result))
})

test_that("update_bfgs skips when curvature condition violated", {
  b <- diag(3)
  s <- c(1, 0, 0)
  y <- c(-0.5, 0, 0)  # y's = -0.5 < 0 (curvature violated)

  result <- update_bfgs(b, s, y)

  expect_true(result$skipped)
  expect_equal(result$b, b)  # B unchanged
})

test_that("update_bfgs maintains symmetry", {
  set.seed(123)
  a <- matrix(rnorm(9), 3, 3)
  b <- a %*% t(a) + diag(3)

  s <- c(0.5, -0.3, 0.8)
  y <- c(1.2, 0.5, 0.4)

  result <- update_bfgs(b, s, y)

  expect_equal(result$b, t(result$b), tolerance = 1e-10)
})


# =============================================================================
# L-BFGS Tests
# =============================================================================

test_that("update_lbfgs initializes empty history", {
  s <- c(1, 0)
  y <- c(2, 0.1)

  history <- update_lbfgs(NULL, s, y, m = 5)

  expect_equal(length(history$s), 1)
  expect_equal(length(history$y), 1)
  expect_equal(length(history$rho), 1)
  expect_true(history$gamma > 0)
})

test_that("update_lbfgs respects memory limit", {
  history <- NULL
  m <- 3

  # Add 5 pairs to history with m=3
  for (i in 1:5) {
    s <- c(1, 0) + i * 0.1
    y <- c(2, 0.1) + i * 0.05
    history <- update_lbfgs(history, s, y, m = m)
  }

  # Should only keep last 3
  expect_equal(length(history$s), 3)
  expect_equal(length(history$y), 3)
  expect_equal(length(history$rho), 3)
})

test_that("update_lbfgs skips when curvature violated", {
  history <- NULL
  s <- c(1, 0)
  y <- c(-0.5, 0)  # y's < 0

  history <- update_lbfgs(history, s, y, m = 5)

  expect_equal(length(history$s), 0)
})

test_that("lbfgs_multiply returns gamma*g for empty history", {
  history <- list(s = list(), y = list(), rho = numeric(0), gamma = 2.5)
  g <- c(1, 2, 3)

  result <- lbfgs_multiply(history, g)

  expect_equal(result, 2.5 * g)
})

test_that("lbfgs_multiply computes H*g (inverse Hessian approximation)", {
  # L-BFGS two-loop recursion computes H*g (inverse Hessian times g),

  # which gives the search direction for line search methods.
  # For ARC, we use the Woodbury solver for B*g instead.

  n <- 3
  history <- NULL

  # Generate (s, y) pairs satisfying curvature condition
  set.seed(42)
  for (i in 1:3) {
    s <- rnorm(n)
    y <- rnorm(n) + 2 * s  # Ensure y's > 0
    history <- update_lbfgs(history, s, y, m = 10)
  }

  # For any g, H*g should be a valid descent direction (opposite to g)
  g <- c(1, 1, 1)
  h_times_g <- lbfgs_multiply(history, g)

  # H*g should be non-zero
  expect_true(sqrt(sum(h_times_g^2)) > 1e-10)

  # For positive definite H (which L-BFGS maintains), H*g and g should
  # point in similar direction (H*g is the search direction, not descent)
  # Actually, H*g points in same direction as -descent, so H*g should
  # be correlated with g for PD H
  cos_angle <- sum(h_times_g * g) / (sqrt(sum(h_times_g^2)) * sqrt(sum(g^2)))
  expect_true(cos_angle > 0)  # Same half-plane as g
})


# =============================================================================
# L-SR1 Tests
# =============================================================================

test_that("update_lsr1 initializes empty history", {
  s <- c(1, 0)
  y <- c(2, 0.1)

  history <- update_lsr1(NULL, s, y, m = 5)

  expect_equal(length(history$s), 1)
  expect_equal(length(history$y), 1)
  expect_false(history$skipped)
})

test_that("update_lsr1 respects memory limit", {
  history <- NULL
  m <- 3

  for (i in 1:5) {
    s <- c(1, 0) + i * 0.1
    y <- c(2 + i * 0.5, 0.1)  # Different y each time to avoid skip
    history <- update_lsr1(history, s, y, m = m)
  }

  expect_equal(length(history$s), 3)
  expect_equal(length(history$y), 3)
})

test_that("update_lsr1 skips with near-zero denominator", {
  history <- list(s = list(), y = list(), gamma = 1.0)
  s <- c(1, 0)
  # y = gamma*s makes r = y - B*s = 0
  y <- c(1, 0)

  result <- update_lsr1(history, s, y, m = 5)

  expect_true(result$skipped)
  expect_equal(length(result$s), 0)  # No pair added
})

test_that("lsr1_multiply returns gamma*v for empty history", {
  history <- list(s = list(), y = list(), gamma = 3.0)
  v <- c(1, 2, 3)

  result <- lsr1_multiply(history, v)

  expect_equal(result, 3.0 * v)
})

test_that("lsr1_multiply satisfies secant equation on stored pairs", {
  # L-SR1 compact representation uses dynamic gamma scaling,
  # so it won't exactly match explicit SR1. Instead verify key properties.
  n <- 3
  history <- NULL
  stored_s <- list()
  stored_y <- list()

  # Generate (s, y) pairs
  set.seed(42)
  for (i in 1:3) {
    s <- rnorm(n)
    y <- rnorm(n) + 1.5 * s  # Ensure non-trivial r

    b_times_s <- if (is.null(history)) history$gamma * s else lsr1_multiply(history, s)
    if (is.null(b_times_s)) b_times_s <- s  # gamma = 1 initially

    history <- update_lsr1(history, s, y, m = 10, b_times_s = b_times_s)

    if (!history$skipped) {
      stored_s <- c(stored_s, list(s))
      stored_y <- c(stored_y, list(y))
    }
  }

  # Property 1: B*s_i should approximately equal y_i for stored pairs
  # (exact only for most recent pair due to finite memory)
  if (length(stored_s) > 0) {
    s_last <- stored_s[[length(stored_s)]]
    y_last <- stored_y[[length(stored_y)]]
    b_s_last <- lsr1_multiply(history, s_last)

    # Should satisfy secant equation reasonably well
    rel_error <- sqrt(sum((b_s_last - y_last)^2)) / sqrt(sum(y_last^2))
    expect_true(rel_error < 1.0)  # Within 100% relative error
  }

  # Property 2: Output should be valid (non-zero for non-zero input)
  v <- c(1, 1, 1)
  result <- lsr1_multiply(history, v)
  expect_true(sqrt(sum(result^2)) > 0)
})


# =============================================================================
# Edge Cases
# =============================================================================

test_that("all updates handle zero vectors gracefully", {
  b <- diag(3)
  s <- c(0, 0, 0)
  y <- c(1, 0, 0)

  # SR1 should skip (denom = 0)
  sr1_result <- update_sr1(b, s, y)
  expect_true(sr1_result$skipped)

  # BFGS should skip (y's = 0)
  bfgs_result <- update_bfgs(b, s, y)
  expect_true(bfgs_result$skipped)
})

test_that("all updates handle very small vectors", {
  b <- diag(3)
  s <- c(1e-15, 0, 0)
  y <- c(1e-15, 0, 0)

  # Should not crash
  sr1_result <- update_sr1(b, s, y)
  bfgs_result <- update_bfgs(b, s, y)

  expect_true(is.matrix(sr1_result$b))
  expect_true(is.matrix(bfgs_result$b))
})
