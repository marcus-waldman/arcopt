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


# =============================================================================
# Powell-Damped BFGS Tests
# =============================================================================

test_that("update_bfgs_powell applies damping when curvature violated", {
  b <- diag(3)
  s <- c(1, 0, 0)
  y <- c(-0.5, 0, 0)  # y's = -0.5 < 0, needs damping

  result <- update_bfgs_powell(b, s, y)

  # Should apply damping (theta < 1)
  expect_true(result$theta < 1.0)
  expect_true(result$theta > 0)
  expect_false(result$skipped)

  # Damped y should satisfy curvature condition
  ys_damped <- sum(result$y_damped * s)
  expect_true(ys_damped > 0)

  # B should be updated (different from input)
  expect_false(all(result$b == b))
})

test_that("update_bfgs_powell skips damping when curvature satisfied", {
  b <- diag(3)
  s <- c(1, 0, 0)
  y <- c(2, 0, 0)  # y's = 2 > 0, no damping needed

  result <- update_bfgs_powell(b, s, y)

  expect_equal(result$theta, 1.0)
  expect_equal(result$y_damped, y)
  expect_false(result$skipped)
})

test_that("update_bfgs_powell ensures minimum curvature", {
  b <- diag(3)
  s <- c(1, 0, 0)
  y <- c(0.1, 0, 0)  # y's = 0.1, but s'Bs = 1, so 0.1 < 0.2 * 1

  result <- update_bfgs_powell(b, s, y, damping_threshold = 0.2)

  # Should apply damping since y's < 0.2 * s'Bs
  expect_true(result$theta < 1.0)

  # Damped curvature should be at least threshold * s'Bs
  ys_damped <- sum(result$y_damped * s)
  sbs <- sum(s * b %*% s)
  expect_true(ys_damped >= 0.2 * sbs - 1e-10)
})

test_that("update_bfgs_powell maintains symmetry", {
  set.seed(42)
  a <- matrix(rnorm(9), 3, 3)
  b <- a %*% t(a) + diag(3)
  s <- c(0.5, -0.3, 0.8)
  y <- c(-0.2, 0.5, -0.4)  # Likely violates curvature

  result <- update_bfgs_powell(b, s, y)

  expect_equal(result$b, t(result$b), tolerance = 1e-10)
})

test_that("update_bfgs_powell handles edge case of zero s'Bs", {
  # Degenerate B with zero in direction of s
  b <- diag(c(0, 1, 1))
  s <- c(1, 0, 0)  # s'Bs = 0
  y <- c(1, 0, 0)

  result <- update_bfgs_powell(b, s, y)

  # Should skip since s'Bs is too small for damping
  expect_true(result$skipped)
})


# =============================================================================
# Hybrid BFGS/SR1 Routing Tests
# =============================================================================

test_that("update_hybrid uses BFGS when curvature satisfied", {
  b <- diag(3)
  s <- c(1, 0, 0)
  y <- c(2, 0, 0)  # y's = 2 > 0

  result <- update_hybrid(b, s, y)

  expect_equal(result$update_type, "bfgs")
  expect_false(result$skipped)
  expect_equal(result$skip_count, 0L)
  expect_false(result$restarted)
})

test_that("update_hybrid falls back to SR1 when BFGS fails", {
  b <- diag(3)
  s <- c(1, 0, 0)
  # y's <= 0 but SR1 denominator is good
  # y = (-0.5, 0.5, 0) gives y's = -0.5
  # r = y - Bs = (-0.5, 0.5, 0) - (1, 0, 0) = (-1.5, 0.5, 0)
  # denom = r's = -1.5, |denom| = 1.5 which is > skip_tol * ||r|| * ||s||
  y <- c(-0.5, 0.5, 0)

  result <- update_hybrid(b, s, y)

  expect_equal(result$update_type, "sr1")
  expect_false(result$skipped)
  expect_equal(result$skip_count, 0L)
})

test_that("update_hybrid falls back to Powell when BFGS and SR1 fail", {
  b <- diag(3)
  s <- c(1, 0, 0)
  # To make both BFGS and SR1 fail:
  # - BFGS fails when y's <= bfgs_tol
  # - SR1 fails when y ≈ Bs (so r ≈ 0)
  # Use high bfgs_tol and y = Bs exactly
  y <- c(1, 0, 0)  # y = Bs exactly, so r = 0

  # With high bfgs_tol, y's = 1 < bfgs_tol = 10, so BFGS fails
  # r = y - Bs = 0, so SR1 also fails (zero residual)
  result <- update_hybrid(b, s, y, bfgs_tol = 10)

  expect_equal(result$update_type, "powell")
  expect_false(result$skipped)
  expect_true(!is.na(result$theta))
})

test_that("update_hybrid produces same result as BFGS for well-conditioned case", {
  b <- diag(3)
  s <- c(1, 0, 0)
  y <- c(2, 0.1, 0)

  hybrid_result <- update_hybrid(b, s, y)
  bfgs_result <- update_bfgs(b, s, y)

  expect_equal(hybrid_result$update_type, "bfgs")
  expect_equal(hybrid_result$b, bfgs_result$b, tolerance = 1e-10)
})

test_that("update_hybrid tracks skip_count correctly", {
  b <- diag(3)
  s <- c(0, 0, 0)  # Zero step - should skip all methods
  y <- c(1, 0, 0)

  result <- update_hybrid(b, s, y, skip_count = 2L)

  expect_true(result$skipped)
  expect_equal(result$update_type, "skipped")
  # skip_count should remain unchanged for zero step (early exit)
  expect_equal(result$skip_count, 2L)
})

test_that("update_hybrid resets skip_count on successful update", {
  b <- diag(3)
  s <- c(1, 0, 0)
  y <- c(2, 0, 0)

  result <- update_hybrid(b, s, y, skip_count = 3L)

  expect_false(result$skipped)
  expect_equal(result$skip_count, 0L)  # Reset on success
})

test_that("update_hybrid maintains symmetry for all routing paths", {
  set.seed(42)
  a <- matrix(rnorm(9), 3, 3)
  b <- a %*% t(a) + diag(3)

  # Test BFGS path
  s1 <- c(1, 0.5, 0.3)
  y1 <- c(2, 1, 0.5)  # Positive curvature
  result1 <- update_hybrid(b, s1, y1)
  expect_equal(result1$b, t(result1$b), tolerance = 1e-10)

  # Test SR1 path
  s2 <- c(1, 0, 0)
  y2 <- c(-0.5, 0.5, 0)  # Negative curvature, good SR1 denominator
  result2 <- update_hybrid(b, s2, y2)
  expect_equal(result2$b, t(result2$b), tolerance = 1e-10)

  # Test Powell path
  s3 <- c(1, 0, 0)
  y3 <- c(1, 1, 0)  # BFGS boundary, SR1 perpendicular
  result3 <- update_hybrid(b, s3, y3, bfgs_tol = 1e-10)
  expect_equal(result3$b, t(result3$b), tolerance = 1e-10)
})

test_that("update_hybrid respects bfgs_tol parameter", {
  b <- diag(3)
  s <- c(1, 0, 0)
  y <- c(0.001, 0.5, 0)  # Very small positive curvature


  # With high bfgs_tol, should fall back to SR1
  result_high <- update_hybrid(b, s, y, bfgs_tol = 0.01)
  expect_equal(result_high$update_type, "sr1")

  # With low bfgs_tol, should use BFGS
  result_low <- update_hybrid(b, s, y, bfgs_tol = 1e-10)
  expect_equal(result_low$update_type, "bfgs")
})
