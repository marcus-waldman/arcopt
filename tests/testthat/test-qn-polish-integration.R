test_that("qn_polish is disabled by default", {
  rosen <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  rosen_gr <- function(x) {
    c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
      200 * (x[2] - x[1]^2))
  }
  rosen_hess <- function(x) {
    matrix(c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
             -400 * x[1], 200), 2, 2)
  }

  res <- arcopt(c(-1.2, 1), rosen, rosen_gr, rosen_hess,
                control = list(maxit = 200, trace = 0))

  expect_true(res$converged)
  expect_equal(res$diagnostics$qn_polish_switches, 0L)
  expect_equal(res$diagnostics$qn_polish_reverts, 0L)
  expect_equal(res$diagnostics$solver_mode_final, "cubic")
})

test_that("arcopt result includes qn_polish fields", {
  fn <- function(x) 0.5 * sum(x^2)
  gr <- function(x) x
  hess_fn <- function(x) diag(length(x))

  res <- arcopt(c(1, 1), fn, gr, hess_fn,
                control = list(maxit = 50, trace = 0))

  for (f in c("qn_polish_switches", "qn_polish_reverts",
              "hess_evals_at_polish_switch")) {
    expect_true(f %in% names(res$diagnostics),
                info = paste("missing diagnostics field:", f))
  }
  expect_equal(res$diagnostics$qn_polish_switches, 0L)
  expect_equal(res$diagnostics$qn_polish_reverts, 0L)
})

test_that("qn_polish fires on smooth convex with quartic tail", {
  # f(x) = sum(0.5 * x^2 + x^4 / 12), PD everywhere, well-conditioned in
  # basin. Newton takes ~8 iterations from x0=10*rep(1,n); cubic/Newton
  # transitional iterations don't occur because H is PD throughout.
  fn <- function(x) sum(0.5 * x^2 + x^4 / 12)
  gr <- function(x) x + x^3 / 3
  hess_fn <- function(x) diag(1 + x^2)

  x0 <- rep(10, 5)
  res <- arcopt(x0, fn, gr, hess_fn,
                control = list(maxit = 50, trace = 0,
                               qn_polish_enabled = TRUE))

  expect_true(res$converged)
  expect_true(max(abs(res$par)) < 1e-4)
  expect_gte(res$diagnostics$qn_polish_switches, 1L)
  expect_equal(res$diagnostics$solver_mode_final, "qn_polish")
})

test_that("qn_polish saves >= 30% Hessian evaluations on smooth convex", {
  fn <- function(x) sum(0.5 * x^2 + x^4 / 12)
  gr <- function(x) x + x^3 / 3
  hess_fn <- function(x) diag(1 + x^2)

  x0 <- rep(10, 5)

  res_off <- arcopt(x0, fn, gr, hess_fn,
                    control = list(maxit = 50, trace = 0,
                                   qn_polish_enabled = FALSE))
  evals_off <- res_off$evaluations$hess

  res_on <- arcopt(x0, fn, gr, hess_fn,
                   control = list(maxit = 50, trace = 0,
                                  qn_polish_enabled = TRUE))
  evals_on <- res_on$evaluations$hess

  expect_true(res_off$converged)
  expect_true(res_on$converged)
  expect_gte(res_on$diagnostics$qn_polish_switches, 1L)
  # Saves >= 30% Hessian evals on this problem with default thresholds
  expect_true(evals_on <= 0.7 * evals_off,
              info = sprintf("evals_off=%d evals_on=%d", evals_off, evals_on))
})

test_that("qn_polish saves >= 50% Hessian evaluations with loose thresholds", {
  fn <- function(x) sum(0.5 * x^2 + x^4 / 12)
  gr <- function(x) x + x^3 / 3
  hess_fn <- function(x) diag(1 + x^2)

  x0 <- rep(10, 5)
  res_off <- arcopt(x0, fn, gr, hess_fn,
                    control = list(maxit = 50, trace = 0,
                                   qn_polish_enabled = FALSE))

  res_on <- arcopt(x0, fn, gr, hess_fn,
                   control = list(maxit = 50, trace = 0,
                                  qn_polish_enabled = TRUE,
                                  qn_polish_window = 3,
                                  qn_polish_rho = 0.3,
                                  qn_polish_g_decay = 0.95))

  expect_true(res_off$converged)
  expect_true(res_on$converged)
  expect_gte(res_on$diagnostics$qn_polish_switches, 1L)
  expect_true(res_on$evaluations$hess <= 0.5 * res_off$evaluations$hess,
              info = sprintf("evals_off=%d evals_on=%d",
                             res_off$evaluations$hess,
                             res_on$evaluations$hess))
})

test_that("qn_polish does not fire on pure quadratic (Newton done in 1 iter)", {
  # Pure quadratic: Newton converges in one step -- detector window never
  # fills. This is a test that the detector doesn't spuriously fire.
  A <- diag(c(1, 2, 4, 8, 16))
  b <- c(1, 1, 1, 1, 1)
  fn <- function(x) 0.5 * as.numeric(t(x) %*% A %*% x) - sum(b * x)
  gr <- function(x) as.vector(A %*% x) - b
  hess_fn <- function(x) A

  res <- arcopt(rep(0, 5), fn, gr, hess_fn,
                control = list(maxit = 50, trace = 0,
                               qn_polish_enabled = TRUE))
  expect_true(res$converged)
  expect_equal(res$diagnostics$qn_polish_switches, 0L)
  expect_equal(res$diagnostics$solver_mode_final, "cubic")
  expect_equal(res$par, solve(A, b), tolerance = 1e-6)
})

test_that("qn_polish and tr_fallback do not collide on well-behaved problems", {
  # Neither should fire on Rosenbrock with default settings
  rosen <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  rosen_gr <- function(x) {
    c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
      200 * (x[2] - x[1]^2))
  }
  rosen_hess <- function(x) {
    matrix(c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
             -400 * x[1], 200), 2, 2)
  }

  res <- arcopt(c(-1.2, 1), rosen, rosen_gr, rosen_hess,
                control = list(maxit = 200, trace = 0,
                               tr_fallback_enabled = TRUE,
                               qn_polish_enabled = TRUE))
  expect_true(res$converged)
  expect_equal(res$par, c(1, 1), tolerance = 1e-4)
  expect_equal(res$diagnostics$ridge_switches, 0L)
})

test_that("qn_polish branch produces finite outputs on near-minimum start", {
  # Starting very close to the minimum: polish detector should fire fast
  fn <- function(x) 0.5 * sum(x^2) + 0.01 * sum(x^4)
  gr <- function(x) x + 0.04 * x^3
  hess_fn <- function(x) diag(1 + 0.12 * x^2)

  x0 <- rep(0.01, 5) # very close to the minimum
  res <- arcopt(x0, fn, gr, hess_fn,
                control = list(maxit = 50, trace = 0,
                               qn_polish_enabled = TRUE))
  expect_true(all(is.finite(res$par)))
  expect_true(is.finite(res$value))
})

# ---- Phase 4: bidirectional logic (revert + cooldown) ----

test_that("qn_polish reverts to cubic when line search fails repeatedly", {
  # Force line search to fail by giving it a single iteration budget plus
  # a very tight curvature condition c2 = 0.01. After
  # qn_polish_max_fail = 2 consecutive failures, polish should revert to
  # cubic and the revert counter should reflect this.
  fn <- function(x) sum(0.5 * x^2 + x^4 / 12)
  gr <- function(x) x + x^3 / 3
  hess_fn <- function(x) diag(1 + x^2)

  x0 <- rep(10, 5)
  res <- arcopt(x0, fn, gr, hess_fn,
                control = list(maxit = 50, trace = 0,
                               qn_polish_enabled = TRUE,
                               qn_polish_window = 3,
                               qn_polish_rho = 0.3,
                               qn_polish_g_decay = 0.95,
                               qn_polish_c2 = 0.01,
                               qn_polish_max_ls_iter = 1,
                               qn_polish_max_fail = 2))

  # Polish should fire, then revert, possibly multiple times.
  # Must still converge.
  expect_true(res$converged)
  expect_gte(res$diagnostics$qn_polish_switches, 1L,
             label = "polish fired at least once")
  expect_gte(res$diagnostics$qn_polish_reverts, 1L,
             label = "revert happened at least once under tight c2")
})

test_that("qn_polish cooldown prevents immediate re-entry after revert", {
  # Same setup as above, but check that reverts don't happen every
  # single iteration (cooldown spreads them out). If the cooldown weren't
  # respected, we'd see reverts == switches - 0 or 1 with each iteration
  # reverting. With cooldown, reverts should be bounded even on a
  # LS-hostile problem.
  fn <- function(x) sum(0.5 * x^2 + x^4 / 12)
  gr <- function(x) x + x^3 / 3
  hess_fn <- function(x) diag(1 + x^2)

  x0 <- rep(10, 5)
  res <- arcopt(x0, fn, gr, hess_fn,
                control = list(maxit = 50, trace = 0,
                               qn_polish_enabled = TRUE,
                               qn_polish_window = 3,
                               qn_polish_rho = 0.3,
                               qn_polish_g_decay = 0.95,
                               qn_polish_c2 = 0.01,
                               qn_polish_max_ls_iter = 1,
                               qn_polish_max_fail = 2,
                               qn_polish_reenter_delay = 5))

  expect_true(res$converged)
  # With reenter_delay = 5, between consecutive polish->cubic->polish
  # cycles there must be at least 5 cubic iterations. Over maxit = 50,
  # this bounds the number of switches.
  expect_true(res$diagnostics$qn_polish_switches <= ceiling(50 / 5),
              info = sprintf("switches = %d", res$diagnostics$qn_polish_switches))
})

test_that("qn_polish revert leaves cubic state usable", {
  # After a revert, cubic mode must resume correctly and converge
  fn <- function(x) sum(0.5 * x^2 + x^4 / 12)
  gr <- function(x) x + x^3 / 3
  hess_fn <- function(x) diag(1 + x^2)

  x0 <- rep(10, 5)
  res <- arcopt(x0, fn, gr, hess_fn,
                control = list(maxit = 100, trace = 0,
                               qn_polish_enabled = TRUE,
                               qn_polish_window = 3,
                               qn_polish_rho = 0.3,
                               qn_polish_g_decay = 0.95,
                               qn_polish_c2 = 0.01,
                               qn_polish_max_ls_iter = 1,
                               qn_polish_max_fail = 2))
  expect_true(res$converged)
  expect_true(max(abs(res$par)) < 1e-3,
              info = "post-revert cubic mode must still find minimum")
  expect_true(all(is.finite(res$gradient)))
})
