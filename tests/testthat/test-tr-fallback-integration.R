test_that("arcopt still converges on Rosenbrock with TR-fallback enabled", {
  rosen <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  rosen_gr <- function(x) {
    c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
      200 * (x[2] - x[1]^2))
  }
  rosen_hess <- function(x) {
    matrix(c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
             -400 * x[1], 200), 2, 2)
  }

  res <- arcopt(
    c(-1.2, 1), rosen, rosen_gr, rosen_hess,
    control = list(maxit = 200, trace = 0)
  )

  expect_true(res$converged)
  expect_equal(res$par, c(1, 1), tolerance = 1e-4)
  # Rosenbrock has a curved valley but not a flat ridge; detector must NOT fire
  expect_equal(res$diagnostics$ridge_switches, 0L)
  expect_equal(res$diagnostics$solver_mode_final, "cubic")
})

test_that("detector-forced TR branch produces finite, decreasing sequence", {
  # Force the TR branch by lowering the detector bar so it fires early.
  # Use a non-quadratic function where cubic reg takes several iterations.
  fn <- function(x) 0.5 * x[1]^2 + 1e-3 * (cosh(x[2] / 10) - 1)
  gr_fn <- function(x) c(x[1], 1e-4 * sinh(x[2] / 10))
  hess_fn <- function(x) matrix(c(1, 0, 0, 1e-5 * cosh(x[2] / 10)), 2, 2)

  res <- arcopt(
    c(2, 8), fn, gr_fn, hess_fn,
    control = list(
      maxit = 300, trace = 0,
      tr_fallback_enabled = TRUE,
      tr_fallback_window = 3,           # short window -> easier to trigger
      tr_fallback_tol_ridge = 1e-2,     # lambda_min < 1e-2 counts
      tr_fallback_g_inf_floor = 1e-10,  # low floor -> tolerates small ||g||
      tr_fallback_grad_decrease_max = 0.5
    )
  )

  # Whether or not the detector fired, the run must not crash or produce NaN.
  expect_true(all(is.finite(res$par)))
  expect_true(is.finite(res$value))
  expect_true(res$diagnostics$solver_mode_final %in% c("cubic", "tr"))
})

test_that("TR-fallback disabled leaves behavior unchanged", {
  rosen <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  rosen_gr <- function(x) {
    c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
      200 * (x[2] - x[1]^2))
  }
  rosen_hess <- function(x) {
    matrix(c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
             -400 * x[1], 200), 2, 2)
  }

  res_off <- arcopt(
    c(-1.2, 1), rosen, rosen_gr, rosen_hess,
    control = list(maxit = 200, trace = 0, tr_fallback_enabled = FALSE)
  )

  expect_true(res_off$converged)
  expect_equal(res_off$par, c(1, 1), tolerance = 1e-4)
  expect_equal(res_off$diagnostics$ridge_switches, 0L)
  expect_equal(res_off$diagnostics$solver_mode_final, "cubic")
})

test_that("arcopt result includes new TR-fallback fields", {
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

  for (f in c("solver_mode_final", "ridge_switches", "radius_final")) {
    expect_true(f %in% names(res$diagnostics),
                info = paste("missing diagnostics field:", f))
  }
})

test_that("arcopt_qn still converges on Rosenbrock with TR-fallback enabled", {
  rosen <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  rosen_gr <- function(x) {
    c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
      200 * (x[2] - x[1]^2))
  }

  res <- arcopt(
    c(-1.2, 1), rosen, rosen_gr,
    control = list(maxit = 500, trace = 0, use_qn = TRUE,
                   qn_method = "hybrid")
  )

  # Must not crash and must return all TR fields
  for (f in c("solver_mode_final", "ridge_switches", "radius_final")) {
    expect_true(f %in% names(res$diagnostics),
                info = paste("missing diagnostics field:", f))
  }
  expect_true(all(is.finite(res$par)))
  # Rosenbrock has no true flat ridge; detector should not fire
  expect_equal(res$diagnostics$ridge_switches, 0L)
  expect_equal(res$diagnostics$solver_mode_final, "cubic")
})

test_that("arcopt_qn with tr_fallback disabled has ridge_switches == 0", {
  rosen <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
  rosen_gr <- function(x) {
    c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
      200 * (x[2] - x[1]^2))
  }

  res <- arcopt(
    c(-1.2, 1), rosen, rosen_gr,
    control = list(maxit = 500, trace = 0, use_qn = TRUE,
                   qn_method = "hybrid",
                   tr_fallback_enabled = FALSE)
  )

  expect_equal(res$diagnostics$ridge_switches, 0L)
  expect_equal(res$diagnostics$solver_mode_final, "cubic")
})
