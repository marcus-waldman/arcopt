# Push tolerances down hard, and force a non-default qn_polish_window
# to see when (if ever) polish can fire on the multivariate-t MLE.

source("benchmarks/multivariate_t/setup.R")
devtools::load_all(".", quiet = TRUE)

cold_start <- function(p) mvt_pack(rep(0, p), diag(p), nu = 30)

run_one <- function(data, x0, ctrl) {
  fns <- make_mvt_fns(data)
  arcopt(x0, fn = fns$fn, gr = fns$gr, hess = fns$hess, control = ctrl)
}

scenarios <- list(
  list(p = 5L,  gtol = 1e-12, ftol = 1e-15, xtol = 1e-15,
       window = 5L, desc = "p=5 default-window=5, very tight tol"),
  list(p = 5L,  gtol = 1e-12, ftol = 1e-15, xtol = 1e-15,
       window = 3L, desc = "p=5 window=3, very tight tol"),
  list(p = 10L, gtol = 1e-10, ftol = 1e-13, xtol = 1e-13,
       window = 5L, desc = "p=10 default-window=5, tight tol"),
  list(p = 10L, gtol = 1e-10, ftol = 1e-13, xtol = 1e-13,
       window = 3L, desc = "p=10 window=3, tight tol"),
  list(p = 20L, gtol = 1e-10, ftol = 1e-13, xtol = 1e-13,
       window = 5L, desc = "p=20 default-window=5, tight tol")
)

for (sc in scenarios) {
  cat(sprintf("\n==== %s ====\n", sc$desc))
  data <- simulate_mvt(p = sc$p, n = 500L, nu = 5, seed = 101L)
  x0 <- cold_start(sc$p)

  for (polish in c(FALSE, TRUE)) {
    ctrl <- list(maxit = 500L,
                 gtol_abs = sc$gtol, ftol_abs = sc$ftol, xtol_abs = sc$xtol)
    if (polish) {
      ctrl$qn_polish_enabled <- TRUE
      ctrl$qn_polish_window <- sc$window
    }
    res <- run_one(data, x0, ctrl)
    cat(sprintf("  polish=%-5s  iter=%3d  hess=%3d  mode=%s",
                polish, res$iterations, res$evaluations$hess,
                res$diagnostics$solver_mode_final))
    if (polish) {
      cat(sprintf("  switches=%d  reverts=%d  hess@switch=%s",
                  res$diagnostics$qn_polish_switches,
                  res$diagnostics$qn_polish_reverts,
                  format(res$diagnostics$hess_evals_at_polish_switch)))
    }
    cat(sprintf("  |g|=%.2e\n", max(abs(res$gradient))))
  }
}
