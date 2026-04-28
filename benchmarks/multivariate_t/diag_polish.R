# Diagnostic: figure out why polish never fires at MoM start.
# Try (a) a cold start (mu=0, Sigma=I, nu=10) and (b) heavier-tail truth.

source("benchmarks/multivariate_t/setup.R")
devtools::load_all(".", quiet = TRUE)

cold_start <- function(p) {
  mvt_pack(rep(0, p), diag(p), nu = 30)
}

run_one <- function(data, x0, polish, gtol = 1e-8) {
  fns <- make_mvt_fns(data)
  ctrl <- list(maxit = 500L, gtol_abs = gtol, trace = 1L)
  if (polish) ctrl$qn_polish_enabled <- TRUE
  res <- arcopt(x0, fn = fns$fn, gr = fns$gr, hess = fns$hess,
                control = ctrl)
  list(res = res, x0 = x0)
}

scenarios <- list(
  list(p = 5L,  nu_true = 5,  start = "MoM",  desc = "p=5  / nu=5  / MoM"),
  list(p = 5L,  nu_true = 5,  start = "cold", desc = "p=5  / nu=5  / cold"),
  list(p = 5L,  nu_true = 3,  start = "MoM",  desc = "p=5  / nu=3  / MoM"),
  list(p = 5L,  nu_true = 3,  start = "cold", desc = "p=5  / nu=3  / cold"),
  list(p = 10L, nu_true = 3,  start = "cold", desc = "p=10 / nu=3  / cold")
)

for (sc in scenarios) {
  cat(sprintf("\n==== %s ====\n", sc$desc))
  data <- simulate_mvt(p = sc$p, n = 200L, nu = sc$nu_true, seed = 101L)
  x0 <- if (sc$start == "MoM") mvt_start(data) else cold_start(sc$p)

  for (polish in c(FALSE, TRUE)) {
    out <- run_one(data, x0, polish)
    res <- out$res
    cat(sprintf("  polish=%-5s  iter=%3d  hess=%3d  mode=%s",
                polish, res$iterations, res$evaluations$hess,
                res$diagnostics$solver_mode_final))
    if (polish) {
      cat(sprintf("  switches=%d  reverts=%d  hess@switch=%s",
                  res$diagnostics$qn_polish_switches,
                  res$diagnostics$qn_polish_reverts,
                  format(res$diagnostics$hess_evals_at_polish_switch)))
    }
    cat("\n")
  }
}
