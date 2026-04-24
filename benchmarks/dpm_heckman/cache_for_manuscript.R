## cache_for_manuscript.R
## ======================
## Runs the trust / arcopt-AD-nohybrid / arcopt-AD-hybrid comparison on the
## seed=101 DPM Heckman MAP problem (n_params=299) and writes a compact
## results object to manuscript/_cache/dpm_heckman_results.rds so that
## arcopt-jss.qmd can render without depending on BridgeStan.
##
## Reproduces the numbers reported as "Example 4: DPM Heckman MAP" in the
## manuscript. Total runtime on a 2026 workstation is roughly 11 minutes,
## dominated by the nohybrid run exhausting its 200-iter budget.
##
## Usage (from repo root):
##   Rscript benchmarks/dpm_heckman/cache_for_manuscript.R

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(arcopt)
  }
  library(bridgestan)
})

source(file.path("benchmarks", "dpm_heckman", "setup.R"))

seed  <- 101L
maxit <- 200L
gtol  <- 1e-5

model_so <- compile_dpm_heckman_model()
data <- simulate_dpm_heckman(n_obs = 500, k_true = 2, p_out = 5,
                             p_sel = 7, n_excl = 2, seed = seed)
fns <- make_dpm_heckman_fns(data, model_so, k_trunc = 20L)
x0 <- fns$x0

run_one <- function(optimizer) {
  t0 <- proc.time()[["elapsed"]]
  fit <- if (identical(optimizer, "trust")) {
    objfun <- function(p) list(value = fns$fn(p),
                               gradient = fns$gr(p),
                               hessian = fns$hess_ad(p))
    r <- trust::trust(objfun, x0, rinit = 1, rmax = 100, iterlim = maxit)
    list(par = r$argument, value = r$value, iterations = r$iterations,
         gradient = fns$gr(r$argument), solver_mode_final = NA_character_,
         ridge_switches = NA_integer_)
  } else if (identical(optimizer, "arcopt-nohybrid")) {
    arcopt::arcopt(x0 = x0, fn = fns$fn, gr = fns$gr, hess = fns$hess_ad,
                   control = list(maxit = maxit, gtol_abs = gtol, trace = 0L,
                                  tr_fallback_enabled = FALSE))
  } else if (identical(optimizer, "arcopt-hybrid")) {
    arcopt::arcopt(x0 = x0, fn = fns$fn, gr = fns$gr, hess = fns$hess_ad,
                   control = list(maxit = maxit, gtol_abs = gtol, trace = 0L,
                                  tr_fallback_enabled = TRUE))
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  eig_min <- min(eigen(fns$hess_ad(fit$par),
                       symmetric = TRUE, only.values = TRUE)$values)
  list(method = optimizer,
       fn = fit$value,
       g_inf = max(abs(fit$gradient)),
       iter = fit$iterations,
       time_s = elapsed,
       eig_min = eig_min,
       mode = fit$solver_mode_final,
       switches = fit$ridge_switches)
}

cat("[cache] running trust...\n")
res_trust    <- run_one("trust")
cat("[cache] running arcopt-nohybrid...\n")
res_nohybrid <- run_one("arcopt-nohybrid")
cat("[cache] running arcopt-hybrid...\n")
res_hybrid   <- run_one("arcopt-hybrid")

results <- list(
  meta = list(seed = seed, maxit = maxit, gtol = gtol,
              n_params = fns$n_params,
              selection_rate = data$selection_rate,
              f_x0 = fns$fn(x0),
              n_obs = data$n_obs, k_trunc = 20L, k_true = data$k_true,
              timestamp = Sys.time()),
  trust    = res_trust,
  nohybrid = res_nohybrid,
  hybrid   = res_hybrid
)

out_dir <- file.path("manuscript", "_cache")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_file <- file.path(out_dir, "dpm_heckman_results.rds")
saveRDS(results, out_file)
cat(sprintf("[cache] wrote %s\n", out_file))
