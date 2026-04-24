# Diagnostic: does trust::trust converge from the iterate where arcopt's
# cubic regularization gets stuck? Tests whether the hybrid's failure on
# DPM Heckman is a (a) property of the midtrajectory starting point or
# (b) a bug in my TR implementation.

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})
source(file.path("benchmarks", "dpm_heckman", "setup.R"))

data <- simulate_dpm_heckman(n_obs = 500L, k_true = 2L,
                             p_out = 5L, p_sel = 7L, n_excl = 2L,
                             seed = 101L)
model_so <- compile_dpm_heckman_model()
fns <- make_dpm_heckman_fns(data, model_so, k_trunc = 20L)

# Step 1: run arcopt cubic-only for ~65 iterations to get the mid-trajectory iterate
cat("=== Step 1: arcopt cubic-only, maxit=65 ===\n")
t0 <- proc.time()[["elapsed"]]
fit_cubic <- arcopt(fns$x0, fns$fn, fns$gr, fns$hess_ad,
                    control = list(maxit = 65L, gtol_abs = 1e-5, trace = 0,
                                   tr_fallback_enabled = FALSE))
cat(sprintf("  f=%.3f |g|_inf=%.2e iters=%d elapsed=%.1fs\n",
            fit_cubic$value, max(abs(fit_cubic$gradient)),
            fit_cubic$iterations, proc.time()[["elapsed"]] - t0))
x_mid <- fit_cubic$par

# Step 2: run trust::trust from that iterate
cat("\n=== Step 2: trust::trust from the midtrajectory iterate ===\n")
objfun <- function(p) list(value = fns$fn(p),
                           gradient = fns$gr(p),
                           hessian = fns$hess_ad(p))
t0 <- proc.time()[["elapsed"]]
fit_trust <- trust::trust(objfun, x_mid, rinit = 1, rmax = 100,
                          iterlim = 150)
g <- fns$gr(fit_trust$argument)
eig_min <- min(eigen(fns$hess_ad(fit_trust$argument), symmetric = TRUE,
                     only.values = TRUE)$values)
cat(sprintf("  f=%.3f |g|_inf=%.2e iter=%d elapsed=%.1fs eig_min=%+.3e converged=%s\n",
            fit_trust$value, max(abs(g)), fit_trust$iterations,
            proc.time()[["elapsed"]] - t0, eig_min, fit_trust$converged))

# Step 3: for reference, trust::trust from x0
cat("\n=== Step 3: trust::trust from x0 (baseline) ===\n")
t0 <- proc.time()[["elapsed"]]
fit_trust_x0 <- trust::trust(objfun, fns$x0, rinit = 1, rmax = 100,
                             iterlim = 150)
g <- fns$gr(fit_trust_x0$argument)
eig_min <- min(eigen(fns$hess_ad(fit_trust_x0$argument), symmetric = TRUE,
                     only.values = TRUE)$values)
cat(sprintf("  f=%.3f |g|_inf=%.2e iter=%d elapsed=%.1fs eig_min=%+.3e converged=%s\n",
            fit_trust_x0$value, max(abs(g)), fit_trust_x0$iterations,
            proc.time()[["elapsed"]] - t0, eig_min, fit_trust_x0$converged))

cat("\n=== Verdict ===\n")
cat("If trust from midtrajectory escapes (eig_min > 0), my hybrid has a bug.\n")
cat("If trust from midtrajectory also gets stuck, the saddle trajectory is\n")
cat("  the problem, not the TR implementation — we need to switch EARLIER.\n")
