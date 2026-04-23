# Diagnostic trace of flat-ridge detector signals on DPM Heckman.
# Runs arcopt with hybrid disabled and dumps the four signals per iteration
# so we can see whether the detector SHOULD have fired but didn't, or
# whether the landscape actually doesn't match the ridge pattern.

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})
source(file.path("benchmarks", "dpm_heckman", "setup.R"))

data <- simulate_dpm_heckman(n_obs = 500L, k_true = 2L,
                             p_out = 5L, p_sel = 7L, n_excl = 2L,
                             seed = 101L)
model_so <- compile_dpm_heckman_model()
fns <- make_dpm_heckman_fns(data, model_so, k_trunc = 20L)

# Use a custom hess wrapper that also caches per-iteration eigenvalues.
# We pipe them into a global via environment.
trace_env <- new.env()
trace_env$lambda_min <- numeric(0)
trace_env$g_inf <- numeric(0)

hess_wrapped <- function(x) {
  h <- fns$hess_ad(x)
  trace_env$lambda_min <- c(
    trace_env$lambda_min,
    min(eigen(h, symmetric = TRUE, only.values = TRUE)$values)
  )
  h
}

cat("Running arcopt-AD with hybrid DISABLED, trace=2, maxit=100\n")
cat("Per-iteration signals logged to trace_ridge_diag.log\n")

t0 <- proc.time()[["elapsed"]]
fit <- arcopt(
  x0 = fns$x0, fn = fns$fn, gr = fns$gr, hess = hess_wrapped,
  control = list(maxit = 100L, gtol_abs = 1e-5, trace = 2L,
                 tr_fallback_enabled = FALSE)
)
elapsed <- proc.time()[["elapsed"]] - t0

cat(sprintf("\nFinal: f=%.3f  iters=%d  |g|_inf=%.2e  elapsed=%.1fs\n",
            fit$value, fit$iterations, max(abs(fit$gradient)), elapsed))

sig <- fit$trace$sigma
rho <- fit$trace$rho
# trace$g_norm is l2; derive |g|_inf from trace if present else skip
g_norm <- fit$trace$g_norm
lam <- trace_env$lambda_min
# lam length may differ from iter count; align from tail
min_len <- min(length(sig), length(rho), length(lam), length(g_norm))
sig <- tail(sig, min_len)
rho <- tail(rho, min_len)
g_norm <- tail(g_norm, min_len)
lam <- tail(lam, min_len)

log_path <- "benchmarks/dpm_heckman/trace_ridge_diag.log"
try(file.remove(log_path), silent = TRUE)
log_w <- function(x) { cat(x, file = log_path, append = TRUE); cat(x) }

log_w(sprintf("iter   sigma        rho      ||g||_2    lambda_min\n"))
for (i in seq_len(min_len)) {
  log_w(sprintf("%4d  %10.4e  %7.4f  %10.4e  %+10.4e\n",
                i, sig[i], rho[i], g_norm[i], lam[i]))
}

# Summary: how often would each signal fire?
sigma_min_default <- 1e-6
sigma_pinned <- sig <= 10 * sigma_min_default
rho_near1 <- abs(rho - 1) < 0.1
pd_and_ridge <- lam > 0 & lam < 1e-3
any_indef <- lam < 0
log_w(sprintf("\nSignal frequencies over %d iterations:\n", min_len))
log_w(sprintf("  sigma <= 10*sigma_min  : %d / %d  (%.0f%%)\n",
              sum(sigma_pinned), min_len, 100 * mean(sigma_pinned)))
log_w(sprintf("  |rho - 1| < 0.1         : %d / %d  (%.0f%%)\n",
              sum(rho_near1, na.rm = TRUE), min_len,
              100 * mean(rho_near1, na.rm = TRUE)))
log_w(sprintf("  0 < lambda_min < 1e-3   : %d / %d  (%.0f%%)\n",
              sum(pd_and_ridge), min_len, 100 * mean(pd_and_ridge)))
log_w(sprintf("  lambda_min < 0 (indef.) : %d / %d  (%.0f%%)\n",
              sum(any_indef), min_len, 100 * mean(any_indef)))
