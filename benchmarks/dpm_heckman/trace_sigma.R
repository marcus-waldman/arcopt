library(arcopt)
source("benchmarks/dpm_heckman/setup.R")

data <- simulate_dpm_heckman(n_obs = 500L, k_true = 2L,
                             p_out = 5L, p_sel = 7L, n_excl = 2L,
                             seed = 101L)
model_so <- compile_dpm_heckman_model()
fns <- make_dpm_heckman_fns(data, model_so, k_trunc = 20L)

cat("Running arcopt-AD with trace=2, maxit=100 to see sigma evolution\n")
t0 <- proc.time()[["elapsed"]]
fit <- arcopt(x0 = fns$x0, fn = fns$fn, gr = fns$gr, hess = fns$hess_ad,
              control = list(maxit = 100L, gtol_abs = 1e-5, trace = 2L))
cat(sprintf("\nFinal: f=%.2f  iters=%d  time=%.1fs  |g|_inf=%.2e\n",
            fit$value, fit$iterations,
            proc.time()[["elapsed"]] - t0,
            max(abs(fit$gradient))))
cat("\nsigma trajectory (first 30 + last 10 iterations):\n")
if (!is.null(fit$trace$sigma)) {
  sig <- fit$trace$sigma
  rho <- fit$trace$rho
  n_iter <- length(sig)
  cat(sprintf("  iter  sigma            rho\n"))
  for (i in c(seq_len(min(30, n_iter)),
              if (n_iter > 40) seq.int(n_iter - 9, n_iter))) {
    cat(sprintf("  %4d  %12.4e  %8.4f\n", i, sig[i], rho[i]))
  }
  cat(sprintf("\n  min(sigma) = %.4e  max(sigma) = %.4e\n",
              min(sig), max(sig)))
  cat(sprintf("  fraction of unsuccessful steps (rho < 0.1): %.2f\n",
              mean(rho < 0.1, na.rm = TRUE)))
} else {
  cat("  trace$sigma not populated\n")
  cat("  names in fit$trace:", names(fit$trace), "\n")
}
