# Trace Rasch JML: see whether polish can fire at default thresholds.

source("benchmarks/rasch_jml/setup.R")
devtools::load_all(".", quiet = TRUE)

data <- simulate_rasch(n_persons = 400L, n_items = 100L, seed = 42L)
fns <- make_rasch_fns(data)
x0 <- rasch_cold_start(400L, 100L)

cat("---- arcopt run with trace, gtol = 1e-10 ----\n")

res <- arcopt(x0, fn = fns$fn, gr = fns$gr, hess = fns$hess,
              control = list(maxit = 200L,
                             gtol_abs = 1e-10, ftol_abs = 1e-13,
                             xtol_abs = 1e-13,
                             qn_polish_enabled = TRUE,
                             trace = 2L))

tr <- res$trace
trace_df <- data.frame(
  iter = seq_along(tr$f) - 1L,
  f = tr$f,
  g_inf = tr$g_norm,
  sigma = tr$sigma,
  rho = tr$rho,
  step_type = tr$step_type,
  lambda_min = if (!is.null(tr$lambda)) tr$lambda else NA_real_
)
print(trace_df, row.names = FALSE, digits = 4)

cat(sprintf("\nFinal: iter=%d  hess=%d  fn=%d  gr=%d  mode=%s\n",
            res$iterations, res$evaluations$hess,
            res$evaluations$fn, res$evaluations$gr,
            res$diagnostics$solver_mode_final))
cat(sprintf("Polish: switches=%d reverts=%d hess@switch=%s\n",
            res$diagnostics$qn_polish_switches,
            res$diagnostics$qn_polish_reverts,
            format(res$diagnostics$hess_evals_at_polish_switch)))

# also count consecutive iterations meeting EACH polish signal at default
# thresholds (window=5, rho>=0.9, lambda_min>=1e-3, g_decay<=0.5,
# g_inf_floor=1e-8 at window start). This tells us which signal blocks
# polish from firing.
cat("\n---- per-iteration polish signals (default thresholds) ----\n")
n <- nrow(trace_df)
sig_rho <- trace_df$rho >= 0.9 & !is.na(trace_df$rho)
sig_lam <- trace_df$lambda_min >= 1e-3 & !is.na(trace_df$lambda_min)
sig_decay <- c(NA, trace_df$g_inf[-1] / trace_df$g_inf[-n] <= 0.5)
sig_floor <- trace_df$g_inf > 1e-8
sig_all <- sig_rho & sig_lam & sig_decay & sig_floor
sig_df <- data.frame(
  iter = trace_df$iter,
  rho_ok = sig_rho, lam_ok = sig_lam,
  decay_ok = sig_decay, floor_ok = sig_floor, all_ok = sig_all
)
print(sig_df, row.names = FALSE)

# longest streak of all_ok
streak <- 0L; max_streak <- 0L
for (k in seq_len(n)) {
  if (isTRUE(sig_all[k])) {
    streak <- streak + 1L
    if (streak > max_streak) max_streak <- streak
  } else {
    streak <- 0L
  }
}
cat(sprintf("\nLongest streak of all-signals-pass: %d (need 5 for polish)\n",
            max_streak))
