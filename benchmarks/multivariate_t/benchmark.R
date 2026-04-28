# §4.4 multivariate-t MLE polish go/no-go
# =======================================
#
# Question: at default `qn_polish_*` thresholds, does enabling polish
# save >= 30% of Hessian evaluations on a multivariate-t MLE in the
# strongly-convex post-warm-up basin?
#
# Headline configuration: p = 5, n = 200 observations, true nu = 5, MoM
# warm start. n_par = 21.
#
# Stress configuration: p = 10, n = 500 observations, true nu = 5,
# MoM warm start. n_par = 66.
#
# Records, per replication and configuration:
#   * iterations to convergence (gtol_abs = 1e-8 default)
#   * fn / gr / hess evaluation counts
#   * solver mode at termination
#   * elapsed time
#   * final |grad|_inf and final fval (cross-check both modes agree)

source("benchmarks/multivariate_t/setup.R")
devtools::load_all(".", quiet = TRUE)

run_one <- function(data, polish_enabled, gtol = 1e-8) {
  fns <- make_mvt_fns(data)
  x0 <- mvt_start(data)

  ctrl <- list(maxit = 500L, gtol_abs = gtol)
  if (polish_enabled) ctrl$qn_polish_enabled <- TRUE

  t0 <- proc.time()[["elapsed"]]
  res <- arcopt(x0, fn = fns$fn, gr = fns$gr, hess = fns$hess,
                control = ctrl)
  t1 <- proc.time()[["elapsed"]]

  data.frame(
    polish = polish_enabled,
    iter = res$iterations,
    converged = res$converged,
    fn_evals = res$evaluations$fn,
    gr_evals = res$evaluations$gr,
    hess_evals = res$evaluations$hess,
    final_mode = res$diagnostics$solver_mode_final,
    polish_switches = res$diagnostics$qn_polish_switches,
    polish_reverts = res$diagnostics$qn_polish_reverts,
    g_inf = max(abs(res$gradient)),
    fval = res$value,
    time_s = t1 - t0,
    stringsAsFactors = FALSE
  )
}

run_config <- function(p_dim, n_obs, nu_true, n_reps,
                       gtol = 1e-8, label = NULL) {
  if (is.null(label)) label <- sprintf("p=%d, n=%d", p_dim, n_obs)
  cat(sprintf("\n==== %s (n_par = %d), %d reps ====\n",
              label, mvt_npar(p_dim), n_reps))

  results <- list()
  for (rep_i in seq_len(n_reps)) {
    data <- simulate_mvt(p = p_dim, n = n_obs, nu = nu_true,
                         seed = 100L + rep_i)
    for (polish in c(FALSE, TRUE)) {
      r <- run_one(data, polish, gtol = gtol)
      r$rep <- rep_i
      results[[length(results) + 1L]] <- r
    }
  }
  do.call(rbind, results)
}

summarize_run <- function(df) {
  agg <- aggregate(
    cbind(iter, hess_evals, gr_evals, fn_evals, time_s) ~ polish,
    data = df, FUN = function(z) c(median = stats::median(z),
                                   mean = mean(z),
                                   sd = stats::sd(z)))

  # Hess-eval saving = (default - polish) / default
  med_default <- stats::median(df$hess_evals[!df$polish])
  med_polish <- stats::median(df$hess_evals[df$polish])
  saving_med <- (med_default - med_polish) / med_default

  mean_default <- mean(df$hess_evals[!df$polish])
  mean_polish <- mean(df$hess_evals[df$polish])
  saving_mean <- (mean_default - mean_polish) / mean_default

  cat("\nPer-mode summary:\n")
  print(agg)
  cat(sprintf("\nHess-eval saving (median): %.1f%%\n", 100 * saving_med))
  cat(sprintf("Hess-eval saving (mean):   %.1f%%\n", 100 * saving_mean))
  cat(sprintf("Polish switches (mean):    %.2f\n",
              mean(df$polish_switches[df$polish])))
  cat(sprintf("Polish reverts  (mean):    %.2f\n",
              mean(df$polish_reverts[df$polish])))
  cat(sprintf("Final mode polish runs:   %s\n",
              paste(table(df$final_mode[df$polish]), collapse = " / ")))

  # Convergence and value cross-check
  fval_default <- df$fval[!df$polish]
  fval_polish  <- df$fval[df$polish]
  cat(sprintf("Max |fval_default - fval_polish|: %.3e\n",
              max(abs(fval_default - fval_polish))))
  cat(sprintf("All converged? default=%s polish=%s\n",
              all(df$converged[!df$polish]),
              all(df$converged[df$polish])))

  invisible(list(saving_median = saving_med, saving_mean = saving_mean))
}

# ---- run ----

cat("Running headline (p = 5)...\n")
res_p5 <- run_config(p_dim = 5L, n_obs = 200L, nu_true = 5,
                     n_reps = 10L, gtol = 1e-8,
                     label = "Headline p=5, n=200, nu_true=5")
sum_p5 <- summarize_run(res_p5)

cat("\n\nRunning stress (p = 10)...\n")
res_p10 <- run_config(p_dim = 10L, n_obs = 500L, nu_true = 5,
                      n_reps = 5L, gtol = 1e-8,
                      label = "Stress p=10, n=500, nu_true=5")
sum_p10 <- summarize_run(res_p10)

# Persist results so the manuscript can reproduce them later.
dir.create("benchmarks/multivariate_t/results", showWarnings = FALSE,
           recursive = TRUE)
utils::write.csv(res_p5,
                 "benchmarks/multivariate_t/results/results_p5_n200.csv",
                 row.names = FALSE)
utils::write.csv(res_p10,
                 "benchmarks/multivariate_t/results/results_p10_n500.csv",
                 row.names = FALSE)

cat("\n\n====================================================\n")
cat("Go/no-go verdict (>= 30% Hessian-eval saving at p = 5)\n")
cat("====================================================\n")
cat(sprintf("p=5  median saving: %.1f%%   (criterion: >= 30%%) -> %s\n",
            100 * sum_p5$saving_median,
            if (sum_p5$saving_median >= 0.30) "PASS" else "FAIL"))
cat(sprintf("p=10 median saving: %.1f%%   (informational)\n",
            100 * sum_p10$saving_median))
