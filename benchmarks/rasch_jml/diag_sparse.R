# Look inside one sparse-Rasch run and confirm robustness across seeds
# and items-per-person choices.

source("benchmarks/rasch_jml/setup.R")
source("benchmarks/rasch_jml/benchmark.R", echo = FALSE)
devtools::load_all(".", quiet = TRUE)

cat("\n---- single sparse run trace (10 items/person, seed=101) ----\n")
set.seed(101L)
data <- simulate_rasch(n_persons = 400L, n_items = 100L, seed = 101L)
mask <- matrix(0L, 400L, 100L)
for (i in seq_len(400L)) mask[i, sample.int(100L, 10L)] <- 1L
data$mask <- mask
fns <- make_rasch_fns_sparse(data)
x0 <- rasch_cold_start(400L, 100L)

res <- arcopt(x0, fn = fns$fn, gr = fns$gr, hess = fns$hess,
              control = list(maxit = 500L, gtol_abs = 1e-8,
                             qn_polish_enabled = TRUE,
                             trace = 2L))
tr <- res$trace
trace_df <- data.frame(
  iter = seq_along(tr$f) - 1L,
  f = tr$f, g_inf = tr$g_norm, sigma = tr$sigma,
  rho = tr$rho, step = tr$step_type,
  lam_min = if (!is.null(tr$lambda)) tr$lambda else NA_real_
)
print(trace_df, row.names = FALSE, digits = 4)
cat(sprintf("\nFinal: iter=%d hess=%d mode=%s switches=%d hess@switch=%s\n",
            res$iterations, res$evaluations$hess,
            res$diagnostics$solver_mode_final,
            res$diagnostics$qn_polish_switches,
            format(res$diagnostics$hess_evals_at_polish_switch)))

# Sweep across items-per-person to see how saving varies
cat("\n---- sweep: items/person in {5, 10, 20, 50, 100} ----\n")
sweep_res <- list()
for (ipp in c(5L, 10L, 20L, 50L, 100L)) {
  for (rep_i in 1L:5L) {
    data_s <- simulate_rasch(n_persons = 400L, n_items = 100L,
                             seed = 200L + rep_i)
    if (ipp < 100L) {
      mask <- matrix(0L, 400L, 100L)
      for (i in seq_len(400L)) mask[i, sample.int(100L, ipp)] <- 1L
      data_s$mask <- mask
      fns_s <- make_rasch_fns_sparse(data_s)
    } else {
      fns_s <- make_rasch_fns(data_s)
    }
    x0_s <- rasch_cold_start(400L, 100L)
    for (polish in c(FALSE, TRUE)) {
      ctrl <- list(maxit = 500L, gtol_abs = 1e-8)
      if (polish) ctrl$qn_polish_enabled <- TRUE
      t0 <- proc.time()[["elapsed"]]
      r <- arcopt(x0_s, fn = fns_s$fn, gr = fns_s$gr, hess = fns_s$hess,
                  control = ctrl)
      t1 <- proc.time()[["elapsed"]]
      sweep_res[[length(sweep_res) + 1L]] <- data.frame(
        ipp = ipp, rep = rep_i, polish = polish,
        iter = r$iterations, hess = r$evaluations$hess,
        time_s = t1 - t0,
        switches = r$diagnostics$qn_polish_switches
      )
    }
  }
}
sw_df <- do.call(rbind, sweep_res)

agg <- aggregate(cbind(iter, hess, time_s) ~ ipp + polish, data = sw_df,
                 FUN = stats::median)
agg <- agg[order(agg$ipp, agg$polish), ]
print(agg, row.names = FALSE)

# Saving by ipp
cat("\nHess saving by items-per-person:\n")
for (ipp in unique(sw_df$ipp)) {
  d_med <- stats::median(sw_df$hess[sw_df$ipp == ipp & !sw_df$polish])
  p_med <- stats::median(sw_df$hess[sw_df$ipp == ipp & sw_df$polish])
  t_d <- stats::median(sw_df$time_s[sw_df$ipp == ipp & !sw_df$polish])
  t_p <- stats::median(sw_df$time_s[sw_df$ipp == ipp & sw_df$polish])
  cat(sprintf("  ipp=%3d   default=%4.1f hess  polish=%4.1f hess  saving=%5.1f%%   time saving=%5.1f%%\n",
              ipp, d_med, p_med, 100 * (d_med - p_med) / d_med,
              100 * (t_d - t_p) / t_d))
}
