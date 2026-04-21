# Deep linear network benchmark (n = 504)
# =======================================
#
# Compares arcopt QN variants against optim BFGS / L-BFGS-B / nlminb on
# a 3-layer linear network with 504 parameters. Saddles are the sole
# non-global critical points (Kawaguchi 2016), so any method's failure
# to reach f* = 0 under a tight gradient tolerance is a saddle-escape
# problem, not a spurious-local-min problem.
#
# Usage (from repo root):
#   Rscript benchmarks/deep_linear/benchmark.R [n_seeds] [init_sd]

suppressPackageStartupMessages({
  library(arcopt)
})

source(file.path("benchmarks", "deep_linear", "setup.R"))

args <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(args) >= 1L) as.integer(args[[1]]) else 5L
init_sd <- if (length(args) >= 2L) as.numeric(args[[2]]) else 0.3

maxit <- 2000L
gtol <- 1e-5

prob <- make_deep_linear_problem(
  d_in = 5L, h1 = 18L, h2 = 18L, d_out = 5L,
  n_samples = 20L, w_star_seed = 1L, data_seed = 2L
)

cat(sprintf("Deep linear network: n = %d\n", prob$dims$n))
cat(sprintf("  d_in=%d  h1=%d  h2=%d  d_out=%d  N=%d\n",
            prob$dims$d_in, prob$dims$h1, prob$dims$h2,
            prob$dims$d_out, prob$dims$n_samples))
cat(sprintf("  ||Y||_F^2/2 (f at origin) = %.2f\n",
            0.5 * sum(prob$y_data^2)))
cat(sprintf("  seeds=%d  init_sd=%.3f  maxit=%d  gtol=%.0e\n\n",
            n_seeds, init_sd, maxit, gtol))

methods <- list(
  list(name = "arcopt-hybrid", runner = "arcopt_qn",
       qn_method = "hybrid"),
  list(name = "arcopt-sr1",    runner = "arcopt_qn",
       qn_method = "sr1"),
  list(name = "arcopt-bfgs",   runner = "arcopt_qn",
       qn_method = "bfgs"),
  list(name = "optim-BFGS",    runner = "optim_bfgs"),
  list(name = "optim-LBFGSB",  runner = "optim_lbfgsb"),
  list(name = "nlminb",        runner = "nlminb")
)

run_method <- function(method, x0) {
  t0 <- proc.time()[["elapsed"]]

  if (method$runner == "optim_bfgs") {
    fit <- tryCatch(stats::optim(
      par = x0, fn = prob$fn, gr = prob$gr, method = "BFGS",
      control = list(maxit = maxit, reltol = 1e-12)
    ), error = function(e) NULL)
    if (is.null(fit)) {
      return(list(f = NA_real_, gnorm_inf = NA_real_, gnorm_2 = NA_real_,
                  iter = NA_integer_,
                  time = proc.time()[["elapsed"]] - t0, conv = FALSE))
    }
    g_final <- prob$gr(fit$par)
    return(list(
      f = fit$value,
      gnorm_inf = max(abs(g_final)),
      gnorm_2 = sqrt(sum(g_final^2)),
      iter = fit$counts[["gradient"]],
      time = proc.time()[["elapsed"]] - t0,
      conv = fit$convergence == 0L && max(abs(g_final)) < gtol
    ))
  }

  if (method$runner == "optim_lbfgsb") {
    fit <- tryCatch(stats::optim(
      par = x0, fn = prob$fn, gr = prob$gr, method = "L-BFGS-B",
      control = list(maxit = maxit, factr = 1e1, pgtol = gtol)
    ), error = function(e) NULL)
    if (is.null(fit)) {
      return(list(f = NA_real_, gnorm_inf = NA_real_, gnorm_2 = NA_real_,
                  iter = NA_integer_,
                  time = proc.time()[["elapsed"]] - t0, conv = FALSE))
    }
    g_final <- prob$gr(fit$par)
    return(list(
      f = fit$value,
      gnorm_inf = max(abs(g_final)),
      gnorm_2 = sqrt(sum(g_final^2)),
      iter = fit$counts[["gradient"]],
      time = proc.time()[["elapsed"]] - t0,
      conv = fit$convergence == 0L && max(abs(g_final)) < gtol
    ))
  }

  if (method$runner == "nlminb") {
    fit <- tryCatch(stats::nlminb(
      start = x0, objective = prob$fn, gradient = prob$gr,
      control = list(iter.max = maxit, eval.max = 4L * maxit,
                     abs.tol = 1e-14, rel.tol = 1e-12)
    ), error = function(e) NULL)
    if (is.null(fit)) {
      return(list(f = NA_real_, gnorm_inf = NA_real_, gnorm_2 = NA_real_,
                  iter = NA_integer_,
                  time = proc.time()[["elapsed"]] - t0, conv = FALSE))
    }
    g_final <- prob$gr(fit$par)
    return(list(
      f = fit$objective,
      gnorm_inf = max(abs(g_final)),
      gnorm_2 = sqrt(sum(g_final^2)),
      iter = fit$iterations,
      time = proc.time()[["elapsed"]] - t0,
      conv = fit$convergence == 0L && max(abs(g_final)) < gtol
    ))
  }

  # arcopt_qn path
  ctrl <- list(
    maxit = maxit, gtol_abs = gtol, trace = 0L,
    use_qn = TRUE, qn_method = method$qn_method
  )
  fit <- tryCatch(
    arcopt::arcopt(x0 = x0, fn = prob$fn, gr = prob$gr, control = ctrl),
    error = function(e) NULL
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  if (is.null(fit)) {
    return(list(f = NA_real_, gnorm_inf = NA_real_, gnorm_2 = NA_real_,
                iter = NA_integer_, time = elapsed, conv = FALSE))
  }
  g_final <- fit$gradient
  list(
    f = fit$value,
    gnorm_inf = max(abs(g_final)),
    gnorm_2 = sqrt(sum(g_final^2)),
    iter = fit$iterations,
    time = elapsed,
    conv = isTRUE(fit$converged) && max(abs(g_final)) < gtol
  )
}

results <- list()
for (seed in seq_len(n_seeds)) {
  x0 <- make_deep_linear_start(prob$dims, seed = 100L + seed, sd = init_sd)
  cat(sprintf("Seed %d (f0 = %.2e):\n", seed, prob$fn(x0)))
  for (method in methods) {
    res <- run_method(method, x0)
    results[[length(results) + 1L]] <- data.frame(
      seed = seed, method = method$name,
      f = res$f, gnorm_inf = res$gnorm_inf, gnorm_2 = res$gnorm_2,
      iter = res$iter, time_s = round(res$time, 2),
      converged = res$conv, stringsAsFactors = FALSE
    )
    cat(sprintf("  %-14s f=%.2e  |g|_inf=%.2e  iter=%5s  time=%6.2fs  conv=%s\n",
                method$name, res$f, res$gnorm_inf,
                ifelse(is.na(res$iter), "NA", format(res$iter)),
                res$time, res$conv))
  }
  cat("\n")
}

all_results <- do.call(rbind, results)

cat("\n=============================================================\n")
cat(sprintf("Summary across %d seeds (converged = |g|_inf < %.0e)\n",
            n_seeds, gtol))
cat("=============================================================\n")
summary_tbl <- do.call(rbind, lapply(split(all_results, all_results$method),
  function(g) {
    data.frame(
      method = g$method[[1]],
      n_conv = sum(g$converged, na.rm = TRUE),
      n_total = nrow(g),
      pct_conv = round(100 * mean(g$converged, na.rm = TRUE), 0),
      median_f = stats::median(g$f, na.rm = TRUE),
      median_iter = stats::median(g$iter, na.rm = TRUE),
      median_time_s = round(stats::median(g$time_s, na.rm = TRUE), 2),
      stringsAsFactors = FALSE
    )
  }))
summary_tbl <- summary_tbl[match(
  vapply(methods, function(m) m$name, character(1)),
  summary_tbl$method
), ]
rownames(summary_tbl) <- NULL
print(summary_tbl, row.names = FALSE)

out_dir <- file.path("benchmarks", "deep_linear", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(
  all_results,
  file = file.path(out_dir,
                   sprintf("results_n%d_sd%s_seeds%d.csv",
                           prob$dims$n, format(init_sd), n_seeds)),
  row.names = FALSE
)
cat(sprintf("\nResults: %s\n",
            file.path(out_dir,
                      sprintf("results_n%d_sd%s_seeds%d.csv",
                              prob$dims$n, format(init_sd), n_seeds))))
