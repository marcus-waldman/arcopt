# DPM Mixture Heckman benchmark
# =============================
#
# Compare arcopt variants against optim-BFGS, nlminb, and trust on the
# MAP estimation problem for a truncated DPM Heckman model.
#
# Usage (from repo root):
#   Rscript benchmarks/dpm_heckman/benchmark.R [n_seeds] [n_obs] [k_trunc] [k_true]

suppressPackageStartupMessages({
  library(arcopt)
})

source(file.path("benchmarks", "dpm_heckman", "setup.R"))

args <- commandArgs(trailingOnly = TRUE)
n_seeds  <- if (length(args) >= 1L) as.integer(args[[1]]) else 3L
n_obs    <- if (length(args) >= 2L) as.integer(args[[2]]) else 500L
k_trunc  <- if (length(args) >= 3L) as.integer(args[[3]]) else 20L
k_true   <- if (length(args) >= 4L) as.integer(args[[4]]) else 2L

p_out <- 5L
p_sel <- 7L
n_excl <- 2L
maxit <- 500L   # tight cap so one slow method doesn't drag wall time
gtol <- 1e-5

# Live progress log (writes bypass R's stdout buffering)
log_path <- file.path("benchmarks", "dpm_heckman", "progress.log")
try(file.remove(log_path), silent = TRUE)
progress <- function(...) {
  msg <- paste0(sprintf(...), "\n")
  cat(msg)
  cat(msg, file = log_path, append = TRUE)
}

# Compile once outside the seed loop (cached via hash on the .so).
cat("Compiling Stan model...\n")
model_so <- compile_dpm_heckman_model()

methods <- list(
  list(name = "arcopt-AD",     runner = "arcopt_exact", hess = "ad"),
  list(name = "arcopt-FD",     runner = "arcopt_exact", hess = "fd"),
  list(name = "arcopt-hybrid", runner = "arcopt_qn",    qn = "hybrid"),
  list(name = "arcopt-sr1",    runner = "arcopt_qn",    qn = "sr1"),
  list(name = "arcopt-bfgs",   runner = "arcopt_qn",    qn = "bfgs"),
  list(name = "optim-BFGS",    runner = "optim_bfgs"),
  list(name = "nlminb",        runner = "nlminb"),
  list(name = "trust",         runner = "trust")
)

run_method <- function(method, fns, x0) {
  t0 <- proc.time()[["elapsed"]]
  if (method$runner == "optim_bfgs") {
    fit <- tryCatch(stats::optim(
      par = x0, fn = fns$fn, gr = fns$gr, method = "BFGS",
      control = list(maxit = maxit, reltol = 1e-12)
    ), error = function(e) NULL)
    elapsed <- proc.time()[["elapsed"]] - t0
    if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                  iter = NA_integer_, time = elapsed,
                                  conv = FALSE, eig_min = NA_real_,
                                  par = NULL))
    g <- fns$gr(fit$par)
    return(list(f = fit$value, gnorm_inf = max(abs(g)),
                iter = unname(fit$counts[["gradient"]]),
                time = elapsed,
                conv = fit$convergence == 0L && max(abs(g)) < gtol,
                eig_min = NA_real_, par = fit$par))
  }
  if (method$runner == "nlminb") {
    fit <- tryCatch(stats::nlminb(
      start = x0, objective = fns$fn, gradient = fns$gr,
      control = list(iter.max = maxit, eval.max = 4L * maxit,
                     abs.tol = 1e-14, rel.tol = 1e-12)
    ), error = function(e) NULL)
    elapsed <- proc.time()[["elapsed"]] - t0
    if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                  iter = NA_integer_, time = elapsed,
                                  conv = FALSE, eig_min = NA_real_,
                                  par = NULL))
    g <- fns$gr(fit$par)
    return(list(f = fit$objective, gnorm_inf = max(abs(g)),
                iter = fit$iterations, time = elapsed,
                conv = fit$convergence == 0L && max(abs(g)) < gtol,
                eig_min = NA_real_, par = fit$par))
  }
  if (method$runner == "trust") {
    if (!requireNamespace("trust", quietly = TRUE)) {
      return(list(f = NA_real_, gnorm_inf = NA_real_,
                  iter = NA_integer_, time = 0,
                  conv = FALSE, eig_min = NA_real_, par = NULL))
    }
    objfun <- function(p) list(value = fns$fn(p),
                               gradient = fns$gr(p),
                               hessian = fns$hess_ad(p))
    fit <- tryCatch(trust::trust(objfun, x0, rinit = 1, rmax = 100,
                                 iterlim = maxit),
                    error = function(e) NULL)
    elapsed <- proc.time()[["elapsed"]] - t0
    if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                  iter = NA_integer_, time = elapsed,
                                  conv = FALSE, eig_min = NA_real_,
                                  par = NULL))
    g <- fns$gr(fit$argument)
    return(list(f = fit$value, gnorm_inf = max(abs(g)),
                iter = fit$iterations, time = elapsed,
                conv = fit$converged && max(abs(g)) < gtol,
                eig_min = NA_real_, par = fit$argument))
  }
  if (method$runner == "arcopt_exact") {
    hess_fn <- if (method$hess == "ad") fns$hess_ad else fns$hess_fd
    ctrl <- list(maxit = maxit, gtol_abs = gtol, trace = 0L)
    fit <- tryCatch(
      arcopt::arcopt(x0 = x0, fn = fns$fn, gr = fns$gr, hess = hess_fn,
                     control = ctrl),
      error = function(e) NULL
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                  iter = NA_integer_, time = elapsed,
                                  conv = FALSE, eig_min = NA_real_,
                                  par = NULL))
    g <- fit$gradient
    return(list(f = fit$value, gnorm_inf = max(abs(g)),
                iter = fit$iterations, time = elapsed,
                conv = isTRUE(fit$converged) && max(abs(g)) < gtol,
                eig_min = NA_real_, par = fit$par))
  }
  # arcopt_qn
  ctrl <- list(maxit = maxit, gtol_abs = gtol, trace = 0L,
               use_qn = TRUE, qn_method = method$qn)
  fit <- tryCatch(
    arcopt::arcopt(x0 = x0, fn = fns$fn, gr = fns$gr, control = ctrl),
    error = function(e) NULL
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                iter = NA_integer_, time = elapsed,
                                conv = FALSE, eig_min = NA_real_,
                                par = NULL))
  g <- fit$gradient
  list(f = fit$value, gnorm_inf = max(abs(g)),
       iter = fit$iterations, time = elapsed,
       conv = isTRUE(fit$converged) && max(abs(g)) < gtol,
       eig_min = NA_real_, par = fit$par)
}

progress("DPM Mixture Heckman: n_obs=%d  K_trunc=%d  K_true=%d",
         n_obs, k_trunc, k_true)
progress("  P_out=%d  P_sel=%d  n_excl=%d  maxit=%d  gtol=%.0e",
         p_out, p_sel, n_excl, maxit, gtol)

results <- list()
for (seed in seq_len(n_seeds)) {
  data <- simulate_dpm_heckman(n_obs = n_obs, k_true = k_true,
                               p_out = p_out, p_sel = p_sel,
                               n_excl = n_excl, seed = 100L + seed)
  fns <- make_dpm_heckman_fns(data, model_so, k_trunc = k_trunc)
  x0 <- fns$x0
  f0 <- fns$fn(x0)
  progress("Seed %d: n_params=%d  sel_rate=%.2f  f(x0)=%.2f",
           seed, fns$n_params, data$selection_rate, f0)
  for (method in methods) {
    progress("  [running  ] %s", method$name)
    res <- run_method(method, fns, x0)
    # Post-hoc Hessian eigenvalue at solution, for saddle check.
    if (!is.null(res$par)) {
      ev <- tryCatch(
        min(eigen(fns$hess_ad(res$par), symmetric = TRUE,
                  only.values = TRUE)$values),
        error = function(e) NA_real_
      )
      res$eig_min <- ev
    }
    status <- if (is.na(res$eig_min)) "?" else
      if (res$eig_min < -1e-6) "saddle" else "local min"
    results[[length(results) + 1L]] <- data.frame(
      seed = seed, method = method$name,
      f = res$f, gnorm_inf = res$gnorm_inf, iter = res$iter,
      time_s = round(res$time, 2),
      converged = res$conv,
      eig_min = signif(res$eig_min, 3),
      status = status,
      stringsAsFactors = FALSE
    )
    progress("  %-14s f=%10.3f  |g|_inf=%7.2e  iter=%5s  time=%6.2fs  eig_min=%8.2e  conv=%s  %s",
             method$name, res$f, res$gnorm_inf,
             ifelse(is.na(res$iter), "NA", format(res$iter)),
             res$time, res$eig_min, res$conv, status)
  }
  cat("\n")
}

all_results <- do.call(rbind, results)

cat("\n=============================================================\n")
cat(sprintf("Summary across %d seed(s)\n", n_seeds))
cat("=============================================================\n")
summary_tbl <- do.call(rbind, lapply(split(all_results, all_results$method),
  function(g) {
    n_sad <- sum(g$status == "saddle", na.rm = TRUE)
    data.frame(
      method = g$method[[1]],
      n_conv = sum(g$converged, na.rm = TRUE),
      n_local_min = sum(g$status == "local min", na.rm = TRUE),
      n_saddle = n_sad,
      n_total = nrow(g),
      median_f = round(stats::median(g$f, na.rm = TRUE), 2),
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

out_dir <- file.path("benchmarks", "dpm_heckman", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir,
                      sprintf("results_N%d_K%d_Ktrue%d_seeds%d.csv",
                              n_obs, k_trunc, k_true, n_seeds))
utils::write.csv(all_results, file = out_file, row.names = FALSE)
cat(sprintf("\nResults: %s\n", out_file))
