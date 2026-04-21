# Chebyshev-Rosenbrock benchmark (n = 500)
# ========================================
#
# Nesterov's worst-case ARC benchmark. f has 2^(n-1) critical points,
# most of which are genuine local minima. A local optimizer's
# convergence to the global minimum is therefore not guaranteed, but
# methods that handle negative curvature properly (ARC) should climb
# out of saddle basins that first-order methods get stuck in.
#
# Starting points:
#   "nesterov":    x_0 = (-1, 1, 1, ..., 1)  (Cartis-Gould-Toint 2010)
#   "all-neg":     x_0 = (-1, -1, ..., -1)
#   "origin":      x_0 = (0, 0, ..., 0)  (saddle-like: Hessian has
#                  n-1 negative eigenvalues of -6 to -7.5 and one +2)
#   "random":      x_0 ~ N(0, sd^2)
#
# Usage (from repo root):
#   Rscript benchmarks/chebyshev_rosenbrock/benchmark.R [n] [kind] [n_seeds]
#
# For deterministic starts ("nesterov" or "all-neg"), n_seeds has no
# effect beyond running the same case multiple times.

suppressPackageStartupMessages({
  library(arcopt)
})

source(file.path("benchmarks", "chebyshev_rosenbrock", "setup.R"))

args <- commandArgs(trailingOnly = TRUE)
n     <- if (length(args) >= 1L) as.integer(args[[1]]) else 500L
kind  <- if (length(args) >= 2L) args[[2]]            else "nesterov"
n_seeds <- if (length(args) >= 3L) as.integer(args[[3]]) else 1L

maxit <- 2000L
gtol <- 1e-5

prob <- make_chebrosen_problem(n = n)

cat(sprintf("Chebyshev-Rosenbrock: n = %d  start = %s  n_seeds = %d\n",
            n, kind, n_seeds))
cat(sprintf("  maxit = %d   gtol = %.0e\n\n", maxit, gtol))

methods <- list(
  list(name = "arcopt-exact",  runner = "arcopt_exact"),
  list(name = "arcopt-hybrid", runner = "arcopt_qn", qn_method = "hybrid"),
  list(name = "arcopt-sr1",    runner = "arcopt_qn", qn_method = "sr1"),
  list(name = "arcopt-bfgs",   runner = "arcopt_qn", qn_method = "bfgs"),
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
    elapsed <- proc.time()[["elapsed"]] - t0
    if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                  iter = NA_integer_, time = elapsed,
                                  conv = FALSE))
    g <- prob$gr(fit$par)
    return(list(f = fit$value, gnorm_inf = max(abs(g)),
                iter = fit$counts[["gradient"]], time = elapsed,
                conv = fit$convergence == 0L && max(abs(g)) < gtol))
  }

  if (method$runner == "optim_lbfgsb") {
    fit <- tryCatch(stats::optim(
      par = x0, fn = prob$fn, gr = prob$gr, method = "L-BFGS-B",
      control = list(maxit = maxit, factr = 1e1, pgtol = gtol)
    ), error = function(e) NULL)
    elapsed <- proc.time()[["elapsed"]] - t0
    if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                  iter = NA_integer_, time = elapsed,
                                  conv = FALSE))
    g <- prob$gr(fit$par)
    return(list(f = fit$value, gnorm_inf = max(abs(g)),
                iter = fit$counts[["gradient"]], time = elapsed,
                conv = fit$convergence == 0L && max(abs(g)) < gtol))
  }

  if (method$runner == "nlminb") {
    fit <- tryCatch(stats::nlminb(
      start = x0, objective = prob$fn, gradient = prob$gr,
      control = list(iter.max = maxit, eval.max = 4L * maxit,
                     abs.tol = 1e-14, rel.tol = 1e-12)
    ), error = function(e) NULL)
    elapsed <- proc.time()[["elapsed"]] - t0
    if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                  iter = NA_integer_, time = elapsed,
                                  conv = FALSE))
    g <- prob$gr(fit$par)
    return(list(f = fit$objective, gnorm_inf = max(abs(g)),
                iter = fit$iterations, time = elapsed,
                conv = fit$convergence == 0L && max(abs(g)) < gtol))
  }

  if (method$runner == "arcopt_exact") {
    ctrl <- list(maxit = maxit, gtol_abs = gtol, trace = 0L)
    fit <- tryCatch(
      arcopt::arcopt(x0 = x0, fn = prob$fn, gr = prob$gr, hess = prob$hess,
                     control = ctrl),
      error = function(e) NULL
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                  iter = NA_integer_, time = elapsed,
                                  conv = FALSE))
    g <- fit$gradient
    return(list(f = fit$value, gnorm_inf = max(abs(g)),
                iter = fit$iterations, time = elapsed,
                conv = isTRUE(fit$converged) && max(abs(g)) < gtol))
  }

  # arcopt_qn
  ctrl <- list(maxit = maxit, gtol_abs = gtol, trace = 0L,
               use_qn = TRUE, qn_method = method$qn_method)
  fit <- tryCatch(
    arcopt::arcopt(x0 = x0, fn = prob$fn, gr = prob$gr, control = ctrl),
    error = function(e) NULL
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                iter = NA_integer_, time = elapsed,
                                conv = FALSE))
  g <- fit$gradient
  list(f = fit$value, gnorm_inf = max(abs(g)),
       iter = fit$iterations, time = elapsed,
       conv = isTRUE(fit$converged) && max(abs(g)) < gtol)
}

results <- list()
for (seed in seq_len(n_seeds)) {
  if (kind == "nesterov") {
    x0 <- c(-1, rep(1, n - 1L))
  } else if (kind == "all-neg") {
    x0 <- rep(-1, n)
  } else if (kind == "origin") {
    x0 <- rep(0, n)
  } else {
    x0 <- make_chebrosen_start(n, kind = "random", seed = 100L + seed)
  }
  f0 <- prob$fn(x0)
  cat(sprintf("Seed %d (f0 = %.2e):\n", seed, f0))
  for (method in methods) {
    res <- run_method(method, x0)
    results[[length(results) + 1L]] <- data.frame(
      seed = seed, method = method$name, f = res$f,
      gnorm_inf = res$gnorm_inf, iter = res$iter,
      time_s = round(res$time, 2), converged = res$conv,
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %-14s f=%10.3e  |g|_inf=%8.2e  iter=%5s  time=%7.2fs  conv=%s\n",
                method$name, res$f, res$gnorm_inf,
                ifelse(is.na(res$iter), "NA", format(res$iter)),
                res$time, res$conv))
  }
  cat("\n")
}

all_results <- do.call(rbind, results)

cat("\n=============================================================\n")
cat(sprintf("Summary across %d runs (converged = |g|_inf < %.0e AND f < 1e-6)\n",
            n_seeds, gtol))
cat("=============================================================\n")
all_results$reached_opt <- all_results$converged &
  !is.na(all_results$f) & all_results$f < 1e-6
summary_tbl <- do.call(rbind, lapply(split(all_results, all_results$method),
  function(g) {
    data.frame(
      method = g$method[[1]],
      n_opt = sum(g$reached_opt, na.rm = TRUE),
      n_total = nrow(g),
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

out_dir <- file.path("benchmarks", "chebyshev_rosenbrock", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir,
                      sprintf("results_n%d_%s_seeds%d.csv",
                              n, kind, n_seeds))
utils::write.csv(all_results, file = out_file, row.names = FALSE)
cat(sprintf("\nResults: %s\n", out_file))
