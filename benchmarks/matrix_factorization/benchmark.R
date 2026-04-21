# Large-scale matrix factorization benchmark
# ==========================================
#
# Compares arcopt (exact + QN variants) and optim BFGS / nlminb on a
# 500-parameter symmetric matrix factorization problem. The instance is
# d = 50, r = 10 so U is 50 x 10 and U U' - M is 50 x 50. Ground truth
# M = U_true U_true' has rank r, so f* = 0.
#
# Usage (from repo root):
#   Rscript benchmarks/matrix_factorization/benchmark.R [n_seeds]

suppressPackageStartupMessages({
  library(arcopt)
})

source(file.path("benchmarks", "matrix_factorization", "setup.R"))

args <- commandArgs(trailingOnly = TRUE)
n_seeds <- if (length(args) >= 1L) as.integer(args[[1]]) else 5L

d <- 50L
r <- 10L
maxit <- 2000L
gtol <- 1e-5

cat(sprintf("Matrix factorization benchmark: d=%d, r=%d, n=%d\n", d, r, d * r))
cat(sprintf("Seeds: %d    maxit: %d    gtol: %.0e\n\n",
            n_seeds, maxit, gtol))

prob <- make_matrix_factorization_problem(d = d, r = r, u_true_seed = 1L)

methods <- list(
  list(name = "arcopt-exact", use_qn = FALSE, qn_method = NA_character_),
  list(name = "arcopt-hybrid", use_qn = TRUE, qn_method = "hybrid"),
  list(name = "arcopt-sr1", use_qn = TRUE, qn_method = "sr1"),
  list(name = "arcopt-bfgs", use_qn = TRUE, qn_method = "bfgs"),
  list(name = "optim-BFGS", use_qn = NA, qn_method = NA_character_)
)

run_method <- function(method, x0) {
  t0 <- proc.time()[["elapsed"]]
  if (method$name == "optim-BFGS") {
    fit <- tryCatch(
      stats::optim(
        par = x0, fn = prob$fn, gr = prob$gr, method = "BFGS",
        control = list(maxit = maxit, reltol = 1e-10)
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(list(f = NA_real_, gnorm = NA_real_, iter = NA_integer_,
                  time = proc.time()[["elapsed"]] - t0, conv = FALSE))
    }
    return(list(
      f = fit$value,
      gnorm = sqrt(sum(prob$gr(fit$par)^2)),
      iter = fit$counts[["function"]],
      time = proc.time()[["elapsed"]] - t0,
      conv = fit$convergence == 0L
    ))
  }

  ctrl <- list(maxit = maxit, gtol_abs = gtol, trace = 0L)
  if (isTRUE(method$use_qn)) {
    ctrl$use_qn <- TRUE
    ctrl$qn_method <- method$qn_method
    fit <- tryCatch(
      arcopt::arcopt(x0 = x0, fn = prob$fn, gr = prob$gr, control = ctrl),
      error = function(e) NULL
    )
  } else {
    fit <- tryCatch(
      arcopt::arcopt(x0 = x0, fn = prob$fn, gr = prob$gr, hess = prob$hess,
                     control = ctrl),
      error = function(e) NULL
    )
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  if (is.null(fit)) {
    return(list(f = NA_real_, gnorm = NA_real_, iter = NA_integer_,
                time = elapsed, conv = FALSE))
  }
  list(
    f = fit$value,
    gnorm = sqrt(sum(fit$gradient^2)),
    iter = fit$iterations,
    time = elapsed,
    conv = isTRUE(fit$converged)
  )
}

results <- list()
for (seed in seq_len(n_seeds)) {
  x0 <- make_start_point(d = d, r = r, seed = 1000L + seed)
  cat(sprintf("Seed %d:\n", seed))
  for (method in methods) {
    res <- run_method(method, x0)
    results[[length(results) + 1L]] <- data.frame(
      seed = seed,
      method = method$name,
      f = res$f,
      gnorm = res$gnorm,
      iter = res$iter,
      time_s = round(res$time, 2),
      converged = res$conv,
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %-14s f=%.2e  ||g||=%.2e  iter=%5s  time=%5.2fs  conv=%s\n",
                method$name, res$f, res$gnorm,
                ifelse(is.na(res$iter), "NA", format(res$iter)),
                res$time, res$conv))
  }
  cat("\n")
}

all_results <- do.call(rbind, results)

cat("\n=============================================================\n")
cat("Summary by method (converged = ||grad|| < gtol)\n")
cat("=============================================================\n")
summary_tbl <- do.call(
  rbind,
  lapply(split(all_results, all_results$method), function(g) {
    data.frame(
      method = g$method[[1]],
      n_conv = sum(g$converged, na.rm = TRUE),
      n_total = nrow(g),
      pct_conv = 100 * mean(g$converged, na.rm = TRUE),
      median_f = stats::median(g$f, na.rm = TRUE),
      median_iter = stats::median(g$iter, na.rm = TRUE),
      median_time_s = round(stats::median(g$time_s, na.rm = TRUE), 2),
      stringsAsFactors = FALSE
    )
  })
)
summary_tbl <- summary_tbl[match(
  vapply(methods, function(m) m$name, character(1)),
  summary_tbl$method
), ]
rownames(summary_tbl) <- NULL
print(summary_tbl, row.names = FALSE)

out_dir <- file.path("benchmarks", "matrix_factorization", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(
  all_results,
  file = file.path(out_dir, sprintf("results_d%d_r%d_seeds%d.csv",
                                    d, r, n_seeds)),
  row.names = FALSE
)
cat(sprintf("\nResults written to %s\n",
            file.path(out_dir, sprintf("results_d%d_r%d_seeds%d.csv",
                                       d, r, n_seeds))))
