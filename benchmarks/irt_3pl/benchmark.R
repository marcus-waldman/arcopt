# 3PL IRT Marginal MLE benchmark
# ==============================
#
# Compare arcopt-exact, arcopt-hybrid, arcopt-sr1, optim-BFGS, and
# nlminb on 3PL MML without a prior on c. Default problem size is
# J = 10 items, N = 1000 examinees -> n = 30 parameters.
#
# Hard items first ("hard"): b_i in [0.5, 2.5], most examinees far
# below b so c_i is well-sampled and estimable.
#
# Usage (from repo root):
#   Rscript benchmarks/irt_3pl/benchmark.R [items] [n_seeds] [j_items] [n_ex]

suppressPackageStartupMessages({
  library(arcopt)
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("package 'numDeriv' is required")
  }
})

source(file.path("benchmarks", "irt_3pl", "setup.R"))

args <- commandArgs(trailingOnly = TRUE)
items    <- if (length(args) >= 1L) args[[1]]             else "hard"
n_seeds  <- if (length(args) >= 2L) as.integer(args[[2]]) else 5L
j_items  <- if (length(args) >= 3L) as.integer(args[[3]]) else 10L
n_ex     <- if (length(args) >= 4L) as.integer(args[[4]]) else 1000L

maxit <- 2000L
gtol <- 1e-5

methods <- list(
  list(name = "arcopt-exact",  runner = "arcopt_exact"),
  list(name = "arcopt-hybrid", runner = "arcopt_qn", qn_method = "hybrid"),
  list(name = "arcopt-sr1",    runner = "arcopt_qn", qn_method = "sr1"),
  list(name = "optim-BFGS",    runner = "optim_bfgs"),
  list(name = "nlminb",        runner = "nlminb")
)

standard_start <- function(data) {
  p_item <- colMeans(data$y)   # proportion correct
  p_item <- pmin(pmax(p_item, 0.05), 0.95)  # clip to avoid +-Inf
  # a_i = 1, b_i = -logit(p_i), c_i = 0.2
  irt3pl_pack(a = rep(1, data$j_items),
              b = -stats::qlogis(p_item),
              c = rep(0.2, data$j_items))
}

run_method <- function(method, fn, gr, hess, x0) {
  t0 <- proc.time()[["elapsed"]]

  if (method$runner == "optim_bfgs") {
    fit <- tryCatch(stats::optim(
      par = x0, fn = fn, gr = gr, method = "BFGS",
      control = list(maxit = maxit, reltol = 1e-12)
    ), error = function(e) NULL)
    elapsed <- proc.time()[["elapsed"]] - t0
    if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                  iter = NA_integer_, time = elapsed,
                                  conv = FALSE, par = NULL))
    g <- gr(fit$par)
    return(list(f = fit$value, gnorm_inf = max(abs(g)),
                iter = fit$counts[["gradient"]], time = elapsed,
                conv = fit$convergence == 0L && max(abs(g)) < gtol,
                par = fit$par))
  }

  if (method$runner == "nlminb") {
    fit <- tryCatch(stats::nlminb(
      start = x0, objective = fn, gradient = gr,
      control = list(iter.max = maxit, eval.max = 4L * maxit,
                     abs.tol = 1e-14, rel.tol = 1e-12)
    ), error = function(e) NULL)
    elapsed <- proc.time()[["elapsed"]] - t0
    if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                  iter = NA_integer_, time = elapsed,
                                  conv = FALSE, par = NULL))
    g <- gr(fit$par)
    return(list(f = fit$objective, gnorm_inf = max(abs(g)),
                iter = fit$iterations, time = elapsed,
                conv = fit$convergence == 0L && max(abs(g)) < gtol,
                par = fit$par))
  }

  if (method$runner == "arcopt_exact") {
    ctrl <- list(maxit = maxit, gtol_abs = gtol, trace = 0L)
    fit <- tryCatch(
      arcopt::arcopt(x0 = x0, fn = fn, gr = gr, hess = hess,
                     control = ctrl),
      error = function(e) NULL
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                  iter = NA_integer_, time = elapsed,
                                  conv = FALSE, par = NULL))
    g <- fit$gradient
    return(list(f = fit$value, gnorm_inf = max(abs(g)),
                iter = fit$iterations, time = elapsed,
                conv = isTRUE(fit$converged) && max(abs(g)) < gtol,
                par = fit$par))
  }

  ctrl <- list(maxit = maxit, gtol_abs = gtol, trace = 0L,
               use_qn = TRUE, qn_method = method$qn_method)
  fit <- tryCatch(
    arcopt::arcopt(x0 = x0, fn = fn, gr = gr, control = ctrl),
    error = function(e) NULL
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  if (is.null(fit)) return(list(f = NA_real_, gnorm_inf = NA_real_,
                                iter = NA_integer_, time = elapsed,
                                conv = FALSE, par = NULL))
  g <- fit$gradient
  list(f = fit$value, gnorm_inf = max(abs(g)),
       iter = fit$iterations, time = elapsed,
       conv = isTRUE(fit$converged) && max(abs(g)) < gtol,
       par = fit$par)
}

cat(sprintf("3PL MML benchmark: items=%s  J=%d  N=%d  seeds=%d  n_params=%d\n",
            items, j_items, n_ex, n_seeds, 3L * j_items))
cat(sprintf("  maxit=%d  gtol=%.0e\n\n", maxit, gtol))

results <- list()
for (seed in seq_len(n_seeds)) {
  data <- simulate_irt3pl(j_items = j_items, n_examinees = n_ex,
                          items = items, seed = 1000L + seed)
  fns <- make_irt3pl_fns(data, q_nodes = 40L)
  fn <- fns$fn
  gr <- fns$gr
  hess <- function(x) numDeriv::jacobian(gr, x, method = "Richardson")

  x0 <- standard_start(data)
  f_true <- fn(data$x_true)
  f_x0 <- fn(x0)
  cat(sprintf("Seed %d: f(truth) = %.2f  f(x0) = %.2f\n",
              seed, f_true, f_x0))
  for (method in methods) {
    res <- run_method(method, fn, gr, hess, x0)
    c_est <- if (!is.null(res$par)) {
      pars <- irt3pl_unpack(res$par, j_items)
      paste(sprintf("%.2f", pars$c), collapse = ",")
    } else NA_character_
    results[[length(results) + 1L]] <- data.frame(
      seed = seed, method = method$name,
      f = res$f, gnorm_inf = res$gnorm_inf, iter = res$iter,
      time_s = round(res$time, 2), converged = res$conv,
      f_minus_truth = res$f - f_true,
      c_est = c_est,
      stringsAsFactors = FALSE
    )
    cat(sprintf("  %-14s f=%10.3f  |g|_inf=%8.2e  iter=%5s  time=%6.2fs  conv=%s\n",
                method$name, res$f, res$gnorm_inf,
                ifelse(is.na(res$iter), "NA", format(res$iter)),
                res$time, res$conv))
  }
  cat("\n")
}

all_results <- do.call(rbind, results)

cat("\n=============================================================\n")
cat(sprintf("Summary across %d seeds\n", n_seeds))
cat("=============================================================\n")
summary_tbl <- do.call(rbind, lapply(split(all_results, all_results$method),
  function(g) {
    data.frame(
      method = g$method[[1]],
      n_conv = sum(g$converged, na.rm = TRUE),
      n_total = nrow(g),
      median_f = round(stats::median(g$f, na.rm = TRUE), 3),
      median_f_minus_truth = round(stats::median(g$f_minus_truth,
                                                 na.rm = TRUE), 3),
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

out_dir <- file.path("benchmarks", "irt_3pl", "results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir,
                      sprintf("results_%s_J%d_N%d_seeds%d.csv",
                              items, j_items, n_ex, n_seeds))
utils::write.csv(all_results, file = out_file, row.names = FALSE)
cat(sprintf("\nResults: %s\n", out_file))
