# §4.4 Rasch-JML polish go/no-go
# ==============================
#
# Compare arcopt default vs arcopt-polish on N=400, J=100 Rasch JML at a
# range of tolerances and starting points, reporting Hessian-eval saving.

source("benchmarks/rasch_jml/setup.R")
devtools::load_all(".", quiet = TRUE)

run_one <- function(data, x0, polish, ctrl_extra = list()) {
  fns <- make_rasch_fns(data)
  ctrl <- c(list(maxit = 500L, gtol_abs = 1e-8), ctrl_extra)
  if (polish) ctrl$qn_polish_enabled <- TRUE

  t0 <- proc.time()[["elapsed"]]
  res <- arcopt(x0, fn = fns$fn, gr = fns$gr, hess = fns$hess,
                control = ctrl)
  t1 <- proc.time()[["elapsed"]]

  data.frame(
    polish = polish,
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

run_config <- function(label, sim_args, ctrl_extra = list(),
                       n_reps = 5L, sparsity = NULL) {
  cat(sprintf("\n==== %s ====\n", label))
  results <- list()
  for (rep_i in seq_len(n_reps)) {
    sim_args$seed <- 100L + rep_i
    data <- do.call(simulate_rasch, sim_args)

    # Optional sparsity: zero-out a fraction of (i,j) responses (treat as
    # missing -- but Rasch JML as written treats Y as 0/1 fully observed.
    # If sparsity is non-NULL, we instead simulate fewer items per person.)
    if (!is.null(sparsity)) {
      # Each person answers `sparsity$items_per_person` items (random subset)
      n_per <- sparsity$items_per_person
      mask <- matrix(0L, data$n_persons, data$n_items)
      for (i in seq_len(data$n_persons)) {
        chosen <- sample.int(data$n_items, n_per)
        mask[i, chosen] <- 1L
      }
      # Set unobserved cells' contribution to zero by zeroing y AND
      # adjusting the fn/gr/hess to skip those cells. Easier: keep mask
      # in a closure and rebuild fns. (See sparse helpers below.)
      data$mask <- mask
    }

    x0 <- rasch_cold_start(data$n_persons, data$n_items)
    for (polish in c(FALSE, TRUE)) {
      r <- if (is.null(sparsity)) {
        run_one(data, x0, polish, ctrl_extra)
      } else {
        run_one_sparse(data, x0, polish, ctrl_extra)
      }
      r$rep <- rep_i
      results[[length(results) + 1L]] <- r
    }
  }
  do.call(rbind, results)
}

# ---- sparse variant: each person sees a subset of items --------------

make_rasch_fns_sparse <- function(data) {
  y_mat <- data$y_mat
  mask <- data$mask                          # 1 if observed, 0 otherwise
  n_persons <- data$n_persons
  n_items <- data$n_items
  yobs <- y_mat * mask
  row_sum_y <- rowSums(yobs)
  col_sum_y <- colSums(yobs)
  row_obs <- rowSums(mask)                   # number of items each person sees
  col_obs <- colSums(mask)                   # number of persons per item

  fn <- function(x) {
    pars <- rasch_unpack(x, n_persons, n_items)
    z_mat <- outer(pars$theta, pars$beta, `-`)
    sum(-pars$theta * row_sum_y) + sum(pars$beta * col_sum_y) +
      sum(mask * (log1p(exp(-abs(z_mat))) + pmax(z_mat, 0)))
  }

  gr <- function(x) {
    pars <- rasch_unpack(x, n_persons, n_items)
    z_mat <- outer(pars$theta, pars$beta, `-`)
    p_mat <- stats::plogis(z_mat) * mask
    g_theta <- rowSums(p_mat) - row_sum_y
    g_beta_full <- col_sum_y - colSums(p_mat)
    c(g_theta, g_beta_full[-1L])
  }

  hess <- function(x) {
    pars <- rasch_unpack(x, n_persons, n_items)
    z_mat <- outer(pars$theta, pars$beta, `-`)
    p_mat <- stats::plogis(z_mat)
    w_mat <- p_mat * (1 - p_mat) * mask
    d_vec <- rowSums(w_mat)
    e_vec <- colSums(w_mat)[-1L]
    n_par <- n_persons + n_items - 1L
    h_mat <- matrix(0, n_par, n_par)
    diag(h_mat)[seq_len(n_persons)] <- d_vec
    beta_idx <- n_persons + seq_len(n_items - 1L)
    diag(h_mat)[beta_idx] <- e_vec
    cross <- -w_mat[, -1L, drop = FALSE]
    h_mat[seq_len(n_persons), beta_idx] <- cross
    h_mat[beta_idx, seq_len(n_persons)] <- t(cross)
    h_mat
  }

  list(fn = fn, gr = gr, hess = hess)
}

run_one_sparse <- function(data, x0, polish, ctrl_extra = list()) {
  fns <- make_rasch_fns_sparse(data)
  ctrl <- c(list(maxit = 500L, gtol_abs = 1e-8), ctrl_extra)
  if (polish) ctrl$qn_polish_enabled <- TRUE

  t0 <- proc.time()[["elapsed"]]
  res <- arcopt(x0, fn = fns$fn, gr = fns$gr, hess = fns$hess,
                control = ctrl)
  t1 <- proc.time()[["elapsed"]]

  data.frame(
    polish = polish,
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

summarize_run <- function(df, name) {
  default_med <- stats::median(df$hess_evals[!df$polish])
  polish_med <- stats::median(df$hess_evals[df$polish])
  saving <- (default_med - polish_med) / default_med
  default_time <- stats::median(df$time_s[!df$polish])
  polish_time <- stats::median(df$time_s[df$polish])
  time_saving <- (default_time - polish_time) / default_time
  cat(sprintf("\n%s\n", name))
  cat(sprintf("  default hess (median): %.1f   polish hess (median): %.1f\n",
              default_med, polish_med))
  cat(sprintf("  median saving (Hess): %.1f%%\n", 100 * saving))
  cat(sprintf("  median saving (time): %.1f%%\n", 100 * time_saving))
  cat(sprintf("  polish switches (mean): %.2f   reverts: %.2f\n",
              mean(df$polish_switches[df$polish]),
              mean(df$polish_reverts[df$polish])))
  cat(sprintf("  fval agreement (max diff): %.3e\n",
              max(abs(df$fval[!df$polish] - df$fval[df$polish]))))
  invisible(saving)
}

# ---- runs ----

cat("Running headline (N=400, J=100, dense, gtol=1e-8)...\n")
res_dense <- run_config(
  "Headline: dense, N=400, J=100, gtol=1e-8",
  sim_args = list(n_persons = 400L, n_items = 100L), n_reps = 5L
)
sav_dense <- summarize_run(res_dense, "dense")

cat("\nRunning tighter tol (gtol=1e-10)...\n")
res_tight <- run_config(
  "Tight tol: dense, N=400, J=100, gtol=1e-10",
  sim_args = list(n_persons = 400L, n_items = 100L),
  ctrl_extra = list(gtol_abs = 1e-10, ftol_abs = 1e-13, xtol_abs = 1e-13),
  n_reps = 5L
)
sav_tight <- summarize_run(res_tight, "tight tol")

cat("\nRunning sparse (each person sees 10/100 items)...\n")
res_sparse <- run_config(
  "Sparse: 10 of 100 items per person",
  sim_args = list(n_persons = 400L, n_items = 100L),
  sparsity = list(items_per_person = 10L),
  n_reps = 5L
)
sav_sparse <- summarize_run(res_sparse, "sparse")

dir.create("benchmarks/rasch_jml/results", showWarnings = FALSE,
           recursive = TRUE)
utils::write.csv(res_dense,
  "benchmarks/rasch_jml/results/dense_n400_j100.csv", row.names = FALSE)
utils::write.csv(res_tight,
  "benchmarks/rasch_jml/results/tight_n400_j100.csv", row.names = FALSE)
utils::write.csv(res_sparse,
  "benchmarks/rasch_jml/results/sparse_n400_j100.csv", row.names = FALSE)

cat("\n\n====================================================\n")
cat("Go/no-go verdict (>= 30% Hessian-eval saving)\n")
cat("====================================================\n")
verdict <- function(sav, lbl) {
  cat(sprintf("%-15s saving=%.1f%%   %s\n",
              lbl, 100 * sav,
              if (sav >= 0.30) "PASS" else "FAIL"))
}
verdict(sav_dense, "dense, gtol=1e-8")
verdict(sav_tight, "dense, gtol=1e-10")
verdict(sav_sparse, "sparse 10/100")
