## cache_for_manuscript.R
## ======================
## Runs the four-method comparison (nlminb, trust, arcopt default,
## arcopt-polish) on 100 random sparse-Rasch JML datasets at the
## §4.4 headline configuration (N = 400 participants, J = 100 items,
## 10 items per participant, gtol_abs = 1e-8, cold start) and writes
## the per-seed records and aggregated medians/IQRs to
## manuscript/_cache/rasch_jml_results.rds.
##
## Total runtime on a 2026 workstation is roughly 8 minutes, dominated
## by the 100 trust calls (~3 s each at n_par = 499).
##
## Usage (from repo root):
##   Rscript benchmarks/rasch_jml/cache_for_manuscript.R

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(arcopt)
  }
})

source(file.path("benchmarks", "rasch_jml", "setup.R"))

# Configuration
n_seeds   <- 100L
n_persons <- 400L
n_items   <- 100L
items_per_person <- 10L
gtol      <- 1e-8
maxit     <- 500L

build_sparse_data <- function(seed) {
  set.seed(seed)
  data <- simulate_rasch(n_persons = n_persons, n_items = n_items,
                         seed = seed)
  mask <- matrix(0L, n_persons, n_items)
  for (i in seq_len(n_persons)) {
    mask[i, sample.int(n_items, items_per_person)] <- 1L
  }
  data$mask <- mask
  data
}

make_sparse_fns <- function(data) {
  y_mat <- data$y_mat
  mask_local <- data$mask
  np <- data$n_persons
  ni <- data$n_items
  yobs <- y_mat * mask_local
  row_sum_y <- rowSums(yobs)
  col_sum_y <- colSums(yobs)

  fn <- function(x) {
    pars <- rasch_unpack(x, np, ni)
    z_mat <- outer(pars$theta, pars$beta, `-`)
    sum(-pars$theta * row_sum_y) + sum(pars$beta * col_sum_y) +
      sum(mask_local * (log1p(exp(-abs(z_mat))) + pmax(z_mat, 0)))
  }
  gr <- function(x) {
    pars <- rasch_unpack(x, np, ni)
    z_mat <- outer(pars$theta, pars$beta, `-`)
    p_mat <- stats::plogis(z_mat) * mask_local
    g_theta <- rowSums(p_mat) - row_sum_y
    g_beta_full <- col_sum_y - colSums(p_mat)
    c(g_theta, g_beta_full[-1L])
  }
  hess <- function(x) {
    pars <- rasch_unpack(x, np, ni)
    z_mat <- outer(pars$theta, pars$beta, `-`)
    p_mat <- stats::plogis(z_mat)
    w_mat <- p_mat * (1 - p_mat) * mask_local
    d_vec <- rowSums(w_mat)
    e_vec <- colSums(w_mat)[-1L]
    n_par <- np + ni - 1L
    h_mat <- matrix(0, n_par, n_par)
    diag(h_mat)[seq_len(np)] <- d_vec
    beta_idx <- np + seq_len(ni - 1L)
    diag(h_mat)[beta_idx] <- e_vec
    cross <- -w_mat[, -1L, drop = FALSE]
    h_mat[seq_len(np), beta_idx] <- cross
    h_mat[beta_idx, seq_len(np)] <- t(cross)
    h_mat
  }
  list(fn = fn, gr = gr, hess = hess)
}

run_methods <- function(fns, x0) {
  results <- list()

  # nlminb (PORT, gradient + Hessian)
  t0 <- proc.time()[["elapsed"]]
  r1 <- stats::nlminb(x0, fns$fn, fns$gr, fns$hess,
                      control = list(iter.max = maxit, eval.max = maxit * 2L,
                                     rel.tol = 1e-10))
  results$nlminb <- list(
    nll = r1$objective, g_inf = max(abs(fns$gr(r1$par))),
    iter = unname(r1$evaluations[1]),
    hess = unname(r1$evaluations[1]),  # nlminb doesn't separate; treat each iter as one Hess call
    time_s = proc.time()[["elapsed"]] - t0,
    converged = r1$convergence == 0L,
    notes = if (r1$convergence == 0L) "" else r1$message
  )

  # trust (Hessian-using trust region)
  t0 <- proc.time()[["elapsed"]]
  obj <- function(p) list(value = fns$fn(p), gradient = fns$gr(p),
                          hessian = fns$hess(p))
  r2 <- trust::trust(obj, x0, rinit = 1, rmax = 100, iterlim = maxit)
  results$trust <- list(
    nll = r2$value, g_inf = max(abs(r2$gradient)),
    iter = r2$iterations, hess = r2$iterations,
    time_s = proc.time()[["elapsed"]] - t0,
    converged = isTRUE(r2$converged), notes = ""
  )

  # arcopt default
  t0 <- proc.time()[["elapsed"]]
  r3 <- arcopt::arcopt(x0, fns$fn, fns$gr, fns$hess,
                       control = list(maxit = maxit, gtol_abs = gtol))
  results$arcopt <- list(
    nll = r3$value, g_inf = max(abs(r3$gradient)),
    iter = r3$iterations, hess = r3$evaluations$hess,
    time_s = proc.time()[["elapsed"]] - t0,
    converged = isTRUE(r3$converged) && r3$message != "max_iter",
    mode = r3$diagnostics$solver_mode_final,
    notes = sprintf("mode=%s", r3$diagnostics$solver_mode_final)
  )

  # arcopt with polish
  t0 <- proc.time()[["elapsed"]]
  r4 <- arcopt::arcopt(x0, fns$fn, fns$gr, fns$hess,
                       control = list(maxit = maxit, gtol_abs = gtol,
                                      qn_polish_enabled = TRUE))
  results$`arcopt-polish` <- list(
    nll = r4$value, g_inf = max(abs(r4$gradient)),
    iter = r4$iterations, hess = r4$evaluations$hess,
    time_s = proc.time()[["elapsed"]] - t0,
    converged = isTRUE(r4$converged) && r4$message != "max_iter",
    mode = r4$diagnostics$solver_mode_final,
    polish_switches = r4$diagnostics$qn_polish_switches,
    polish_reverts = r4$diagnostics$qn_polish_reverts,
    notes = sprintf("mode=%s, switches=%d",
                    r4$diagnostics$solver_mode_final,
                    r4$diagnostics$qn_polish_switches)
  )

  results
}

# --- main loop ---

t_start <- proc.time()[["elapsed"]]
seeds <- 100L + seq_len(n_seeds)
records <- vector("list", length(seeds))

for (i in seq_along(seeds)) {
  seed <- seeds[i]
  data <- build_sparse_data(seed)
  fns <- make_sparse_fns(data)
  x0 <- rasch_cold_start(n_persons, n_items)
  res <- run_methods(fns, x0)
  for (m in names(res)) {
    res[[m]]$seed <- seed
    res[[m]]$method <- m
  }
  records[[i]] <- res
  if (i %% 10L == 0L) {
    cat(sprintf("[cache] %d / %d seeds done (%.1f min elapsed)\n",
                i, length(seeds),
                (proc.time()[["elapsed"]] - t_start) / 60))
  }
}

# Flatten to a long data frame
flat_rows <- list()
for (rec in records) {
  for (m in names(rec)) {
    r <- rec[[m]]
    flat_rows[[length(flat_rows) + 1L]] <- data.frame(
      seed = r$seed, method = r$method,
      iter = r$iter, hess = r$hess,
      time_s = r$time_s, nll = r$nll, g_inf = r$g_inf,
      converged = r$converged,
      polish_switches = if (!is.null(r$polish_switches)) r$polish_switches
                        else NA_integer_,
      stringsAsFactors = FALSE
    )
  }
}
df_long <- do.call(rbind, flat_rows)

# Summary stats per method
summary_method <- function(m, col) {
  v <- df_long[[col]][df_long$method == m]
  c(median = stats::median(v, na.rm = TRUE),
    q25 = stats::quantile(v, 0.25, na.rm = TRUE, names = FALSE),
    q75 = stats::quantile(v, 0.75, na.rm = TRUE, names = FALSE),
    mean = mean(v, na.rm = TRUE),
    sd = stats::sd(v, na.rm = TRUE))
}
methods <- unique(df_long$method)
summary_tbl <- list()
for (m in methods) {
  summary_tbl[[m]] <- list(
    iter = summary_method(m, "iter"),
    hess = summary_method(m, "hess"),
    time_s = summary_method(m, "time_s"),
    nll = summary_method(m, "nll"),
    g_inf = summary_method(m, "g_inf"),
    converged_n = sum(df_long$converged[df_long$method == m], na.rm = TRUE),
    n_seeds = sum(df_long$method == m)
  )
}

# Polish-saving distribution: per-seed (default - polish) / default
saving <- merge(
  df_long[df_long$method == "arcopt", c("seed", "hess")],
  df_long[df_long$method == "arcopt-polish", c("seed", "hess")],
  by = "seed", suffixes = c("_default", "_polish"))
saving$rel_saving <- (saving$hess_default - saving$hess_polish) /
                     saving$hess_default

# Time saving similarly
time_save <- merge(
  df_long[df_long$method == "arcopt", c("seed", "time_s")],
  df_long[df_long$method == "arcopt-polish", c("seed", "time_s")],
  by = "seed", suffixes = c("_default", "_polish"))
time_save$rel_saving <- (time_save$time_s_default -
                         time_save$time_s_polish) /
                        time_save$time_s_default

# NLL agreement
nll_diff <- merge(
  df_long[df_long$method == "arcopt", c("seed", "nll")],
  df_long[df_long$method == "arcopt-polish", c("seed", "nll")],
  by = "seed", suffixes = c("_default", "_polish"))
nll_diff$abs_diff <- abs(nll_diff$nll_default - nll_diff$nll_polish)

results <- list(
  meta = list(
    n_seeds = n_seeds, n_persons = n_persons, n_items = n_items,
    items_per_person = items_per_person,
    gtol = gtol, maxit = maxit,
    n_par = n_persons + n_items - 1L,
    timestamp = Sys.time(),
    elapsed_total_s = proc.time()[["elapsed"]] - t_start
  ),
  records = records,
  long = df_long,
  summary = summary_tbl,
  hess_saving = saving,
  time_saving = time_save,
  nll_agreement = nll_diff
)

out_dir <- file.path("manuscript", "_cache")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_file <- file.path(out_dir, "rasch_jml_results.rds")
saveRDS(results, out_file)

cat(sprintf("\n[cache] wrote %s\n", out_file))
cat(sprintf("[cache] total elapsed: %.1f min\n",
            (proc.time()[["elapsed"]] - t_start) / 60))
cat(sprintf("[cache] median Hess saving: %.1f%%   (q25-q75: %.1f%% to %.1f%%)\n",
            100 * stats::median(saving$rel_saving),
            100 * stats::quantile(saving$rel_saving, 0.25, names = FALSE),
            100 * stats::quantile(saving$rel_saving, 0.75, names = FALSE)))
cat(sprintf("[cache] pass rate (Hess saving > 30%%): %d / %d\n",
            sum(saving$rel_saving > 0.30), nrow(saving)))
