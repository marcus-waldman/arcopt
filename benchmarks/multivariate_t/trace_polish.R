# Look inside one run: print per-iteration ||g||_inf, rho, lambda_min,
# sigma to see exactly which polish-detector signal is failing.

source("benchmarks/multivariate_t/setup.R")
devtools::load_all(".", quiet = TRUE)

data <- simulate_mvt(p = 5L, n = 200L, nu = 5, seed = 101L)
fns <- make_mvt_fns(data)
x0 <- mvt_pack(rep(0, 5L), diag(5L), nu = 30)

res <- arcopt(x0, fn = fns$fn, gr = fns$gr, hess = fns$hess,
              control = list(maxit = 500L, gtol_abs = 1e-8,
                             qn_polish_enabled = TRUE,
                             trace = 3L))

cat("\n--- iteration trace ---\n")
tr <- res$trace
cat(sprintf("Iterations recorded: %d\n", length(tr$f)))
trace_df <- data.frame(
  iter = seq_along(tr$f) - 1L,
  f = tr$f,
  g_inf = tr$g_norm,
  sigma = tr$sigma,
  rho = tr$rho,
  step_type = tr$step_type,
  lambda_min = if (!is.null(tr$lambda)) tr$lambda else NA
)
print(trace_df, row.names = FALSE)

cat("\n--- diagnostics ---\n")
print(res$diagnostics)
cat(sprintf("\nIterations=%d  hess_evals=%d  fn_evals=%d  gr_evals=%d\n",
            res$iterations, res$evaluations$hess,
            res$evaluations$fn, res$evaluations$gr))
