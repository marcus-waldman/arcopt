# Sanity-check the analytic gradient and Hessian for Rasch JML.

source("benchmarks/utils.R")
source("benchmarks/rasch_jml/setup.R")
devtools::load_all(".", quiet = TRUE)

cat("====================================================\n")
cat("Rasch JML: derivative checks\n")
cat("====================================================\n\n")

# Use a small problem for FD checks (full N=400 J=100 is too many params
# for finite-difference Hessian via numDeriv).
data_small <- simulate_rasch(n_persons = 20L, n_items = 10L, seed = 7L)
fns <- make_rasch_fns(data_small)
n_par <- rasch_npar(20L, 10L)
cat(sprintf("Small problem: N=20, J=10, n_par=%d\n", n_par))

# Gradient checks at three points
x_true <- data_small$x_true
x_cold <- rasch_cold_start(20L, 10L)
set.seed(13L); x_pert <- x_true + stats::rnorm(n_par, sd = 0.3)

for (lbl in c("truth", "cold", "perturbed")) {
  xx <- switch(lbl, truth = x_true, cold = x_cold, perturbed = x_pert)
  g_chk <- check_gradient(fns$fn, fns$gr, xx, tol = 1e-5, h = 1e-7)
  cat(sprintf("  grad @ %-9s  max|num-ana| = %.3e   passed=%s\n",
              lbl, g_chk$max_diff, g_chk$passed))
}

# Hessian check at truth (FD of analytic gradient)
h_chk <- check_hessian(fns$gr, fns$hess, x_true, tol = 1e-4, h = 1e-6)
cat(sprintf("  hess @ truth     max|num-ana| = %.3e   passed=%s\n",
            h_chk$max_diff, h_chk$passed))

# Eigenvalue check on full-size Hessian
cat("\n---- full-size N=400, J=100 ----\n")
data_full <- simulate_rasch(n_persons = 400L, n_items = 100L, seed = 7L)
fns_full <- make_rasch_fns(data_full)
n_par_full <- rasch_npar(400L, 100L)
cat(sprintf("Full problem: N=400, J=100, n_par=%d\n", n_par_full))

x0 <- rasch_cold_start(400L, 100L)
g0 <- fns_full$gr(x0)
cat(sprintf("  ||g||_inf at cold start: %.3e\n", max(abs(g0))))
h0 <- fns_full$hess(x0)
asym <- max(abs(h0 - t(h0)))
cat(sprintf("  |H - H^T| at cold start: %.3e\n", asym))
ev_h0 <- eigen(h0, symmetric = TRUE, only.values = TRUE)$values
cat(sprintf("  Eigenvalues: lambda_min=%.3e   lambda_max=%.3e   cond=%.2e\n",
            min(ev_h0), max(ev_h0), max(ev_h0) / min(ev_h0)))

h_truth <- fns_full$hess(data_full$x_true)
ev_truth <- eigen(h_truth, symmetric = TRUE, only.values = TRUE)$values
cat(sprintf("  At truth: lambda_min=%.3e   lambda_max=%.3e   cond=%.2e\n",
            min(ev_truth), max(ev_truth),
            max(ev_truth) / min(ev_truth)))

# Quick solve to confirm convergence
cat("\n---- quick solve from cold start ----\n")
t0 <- proc.time()[["elapsed"]]
res <- arcopt(x0, fn = fns_full$fn, gr = fns_full$gr, hess = fns_full$hess,
              control = list(maxit = 200L, gtol_abs = 1e-6))
t1 <- proc.time()[["elapsed"]]
cat(sprintf("  iter=%d  hess_evals=%d  fn=%d  gr=%d  conv=%s\n",
            res$iterations, res$evaluations$hess,
            res$evaluations$fn, res$evaluations$gr, res$converged))
cat(sprintf("  |g|_inf=%.3e  fval=%.3f  mode=%s  time=%.2fs\n",
            max(abs(res$gradient)), res$value,
            res$diagnostics$solver_mode_final, t1 - t0))

# Estimated parameter recovery
err_theta <- max(abs(res$par[1:400] - data_full$theta_true))
err_beta_free <- max(abs(res$par[401:499] - data_full$beta_true[-1L]))
cat(sprintf("  max|theta_hat - theta*| = %.3e\n", err_theta))
cat(sprintf("  max|beta_hat - beta*|   = %.3e\n", err_beta_free))
