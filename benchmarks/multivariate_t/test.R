# FD-check the analytic gradient on a few non-degenerate points and run
# a sanity solve to confirm we land near the truth.

source("benchmarks/utils.R")
source("benchmarks/multivariate_t/setup.R")

devtools::load_all(".", quiet = TRUE)

cat("====================================================\n")
cat("Multivariate-t MLE: setup sanity checks\n")
cat("====================================================\n\n")

for (p_dim in c(3L, 5L)) {
  cat(sprintf("---- p = %d (n_par = %d) ----\n", p_dim, mvt_npar(p_dim)))

  data <- simulate_mvt(p = p_dim, n = 200L, nu = 5, seed = 7L)
  fns <- make_mvt_fns(data)

  # Check 1: gradient at the truth
  g_at_truth <- check_gradient(fns$fn, fns$gr, data$x_true,
                               tol = 1e-4, h = 1e-7)
  cat(sprintf("  grad @ truth      max|num - ana| = %.3e   passed=%s\n",
              g_at_truth$max_diff, g_at_truth$passed))

  # Check 2: gradient at the MoM start
  x0 <- mvt_start(data)
  g_at_start <- check_gradient(fns$fn, fns$gr, x0,
                               tol = 1e-4, h = 1e-7)
  cat(sprintf("  grad @ start      max|num - ana| = %.3e   passed=%s\n",
              g_at_start$max_diff, g_at_start$passed))

  # Check 3: gradient at a perturbed point
  set.seed(13L)
  x_pert <- data$x_true + stats::rnorm(length(data$x_true), sd = 0.2)
  g_at_pert <- check_gradient(fns$fn, fns$gr, x_pert,
                              tol = 1e-4, h = 1e-7)
  cat(sprintf("  grad @ perturbed  max|num - ana| = %.3e   passed=%s\n",
              g_at_pert$max_diff, g_at_pert$passed))

  # Check 4: Hessian symmetric and PD at the truth (large-sample)
  h_at_truth <- fns$hess(data$x_true)
  asym <- max(abs(h_at_truth - t(h_at_truth)))
  eig_min <- min(eigen(h_at_truth, symmetric = TRUE,
                       only.values = TRUE)$values)
  cat(sprintf("  hess @ truth      |H - H^T| = %.3e   lambda_min = %.3e\n",
              asym, eig_min))

  # Check 5: solve from MoM start; check final |grad|, distance to truth
  res <- arcopt(x0, fn = fns$fn, gr = fns$gr, hess = fns$hess,
                control = list(maxit = 200L, gtol_abs = 1e-6))
  err_mu <- max(abs(res$par[seq_len(p_dim)] - data$mu_true))
  err_nu <- abs(exp(res$par[length(res$par)]) - data$nu_true)
  cat(sprintf("  solve: %3d iter, conv=%s, |g|_inf=%.2e, max|mu-mu*|=%.2e, |nu-nu*|=%.2e\n",
              res$iterations, res$converged,
              max(abs(res$gradient)), err_mu, err_nu))
  cat(sprintf("         hess_evals=%d  mode=%s\n",
              res$evaluations$hess, res$diagnostics$solver_mode_final))
  cat("\n")
}

cat("Done.\n")
