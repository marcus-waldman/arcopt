# Validate Stan-derived fn/gr/hess for DPM mixture Heckman
# ========================================================

source("setup.R")

cat("Simulating small instance for gradient check\n")
data <- simulate_dpm_heckman(n_obs = 100L, k_true = 2L,
                             p_out = 3L, p_sel = 5L, n_excl = 2L,
                             seed = 1L)
cat("  N =", data$n_obs,
    "  P_out =", data$p_out,
    "  P_sel =", data$p_sel,
    "  selection rate =", round(data$selection_rate, 3), "\n")

cat("\nCompiling Stan model...\n")
model_so <- compile_dpm_heckman_model(
  stan_file = file.path("..", "dpm_heckman", "dpm_heckman.stan")
)

cat("\nBuilding fn/gr/hess...\n")
fns <- make_dpm_heckman_fns(data, model_so, k_trunc = 5L,
                            sigma_beta = 5.0, sigma_gamma = 5.0,
                            alpha = 1.0)
cat("  n_params (unconstrained) =", fns$n_params, "\n")

cat("\nAt x0:\n")
cat("  f(x0) =", fns$fn(fns$x0), "\n")
g0 <- fns$gr(fns$x0)
cat("  ||g(x0)|| =", sqrt(sum(g0^2)),
    "  max|g| =", max(abs(g0)), "\n")

cat("\nGradient check via numDeriv::grad (slow-but-sure)...\n")
g_fd <- numDeriv::grad(fns$fn, fns$x0, method = "Richardson")
max_diff <- max(abs(g0 - g_fd))
cat("  max |grad_AD - grad_FD| =", format(max_diff, digits = 3), "\n")
cat("  passed:", max_diff < 1e-5, "\n")

cat("\nHessian checks:\n")
h_ad <- fns$hess_ad(fns$x0)
h_fd <- fns$hess_fd(fns$x0)
cat("  dim(H_AD) =", dim(h_ad), "\n")
cat("  H_AD symmetric (tol 1e-8):", isSymmetric(h_ad, tol = 1e-8), "\n")
max_h_diff <- max(abs(h_ad - h_fd))
cat("  max |H_AD - H_FD| =", format(max_h_diff, digits = 3), "\n")

eig_ad <- eigen(h_ad, symmetric = TRUE, only.values = TRUE)$values
cat("  eigenvalue range(H_AD): [",
    format(min(eig_ad), digits = 3), ",",
    format(max(eig_ad), digits = 3), "]\n")
cat("  n_neg(H_AD):", sum(eig_ad < -1e-8), "\n")

if (max_diff >= 1e-5) stop("Gradient validation failed")
cat("\nAll checks passed.\n")
