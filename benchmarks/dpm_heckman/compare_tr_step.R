# Compare one step of my solve_tr_eigen against trust::trust's first step
# on the DPM Heckman midtrajectory iterate.

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})
source(file.path("benchmarks", "dpm_heckman", "setup.R"))

data <- simulate_dpm_heckman(n_obs = 500L, k_true = 2L,
                             p_out = 5L, p_sel = 7L, n_excl = 2L,
                             seed = 101L)
model_so <- compile_dpm_heckman_model()
fns <- make_dpm_heckman_fns(data, model_so, k_trunc = 20L)

# Cubic-only for 65 iters to get the midtrajectory iterate
cat("Getting midtrajectory iterate (arcopt cubic, maxit=65)...\n")
fit_cubic <- arcopt(fns$x0, fns$fn, fns$gr, fns$hess_ad,
                    control = list(maxit = 65L, gtol_abs = 1e-5, trace = 0,
                                   tr_fallback_enabled = FALSE))
x_mid <- fit_cubic$par
g_mid <- fns$gr(x_mid)
H_mid <- fns$hess_ad(x_mid)
f_mid <- fns$fn(x_mid)

cat(sprintf("  f=%.3f  |g|=%.3f  lambda_min=%.3f  lambda_max=%.3f\n",
            f_mid, max(abs(g_mid)),
            min(eigen(H_mid, symmetric=TRUE, only.values=TRUE)$values),
            max(eigen(H_mid, symmetric=TRUE, only.values=TRUE)$values)))

# My TR step
my_result <- solve_tr_eigen(g_mid, H_mid, radius = 1.0)

# Trust's TR step (force 1 iteration)
objfun <- function(p) list(value = fns$fn(p), gradient = fns$gr(p),
                           hessian = fns$hess_ad(p))
# iterlim=1 to force only one step
t_result <- trust::trust(objfun, x_mid, rinit = 1, rmax = 1,
                         iterlim = 1, fterm = 0, mterm = 0)
# trust's step: x_new - x_mid
trust_step <- t_result$argument - x_mid

cat("\n=== My solve_tr_eigen (radius=1) ===\n")
cat(sprintf("  case_type: %s\n", my_result$case_type))
cat(sprintf("  on_boundary: %s\n", my_result$on_boundary))
cat(sprintf("  lambda: %.3e\n", my_result$lambda))
cat(sprintf("  ||s||: %.6f\n", sqrt(sum(my_result$s^2))))
cat(sprintf("  pred_reduction: %.6e\n", my_result$pred_reduction))

# Evaluate actual reduction at my step
x_trial_mine <- x_mid + my_result$s
f_trial_mine <- fns$fn(x_trial_mine)
cat(sprintf("  actual_reduction: %.6e\n", f_mid - f_trial_mine))
cat(sprintf("  rho: %.4f\n", (f_mid - f_trial_mine) / my_result$pred_reduction))

cat("\n=== trust::trust (rinit=1, rmax=1, iterlim=1) ===\n")
cat(sprintf("  ||step||: %.6f\n", sqrt(sum(trust_step^2))))
cat(sprintf("  value at step: %.3f\n", t_result$value))
cat(sprintf("  trust rho tracking: accepted=%s\n", t_result$accept))

cat("\n=== Similarity ===\n")
my_step_norm <- sqrt(sum(my_result$s^2))
trust_step_norm <- sqrt(sum(trust_step^2))
if (my_step_norm > 0 && trust_step_norm > 0) {
  cos_sim <- sum(my_result$s * trust_step) / (my_step_norm * trust_step_norm)
  cat(sprintf("  step norm ratio (mine / trust): %.3f\n",
              my_step_norm / trust_step_norm))
  cat(sprintf("  step direction cosine: %.6f\n", cos_sim))
}
cat(sprintf("  ||my_step - trust_step|| / ||trust_step||: %.3e\n",
            sqrt(sum((my_result$s - trust_step)^2)) / trust_step_norm))
