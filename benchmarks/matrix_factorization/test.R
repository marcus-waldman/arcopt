# Validate matrix factorization derivatives
# =========================================
#
# Uses central finite differences on a small instance (d=6, r=2 -> n=12)
# to confirm gradient and Hessian are correct.

source("../utils.R")
source("setup.R")

cat("Testing matrix factorization derivatives (d=6, r=2, n=12)\n")

prob <- make_matrix_factorization_problem(d = 6L, r = 2L, u_true_seed = 7L)
set.seed(42)
x_test <- stats::rnorm(prob$n, sd = 0.5)

grad_check <- check_gradient(prob$fn, prob$gr, x_test, tol = 1e-6)
cat(sprintf("  gradient max diff: %.2e  passed: %s\n",
            grad_check$max_diff, grad_check$passed))

hess_check <- check_hessian(prob$gr, prob$hess, x_test, tol = 1e-4)
cat(sprintf("  hessian  max diff: %.2e  passed: %s\n",
            hess_check$max_diff, hess_check$passed))

# At a global minimum U = u_true the gradient should be ~0 and the
# residual should be ~0.
x_min <- as.vector(prob$u_true)
cat(sprintf("\n  at u_true: f=%.2e  ||grad||=%.2e\n",
            prob$fn(x_min), sqrt(sum(prob$gr(x_min)^2))))

# At origin: gradient 0, Hessian indefinite (symmetric saddle).
x_origin <- rep(0, prob$n)
cat(sprintf("  at origin: f=%.2e  ||grad||=%.2e\n",
            prob$fn(x_origin), sqrt(sum(prob$gr(x_origin)^2))))
eig_origin <- eigen(prob$hess(x_origin), symmetric = TRUE, only.values = TRUE)$values
cat(sprintf("  origin hessian eigenvalue range: [%.2f, %.2f]  (n_neg=%d)\n",
            min(eig_origin), max(eig_origin), sum(eig_origin < -1e-8)))

all_passed <- grad_check$passed && hess_check$passed
if (!all_passed) {
  stop("Derivative validation failed")
}
cat("\nAll checks passed.\n")
