# Validate Brown Almost-Linear derivatives
# ========================================

source("../utils.R")
source("setup.R")

cat("Testing Brown Almost-Linear derivatives (n = 5)\n")

prob <- make_brown_problem(n = 5L)

set.seed(17)
x_test <- stats::runif(prob$n, min = 0.3, max = 1.5)

grad_check <- check_gradient(prob$fn, prob$gr, x_test, tol = 1e-6)
cat(sprintf("  gradient max diff: %.2e  passed: %s\n",
            grad_check$max_diff, grad_check$passed))

hess_check <- check_hessian(prob$gr, prob$hess, x_test, tol = 1e-4)
cat(sprintf("  hessian  max diff: %.2e  passed: %s\n",
            hess_check$max_diff, hess_check$passed))

# At global min x* = (1,...,1)
x_star <- rep(1, prob$n)
cat(sprintf("\n  at x* = (1,...,1): F = %.2e  ||grad|| = %.2e\n",
            prob$fn(x_star), sqrt(sum(prob$gr(x_star)^2))))
eig_star <- eigen(prob$hess(x_star), symmetric = TRUE,
                  only.values = TRUE)$values
cat(sprintf("  hessian at x*: eigenvalues = [%.2f, %.2f]  (n_neg = %d)\n",
            min(eig_star), max(eig_star), sum(eig_star < -1e-8)))

# At standard start x_0 = (0.5,...,0.5)
x0 <- rep(0.5, prob$n)
cat(sprintf("\n  at x_0 = (0.5,...): F = %.2f  ||grad|| = %.2e\n",
            prob$fn(x0), sqrt(sum(prob$gr(x0)^2))))

if (!grad_check$passed || !hess_check$passed) {
  stop("Derivative validation failed")
}
cat("\nAll checks passed.\n")
