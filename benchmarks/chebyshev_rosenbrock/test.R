# Validate Chebyshev-Rosenbrock derivatives
# =========================================

source("../utils.R")
source("setup.R")

cat("Testing Chebyshev-Rosenbrock derivatives (n = 8)\n")

prob <- make_chebrosen_problem(n = 8L)

set.seed(123)
x_test <- stats::rnorm(prob$n, sd = 0.7)

grad_check <- check_gradient(prob$fn, prob$gr, x_test, tol = 1e-6)
cat(sprintf("  gradient max diff: %.2e  passed: %s\n",
            grad_check$max_diff, grad_check$passed))

hess_check <- check_hessian(prob$gr, prob$hess, x_test, tol = 1e-4)
cat(sprintf("  hessian  max diff: %.2e  passed: %s\n",
            hess_check$max_diff, hess_check$passed))

# At global min
x_star <- rep(1, prob$n)
cat(sprintf("\n  at x* = (1,...,1): f = %.2e  ||grad|| = %.2e\n",
            prob$fn(x_star), sqrt(sum(prob$gr(x_star)^2))))

# At adversarial start
x_adv <- rep(-1, prob$n)
cat(sprintf("  at x0 = (-1,...,-1): f = %.2f (expect %d)  ||grad|| = %.2e\n",
            prob$fn(x_adv), 4 * prob$n - 3, sqrt(sum(prob$gr(x_adv)^2))))

eig_adv <- eigen(prob$hess(x_adv), symmetric = TRUE, only.values = TRUE)$values
cat(sprintf("  hessian eigenvalue range at x0: [%.2f, %.2f]  (n_neg = %d)\n",
            min(eig_adv), max(eig_adv), sum(eig_adv < -1e-8)))

if (!grad_check$passed || !hess_check$passed) {
  stop("Derivative validation failed")
}
cat("\nAll checks passed.\n")
