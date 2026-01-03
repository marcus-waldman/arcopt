# Test Rosenbrock Function Derivatives
# ======================================
#
# Validates analytic gradient and Hessian against finite differences

source("../utils.R")
source("fn.R")
source("gr.R")
source("hess.R")

cat("="[rep(1, 60)], "\n", sep = "")
cat("Testing Rosenbrock Function\n")
cat("="[rep(1, 60)], "\n\n", sep = "")

# Test point 1: Random point in R^5
set.seed(333)
x1 <- rnorm(5)

cat("Test 1: Random point in R^5\n")
cat("x =", sprintf("%.4f", x1), "\n\n")

# Check gradient (use looser tolerance for Rosenbrock due to high curvature)
cat("Gradient check:\n")
grad_check <- check_gradient(rosenbrock_fn, rosenbrock_gr, x1, tol = 1e-4)
cat("  Max difference:", sprintf("%.2e", grad_check$max_diff), "\n")
cat("  Passed:        ", grad_check$passed, "\n\n")

# Check Hessian
cat("Hessian check:\n")
hess_check <- check_hessian(rosenbrock_gr, rosenbrock_hess, x1)
cat("  Max difference:", sprintf("%.2e", hess_check$max_diff), "\n")
cat("  Passed:        ", hess_check$passed, "\n\n")

cat("-"[rep(1, 60)], "\n\n", sep = "")

# Test point 2: At optimum
x2 <- rep(1, 5)

cat("Test 2: At optimum (x = 1, ..., 1)\n")
cat("x =", x2, "\n")
cat("f(x) =", sprintf("%.4f", rosenbrock_fn(x2)), "(should be 0)\n\n")

grad_check2 <- check_gradient(rosenbrock_fn, rosenbrock_gr, x2, tol = 1e-4)
cat("Gradient check:\n")
cat("  Max difference:", sprintf("%.2e", grad_check2$max_diff), "\n")
cat("  Passed:        ", grad_check2$passed, "\n\n")

hess_check2 <- check_hessian(rosenbrock_gr, rosenbrock_hess, x2)
cat("Hessian check:\n")
cat("  Max difference:", sprintf("%.2e", hess_check2$max_diff), "\n")
cat("  Passed:        ", hess_check2$passed, "\n\n")

cat("-"[rep(1, 60)], "\n\n", sep = "")

# Test point 3: Classic starting point (-1.2, 1, ...)
x3 <- c(-1.2, 1, -1.2, 1, -1.2)
cat("Test 3: Classic starting point\n")
cat("x =", sprintf("%.4f", x3), "\n\n")

grad_check3 <- check_gradient(rosenbrock_fn, rosenbrock_gr, x3, tol = 1e-4)
cat("Gradient check:\n")
cat("  Max difference:", sprintf("%.2e", grad_check3$max_diff), "\n")
cat("  Passed:        ", grad_check3$passed, "\n\n")

hess_check3 <- check_hessian(rosenbrock_gr, rosenbrock_hess, x3)
cat("Hessian check:\n")
cat("  Max difference:", sprintf("%.2e", hess_check3$max_diff), "\n")
cat("  Passed:        ", hess_check3$passed, "\n\n")

cat("="[rep(1, 60)], "\n", sep = "")

# Summary
all_passed <- grad_check$passed && grad_check2$passed && grad_check3$passed &&
              hess_check$passed && hess_check2$passed && hess_check3$passed

if (all_passed) {
  cat("✓ All checks passed!\n")
} else {
  cat("✗ Some checks failed!\n")
  stop("Derivative validation failed")
}

cat("="[rep(1, 60)], "\n", sep = "")
