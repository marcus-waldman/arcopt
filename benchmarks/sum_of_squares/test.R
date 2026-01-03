# Test Sum of Squares Function Derivatives
# ==========================================
#
# Validates analytic gradient and Hessian against finite differences

source("../utils.R")
source("fn.R")
source("gr.R")
source("hess.R")

cat("="[rep(1, 60)], "\n", sep = "")
cat("Testing Sum of Squares Function\n")
cat("="[rep(1, 60)], "\n\n", sep = "")

# Test point 1: Random point in R^5
set.seed(456)
x1 <- rnorm(5)

cat("Test 1: Random point in R^5\n")
cat("x =", sprintf("%.4f", x1), "\n\n")

# Check gradient
cat("Gradient check:\n")
grad_check <- check_gradient(sum_of_squares_fn, sum_of_squares_gr, x1)
cat("  Max difference:", sprintf("%.2e", grad_check$max_diff), "\n")
cat("  Passed:        ", grad_check$passed, "\n\n")

# Check Hessian
cat("Hessian check:\n")
hess_check <- check_hessian(sum_of_squares_gr, sum_of_squares_hess, x1)
cat("  Max difference:", sprintf("%.2e", hess_check$max_diff), "\n")
cat("  Passed:        ", hess_check$passed, "\n\n")

cat("-"[rep(1, 60)], "\n\n", sep = "")

# Test point 2: At optimum
x2 <- rep(0, 5)

cat("Test 2: At optimum (x = 0)\n")
cat("x =", x2, "\n\n")

grad_check2 <- check_gradient(sum_of_squares_fn, sum_of_squares_gr, x2)
cat("Gradient check:\n")
cat("  Max difference:", sprintf("%.2e", grad_check2$max_diff), "\n")
cat("  Passed:        ", grad_check2$passed, "\n\n")

hess_check2 <- check_hessian(sum_of_squares_gr, sum_of_squares_hess, x2)
cat("Hessian check:\n")
cat("  Max difference:", sprintf("%.2e", hess_check2$max_diff), "\n")
cat("  Passed:        ", hess_check2$passed, "\n\n")

cat("-"[rep(1, 60)], "\n\n", sep = "")

# Test point 3: Test ill-conditioning with higher dimension
x3 <- rnorm(20)
cat("Test 3: High dimension (R^20, ill-conditioned)\n")

grad_check3 <- check_gradient(sum_of_squares_fn, sum_of_squares_gr, x3)
cat("Gradient check:\n")
cat("  Max difference:", sprintf("%.2e", grad_check3$max_diff), "\n")
cat("  Passed:        ", grad_check3$passed, "\n\n")

hess_check3 <- check_hessian(sum_of_squares_gr, sum_of_squares_hess, x3)
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
