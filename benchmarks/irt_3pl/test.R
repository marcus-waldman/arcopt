# Validate 3PL MML fn and gr via finite differences
# =================================================

source("../utils.R")
source("setup.R")

data <- simulate_irt3pl(j_items = 10L, n_examinees = 1000L,
                        items = "hard", seed = 42L)
fns <- make_irt3pl_fns(data, q_nodes = 40L)

cat(sprintf("3PL (MML): J=%d items (hard), N=%d examinees, Q=40 nodes, n=%d params\n",
            data$j_items, data$n_examinees, length(data$x_true)))

set.seed(7)
x_test <- data$x_true + stats::rnorm(length(data$x_true), sd = 0.2)

cat(sprintf("\n  f(x_true)      = %.4f\n", fns$fn(data$x_true)))
cat(sprintf("  f(x_test)      = %.4f\n", fns$fn(x_test)))
cat(sprintf("  ||grad(x_true)|| = %.3e  (not zero due to finite N)\n",
            sqrt(sum(fns$gr(data$x_true)^2))))
cat(sprintf("  ||grad(x_test)|| = %.3e\n\n",
            sqrt(sum(fns$gr(x_test)^2))))

cat("Gradient check at x_test:\n")
grad_check <- check_gradient(fns$fn, fns$gr, x_test, tol = 1e-4, h = 1e-6)
cat(sprintf("  max diff: %.2e  passed: %s\n",
            grad_check$max_diff, grad_check$passed))

if (!grad_check$passed) stop("Gradient validation failed")
cat("\nGradient OK.\n")
