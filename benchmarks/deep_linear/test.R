# Validate deep linear network gradient via FD
# ============================================

source("../utils.R")
source("setup.R")

cat("Testing deep linear network gradient (d_in=3, h1=4, h2=4, d_out=3)\n")

prob <- make_deep_linear_problem(d_in = 3L, h1 = 4L, h2 = 4L,
                                 d_out = 3L, n_samples = 6L,
                                 w_star_seed = 11L, data_seed = 13L)
cat(sprintf("  n = %d\n", prob$dims$n))

set.seed(42)
x_test <- stats::rnorm(prob$dims$n, sd = 0.3)

grad_check <- check_gradient(prob$fn, prob$gr, x_test, tol = 1e-6)
cat(sprintf("  gradient max diff: %.2e  passed: %s\n",
            grad_check$max_diff, grad_check$passed))

x_star <- c(
  as.vector(rbind(diag(3), matrix(0, 1, 3))),
  as.vector(rbind(diag(3), matrix(0, 1, 3)) %*% cbind(diag(3), rep(0, 3)) +
              0),
  as.vector(prob$w_star %*% cbind(diag(3), rep(0, 3)))
)

cat("\n  At origin: f(0) =", prob$fn(rep(0, prob$dims$n)),
    " (should equal 0.5 * ||Y||_F^2 =", 0.5 * sum(prob$y_data^2), ")\n")

g_origin <- prob$gr(rep(0, prob$dims$n))
cat(sprintf("  ||grad(0)||_2 = %.2e   (origin is a critical point)\n",
            sqrt(sum(g_origin^2))))

if (!grad_check$passed) stop("Gradient validation failed")
cat("\nAll checks passed.\n")
