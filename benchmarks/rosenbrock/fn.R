# Rosenbrock Function
# ====================
#
# Functional Form:
#   f(x) = sum_{i=1}^{d-1} [100*(x_{i+1} - x_i^2)^2 + (1 - x_i)^2]
#
# Global Minimum:
#   x* = (1, 1, ..., 1)
#   f* = 0
#
# Properties:
#   - Non-convex for d > 2
#   - Non-separable
#   - Unimodal
#   - Famous narrow curved valley

rosenbrock_fn <- function(x) {
  d <- length(x)

  if (d < 2) {
    stop("Rosenbrock function requires at least 2 dimensions")
  }

  # Sum over pairs of adjacent variables
  sum_val <- 0
  for (i in 1:(d - 1)) {
    sum_val <- sum_val + 100 * (x[i + 1] - x[i]^2)^2 + (1 - x[i])^2
  }

  sum_val
}
