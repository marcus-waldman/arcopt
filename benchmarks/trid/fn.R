# Trid Function
# ==============
#
# Functional Form:
#   f(x) = sum_{i=1}^d (x_i - 1)^2 - sum_{i=2}^d x_i * x_{i-1}
#
# Global Minimum:
#   x_i^* = i*(d + 1 - i)
#   f^* = -d*(d + 4)*(d - 1)/6
#
# For d=5: x* = (5, 8, 9, 8, 5), f* = -50
#
# Properties:
#   - Convex
#   - Non-separable (coupling between adjacent variables)
#   - Unimodal
#   - Origin is a saddle point

trid_fn <- function(x) {
  d <- length(x)

  # First term: sum of squared deviations from 1
  term1 <- sum((x - 1)^2)

  # Second term: negative sum of products of adjacent variables
  if (d > 1) {
    term2 <- sum(x[2:d] * x[1:(d - 1)])
  } else {
    term2 <- 0
  }

  term1 - term2
}
