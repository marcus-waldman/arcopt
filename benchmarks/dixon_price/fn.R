# Dixon-Price Function
# =====================
#
# Functional Form:
#   f(x) = (x_1 - 1)^2 + sum_{i=2}^d i*(2*x_i^2 - x_{i-1})^2
#
# Global Minimum:
#   x_i^* = 2^{-(2^i - 2)/2^i}
#   f^* = 0
#
# Properties:
#   - Unimodal
#   - Non-separable
#   - Narrow valley

dixon_price_fn <- function(x) {
  d <- length(x)

  if (d < 2) {
    stop("Dixon-Price function requires at least 2 dimensions")
  }

  # First term
  term1 <- (x[1] - 1)^2

  # Sum term
  sum_val <- 0
  for (i in 2:d) {
    sum_val <- sum_val + i * (2 * x[i]^2 - x[i - 1])^2
  }

  term1 + sum_val
}
