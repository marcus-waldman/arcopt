# Gradient of Trid Function
# ===========================
#
# Gradient:
#   ∂f/∂x_1 = 2*(x_1 - 1) - x_2
#   ∂f/∂x_i = 2*(x_i - 1) - x_{i-1} - x_{i+1}  for 1 < i < d
#   ∂f/∂x_d = 2*(x_d - 1) - x_{d-1}
#
# Derivation:
#   f = sum((x_i - 1)^2) - sum(x_i * x_{i-1})
#   Each x_i appears in:
#     - (x_i - 1)^2
#     - -x_i * x_{i-1} (if i > 1)
#     - -x_i * x_{i+1} (if i < d)

trid_gr <- function(x) {
  d <- length(x)
  g <- numeric(d)

  if (d == 1) {
    # Single variable case
    g[1] <- 2 * (x[1] - 1)
  } else {
    # First element
    g[1] <- 2 * (x[1] - 1) - x[2]

    # Middle elements
    if (d > 2) {
      for (i in 2:(d - 1)) {
        g[i] <- 2 * (x[i] - 1) - x[i - 1] - x[i + 1]
      }
    }

    # Last element
    g[d] <- 2 * (x[d] - 1) - x[d - 1]
  }

  g
}
