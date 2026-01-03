# Gradient of Rosenbrock Function
# =================================
#
# Gradient:
#   ∂f/∂x_1 = -400*x_1*(x_2 - x_1^2) - 2*(1 - x_1)
#   ∂f/∂x_i = 200*(x_i - x_{i-1}^2) - 400*x_i*(x_{i+1} - x_i^2) - 2*(1 - x_i)
#             for 1 < i < d
#   ∂f/∂x_d = 200*(x_d - x_{d-1}^2)
#
# Derivation:
#   Each x_i appears in:
#     - 100*(x_{i+1} - x_i^2)^2 (if i < d)
#     - 100*(x_i - x_{i-1}^2)^2 (if i > 1)
#     - (1 - x_i)^2 (if i < d)

rosenbrock_gr <- function(x) {
  d <- length(x)

  if (d < 2) {
    stop("Rosenbrock function requires at least 2 dimensions")
  }

  g <- numeric(d)

  # First element
  g[1] <- -400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1])

  # Middle elements
  if (d > 2) {
    for (i in 2:(d - 1)) {
      g[i] <- 200 * (x[i] - x[i - 1]^2) - 400 * x[i] * (x[i + 1] - x[i]^2) - 2 * (1 - x[i])
    }
  }

  # Last element
  g[d] <- 200 * (x[d] - x[d - 1]^2)

  g
}
