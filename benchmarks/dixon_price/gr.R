# Gradient of Dixon-Price Function
# ==================================
#
# Gradient (let g_i = 2*x_i^2 - x_{i-1}):
#   ∂f/∂x_1 = 2*(x_1 - 1) - 2*2*g_2
#           = 2*(x_1 - 1) - 4*(2*x_2^2 - x_1)
#
#   ∂f/∂x_i = 8*i*x_i*g_i - 2*(i+1)*g_{i+1}  for 1 < i < d
#           = 8*i*x_i*(2*x_i^2 - x_{i-1}) - 2*(i+1)*(2*x_{i+1}^2 - x_i)
#
#   ∂f/∂x_d = 8*d*x_d*g_d
#           = 8*d*x_d*(2*x_d^2 - x_{d-1})

dixon_price_gr <- function(x) {
  d <- length(x)

  if (d < 2) {
    stop("Dixon-Price function requires at least 2 dimensions")
  }

  g <- numeric(d)

  # First element
  g[1] <- 2 * (x[1] - 1) - 4 * (2 * x[2]^2 - x[1])

  # Middle elements
  if (d > 2) {
    for (i in 2:(d - 1)) {
      g_i <- 2 * x[i]^2 - x[i - 1]
      g_i_plus_1 <- 2 * x[i + 1]^2 - x[i]
      g[i] <- 8 * i * x[i] * g_i - 2 * (i + 1) * g_i_plus_1
    }
  }

  # Last element
  g_d <- 2 * x[d]^2 - x[d - 1]
  g[d] <- 8 * d * x[d] * g_d

  g
}
