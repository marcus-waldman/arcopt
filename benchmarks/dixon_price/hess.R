# Hessian of Dixon-Price Function
# =================================
#
# Hessian:
#   Banded structure with bandwidth 2 (depends on x)
#
#   Diagonal elements:
#     H_{1,1} = 2 + 4 = 6
#
#     H_{i,i} = 48*i*x_i^2 - 8*i*x_{i-1} + 2*(i+1)  for 1 < i < d
#
#     H_{d,d} = 48*d*x_d^2 - 8*d*x_{d-1}
#
#   Off-diagonal elements:
#     H_{i-1,i} = H_{i,i-1} = -8*i*x_i  for i = 2, ..., d
#     H_{i,i+1} = H_{i+1,i} = -8*(i+1)*x_{i+1}  for i = 1, ..., d-1

dixon_price_hess <- function(x) {
  d <- length(x)

  if (d < 2) {
    stop("Dixon-Price function requires at least 2 dimensions")
  }

  # Initialize Hessian
  H <- matrix(0, d, d)

  # First diagonal element
  H[1, 1] <- 6

  # Middle diagonal elements
  if (d > 2) {
    for (i in 2:(d - 1)) {
      H[i, i] <- 48 * i * x[i]^2 - 8 * i * x[i - 1] + 2 * (i + 1)
    }
  }

  # Last diagonal element
  H[d, d] <- 48 * d * x[d]^2 - 8 * d * x[d - 1]

  # Off-diagonal elements
  for (i in 2:d) {
    # H[i-1, i] from derivative of (i-1) row with respect to x_i
    H[i - 1, i] <- -8 * i * x[i]
    H[i, i - 1] <- -8 * i * x[i]
  }

  H
}
