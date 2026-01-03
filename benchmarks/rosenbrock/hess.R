# Hessian of Rosenbrock Function
# ================================
#
# Hessian:
#   Band-diagonal structure (tridiagonal-like, but depends on x)
#
#   Let r_i = x_{i+1} - x_i^2
#
#   Diagonal elements:
#     H_{1,1} = -400*(x_2 - x_1^2) + 800*x_1^2 + 2
#             = -400*r_1 + 800*x_1^2 + 2
#
#     H_{i,i} = 200 + 1200*x_i^2 - 400*x_{i+1} + 2 (for 1 < i < d)
#             = 200 - 400*r_i + 800*x_i^2 + 2
#
#     H_{d,d} = 200
#
#   Off-diagonal elements (only adjacent coupling):
#     H_{i,i+1} = H_{i+1,i} = -400*x_i  (for i = 1, ..., d-1)
#
#   All other entries are zero

rosenbrock_hess <- function(x) {
  d <- length(x)

  if (d < 2) {
    stop("Rosenbrock function requires at least 2 dimensions")
  }

  # Initialize Hessian
  H <- matrix(0, d, d)

  # First diagonal element
  H[1, 1] <- -400 * (x[2] - x[1]^2) + 800 * x[1]^2 + 2

  # Middle diagonal elements
  if (d > 2) {
    for (i in 2:(d - 1)) {
      H[i, i] <- 200 - 400 * (x[i + 1] - x[i]^2) + 800 * x[i]^2 + 2
    }
  }

  # Last diagonal element
  H[d, d] <- 200

  # Off-diagonal elements
  for (i in 1:(d - 1)) {
    H[i, i + 1] <- -400 * x[i]
    H[i + 1, i] <- -400 * x[i]
  }

  H
}
