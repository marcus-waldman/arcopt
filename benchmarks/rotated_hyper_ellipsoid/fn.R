# Rotated Hyper-Ellipsoid Function
# ==================================
#
# Functional Form:
#   f(x) = sum_{j=1}^d (d - j + 1) * x_j^2
#
# Expanded:
#   f(x) = d*x_1^2 + (d-1)*x_2^2 + ... + 1*x_d^2
#
# Global Minimum:
#   x* = (0, 0, ..., 0)
#   f* = 0
#
# Properties:
#   - Convex
#   - Separable
#   - Unimodal
#   - Highly ill-conditioned (condition number = d)

rotated_hyper_ellipsoid_fn <- function(x) {
  n <- length(x)
  sum((n:1) * x^2)
}
