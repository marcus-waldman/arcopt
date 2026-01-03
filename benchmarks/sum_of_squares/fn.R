# Sum of Squares Function
# ========================
#
# Functional Form:
#   f(x) = sum(i * x_i^2) for i = 1, ..., d
#
# Global Minimum:
#   x* = (0, 0, ..., 0)
#   f* = 0
#
# Properties:
#   - Convex
#   - Separable
#   - Unimodal
#   - Diagonal Hessian with increasing conditioning

sum_of_squares_fn <- function(x) {
  n <- length(x)
  sum((1:n) * x^2)
}
