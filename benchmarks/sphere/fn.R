# Sphere Function
# ================
#
# Functional Form:
#   f(x) = sum(x_i^2) for i = 1, ..., d
#
# Global Minimum:
#   x* = (0, 0, ..., 0)
#   f* = 0
#
# Properties:
#   - Convex
#   - Separable
#   - Unimodal
#   - Constant Hessian (H = 2I)

sphere_fn <- function(x) {
  sum(x^2)
}
