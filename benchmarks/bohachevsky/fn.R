# Bohachevsky Function (Function 1)
# ===================================
#
# Functional Form (2 dimensions):
#   f(x1, x2) = x1^2 + 2*x2^2 - 0.3*cos(3π*x1) - 0.4*cos(4π*x2) + 0.7
#
# Global Minimum:
#   x* = (0, 0)
#   f* = 0
#
# Properties:
#   - Multimodal (due to cosine terms)
#   - Separable
#   - 2 dimensions only

bohachevsky_fn <- function(x) {
  if (length(x) != 2) {
    stop("Bohachevsky function is defined for 2 dimensions only")
  }

  x[1]^2 + 2 * x[2]^2 - 0.3 * cos(3 * pi * x[1]) - 0.4 * cos(4 * pi * x[2]) + 0.7
}
