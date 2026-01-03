# Hessian of Bohachevsky Function
# =================================
#
# Hessian (2x2 matrix):
#   H11 = 2 + 2.7*π^2*cos(3*π*x1)
#   H22 = 4 + 6.4*π^2*cos(4*π*x2)
#   H12 = H21 = 0 (separable function)
#
# Derivation:
#   ∂²f/∂x1² = 2 + 0.9*π*3*π*cos(3*π*x1) = 2 + 2.7*π^2*cos(3*π*x1)
#   ∂²f/∂x2² = 4 + 1.6*π*4*π*cos(4*π*x2) = 4 + 6.4*π^2*cos(4*π*x2)
#   ∂²f/∂x1∂x2 = 0 (no cross terms)
#
# Note:
#   Diagonal Hessian (separable function)
#   Hessian depends on x due to cosine terms

bohachevsky_hess <- function(x) {
  if (length(x) != 2) {
    stop("Bohachevsky function is defined for 2 dimensions only")
  }

  matrix(c(
    2 + 2.7 * pi^2 * cos(3 * pi * x[1]), 0,
    0, 4 + 6.4 * pi^2 * cos(4 * pi * x[2])
  ), 2, 2)
}
