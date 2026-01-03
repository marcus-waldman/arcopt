# Gradient of Bohachevsky Function
# ==================================
#
# Gradient (2 dimensions):
#   ∂f/∂x1 = 2*x1 + 0.9*π*sin(3*π*x1)
#   ∂f/∂x2 = 4*x2 + 1.6*π*sin(4*π*x2)
#
# Derivation:
#   ∂/∂x1 [x1^2 - 0.3*cos(3π*x1)] = 2*x1 + 0.3*3π*sin(3π*x1)
#   ∂/∂x2 [2*x2^2 - 0.4*cos(4π*x2)] = 4*x2 + 0.4*4π*sin(4π*x2)

bohachevsky_gr <- function(x) {
  if (length(x) != 2) {
    stop("Bohachevsky function is defined for 2 dimensions only")
  }

  c(
    2 * x[1] + 0.9 * pi * sin(3 * pi * x[1]),
    4 * x[2] + 1.6 * pi * sin(4 * pi * x[2])
  )
}
