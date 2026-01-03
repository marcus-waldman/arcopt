# Gradient of Powell Singular Function
# ======================================
#
# Gradient (4 dimensions):
#   ∂f/∂x_1 = 2*(x_1 + 10*x_2) + 40*(x_1 - x_4)^3
#   ∂f/∂x_2 = 20*(x_1 + 10*x_2) + 4*(x_2 - 2*x_3)^3
#   ∂f/∂x_3 = 10*(x_3 - x_4) - 8*(x_2 - 2*x_3)^3
#   ∂f/∂x_4 = -10*(x_3 - x_4) - 40*(x_1 - x_4)^3
#
# Derivation:
#   From (x_1 + 10*x_2)^2: 2*(x_1 + 10*x_2)*{1, 10, 0, 0}
#   From 5*(x_3 - x_4)^2: 10*(x_3 - x_4)*{0, 0, 1, -1}
#   From (x_2 - 2*x_3)^4: 4*(x_2 - 2*x_3)^3*{0, 1, -2, 0}
#   From 10*(x_1 - x_4)^4: 40*(x_1 - x_4)^3*{1, 0, 0, -1}

powell_singular_gr <- function(x) {
  if (length(x) != 4) {
    stop("Powell Singular function is defined for 4 dimensions only")
  }

  c(
    2 * (x[1] + 10 * x[2]) + 40 * (x[1] - x[4])^3,
    20 * (x[1] + 10 * x[2]) + 4 * (x[2] - 2 * x[3])^3,
    10 * (x[3] - x[4]) - 8 * (x[2] - 2 * x[3])^3,
    -10 * (x[3] - x[4]) - 40 * (x[1] - x[4])^3
  )
}
