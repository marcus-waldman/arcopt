# Hessian of Sphere Function
# ===========================
#
# Hessian:
#   H = 2I_d
#   (2 times the d×d identity matrix)
#
# Derivation:
#   ∂²f/∂x_i∂x_j = 2 if i == j, 0 otherwise
#
# Note:
#   Constant Hessian - does not depend on x

sphere_hess <- function(x) {
  n <- length(x)
  2 * diag(n)
}
