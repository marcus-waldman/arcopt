# Gradient of Rotated Hyper-Ellipsoid Function
# ==============================================
#
# Gradient:
#   ∇f(x) = [2*d*x_1, 2*(d-1)*x_2, ..., 2*1*x_d]^T
#         = 2 * [d:1] .* x
#
# Derivation:
#   ∂f/∂x_j = ∂/∂x_j ((d - j + 1) * x_j^2) = 2*(d - j + 1)*x_j

rotated_hyper_ellipsoid_gr <- function(x) {
  n <- length(x)
  2 * (n:1) * x
}
