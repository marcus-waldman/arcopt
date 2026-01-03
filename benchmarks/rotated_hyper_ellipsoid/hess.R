# Hessian of Rotated Hyper-Ellipsoid Function
# =============================================
#
# Hessian:
#   H = diag(2*d, 2*(d-1), ..., 2*1)
#     = diag(2 * [d:1])
#
# Derivation:
#   ∂²f/∂x_j∂x_k = 2*(d - j + 1) if j == k, 0 otherwise
#
# Note:
#   Constant diagonal Hessian - does not depend on x
#   Highly ill-conditioned: cond(H) ≈ d
#   First variable has d times larger curvature than last

rotated_hyper_ellipsoid_hess <- function(x) {
  n <- length(x)
  diag(2 * (n:1))
}
