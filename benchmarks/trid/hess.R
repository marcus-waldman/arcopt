# Hessian of Trid Function
# ==========================
#
# Hessian:
#   Tridiagonal matrix:
#     H[i,i] = 2 for all i
#     H[i,i+1] = H[i+1,i] = -1 for i = 1, ..., d-1
#     All other entries = 0
#
# Matrix form:
#   H = [ 2  -1   0   0  ...  0 ]
#       [-1   2  -1   0  ...  0 ]
#       [ 0  -1   2  -1  ...  0 ]
#       [ :   :   :   :   :   : ]
#       [ 0   0   0   0  -1   2 ]
#
# Derivation:
#   ∂²f/∂x_i² = 2 for all i
#   ∂²f/∂x_i∂x_{i+1} = -1 (from -x_i*x_{i+1} term)
#   All other second derivatives = 0
#
# Note:
#   Constant Hessian - does not depend on x
#   Symmetric positive definite

trid_hess <- function(x) {
  d <- length(x)

  # Initialize matrix
  H <- matrix(0, d, d)

  # Fill diagonal with 2
  diag(H) <- 2

  # Fill off-diagonals with -1
  if (d > 1) {
    for (i in 1:(d - 1)) {
      H[i, i + 1] <- -1
      H[i + 1, i] <- -1
    }
  }

  H
}
