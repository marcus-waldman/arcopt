# Hessian of Powell Singular Function
# =====================================
#
# Hessian (4x4 matrix):
#   Let u = x_1 - x_4 and v = x_2 - 2*x_3
#
#   H[1,1] = 2 + 120*u^2
#   H[1,2] = H[2,1] = 20
#   H[1,3] = H[3,1] = 0
#   H[1,4] = H[4,1] = -120*u^2
#
#   H[2,2] = 200 + 12*v^2
#   H[2,3] = H[3,2] = -24*v^2
#   H[2,4] = H[4,2] = 0
#
#   H[3,3] = 10 + 48*v^2
#   H[3,4] = H[4,3] = -10
#
#   H[4,4] = 10 + 120*u^2
#
# Note:
#   At the solution (0,0,0,0), u=0 and v=0, giving:
#     H = [  2  20   0   0 ]
#         [ 20 200   0   0 ]
#         [  0   0  10 -10 ]
#         [  0   0 -10  10 ]
#   This matrix is singular (rank 2)

powell_singular_hess <- function(x) {
  if (length(x) != 4) {
    stop("Powell Singular function is defined for 4 dimensions only")
  }

  u <- x[1] - x[4]
  v <- x[2] - 2 * x[3]

  matrix(c(
    2 + 120 * u^2, 20, 0, -120 * u^2,
    20, 200 + 12 * v^2, -24 * v^2, 0,
    0, -24 * v^2, 10 + 48 * v^2, -10,
    -120 * u^2, 0, -10, 10 + 120 * u^2
  ), 4, 4, byrow = TRUE)
}
