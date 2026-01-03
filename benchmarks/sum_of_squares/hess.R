# Hessian of Sum of Squares Function
# ====================================
#
# Hessian:
#   H = diag(2, 4, 6, ..., 2*d)
#     = diag(2 * [1:d])
#
# Derivation:
#   ∂²f/∂x_i∂x_j = 2*i if i == j, 0 otherwise
#
# Note:
#   Constant diagonal Hessian - does not depend on x
#   Becomes increasingly ill-conditioned as dimension grows
#   Condition number = d

sum_of_squares_hess <- function(x) {
  n <- length(x)
  diag(2 * (1:n))
}
