# Gradient of Sum of Squares Function
# =====================================
#
# Gradient:
#   ∇f(x) = [2*1*x_1, 2*2*x_2, ..., 2*d*x_d]^T
#         = 2 * [1:d] .* x
#
# Derivation:
#   ∂f/∂x_i = ∂/∂x_i (i * x_i^2) = 2*i*x_i

sum_of_squares_gr <- function(x) {
  n <- length(x)
  2 * (1:n) * x
}
