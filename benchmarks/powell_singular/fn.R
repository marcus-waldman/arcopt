# Powell Singular Function
# =========================
#
# Functional Form (4 dimensions only):
#   f(x) = (x_1 + 10*x_2)^2 + 5*(x_3 - x_4)^2 +
#          (x_2 - 2*x_3)^4 + 10*(x_1 - x_4)^4
#
# Global Minimum:
#   x* = (0, 0, 0, 0)
#   f* = 0
#
# Properties:
#   - Non-separable
#   - Hessian is singular at the solution
#   - Challenging for optimization methods

powell_singular_fn <- function(x) {
  if (length(x) != 4) {
    stop("Powell Singular function is defined for 4 dimensions only")
  }

  (x[1] + 10 * x[2])^2 + 5 * (x[3] - x[4])^2 +
    (x[2] - 2 * x[3])^4 + 10 * (x[1] - x[4])^4
}
