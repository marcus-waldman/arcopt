# Saddle Point Test Function for Hex Sticker
# Double-well potential with saddle at origin

#' Saddle point test function
#'
#' f(x, y) = x^4/4 - x^2/2 + y^2/2
#'
#' Properties:
#' - Saddle point at (0, 0) with Hessian diag(-1, 1) (indefinite)
#' - Two minima at (-1, 0) and (1, 0)
#' - Contour structure: two basins separated by saddle ridge
#'
#' @param x Numeric vector of length 2
#' @return Function value (scalar)
saddle_fn <- function(x) {
  x[1]^4 / 4 - x[1]^2 / 2 + x[2]^2 / 2
}

#' Gradient of saddle point function
#'
#' @param x Numeric vector of length 2
#' @return Gradient vector of length 2
saddle_gr <- function(x) {
  c(x[1]^3 - x[1], x[2])
}

#' Hessian of saddle point function
#'
#' @param x Numeric vector of length 2
#' @return 2x2 Hessian matrix
saddle_hess <- function(x) {
  matrix(c(3 * x[1]^2 - 1, 0,
           0, 1),
         nrow = 2, ncol = 2, byrow = TRUE)
}

#' Vectorized version for grid evaluation
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @return Function values
saddle_fn_vectorized <- function(x, y) {
  x^4 / 4 - x^2 / 2 + y^2 / 2
}
