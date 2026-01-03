# arcopt

Adaptive Regularization using Cubics (ARC) optimizer for R.

Handles indefinite Hessians, ill-conditioned problems, and saddle points automatically. Designed for robust optimization in statistical applications (MLE, posterior modes, nonlinear regression).

## Installation

```r
pak::pak("marcus-waldman/arcopt")
```

## Usage

```r
library(arcopt)

# Rosenbrock function (classic test problem)
rosenbrock <- function(x) {
  (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
}

rosenbrock_gr <- function(x) {
  c(
    -2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
    200 * (x[2] - x[1]^2)
  )
}

rosenbrock_hess <- function(x) {
  matrix(c(
    1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
    -400 * x[1], 200
  ), 2, 2)
}

# Optimize
result <- arcopt(
  x0 = c(-1.2, 1),
  fn = rosenbrock,
  gr = rosenbrock_gr,
  hess = rosenbrock_hess
)

result$par       # Optimal parameters
result$value     # Optimal function value
result$converged # Convergence status
```

## Features

- **Robust**: Handles indefinite and ill-conditioned Hessians via cubic regularization
- **Box constraints**: Optional parameter bounds
- **Safeguards**: Automatic stagnation detection and NaN/Inf protection
- **Newton-first**: Tries pure Newton when Hessian is positive definite
- **Adaptive**: Two sigma update strategies (classical CGT and interpolation)

## Requirements

- R >= 4.0
- Analytic Hessian function required

## Reference

Cartis, C., Gould, N. I. M., & Toint, P. L. (2011). Adaptive cubic regularisation methods for unconstrained optimization. *Mathematical Programming*, 127(2), 245-295.

## License

MIT © 2026 Marcus Waldman
