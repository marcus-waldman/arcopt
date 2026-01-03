# arcopt: Adaptive Regularization using Cubics for R

<!-- badges: start -->
[![R-CMD-check](https://github.com/marcus-waldman/arcopt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/marcus-waldman/arcopt/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/marcus-waldman/arcopt/branch/main/graph/badge.svg)](https://app.codecov.io/gh/marcus-waldman/arcopt)
[![lint](https://github.com/marcus-waldman/arcopt/actions/workflows/lint.yaml/badge.svg)](https://github.com/marcus-waldman/arcopt/actions/workflows/lint.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/arcopt)](https://CRAN.R-project.org/package=arcopt)
<!-- badges: end -->

## Overview

**arcopt** implements cubic regularization methods (ARC) for robust local optimization in statistical applications. Designed for statisticians and applied researchers working with challenging optimization landscapes common in maximum likelihood estimation, posterior mode finding, and nonlinear regression.

### Key Features

- **Hessian-Centric**: Leverages second-order information for faster convergence on ill-conditioned problems
- **Robust by Default**: Handles indefinite Hessians and saddle points automatically without manual tuning
- **Constraint Support**: Box constraints and linear equality constraints built-in
- **Diagnostic Rich**: Comprehensive convergence diagnostics and algorithm insights

## Installation

You can install the development version of arcopt from [GitHub](https://github.com/marcus-waldman/arcopt) with:

```r
# install.packages("pak")
pak::pak("marcus-waldman/arcopt")
```

## Quick Start

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

result$par      # Optimal parameters
result$value    # Optimal function value
result$converged  # Convergence status
```

## Why arcopt?

Traditional optimization methods struggle with:
- **Ill-conditioned Hessians**: Poor scaling, high condition numbers
- **Indefinite curvature**: Negative eigenvalues, saddle points
- **Nonconvex landscapes**: Local minima, ridges

**arcopt** handles these challenges through:
1. **Cubic regularization**: Adaptively regularizes indefinite Hessians
2. **Modified Cholesky**: Ensures numerical stability even with negative curvature
3. **Intelligent safeguards**: Automatic stagnation detection, NaN/Inf handling

## Documentation

- [Function Reference](https://marcus-waldman.github.io/arcopt/reference/): Complete API documentation
- [Introduction Vignette](https://marcus-waldman.github.io/arcopt/articles/arcopt-introduction.html): Getting started guide
- [Algorithm Details](https://github.com/marcus-waldman/arcopt/blob/main/design/pseudocode.qmd): Implementation specifications
- [Design Principles](https://github.com/marcus-waldman/arcopt/blob/main/design/design-principles.qmd): Design philosophy

## Comparison with Other Methods

| Method | Handles Indefinite H | Escapes Saddles | Requires Tuning |
|--------|---------------------|-----------------|-----------------|
| Newton | ❌ (fails) | ❌ | N/A |
| Trust-region Newton | ⚠️ (with modifications) | ⚠️ | ✅ (radius) |
| BFGS/L-BFGS | ⚠️ (maintains PD approx) | ❌ | ✅ (initial Hessian) |
| **ARC (arcopt)** | ✅ (native) | ✅ (automatic) | ❌ (adaptive σ) |

## Citation

If you use arcopt in your research, please cite:

```r
citation("arcopt")
```

And consider citing the foundational ARC paper:

> Cartis, C., Gould, N. I. M., & Toint, P. L. (2011). Adaptive cubic regularisation methods for unconstrained optimization. Part I: Motivation, convergence and numerical results. *Mathematical Programming*, 127(2), 245-295. https://doi.org/10.1007/s10107-009-0286-5

## Contributing

Contributions are welcome! Please see the [development guidelines](https://marcus-waldman.github.io/arcopt/articles/contributing.html) for details.

## License

MIT © 2026 Marcus Waldman
