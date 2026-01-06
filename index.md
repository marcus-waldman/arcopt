# arcopt ![](reference/figures/logo.png)

> **Robust nonlinear optimization for R**

When [`optim()`](https://rdrr.io/r/stats/optim.html) fails on your
ill-conditioned model, **arcopt** is designed to succeed. It uses
Adaptive Regularization with Cubics (ARC) to handle indefinite Hessians
and escape saddle points automatically.

## Installation

``` r
# Install from GitHub
pak::pak("marcus-waldman/arcopt")

# Or with devtools
devtools::install_github("marcus-waldman/arcopt")
```

## Quick Start

``` r
library(arcopt)

# Rosenbrock function - a classic difficult optimization problem
result <- arcopt(
  x0 = c(-1.2, 1),
  fn = function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2,
  gr = function(x) c(
    -2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
    200 * (x[2] - x[1]^2)
  ),
  hess = function(x) matrix(c(
    1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
    -400 * x[1], 200
  ), 2, 2)
)

result$par
#> [1] 1 1
```

## Two Modes: Exact Hessian vs Quasi-Newton

**arcopt** offers two optimization strategies:

### Mode 1: Exact Hessian (Default)

Best when you can compute the Hessian analytically or via automatic
differentiation.

``` r
# Provide fn, gr, and hess
result <- arcopt(x0, fn, gr, hess)
```

### Mode 2: Quasi-Newton (No Hessian Required)

Best when computing the Hessian is expensive or unavailable. Uses
BFGS/SR1 approximations.

``` r
# Just provide fn and gr - no Hessian needed
result <- arcopt(x0, fn, gr, control = list(use_qn = TRUE))
```

### Which Mode Should I Use?

| Scenario                           | Recommendation              |
|------------------------------------|-----------------------------|
| Analytic Hessian available         | Use exact Hessian (default) |
| Hessian via autodiff (torch, Stan) | Use exact Hessian           |
| Hessian expensive (n \> 100)       | Use `use_qn = TRUE`         |
| No Hessian function                | Use `use_qn = TRUE`         |

## Benchmark Results

Tested on 8 standard functions with 10 starting points each (80 test
cases):

| Method            | Convergence | Avg Iterations | Hessian Evals |
|-------------------|-------------|----------------|---------------|
| **Exact Hessian** | 100%        | 22             | 1,443         |
| **BFGS QN**       | 100%        | 28             | 0             |
| **SR1 QN**        | 99%         | 61             | 0             |

**BFGS Quasi-Newton** achieves the same 100% convergence as exact
Hessian with only 26% more iterations and zero Hessian evaluations.

## Box Constraints

``` r
# Optimize with bounds: 0 <= x <= 0.5
result <- arcopt(
  x0 = c(0.1, 0.1),
  fn = fn, gr = gr, hess = hess,
  lower = c(0, 0),
  upper = c(0.5, 0.5)
)
```

## Why ARC?

Standard optimizers struggle when the Hessian is indefinite or nearly
singular—common in:

- Maximum likelihood with complex likelihoods
- Bayesian posterior modes (MAP estimation)
- Nonlinear mixed effects models
- Any model with poor identifiability

**Cubic regularization** automatically handles negative curvature
without line searches or trust-region adjustments. It provably escapes
saddle points and converges to local minima.

### The Key Idea

Instead of the quadratic model used by Newton’s method:

    m(s) = f + g's + (1/2) s'Hs

ARC adds a cubic term:

    m(s) = f + g's + (1/2) s'Hs + (σ/3)||s||³

This cubic term prevents unbounded steps when the Hessian has negative
eigenvalues.

## Common Options

``` r
arcopt(x0, fn, gr, hess,
  control = list(
    # Quasi-Newton mode (no Hessian needed)
    use_qn = TRUE,
    qn_method = "bfgs",      # or "sr1", "lbfgs", "lsr1"

    # Convergence tolerances
    gtol_abs = 1e-5,         # Gradient norm
    maxit = 1000,            # Max iterations

    # Momentum (enabled by default, helps ill-conditioned problems)
    use_momentum = TRUE,

    # Diagnostics
    trace = 1                # Print progress
  )
)
```

## Comparison with optim()

| Feature                    | arcopt              | optim(BFGS)   | optim(L-BFGS-B)         |
|----------------------------|---------------------|---------------|-------------------------|
| Handles indefinite Hessian | ✓                   | ✗             | ✗                       |
| Escapes saddle points      | ✓                   | ✗             | ✗                       |
| Box constraints            | ✓                   | ✗             | ✓                       |
| Hessian-free mode          | ✓                   | ✓             | ✓                       |
| Best for                   | Complex likelihoods | Smooth convex | Large bound-constrained |

## When NOT to Use arcopt

**arcopt** is a **local optimizer**. Use something else for:

- **Global optimization**: Use `DEoptim`, `GenSA`, or `nloptr` with
  global algorithms
- **Derivative-free**: Use `nloptr` with COBYLA or `dfoptim`
- **Very large scale** (n \> 10,000): Use `lbfgs` or sparse methods
- **Linear/Quadratic programming**: Use `lpSolve`, `quadprog`

**Best for**: 2-500 parameters, nonconvex objectives, need robust
convergence.

## Return Value

``` r
result <- arcopt(x0, fn, gr, hess)

result$par          # Optimal parameters
result$value        # Optimal function value
result$gradient     # Gradient at solution
result$hessian      # Hessian at solution (if provided)
result$converged    # TRUE if converged
result$iterations   # Number of iterations
result$evaluations  # List: fn, gr, hess counts
result$message      # Why it stopped
```

## References

**Adaptive Regularization with Cubics (ARC)**

Cartis, C., Gould, N. I. M., & Toint, P. L. (2011). Adaptive cubic
regularisation methods for unconstrained optimization. Part I:
motivation, convergence and numerical results. *Mathematical
Programming*, 127(2), 245-295.

Nesterov, Y., & Polyak, B. T. (2006). Cubic regularization of Newton
method and its global performance. *Mathematical Programming*, 108(1),
177-205.

**Quasi-Newton with Cubic Regularization**

Kamzolov, D., Ziu, K., Agafonov, A., & Takáč, M. (2025). Accelerated
Adaptive Cubic Regularized Quasi-Newton Methods. *Journal of
Optimization Theory and Applications*, 208.
<https://doi.org/10.1007/s10957-025-02804-3>

Kamzolov, D., Ziu, K., Agafonov, A., & Takáč, M. (2023). Cubic
Regularization is the Key! The First Accelerated Quasi-Newton Method
with a Global Convergence Rate of O(k^-2) for Convex Functions.
*arXiv:2302.04987*.

Benson, H. Y., & Shanno, D. F. (2018). Cubic regularization in symmetric
rank-1 quasi-Newton methods. *Mathematical Programming Computation*,
10(4), 457-494.

## License

MIT © 2025 Marcus Waldman

## Contributing

Issues and PRs welcome at:
<https://github.com/marcus-waldman/arcopt/issues>
