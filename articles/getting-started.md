# Getting Started with arcopt

**arcopt** implements Adaptive Regularization with Cubics (ARC) for
nonlinear optimization in R. Its target users are statisticians and
applied researchers solving maximum likelihood, MAP, and penalized
regression problems with 2-500 parameters where indefinite Hessians,
saddle points, and ill-conditioning are common.

This vignette walks through a single end-to-end optimization. For an
overview of arcopt’s three adaptive subproblem modes (cubic,
trust-region fallback, and the optional quasi-Newton polish), see
[`vignette("solver-modes", package = "arcopt")`](https://marcus-waldman.github.io/arcopt/articles/solver-modes.md).

## Installation

``` r
# CRAN release
install.packages("arcopt")

# Development version from GitHub (requires the devtools package)
devtools::install_github("marcus-waldman/arcopt")
```

## The interface

[`arcopt()`](https://marcus-waldman.github.io/arcopt/reference/arcopt.md)
has a single entry point:

``` r
arcopt(x0, fn, gr, hess, ..., lower, upper, control)
```

The required arguments are the starting parameter vector `x0`, the
objective function `fn`, its gradient `gr`, and its Hessian `hess`. Box
constraints are supplied via `lower` and `upper`; the `control` list
holds optional tolerances and solver knobs documented in
[`?arcopt`](https://marcus-waldman.github.io/arcopt/reference/arcopt.md).
For problems where only a gradient is available, set
`control = list(use_qn = TRUE)` to dispatch to the quasi-Newton variant,
which seeds its Hessian approximation from a one-time finite-difference
computation; see
[`?arcopt`](https://marcus-waldman.github.io/arcopt/reference/arcopt.md)
for details.

## A first solve: Rosenbrock

``` r
library(arcopt)

rosenbrock <- function(x) {
  (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
}

rosenbrock_gr <- function(x) {
  c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
    200 * (x[2] - x[1]^2))
}

rosenbrock_hess <- function(x) {
  matrix(c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
           -400 * x[1],                     200),
         2, 2)
}

result <- arcopt(
  x0 = c(-1.2, 1),
  fn = rosenbrock,
  gr = rosenbrock_gr,
  hess = rosenbrock_hess
)
```

## Inspecting the result

The returned list reports the optimal parameters, the objective value,
gradient and Hessian at the solution, and convergence diagnostics:

``` r
result$par
#> [1] 0.9999994 0.9999989
result$value
#> [1] 3.138905e-13
result$converged
#> [1] TRUE
result$iterations
#> [1] 28
```

Two diagnostics are worth knowing about for nonconvex problems. The
minimum eigenvalue of the Hessian at the reported solution distinguishes
a local minimum (positive) from a saddle point (negative):

``` r
min(eigen(result$hessian, only.values = TRUE)$values)
#> [1] 0.3993613
```

And the `diagnostics` sublist records which solver mode was active at
termination, plus how many times each adaptive transition fired during
the run:

``` r
str(result$diagnostics)
#> List of 6
#>  $ solver_mode_final          : chr "cubic"
#>  $ ridge_switches             : int 0
#>  $ radius_final               : num NA
#>  $ qn_polish_switches         : int 0
#>  $ qn_polish_reverts          : int 0
#>  $ hess_evals_at_polish_switch: int NA
```

For this well-conditioned quadratic-near-optimum Rosenbrock the run
terminates in cubic mode with no transitions; the trust-region fallback
and quasi-Newton polish modes are dormant safety nets that fire only
when their detectors trigger. The companion vignette (`solver-modes`)
shows examples where each mode is load-bearing.

## Where to go next

- [`?arcopt`](https://marcus-waldman.github.io/arcopt/reference/arcopt.md)
  – full reference for the user-facing controls.
- [`?arcopt_advanced_controls`](https://marcus-waldman.github.io/arcopt/reference/arcopt_advanced_controls.md)
  – the cubic-regularization, trust-region, and polish-mode tuning
  parameters for power users.
- [`vignette("solver-modes", package = "arcopt")`](https://marcus-waldman.github.io/arcopt/articles/solver-modes.md)
  – the three subproblem modes and when each one matters.
