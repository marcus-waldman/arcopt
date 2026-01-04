# Solve Cubic Subproblem with Automatic Solver Selection

Unified dispatcher for cubic subproblem solvers with automatic solver
selection based on problem size and Hessian representation.

## Usage

``` r
solve_cubic_subproblem_dispatch(g, H = NULL, sigma, solver = "auto", ...)
```

## Arguments

- g:

  Gradient vector (length n)

- H:

  Hessian matrix (n x n symmetric, may be indefinite)

- sigma:

  Regularization parameter (positive scalar)

- solver:

  Solver to use: "auto" (default) or "eigen". Auto-selection uses
  eigendecomposition (Algorithm 5a).

- ...:

  Additional arguments passed to specific solvers

## Value

List with:

- s:

  Solution vector

- lambda:

  Lagrange multiplier

- pred_reduction:

  Predicted reduction

- converged:

  TRUE if solver converged

- solver_used:

  Name of solver that was used

## Details

Uses eigendecomposition-based solver (Algorithm 5a) which provides
robust handling of indefinite Hessians and hard cases.

For large-scale matrix-free optimization (n \> 500), see
design/scalable-arcs.qmd for deferred large-scale solver
implementations.
