# Simple Cubic Solver for Scaled Identity Hessian

Solves the cubic subproblem when B = gamma \* I (no QN updates). Uses
closed-form solution: s = -g / (gamma + lambda) with \|\|s\|\| = lambda
/ sigma.

## Usage

``` r
solve_cubic_simple(g, gamma, sigma, tol = 1e-10, max_iter = 50L)
```

## Arguments

- g:

  Gradient vector.

- gamma:

  Positive scalar (B = gamma \* I).

- sigma:

  Cubic regularization parameter.

- tol:

  Convergence tolerance.

- max_iter:

  Maximum iterations.

## Value

List with s, lambda, iterations, converged.
