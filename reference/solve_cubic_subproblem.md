# Solve Cubic Regularization Subproblem

Solves the cubic subproblem: min_s m(s) = g^T s + 1/2 s^T H s + sigma/3
\|\|s\|\|^3 using a modified Cholesky factorization approach with
secular equation iteration.

## Usage

``` r
solve_cubic_subproblem(g, H, sigma, max_lambda_iter = 50, lambda_tol = 1e-10)
```

## Arguments

- g:

  Gradient vector (length n)

- H:

  Hessian matrix (n x n symmetric matrix, may be indefinite)

- sigma:

  Regularization parameter (positive scalar)

- max_lambda_iter:

  Maximum iterations for lambda secular equation solver (default: 50)

- lambda_tol:

  Convergence tolerance for secular equation (default: 1e-10)

## Value

List with components:

- s:

  Solution vector (step)

- lambda:

  Final Lagrange multiplier

- pred_reduction:

  Predicted reduction in objective: -g^T s - 1/2 s^T H s - sigma/3
  \|\|s\|\|^3

- converged:

  Logical; TRUE if secular equation converged

## Details

The algorithm solves the secular equation norm(s(lambda)) = lambda/sigma
where s(lambda) = -(H + lambda\*I)^(-1) \* g. This is done using:

1.  Cholesky factorization with iterative diagonal shift for indefinite
    H

2.  Newton-Raphson iteration on the secular equation phi(lambda) =
    \|\|s\|\| - lambda/sigma

When H is positive definite, the algorithm starts with lambda = 0. When
H is indefinite, it automatically adds a diagonal shift until
factorization succeeds.
