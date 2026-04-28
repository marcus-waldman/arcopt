# Solve Cubic Subproblem via Eigendecomposition (Algorithm 5a)

Solves min_s m(s) = g^T s + 1/2 s^T H s + sigma/3 \|\|s\|\|^3 using full
eigendecomposition with explicit hard-case detection.

## Usage

``` r
solve_cubic_eigen(
  g,
  H,
  sigma,
  max_lambda_iter = 50,
  lambda_tol = 1e-10,
  hard_case_tol = 1e-08
)
```

## Arguments

- g:

  Gradient vector (length n)

- H:

  Hessian matrix (n x n symmetric, may be indefinite)

- sigma:

  Regularization parameter (positive scalar)

- max_lambda_iter:

  Maximum secular equation iterations (default: 50)

- lambda_tol:

  Secular equation convergence tolerance (default: 1e-10)

- hard_case_tol:

  Hard case detection threshold (default: 1e-8)

## Value

List with:

- s:

  Solution vector

- lambda:

  Lagrange multiplier

- pred_reduction:

  Predicted reduction: -g^T s - 1/2 s^T H s - sigma/3 \|\|s\|\|^3

- converged:

  TRUE if secular equation converged

- case_type:

  "easy", "hard", or "zero_gradient"

- eigenvalues:

  Eigenvalues of `H` sorted increasing, or `NULL` for the zero-gradient
  early-return path. Exposed so callers can reuse the spectrum (e.g. for
  `lambda_min` in end-of-iteration detectors) instead of re-decomposing
  `H`.

## Details

Uses eigen() for O(n^3) eigendecomposition. Hard case occurs when
gradient is nearly orthogonal to smallest eigenvector (\|g_1\| \<
tol_hard \* \|\|g\|\|). Recommended for n \<= 500 due to cubic
computational cost.
