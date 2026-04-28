# Solve Trust-Region Subproblem via Eigendecomposition (Algorithm 5c)

Solves min_s m(s) = g^T s + 1/2 s^T H s subject to \|\|s\|\| \<= radius
using full eigendecomposition with explicit hard-case detection.

## Usage

``` r
solve_tr_eigen(
  g,
  H,
  radius,
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

- radius:

  Trust-region radius (positive scalar)

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

  Lagrange multiplier (0 if interior Newton step)

- pred_reduction:

  Predicted reduction: -g^T s - 1/2 s^T H s

- converged:

  TRUE if solver returned a valid step

- case_type:

  "interior", "easy", "hard", or "zero_gradient"

- on_boundary:

  TRUE if \|\|s\|\| == radius (active constraint)

## Details

This is the trust-region counterpart to `solve_cubic_eigen`. The only
structural difference is the constraint: the cubic solver enforces
`lambda = sigma * ||s||`, while the trust-region solver enforces
`||s|| <= radius`. Both reduce to `(H + lambda*I) s = -g` with
`lambda >= max(0, -lambda_min(H))`.

Three cases are handled:

- **Interior**: H is positive definite and the unconstrained Newton step
  `-H^{-1} g` lies inside the trust region. Returned with `lambda = 0`
  and `on_boundary = FALSE`.

- **Easy boundary**: A unique `lambda > max(0, -lambda_min)` places
  `||s(lambda)|| = radius`. Found by Newton on the secular equation.

- **Hard case**: Hessian is indefinite and the gradient is nearly
  orthogonal to the smallest-eigenvalue eigenvector. Resolved per
  More-Sorensen by setting `lambda = -lambda_min(H)` and adding a
  component along the smallest eigenvector to reach the boundary.

Recommended for `n <= 500` due to O(n^3) eigendecomposition cost.
