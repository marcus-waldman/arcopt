# Woodbury-Based Cubic Subproblem Solver

Solves the cubic regularization subproblem when the Hessian
approximation has low-rank structure from quasi-Newton updates: B =
gamma \* I + w' c_mat w.

## Usage

``` r
solve_cubic_woodbury(g, gamma, w, c_mat, sigma, tol = 1e-10, max_iter = 50L)
```

## Arguments

- g:

  Gradient vector (length n).

- gamma:

  Positive scalar for initial Hessian scaling (B0 = gamma \* I).

- w:

  Low-rank update matrix (m x n, where m \<\< n). If NULL, uses
  eigendecomposition fallback.

- c_mat:

  Small coefficient matrix (m x m, symmetric). May be indefinite for SR1
  updates.

- sigma:

  Positive cubic regularization parameter.

- tol:

  Convergence tolerance for secular equation (default: 1e-10).

- max_iter:

  Maximum Newton/bisection iterations (default: 50).

## Value

A list with components:

- s:

  Step vector (length n).

- lambda:

  Regularization parameter satisfying secular equation.

- iterations:

  Number of iterations used.

- converged:

  Logical indicating convergence.

## Details

For limited-memory quasi-Newton methods (L-BFGS, L-SR1), the Hessian
approximation has the form B = gamma \* I + w' c_mat w where w contains
the secant vectors and c_mat is a small m x m matrix.

This solver uses the Woodbury identity to compute: (B + lambda \* I)^-1
g = (1/(gamma + lambda)) g - (1/(gamma + lambda)^2) w' M^-1 w g where M
= c_mat^-1 + (1/(gamma + lambda)) w w'.

Complexity is O(m^2 n + m^3 log(1/tol)) vs O(n^3) for
eigendecomposition.

## References

Algorithm 5b in design/pseudocode.qmd
