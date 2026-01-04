# Solve Cubic Subproblem with Automatic Solver Selection

Unified dispatcher for cubic subproblem solvers with automatic solver
selection based on problem size and Hessian representation.

## Usage

``` r
solve_cubic_subproblem_dispatch(
  g,
  H = NULL,
  hess_vec = NULL,
  sigma,
  solver = "auto",
  solver_threshold = 500,
  ...
)
```

## Arguments

- g:

  Gradient vector (length n)

- H:

  Hessian matrix (n x n symmetric, may be indefinite). Optional if
  hess_vec is provided.

- hess_vec:

  Function computing Hessian-vector products: hess_vec(v) -\> Hv.
  Optional if H is provided.

- sigma:

  Regularization parameter (positive scalar)

- solver:

  Solver to use: "auto" (default), "ldl", "eigen", or "cg"

- solver_threshold:

  Threshold for auto-selection (default: 500). Problems with n \<=
  threshold use eigendecomposition, larger problems use CG.

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

Auto-selection logic:

- If hess_vec is provided: use "cg" (matrix-free optimization)

- Else if n \<= solver_threshold: use "eigen" (full eigendecomposition)

- Else: construct hess_vec from H and use "cg" (large-scale)

Manual solver options:

- "ldl": LDL factorization with modified Cholesky (legacy solver)

- "eigen": Eigendecomposition with explicit hard-case handling
  (Algorithm 5a) **Recommended for most applications**

- "cg": ARCqK multi-shift CG-Lanczos (Algorithm 5b) **EXPERIMENTAL - may
  perform poorly on small problems (n \< 500)**

## Performance Notes

The CG solver is designed for large-scale problems (n \> 500) and
matrix-free optimization. Benchmarks show it performs poorly on small
problems due to discrete shift selection and Lanczos breakdown issues.
For best results:

- Use "auto" (default) to let the dispatcher choose

- Use "eigen" explicitly for problems with n \<= 500

- Use "cg" only for large-scale or matrix-free problems
