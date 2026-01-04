# Solve Cubic Subproblem via ARCqK Multi-Shift CG-Lanczos (Algorithm 5b)

Solves min_s m(s) = g^T s + 1/2 s^T H s + sigma/3 \|\|s\|\|^3 using
multi-shift CG-Lanczos with automatic negative curvature detection.

## Usage

``` r
solve_cubic_cg(
  g,
  hess_vec,
  sigma,
  shift_psi = sqrt(10),
  shift_e_lower = -10,
  shift_e_upper = 20,
  max_cg_iter = NULL,
  residual_zeta = 0.5,
  residual_xi = 1e-06
)
```

## Arguments

- g:

  Gradient vector (length n)

- hess_vec:

  Function computing Hessian-vector products: hess_vec(v) -\> Hv

- sigma:

  Regularization parameter (positive scalar)

- shift_psi:

  Shift spacing ratio (default: sqrt(10) ≈ 3.16)

- shift_e_lower:

  Min shift exponent (default: -10 → 10^-10)

- shift_e_upper:

  Max shift exponent (default: 20 → 10^20)

- max_cg_iter:

  Maximum CG iterations (default: min(100, 2\*n))

- residual_zeta:

  Residual tolerance exponent (default: 0.5)

- residual_xi:

  Residual tolerance constant (default: 1e-6)

## Value

List with:

- s:

  Solution vector

- lambda:

  Lagrange multiplier (selected shift)

- pred_reduction:

  Predicted reduction: -g^T s - 1/2 s^T H s - sigma/3 \|\|s\|\|^3

- converged:

  TRUE if CG converged for selected shift

- cg_iters:

  Number of CG iterations performed

## Details

**EXPERIMENTAL SOLVER** - This solver is designed for large-scale
problems (n \> 500) and may perform poorly on small problems. Benchmarks
show it can produce large errors and numerical instability on problems
with n \< 100. Use solve_cubic_eigen() for most applications.

Uses CG-Lanczos with multiple shifts (typically 61) geometrically spaced
from 10^-10 to 10^20. Single Hessian-vector product per iteration shared
across all shifts. Negative curvature detected via gamma_j \< 0.
Recommended only for large-scale problems (n \> 500) or matrix-free
optimization.
