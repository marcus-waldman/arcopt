# BFGS (Broyden-Fletcher-Goldfarb-Shanno) Hessian Update

Updates the Hessian approximation using the BFGS formula. BFGS maintains
positive definiteness when the curvature condition y'\*s \> 0 is
satisfied.

## Usage

``` r
update_bfgs(b, s, y)
```

## Arguments

- b:

  Current Hessian approximation (n x n symmetric positive definite
  matrix).

- s:

  Step vector: s = x_new - x_old.

- y:

  Gradient difference: y = g_new - g_old.

## Value

A list with components:

- b:

  Updated Hessian approximation (n x n matrix).

- skipped:

  Logical indicating if update was skipped due to curvature condition
  violation.

## Details

The BFGS update formula is: \$\$B\_{k+1} = B_k - \frac{B_k s s^T
B_k}{s^T B_k s} + \frac{y y^T}{y^T s}\$\$

The update is skipped when the curvature condition \\y^T s \> 0\\ is
violated. This can happen on nonconvex problems. When skipped, the
previous approximation is retained.

Note: BFGS maintains positive definiteness, which wastes ARC's ability
to handle indefinite Hessians. Consider using SR1 for nonconvex
problems.

## References

Algorithm 4a in design/pseudocode.qmd

Nocedal, J., & Wright, S. J. (2006). Numerical Optimization (2nd ed.).
Springer. Section 6.1.
