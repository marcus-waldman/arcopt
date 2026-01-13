# Powell-Damped BFGS Update

Applies Powell damping to ensure the curvature condition is satisfied,
then performs a BFGS update. This guarantees an update even when the
standard curvature condition y's \> 0 fails.

## Usage

``` r
update_bfgs_powell(b, s, y, damping_threshold = 0.2)
```

## Arguments

- b:

  Current Hessian approximation (n x n symmetric matrix).

- s:

  Step vector: s = x_new - x_old.

- y:

  Gradient difference: y = g_new - g_old.

- damping_threshold:

  Threshold for damping (default: 0.2). Damping applied if y's \<
  damping_threshold \* s'Bs.

## Value

A list with components:

- b:

  Updated Hessian approximation (n x n matrix).

- theta:

  Damping parameter used (1.0 = no damping).

- y_damped:

  The damped y vector used in update.

- skipped:

  Logical indicating if update was skipped (only if s'Bs is too small to
  apply damping).

## Details

Powell damping modifies y to y_damped = theta\*y + (1-theta)\*Bs such
that s'\*y_damped \>= damping_threshold \* s'\*Bs, ensuring positive
curvature.

The formula for theta when y's \< damping_threshold \* s'Bs is: theta =
(1 - damping_threshold) \* s'Bs / (s'Bs - y's)

This ensures s'\*y_damped = damping_threshold \* s'Bs exactly.

## References

Powell, M. J. D. (1978). A fast algorithm for nonlinearly constrained
optimization calculations. Numerical Analysis, 144-157.

Algorithm 4a-hybrid in design/pseudocode.qmd
