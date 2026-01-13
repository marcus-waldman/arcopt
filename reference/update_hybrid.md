# Hybrid BFGS/SR1 Update with Automatic Routing

Attempts BFGS first, falls back to SR1, then to Powell-damped BFGS.
Provides robust quasi-Newton updates for both convex and nonconvex
regions.

## Usage

``` r
update_hybrid(
  b,
  s,
  y,
  bfgs_tol = 1e-10,
  sr1_skip_tol = 1e-08,
  skip_count = 0L,
  restart_threshold = 5L
)
```

## Arguments

- b:

  Current Hessian approximation (n x n symmetric matrix).

- s:

  Step vector: s = x_new - x_old.

- y:

  Gradient difference: y = g_new - g_old.

- bfgs_tol:

  Tolerance for BFGS curvature condition (default: 1e-10).

- sr1_skip_tol:

  Tolerance for SR1 skip test (default: 1e-8).

- skip_count:

  Current count of consecutive skipped updates.

- restart_threshold:

  Number of consecutive skips before restart (default: 5).

## Value

A list with components:

- b:

  Updated Hessian approximation (n x n matrix).

- update_type:

  Character: "bfgs", "sr1", "powell", or "skipped".

- skipped:

  Logical indicating if no update was applied.

- skip_count:

  Updated skip counter (reset on successful update).

- restarted:

  Logical indicating if b was reset to scaled identity.

- theta:

  Powell damping parameter (only if update_type = "powell").

## Details

Routing priority:

1.  BFGS if y's \> bfgs_tol (stable, positive definite)

2.  SR1 if \|r's\| \>= sr1_skip_tol \* \|\|r\|\| \* \|\|s\|\| (allows
    indefiniteness)

3.  Powell-damped BFGS otherwise (guaranteed update)

The hybrid approach is ideal for problems that transition between convex
and nonconvex regions.

## References

Algorithm 4a-hybrid in design/pseudocode.qmd
