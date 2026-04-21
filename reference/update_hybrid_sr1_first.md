# Hybrid Update with SR1-First Routing

Reverse priority of update_hybrid: tries SR1 first (to preserve
indefinite curvature information), then BFGS, then Powell-damped BFGS.
Used by arcopt_qn when the current Hessian approximation is indefinite
or when recent BFGS predictions have been inaccurate.

## Usage

``` r
update_hybrid_sr1_first(
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

Same structure as update_hybrid.
