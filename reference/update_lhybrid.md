# Limited-Memory Hybrid BFGS/SR1 Update

Limited-memory version of
[`update_hybrid`](https://marcus-waldman.github.io/arcopt/reference/update_hybrid.md).
Routes between L-BFGS, L-SR1, and Powell-damped L-BFGS based on
curvature conditions.

## Usage

``` r
update_lhybrid(
  history,
  s,
  y,
  m = 10L,
  bfgs_tol = 1e-10,
  sr1_skip_tol = 1e-08,
  damping_threshold = 0.2
)
```

## Arguments

- history:

  L-Hybrid history list with s, y, gamma. Can be NULL.

- s:

  Step vector: s = x_new - x_old.

- y:

  Gradient difference: y = g_new - g_old.

- m:

  Maximum number of pairs to store (default: 10).

- bfgs_tol:

  BFGS curvature tolerance (default: 1e-10).

- sr1_skip_tol:

  SR1 skip threshold (default: 1e-8).

- damping_threshold:

  Powell damping threshold (default: 0.2).

## Value

A list with components:

- history:

  Updated history structure with s, y, gamma.

- update_type:

  "lbfgs", "lsr1", or "powell".

- theta:

  Damping parameter (1.0 if no damping applied).

## Details

Routing priority:

1.  L-BFGS if y'\*s \> bfgs_tol (stable, positive definite)

2.  L-SR1 if SR1 denominator OK (allows indefiniteness)

3.  Powell-damped L-BFGS (guaranteed update)

Uses a unified history format (s, y, gamma) compatible with both
lbfgs_compact_form and lbfgs_multiply_b. Does not use rho since we don't
need the two-loop recursion (lbfgs_multiply).

This is Algorithm 4a-lhybrid from design/pseudocode.qmd.
