# Detect Stagnation and Recommend Action

Monitors optimization progress and detects stagnation based on small
consecutive steps or function changes. Returns recommended action:
continue, refresh Hessian, or stop.

## Usage

``` r
detect_stagnation(
  step_norms,
  f_values,
  tol_step = 1e-12,
  tol_f = 1e-14,
  max_stagnant = 5,
  is_qn = FALSE,
  already_refreshed = FALSE
)
```

## Arguments

- step_norms:

  Vector of recent step norms (most recent last)

- f_values:

  Vector of recent function values (most recent last)

- tol_step:

  Step size tolerance for stagnation (default: 1e-12)

- tol_f:

  Function change tolerance for stagnation (default: 1e-14)

- max_stagnant:

  Maximum consecutive stagnant iterations before action (default: 5)

- is_qn:

  Logical; TRUE if using quasi-Newton (enables Hessian refresh option)

- already_refreshed:

  Logical; TRUE if Hessian was already refreshed once

## Value

Character; one of "continue", "refresh_hessian", or "stop"

## Details

Stagnation is detected when either:

- Step norms are consecutively below tol_step for max_stagnant
  iterations

- Function changes are consecutively below tol_f for max_stagnant
  iterations

When stagnation is detected:

- If using quasi-Newton and not already refreshed: return
  "refresh_hessian"

- Otherwise: return "stop" (declare convergence failure)
