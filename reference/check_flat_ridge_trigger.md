# Check Flat-Ridge Trigger

Evaluates the four-signal rule against the current sliding window.

## Usage

``` r
check_flat_ridge_trigger(
  state,
  sigma_min,
  tol_ridge = 0.001,
  rho_tol = 0.1,
  grad_decrease_max = 0.9,
  g_inf_floor = 1e-06
)
```

## Arguments

- state:

  List returned by
  [`init_flat_ridge_state`](https://marcus-waldman.github.io/arcopt/reference/init_flat_ridge_state.md)
  and maintained by
  [`update_flat_ridge_state`](https://marcus-waldman.github.io/arcopt/reference/update_flat_ridge_state.md).

- sigma_min:

  Minimum regularization floor used by the orchestrator.

- tol_ridge:

  Upper bound on `lambda_min` that counts as "poorly conditioned"
  (default `1e-3`).

- rho_tol:

  Maximum `|rho - 1|` that counts as "near-perfect model" (default
  `0.1`).

- grad_decrease_max:

  Ratio of latest to oldest `||g||_inf` above which the gradient is
  deemed stagnant (default `0.9`).

- g_inf_floor:

  Absolute lower bound on `||g||_inf`; the trigger will not fire below
  this value to avoid firing at true local minima (default `1e-6`).

## Value

`TRUE` if the trust-region fallback should be activated, `FALSE`
otherwise. Returns `FALSE` whenever the window is not yet full or any
diagnostic is non-finite.
