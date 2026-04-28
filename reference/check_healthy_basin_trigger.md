# Check Healthy-Basin Trigger

Evaluates the five-signal rule against the current sliding window.

## Usage

``` r
check_healthy_basin_trigger(
  state,
  rho_polish = 0.9,
  lambda_min_polish = 0.001,
  g_decay_polish = 0.5,
  g_inf_floor_polish = 1e-08
)
```

## Arguments

- state:

  List returned by
  [`init_healthy_basin_state`](https://marcus-waldman.github.io/arcopt/reference/init_healthy_basin_state.md)
  and maintained by
  [`update_healthy_basin_state`](https://marcus-waldman.github.io/arcopt/reference/update_healthy_basin_state.md).

- rho_polish:

  Minimum acceptable `rho` throughout the window (default `0.9`).

- lambda_min_polish:

  Minimum acceptable `lambda_min(H)` throughout the window (default
  `1e-3`).

- g_decay_polish:

  Maximum acceptable ratio of consecutive `||g||_inf` values (default
  `0.5` = 2x-per-iter contraction).

- g_inf_floor_polish:

  Absolute lower bound on `||g||_inf` at window start (default `1e-8`);
  prevents firing at convergence.

## Value

`TRUE` if the qn_polish transition should be activated, `FALSE`
otherwise. Returns `FALSE` whenever the window is not yet full or any
diagnostic is non-finite.
