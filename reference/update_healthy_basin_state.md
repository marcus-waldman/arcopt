# Update Healthy-Basin Detector State with One Iteration's Diagnostics

Pushes a new record onto the sliding window and trims the oldest entry
once the window is full.

## Usage

``` r
update_healthy_basin_state(state, used_newton, rho, lambda_min, g_inf)
```

## Arguments

- state:

  List returned by
  [`init_healthy_basin_state`](https://marcus-waldman.github.io/arcopt/reference/init_healthy_basin_state.md).

- used_newton:

  Logical; `TRUE` if the current iteration accepted a Newton step (i.e.,
  `step_type == "newton"`).

- rho:

  Current step-acceptance ratio (NA if step was rejected).

- lambda_min:

  Current smallest Hessian eigenvalue.

- g_inf:

  Current `||g||_inf`.

## Value

The updated state list.
