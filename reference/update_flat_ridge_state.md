# Update Flat-Ridge Detector State with One Iteration's Diagnostics

Pushes a new record onto the sliding window and trims the oldest entry
once the window is full.

## Usage

``` r
update_flat_ridge_state(state, sigma, rho, g_inf, lambda_min)
```

## Arguments

- state:

  List returned by
  [`init_flat_ridge_state`](https://marcus-waldman.github.io/arcopt/reference/init_flat_ridge_state.md).

- sigma:

  Current regularization parameter.

- rho:

  Current step-acceptance ratio.

- g_inf:

  Current `||g||_inf`.

- lambda_min:

  Current smallest Hessian eigenvalue (may be negative for indefinite
  Hessians; pass `NA_real_` if unavailable).

## Value

The updated state list.
