# Initialize Healthy-Basin Detector State

Initialize Healthy-Basin Detector State

## Usage

``` r
init_healthy_basin_state(window = 5)
```

## Arguments

- window:

  Integer sliding-window length (default 5; shorter than the flat-ridge
  detector's 10 because the basin signal is stronger and
  revert-on-false-positive is cheap).

## Value

A list with empty diagnostic vectors and the window size.

## See also

[`update_healthy_basin_state`](https://marcus-waldman.github.io/arcopt/reference/update_healthy_basin_state.md),
[`check_healthy_basin_trigger`](https://marcus-waldman.github.io/arcopt/reference/check_healthy_basin_trigger.md)
