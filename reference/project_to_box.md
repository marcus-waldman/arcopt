# Project Point to Box Constraints

Projects a point to the feasible region defined by box constraints.

## Usage

``` r
project_to_box(x, lower, upper)
```

## Arguments

- x:

  Point to project

- lower:

  Lower bounds

- upper:

  Upper bounds

## Value

Projected point within bounds

## Details

Component-wise projection: each component is clamped to its bounds
