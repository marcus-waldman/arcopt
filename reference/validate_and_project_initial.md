# Validate and Project Initial Point

Validates bounds and projects initial point to feasible region if
needed.

## Usage

``` r
validate_and_project_initial(x0, lower, upper)
```

## Arguments

- x0:

  User-provided initial point

- lower:

  Lower bounds

- upper:

  Upper bounds

## Value

Feasible initial point

## Details

This implements Algorithm 7b from the design specification:

1.  Validate bounds (lower \<= upper, correct dimensions)

2.  Check if x0 is feasible

3.  Project to box if infeasible (with warning)
