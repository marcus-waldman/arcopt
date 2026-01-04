# Create Hessian-Vector Product Function from Full Matrix

Wraps a full Hessian matrix in a function interface for use with
solve_cubic_cg.

## Usage

``` r
hess_vec_fd(H)
```

## Arguments

- H:

  Hessian matrix (n x n)

## Value

Function that computes H %\*% v for any vector v
