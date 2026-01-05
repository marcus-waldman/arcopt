# L-BFGS Two-Loop Recursion

Computes the product B\*g using the two-loop recursion algorithm, where
B is the L-BFGS approximation represented by the history.

## Usage

``` r
lbfgs_multiply(history, g)
```

## Arguments

- history:

  L-BFGS history list from `update_lbfgs`.

- g:

  Vector to multiply by B (typically the gradient).

## Value

The product B\*g as a numeric vector.

## Details

The two-loop recursion computes B\*g in O(mn) operations where m is the
history length and n is the dimension, avoiding explicit formation of
the n x n matrix B.

The algorithm:

1.  Backward loop: compute alpha_i and update q

2.  Initial scaling: r = gamma \* q

3.  Forward loop: update r using beta_i

## References

Nocedal, J., & Wright, S. J. (2006). Numerical Optimization (2nd ed.).
Springer. Algorithm 7.4.
