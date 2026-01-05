# L-SR1 History Update

Updates the L-SR1 history by adding a new (s, y) pair if it passes the
skip test.

## Usage

``` r
update_lsr1(history, s, y, m = 10L, skip_tol = 1e-08, b_times_s = NULL)
```

## Arguments

- history:

  List with components `s`, `y`, and `gamma`. Can be NULL.

- s:

  New step vector: s = x_new - x_old.

- y:

  New gradient difference: y = g_new - g_old.

- m:

  Maximum number of pairs to store (default: 10).

- skip_tol:

  Tolerance for skip test (default: 1e-8).

- b_times_s:

  Product B*s where B is current approximation. If NULL, uses gamma*s.

## Value

Updated history list with components:

- s:

  List of step vectors.

- y:

  List of gradient differences.

- gamma:

  Scaling factor for initial Hessian.

- skipped:

  Logical indicating if update was skipped.

## Details

Unlike L-BFGS, L-SR1 can represent indefinite Hessians and uses a
different skip test based on the denominator of the SR1 formula.
