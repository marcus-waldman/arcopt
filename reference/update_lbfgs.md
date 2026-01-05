# L-BFGS History Update

Updates the L-BFGS history by adding a new (s, y) pair and dropping the
oldest if the memory limit is exceeded.

## Usage

``` r
update_lbfgs(history, s, y, m = 10L)
```

## Arguments

- history:

  List with components `s` (list of step vectors) and `y` (list of
  gradient differences). Can be NULL for first iteration.

- s:

  New step vector: s = x_new - x_old.

- y:

  New gradient difference: y = g_new - g_old.

- m:

  Maximum number of pairs to store (default: 10).

## Value

Updated history list with components:

- s:

  List of step vectors (most recent last).

- y:

  List of gradient differences.

- rho:

  Vector of 1/(y'\*s) values for each pair.

- gamma:

  Scaling factor for initial Hessian: (y'\*s)/(y'\*y).

## Details

L-BFGS stores the m most recent (s, y) pairs and represents the Hessian
implicitly. The scaling gamma is the Barzilai-Borwein estimate from the
most recent pair.
