# Protect Against NaN/Inf in Function Evaluations

Detects NaN or Inf values in function and gradient evaluations and
attempts recovery by increasing regularization to produce smaller steps.

## Usage

``` r
check_finite(f_value, g_value)
```

## Arguments

- f_value:

  Function value to check

- g_value:

  Gradient vector to check

## Value

Logical; TRUE if values are finite and valid, FALSE if NaN/Inf detected

## Details

This is a simple check function. The recovery logic (re-solving with
larger sigma) should be implemented in the main optimization loop.

Checks that:

- Function value is finite (not NaN or Inf)

- All gradient components are finite
