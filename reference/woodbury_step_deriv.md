# Compute Step Derivative for Newton Update

Computes the second-order inverse (B + lambda*I)^-2 times g, which
equals (B + lambda*I)^-1 times s, for the Newton derivative of the
secular equation.

## Usage

``` r
woodbury_step_deriv(lambda, s, gamma, w_wt, c_inv, w)
```

## Arguments

- lambda:

  Current regularization parameter.

- s:

  Current step vector.

- gamma:

  Scaling parameter.

- w_wt:

  Precomputed w times t(w).

- c_inv:

  Precomputed inverse of c_mat.

- w:

  Low-rank update matrix.

## Value

Vector equal to (B + lambda\*I) inverse times s.
