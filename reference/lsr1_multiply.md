# L-SR1 Hessian-Vector Product

Computes the product B\*v using the L-SR1 compact representation.

## Usage

``` r
lsr1_multiply(history, v)
```

## Arguments

- history:

  L-SR1 history list from `update_lsr1`.

- v:

  Vector to multiply by B.

## Value

The product B\*v as a numeric vector.

## Details

The L-SR1 Hessian has the compact form: \$\$B = \gamma I + \Psi (\Psi^T
S)^{-1} \Psi^T\$\$ where \\\Psi = Y - \gamma S\\ and S, Y are matrices
of stored vectors.

This function computes B\*v in O(m^2 + mn) operations.
