# Build L-SR1 Compact Form Matrices

Constructs the w and c_mat matrices for the compact representation B =
gamma \* I + w' c_mat w from L-SR1 history.

## Usage

``` r
lsr1_compact_form(history)
```

## Arguments

- history:

  L-SR1 history from update_lsr1.

## Value

List with w (m x n matrix) and c_mat (m x m matrix), or NULL if history
is empty.

## Details

For L-SR1, the compact form is: B = gamma \* I + Psi \* (Psi' S)^-1 \*
Psi'

where Psi = Y - gamma \* S. This can be written as B = gamma\*I + w'
c_mat w with w = Psi' (so w is m x n) and c_mat = (Psi' S)^-1.
