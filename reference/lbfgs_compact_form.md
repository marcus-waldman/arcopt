# Build L-BFGS Compact Form Matrices

Constructs the w and c_mat matrices for the compact representation B =
gamma \* I + w' c_mat w from L-BFGS history.

## Usage

``` r
lbfgs_compact_form(history)
```

## Arguments

- history:

  L-BFGS history from update_lbfgs.

## Value

List with w (m x n matrix) and c_mat (m x m matrix), or NULL if history
is empty.

## Details

The compact L-BFGS form uses a simplified SR1-like formulation: B =
gamma \* I + Psi (Psi' S)^-1 Psi' where Psi = Y - gamma \* S.

## References

Nocedal, J., & Wright, S. J. (2006). Numerical Optimization (2nd ed.).
Springer. Section 7.2.
