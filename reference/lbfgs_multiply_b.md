# L-BFGS Hessian-Vector Product (B\*v)

Computes B*v (Hessian approximation times vector) from L-BFGS history.
Note: This is different from
[`lbfgs_multiply`](https://marcus-waldman.github.io/arcopt/reference/lbfgs_multiply.md)
which computes H*v (inverse Hessian times vector).

## Usage

``` r
lbfgs_multiply_b(history, v, gamma = NULL)
```

## Arguments

- history:

  L-BFGS history list with components `s`, `y`, `rho`, `gamma`. Can be
  NULL (returns gamma\*v or just v if gamma not available).

- v:

  Vector to multiply (length n).

- gamma:

  Initial scaling. If NULL, uses history\$gamma or 1.0.

## Value

B\*v where B is the L-BFGS Hessian approximation.

## Details

Uses compact form representation: B = gamma*I + Psi \* inv(M) \* Psi'
where Psi = Y - gamma*S and M = Psi' \* S.

This is needed for Powell damping in L-Hybrid, which requires s^T B s.
