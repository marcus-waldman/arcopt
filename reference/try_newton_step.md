# Try Newton Step with Gershgorin Check

Attempts to compute and validate a Newton step when the Hessian is
positive definite. Uses Gershgorin circle theorem for a cheap O(n^2)
check before attempting Cholesky factorization.

## Usage

``` r
try_newton_step(g, H)
```

## Arguments

- g:

  Gradient vector (length n)

- H:

  Hessian matrix (n x n symmetric matrix)

## Value

List with components:

- success:

  Logical; TRUE if Newton step computed successfully

- s:

  Newton step vector if successful, NULL otherwise

- pred_reduction:

  Predicted reduction from quadratic model if successful, NULL otherwise

## Details

The algorithm:

1.  Computes Gershgorin lower bound on minimum eigenvalue

2.  If bound \<= 0, returns failure (H may be indefinite)

3.  Attempts Cholesky factorization of H

4.  If successful, solves H \* s = -g for Newton step

The Gershgorin theorem states that for symmetric matrix H, the minimum
eigenvalue satisfies: lambda_min \>= min over i of (H_ii - sum of
abs(H_ij) for j != i). If this lower bound is positive, H is guaranteed
positive definite.
