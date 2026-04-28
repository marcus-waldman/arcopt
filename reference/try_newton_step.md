# Try Newton Step

Attempts to compute and validate a Newton step when the Hessian is
positive definite. Uses Cholesky factorization to check positive
definiteness and solve the Newton system.

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

1.  Attempts Cholesky factorization of H (fails if H is indefinite)

2.  If successful, solves H \* s = -g for Newton step

3.  Computes predicted reduction from quadratic model

No pre-filtering based on sigma is required - the ratio test in the main
ARC loop naturally determines whether Newton steps are appropriate. See
the design notes in `literature/consensus_reviews/` for the rationale
behind not imposing an explicit sigma threshold on the Newton step.
