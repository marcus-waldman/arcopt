# Compute Step Using Woodbury Identity

Computes s(lambda) = -(B + lambda \* I) inverse times g using the
Woodbury formula.

## Usage

``` r
woodbury_step(lambda, g, gamma, w_g, w_wt, c_inv, w)
```

## Arguments

- lambda:

  Current regularization parameter.

- g:

  Gradient vector.

- gamma:

  Scaling parameter (B0 = gamma \* I).

- w_g:

  Precomputed w times g.

- w_wt:

  Precomputed w times t(w).

- c_inv:

  Precomputed inverse of c_mat.

- w:

  Low-rank update matrix.

## Value

List with s (step vector) and m_mat (for potential reuse).
