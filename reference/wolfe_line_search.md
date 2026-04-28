# Wolfe Line Search

Performs a line search along a descent direction `d` from point `x`,
finding a step length `alpha` that satisfies the strong Wolfe
conditions. Implements the bracketing + zoom scheme from Nocedal &
Wright (2006), Algorithms 3.5 (bracketing) and 3.6 (zoom).

## Usage

``` r
wolfe_line_search(
  fn,
  gr,
  x,
  d,
  f_x,
  g_x,
  c1 = 1e-04,
  c2 = 0.9,
  alpha_max = 1,
  max_iter = 20,
  ...
)
```

## Arguments

- fn:

  Objective function. Called as `fn(x_new, ...)`.

- gr:

  Gradient function. Called as `gr(x_new, ...)`.

- x:

  Numeric vector; current iterate.

- d:

  Numeric vector; search direction (must be a descent direction, i.e.,
  `sum(g_x * d) < 0`).

- f_x:

  Current objective value `fn(x)`.

- g_x:

  Current gradient `gr(x)`.

- c1:

  Armijo (sufficient decrease) constant; default `1e-4`.

- c2:

  Curvature constant (strong Wolfe); default `0.9` (standard for
  quasi-Newton). Smaller values give stricter curvature conditions; for
  linear-CG use `0.1`.

- alpha_max:

  Maximum step length; default `1.0`. Line search will try
  `alpha_curr = alpha_max` first (appropriate for Newton-like
  directions).

- max_iter:

  Maximum total evaluations (bracketing + zoom); default `20`.

- ...:

  Extra arguments forwarded to `fn` and `gr`.

## Value

A list with components:

- success:

  TRUE if strong Wolfe was satisfied; FALSE otherwise.

- alpha:

  Final step length (best available even on failure).

- x_new:

  `x + alpha * d`.

- f_new:

  `fn(x_new)`.

- g_new:

  `gr(x_new)`.

- evals_f:

  Number of objective evaluations performed.

- evals_g:

  Number of gradient evaluations performed.

- reason:

  String describing success or failure mode.

## Details

Strong Wolfe conditions:

- Armijo (sufficient decrease): \\f(x + \alpha d) \le f(x) + c_1 \alpha
  g_x^\top d\\

- Curvature: \\\|g(x + \alpha d)^\top d\| \le c_2 \|g_x^\top d\|\\

The algorithm first brackets an interval containing a Wolfe point, then
bisects within the bracket to locate one. Bisection is used for
simplicity; quadratic or cubic interpolation as in Nocedal & Wright Sec.
3.5 can be added as a future refinement.

Returns `success = FALSE` when no strong Wolfe alpha is found within
`max_iter` evaluations. The caller can still consult `alpha` (and the
associated trial point) to decide whether to take the best-known step or
reject the iteration.
