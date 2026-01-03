# Check Convergence Criteria (Stan-style)

Checks multiple convergence criteria following Stan's L-BFGS
implementation. Termination occurs when **any** criterion is satisfied,
preventing premature termination on ill-scaled problems while ensuring
robust detection.

## Usage

``` r
check_convergence(
  g_current,
  f_current,
  f_previous,
  x_current,
  x_previous,
  iter,
  tol_grad = 1e-08,
  tol_rel_grad = 1e-06,
  tol_obj = 1e-12,
  tol_rel_obj = 1e-08,
  tol_param = 1e-08,
  max_iter = 1000
)
```

## Arguments

- g_current:

  Current gradient vector (length Q)

- f_current:

  Current objective function value (scalar)

- f_previous:

  Previous objective function value (scalar), NA for first iteration

- x_current:

  Current parameter vector (length Q)

- x_previous:

  Previous parameter vector (length Q), NULL for first iteration

- iter:

  Current iteration number (0-indexed)

- tol_grad:

  Absolute gradient tolerance

- tol_rel_grad:

  Relative gradient tolerance

- tol_obj:

  Absolute objective change tolerance

- tol_rel_obj:

  Relative objective change tolerance

- tol_param:

  Parameter change tolerance (infinity norm)

- max_iter:

  Maximum iterations allowed

## Value

List with two components:

- converged:

  Logical; TRUE if any criterion satisfied

- reason:

  Character; description of which criterion triggered, or "" if not
  converged

## Details

The function checks 6 convergence criteria in order:

1.  **max_iter**: Iteration limit reached

2.  **gradient_abs**: max(abs(g)) \< tol_grad

3.  **gradient_rel**: max(abs(g)) \< tol_rel_grad \* max(1, abs(f))

4.  **objective_abs**: abs(f_current - f_previous) \< tol_obj

5.  **objective_rel**: abs(f_current - f_previous) \< tol_rel_obj \*
    max(1, abs(f_current))

6.  **parameter**: max(abs(x_current - x_previous)) \< tol_param
