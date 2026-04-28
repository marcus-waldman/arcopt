# Adaptive Regularization using Cubics Optimizer

Minimizes a nonlinear objective function using Adaptive Regularization
with Cubics (ARC). Designed for robust optimization of ill-conditioned,
nonconvex, and indefinite Hessian problems common in statistical
applications.

## Usage

``` r
arcopt(
  x0,
  fn,
  gr,
  hess = NULL,
  ...,
  lower = rep(-Inf, length(x0)),
  upper = rep(Inf, length(x0)),
  control = list()
)
```

## Arguments

- x0:

  Numeric vector of initial parameter values (length Q).

- fn:

  Function that computes the objective function value. Should take a
  numeric vector of length Q and return a scalar.

- gr:

  Function that computes the gradient. Should take a numeric vector of
  length Q and return a numeric vector of length Q. Required.

- hess:

  Function that computes the Hessian matrix. Should take a numeric
  vector of length Q and return a Q-by-Q symmetric matrix. Required
  (unless `control$use_qn = TRUE`; see Details).

- ...:

  Additional arguments passed to `fn`, `gr`, and `hess`.

- lower:

  Numeric vector of lower bounds (length Q). Use `-Inf` for unbounded
  parameters. Default: all `-Inf`.

- upper:

  Numeric vector of upper bounds (length Q). Use `Inf` for unbounded
  parameters. Default: all `Inf`.

- control:

  A named list of control parameters. The user-facing tolerances and
  switches are documented below; advanced regularization tuning, the
  trust-region fallback, and the quasi-Newton polish mode live on a
  separate help page (see `\link{arcopt_advanced_controls}`). Recognized
  entries:

  `maxit`

  :   Maximum number of iterations (default `1000`).

  `gtol_abs`

  :   Absolute gradient-norm tolerance for convergence (default `1e-5`).

  `ftol_abs`

  :   Absolute objective-value tolerance (default `1e-8`).

  `xtol_abs`

  :   Absolute step-size tolerance (default `1e-8`).

  `trace`

  :   Integer in `0:3`. Depth of per-iteration data captured into
      `result$trace`: `0` collects nothing, `1` (the default) collects
      function value and gradient norm, `2` adds sigma, rho, step type
      and reciprocal Hessian condition number, `3` adds the full
      iterate, step, Hessian and convergence-criterion record. This flag
      controls *saved* data only – for live console output see
      `verbose`.

  `verbose`

  :   Logical. If `TRUE`, prints one line per iteration to the console
      showing iteration number, objective value, `||g||_inf`, ratio rho,
      regularization scale, and active solver mode (default `FALSE`).
      Orthogonal to `trace`.

  `use_qn`

  :   Logical. If `TRUE`, route to the quasi-Newton ARC variant, which
      approximates the Hessian via SR1/BFGS updates and does not require
      an analytic `hess` function. See the advanced-controls page for
      QN-specific parameters (default `FALSE`).

  See `\link{arcopt_advanced_controls}` for the full set of advanced
  tuning parameters governing the cubic regularization, the trust-region
  fallback, and the quasi-Newton polish mode.

## Value

A list with components:

- `par`: Optimal parameter vector.

- `value`: Objective value at `par`.

- `gradient`: Gradient at `par`.

- `hessian`: Hessian at `par` (or the final BFGS approximation if the
  run ended in qn_polish mode).

- `sigma`: Final cubic regularization parameter.

- `converged`: Logical; whether convergence criteria were met.

- `iterations`: Number of iterations performed.

- `evaluations`: Named list of `fn`, `gr`, and `hess` evaluation counts.

- `message`: Convergence reason.

- `trace`: Per-iteration trace data (depth controlled by
  `control$trace`); `NULL` when `trace = 0`.

- `diagnostics`: Sublist of internal mode-dispatch diagnostics –
  `solver_mode_final`, `ridge_switches`, `radius_final`,
  `qn_polish_switches`, `qn_polish_reverts`, and
  `hess_evals_at_polish_switch`. See `\link{arcopt_advanced_controls}`
  for the meaning of each field. Most users do not need to inspect this;
  it is preserved for diagnostic and benchmarking use.

## Details

The ARC algorithm iteratively minimizes a cubic regularization model:
\$\$m_k(s) = f_k + g_k^\top s + \frac{1}{2} s^\top H_k s +
\frac{\sigma_k}{3} \\s\\^3\$\$ where \\\sigma_k\\ is adapted from
observed model accuracy. arcopt may transparently fall back to a
trust-region subproblem in flat-ridge regimes and (optionally, opt-in)
to a line-search BFGS polish in the quadratic attraction basin. The
transitions are observable via `result$diagnostics`; the algorithmic
details and tunable thresholds are documented under
`\link{arcopt_advanced_controls}`.

### Hessian Requirement

arcopt is Hessian-centric: an analytic `hess` function is strongly
recommended. If the analytic form is unavailable, set
`control$use_qn = TRUE` to obtain Hessian-free quasi-Newton updates (see
the advanced-controls page).

## See also

[`arcopt_advanced_controls`](https://marcus-waldman.github.io/arcopt/reference/arcopt_advanced_controls.md)
for advanced tuning of the cubic regularization, trust-region fallback,
and quasi-Newton polish mode.

## Examples

``` r
# \donttest{
# Rosenbrock function
rosenbrock <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2

rosenbrock_gr <- function(x) {
  c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
    200 * (x[2] - x[1]^2))
}

rosenbrock_hess <- function(x) {
  matrix(c(
    1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
    -400 * x[1], 200
  ), 2, 2)
}

result <- arcopt(
  x0 = c(-1.2, 1),
  fn = rosenbrock,
  gr = rosenbrock_gr,
  hess = rosenbrock_hess
)

print(result$par)      # Should be near c(1, 1)
#> [1] 0.9999994 0.9999989
print(result$value)    # Should be near 0
#> [1] 3.138905e-13
# }
```
