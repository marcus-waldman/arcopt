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
  vector of length Q and return a Q×Q symmetric matrix. Required.

- lower:

  Numeric vector of lower bounds (length Q). Use `-Inf` for unbounded
  parameters. Default: all `-Inf`.

- upper:

  Numeric vector of upper bounds (length Q). Use `Inf` for unbounded
  parameters. Default: all `Inf`.

- control:

  List of control parameters (see Details).

## Value

A list with components:

- `par`: Optimal parameter vector

- `value`: Optimal function value

- `gradient`: Gradient at optimum

- `hessian`: Hessian at optimum (if `hess` provided)

- `converged`: Logical, whether convergence criteria met

- `iterations`: Number of iterations performed

- `evaluations`: List with `fn`, `gr`, and `hess` evaluation counts

- `message`: Convergence message

## Details

The ARC algorithm iteratively minimizes a cubic regularization model:
\$\$m_k(s) = f_k + g_k^T s + \frac{1}{2} s^T H_k s + \frac{\sigma_k}{3}
\\s\\^3\$\$

where \\\sigma_k\\ is adaptively adjusted based on model accuracy.

### Hessian Requirement

ARC critically depends on accurate curvature information. The `hess`
argument is required. For convenience, finite-difference Hessians can be
computed automatically with `control$use_fd = TRUE`, but analytic
Hessians are strongly recommended for best performance.

### Control Parameters

The `control` list accepts:

- `maxit`: Maximum iterations (default: 1000)

- `ftol_abs`: Absolute function tolerance (default: 1e-8)

- `gtol_abs`: Absolute gradient norm tolerance (default: 1e-5)

- `xtol_abs`: Absolute step size tolerance (default: 1e-8)

- `sigma0`: Initial regularization parameter (default: 1.0)

- `eta1`: Acceptance threshold for step (default: 0.1)

- `eta2`: Very successful step threshold (default: 0.9)

- `gamma1`: Regularization decrease factor (default: 0.5)

- `gamma2`: Regularization increase factor (default: 2.0)

- `cubic_solver`: Solver selection: "auto" (recommended) or "eigen"
  (default: "auto"). Auto-selection uses eigendecomposition (Algorithm
  5a) for robust handling of indefinite Hessians and hard cases.

- `use_momentum`: Enable momentum acceleration (default: FALSE).
  Implements Gao et al. (2022) ARCm with recursive momentum and
  bisection search for monotonicity. **Only recommended for known
  ill-conditioned problems** - on well-conditioned problems, the
  bisection overhead may negate iteration savings. Empirically shows
  mixed results: can dramatically reduce iterations on some problems
  while increasing them on others.

- `momentum_tau`: Maximum momentum parameter (default: 0.5, Gao's τ)

- `momentum_alpha1`: Linear step scaling constant (default: 0.1, Gao's
  α₁)

- `momentum_alpha2`: Quadratic step scaling constant (default: 1.0,
  Gao's α₂)

- `trace`: Print iteration progress (default: FALSE)

- `use_qn`: Use quasi-Newton Hessian approximation (default: FALSE).
  When TRUE, routes to arcopt_qn which uses QN updates instead of exact
  Hessians. See `qn_method` for available QN methods.

- `qn_method`: QN update method - "sr1", "bfgs", "lbfgs", "lsr1"
  (default: "sr1"). Only used when `use_qn = TRUE`.

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
