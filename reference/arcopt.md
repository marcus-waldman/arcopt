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
  hess_vec = NULL,
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
  vector of length Q and return a Q×Q symmetric matrix. Either `hess` or
  `hess_vec` must be provided.

- hess_vec:

  Function that computes Hessian-vector products. Should take a numeric
  vector v of length Q and return H\*v (length Q). Allows matrix-free
  optimization for large-scale problems. Optional if `hess` is provided.

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
argument is strongly recommended. If `hess = NULL`, the algorithm falls
back to SR1 quasi-Newton approximation, which is less robust but avoids
Hessian computation.

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

- `use_sr1`: Use SR1 quasi-Newton if `hess = NULL` (default: TRUE)

- `cubic_solver`: Solver selection: "auto" (recommended), "ldl",
  "eigen", "cg" (default: "auto"). The "cg" solver is **experimental**
  and may perform poorly on small problems; use "eigen" for most
  applications.

- `cubic_solver_threshold`: Problem size threshold for auto-selection
  (default: 500)

- `trace`: Print iteration progress (default: FALSE)

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
#> [1] 1 1
print(result$value)    # Should be near 0
#> [1] 7.067483e-18
# }
```
