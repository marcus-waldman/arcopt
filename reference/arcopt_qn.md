# Quasi-Newton ARC Optimizer

Minimizes a nonlinear objective function using Adaptive Regularization
with Cubics and quasi-Newton Hessian approximations. Does not require an
exact Hessian function.

## Usage

``` r
arcopt_qn(
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

  Numeric vector of initial parameter values (length n).

- fn:

  Function that computes the objective function value.

- gr:

  Function that computes the gradient.

- hess:

  Optional Hessian function for hybrid mode. If provided, initializes
  B_0 = H(x_0) and can refresh when approximation degrades.

- ...:

  Additional arguments passed to `fn`, `gr`, and `hess`.

- lower:

  Numeric vector of lower bounds (length n). Default: all -Inf.

- upper:

  Numeric vector of upper bounds (length n). Default: all Inf.

- control:

  List of control parameters (see Details).

## Value

Same structure as arcopt, plus:

- `qn_updates`: Number of successful QN updates

- `qn_skips`: Number of skipped updates

- `qn_restarts`: Number of approximation restarts

## Details

The QN-ARC algorithm maintains a quasi-Newton Hessian approximation B_k
and solves the cubic regularization subproblem: \$\$m_k(s) = f_k + g_k^T
s + 1/2 s^T B_k s + sigma_k/3 \|\|s\|\|^3\$\$

### Control Parameters

In addition to standard arcopt controls:

- `qn_method`: Update method (default: "hybrid"):

  - "hybrid": Tries BFGS, falls back to SR1, then Powell-damped BFGS

  - "sr1": Symmetric Rank-1 (allows indefinite Hessians)

  - "bfgs": Standard BFGS (maintains positive definiteness)

  - "lbfgs": Limited-memory BFGS

  - "lsr1": Limited-memory SR1

  - "lhybrid": Limited-memory hybrid (L-BFGS -\> L-SR1 -\> Powell)

- `bfgs_tol`: Curvature tolerance for BFGS in hybrid mode (default:
  1e-10)

- `qn_memory`: History size for limited-memory methods (default: 10)

- `sr1_skip_tol`: SR1 skip test tolerance (default: 1e-8)

- `sr1_restart_threshold`: Consecutive skips before restart (default: 5)

- `use_accel_qn`: **EXPERIMENTAL** Enable Nesterov acceleration
  (default: FALSE). May improve convergence on strongly convex problems
  but can hurt performance on nonconvex problems. Use with caution.
