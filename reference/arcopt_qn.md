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

- `qn_fd_refreshes`: Number of FD Hessian refreshes performed (hybrid
  mode only; 0 for other qn_method values)

## Details

The QN-ARC algorithm maintains a quasi-Newton Hessian approximation B_k
and solves the cubic regularization subproblem: \$\$m_k(s) = f_k + g_k^T
s + 1/2 s^T B_k s + sigma_k/3 \|\|s\|\|^3\$\$

### Initialization of B_0

When `hess` is not supplied, `arcopt_qn()` seeds the initial Hessian
approximation with a one-time finite-difference Hessian computed from
the supplied gradient at `x0` (cost: `2 * length(x0)` gradient
evaluations once at startup). This provides negative-curvature
information that the pure-identity initialization used by classical
quasi-Newton methods would miss on saddle-prone problems. To opt out and
recover the identity-initialization convention, pass
`hess = function(x) diag(length(x))` explicitly.

### Control Parameters

In addition to standard arcopt controls:

- `qn_method`: Update method (default: "hybrid"):

  - "hybrid": State-aware routing between SR1-first and BFGS-first
    orderings based on the current B_k's eigenstructure and recent rho_k
    values (see below). Includes automatic FD refresh of B_k when SR1 is
    stuck. Recommended for saddle-prone problems.

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

### Hybrid Routing Parameters

The `"hybrid"` method maintains an internal mode flag that starts in
`"indefinite"` (SR1-first priority) and promotes to `"pd"` (BFGS-first)
once the B_k approximation has been reliably positive definite and
cubic-model predictions have tracked the true objective:

- `qn_route_demote_rho` (default 0.25): rho_k below this counts as a
  "bad" step

- `qn_route_promote_rho` (default 0.5): rho_k above this counts as a
  "good" step

- `qn_route_demote_k` (default 2): consecutive bad steps in "pd" mode
  demote back to "indefinite"

- `qn_route_promote_k` (default 3): consecutive good PD steps in
  "indefinite" mode promote to "pd"

- `qn_fd_refresh_k` (default 3): while in "indefinite" mode, this many
  consecutive bad rho steps rebuild B_k from a fresh FD Hessian

- `qn_stuck_refresh_k` (default 100): while in "indefinite" mode, this
  many iterations without promoting also triggers an FD refresh (safety
  net for secondary-saddle stalls)
