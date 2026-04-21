# Changelog

## arcopt 0.1.1 (development)

### Quasi-Newton improvements

- [`arcopt_qn()`](https://marcus-waldman.github.io/arcopt/reference/arcopt_qn.md)
  now seeds the initial Hessian approximation $B_{0}$ with a one-time
  finite-difference Hessian computed from the supplied gradient when no
  `hess` function is provided (previously seeded with the identity
  matrix). Costs $2n$ gradient evaluations at startup but dramatically
  improves saddle-escape behavior. To recover the classical
  identity-init BFGS convention, pass
  `hess = function(x) diag(length(x))`.

- `qn_method = "hybrid"` now uses state-aware routing between SR1-first
  and BFGS-first update priorities. The mode switches based on (i) the
  smallest eigenvalue of the current $B_{k}$ (Cholesky inertia) and

  2.  the cubic-subproblem ratio $\rho_{k}$ over recent steps. New
      control parameters:

  - `qn_route_demote_rho` (default 0.25): $\rho$ below this counts as a
    bad step.
  - `qn_route_promote_rho` (default 0.5): $\rho$ above this counts as a
    good step.
  - `qn_route_demote_k` (default 2): consecutive bad steps in `"pd"`
    mode trigger demotion to `"indefinite"`.
  - `qn_route_promote_k` (default 3): consecutive good PD steps in
    `"indefinite"` mode trigger promotion to `"pd"`.

- Two FD-Hessian refresh triggers rebuild $B_{k}$ from a fresh
  finite-difference Hessian when the iteration is stuck in
  `"indefinite"` mode:

  - `qn_fd_refresh_k` (default 3): $\rho_{k}$ below
    `qn_route_demote_rho` for this many consecutive steps signals SR1
    drift.
  - `qn_stuck_refresh_k` (default 100): iteration has spent this many
    steps in `"indefinite"` mode without promoting, indicating a
    secondary-saddle stall where per-step $\rho$ is acceptable but
    global progress has halted.

- The result list from
  [`arcopt_qn()`](https://marcus-waldman.github.io/arcopt/reference/arcopt_qn.md)
  now includes `qn_fd_refreshes`, the number of FD Hessian refreshes
  performed during optimization.

- Empirical: on a 6-parameter growth-mixture-model symmetric saddle,
  `qn_method = "hybrid"` now reaches a proper local minimum on 50/50
  (100%) of seeds, matching the full-Hessian solvers `nlminb`, `trust`,
  and arcopt with analytic Hessian. Identity-seeded BFGS in
  [`optim()`](https://rdrr.io/r/stats/optim.html) reaches 29/50 (58%) on
  the same benchmark.

### Interface changes

- [`arcopt()`](https://marcus-waldman.github.io/arcopt/reference/arcopt.md)
  and
  [`arcopt_qn()`](https://marcus-waldman.github.io/arcopt/reference/arcopt_qn.md)
  now accept `...` arguments that are passed through to `fn`, `gr`, and
  `hess`, matching the [`optim()`](https://rdrr.io/r/stats/optim.html)
  convention. Tests cover both direct calls and the
  `control = list(use_qn = TRUE)` dispatch.

- `DESCRIPTION`: removed the “linear equality constraints” claim; these
  are now listed in the future-work roadmap. Added `trust` and
  `marqLevAlg` to `Suggests:` to support the comparative benchmarks in
  the manuscript and package vignettes.

### Scope narrowing

- Removed the limited-memory quasi-Newton methods (`"lbfgs"`, `"lsr1"`,
  `"lhybrid"`) and the Woodbury-identity cubic subproblem solver
  (`R/cubic_woodbury.R`). These variants were intended to scale to
  larger problems, but the cubic subproblem solver is still `O(n^3)`
  through the eigendecomposition path, so the limited-memory `B_k` did
  not deliver a scalability advantage in practice. The removal keeps
  [`arcopt_qn()`](https://marcus-waldman.github.io/arcopt/reference/arcopt_qn.md)
  focused on the 2-500 parameter regime documented in the package
  philosophy.

- Valid `qn_method` values are now `"bfgs"`, `"sr1"`, and `"hybrid"`.

- The full limited-memory implementation, including the Woodbury solver
  and matching tests, is preserved on the `scalable-arc` branch for
  future work on matrix-free / large-scale ARC (see
  `design/scalable-arcs.qmd`).

## arcopt 0.1.0

- Initial development version
- CRAN-compliant package structure with Rcpp integration
- Tidyverse style enforcement via lintr
- GitHub Actions CI/CD for automated testing and checks
