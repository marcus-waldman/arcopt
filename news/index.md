# Changelog

## arcopt 0.3.0 — first CRAN release

This is the first CRAN release of **arcopt**. The headline feature is
the tri-modal solver shipped in v0.2.0 (cubic / trust-region fallback /
quasi-Newton polish), wrapped with the tiered control surface documented
under
[`?arcopt`](https://marcus-waldman.github.io/arcopt/reference/arcopt.md)
(Tier 1) and
[`?arcopt_advanced_controls`](https://marcus-waldman.github.io/arcopt/reference/arcopt_advanced_controls.md)
(Tier 2). See the v0.2.0 entry below for the full feature description.

### Other changes since v0.2.0

- Cubic solver now reuses its eigendecomposition at end-of-iteration
  detector checks, removing redundant work and trimming roughly 0.3% of
  total runtime on large-scale benchmarks.

- Added benchmarking and profiling scaffolding under `benchmarks/` for
  identifying C++ port candidates. Excluded from the package build;
  available in the GitHub source tree.

## arcopt 0.2.0 — user-experience overhaul

This release reduces the user-facing control surface and the return-list
clutter without changing default optimization behavior on any tested
problem. The `arcopt(x0, fn, gr, hess)` core signature is unchanged.

### Tiered control surface

- The
  [`?arcopt`](https://marcus-waldman.github.io/arcopt/reference/arcopt.md)
  help page now documents only the seven Tier 1 controls: `maxit`,
  `gtol_abs`, `ftol_abs`, `xtol_abs`, `trace`, `verbose`, and `use_qn`.
  All other tunable parameters are documented on a new
  [`?arcopt_advanced_controls`](https://marcus-waldman.github.io/arcopt/reference/arcopt_advanced_controls.md)
  topic, organized by subsystem (cubic regularization, trust-region
  fallback, quasi-Newton polish, quasi-Newton variant). They remain
  settable via `control = list(...)` exactly as before.

- New `verbose = FALSE` Tier 1 control: when `TRUE`, prints one line per
  iteration to the console (`iter`, `f`, `||g||_inf`, `rho`,
  regularization scale, mode). Orthogonal to `trace`, which still
  controls the depth of data captured into `result$trace`.

- [`arcopt_qn()`](https://marcus-waldman.github.io/arcopt/reference/arcopt_qn.md)
  is no longer publicly exported. It remains an internal function
  reachable via `control$use_qn = TRUE` from the main
  [`arcopt()`](https://marcus-waldman.github.io/arcopt/reference/arcopt.md)
  entry point. The internal documentation is still browsable via
  [`?arcopt_qn`](https://marcus-waldman.github.io/arcopt/reference/arcopt_qn.md).

### Diagnostics nesting (breaking change)

- The mode-dispatch diagnostic fields previously returned at the top
  level of the result list — `solver_mode_final`, `ridge_switches`,
  `radius_final`, `qn_polish_switches`, `qn_polish_reverts`, and
  `hess_evals_at_polish_switch` — are now nested under
  `result$diagnostics`. Update downstream code from
  `result$solver_mode_final` to `result$diagnostics$solver_mode_final`,
  etc.

- [`arcopt_qn()`](https://marcus-waldman.github.io/arcopt/reference/arcopt_qn.md)
  runs additionally place `qn_updates`, `qn_skips`, `qn_restarts`, and
  `qn_fd_refreshes` inside the same `diagnostics` sublist (previously
  top-level).

### Removed controls

- `use_momentum`, `momentum_tau`, `momentum_alpha1`, `momentum_alpha2`
  and the associated bisection-search code (~70 LOC of Gao et al. 2022
  ARCm) have been removed. The mechanism showed mixed empirical results
  in our test suite and was off by default. If you need it, recover the
  v0.1.1 implementation from git history.

- `cubic_solver` is no longer accepted as a documented control. The
  internal dispatcher always selects the eigendecomposition solver
  (Algorithm 5a). The dispatcher itself is preserved as an internal
  extensibility hook for the deferred large-scale Algorithm 5b.

## arcopt 0.1.1 (development)

### QN-polish fallback for the strongly-convex basin

- [`arcopt()`](https://marcus-waldman.github.io/arcopt/reference/arcopt.md)
  now supports an opt-in bidirectional cubic \<-\> qn_polish transition
  (`control$qn_polish_enabled`, default `FALSE`) that replaces the cubic
  subproblem with a Wolfe line search along the BFGS-approximated Newton
  direction once the iterate has entered the quadratic attraction basin
  of a strict local minimum. At the switch, the current exact Hessian
  $H_{k}$ warmstarts the BFGS approximation $B_{0}$; subsequent
  iterations do not call `hess()` unless we revert.

- Motivation: near a strict minimum where $\sigma$ has decayed to its
  floor, the cubic penalty is numerically inactive yet arcopt keeps
  evaluating `hess()`. For expensive Hessians (AD via Stan, FD) the
  per-iteration cost dominates the polish phase. QN-polish skips these
  evaluations and preserves Newton’s local quadratic convergence via the
  exact-H warm start plus BFGS’s superlinear updates.

- Healthy-basin detector fires when all five signals hold over a sliding
  window (default width 5): Newton step accepted every iteration,
  $\rho_{k} \geq 0.9$,
  $\lambda_{\min}\left( H_{k} \right) \geq 10^{- 3}$ (strict
  well-conditioned PD), superlinear gradient decay
  ($\parallel g_{k} \parallel_{\infty}/ \parallel g_{k - 1} \parallel_{\infty} \leq 0.5$),
  and $\parallel g_{0} \parallel_{\infty}$ above a convergence floor
  ($10^{- 8}$).

- Bidirectional: if BFGS drifts non-PD (Cholesky fails), the search
  direction is not descent, or the Wolfe line search fails on
  `qn_polish_max_fail` consecutive iterations (default 3), the solver
  reverts to cubic mode, recomputes $H\left( x_{k} \right)$, and enters
  a cooldown period (`qn_polish_reenter_delay`, default 5 cubic
  iterations) before qn_polish may re-fire. Prevents ping-ponging while
  letting the detector catch the basin signal if cubic re-stabilizes.

- New control parameters: `qn_polish_enabled`, `qn_polish_window`,
  `qn_polish_rho`, `qn_polish_lambda_min`, `qn_polish_g_decay`,
  `qn_polish_g_inf_floor`, `qn_polish_c1`, `qn_polish_c2`,
  `qn_polish_alpha_max`, `qn_polish_max_ls_iter`, `qn_polish_max_fail`,
  `qn_polish_reenter_delay`, `qn_polish_curv_eps`.

- New return fields: `qn_polish_switches` (count of cubic-\>qn_polish
  transitions), `qn_polish_reverts` (count of reversions),
  `hess_evals_at_polish_switch` (Hessian-eval count at first switch;
  difference with final `evaluations$hess` quantifies savings).
  `solver_mode_final` now includes `"qn_polish"` as a possible value.

- New internal functions:
  [`wolfe_line_search()`](https://marcus-waldman.github.io/arcopt/reference/wolfe_line_search.md)
  (Nocedal & Wright 2006, Algorithms 3.5 + 3.6),
  [`init_healthy_basin_state()`](https://marcus-waldman.github.io/arcopt/reference/init_healthy_basin_state.md),
  [`update_healthy_basin_state()`](https://marcus-waldman.github.io/arcopt/reference/update_healthy_basin_state.md),
  [`check_healthy_basin_trigger()`](https://marcus-waldman.github.io/arcopt/reference/check_healthy_basin_trigger.md).

- Empirical: on a smooth non-quadratic convex problem
  ($f(x) = \sum\left( 0.5x_{i}^{2} + x_{i}^{4}/12 \right)$, $n = 5$,
  $x_{0} = 10 \cdot \mathbf{1}$) polish mode with default thresholds
  reduces Hessian evaluations from 9 to 6 (33%); with loose thresholds
  ($W = 3$, $\rho_{\text{polish}} = 0.3$) from 9 to 4 (56%). On
  Rosenbrock the detector does not fire under default settings because
  the banana-valley iterations produce intermittent $\rho$ dips;
  Rosenbrock converges as before via the cubic-only path.

### Trust-region fallback for flat-ridge stagnation

- [`arcopt()`](https://marcus-waldman.github.io/arcopt/reference/arcopt.md)
  and
  [`arcopt_qn()`](https://marcus-waldman.github.io/arcopt/reference/arcopt_qn.md)
  now include a one-way cubic-to-trust-region fallback that activates
  when the cubic-regularization subproblem stalls in a near-singular
  positive-definite regime (“flat ridge”). The switch fires when all
  four runtime signals hold for a sliding window: (i) $\sigma_{k}$
  pinned at its floor, (ii) $\left| \rho_{k} - 1 \right|$ small (model
  is accurate), (iii) $\parallel g_{k} \parallel_{\infty}$ stagnant
  across the window, and (iv) the smallest Hessian eigenvalue strictly
  positive but below `tr_fallback_tol_ridge`. The detector is an
  empirical proxy for violation of the local-error-bound (EB) condition
  under which cubic regularization is guaranteed to converge
  quadratically (Yue, Zhou & So 2018); when all four fire, convergence
  theory for pure cubic regularization has no teeth, and a trust-region
  step-norm constraint is a principled alternative.

- New control parameters (all tunable, sensible defaults):

  - `tr_fallback_enabled` (default `TRUE`) — one-way cubic-\>TR switch.
  - `tr_fallback_window` (default `10`), `tr_fallback_tol_ridge`
    (default `1e-3`), `tr_fallback_rho_tol` (default `0.1`),
    `tr_fallback_grad_decrease_max` (default `0.9`),
    `tr_fallback_g_inf_floor` (default `1e-6`) — detector parameters.
  - `tr_rmax` (default `1e6`), `tr_eta1` (default `0.25`), `tr_eta2`
    (default `0.75`), `tr_gamma_shrink` (default `0.25`),
    `tr_gamma_grow` (default `2.0`) — trust-region update parameters.

- New return fields: `solver_mode_final` (`"cubic"` or `"tr"`),
  `ridge_switches` (count of cubic-\>TR transitions, 0 or 1 in v1),
  `radius_final` (final trust-region radius if the solver switched to TR
  mode).

- New internal functions:
  [`solve_tr_eigen()`](https://marcus-waldman.github.io/arcopt/reference/solve_tr_eigen.md)
  (trust-region subproblem via eigendecomposition, mirrors
  [`solve_cubic_eigen()`](https://marcus-waldman.github.io/arcopt/reference/solve_cubic_eigen.md)),
  [`init_flat_ridge_state()`](https://marcus-waldman.github.io/arcopt/reference/init_flat_ridge_state.md),
  [`update_flat_ridge_state()`](https://marcus-waldman.github.io/arcopt/reference/update_flat_ridge_state.md),
  [`check_flat_ridge_trigger()`](https://marcus-waldman.github.io/arcopt/reference/check_flat_ridge_trigger.md)
  (detector state machine).

- [`solve_tr_eigen()`](https://marcus-waldman.github.io/arcopt/reference/solve_tr_eigen.md)
  uses [`stats::uniroot`](https://rdrr.io/r/stats/uniroot.html) on
  $\phi(\lambda) = \parallel s(\lambda) \parallel - r$ for the easy-case
  secular equation rather than plain Newton. Plain Newton was
  numerically unstable on ill-conditioned indefinite Hessians
  (overshoots into the lower-bound neighborhood where
  $\parallel s(\lambda) \parallel$ blows up); `uniroot` is
  unconditionally stable given a valid bracket.

- Signal 4 of the flat-ridge detector was broadened from
  $0 < \lambda_{\min}(H) < \texttt{𝚝𝚘𝚕\_𝚛𝚒𝚍𝚐𝚎}$ to
  $\lambda_{\min}(H) < \texttt{𝚝𝚘𝚕\_𝚛𝚒𝚍𝚐𝚎}$ — the detector now fires on
  both the classical flat-ridge (small positive $\lambda_{\min}$) and
  stuck-at-indefinite-saddle (negative $\lambda_{\min}$) regimes, since
  cubic regularization loses its grip in both.

- The $\left. \sigma\leftrightarrow r \right.$ duality between cubic
  regularization and trust-region step bounds is formally established in
  Dussault (2018; ARCq) and Martínez & Raydan (2017); this release’s
  adaptive switch between the two subproblem formulations appears to be
  novel. See `todo/tr-fallback-hybrid-briefing.md` for design and
  literature context.

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
