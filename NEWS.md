# arcopt 0.1.1 (development)

## QN-polish fallback for the strongly-convex basin

* `arcopt()` now supports an opt-in bidirectional cubic <-> qn_polish
  transition (`control$qn_polish_enabled`, default `FALSE`) that
  replaces the cubic subproblem with a Wolfe line search along the
  BFGS-approximated Newton direction once the iterate has entered the
  quadratic attraction basin of a strict local minimum. At the switch,
  the current exact Hessian $H_k$ warmstarts the BFGS approximation
  $B_0$; subsequent iterations do not call `hess()` unless we revert.

* Motivation: near a strict minimum where $\sigma$ has decayed to its
  floor, the cubic penalty is numerically inactive yet arcopt keeps
  evaluating `hess()`. For expensive Hessians (AD via Stan, FD) the
  per-iteration cost dominates the polish phase. QN-polish skips
  these evaluations and preserves Newton's local quadratic convergence
  via the exact-H warm start plus BFGS's superlinear updates.

* Healthy-basin detector fires when all five signals hold over a
  sliding window (default width 5): Newton step accepted every
  iteration, $\rho_k \ge 0.9$, $\lambda_{\min}(H_k) \ge 10^{-3}$ (strict
  well-conditioned PD), superlinear gradient decay ($\|g_k\|_\infty /
  \|g_{k-1}\|_\infty \le 0.5$), and $\|g_0\|_\infty$ above a
  convergence floor ($10^{-8}$).

* Bidirectional: if BFGS drifts non-PD (Cholesky fails), the search
  direction is not descent, or the Wolfe line search fails on
  `qn_polish_max_fail` consecutive iterations (default 3), the solver
  reverts to cubic mode, recomputes $H(x_k)$, and enters a cooldown
  period (`qn_polish_reenter_delay`, default 5 cubic iterations)
  before qn_polish may re-fire. Prevents ping-ponging while letting
  the detector catch the basin signal if cubic re-stabilizes.

* New control parameters: `qn_polish_enabled`, `qn_polish_window`,
  `qn_polish_rho`, `qn_polish_lambda_min`, `qn_polish_g_decay`,
  `qn_polish_g_inf_floor`, `qn_polish_c1`, `qn_polish_c2`,
  `qn_polish_alpha_max`, `qn_polish_max_ls_iter`, `qn_polish_max_fail`,
  `qn_polish_reenter_delay`, `qn_polish_curv_eps`.

* New return fields: `qn_polish_switches` (count of cubic->qn_polish
  transitions), `qn_polish_reverts` (count of reversions),
  `hess_evals_at_polish_switch` (Hessian-eval count at first switch;
  difference with final `evaluations$hess` quantifies savings).
  `solver_mode_final` now includes `"qn_polish"` as a possible value.

* New internal functions: `wolfe_line_search()` (Nocedal & Wright
  2006, Algorithms 3.5 + 3.6), `init_healthy_basin_state()`,
  `update_healthy_basin_state()`, `check_healthy_basin_trigger()`.

* Empirical: on a smooth non-quadratic convex problem
  ($f(x) = \sum (0.5 x_i^2 + x_i^4/12)$, $n=5$, $x_0 = 10 \cdot \mathbf{1}$)
  polish mode with default thresholds reduces Hessian evaluations from
  9 to 6 (33%); with loose thresholds ($W=3$, $\rho_{\text{polish}}=0.3$)
  from 9 to 4 (56%). On Rosenbrock the detector does not fire under
  default settings because the banana-valley iterations produce
  intermittent $\rho$ dips; Rosenbrock converges as before via the
  cubic-only path.

## Trust-region fallback for flat-ridge stagnation

* `arcopt()` and `arcopt_qn()` now include a one-way cubic-to-trust-region
  fallback that activates when the cubic-regularization subproblem stalls
  in a near-singular positive-definite regime ("flat ridge"). The switch
  fires when all four runtime signals hold for a sliding window: (i)
  $\sigma_k$ pinned at its floor, (ii) $|\rho_k - 1|$ small (model is
  accurate), (iii) $\|g_k\|_\infty$ stagnant across the window, and (iv)
  the smallest Hessian eigenvalue strictly positive but below
  `tr_fallback_tol_ridge`. The detector is an empirical proxy for
  violation of the local-error-bound (EB) condition under which cubic
  regularization is guaranteed to converge quadratically
  (Yue, Zhou & So 2018); when all four fire, convergence theory for
  pure cubic regularization has no teeth, and a trust-region step-norm
  constraint is a principled alternative.

* New control parameters (all tunable, sensible defaults):
  - `tr_fallback_enabled` (default `TRUE`) — one-way cubic->TR switch.
  - `tr_fallback_window` (default `10`), `tr_fallback_tol_ridge` (default
    `1e-3`), `tr_fallback_rho_tol` (default `0.1`),
    `tr_fallback_grad_decrease_max` (default `0.9`),
    `tr_fallback_g_inf_floor` (default `1e-6`) — detector parameters.
  - `tr_rmax` (default `1e6`), `tr_eta1` (default `0.25`),
    `tr_eta2` (default `0.75`), `tr_gamma_shrink` (default `0.25`),
    `tr_gamma_grow` (default `2.0`) — trust-region update parameters.

* New return fields: `solver_mode_final` (`"cubic"` or `"tr"`),
  `ridge_switches` (count of cubic->TR transitions, 0 or 1 in v1),
  `radius_final` (final trust-region radius if the solver switched to
  TR mode).

* New internal functions: `solve_tr_eigen()` (trust-region subproblem
  via eigendecomposition, mirrors `solve_cubic_eigen()`),
  `init_flat_ridge_state()`, `update_flat_ridge_state()`,
  `check_flat_ridge_trigger()` (detector state machine).

* `solve_tr_eigen()` uses `stats::uniroot` on
  $\phi(\lambda) = \|s(\lambda)\| - r$ for the easy-case secular
  equation rather than plain Newton. Plain Newton was numerically
  unstable on ill-conditioned indefinite Hessians (overshoots into
  the lower-bound neighborhood where $\|s(\lambda)\|$ blows up);
  `uniroot` is unconditionally stable given a valid bracket.

* Signal 4 of the flat-ridge detector was broadened from
  $0 < \lambda_{\min}(H) < \texttt{tol\_ridge}$ to
  $\lambda_{\min}(H) < \texttt{tol\_ridge}$ — the detector now fires
  on both the classical flat-ridge (small positive $\lambda_{\min}$)
  and stuck-at-indefinite-saddle (negative $\lambda_{\min}$) regimes,
  since cubic regularization loses its grip in both.

* The $\sigma \leftrightarrow r$ duality between cubic regularization
  and trust-region step bounds is formally established in
  Dussault (2018; ARCq) and Martínez & Raydan (2017); this release's
  adaptive switch between the two subproblem formulations appears to be
  novel. See `todo/tr-fallback-hybrid-briefing.md` for design and
  literature context.

## Quasi-Newton improvements

* `arcopt_qn()` now seeds the initial Hessian approximation $B_0$ with
  a one-time finite-difference Hessian computed from the supplied
  gradient when no `hess` function is provided (previously seeded with
  the identity matrix). Costs $2n$ gradient evaluations at startup
  but dramatically improves saddle-escape behavior. To recover the
  classical identity-init BFGS convention, pass
  `hess = function(x) diag(length(x))`.

* `qn_method = "hybrid"` now uses state-aware routing between SR1-first
  and BFGS-first update priorities. The mode switches based on (i) the
  smallest eigenvalue of the current $B_k$ (Cholesky inertia) and
  (ii) the cubic-subproblem ratio $\rho_k$ over recent steps. New
  control parameters:
  - `qn_route_demote_rho` (default 0.25): $\rho$ below this counts as a
    bad step.
  - `qn_route_promote_rho` (default 0.5): $\rho$ above this counts as
    a good step.
  - `qn_route_demote_k` (default 2): consecutive bad steps in `"pd"`
    mode trigger demotion to `"indefinite"`.
  - `qn_route_promote_k` (default 3): consecutive good PD steps in
    `"indefinite"` mode trigger promotion to `"pd"`.

* Two FD-Hessian refresh triggers rebuild $B_k$ from a fresh
  finite-difference Hessian when the iteration is stuck in
  `"indefinite"` mode:
  - `qn_fd_refresh_k` (default 3): $\rho_k$ below `qn_route_demote_rho`
    for this many consecutive steps signals SR1 drift.
  - `qn_stuck_refresh_k` (default 100): iteration has spent this many
    steps in `"indefinite"` mode without promoting, indicating a
    secondary-saddle stall where per-step $\rho$ is acceptable but
    global progress has halted.

* The result list from `arcopt_qn()` now includes `qn_fd_refreshes`,
  the number of FD Hessian refreshes performed during optimization.

* Empirical: on a 6-parameter growth-mixture-model symmetric saddle,
  `qn_method = "hybrid"` now reaches a proper local minimum on 50/50
  (100%) of seeds, matching the full-Hessian solvers `nlminb`, `trust`,
  and arcopt with analytic Hessian. Identity-seeded BFGS in
  `optim()` reaches 29/50 (58%) on the same benchmark.

## Interface changes

* `arcopt()` and `arcopt_qn()` now accept `...` arguments that are
  passed through to `fn`, `gr`, and `hess`, matching the `optim()`
  convention. Tests cover both direct calls and the
  `control = list(use_qn = TRUE)` dispatch.

* `DESCRIPTION`: removed the "linear equality constraints" claim;
  these are now listed in the future-work roadmap. Added `trust` and
  `marqLevAlg` to `Suggests:` to support the comparative benchmarks
  in the manuscript and package vignettes.

## Scope narrowing

* Removed the limited-memory quasi-Newton methods (`"lbfgs"`,
  `"lsr1"`, `"lhybrid"`) and the Woodbury-identity cubic subproblem
  solver (`R/cubic_woodbury.R`). These variants were intended to
  scale to larger problems, but the cubic subproblem solver is
  still `O(n^3)` through the eigendecomposition path, so the
  limited-memory `B_k` did not deliver a scalability advantage in
  practice. The removal keeps `arcopt_qn()` focused on the 2-500
  parameter regime documented in the package philosophy.

* Valid `qn_method` values are now `"bfgs"`, `"sr1"`, and `"hybrid"`.

* The full limited-memory implementation, including the Woodbury
  solver and matching tests, is preserved on the `scalable-arc`
  branch for future work on matrix-free / large-scale ARC (see
  `design/scalable-arcs.qmd`).

# arcopt 0.1.0

* Initial development version
* CRAN-compliant package structure with Rcpp integration
* Tidyverse style enforcement via lintr
* GitHub Actions CI/CD for automated testing and checks
