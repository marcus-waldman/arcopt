# QN-Polish Mode for arcopt — Implementation Briefing

**Purpose of this document**: Self-contained context to start a fresh session
and implement the `"qn_polish"` mode — a third branch in the arcopt solver
state machine that fires when the iterate has entered the quadratic
attraction basin of a local minimum, switches to line-search BFGS with the
exact Hessian as warm-started B₀, and stops evaluating the user's `hess()`
function for the remainder of the run.

**Status**: design agreed; implementation not started. Extends the
cubic↔TR hybrid shipped in commit `a3fed60`
(`feature/tr-fallback-hybrid`).

**Parent branch**: `feature/tr-fallback-hybrid` (merge to `main` first if
desired, or branch off directly).

---

## 1. Why this exists

arcopt is Hessian-centric: every iteration evaluates the user's `hess()`
function (analytic, AD, or FD) and optionally performs an O(n³) Cholesky
for the Newton-first path. Near a strict local minimum where H is
well-conditioned PD, this is overkill — **a line-search BFGS step converges
superlinearly without re-evaluating H at all.**

For expensive Hessians (AD via BridgeStan/Stan takes seconds per call; FD
takes n gradient evaluations) the wasted work in the final polish phase
can dominate wall-clock time. Empirical example: the DPM Heckman MAP
(n = 299) takes ~2 seconds per arcopt iteration, of which ~1 second is
`hess()`. A 50-iteration polish phase spends ~50 seconds computing
Hessians that are barely changing.

Classical QN methods (optim-BFGS, nlminb) skip this cost but can't escape
saddles. arcopt's cubic regularization excels at getting into the basin;
QN excels at converging once there. A mode that dispatches each regime
to its best tool is the natural extension of the cubic↔TR hybrid.

Conceptual framing (extends §3a of `design/design-principles.qmd`):

| Regime | Best subproblem | Rationale |
|---|---|---|
| Saddle / indefinite / far from minimum | Cubic | Intrinsic negative-curvature exploration; optimal complexity |
| Flat ridge / stuck indefinite | Trust-region | σ at floor has no grip; step-norm constraint bounds exploration |
| **Strongly convex basin** | **Line-search BFGS** | **H well-conditioned; exact H is overkill; superlinear convergence without Hessian evals** |

The cubic→TR switch (v1) handles the second regime. The cubic→qn_polish
switch (this briefing, v2) handles the third.

---

## 2. The idea — three-mode solver state

Extend `solver_mode` from `{cubic, tr}` to `{cubic, tr, qn_polish}`:

```
            healthy-basin detector
cubic ────────────────────────────► qn_polish
  │                                      │
  │                                      │ line-search failure
  │                                      │ or BFGS drift
  │                                      ▼
  │◄─────────────────────────────────────┘
  │
  │   flat-ridge / stuck-saddle detector
  └────────────────────────────────► tr
                                       │
                                       │ (optional v3: tr → qn_polish
                                       │  if TR finds a basin)
                                       ▼
                                   [not in this briefing]
```

**Key decisions**:

1. **Bidirectional for cubic ↔ qn_polish**. If the line-search phase
   encounters a regime the BFGS approximation can't handle (e.g., BFGS
   predicts decrease but f increases; line search fails to find a Wolfe
   α within max iterations; or several consecutive skipped curvature
   conditions), revert to cubic and recompute a fresh exact Hessian. This
   is *different* from cubic ↔ tr which is one-way.

2. **Not triggered from tr mode**. If tr has fired, we're already in a
   pathological regime; staying in tr is safer than adding another
   transition. v3 can add tr → qn_polish if empirical evidence warrants.

3. **Not enabled by default in v1 of this feature**. Introduce as
   `qn_polish_enabled = FALSE`; set to `TRUE` as default after benchmarks
   prove it doesn't regress existing problems.

---

## 3. Design — agreed decisions

### 3.1 Healthy-basin detector (fires cubic → qn_polish)

Sliding window of W iterations; fires when **all** five signals hold
throughout the window:

1. **Newton step accepted** at every iteration in the window
   (i.e., `step_type == "newton"` for W consecutive iters; the cubic
   subproblem solver hasn't been invoked because `try_newton_step`
   succeeded and its ρ passed the acceptance test).

2. **High ρ**: `ρ_k ≥ rho_polish` (default `0.9`) throughout the window.
   Quadratic model is predictive.

3. **Hessian well-conditioned PD**: `λ_min(H_k) ≥ lambda_min_polish`
   (default `1e-3`) throughout the window. Not just "chol succeeded"
   (that only guarantees PD, not well-conditioning).

4. **Gradient decreasing superlinearly**:
   `‖g_k‖_∞ / ‖g_{k−1}‖_∞ ≤ g_decay_polish` (default `0.5`) — geometric
   contraction at rate ≥ 2× per iteration. This is the strongest
   signal that we are in the quadratic attraction basin.

5. **Absolute gradient floor crossed**: `‖g_0‖_∞ > g_inf_floor_polish`
   at window start (default `1e-8`). Prevents firing at convergence
   where both grad_decay and Newton-accept are trivially satisfied.

Parameters:

- `qn_polish_window = 5` (shorter than TR-fallback's 10 because the basin
  signal is stronger; false positives are cheaper because we can revert)
- `qn_polish_rho = 0.9`
- `qn_polish_lambda_min = 1e-3`
- `qn_polish_g_decay = 0.5`
- `qn_polish_g_inf_floor = 1e-8`

### 3.2 Hand-off mechanics at cubic → qn_polish transition

- **B₀ initialization**: `B_0 = H(x_current)`. Reuse the most recent
  `h_current` available — no extra Hessian eval.
- **Line-search state**: initialize `alpha_prev = 1.0` (try full Newton
  step first every iteration).
- **No more calls to `hess()`** after the transition — this is the whole
  point of the mode.
- **Track `hess_evals` at transition**: useful diagnostic for the
  "Hessian evaluations saved" metric.

### 3.3 qn_polish main loop

For each iteration in qn_polish mode:

1. **Compute direction** `d_k = -B_k⁻¹ g_k`. Use Cholesky of `B_k` (since
   BFGS preserves PD under Wolfe curvature). If Cholesky fails (drift
   accumulated indefiniteness), **revert to cubic** immediately.
2. **Sanity check**: require `g_k' d_k < 0` (descent direction). If not,
   **revert to cubic**.
3. **Line search** along `d_k` satisfying Wolfe conditions
   (see §3.4 for algorithm).
4. If line search succeeds with some `α`:
   - `x_{k+1} = x_k + α d_k`
   - Evaluate `g_{k+1}`; compute `s_k = α d_k`, `y_k = g_{k+1} − g_k`.
   - BFGS update `B_{k+1}`:
     - If `y_k' s_k > ε_curv · ‖s_k‖ ‖y_k‖` (default `ε_curv = 1e-10`):
       apply standard BFGS formula
       $B_{k+1} = B_k - \frac{B_k s_k s_k^\top B_k}{s_k^\top B_k s_k}
                      + \frac{y_k y_k^\top}{y_k^\top s_k}$
     - Else: skip update (preserves PD of B).
   - Check convergence (standard `check_convergence`).
5. If line search fails (no Wolfe α within `qn_polish_max_ls_iter` =
   default `20`):
   - Increment `qn_polish_fail_count`.
   - If `qn_polish_fail_count ≥ qn_polish_max_fail` (default `3`):
     **revert to cubic**, recompute `H`, zero the counter.
   - Else: reject this iteration; `x` unchanged; try again next iter
     (B_k unchanged).

### 3.4 Wolfe line search (new `R/line_search.R`)

Standard backtracking-then-zoom algorithm (Nocedal & Wright, Alg 3.5–3.6):

```
wolfe_line_search(fn, gr, x, d, f_x, g_x,
                  c1 = 1e-4, c2 = 0.9,
                  alpha_max = 1.0, max_iter = 20, ...)
  returns list(alpha, f_new, g_new, evals_f, evals_g, success)
```

Implements:
- **Armijo (sufficient decrease)**: `f(x + αd) ≤ f(x) + c₁·α·g_x'd`
- **Strong Wolfe (curvature)**: `|g(x + αd)'d| ≤ c₂·|g_x'd|`

Algorithm:
1. Bracketing phase: try α = α_max, α_max/2, α_max/4, … until Armijo
   fails OR curvature satisfied; bracket [α_lo, α_hi] where Armijo breaks.
2. Zoom phase: bisection within bracket until Wolfe found or max_iter
   hit.

~100 LOC. No external dependencies; pure base R.

### 3.5 qn_polish → cubic revert

Triggered by:
- Cholesky of B_k fails (B lost PD despite skip-on-bad-curvature — rare
  but possible via numerical accumulation).
- `g_k' d_k ≥ 0` (direction isn't descent, shouldn't happen with PD B
  but guard against numerical weirdness).
- `qn_polish_fail_count ≥ qn_polish_max_fail` consecutive line-search
  failures.

Revert action:
- `solver_mode ← "cubic"`
- `h_current ← hess(x_current, ...)`  (recompute exact H; `hess_evals++`)
- `qn_polish_revert_count ← qn_polish_revert_count + 1`
- Discard `B_k` (no longer needed; cubic mode uses `h_current`).
- Healthy-basin detector state is NOT reset — the window keeps
  accumulating, so if cubic regains stability quickly we can re-enter
  polish on the same window data. BUT add hysteresis: require
  `qn_polish_reenter_delay` iterations (default `5`) before the
  detector can fire again. Prevents ping-ponging.

### 3.6 Interaction with arcopt_qn

v1 of this mode is **arcopt()-only**. `arcopt_qn()` already operates
without Hessian evals after startup — adding polish mode there would be
redundant (the QN path already has BFGS updates in cubic mode). A later
version could add polish mode to arcopt_qn to replace the cubic
subproblem with a pure line search once basin detected, but this is out
of scope for v1.

### 3.7 Interaction with momentum

Momentum (`use_momentum = TRUE`) is disabled in qn_polish mode. The
momentum buffer is reset at transition (`v_prev ← 0`). Rationale:
momentum is tied to cubic's ρ/σ dynamics; it has no analog in BFGS line
search.

---

## 4. Critical files

### Existing (read to understand)

- `R/arcopt.R` — main loop; add qn_polish branch after the existing
  cubic and tr branches. Roughly: between the stagnation check (STEP 1b)
  and STEP 2.
- `R/flat_ridge.R` — pattern to mirror for healthy-basin detector
  (state machine idiom: `init_*`, `update_*`, `check_*`).
- `R/tr_eigen.R` — neighbor in spirit; a new `R/healthy_basin.R`
  should follow the same style.
- `R/qn_updates.R` — contains `update_bfgs()`; can be reused directly
  for the BFGS update in qn_polish. Verify it exposes the skip-on-bad-
  curvature behavior or wrap with a skip guard in the caller.
- `R/newton_step.R` — `try_newton_step()`; polish mode's acceptance
  predicate (step_type == "newton" && ρ ≥ 0.9) reads from arcopt's
  per-iteration state.
- `todo/tr-fallback-hybrid-briefing.md` — style and structural reference
  for this document; the completed work it describes is the
  architectural precedent.

### New files to create

- `R/healthy_basin.R` — detector state machine mirroring
  `R/flat_ridge.R`. ~80 LOC.
- `R/line_search.R` — `wolfe_line_search()`. ~100 LOC.
- `tests/testthat/test-healthy-basin.R` — unit tests for the detector,
  mirroring `test-flat-ridge.R`.
- `tests/testthat/test-line-search.R` — unit tests for Wolfe line
  search on known-answer quadratic problems.
- `tests/testthat/test-qn-polish-integration.R` — end-to-end tests:
  Rosenbrock reaches polish mode, DPM Heckman probably does NOT (TR
  fires first and stays), strongly convex sphere reaches polish mode
  and converges with minimal Hessian evals.

### Modified files

- `R/arcopt.R`:
  - Add `qn_polish_*` controls to `control_defaults`.
  - Add state variables at init: `solver_mode` already exists; add
    `basin_state`, `b_current_polish`, `qn_polish_switches`,
    `qn_polish_reverts`, `qn_polish_fail_count`,
    `qn_polish_reenter_cooldown`, `hess_evals_at_polish_switch`.
  - Add `if (solver_mode == "qn_polish") { ... }` branch in main loop.
  - After cubic block: update basin_state, check trigger, switch if
    ready.
  - Extend return list: `solver_mode_final` already covers the new
    mode; add `qn_polish_switches`, `qn_polish_reverts`, and
    `hess_evals_saved` (= `hess_evals_at_polish_switch` - whatever a
    full cubic run would have cost, estimable from
    `iterations - switch_iteration`).

- `NEWS.md`: new section for qn_polish mode.
- `design/pseudocode.qmd`: Algorithm 5d (Wolfe line search),
  Algorithm 6d (healthy-basin detector), and extended Algorithm 1 main
  loop showing the three-way dispatch. Update the Algorithm Index.
- `design/design-principles.qmd`: extend §3a (σ↔r duality) into a
  tri-modal framing. The section after it (3. Escape Saddles Naturally)
  remains valid; add a paragraph that polish mode preserves
  Newton-local quadratic convergence while avoiding Hessian recomputation.

---

## 5. Implementation phases

Estimate: 4–6 working days for core + benchmarks. Can split into
smaller branches (`feature/qn-polish-detector`,
`feature/qn-polish-line-search`, `feature/qn-polish-integration`) or
land as one.

### Phase 1 — Wolfe line search

**New file**: `R/line_search.R`

```r
wolfe_line_search(fn, gr, x, d, f_x, g_x,
                  c1 = 1e-4, c2 = 0.9,
                  alpha_max = 1.0, max_iter = 20, ...)
  -> list(alpha, f_new, g_new, x_new, evals_f, evals_g,
          success, reason)
```

Implementation from Nocedal & Wright Alg 3.5 (bracketing) + 3.6 (zoom).
~100 LOC.

**New tests**: `tests/testthat/test-line-search.R`
- Quadratic f(x) = 0.5 x'Ax with PD A: succeeds with α=1
- Nonconvex along d: α ≤ 1 found
- No Wolfe α exists (e.g., infinite descent): returns `success=FALSE`
  within max_iter
- Agreement with optim-BFGS's internal line search on a few configs
  (loose tolerance)

### Phase 2 — Healthy-basin detector

**New file**: `R/healthy_basin.R`

```r
init_healthy_basin_state(window = 5) -> state
update_healthy_basin_state(state, used_newton, rho, g_inf, lambda_min,
                           chol_success) -> state
check_healthy_basin_trigger(state, rho_polish = 0.9,
                            lambda_min_polish = 1e-3,
                            g_decay_polish = 0.5,
                            g_inf_floor_polish = 1e-8) -> logical
```

~80 LOC. Tests mirror `test-flat-ridge.R`:
- Trigger fires when all 5 signals align for full window
- Trigger does NOT fire during cubic-dominated iterations
- Trigger does NOT fire at convergence (grad below floor)
- Trigger does NOT fire during TR mode (Newton not used)
- Hysteresis (cooldown) works after a revert

### Phase 3 — arcopt.R integration

**Edit**: `R/arcopt.R`.

Add to iteration state (before main loop):
- `basin_state <- init_healthy_basin_state(...)`
- `b_current_polish <- NULL`
- `qn_polish_switches <- 0L`
- `qn_polish_reverts <- 0L`
- `qn_polish_fail_count <- 0L`
- `qn_polish_cooldown <- 0L`
- `hess_evals_at_polish_switch <- NA_integer_`

Main-loop logic (pseudocode):

```
if (solver_mode == "qn_polish") {
  # Compute direction
  chol_B <- try(chol(b_current_polish))
  if (Cholesky failed OR descent check fails) {
    revert_to_cubic(); next
  }
  d <- solve(chol_B, -g_current via backsolve)

  # Line search
  ls <- wolfe_line_search(fn, gr, x_current, d, f_current, g_current, ...)
  fn_evals += ls$evals_f; gr_evals += ls$evals_g

  if (!ls$success) {
    qn_polish_fail_count++
    if (qn_polish_fail_count >= qn_polish_max_fail) revert_to_cubic()
    iter++; next
  }
  qn_polish_fail_count <- 0L

  # Accept
  x_previous <- x_current; f_previous <- f_current
  x_current <- ls$x_new; f_current <- ls$f_new; g_current <- ls$g_new

  # BFGS update
  s_vec <- ls$alpha * d
  y_vec <- g_current - g_previous
  if (curvature OK) b_current_polish <- update_bfgs(b_current_polish, s_vec, y_vec)$b
  # else skip (preserves PD)

  step_norms <- c(step_norms, ||s_vec||); f_values <- c(f_values, f_current)
}

else if (solver_mode == "cubic") {
  # Existing code unchanged

  # AFTER sigma update, update detector
  basin_state <- update_healthy_basin_state(basin_state, used_newton, rho, ...,
                                            lambda_min(H), chol_success)
  if (qn_polish_cooldown > 0L) qn_polish_cooldown <- qn_polish_cooldown - 1L

  if (control$qn_polish_enabled &&
      qn_polish_cooldown == 0L &&
      check_healthy_basin_trigger(basin_state, ...)) {
    solver_mode <- "qn_polish"
    b_current_polish <- h_current  # warmstart
    qn_polish_switches <- qn_polish_switches + 1L
    hess_evals_at_polish_switch <- hess_evals
  }

  # Existing tr-fallback detector logic stays the same
}
```

revert_to_cubic() helper:
```
solver_mode <- "cubic"
h_current <- hess(x_current, ...); hess_evals++
qn_polish_reverts <- qn_polish_reverts + 1L
b_current_polish <- NULL
qn_polish_fail_count <- 0L
qn_polish_cooldown <- control$qn_polish_reenter_delay
```

### Phase 4 — Bidirectional + hysteresis

Validate that revert cleanly restores cubic state and that the cooldown
prevents immediate re-entry. Unit test this specifically by constructing
a problem where qn_polish revers and verifying the cooldown.

### Phase 5 — Benchmark validation

**Regression** (polish must NOT degrade these; ideally neutral or
faster):
- Rosenbrock — expected: polish fires by iter ~10, converges 20–30% faster
- Mixture saddle, two-exp, GMM (manuscript examples) — expected: polish
  fires after saddle escape, net faster
- DPM Heckman — expected: TR fires FIRST, stays in TR (polish never
  activates because we're not in basin); same result as
  v1 hybrid (f=675, 125 iters)

**Win cases** (polish should demonstrate meaningful speedup):
- Sphere / sum-of-squares / trid (well-conditioned convex) — polish
  should fire almost immediately, converge in O(log log ε) iters
- Rotated-hyper-ellipsoid (moderately conditioned) — polish should
  fire after 5–10 cubic iters, converge quickly
- Any large-n problem with expensive AD Hessian — measure wall-clock
  Hessian-eval savings

**New synthetic test case**: a problem where cubic runs for ~15 iters
to reach basin, then polish saves ~30 Hessian evals. Good candidate:
smoothed versions of mildly-non-convex objectives.

### Phase 6 — Documentation

- `R/arcopt.R` roxygen: new control params, tri-modal state description
- `NEWS.md`: new 0.3.0 section
- `design/pseudocode.qmd`: Algorithm 5d (Wolfe LS), 6d (healthy-basin
  detector), updated Algorithm 1
- `design/design-principles.qmd`: extend §3a into tri-modal framing
- Manuscript: optional — could extend `§sec-tr-fallback` into
  `§sec-adaptive-mode-selection` describing the full three-mode
  architecture

---

## 6. Starting instructions for a new session

1. **Read this briefing file first**.

2. **Read the tr-fallback briefing** at
   `todo/tr-fallback-hybrid-briefing.md` for architectural precedent and
   Appendix A's Yue-So (2018) review (the EB-condition framing still
   applies — polish mode is the "we are inside EB and it's tight" case).

3. **Read the 2026-04-23 memory** at
   `C:\Users\marcu\.claude\projects\C--Users-marcu-git-repositories-arcopt\memory\tr_fallback_hybrid_implemented_2026-04-23.md`
   for the validation methodology and debugging lessons. In particular,
   note:
   - Use `devtools::load_all(".")` in benchmark scripts, not
     `library(arcopt)` (otherwise you'll run installed code, not dev)
   - The `uniroot` pattern used in `solve_tr_eigen` was chosen after
     plain Newton on the secular equation proved unstable; the
     healthy-basin detector and Wolfe line search don't have that
     problem but should still use robust primitives (Brent, bisection)
     not bespoke Newton iterators where possible.

4. **Ensure clean branch state**:
   ```
   git checkout feature/tr-fallback-hybrid    # or main after merge
   git checkout -b feature/qn-polish-mode
   ```

5. **Baseline tests must pass before any changes**:
   ```
   devtools::test()   # expect 537 passing
   ```

6. **Start with Phase 1** (Wolfe line search) — fewest dependencies, most
   testable in isolation. Then Phase 2 (detector), Phase 3 (arcopt
   integration), Phase 4 (bidirectional logic), Phase 5 (benchmarks),
   Phase 6 (docs).

7. **Verify after each phase** with
   `devtools::test() + lintr::lint_package() + targeted benchmark`.

---

## 7. Open questions / risks

- **BFGS drift vs. occasional Hessian refresh**. In very long polish
  phases (hundreds of BFGS updates) the approximation can accumulate
  error. Possible remedy: refresh B_k from `hess(x_current)` every
  K iterations (default K = ∞, i.e., never, in v1). Add
  `qn_polish_hess_refresh_k` as a control parameter. v1 default: no
  refresh. Monitor in benchmarks; enable if drift visibly hurts.

- **Line search c₂**. Standard BFGS uses c₂ = 0.9 (loose curvature).
  Strict Wolfe would be c₂ = 0.1. For arcopt's regime (transitioning
  from cubic, B₀ is exact H), either should work. Default to 0.9;
  tune only if benchmarks motivate it.

- **Interaction with box constraints**. arcopt supports `lower`/`upper`
  via `apply_box_constraints`. Line search + box constraints is a
  known-hard combination. For v1, apply box truncation to `α·d` the
  same way cubic path does. For v2, consider projected BFGS (reduced-
  space) which is more principled but substantially more code.

- **Detector false positives / negatives**. The 5-signal rule is a
  heuristic. Calibration will come from benchmark runs. If positives
  cause unnecessary reverts (measured by `qn_polish_reverts` > 0 on
  benign problems), tighten signals. If negatives prevent useful polish
  (measured by "polish never fires on sphere problem"), loosen.

- **Linear equality constraints**. Not yet implemented in arcopt
  anywhere; design/pseudocode.qmd §9.2 has the null-space transformation
  planned. When added, polish mode operates in the reduced z-space
  transparently (it sees only fn/gr in z-space as far as it's concerned).

- **`hess_evals` semantics in return list**. Currently counts every
  `hess()` call. Polish mode only calls `hess()` at init + on
  reverts. A new `hess_evals_saved` field reports `(total_iters_in_polish) -
  (reverts)`, so users can see the win. Alternatively, just expose
  `hess_evals` honestly and let users compute the savings from
  `iterations` − `hess_evals`.

---

## 8. Success criteria

Done when:

1. `devtools::check()` passes 0 errors / 0 warnings / 0 notes.
2. `lintr::lint_package()` shows no new lints on the new files.
3. `devtools::test()` passes all old tests + ~30 new tests for
   line search, healthy-basin detector, and qn_polish integration.
4. Regression: Rosenbrock, mixture saddle, two-exp, GMM, DPM Heckman
   all converge to the same local minimum (or better) as the v1
   cubic-only / cubic+TR run. Iteration count **may be smaller**
   because polish is faster per iter in the basin, but should not be
   larger.
5. Win: on a well-conditioned convex problem (sphere, n ≥ 100), polish
   mode fires within 5–10 iterations and reduces `hess_evals` by ≥ 50%
   compared to `qn_polish_enabled = FALSE`.
6. `solver_mode_final` correctly reports `"qn_polish"` when the run
   terminates in polish mode; `qn_polish_switches` counts each
   cubic→polish transition; `qn_polish_reverts` counts each
   polish→cubic transition.
7. Documentation updated across `NEWS.md`, `R/arcopt.R` roxygen,
   `design/pseudocode.qmd`, `design/design-principles.qmd`.

---

## 9. Not in scope (deferred to later iterations)

- **tr → qn_polish transition**. If the TR fallback reaches a basin, we
  could polish there too. v1 keeps TR as a terminal state.
- **arcopt_qn polish mode**. arcopt_qn already skips `hess()` after
  startup; the polish analog there would be "skip the cubic subproblem
  and just do line search". Useful but orthogonal.
- **Multi-mode telemetry**. Richer diagnostic output (per-mode timing
  breakdown, transition iteration numbers) — nice to have, not
  required for v1.
- **Automatic detector threshold tuning**. Adaptive W, rho_polish, etc.
  based on problem-specific signals. Research-grade; likely premature.

---

**End of briefing.** Start with commit (branch) → Phase 1 (line search)
→ Phase 2 (detector) → Phase 3 (integration) → Phase 4 (bidirectional)
→ Phase 5 (benchmarks) → Phase 6 (docs). Each phase is independently
testable and reviewable.
