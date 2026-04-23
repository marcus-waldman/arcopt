# TR-Fallback Hybrid for arcopt — Implementation Briefing

**Purpose of this document**: Self-contained context to start a fresh session and implement the cubic+trust-region hybrid in arcopt. Covers what prompted this work, the design, phases, and verification plan.

**Status**: design agreed; implementation not started. JSS deadline set aside — do this right, not fast.

---

## 1. Why this exists

On 2026-04-21 we built a DPM mixture Heckman MAP benchmark (`benchmarks/dpm_heckman/`) expecting it to showcase arcopt's saddle-escape advantage. Instead, arcopt **lost decisively** to `trust`, `optim-BFGS`, and `nlminb`:

| Method | f | \|g\|_inf | iter | time | eig_min | status |
|---|---|---|---|---|---|---|
| **trust** | 667.0 | **1.5e-9** | **97** | 164s | +0.009 | local min ✓ |
| optim-BFGS | **616.5** | 3.0e-5 | 407 | 3.5s | +0.027 | local min |
| nlminb | 662.5 | 2.9e-3 | 271 | 1.6s | +0.018 | local min |
| arcopt-AD | 683.7 | 1.61 | 500 | 840s | −276 | saddle (stuck) |
| arcopt-hybrid | 684.0 | 3.65 | 500 | 22s | −0.32 | borderline |
| arcopt-sr1 | 985.3 | 301 | 500 | 17.8s | −120 | saddle |
| arcopt-bfgs | 716.6 | 4.1 | 500 | 17.3s | −3.1 | saddle |
| arcopt-FD | 704.2 | 34.6 | 500 | 3202s | −9440 | saddle |

n=500 observations, K_max=20, K_true=2, dim(X)=5, dim(W)=7, 299 parameters.

Critically, `trust` uses the **same analytic Hessian** as arcopt-AD but converges in 97 iterations while arcopt-AD runs 500 without converging. The landscape isn't stumping Hessian-based methods in general — arcopt's cubic regularization specifically gets trapped.

A σ trace (see `benchmarks/dpm_heckman/trace_sigma.R`) showed arcopt reaching late iterations with:
- σ pinned at `σ_min = 1e-6`
- ρ_k ≈ 1.0 consistently (model predictions perfect)
- ‖g‖_inf stuck at ~1.6 (no meaningful decrease)
- Hessian has **positive eigenvalues** (no saddle at end state)

This is a **flat-ridge** signature, not a saddle. arcopt's cubic regularization `(σ/3)‖s‖³` is tuned to dampen steps in high-curvature / negative-curvature regimes — it has no advantage when the Hessian is near-singular-PD along a ridge. The cubic term vanishes (σ hits floor) and what remains is essentially Newton's method on an ill-conditioned quadratic, which takes many iterations to walk off the ridge.

`trust`, using the same eigendecomposition but with a **hard step-norm constraint** instead of the cubic penalty, takes more aggressive Newton-like moves along the ridge and converges fast.

---

## 2. The idea — hybrid cubic + trust-region

### σ ↔ r duality

Both solvers use the same eigendecomposition of H_k and solve for a step `s`:

- **Cubic (ARC)**: `min g'·s + (1/2) s'Hs + (σ/3)‖s‖³`. First-order optimality gives `(H + λI) s = −g`, `λ = σ‖s‖`.
- **Trust-region**: `min g'·s + (1/2) s'Hs` subject to `‖s‖ ≤ r`. First-order optimality gives `(H + λI) s = −g`, `λ ≥ 0`, with `‖s‖ = r` when constraint is active.

Both reduce to the same eigendecomposition machinery (rotate to eigenbasis of H; reduce to 1-D root-finding for λ). The secular equations differ only in the constraint: `λ = σ‖s(λ)‖` vs `‖s(λ)‖ = r`.

### Regime appropriateness

- **Indefinite / highly non-convex**: cubic wins. `(σ/3)‖s‖³` naturally dampens steps when λ_min(H) < 0; the step along the negative-curvature eigendirection is bounded.
- **Flat-ridge / near-singular-PD**: TR wins. The step-norm constraint keeps us within a controlled region, but allows the maximum ‖s‖ = r in the small-eigenvalue direction — which is exactly what you want to walk off a ridge in one step.

Neither is universal. A method that detects which regime it's in and dispatches accordingly gets the best of both.

### Precedent

- Cartis-Gould-Toint's original ARC papers discuss cubic regularization's complexity advantage on nonconvex problems but are silent on ridge degeneracy.
- Trust-region methods (Conn-Gould-Toint book) are well-studied but struggle at saddles without special modifications (e.g., Moré-Sorensen's "hard case").
- We're not aware of prior work on adaptive selection between cubic and TR subproblems based on local landscape diagnostics. If it exists, cite it; if not, this is a real contribution.

---

## 3. Design — agreed decisions

1. **One-way switch (cubic → TR, no going back)** in v1. Rationale: simpler, matches the observation that once we detect ridge behavior, staying in TR is robust. Revisit bidirectional if a real problem demands it.

2. **Default on (`tr_fallback_enabled = TRUE`)**. This is new behavior that genuinely helps; users should get it by default. Keep an opt-out for research reproducibility.

3. **Detect on sliding window, all four signals required**:
   - σ pinned at floor: `σ_k ≤ 10·σ_min` for W consecutive iterations
   - Near-perfect model: `|ρ_k − 1| < 0.1` for W iterations
   - Gradient stagnant: `‖g‖_∞(k) / ‖g‖_∞(k−W) > 0.9` (less than 10% decrease over window)
   - Hessian positive-definite but poorly conditioned: `0 < λ_min(H_k) < tol_ridge`

   Default W=10, tol_ridge=1e-3. All tunable.

4. **TR subproblem solver**: port `cubic_eigen.R` structure; change secular equation. Reuse eigendecomposition cache if possible (the Hessian hasn't changed between subproblem solves within one iteration).

5. **Initial TR radius at switch**: `r_0 ← max(last_step_norm, 1/σ_current)`. Then standard TR adaptation.

6. **Same ρ ratio and acceptance machinery as cubic**. Just swap the step computation and the σ ↔ r update rules.

---

## 4. Critical files

### Existing (read to understand)

- `C:\Users\marcu\git-repositories\arcopt\R\arcopt.R` — main orchestration. Iteration loop calls `solve_newton` then `solve_cubic_eigen`, with σ update via `update_sigma_cgt`. Will need a mode state (`"cubic"` / `"tr"`) and dispatch.
- `C:\Users\marcu\git-repositories\arcopt\R\arcopt_qn.R` — QN variant. Same structure; needs the same mode-state addition.
- `C:\Users\marcu\git-repositories\arcopt\R\cubic_eigen.R` — eigendecomposition solver for cubic subproblem. Reference architecture for the TR solver.
- `C:\Users\marcu\git-repositories\arcopt\R\sigma_update.R` — σ adaptation (η1/η2/γ1/γ2 logic). The TR radius update mirrors this; consider a shared helper.
- `C:\Users\marcu\git-repositories\arcopt\R\convergence.R` — convergence check. Probably unchanged; both solvers feed into the same criteria.
- `C:\Users\marcu\git-repositories\arcopt\design\pseudocode.qmd` — algorithmic specifications. Will need an Algorithm for the TR subproblem + Algorithm for the flat-ridge detector.

### Benchmarks / evidence

- `C:\Users\marcu\git-repositories\arcopt\benchmarks\dpm_heckman\` — motivating failure case. **Not yet committed** as of briefing time. First commit this as-is (prerequisite to the work).
  - `setup.R`, `dpm_heckman.stan`, `benchmark.R`, `trace_sigma.R`, `test.R`, `results/`
- `C:\Users\marcu\git-repositories\arcopt\benchmarks\matrix_factorization\`, `deep_linear\`, `chebyshev_rosenbrock\`, `brown_almost_linear\`, `irt_3pl\` — scope-boundary benchmarks committed on 2026-04-21 (commit `35b637c`). These must not regress with the new hybrid.

### External reference

- `trust` package source on GitHub: `https://github.com/cjgeyer/trust/blob/master/package/trust/R/trust.R`. Brief summary: uses eigendecomposition (not dogleg or Steihaug-CG), solves secular equation via `uniroot`, handles easy/hard-hard cases, updates radius on rho thresholds 1/4 and 3/4 with factor-4 shrink / factor-2 grow up to `rmax`. Its subproblem logic is close to a drop-in template for our `solve_tr_eigen`.

---

## 5. Implementation phases

Estimate: 3–5 working days for core + benchmarks. Additional time for manuscript.

### Phase 1 — TR subproblem solver

**New file**: `R/tr_eigen.R`

```r
# solve_tr_eigen(g, H, radius) -> list(s, pred_reduction, lambda, on_boundary)
# - eigendecomposition of H (symmetric)
# - if Newton step ||H^-1 g|| <= radius AND H is PD: return Newton step
# - else: uniroot for lambda in (max(0, -lambda_min(H)), large) such that
#   ||s(lambda)|| = radius, where s(lambda) = -(H + lambda I)^-1 g
# - compute pred_reduction = -g'·s - (1/2) s'Hs
# - handle hard case (lambda = -lambda_min and direction aligns with
#   eigenvector of lambda_min) per Moré-Sorensen
```

~80 LOC. Mirrors `cubic_eigen.R` tightly; differs only in the secular equation.

**New tests**: `tests/testthat/test-tr-eigen.R`
- Trivial PD quadratic: Newton step inside radius, returned unchanged
- Newton step outside radius: boundary solution
- Indefinite Hessian: lambda > |lambda_min|, step on boundary
- Hard case: Moré-Sorensen corner, step includes eigenvector component
- Agreement with `trust::trust()` single-step solve on a few configs (loose tolerance; different step norm conventions)

### Phase 2 — Flat-ridge detector

**New file**: `R/flat_ridge.R`

```r
# init_flat_ridge_state(window = 10) -> state
# update_flat_ridge_state(state, sigma, rho, g_inf, H_eig_min, sigma_min, tol_ridge) -> state
#   pushes signals onto sliding window of length `window`
# check_flat_ridge_trigger(state, window, tol_ridge) -> logical
#   returns TRUE if ALL four conditions hold throughout the window
```

~80 LOC including state management.

**New tests**: `tests/testthat/test-flat-ridge.R`
- State accumulates correctly over iterations
- Trigger fires only when all four signals align for full window
- Trigger does not fire at a normal local minimum (gradient small)
- Trigger does not fire at a saddle (lambda_min negative)
- Trigger does not fire when sigma is still oscillating

### Phase 3 — Orchestration integration

**Edit**: `R/arcopt.R`, `R/arcopt_qn.R`

Add to iteration state:
- `solver_mode` in `{"cubic", "tr"}`, starts at `"cubic"`
- `ridge_state` from `init_flat_ridge_state()`
- `radius_current` — initialized at switch time

Per-iteration logic:
```
if (solver_mode == "cubic") {
  update ridge_state with current (sigma, rho, g_inf, lambda_min)
  step <- solve_cubic_eigen(g, H, sigma)
  ... (existing)
  if (check_flat_ridge_trigger(ridge_state, ...)) {
    solver_mode <- "tr"
    radius_current <- max(||step||, 1/sigma_current)
    ridge_switches <- ridge_switches + 1  # for diagnostics
  }
} else {  # "tr" mode
  step <- solve_tr_eigen(g, H, radius_current)
  ... compute rho ...
  radius_current <- update_radius_tr(radius_current, rho, step_on_boundary, ...)
}
```

Both modes feed into the same convergence check and Hessian evaluation. σ continues to be updated in cubic mode; radius in TR mode.

**New control parameters** (all tunable, with conservative defaults):
- `tr_fallback_enabled` (default `TRUE`)
- `tr_fallback_window` (default `10`)
- `tr_fallback_tol_ridge` (default `1e-3`)
- `tr_rmax` (default `1e6`)
- `tr_eta1` (default `0.25`) — TR acceptance threshold
- `tr_eta2` (default `0.75`) — TR expansion threshold
- `tr_gamma_shrink` (default `0.25`) — radius shrink factor on bad step
- `tr_gamma_grow` (default `2.0`) — radius grow factor on good boundary step

**New return fields**:
- `solver_mode_final` — `"cubic"` or `"tr"`
- `ridge_switches` — integer count of cubic→tr transitions (0 or 1 in v1)
- `radius_history` / `sigma_history` — for trace mode

### Phase 4 — Benchmark validation

**Motivating win**:
- `benchmarks/dpm_heckman/benchmark.R` — rerun with hybrid enabled, confirm arcopt reaches trust's local min (or better) within similar iter budget.

**Regression checks** (must not worsen materially):
- Existing mixture saddle example (manuscript §Example 1) — arcopt-hybrid should still hit 100% saddle escape.
- Two-exponential saddle (§Example 2)
- GMM saddle (§Example 3)
- `benchmarks/matrix_factorization/`, `deep_linear/`, `chebyshev_rosenbrock/`, `brown_almost_linear/`, `irt_3pl/` — rerun and compare.

**New synthetic flat-ridge test problem** (not ARC-standard):
- Simple 2D quadratic with eigenvalues (1, 1e-4): test that detection fires and TR converges. Add to `benchmarks/` or `tests/testthat/`.

### Phase 5 — Documentation

- `R/arcopt.R` roxygen: new control params, solver modes, return fields
- `NEWS.md`: 0.2.0 section introducing the hybrid
- `README.md`: brief mention of the fallback capability
- `design/pseudocode.qmd`:
  - New Algorithm 5c (TR subproblem via eigendecomposition) alongside existing 5a
  - New Algorithm (flat-ridge detector)
  - Updated Algorithm 1 (main loop) showing mode dispatch
- `design/design-principles.qmd`: add a paragraph on the σ ↔ r duality and why a hybrid is principled, not ad hoc

### Phase 6 — Manuscript (JSS or follow-up)

With the deadline lifted, there's no pressure to cram this into v1. Options:
- **Expanded v1**: add new §3.x subsection on method selection; use DPM MAP as Example 4; extend §2 with duality framing.
- **Follow-up paper**: submit current 3-example manuscript; write "Adaptive solver selection for cubic-regularized optimization" as a second paper.

Deferred. Decide after implementation and benchmarking.

---

## 6. Starting instructions for a new session

1. **Read this briefing file first**.

2. **Read the plan file** at `C:\Users\marcu\.claude\plans\review-the-manuscript-and-happy-bumblebee.md` for context on the DPM Heckman work that motivated this.

3. **Read key memories** under `C:\Users\marcu\.claude\projects\C--Users-marcu-git-repositories-arcopt\memory\`:
   - `scaling_boundary_evidence_2026-04-21.md` — scope-boundary benchmarks context
   - `fd_init_qn_default.md` — current arcopt-hybrid routing state
   - `manuscript_state_2026-04-20.md` — manuscript structure

4. **First action — commit the DPM Heckman benchmark** as-is. This is the motivating evidence and should be preserved before any refactor:
   ```
   git add benchmarks/dpm_heckman
   git commit -m "Add DPM mixture Heckman benchmark: motivating failure for TR-fallback hybrid

   n=299 params MAP. arcopt variants all end at saddles; trust converges cleanly.
   Motivates implementing a cubic->TR hybrid (see todo/tr-fallback-hybrid-briefing.md)."
   ```

5. **Second action — create a branch**:
   ```
   git checkout -b feature/tr-fallback-hybrid
   ```

6. **Then start Phase 1**: implement `R/tr_eigen.R` and its tests.

7. **Verify after each phase** with `devtools::test()` + `lintr::lint_package()` + the relevant benchmark(s).

---

## 7. Open questions / risks

- **Numerical stability of secular-equation root-finding** when Hessian eigenvalues span many orders of magnitude. `uniroot` should cope; if not, consider `nleqslv` or a hand-rolled Newton-on-the-secular-equation per Moré-Sorensen.
- **Eigendecomposition caching**: the Hessian hasn't changed between cubic and TR subproblem solves within one iteration (they're both invoked on H_k). Consider computing the eigendecomposition once per outer iteration and passing it to whichever subproblem solver fires. Small refactor of `cubic_eigen.R` to accept a pre-computed decomposition.
- **Detection false positives**: the detector should not fire at a true local minimum where ‖g‖ is just small. The "gradient stagnant" condition is `‖g‖_∞(k) / ‖g‖_∞(k−W) > 0.9`, which is satisfied at a local min where both are near zero — need an absolute floor, e.g., also require `‖g‖_∞ > some threshold` to avoid firing at convergence.
- **Interaction with arcopt_qn's FD refresh triggers**: the FD-refresh logic in the QN path already handles "stuck" detection differently. Need to ensure the two don't conflict — possibly only one of (QN refresh, TR fallback) fires per regime.
- **`trust` as a dependency**: we reuse its algorithm, not its code. Native implementation means `trust` stays in Suggests, not Imports. If we can't get native to work, drop-in calls to `trust::trust()` would upgrade it to Imports.
- **One-way switch regret**: if a problem has ridge→indefinite→ridge structure, one-way switch loses saddle-escape after the first switch. Monitor via `ridge_switches` count; if this happens in practice, generalize to bidirectional in v2.

---

## 8. Success criteria

The feature is done when:

1. `devtools::check()` passes 0 errors / 0 warnings / 0 notes with the new code.
2. `lintr::lint_package()` is clean on all new/changed files.
3. `devtools::test()` passes all old tests + new tests for TR subproblem and flat-ridge detection.
4. DPM mixture Heckman benchmark: arcopt (any variant) reaches a true local minimum (eig_min > 0) on at least one seed in ≤ 2× trust's iteration count.
5. All existing manuscript examples (mixture, GMM, two-exp) are unchanged in outcome — hybrid doesn't fire on saddle problems.
6. New synthetic flat-ridge test problem: TR fallback triggers and converges to `‖g‖_inf < 1e-8`.
7. Documentation updated across `NEWS.md`, `R/arcopt.R` roxygen, `design/pseudocode.qmd`.

---

**End of briefing.** Reach out with any questions that weren't anticipated here; otherwise start with commit → branch → Phase 1.

---

## Appendix A — Yue, Zhou, So (2018) review (added 2026-04-23)

Reference: *On the Quadratic Convergence of the Cubic Regularization Method under a Local Error Bound Condition*, SIAM J. Optim. 29(1):904–932.

**What they prove.** Under a local error bound (EB) condition `dist(x, X) ≤ κ ||∇f(x)||` in a neighborhood of the second-order critical set `X` (their Definition 1 / Assumption 2), CR attains Q-quadratic convergence of iterates `{x^k}` to a second-order critical point — even when minimizers are *non-isolated* (degenerate). Their Theorem 1 establishes that EB is equivalent to quadratic growth `f(x) ≥ f(x̂) + (α/2) dist²(x, X)` under a mild separation-of-isocosts assumption.

**What they don't cover.** The EB condition fails (or becomes vacuous with κ → ∞) when the landscape is flatter than quadratic-growth — e.g., a stretched ridge where gradient shrinks faster than distance to the true minimum. Yue-So's own example `f(x) = (||x||² − 1)²` satisfies EB (X is the unit sphere), because the gradient `4(||x||² − 1)x` is small only when close to X. A ridge with O(‖x‖⁴) decay along the ridge direction would violate EB: gradient vanishes along the ridge but distance to the global minimum does not.

**Why this sharpens the hybrid contribution.**
- Our DPM Heckman stall is in a regime **weaker than EB** (or equivalently, weaker than Łojasiewicz-1/2 / gradient-dominance, Yue-So Definition 4).
- No published convergence result for CR/ARC covers this regime.
- The four-signal detector (σ floor, ρ≈1, ‖g‖ stagnant, 0 < λ_min(H) < tol_ridge) is an **empirical proxy for EB violation**: when all four fire, the effective κ in `dist ≤ κ ||∇f||` would have to be very large (≥ 1/tol_ridge) — EB is practically vacuous.
- The hybrid is therefore not a workaround for a theoretical gap that Yue-So plugged — it's a solution for landscapes beyond Yue-So's assumption.

**Design confirmations.**
1. The σ floor (`σ_min = 1e-6`) being pinned is equivalent to saying the cubic model trusts local quadratic structure — so EB-dominated convergence should already be in effect. If we're not converging despite that, we're outside EB.
2. `tol_ridge = 1e-3` default is consistent: EB with κ ~ 10³ is practically useless as a convergence guarantee, so using λ_min(H) < 1e-3 as "ridge" is calibrated to where Yue-So's guarantee no longer has teeth.
3. Gradient-stagnation check (‖g‖_k / ‖g‖_{k−W} > 0.9 over window W) directly tests local gradient-dominance violation.

**Citations to add.**
- Manuscript §Related Work / §Background: Yue-So (2018) as the strongest known CR convergence result for degenerate minima, and the paragraph explaining our ridge regime lies *outside* their assumption.
- Detector description in §Algorithm: frame the four signals as EB-violation diagnostics.
- Still open: whether EB is testable in practice (no — κ is not computable without knowing X). Our detector substitutes observable runtime signals.

**No changes needed to briefing body** — the detector design is validated. Proceed to Phase 1.
