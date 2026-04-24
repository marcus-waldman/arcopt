# arcopt simplification briefing — package and manuscript UX

**Purpose**: Self-contained context for a future session focused on
simplifying the user-facing surface of the arcopt package *and* the JSS
manuscript so that the package is actually adoptable. The current state
is feature-complete but over-exposed — too many knobs, too many
first-class concepts, too much internal plumbing in the user-facing
narrative.

**Status**: design + scoping only; no implementation started.

**Parent branch**: `main` (post commits `1b5d928`,
`531acf3` — v0.1.1.9000 with the three-mode hybrid, all tests passing,
manuscript renders cleanly).

---

## 1. Why this exists

The three-mode solver (cubic / TR-fallback / qn_polish) was a sound
engineering exercise, but it has been **shipped as if it were a user
concept**. As of 2026-04-24:

- `control_defaults` in `R/arcopt.R` exposes **~30 control parameters**
  across nine control families (core tolerances, regularization,
  Newton-first, QN updates + method selection, momentum,
  `tr_fallback_*` × 7, `qn_polish_*` × 13).
- The manuscript has **three dedicated subsections** on mode machinery
  (`sec-tr-fallback`, `sec-qn-polish-mode`, `sec-hybrid-practice`) plus
  four full examples (mixture, two-exp, GMM, DPM Heckman) plus a
  six-problem benchmark suite with eight compared methods.
- The return list has **seven dispatch-related diagnostic fields**
  (`solver_mode_final`, `ridge_switches`, `radius_final`,
  `qn_polish_switches`, `qn_polish_reverts`,
  `hess_evals_at_polish_switch`, plus implicit mode in `hessian`).
- Labels are **overloaded**: "arcopt-Hybrid" in the §5 benchmarks means
  `qn_method = "hybrid"` (BFGS↔SR1↔Powell QN updates), while
  "arcopt-hybrid" in the Example 4 Heckman table means
  `tr_fallback_enabled = TRUE` (three-mode solver dispatch). A
  reader landing in §5 cannot reliably guess what `arcopt-hybrid`
  means in the Heckman table.

This is a *product* problem, not a correctness problem. Users whose
only question is "does `arcopt(x0, fn, gr, hess)` converge?" should
not have to read 35 pages and decipher 30 parameters to find out. The
current surface will lose readers in the first five minutes.

---

## 2. Goals

1. **User-facing API**: reduce the documented `control` surface to the
   essential tolerances + one or two explicitly-offered switches.
   Relegate everything else to an advanced / internal help page that
   most users never open.
2. **Return list**: promote `par`, `value`, `gradient`, `converged`,
   `iterations`, `evaluations`, `message` as the *primary* result, and
   collect every mode-dispatch diagnostic into a single `diagnostics`
   sublist that stays out of printed output unless asked for.
3. **Manuscript**: reduce page count and conceptual load. The three
   mode subsections should collapse into a single compact "adaptive
   regularization" subsection; examples and benchmarks should cut
   redundancy.
4. **Defaults**: the defaults should be the *recommendation*, not the
   opt-out baseline. If `qn_polish_enabled = FALSE` is our current
   recommendation, say so and move on; if we plan to flip it, flip it
   and retire the toggle.
5. **Naming**: fix the `arcopt-hybrid` / `arcopt-Hybrid` overload.

Non-goals for this pass: algorithmic behaviour changes, removing
implemented modes, breaking the core `arcopt(x0, fn, gr, hess)`
signature, CRAN submission.

---

## 3. Proposed scope

### 3.1 Package — `R/arcopt.R` and documentation

**Tier 1: promoted (user-visible) control parameters.** Documented at
the top of the `?arcopt` help page with descriptions a statistician
understands without context:

- `maxit`
- `gtol_abs`
- `ftol_abs`, `xtol_abs`
- `trace`
- `use_qn`

That's it for the default help page. Every other knob moves to a
separate, clearly-marked "advanced controls" section (see Tier 2) or
to an internal vignette.

**Tier 2: advanced (hidden-but-documented) controls.** Still
reachable via `control = list(...)`, still tested, but documented in
a separate `?arcopt_advanced_controls` topic (or equivalent) rather
than in the main help page. These include:

- `sigma0`, `sigma_min`, `sigma_max`, `eta1`, `eta2`, `gamma1`, `gamma2`
  (regularization tuning)
- `qn_method` ("sr1" / "bfgs" / "hybrid")
- All `tr_fallback_*` and `qn_polish_*` thresholds
- `use_momentum` (if we keep it — see 3.2)

**Tier 3: internal-only.** Remove from the documented control surface
entirely and/or rename with a `.` prefix to signal experimental:

- Anything added only for unit tests
- Anything the implementation relies on but users should never set

**Return list trimming.** Keep the primary fields unchanged. Collect
every mode-dispatch diagnostic into one sublist:

```r
result$diagnostics$solver_mode_final
result$diagnostics$ridge_switches
result$diagnostics$radius_final
result$diagnostics$qn_polish_switches
result$diagnostics$qn_polish_reverts
result$diagnostics$hess_evals_at_polish_switch
```

so that `str(result)` at the REPL shows a clean primary object and the
diagnostics stay out of the way unless explicitly inspected. This is a
**breaking change** for downstream code that reads those fields at the
top level, but the only such code in-tree is the cache producer and
the integration tests — cheap to update.

### 3.2 Decisions to pin before implementation

1. **Flip `qn_polish_enabled = TRUE` or keep opt-in?** The integration
   tests pass at defaults; the only reason it is opt-in today is
   briefing §9's "pending broader benchmark evidence." If we flip it,
   the control parameter can disappear from the help entirely, and
   one whole subsection of the manuscript goes away. If we keep it
   opt-in, we should commit to a calendar milestone for flipping
   (e.g., after the JSS revision cycle).
2. **Remove or retain `use_momentum`?** The manuscript mentions it in
   one paragraph (lines 461–464). If it has empirical value, elevate
   it; if not, remove it and any tests. No middle ground — having it
   present but undocumented is the worst UX.
3. **Keep `arcopt_qn` as a separate entry point, or fold into
   `arcopt()` with `use_qn = TRUE` only?** The current situation is
   *both*: `arcopt()` accepts `use_qn = TRUE` and there is a separate
   `arcopt_qn()` function. Pick one.
4. **Do we want a `verbose = TRUE` flag** that prints a one-line
   per-iteration summary at the R console when enabled? Current `trace
   = 1` output is in C++ and less useful. A Python-style progress
   indicator would improve UX substantially.

### 3.3 Manuscript — `manuscript/arcopt-jss.qmd`

**Collapse the three mode subsections.** Merge `sec-tr-fallback`
(lines 297–349), `sec-qn-polish-mode` (lines 351–425), and
`sec-hybrid-practice` (lines 1929–2019) into one compact subsection
titled roughly "Adaptive regularization and fallback modes" that:

- Opens by framing cubic regularization as the default.
- Names the two fallbacks (trust-region fallback for flat-ridge
  stagnation; quasi-Newton polish for the quadratic basin) in one
  paragraph each.
- Defers all detector thresholds, signal counts, and state-machine
  details to a single appendix or to the package documentation.
- Target length: 30 lines of prose, no algorithm boxes.

**Cut one example.** The mixture saddle (sec-mixture-example), the
two-exponential saddle (sec-two-exp-example), and the growth mixture
model (sec-gmm-example) all make the *same* point: QN without the
true Hessian escapes symmetric saddles poorly, and the full-Hessian
**arcopt** succeeds. Two of the three is enough. Propose cutting
**two-exp** because (a) it is the most synthetic, (b) it shares the
"identifiability saddle" framing with GMM, and (c) GMM better
motivates the FD-init contribution. Keep mixture (pedagogical),
GMM (statistically realistic saddle), and DPM Heckman (flat-ridge
motivation).

**Tighten the FD-init and QN-hybrid routing prose.** The
"state-aware routing scheme" passage (roughly lines 518–532) spends
15 lines on a three-mode FSM inside the QN path. Compress to a
single paragraph that names the BFGS → SR1 → Powell fallback order
and points to `?qn_method` for details.

**Re-label Example 4's table** to disambiguate from §5's
`arcopt-Hybrid`. Replace `arcopt-nohybrid` / `arcopt-hybrid` with
something like `arcopt (cubic only, v0.1.0)` and
`arcopt (v0.1.1 default)` — or, if the §5 `arcopt-Hybrid` column is
renamed to `arcopt-QN-hybrid`, the Example 4 table can keep
`arcopt-hybrid` without ambiguity.

**Cut `sec-hybrid-practice`'s polish demo** if we're simplifying.
The 33%/56% hess-eval savings on a smooth convex quartic-tail are a
nice flourish but they are buying us zero new readers — a
statistician reading the paper already trusts that BFGS polishes
faster than a second-order method in the basin. Move the empirical
observation to a one-sentence remark inside the collapsed mode
subsection; drop the demo chunk.

**Review the benchmark prose.** Eight methods × six problems is a
lot of ground. Consider cutting SR1 from the benchmark method set
(it is strictly worse than Hybrid in every row) to reduce the
table's column count.

### 3.4 Rough page count targets

Today (v0.1.1, post-Example 4 + hybrid-practice): ~35 printed pages
in the apaquarto-pdf layout.

Target after simplification: ~28 pages. Cuts:

- Merged mode subsection: −3 pages
- Two-exp example cut: −3 pages
- Hybrid-practice polish demo cut: −2 pages
- QN routing prose tighten: −1 page
- Net: ~9 pages saved, yielding ~26–28 pages.

### 3.5 NEWS.md entry (to draft)

```
## v0.2.0 — user experience overhaul

- Renamed mode-dispatch diagnostics into `result$diagnostics` sublist
  (breaking change for downstream readers of `result$solver_mode_final`
  etc.; update to `result$diagnostics$solver_mode_final`).
- Promoted {maxit, gtol_abs, ftol_abs, xtol_abs, trace, use_qn} as the
  user-facing control surface; all other parameters documented under
  `?arcopt_advanced_controls`.
- [optional] Flipped `qn_polish_enabled` default to TRUE based on
  benchmark evidence.
- [optional] Removed / promoted `use_momentum`.
```

---

## 4. Critical files

- **Edit**: `R/arcopt.R` — tier the `control_defaults`; restructure
  roxygen into main + advanced; nest `diagnostics` sublist in return
  list.
- **Edit**: `R/arcopt_qn.R` — same treatment if we keep the separate
  entry point; decide whether to remove.
- **Edit**: `manuscript/arcopt-jss.qmd` — merge mode subsections, cut
  two-exp, cut polish demo, disambiguate hybrid labels, tighten QN
  prose.
- **Edit**: `NEWS.md` — v0.2.0 entry.
- **Edit**: `DESCRIPTION` — bump to `0.2.0.9000`.
- **Edit**: integration tests and cache scripts that read the
  top-level diagnostic fields (nest them).
- **New**: `man/arcopt_advanced_controls.Rd` (or equivalent roxygen
  source) for the tier-2 controls page.

---

## 5. Out of scope

- Any algorithmic change in how cubic / TR / polish modes actually
  behave.
- Removing any currently-implemented detector or mode.
- Breaking the `arcopt(x0, fn, gr, hess, ...)` core signature.
- CRAN submission or compliance work beyond what the simplification
  touches.
- Manuscript figures and typography beyond the cuts listed above.
- The v2-future items from the original briefings (tr → qn_polish
  transition; polish in `arcopt_qn`).

---

## 6. Starting instructions for a fresh session

1. **Read this briefing top to bottom.**
2. **Read the two implementation briefings**
   `todo/tr-fallback-hybrid-briefing.md` and
   `todo/qn-polish-mode-briefing.md` for the algorithmic context —
   this briefing is *only* about UX; algorithms are unchanged.
3. **Read the current `control_defaults` in `R/arcopt.R`** and the
   roxygen header to see the surface we are pruning.
4. **Read `manuscript/arcopt-jss.qmd`** sections 2, 3, 4, 5.3 to see
   the current shape.
5. **Propose the tier-1/2/3 partition** as the first concrete step and
   get the user's sign-off before editing `R/arcopt.R`. The
   partitioning decision is the fulcrum of the whole simplification;
   do not barrel through it.
6. **Decide the four pin-downs in §3.2** with the user before code
   changes. These are product decisions, not implementation details.
7. **Then edit in the order**: `R/arcopt.R` → `NEWS.md` →
   `manuscript/arcopt-jss.qmd` → integration tests / cache scripts →
   render manuscript → commit and push.

---

## 7. Open questions / risks

- **Breaking-change magnitude.** Nesting diagnostics under
  `result$diagnostics$*` breaks any external code reading the fields
  at the top level. No such code is known to exist outside this repo,
  but we should search before committing.
- **Losing the mode-transparency message.** Today's manuscript
  explicitly makes the hybrid's structure visible ("the three-mode
  dispatch..."). If we collapse that into a compact paragraph, we
  trade narrative clarity for reader bandwidth. That trade is
  deliberate, but a reviewer might push back asking "what do these
  other modes do?" Keep the algorithm box / appendix available so
  the answer is one-page-flip away.
- **Cutting two-exp** removes a specific pedagogical point about
  multi-dimensional saddles in nonlinear regression. Confirm the GMM
  example covers that territory well enough.
- **`qn_polish_enabled` flip**. If we flip the default to `TRUE`, we
  commit to the polish mode's correctness on *all* downstream
  statistical problems, not just the current integration test suite.
  Before flipping, run polish with defaults on every manuscript
  example and every benchmark problem and confirm zero regressions
  (zero polish_reverts, zero divergences from the cubic-only answer).
  The probe script at
  `manuscript/todo/probe_hybrid_on_all_problems.R` is the starting
  point for this check.

---

## 8. Success criteria

Done when:

1. `?arcopt` help page fits on one screen (≤60 lines of useful prose
   + the args block) and documents only the tier-1 parameters in the
   main `Arguments` section.
2. `str(arcopt(x0, fn, gr, hess))` on a simple problem shows a
   primary return with a single tidy `diagnostics` sublist.
3. Manuscript page count drops by ≥5 vs. the 2026-04-24 state.
4. The `arcopt-hybrid` / `arcopt-Hybrid` label collision is gone from
   the manuscript.
5. `devtools::check()` passes with 0 errors / 0 warnings / 0 notes.
6. `devtools::test()` passes after diagnostic-path nesting is
   applied.
7. `benchmarks/dpm_heckman/cache_for_manuscript.R` re-runs
   successfully and the manuscript still renders Example 4 with the
   nested diagnostics.
8. A user who has never read the ARC literature can open `?arcopt`,
   copy the first example, run it on their problem, and understand
   the output without reading anything else.

---

## 9. Context for "where we left off" on 2026-04-24

- **Code**: `main` at commit `1b5d928`. v0.1.1.9000 with three-mode
  solver, all tests passing, manuscript renders.
- **Manuscript**: `manuscript/arcopt-jss.qmd` is gitignored
  (local-only). Covers cubic + TR-fallback + qn_polish, with
  Example 4 (DPM Heckman MAP, 299 params, cached at
  `manuscript/_cache/dpm_heckman_results.rds`). PDF at 308 KB.
- **Cache producer**: `benchmarks/dpm_heckman/cache_for_manuscript.R`
  tracked on `main`; full run ~11 minutes.
- **Probe**: `manuscript/todo/probe_hybrid_on_all_problems.R`
  (gitignored) runs every manuscript problem with legacy / default /
  full-hybrid variants; useful as the regression baseline for any
  default change.
- **Installed arcopt**: 0.1.1.9000 in the user's R 4.5.1 library, so
  `library(arcopt)` in the manuscript picks up the new return
  fields. Any `diagnostics` nesting will require a re-install
  (`devtools::install(quick = TRUE)`) before the manuscript renders.
- **Outstanding user comment** that motivated this briefing: the v1
  manuscript exposes too much internal machinery to be a readable
  package paper, and the control surface will deter adoption. This
  briefing is the response.
