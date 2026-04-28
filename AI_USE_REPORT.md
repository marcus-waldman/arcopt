# AI use report — arcopt

This document reports the use of Anthropic's Claude (Sonnet and Opus
4.x model families) in the development of the **arcopt** R package
and its accompanying JSS manuscript. It is intended to satisfy the
disclosure expectations of CRAN reviewers, JSS reviewers, and any
downstream user who wants to understand which components of the
package were AI-assisted and how each was verified.

The package is authored and maintained by Marcus Waldman
(University of Colorado Anschutz). All design decisions, framing
choices, scope decisions, and final acceptances of code and prose are
the author's. Claude's role is described component-by-component
below.

---

## 1. Scope of disclosure

This report covers the v0.3.0 release of arcopt (first CRAN
submission, dated 2026-04-27). Earlier releases (v0.1.x, v0.2.0)
involved Claude assistance to varying degrees across iterative
sessions over the preceding months. The discussion below describes
the *state of v0.3.0* rather than attempting to reconstruct
session-by-session contributions, which were not formally tracked.

The disclosure applies to:

* The R package source tree at the root of the repository
  (`R/`, `src/`, `tests/`, `man/`, `vignettes/`, `NAMESPACE`,
  `DESCRIPTION`, `NEWS.md`, `README.md`, `_pkgdown.yml`,
  `.Rbuildignore`, `cran-comments.md`).
* The JSS manuscript at `manuscript/arcopt-jss.qmd` and its
  bibliography at `manuscript/references.bib`.
* The benchmark and validation scripts under `benchmarks/`.

It does *not* cover external dependencies (Rcpp, RcppEigen, trust,
marqLevAlg, etc.), the algorithm specifications themselves
(Cartis-Gould-Toint ARC, Kamzolov-class quasi-Newton ARC, Yue-Zhou-So
local error-bound theory, Liu-Rubin / Andersen psychometric
references), or the underlying R language and ecosystem.

---

## 2. Components and their AI involvement

The components below are listed roughly in order of the *strength* of
verification applied. Components near the top have automated tests
that catch regressions; components near the bottom rely more on
author review.

### 2.1 Algorithm implementation (`R/`, `src/`)

The cubic-regularization solver, Newton-first dispatcher, sigma
adaptation rules, indefiniteness handling, BFGS / SR1 / hybrid
quasi-Newton updates, the trust-region fallback's four-signal
detector, and the quasi-Newton polish mode's five-signal detector
were all implemented in R with Claude assistance over multiple
sessions. The C++ surface in `src/` is currently a single
placeholder; the algorithmic work is in pure R.

**Algorithmic sources (not Claude):** Cartis, Gould & Toint (2011a/b);
Nesterov & Polyak (2006); Kamzolov et al. (2023, 2025); Yue, Zhou &
So (2019); Cauchy / Newton step composition theory.

**Verification:**
* 594 unit tests in `tests/testthat/` with 0 failures, 0 warnings,
  0 skips on R 4.5.1 (Windows). The tests cover: cubic-subproblem
  correctness against eigendecomposition reference solutions, BFGS /
  SR1 / hybrid update rules under PD and indefinite states, both
  detectors' state-machine transitions, box-constraint projection,
  convergence criteria, and end-to-end integration on Rosenbrock,
  saddle-escape, and flat-ridge problems.
* `R CMD check --as-cran` on the local system: 0/0/0 modulo a known
  network-environment NOTE.
* Multi-platform CI via GitHub Actions:
  `.github/workflows/R-CMD-check.yaml` runs the test suite on
  macOS-release, Windows-release, Ubuntu-devel, Ubuntu-release, and
  Ubuntu-oldrel-1.

### 2.2 Tests (`tests/testthat/`)

Most tests were drafted with Claude assistance, then reviewed and
adjusted by the author. The test suite is the package's primary
verification mechanism: any regression introduced by AI-generated
code is caught here.

**Verification:**
* The tests test the algorithm code; the algorithm code's outputs are
  cross-checked against analytic results where applicable
  (closed-form Hessians, finite-difference Jacobians via
  `benchmarks/utils.R`, hand-verified eigendecompositions on small
  examples).
* No tests are skipped on CRAN (`tests/testthat/` contains no
  `skip_on_cran()` calls), so all 594 tests run on the CRAN
  pre-flight machines.

### 2.3 Benchmark and validation scripts (`benchmarks/`)

The setup files (problem definitions, simulators, fn / gr / hess
wrappers) and the manuscript-cache scripts were drafted with Claude
assistance. The numbers they produce flow into the JSS manuscript via
`manuscript/_cache/*.rds` files and inline `r ...` expressions in the
qmd, so any silent error in these scripts would surface as a
manuscript-table inconsistency. The directory is excluded from the
CRAN tarball via `.Rbuildignore`.

**Verification:**
* `benchmarks/utils.R` provides `check_gradient()` and
  `check_hessian()` finite-difference utilities. New problem setups
  are validated against these before being used in benchmarks. For
  example, the sparse-Rasch JML setup in
  `benchmarks/rasch_jml/setup.R` passes FD gradient and Hessian
  checks at <1e-6 maximum absolute error (test transcript:
  `benchmarks/rasch_jml/test.R`).
* The §4.4 manuscript example uses 100 random datasets cached via
  `benchmarks/rasch_jml/cache_for_manuscript.R`; the qmd loads the
  cache and the displayed numbers are computed at render time, not
  hardcoded. This catches drift between the package and the
  manuscript.
* A separate Phase-1 validation gate caught a mistaken initial choice
  of headline example. The original §4.4 candidate (multivariate-t
  MLE) was implemented and tested, and the polish-saving go/no-go
  criterion failed at default thresholds. This was reported, the
  example was replaced with sparse-Rasch JML, and only then was
  prose written. Files retained: `benchmarks/multivariate_t/`.

### 2.4 Vignettes (`vignettes/`)

Both vignettes (`getting-started.Rmd`, `solver-modes.Rmd`) were
drafted by Claude in a single session under author direction. Code
chunks are run on every package build via `knitr`.

**Verification:**
* `devtools::build_vignettes()` builds both in 7.5 seconds with no
  errors.
* The Rosenbrock and quartic-tail examples used in the vignettes
  match unit tests in `tests/testthat/test-arcopt.R`,
  `tests/testthat/test-tr-fallback-integration.R`, and
  `tests/testthat/test-qn-polish-integration.R`, so any divergence
  between vignette behavior and tested behavior would surface as a
  vignette-knit error.
* The vignettes consume only the public API and exhibit no
  dependencies on `benchmarks/` or external packages beyond
  `Suggests` (no Stan, no BridgeStan, no `numDeriv` at runtime).

### 2.5 JSS manuscript (`manuscript/arcopt-jss.qmd`)

The manuscript is the result of an iterative restructure (v2,
2026-04-26 to 2026-04-27) that switched from a benchmark-suite
structure to a mode-organized structure with one example per
subproblem mode. Substantial Claude assistance went into:

* Restructuring §4 from 6 generic test problems to 4 mode-organized
  examples (§4.1 mixture saddle, §4.2 GMM, §4.3 DPM Heckman, §4.4
  sparse Rasch JML).
* Writing §4.4 from scratch based on validated benchmarks.
* Building the §4.5 master synthesis table.
* Restructuring §5 Discussion into §5.1 (when each mode wins),
  §5.2 (limitations), §5.3 (when arcopt is the wrong choice — with
  an honest NIST head-to-head table showing where LM solvers
  outperform arcopt), and §5.4 (future directions).
* Updating the abstract to reflect the mode-organized examples.

The §4.1 (mixture saddle), §4.2 GMM main table and discussion, and
§4.3 DPM Heckman example were primarily authored prior to the v2
restructure with light editing during it.

**Verification:**
* The manuscript renders cleanly via Quarto (`quarto render
  arcopt-jss.qmd --to pdf`) with no warnings or unresolved
  cross-references; the rendered PDF is currently 290 KB / 40 pages.
* All numerical claims in §4.4 derive from
  `manuscript/_cache/rasch_jml_results.rds`, computed by a runnable
  benchmark script (see §2.3 above).
* The §5.3 NIST head-to-head data derive from
  `benchmarks/nist_strd/results/headtohead.csv`. arcopt's
  underperformance on the NIST higher-difficulty far-start regime is
  documented honestly rather than rationalized — Claude's role
  included surfacing this finding to the author rather than burying
  it.
* Bibliography entries added during the restructure (notably
  Andersen 1973 for the JML inconsistency caveat) were
  hand-verified against the cited journal volumes; speculative
  citations the author was unsure of were removed rather than
  retained (e.g., a candidate Computer-Adaptive Testing reference
  for §4.4 was removed because the author could not confirm the
  citekey).
* The manuscript prose is the author's responsibility to review
  before JSS submission. Acknowledged limitation: peer-review
  confirmation is *not yet available* at v0.3.0; this report will be
  revised after JSS reviewer feedback.

### 2.6 Documentation (`man/`, `README.md`, `NEWS.md`,
`?arcopt_advanced_controls`)

Roxygen comments in `R/` are mostly co-authored across sessions.
`README.md` is primarily author-written. `NEWS.md`'s v0.2.0 entry
predates this session; the v0.3.0 entry was drafted in this session.
`arcopt_advanced_controls` topic page is structured as a
documentation-only Rd page with the author defining the tier
hierarchy (Tier 1 in `?arcopt`, Tier 2 here).

**Verification:**
* `devtools::check()`'s `checking Rd cross-references` check passes
  (all `\link{...}` references resolve).
* `devtools::document()` regenerates `man/*.Rd` from `R/*.R` roxygen
  blocks; any inconsistency between source and Rd is caught at
  build.
* All 9 cross-references checked manually in the documentation audit
  resolve to documented functions.

### 2.7 CRAN-submission infrastructure

The `cran-comments.md`, `vignettes/.gitignore`, the `.Rbuildignore`
additions for `benchmarks/`, `hexsticker/`, `manuscript/`, `lib/`, and
the `DESCRIPTION` email and version-bump edits were drafted in this
session and reviewed by the author.

**Verification:**
* All checks rerun after each edit: `devtools::check()`,
  `devtools::test()`, `devtools::build_vignettes()`, and finally
  `devtools::check(args = "--as-cran")`. All produced 0/0/0 modulo
  the known offline-environment NOTE.

---

## 3. Verification mechanisms summary

The verification stack, ordered from most-automated to most-manual:

| Layer | Catches | Run when |
|---|---|---|
| `tests/testthat/` (594 tests, 0 skips on CRAN) | Algorithm regressions | Every `devtools::test()` and CI run |
| `R CMD check --as-cran` | Doc drift, NAMESPACE inconsistency, example failures, Rd cross-ref breakage, package-build sanity | Before each release; on every CI push |
| Multi-platform CI (`.github/workflows/`) | Platform-specific issues (R-devel, R-oldrel, OS differences) | Every push to main |
| FD derivative checks (`benchmarks/utils.R`) | Gradient / Hessian errors in new problem setups | When introducing a new benchmark or example |
| Vignette knit | Drift between documented behavior and actual API | Every package build |
| Cache-driven manuscript inline numbers | Drift between manuscript claims and benchmark results | Every manuscript render |
| Author manual review | Framing, scope, judgment calls, mathematical correctness Claude can't verify by itself | Before each commit, release, or manuscript revision |

---

## 4. Known limitations of this report

* **No session-level attribution.** This report describes the state
  of v0.3.0 and not which Claude session contributed which line. The
  author's coding history is in `git log`; AI-conversation transcripts
  are not currently archived in the repository.
* **JSS peer review pending.** The manuscript described in §2.5 has
  not yet been through JSS review at the time of this report; reviewer
  feedback may surface issues that require revision and a follow-up
  acknowledgment.
* **Verification depth varies.** §2.1 algorithm code is verified by
  594 tests; §2.5 manuscript prose is verified by author review only.
  Readers should weight claims accordingly.
* **No formal review of Claude prompts.** The prompt history is not
  audited here. The author affirms that no prompts attempted to
  bypass standard verification mechanisms or to misrepresent results.

---

## 5. Contact

For questions about specific components or verification details:
**Marcus Waldman** — marcus.waldman@cuanschutz.edu —
University of Colorado Anschutz Medical Campus, Department of
Biostatistics and Informatics.

This report is versioned in the package source repository at
[https://github.com/marcus-waldman/arcopt](https://github.com/marcus-waldman/arcopt)
and will be revised on subsequent releases.

*Last updated: 2026-04-27 (v0.3.0 first CRAN submission).*
