# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

**arcopt** (Adaptive Regularization using Cubics for R) is an
optimization package implementing cubic regularization methods (ARC) for
local optimization problems common in statistics and applied research.

- **Primary Users**: R users (statisticians, applied researchers)
- **Problem Scale**: Small to medium (typically 2-500 parameters)
- **GitHub**: <https://github.com/marcus-waldman/arcopt>

## Core Architecture

The optimizer implements a hierarchical 4-layer system:

1.  **Orchestration Layer**: Main loop (Algorithm 1) that controls
    iteration flow and dispatches to sub-algorithms
2.  **Core Solvers Layer**: Newton step computation (Algorithm 1a) and
    eigendecomposition cubic solver (Algorithm 5a)
3.  **Adaptation Layer**: Convergence checks (Algorithm 0) and σ (sigma)
    regularization updates (Algorithms 2a/2b)
4.  **Safeguards Layer**: Indefiniteness handling (Algorithms 6a/6b),
    box constraint truncation (Algorithm 7), and linear equality
    constraint handling (Algorithm 8)

### Key Design Philosophy

- **Hessian-Centric**: Accurate second derivative information is
  prioritized. Hessian sources in order of preference: analytic \>
  automatic differentiation \> finite differences \> quasi-Newton
  approximation. QN-ARC mode available via `use_qn = TRUE`.
- **Robust by Default**: Single entry point `arcopt(x0, fn, gr, hess)`
  with sensible defaults; users should not need to manually tune
  regularization or select solvers
- **Saddle Point Escape**: Cubic regularization naturally handles
  negative curvature without explicit eigenvector computation
- **Constraint Handling**: Box constraints via step truncation; linear
  equality constraints via null-space transformation

## Documentation

- **NEWS.md**: Source-of-truth changelog. Read first when picking up a
  new session — it lists what changed in each release.
- **design/design-principles.qmd**: Design philosophy, target user
  profile, problem characteristics, and core principles
- **design/pseudocode.qmd**: Complete algorithmic specifications (8
  algorithm groups) with system flowchart and detailed pseudocode. Note:
  Algorithm 3 (momentum) and Algorithm 4a-lhybrid (limited-memory)
  sections are flagged as removed; preserved for historical reference.
- **design/literature-review.qmd**: Academic references and theoretical
  foundations
- **literature/consensus_reviews/**: Topic-specific reviews on
  regularization parameters, failure modes, and comparisons with other
  ARC implementations

## v0.2.0 user surface (post-simplification)

The user-facing control surface was tiered in v0.2.0. When advising
users:

- **[`?arcopt`](https://marcus-waldman.github.io/arcopt/reference/arcopt.md)**
  documents only seven Tier 1 controls: `maxit`, `gtol_abs`, `ftol_abs`,
  `xtol_abs`, `trace`, `verbose`, `use_qn`. The first four are
  tolerances; `trace` (0/1/2/3) controls saved per-iteration data depth;
  `verbose` (TRUE/FALSE) prints one line per iteration to the console
  (orthogonal to `trace`); `use_qn = TRUE` routes to the (non-exported)
  quasi-Newton variant.
- **[`?arcopt_advanced_controls`](https://marcus-waldman.github.io/arcopt/reference/arcopt_advanced_controls.md)**
  documents everything else — cubic regularization tuning (`sigma0`,
  `eta1/2`, `gamma1/2`, …), the trust-region fallback (`tr_fallback_*`,
  `tr_*`), the optional polish mode (`qn_polish_*`, default off), and
  the QN routing parameters (`qn_method`, `qn_route_*`).
- **`result$diagnostics`** is a sublist holding mode-dispatch
  diagnostics: `solver_mode_final`, `ridge_switches`, `radius_final`,
  `qn_polish_switches`, `qn_polish_reverts`,
  `hess_evals_at_polish_switch`. QN runs additionally include
  `qn_updates`, `qn_skips`, `qn_restarts`, `qn_fd_refreshes`. Most users
  do not need to inspect these.
- **[`arcopt_qn()`](https://marcus-waldman.github.io/arcopt/reference/arcopt_qn.md)
  is internal** (not in NAMESPACE). Users go through
  `arcopt(use_qn = TRUE)`. The internal function still has its own
  [`?arcopt_qn`](https://marcus-waldman.github.io/arcopt/reference/arcopt_qn.md)
  help page for power users.

## Development Guidelines

### Project Stage

Current state (v0.2.0.9000): tri-modal solver shipped (cubic /
TR-fallback / qn_polish-opt-in), 594 unit tests passing,
`devtools::check()` clean modulo four pre-existing notes about
non-package directories (`benchmarks/`, `manuscript/`, `hexsticker/`).
JSS manuscript at `manuscript/arcopt-jss.qmd` (gitignored, local-only);
rendered PDF ~282 KB / 57 pages. CRAN submission not yet attempted.

### When Implementing Features

1.  **Follow the pseudocode specifications** in design/pseudocode.qmd
    exactly—these are implementation-ready
2.  **Maintain Hessian-centric philosophy**: Default interface requires
    Hessian function; provide finite-difference wrapper as convenience
    option
3.  **Prioritize robustness**: Include automatic indefiniteness
    detection, adaptive σ adjustment, and numerical safeguards
4.  **Test on pathological problems**: Ill-conditioned, nonconvex,
    indefinite Hessian, saddle point, and ridge maximum cases
5.  **Reference design principles** when making architectural
    decisions—existing decisions are well-justified

### Implemented

- Core cubic regularization solver (Algorithm 5a) with
  eigendecomposition
- Adaptive regularization updates (Algorithms 2a/2b) for σ adjustment
- Indefiniteness handling for negative curvature
- Box-constraint handling (Algorithm 7) via step truncation
- Quasi-Newton variants (Algorithm 4): SR1, BFGS, hybrid (BFGS → SR1 →
  Powell-damped) with FD-Hessian seeding of B₀
- Trust-region fallback (cubic → TR on flat-ridge stagnation, default
  on)
- QN-polish mode (cubic ↔︎ Wolfe line-search BFGS in healthy basin,
  opt-in)

### Deferred Features

To maintain focus on the core use case (statisticians, 2–500 parameters,
analytic Hessians):

**Removed in v0.2.0 (2026-04-24):** - Momentum acceleration
(`use_momentum`, `momentum_tau`, `momentum_alpha1`, `momentum_alpha2`
controls + the Gao et al. 2022 ARCm bisection block); see Algorithm 3 in
`design/pseudocode.qmd` (banner-marked). - `cubic_solver` user-facing
knob (the internal dispatcher always selects the eigendecomposition
solver; preserved as an extensibility hook for the deferred Algorithm
5b).

**Removed before v0.2.0:** - L-BFGS / L-SR1 / L-Hybrid quasi-Newton
variants (commit `d32abde`, 2026-04-20). Preserved on the `scalable-arc`
branch pending Algorithm 5b. - LDL-based cubic solver (see
`design/historical/ldl-solver.qmd`); replaced by eigendecomposition. -
Linear equality constraint handling (Algorithm 8) — design retained,
implementation deferred.

**Deferred to future releases:** - ARCqK multi-shift CG-Lanczos solver
(Algorithm 5b) for n \> 500 - Matrix-free optimization via `hess_vec`
interface - `tr → qn_polish` transition (v3 feature; today TR is a
terminal mode) - `qn_polish_enabled = TRUE` as default (currently opt-in
pending broader benchmark evidence)

See `design/scalable-arcs.qmd` for preserved L-\* implementations and
design rationale.

## Development Workflow & CRAN Compliance

### Required Tools

``` r
install.packages(c("devtools", "usethis", "roxygen2", "testthat",
                   "pkgdown", "lintr", "styler", "covr", "Rcpp"))
```

### Daily Development Workflow

**Before committing any code**, always run this sequence:

1.  `devtools::load_all()` - Load package for testing
2.  `devtools::document()` - Update documentation from roxygen2 comments
3.  `lintr::lint_package()` - **CRITICAL**: Check Tidyverse style
    compliance
4.  `styler::style_pkg()` - Auto-fix formatting issues
5.  `devtools::test()` - Run all tests
6.  `devtools::check()` - Full R CMD check (must pass with 0 errors, 0
    warnings, 0 notes)

**After modifying C++ code in src/:** - Run
[`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
to regenerate RcppExports.R/cpp - Then run `devtools::document()` to
update NAMESPACE

### Style Requirements

**This package strictly follows the [Tidyverse Style
Guide](https://style.tidyverse.org/)**

Enforced by `.lintr` configuration: - **snake_case** naming for
functions and variables (e.g., `solve_cubic`, not `solveCubic`) - **120
character** line length maximum - **2-space indentation** (never tabs) -
**roxygen2 documentation** required for all exported functions - **No
browser(), print() debugging** statements in production code

**Common violations to avoid:** - CamelCase or PascalCase names (use
snake_case) - Lines over 120 characters (break long lines) - Missing or
incomplete roxygen2 documentation - Trailing whitespace

### Pre-Commit Checklist

Before every commit: - \[ \] `lintr::lint_package()` shows no errors -
\[ \] `devtools::check()` passes (0 errors, 0 warnings, 0 notes) - \[ \]
All tests pass with `devtools::test()` - \[ \] roxygen2 documentation is
current - \[ \] No debugging code (browser(), print()) in commits

### CRAN Compliance

The package must always be ready for CRAN submission:

**DESCRIPTION file:** - Title in title case - Description with 2+
sentences, doesn’t start with package name - All dependencies declared
in Imports/Suggests - Rcpp in both Imports and LinkingTo

**Code requirements:** - All exported functions have complete roxygen2
documentation - Examples in documentation actually run (or use
`\donttest{}`) - No absolute file paths (use
[`system.file()`](https://rdrr.io/r/base/system.file.html) for package
files) - Proper `importFrom()` or `::` for external functions

**Rcpp-specific:** - After editing C++ files, run
[`Rcpp::compileAttributes()`](https://rdrr.io/pkg/Rcpp/man/compileAttributes.html)
before committing - RcppExports.R and RcppExports.cpp should be
committed (tracked by Git) - Makevars and Makevars.win configure
compilation properly

### CI/CD Integration

GitHub Actions automatically enforce standards on every push: -
**R-CMD-check.yaml**: Tests on 5 platforms (macOS, Windows, Linux ×
multiple R versions) - **lint.yaml**: Fails if lintr detects style
violations - **test-coverage.yaml**: Reports code coverage to
codecov.io - **pkgdown.yaml**: Builds and deploys package website

**All workflows must pass before merging PRs** (configure GitHub branch
protection).

### Implementation Guidelines

When implementing algorithms from design/pseudocode.qmd: 1. Start with
pure R implementation (easier to debug) 2. Write comprehensive tests
first (test-driven development) 3. Profile code with
`profvis::profvis()` to identify bottlenecks 4. Move hot paths to C++
only after profiling confirms need 5. Maintain equivalence tests between
R and C++ implementations

### Reference

- [R Packages (2e)](https://r-pkgs.org/) - Comprehensive R package
  development guide
- [Tidyverse Style Guide](https://style.tidyverse.org/) - Required
  coding style
- [lintr documentation](https://lintr.r-lib.org/) - Style checker
  reference
- [Rcpp for R Package
  Development](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-package.pdf)
