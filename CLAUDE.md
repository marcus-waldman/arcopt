# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**arcopt** (Adaptive Regularization using Cubics for R) is an optimization package implementing cubic regularization methods (ARC) for local optimization problems common in statistics and applied research.

- **Primary Users**: R users (statisticians, applied researchers)
- **Problem Scale**: Small to medium (typically 2-500 parameters)
- **GitHub**: https://github.com/marcus-waldman/arcopt

## Core Architecture

The optimizer implements a hierarchical 4-layer system:

1. **Orchestration Layer**: Main loop (Algorithm 1) that controls iteration flow and dispatches to sub-algorithms
2. **Core Solvers Layer**: Newton step computation (Algorithm 1a) and cubic subproblem solver (Algorithm 5)
3. **Adaptation Layer**: Convergence checks (Algorithm 0), σ (sigma) regularization updates (Algorithms 2a/2b), momentum (Algorithm 3), and SR1 quasi-Newton updates (Algorithm 4)
4. **Safeguards Layer**: Indefiniteness handling (Algorithms 6a/6b), box constraint truncation (Algorithm 7), and linear equality constraint handling (Algorithm 8)

### Key Design Philosophy

- **Hessian-Centric**: Accurate second derivative information is prioritized. Hessian sources in order of preference: analytic > automatic differentiation > finite differences > SR1 quasi-Newton
- **Robust by Default**: Single entry point `arcopt(x0, fn, gr, hess)` with sensible defaults; users should not need to manually tune regularization or select solvers
- **Saddle Point Escape**: Cubic regularization naturally handles negative curvature without explicit eigenvector computation
- **Constraint Handling**: Box constraints via step truncation; linear equality constraints via null-space transformation

## Documentation

- **design/design-principles.qmd**: Design philosophy, target user profile, problem characteristics, and core principles
- **design/pseudocode.qmd**: Complete algorithmic specifications (9 algorithm groups) with system flowchart and detailed pseudocode
- **design/literature-review.qmd**: Academic references and theoretical foundations
- **literature/consensus_reviews/**: Topic-specific reviews on regularization parameters, quasi-Newton methods, failure modes, and comparisons with other ARC implementations

## Development Guidelines

### Project Stage

This is a research project in early development stages. The repository currently contains:
- Documentation and design specifications (no implementation yet)
- Git history not yet initialized on main branch
- Focus on algorithm specification before coding

### When Implementing Features

1. **Follow the pseudocode specifications** in design/pseudocode.qmd exactly—these are implementation-ready
2. **Maintain Hessian-centric philosophy**: Default interface requires Hessian function; provide convenience wrappers (FD, SR1) as secondary options
3. **Prioritize robustness**: Include automatic indefiniteness detection, adaptive σ adjustment, and numerical safeguards
4. **Test on pathological problems**: Ill-conditioned, nonconvex, indefinite Hessian, saddle point, and ridge maximum cases
5. **Reference design principles** when making architectural decisions—existing decisions are well-justified

### Future Development Priorities

Based on literature review and design documents:
- Core cubic regularization solver (Algorithm 5) with modified Cholesky factorization
- Adaptive regularization updates (Algorithms 2a/2b) for σ adjustment
- Indefiniteness handling (Algorithms 6a/6b) for negative curvature
- Constraint handling (Algorithms 7-8) for box and linear equality constraints
- SR1 quasi-Newton updates (Algorithm 4) as approximate Hessian option

## Development Workflow & CRAN Compliance

### Required Tools

```r
install.packages(c("devtools", "usethis", "roxygen2", "testthat",
                   "pkgdown", "lintr", "styler", "covr", "Rcpp"))
```

### Daily Development Workflow

**Before committing any code**, always run this sequence:

1. `devtools::load_all()` - Load package for testing
2. `devtools::document()` - Update documentation from roxygen2 comments
3. `lintr::lint_package()` - **CRITICAL**: Check Tidyverse style compliance
4. `styler::style_pkg()` - Auto-fix formatting issues
5. `devtools::test()` - Run all tests
6. `devtools::check()` - Full R CMD check (must pass with 0 errors, 0 warnings, 0 notes)

**After modifying C++ code in src/:**
- Run `Rcpp::compileAttributes()` to regenerate RcppExports.R/cpp
- Then run `devtools::document()` to update NAMESPACE

### Style Requirements

**This package strictly follows the [Tidyverse Style Guide](https://style.tidyverse.org/)**

Enforced by `.lintr` configuration:
- **snake_case** naming for functions and variables (e.g., `solve_cubic`, not `solveCubic`)
- **120 character** line length maximum
- **2-space indentation** (never tabs)
- **roxygen2 documentation** required for all exported functions
- **No browser(), print() debugging** statements in production code

**Common violations to avoid:**
- CamelCase or PascalCase names (use snake_case)
- Lines over 120 characters (break long lines)
- Missing or incomplete roxygen2 documentation
- Trailing whitespace

### Pre-Commit Checklist

Before every commit:
- [ ] `lintr::lint_package()` shows no errors
- [ ] `devtools::check()` passes (0 errors, 0 warnings, 0 notes)
- [ ] All tests pass with `devtools::test()`
- [ ] roxygen2 documentation is current
- [ ] No debugging code (browser(), print()) in commits

### CRAN Compliance

The package must always be ready for CRAN submission:

**DESCRIPTION file:**
- Title in title case
- Description with 2+ sentences, doesn't start with package name
- All dependencies declared in Imports/Suggests
- Rcpp in both Imports and LinkingTo

**Code requirements:**
- All exported functions have complete roxygen2 documentation
- Examples in documentation actually run (or use `\donttest{}`)
- No absolute file paths (use `system.file()` for package files)
- Proper `importFrom()` or `::` for external functions

**Rcpp-specific:**
- After editing C++ files, run `Rcpp::compileAttributes()` before committing
- RcppExports.R and RcppExports.cpp should be committed (tracked by Git)
- Makevars and Makevars.win configure compilation properly

### CI/CD Integration

GitHub Actions automatically enforce standards on every push:
- **R-CMD-check.yaml**: Tests on 5 platforms (macOS, Windows, Linux × multiple R versions)
- **lint.yaml**: Fails if lintr detects style violations
- **test-coverage.yaml**: Reports code coverage to codecov.io
- **pkgdown.yaml**: Builds and deploys package website

**All workflows must pass before merging PRs** (configure GitHub branch protection).

### Implementation Guidelines

When implementing algorithms from design/pseudocode.qmd:
1. Start with pure R implementation (easier to debug)
2. Write comprehensive tests first (test-driven development)
3. Profile code with `profvis::profvis()` to identify bottlenecks
4. Move hot paths to C++ only after profiling confirms need
5. Maintain equivalence tests between R and C++ implementations

### Reference

- [R Packages (2e)](https://r-pkgs.org/) - Comprehensive R package development guide
- [Tidyverse Style Guide](https://style.tidyverse.org/) - Required coding style
- [lintr documentation](https://lintr.r-lib.org/) - Style checker reference
- [Rcpp for R Package Development](https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-package.pdf)
