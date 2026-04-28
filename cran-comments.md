# cran-comments.md

## Submission

This is the first CRAN submission of **arcopt**.

## Test environments

Local:

* Windows 11, R 4.5.1 (release)

GitHub Actions (`.github/workflows/R-CMD-check.yaml`):

* macOS-latest, R-release
* windows-latest, R-release
* ubuntu-latest, R-devel
* ubuntu-latest, R-release
* ubuntu-latest, R-oldrel-1

Pre-submission checks (run 2026-04-27):

* `R CMD check --as-cran` (local, Windows R 4.5.1): 0 errors,
  0 warnings, 0 notes (a "checking for future file timestamps --
  unable to verify current time" NOTE appears in offline sandboxes
  due to lack of NTP access; this NOTE does not reproduce on a
  network-connected machine and does not appear on CRAN's check
  infrastructure).
* `devtools::check_win_devel()` (Windows R-devel via win-builder):
  0 errors, 0 warnings, 1 NOTE -- discussed below.
* `devtools::check_mac_release()` (macOS R-release via mac-builder):
  0 errors, 0 warnings, 0 notes.

## R CMD check results

0 errors | 0 warnings | 1 note (win-devel only).

The single NOTE on win-devel is the standard CRAN first-submission
NOTE plus a flag on two technical terms in the package DESCRIPTION:

```
Maintainer: 'Marcus Waldman <marcus.waldman@cuanschutz.edu>'

New submission

Possibly misspelled words in DESCRIPTION:
  Cubics (3:38)
  nonconvex (10:41)
```

* "Cubics" appears in the package title and refers to the cubic
  polynomial regularizer of Adaptive Regularization with Cubics
  (Cartis, Gould & Toint 2011); it is the established name of the
  algorithm class in the optimization literature.
* "nonconvex" is the standard adjective form used throughout the
  optimization and machine-learning literature (e.g., Nesterov &
  Polyak 2006; Cartis, Gould & Toint 2011).

Both are intentional. The "New submission" portion of the NOTE will
not appear on subsequent releases.

## Reverse dependencies

This is a new submission; there are no reverse dependencies.

## Notes for the reviewer

* **arcopt** implements Adaptive Regularization with Cubics (ARC) for
  nonlinear optimization, the first such implementation on CRAN. It
  is intended for statisticians and applied researchers solving
  maximum likelihood, MAP, and penalized regression problems with 2-500
  parameters where indefinite Hessians, saddle points, and
  ill-conditioning are common.

* The single user-facing entry point is `arcopt(x0, fn, gr, hess)`,
  documented under `?arcopt`. Advanced controls are documented
  separately under `?arcopt_advanced_controls`.

* Two short vignettes are included: `getting-started` walks through a
  complete optimization, and `solver-modes` describes the three
  adaptive subproblem modes (cubic regularization,
  trust-region fallback, optional quasi-Newton polish) with one
  small example each.

* All examples in help pages and vignettes run quickly and use no
  external resources or data.

* The package has 594 unit tests (in `tests/testthat/`) and runs all
  tests on every CI build across the platforms listed above.

## Disclosure: AI-assisted development

Portions of the package code, tests, vignettes, and accompanying
manuscript were drafted with the assistance of Anthropic's Claude
(Sonnet and Opus 4.x). All Claude-generated code was reviewed by the
author and verified by the package's 594-test suite, `R CMD check
--as-cran` on multiple platforms, and finite-difference derivative
checks where applicable. The author retains full responsibility for
the code and its correctness. A more detailed report on which
components were Claude-assisted and how each was verified is
available at `AI_USE_REPORT.md` in the source repository (excluded
from the package tarball via `.Rbuildignore`).
