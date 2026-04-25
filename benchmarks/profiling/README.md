# Profiling: identifying C++ port candidates

This directory profiles `arcopt()` on the three manuscript examples
(mixture saddle, GMM, DPM Heckman MAP) across all three solver modes
(`cubic-default`, `qn`, `qn_polish`) to identify where wall time is
spent and which functions are worth porting to C++.

## Files

- `setup_examples.R` — problem builders. Each returns
  `list(label, n, x0, fn, gr, hess)`. Mixture and GMM are lifted from
  the JSS manuscript chunks; DPM Heckman wraps the existing
  `benchmarks/dpm_heckman/setup.R` (requires `bridgestan`).
- `run_profiling.R` — runs `profvis()` on each (problem, mode) pair,
  writes `<problem>_<mode>.html` (interactive viewer) and
  `<problem>_<mode>.Rprof` (raw samples) to `results/`.
- `summarize_hotspots.R` — reads all 9 Rprof files, classifies each
  function as **arcopt-internal**, **user-callback**, or **other**,
  and emits `results/summary.md` with cross-run rankings and
  line-level drill-down for the top-3 arcopt functions.

## Running the sweep

From the repo root:

```bash
Rscript benchmarks/profiling/run_profiling.R
Rscript benchmarks/profiling/summarize_hotspots.R
```

Selective runs (during iteration) via `--problem` and `--mode` flags
(comma-separated):

```bash
Rscript benchmarks/profiling/run_profiling.R --problem mixture,gmm
Rscript benchmarks/profiling/run_profiling.R --mode cubic-default
```

Expected total wall time: roughly 30 minutes, dominated by the three
DPM Heckman runs (each ~5–10 min). Mixture and GMM together finish in
under a minute even with their repetition counts.

## Interpreting the buckets

The three buckets exist because **only arcopt-internal time is a port
candidate**:

| Bucket | What it is | C++ candidate? |
|---|---|---|
| `arcopt`  | arcopt's own R code (e.g., `solve_cubic`, `update_sigma`) | **Yes** |
| `user`    | the user's `fn`/`gr`/`hess` callbacks                     | No — user's job |
| `other`   | base R, stats, BridgeStan AD, etc.                        | No |

Note that the mixture and GMM examples use **finite-difference
Hessians written in pure R**, which will dominate their profiles.
Those samples land in the `user` bucket and should be ignored when
choosing what to port. The DPM Heckman profile uses Stan AD (already
C++), so its `arcopt` bucket is the most diagnostically meaningful for
package-internal hotspots.

## Reproducibility

Each problem builder uses a fixed seed (`42` for mixture/GMM, `101`
for DPM Heckman) so reruns produce the same trajectory. `results/` is
gitignored — re-run the scripts to regenerate.

## Output: results/

After running both scripts:

```
results/
  mixture_cubic-default.{html,Rprof}
  mixture_qn.{html,Rprof}
  mixture_qn_polish.{html,Rprof}
  gmm_cubic-default.{html,Rprof}
  gmm_qn.{html,Rprof}
  gmm_qn_polish.{html,Rprof}
  dpm_heckman_cubic-default.{html,Rprof}
  dpm_heckman_qn.{html,Rprof}
  dpm_heckman_qn_polish.{html,Rprof}
  summary.md            <- start here
```

Open the HTML files in a browser for the interactive flame graph view;
read `summary.md` for the textual ranking and line-level drill-down.
