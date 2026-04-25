# benchmarks/profiling/compare_baseline.R
# =======================================
# Quick A/B comparison of two Rprof files. Saves a snapshot of an
# Rprof file as `<stem>.baseline.Rprof` BEFORE making a change, then
# compares it with the post-change version of the same name.
#
# Usage (from repo root):
#   # 1. Save baseline before the change:
#   #    cp results/<stem>.Rprof results/<stem>.baseline.Rprof
#   # 2. Re-run profiling.
#   # 3. Compare:
#   #    Rscript benchmarks/profiling/compare_baseline.R <stem>
#
# Default stem: dpm_heckman_cubic-default

args <- commandArgs(trailingOnly = TRUE)
stem <- if (length(args) >= 1) args[1] else "dpm_heckman_cubic-default"

baseline_path <- file.path("benchmarks", "profiling", "results",
                           paste0(stem, ".baseline.Rprof"))
current_path <- file.path("benchmarks", "profiling", "results",
                          paste0(stem, ".Rprof"))

for (p in c(baseline_path, current_path)) {
  if (!file.exists(p)) stop("missing Rprof file: ", p)
}

fmt_self <- function(p, label) {
  s <- summaryRprof(p)
  bs <- s$by.self
  cat(sprintf("\n=== %s (%s, sampling=%.2fs) ===\n",
              label, basename(p), s$sampling.time))
  targets <- c('"eigen"', '".C"', '"solve_cubic_eigen"', '"hess"',
               '"chol.default"', '"solve_cubic_subproblem_dispatch"',
               '"solve_tr_eigen"')
  for (fn in targets) {
    if (fn %in% rownames(bs)) {
      cat(sprintf("  %-35s self=%5.2f%% total=%5.2f%%\n", fn,
                  bs[fn, "self.pct"], bs[fn, "total.pct"]))
    }
  }
}

fmt_lines <- function(p, label, pattern) {
  s <- summaryRprof(p, lines = "show")
  bl <- s$by.line
  bl <- bl[order(-bl$self.pct), ]
  hits <- grep(pattern, rownames(bl))
  cat(sprintf("\n--- %s line-level (matching '%s') ---\n",
              label, pattern))
  print(head(bl[hits, c("self.pct", "total.pct")], 8))
}

fmt_self(baseline_path, "BASELINE")
fmt_self(current_path, "CURRENT")
fmt_lines(baseline_path, "BASELINE",
          "(cubic_eigen|tr_eigen|arcopt\\.R#7[01][0-9])")
fmt_lines(current_path, "CURRENT",
          "(cubic_eigen|tr_eigen|arcopt\\.R#7[01][0-9])")
