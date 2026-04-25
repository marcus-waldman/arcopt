# benchmarks/profiling/run_profiling.R
# =====================================
# Runs profvis() on each (problem in {mixture, gmm, dpm_heckman}) x
# (mode in {cubic-default, qn, qn_polish}) and saves both an interactive
# HTML report and the raw Rprof samples for downstream aggregation.
#
# Total expected wall time: ~30 minutes, dominated by the three DPM
# Heckman runs.
#
# Usage (from repo root):
#   Rscript benchmarks/profiling/run_profiling.R
#
# Selective runs (by problem and/or mode):
#   Rscript benchmarks/profiling/run_profiling.R --problem mixture,gmm
#   Rscript benchmarks/profiling/run_profiling.R --mode cubic-default

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(arcopt)
  }
  if (!requireNamespace("profvis", quietly = TRUE)) {
    stop("package 'profvis' is required (install.packages('profvis'))")
  }
})

source(file.path("benchmarks", "profiling", "setup_examples.R"))

results_dir <- file.path("benchmarks", "profiling", "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# --------------------------------------------------------------------
# CLI arg parsing (very simple: --problem and --mode each take a comma
# separated list of names)
# --------------------------------------------------------------------
parse_csv_arg <- function(args, flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  strsplit(args[i + 1], ",", fixed = TRUE)[[1]]
}

cli <- commandArgs(trailingOnly = TRUE)
problems <- parse_csv_arg(cli, "--problem",
                          c("mixture", "gmm", "dpm_heckman"))
modes <- parse_csv_arg(cli, "--mode",
                       c("cubic-default", "qn", "qn_polish"))

# --------------------------------------------------------------------
# Mode -> control list (passed to arcopt()'s `control` argument).
# - cubic-default: package defaults (TR-fallback on, qn_polish off,
#   use_qn off)
# - qn: gradient-only quasi-Newton variant (does not call hess)
# - qn_polish: cubic + TR + opt-in BFGS polish in healthy basins
# --------------------------------------------------------------------
mode_controls <- list(
  "cubic-default" = list(maxit = 200L, gtol_abs = 1e-5, trace = 0L),
  "qn"            = list(maxit = 200L, gtol_abs = 1e-5, trace = 0L,
                         use_qn = TRUE),
  "qn_polish"     = list(maxit = 200L, gtol_abs = 1e-5, trace = 0L,
                         qn_polish_enabled = TRUE)
)

# --------------------------------------------------------------------
# Repetition counts per problem. Mixture and GMM finish in <1 s, so we
# loop them to accumulate enough samples (Rprof default interval is
# 20 ms). DPM Heckman is multi-minute on its own, so a single run is
# plenty.
# --------------------------------------------------------------------
problem_reps <- list(
  mixture     = 50L,
  gmm         = 20L,
  dpm_heckman = 1L
)

# --------------------------------------------------------------------
# Build problem objects up front so we don't re-simulate inside profvis
# (we only want optimizer time profiled, not data simulation).
# --------------------------------------------------------------------
build_problem <- function(name) {
  switch(name,
    mixture     = make_mixture_problem(),
    gmm         = make_gmm_problem(),
    dpm_heckman = make_dpm_heckman_problem(),
    stop("unknown problem: ", name)
  )
}

run_one <- function(problem, mode_name, reps) {
  ctrl <- mode_controls[[mode_name]]
  out_stem <- file.path(results_dir,
                        sprintf("%s_%s", problem$label, mode_name))
  prof_file <- paste0(out_stem, ".Rprof")
  html_file <- paste0(out_stem, ".html")

  cat(sprintf("[profile] %s / %s : %d rep(s)...\n",
              problem$label, mode_name, reps))

  # Warmup: a single un-profiled call so R's bytecode compiler finishes
  # JIT-compiling every closure before sampling begins. Otherwise
  # findCenvVar / cmp / cmpCall etc. pollute the early samples on
  # short-running problems (mixture, GMM).
  arcopt::arcopt(x0 = problem$x0, fn = problem$fn,
                 gr = problem$gr, hess = problem$hess,
                 control = ctrl)

  t0 <- proc.time()[["elapsed"]]

  # Collect samples with line.profiling so the summary can drill down.
  Rprof(prof_file, line.profiling = TRUE, interval = 0.01,
        memory.profiling = FALSE)
  on.exit(Rprof(NULL), add = TRUE)
  for (i in seq_len(reps)) {
    arcopt::arcopt(x0 = problem$x0, fn = problem$fn,
                   gr = problem$gr, hess = problem$hess,
                   control = ctrl)
  }
  Rprof(NULL)
  on.exit()

  elapsed <- proc.time()[["elapsed"]] - t0
  cat(sprintf("[profile]   wall = %.2fs\n", elapsed))

  # Render HTML viewer from the same Rprof data.
  pv <- profvis::profvis(prof_input = prof_file)
  htmlwidgets::saveWidget(pv, html_file, selfcontained = TRUE)
  cat(sprintf("[profile]   wrote %s\n", html_file))
  cat(sprintf("[profile]   wrote %s\n", prof_file))

  invisible(list(elapsed = elapsed, prof_file = prof_file,
                 html_file = html_file))
}

# --------------------------------------------------------------------
# Main sweep
# --------------------------------------------------------------------
total_t0 <- proc.time()[["elapsed"]]
for (pname in problems) {
  cat(sprintf("\n=== problem: %s ===\n", pname))
  problem <- build_problem(pname)
  reps <- problem_reps[[pname]]
  for (mname in modes) {
    run_one(problem, mname, reps)
  }
}
total_elapsed <- proc.time()[["elapsed"]] - total_t0
cat(sprintf("\n[profile] sweep complete: %.1fs total\n", total_elapsed))
