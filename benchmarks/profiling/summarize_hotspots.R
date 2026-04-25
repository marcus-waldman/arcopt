# benchmarks/profiling/summarize_hotspots.R
# ==========================================
# Aggregates the 9 Rprof samples produced by run_profiling.R into a
# single markdown report identifying C++ port candidates inside arcopt.
#
# The key analytical step is splitting samples into three buckets:
#   - arcopt-internal: time inside arcopt's own R code (genuine ports
#     candidates)
#   - user-callback : time inside the user's fn / gr / hess (the user's
#     job to optimize, NOT a package port candidate)
#   - other         : base R, stats, BridgeStan AD calls, etc.
#
# Output: benchmarks/profiling/results/summary.md
#
# Usage (from repo root):
#   Rscript benchmarks/profiling/summarize_hotspots.R

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(arcopt)
  }
})

results_dir <- file.path("benchmarks", "profiling", "results")
out_file <- file.path(results_dir, "summary.md")

# --------------------------------------------------------------------
# Discover arcopt's function names (exports + internals) and the set
# of source-file basenames so we can classify samples without relying
# on namespace prefixes appearing in the trace.
# --------------------------------------------------------------------
arcopt_fns <- ls(asNamespace("arcopt"), all.names = TRUE)
arcopt_source_files <- basename(list.files("R", pattern = "\\.R$",
                                           full.names = TRUE))

# Names of user-supplied callbacks. These come from setup_examples.R
# (closure variable names) plus the formal parameter names arcopt uses
# when invoking them ("fn", "gr", "hess").
user_callback_fns <- c(
  "fn", "gr", "hess", "hess_ad", "hess_fd",
  "mixture_nll", "mixture_grad", "mixture_hess",
  "gmm_nll", "gmm_grad", "gmm_hess"
)

classify_fn <- function(name) {
  unquoted <- gsub('^"|"$', "", name)
  bare <- sub("^.*::", "", unquoted)
  if (bare %in% arcopt_fns) return("arcopt")
  if (bare %in% user_callback_fns) return("user")
  if (grepl("^arcopt::", unquoted)) return("arcopt")
  "other"
}

# --------------------------------------------------------------------
# Per-Rprof-file summary. Returns a list with:
#   by_self    : data.frame(fn, self_pct, total_pct, bucket)
#   by_line    : data.frame(line_id, self_pct, total_pct) (or NULL)
#   sampling_s : total elapsed sampling time (seconds)
# --------------------------------------------------------------------
summarize_rprof <- function(path) {
  s_fn <- tryCatch(
    summaryRprof(path, lines = "hide"),
    error = function(e) {
      warning("summaryRprof (function-level) failed for ", path, ": ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(s_fn)) return(NULL)

  bs <- s_fn$by.self
  bs$fn <- rownames(bs)
  bs$bucket <- vapply(bs$fn, classify_fn, character(1))
  rownames(bs) <- NULL

  s_line <- tryCatch(
    summaryRprof(path, lines = "show"),
    error = function(e) NULL
  )
  by_line <- NULL
  if (!is.null(s_line) && !is.null(s_line$by.line) &&
      nrow(s_line$by.line) > 0) {
    bl <- s_line$by.line
    bl$line_id <- rownames(bl)
    rownames(bl) <- NULL
    by_line <- bl
  }

  list(
    by_self    = bs,
    by_line    = by_line,
    sampling_s = if (is.null(s_fn$sampling.time)) NA_real_ else
      s_fn$sampling.time
  )
}

# --------------------------------------------------------------------
# Discover all Rprof files in results/
# --------------------------------------------------------------------
prof_files <- list.files(results_dir, pattern = "\\.Rprof$",
                         full.names = TRUE)
if (length(prof_files) == 0L) {
  stop("no .Rprof files found in ", results_dir,
       " -- run run_profiling.R first")
}

# Parse <problem>_<mode>.Rprof filenames (mode may contain a hyphen,
# e.g. cubic-default).
parse_label <- function(path) {
  stem <- sub("\\.Rprof$", "", basename(path))
  problem <- sub("_.*$", "", stem)
  mode <- sub("^[^_]+_", "", stem)
  list(problem = problem, mode = mode, label = stem)
}

per_run <- lapply(prof_files, function(p) {
  meta <- parse_label(p)
  s <- summarize_rprof(p)
  if (is.null(s)) return(NULL)
  c(meta, s)
})
per_run <- per_run[!vapply(per_run, is.null, logical(1))]

# --------------------------------------------------------------------
# Cross-run aggregation: total self-time per arcopt-internal function,
# weighted by sample count from each run. We use self.pct * sampling_s
# / 100 as a seconds-equivalent contribution, which works as long as
# sampling_s is recorded.
# --------------------------------------------------------------------
agg <- list()
for (r in per_run) {
  bs <- r$by_self
  bs <- bs[bs$bucket == "arcopt", , drop = FALSE]
  if (nrow(bs) == 0) next
  weight <- if (is.finite(r$sampling_s) && r$sampling_s > 0) {
    r$sampling_s
  } else 1
  for (i in seq_len(nrow(bs))) {
    fn <- bs$fn[i]
    contribution <- bs$self.pct[i] / 100 * weight
    prior <- if (is.null(agg[[fn]])) 0 else agg[[fn]]
    agg[[fn]] <- prior + contribution
  }
}
if (length(agg) == 0L) {
  arcopt_global_top <- character(0)
  agg_df <- data.frame(fn = character(0), seconds = numeric(0))
} else {
  agg_df <- data.frame(
    fn = names(agg),
    seconds = unlist(agg, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  agg_df <- agg_df[order(-agg_df$seconds), ]
  arcopt_global_top <- head(agg_df$fn, 3)
}

# --------------------------------------------------------------------
# Markdown emit helpers
# --------------------------------------------------------------------
fmt_pct <- function(x) sprintf("%5.1f%%", x)
fmt_sec <- function(x) sprintf("%.2fs", x)

emit_table <- function(df, cols, headers, top = 10) {
  if (is.null(df) || nrow(df) == 0) {
    return("_(no samples)_\n")
  }
  df <- head(df, top)
  hdr <- paste("|", paste(headers, collapse = " | "), "|")
  sep <- paste("|", paste(rep("---", length(headers)), collapse = " | "),
               "|")
  rows <- vapply(seq_len(nrow(df)), function(i) {
    vals <- vapply(cols, function(col) {
      v <- df[[col]][i]
      if (is.numeric(v)) sprintf("%.2f", v) else as.character(v)
    }, character(1))
    paste("|", paste(vals, collapse = " | "), "|")
  }, character(1))
  paste(c(hdr, sep, rows), collapse = "\n")
}

# --------------------------------------------------------------------
# Build the markdown report
# --------------------------------------------------------------------
md <- c(
  "# arcopt profiling: C++ port candidates",
  "",
  sprintf("_Generated %s_", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "",
  "Each profile is split into three buckets:",
  "",
  "- **arcopt** — time inside arcopt's R code (port candidates).",
  "- **user** — time inside the user-supplied `fn`/`gr`/`hess` (the",
  "  user's responsibility, not the package's).",
  "- **other** — base R, stats, BridgeStan, etc.",
  "",
  "## Cross-run ranking of arcopt-internal functions",
  "",
  "Total seconds-equivalent across all (problem x mode) runs.",
  ""
)

if (nrow(agg_df) == 0) {
  md <- c(md, "_(no arcopt-internal samples in any run)_", "")
} else {
  agg_top <- head(agg_df, 15)
  agg_top$rank <- seq_len(nrow(agg_top))
  md <- c(md, emit_table(
    agg_top[, c("rank", "fn", "seconds")],
    cols = c("rank", "fn", "seconds"),
    headers = c("rank", "function", "seconds"),
    top = 15
  ), "")
  md <- c(md,
    sprintf("**Global top 3 (line-level drill-down below):** %s",
            paste(arcopt_global_top, collapse = ", ")),
    "")
}

md <- c(md, "## Per-run function-level top 10", "")

for (r in per_run) {
  md <- c(md,
    sprintf("### %s / %s  (sampling = %s)", r$problem, r$mode,
            fmt_sec(r$sampling_s)),
    "")
  bs <- r$by_self
  bs <- bs[order(-bs$self.pct), ]
  cols <- c("fn", "self.pct", "total.pct", "bucket")
  md <- c(md, emit_table(bs[, cols],
                         cols = cols,
                         headers = c("function", "self %", "total %",
                                     "bucket"),
                         top = 10), "")

  bucket_totals <- tapply(bs$self.pct, bs$bucket, sum)
  if (!is.null(bucket_totals)) {
    bt_lines <- sprintf("- **%s**: %s",
                        names(bucket_totals),
                        fmt_pct(as.numeric(bucket_totals)))
    md <- c(md, "Bucket totals (self %):", "", bt_lines, "")
  }
}

# --------------------------------------------------------------------
# Line-level drill-down: top hot lines that fall inside any arcopt R
# source file. (Mapping line_id back to enclosing function would need
# source-ref parsing; the file+line is enough to navigate to the code.)
# --------------------------------------------------------------------
md <- c(md, "## Line-level drill-down (hot lines in arcopt R source)",
        "")

is_arcopt_line <- function(line_id) {
  # line_id looks like "arcopt.R#649" or "<no location>"
  file_part <- sub("#.*$", "", line_id)
  file_part %in% arcopt_source_files
}

any_lines_emitted <- FALSE
for (r in per_run) {
  bl <- r$by_line
  if (is.null(bl)) next
  bl <- bl[is_arcopt_line(bl$line_id), , drop = FALSE]
  if (nrow(bl) == 0) next
  bl <- bl[order(-bl$self.pct), ]
  md <- c(md, sprintf("### %s / %s", r$problem, r$mode), "")
  md <- c(md, emit_table(bl,
                         cols = c("line_id", "self.pct", "total.pct"),
                         headers = c("file#line", "self %", "total %"),
                         top = 15), "")
  any_lines_emitted <- TRUE
}
if (!any_lines_emitted) {
  md <- c(md, "_(no samples landed inside arcopt R sources)_", "")
}

writeLines(md, out_file)
cat(sprintf("[summarize] wrote %s\n", out_file))
