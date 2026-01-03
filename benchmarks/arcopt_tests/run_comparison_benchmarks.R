# Comparative Benchmark: arcopt vs trust
# ======================================
#
# Compares arcopt and trust::trust on 8 benchmark functions.
# Goal: Identify relative performance and correctness.

# Setup -------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
  script_dir <- dirname(normalizePath(script_path))
} else {
  script_dir <- getwd()
}
benchmarks_dir <- dirname(script_dir)
pkg_dir <- dirname(benchmarks_dir)

# Load packages
cat("Loading packages...\n")
devtools::load_all(pkg_dir)
library(trust)

# Source helper files
source(file.path(script_dir, "starting_points.R"))

# Configuration -----------------------------------------------------------

test_functions <- c(
  "sphere",
  "sum_of_squares",
  "bohachevsky",
  "rotated_hyper_ellipsoid",
  "trid",
  "powell_singular",
  "rosenbrock",
  "dixon_price"
)

# Control parameters
arcopt_control <- list(
  trace = FALSE,
  maxit = 1000,
  gtol_abs = 1e-6
)

trust_iterlim <- 1000

# Results storage
results <- data.frame(
 func_name = character(),
  start_id = integer(),
  difficulty = character(),
  description = character(),
  optimizer = character(),
  converged = logical(),
  iterations = integer(),
  final_grad_norm = numeric(),
  solution_error = numeric(),
  function_error = numeric(),
  time_sec = numeric(),
  error_msg = character(),
  stringsAsFactors = FALSE
)

# Load benchmark function -------------------------------------------------

load_benchmark_function <- function(func_name, benchmarks_dir) {
  func_dir <- file.path(benchmarks_dir, func_name)

  source(file.path(func_dir, "fn.R"))
  source(file.path(func_dir, "gr.R"))
  source(file.path(func_dir, "hess.R"))

  fn_name <- paste0(func_name, "_fn")
  gr_name <- paste0(func_name, "_gr")
  hess_name <- paste0(func_name, "_hess")

  list(
    fn = get(fn_name),
    gr = get(gr_name),
    hess = get(hess_name)
  )
}

# Create trust-compatible objfun wrapper
make_trust_objfun <- function(fn, gr, hess) {
  function(x) {
    list(
      value = fn(x),
      gradient = gr(x),
      hessian = hess(x)
    )
  }
}

# Run a single optimizer --------------------------------------------------

run_arcopt <- function(x0, funcs, control) {
  start_time <- Sys.time()
  result <- tryCatch(
    {
      res <- arcopt(
        x0 = x0,
        fn = funcs$fn,
        gr = funcs$gr,
        hess = funcs$hess,
        control = control
      )
      list(
        par = res$par,
        value = res$value,
        gradient = res$gradient,
        converged = res$converged,
        iterations = res$iterations,
        error = FALSE
      )
    },
    error = function(e) {
      list(error = TRUE, error_msg = conditionMessage(e))
    }
  )
  result$time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  result
}

run_trust <- function(x0, funcs, iterlim) {
  objfun <- make_trust_objfun(funcs$fn, funcs$gr, funcs$hess)

  start_time <- Sys.time()
  result <- tryCatch(
    {
      res <- trust(
        objfun = objfun,
        parinit = x0,
        rinit = 1,
        rmax = 100,
        iterlim = iterlim,
        minimize = TRUE
      )
      # Compute gradient at solution
      final_grad <- funcs$gr(res$argument)
      list(
        par = res$argument,
        value = res$value,
        gradient = final_grad,
        converged = res$converged,
        iterations = res$iterations,
        error = FALSE
      )
    },
    error = function(e) {
      list(error = TRUE, error_msg = conditionMessage(e))
    }
  )
  result$time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  result
}

# Run tests ---------------------------------------------------------------

cat("\n")
cat("=" |> rep(78) |> paste(collapse = ""), "\n")
cat("  Comparative Benchmark: arcopt vs trust\n")
cat("=" |> rep(78) |> paste(collapse = ""), "\n\n")

cat("Testing", length(test_functions), "functions with 10 starting points each\n")
cat("Total test cases:", length(test_functions) * 10 * 2, "(2 optimizers)\n\n")

for (func_name in test_functions) {
  cat("-" |> rep(78) |> paste(collapse = ""), "\n")
  cat(" ", func_name, "\n")
  cat("-" |> rep(78) |> paste(collapse = ""), "\n")

  # Load function
  funcs <- tryCatch(
    load_benchmark_function(func_name, benchmarks_dir),
    error = function(e) {
      cat("  ERROR loading function:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(funcs)) next

  # Get starting points and optimum
  starts <- generate_starting_points(func_name)
  opt <- get_optimum(func_name)

  for (i in seq_along(starts)) {
    start <- starts[[i]]
    x0 <- start$x0

    cat(sprintf("  [%2d/10] %-22s", i, start$description))

    # Run both optimizers
    arc_result <- run_arcopt(x0, funcs, arcopt_control)
    tru_result <- run_trust(x0, funcs, trust_iterlim)

    # Process arcopt result
    if (arc_result$error) {
      arc_status <- "ERR"
      arc_iter <- NA
      arc_err <- NA
      results <- rbind(results, data.frame(
        func_name = func_name, start_id = i, difficulty = start$difficulty,
        description = start$description, optimizer = "arcopt",
        converged = FALSE, iterations = NA_integer_,
        final_grad_norm = NA_real_, solution_error = NA_real_,
        function_error = NA_real_, time_sec = arc_result$time,
        error_msg = arc_result$error_msg, stringsAsFactors = FALSE
      ))
    } else {
      arc_grad <- max(abs(arc_result$gradient))
      arc_err <- sqrt(sum((arc_result$par - opt$x)^2))
      arc_ferr <- abs(arc_result$value - opt$f)
      arc_status <- if (arc_result$converged && arc_err < 1e-4) "OK" else if (arc_result$converged) "WR" else "FL"
      arc_iter <- arc_result$iterations
      results <- rbind(results, data.frame(
        func_name = func_name, start_id = i, difficulty = start$difficulty,
        description = start$description, optimizer = "arcopt",
        converged = arc_result$converged, iterations = arc_iter,
        final_grad_norm = arc_grad, solution_error = arc_err,
        function_error = arc_ferr, time_sec = arc_result$time,
        error_msg = "", stringsAsFactors = FALSE
      ))
    }

    # Process trust result
    if (tru_result$error) {
      tru_status <- "ERR"
      tru_iter <- NA
      tru_err <- NA
      results <- rbind(results, data.frame(
        func_name = func_name, start_id = i, difficulty = start$difficulty,
        description = start$description, optimizer = "trust",
        converged = FALSE, iterations = NA_integer_,
        final_grad_norm = NA_real_, solution_error = NA_real_,
        function_error = NA_real_, time_sec = tru_result$time,
        error_msg = tru_result$error_msg, stringsAsFactors = FALSE
      ))
    } else {
      tru_grad <- max(abs(tru_result$gradient))
      tru_err <- sqrt(sum((tru_result$par - opt$x)^2))
      tru_ferr <- abs(tru_result$value - opt$f)
      tru_status <- if (tru_result$converged && tru_err < 1e-4) "OK" else if (tru_result$converged) "WR" else "FL"
      tru_iter <- tru_result$iterations
      results <- rbind(results, data.frame(
        func_name = func_name, start_id = i, difficulty = start$difficulty,
        description = start$description, optimizer = "trust",
        converged = tru_result$converged, iterations = tru_iter,
        final_grad_norm = tru_grad, solution_error = tru_err,
        function_error = tru_ferr, time_sec = tru_result$time,
        error_msg = "", stringsAsFactors = FALSE
      ))
    }

    # Print comparison
    arc_iter_str <- if (is.na(arc_iter)) "ERR" else sprintf("%4d", arc_iter)
    tru_iter_str <- if (is.na(tru_iter)) "ERR" else sprintf("%4d", tru_iter)
    arc_err_str <- if (is.na(arc_err)) "   ERR" else sprintf("%.1e", arc_err)
    tru_err_str <- if (is.na(tru_err)) "   ERR" else sprintf("%.1e", tru_err)

    cat(sprintf(
      " arc:%s/%s  tru:%s/%s  err: %s / %s\n",
      arc_status, arc_iter_str, tru_status, tru_iter_str, arc_err_str, tru_err_str
    ))
  }

  cat("\n")
}

# Summary -----------------------------------------------------------------

cat("=" |> rep(78) |> paste(collapse = ""), "\n")
cat("  Summary\n")
cat("=" |> rep(78) |> paste(collapse = ""), "\n\n")

# Overall stats by optimizer
for (optim in c("arcopt", "trust")) {
  opt_results <- results[results$optimizer == optim, ]
  n_total <- nrow(opt_results)
  n_converged <- sum(opt_results$converged, na.rm = TRUE)
  n_correct <- sum(opt_results$converged & opt_results$solution_error < 1e-4, na.rm = TRUE)
  n_errors <- sum(is.na(opt_results$converged) | opt_results$error_msg != "", na.rm = TRUE)

  cat(sprintf("%s:\n", toupper(optim)))
  cat(sprintf("  Converged:          %d/%d (%.1f%%)\n", n_converged, n_total, 100 * n_converged / n_total))
  cat(sprintf("  Correct (err<1e-4): %d/%d (%.1f%%)\n", n_correct, n_total, 100 * n_correct / n_total))
  cat(sprintf("  Errors:             %d\n", n_errors))
  cat("\n")
}

# Per-function comparison
cat("Per-function comparison (Correct / Total):\n")
cat("-" |> rep(78) |> paste(collapse = ""), "\n")
cat(sprintf("%-30s %15s %15s %15s\n", "Function", "arcopt", "trust", "Winner"))
cat("-" |> rep(78) |> paste(collapse = ""), "\n")

for (func_name in test_functions) {
  arc_res <- results[results$func_name == func_name & results$optimizer == "arcopt", ]
  tru_res <- results[results$func_name == func_name & results$optimizer == "trust", ]

  arc_corr <- sum(arc_res$converged & arc_res$solution_error < 1e-4, na.rm = TRUE)
  tru_corr <- sum(tru_res$converged & tru_res$solution_error < 1e-4, na.rm = TRUE)
  n <- nrow(arc_res)

  winner <- if (arc_corr > tru_corr) "arcopt" else if (tru_corr > arc_corr) "trust" else "tie"

  cat(sprintf(
    "%-30s %7d/%-7d %7d/%-7d %15s\n",
    func_name, arc_corr, n, tru_corr, n, winner
  ))
}

cat("-" |> rep(78) |> paste(collapse = ""), "\n")

# Iteration comparison for cases where both found correct solution
cat("\nIteration comparison (where both correct):\n")
cat("-" |> rep(78) |> paste(collapse = ""), "\n")

arc_all <- results[results$optimizer == "arcopt", ]
tru_all <- results[results$optimizer == "trust", ]

# Merge by func_name and start_id
merged <- merge(
  arc_all[, c("func_name", "start_id", "converged", "solution_error", "iterations")],
  tru_all[, c("func_name", "start_id", "converged", "solution_error", "iterations")],
  by = c("func_name", "start_id"),
  suffixes = c("_arc", "_tru")
)

both_correct <- merged[
  merged$converged_arc & merged$converged_tru &
  merged$solution_error_arc < 1e-4 & merged$solution_error_tru < 1e-4,
]

if (nrow(both_correct) > 0) {
  cat(sprintf("Cases where both found correct solution: %d\n", nrow(both_correct)))
  cat(sprintf("  arcopt avg iterations: %.1f\n", mean(both_correct$iterations_arc, na.rm = TRUE)))
  cat(sprintf("  trust  avg iterations: %.1f\n", mean(both_correct$iterations_tru, na.rm = TRUE)))
  cat(sprintf("  arcopt wins (fewer iter): %d\n", sum(both_correct$iterations_arc < both_correct$iterations_tru)))
  cat(sprintf("  trust  wins (fewer iter): %d\n", sum(both_correct$iterations_tru < both_correct$iterations_arc)))
  cat(sprintf("  ties:                     %d\n", sum(both_correct$iterations_arc == both_correct$iterations_tru)))
}

# Save results
results_file <- file.path(script_dir, "results", "comparison_results.csv")
write.csv(results, results_file, row.names = FALSE)
cat("\n")
cat("Results saved to:", results_file, "\n")

cat("\n")
cat("=" |> rep(78) |> paste(collapse = ""), "\n")

invisible(results)
