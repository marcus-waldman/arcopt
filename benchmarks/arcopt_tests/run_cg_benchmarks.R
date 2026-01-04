# Run arcopt Benchmark Tests with CG Solver
# ==========================================
#
# Tests arcopt on 8 benchmark functions using the CG-Lanczos solver

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

# Load arcopt package
cat("Loading arcopt package...\n")
devtools::load_all(pkg_dir)

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

# Control parameters for arcopt - FORCE CG SOLVER
arcopt_control <- list(
  trace = FALSE,
  maxit = 1000,
  gtol_abs = 1e-6,
  cubic_solver = "cg"  # Force CG solver
)

# Results storage
results <- data.frame(
  func_name = character(),
  start_id = integer(),
  difficulty = character(),
  description = character(),
  converged = logical(),
  iterations = integer(),
  fn_evals = integer(),
  gr_evals = integer(),
  hess_evals = integer(),
  final_grad_norm = numeric(),
  solution_error = numeric(),
  function_error = numeric(),
  message = character(),
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

# Run tests ---------------------------------------------------------------

cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("  arcopt Benchmark Tests - CG SOLVER\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("Testing", length(test_functions), "functions with 10 starting points each\n")
cat("Total test cases:", length(test_functions) * 10, "\n")
cat("Solver: CG-Lanczos (ARCqK)\n\n")

for (func_name in test_functions) {
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")
  cat(" ", func_name, "\n")
  cat("-" |> rep(70) |> paste(collapse = ""), "\n")

  funcs <- tryCatch(
    load_benchmark_function(func_name, benchmarks_dir),
    error = function(e) {
      cat("  ERROR loading function:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(funcs)) next

  starts <- generate_starting_points(func_name)
  opt <- get_optimum(func_name)

  for (i in seq_along(starts)) {
    start <- starts[[i]]
    x0 <- start$x0

    cat(sprintf("  [%2d/10] %-25s", i, start$description))

    start_time <- Sys.time()
    result <- tryCatch(
      {
        arcopt(
          x0 = x0,
          fn = funcs$fn,
          gr = funcs$gr,
          hess = funcs$hess,
          control = arcopt_control
        )
      },
      error = function(e) {
        list(
          error = TRUE,
          error_msg = conditionMessage(e)
        )
      }
    )
    end_time <- Sys.time()
    elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

    if (!is.null(result$error) && result$error) {
      cat(" ERROR:", substr(result$error_msg, 1, 30), "\n")

      results <- rbind(results, data.frame(
        func_name = func_name,
        start_id = i,
        difficulty = start$difficulty,
        description = start$description,
        converged = FALSE,
        iterations = NA_integer_,
        fn_evals = NA_integer_,
        gr_evals = NA_integer_,
        hess_evals = NA_integer_,
        final_grad_norm = NA_real_,
        solution_error = NA_real_,
        function_error = NA_real_,
        message = "ERROR",
        time_sec = elapsed,
        error_msg = result$error_msg,
        stringsAsFactors = FALSE
      ))
    } else {
      final_grad_norm <- max(abs(result$gradient))
      solution_error <- sqrt(sum((result$par - opt$x)^2))
      function_error <- abs(result$value - opt$f)

      status <- if (result$converged) {
        if (solution_error < 1e-4) "OK" else "WRONG"
      } else {
        "FAIL"
      }

      cat(sprintf(
        " %4s  iter=%4d  |g|=%.1e  err=%.1e\n",
        status, result$iterations, final_grad_norm, solution_error
      ))

      results <- rbind(results, data.frame(
        func_name = func_name,
        start_id = i,
        difficulty = start$difficulty,
        description = start$description,
        converged = result$converged,
        iterations = result$iterations,
        fn_evals = result$evaluations$fn,
        gr_evals = result$evaluations$gr,
        hess_evals = result$evaluations$hess,
        final_grad_norm = final_grad_norm,
        solution_error = solution_error,
        function_error = function_error,
        message = result$message,
        time_sec = elapsed,
        error_msg = "",
        stringsAsFactors = FALSE
      ))
    }
  }

  cat("\n")
}

# Summary -----------------------------------------------------------------

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("  Summary - CG Solver\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

n_total <- nrow(results)
n_converged <- sum(results$converged, na.rm = TRUE)
n_correct <- sum(results$converged & results$solution_error < 1e-4, na.rm = TRUE)
n_errors <- sum(results$message == "ERROR", na.rm = TRUE)

cat(sprintf("Total test cases:    %d\n", n_total))
cat(sprintf("Converged:           %d (%.1f%%)\n", n_converged, 100 * n_converged / n_total))
cat(sprintf("Correct (err<1e-4):  %d (%.1f%%)\n", n_correct, 100 * n_correct / n_total))
cat(sprintf("Errors:              %d\n", n_errors))

cat("\n")

cat("Per-function results:\n")
cat("-" |> rep(70) |> paste(collapse = ""), "\n")
cat(sprintf("%-30s %10s %10s %10s\n", "Function", "Converged", "Correct", "Avg Iter"))
cat("-" |> rep(70) |> paste(collapse = ""), "\n")

for (func_name in test_functions) {
  func_results <- results[results$func_name == func_name, ]
  n <- nrow(func_results)
  n_conv <- sum(func_results$converged, na.rm = TRUE)
  n_corr <- sum(func_results$converged & func_results$solution_error < 1e-4, na.rm = TRUE)
  avg_iter <- mean(func_results$iterations, na.rm = TRUE)

  cat(sprintf(
    "%-30s %5d/%d    %5d/%d    %8.1f\n",
    func_name, n_conv, n, n_corr, n, avg_iter
  ))
}

cat("-" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("=" |> rep(70) |> paste(collapse = ""), "\n")

invisible(results)
