# Run All Benchmark Function Tests
# ==================================
#
# Runs all test scripts for the 8 benchmark functions

cat("\n")
cat("="[rep(1, 70)], "\n", sep = "")
cat("  Running All Benchmark Function Tests\n")
cat("="[rep(1, 70)], "\n\n", sep = "")

test_functions <- c(
  "sphere",
  "bohachevsky",
  "sum_of_squares",
  "rotated_hyper_ellipsoid",
  "trid",
  "rosenbrock",
  "dixon_price",
  "powell_singular"
)

results <- list()

for (func_name in test_functions) {
  cat("\n")
  cat("-"[rep(1, 70)], "\n", sep = "")
  cat("  ", func_name, "\n", sep = "")
  cat("-"[rep(1, 70)], "\n", sep = "")

  test_file <- file.path(func_name, "test.R")

  tryCatch(
    {
      # Change to the function directory and run the test
      old_wd <- getwd()
      setwd(func_name)
      source("test.R")
      setwd(old_wd)
      results[[func_name]] <- "PASS"
    },
    error = function(e) {
      setwd(old_wd)  # Restore working directory on error
      cat("\nERROR:", conditionMessage(e), "\n")
      results[[func_name]] <- "FAIL"
    }
  )
}

# Summary
cat("\n\n")
cat("="[rep(1, 70)], "\n", sep = "")
cat("  Summary\n")
cat("="[rep(1, 70)], "\n\n", sep = "")

pass_count <- sum(sapply(results, function(x) x == "PASS"))
total_count <- length(results)

for (func_name in test_functions) {
  status <- results[[func_name]]
  symbol <- if (status == "PASS") "✓" else "✗"
  cat(sprintf("  %s %-30s %s\n", symbol, func_name, status))
}

cat("\n")
cat(sprintf("  Total: %d/%d tests passed\n", pass_count, total_count))
cat("\n")
cat("="[rep(1, 70)], "\n", sep = "")

if (pass_count == total_count) {
  cat("  ✓ All benchmark tests passed!\n")
} else {
  cat("  ✗ Some tests failed!\n")
  stop("Test suite failed")
}

cat("="[rep(1, 70)], "\n\n", sep = "")
