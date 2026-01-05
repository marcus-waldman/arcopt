devtools::load_all()

# Load exact benchmark functions
source("benchmarks/dixon_price/fn.R")
source("benchmarks/dixon_price/gr.R")
source("benchmarks/dixon_price/hess.R")
source("benchmarks/arcopt_tests/starting_points.R")

# Get starting points
starts <- generate_dixon_price_starts(seed = 49)
x0 <- starts[[6]]$x0

cat("Testing Sigma Reset Fix: After successful Newton, if next Newton fails\n")
cat("(Cholesky failure = indefinite H), reset sigma to sigma0\n")
cat("=======================================================================\n\n")

cat("Starting point:", paste(round(x0, 3), collapse=", "), "\n")
cat("Initial f:", dixon_price_fn(x0), "\n\n")

# Run with trace=2 to see sigma evolution
result <- arcopt(x0, dixon_price_fn, dixon_price_gr, dixon_price_hess,
  control = list(maxit = 100, use_momentum = FALSE, trace = 2))

cat("Final results:\n")
cat("Converged:", result$converged, "\n")
cat("Iterations:", result$iterations, "\n")
cat("Final f:", result$value, "\n")
cat("Gradient norm:", sqrt(sum(result$gradient^2)), "\n")
cat("Final sigma:", format(result$sigma, scientific=TRUE, digits=3), "\n\n")

cat("Sigma evolution (looking for reset to sigma0=1.0):\n")
n_iter <- length(result$trace$sigma)
if (n_iter <= 40) {
  for (i in 1:n_iter) {
    cat(sprintf("Iter %2d: sigma=%8.2e, step_type=%-6s, rho=%7.3f\n",
                i-1, result$trace$sigma[i], result$trace$step_type[i],
                ifelse(is.na(result$trace$rho[i]), 0, result$trace$rho[i])))
  }
} else {
  cat("First 15:\n")
  for (i in 1:15) {
    cat(sprintf("Iter %2d: sigma=%8.2e, step_type=%-6s, rho=%7.3f\n",
                i-1, result$trace$sigma[i], result$trace$step_type[i],
                ifelse(is.na(result$trace$rho[i]), 0, result$trace$rho[i])))
  }
}

# Check if sigma was reset
sigma_vals <- result$trace$sigma
sigma_resets <- which(diff(sigma_vals) > 0.5)  # Sigma increased significantly

if (length(sigma_resets) > 0) {
  cat("\nSIGMA RESETS DETECTED at iterations:", paste(sigma_resets, collapse=", "), "\n")
  cat("Reset values:\n")
  for (idx in sigma_resets) {
    cat(sprintf("  Iter %d: sigma %.2e -> %.2e (reset to sigma0!)\n",
                idx, sigma_vals[idx], sigma_vals[idx+1]))
  }
} else {
  cat("\nNO SIGMA RESETS - sigma only decreased\n")
}

cat("\nStep type distribution:\n")
print(table(result$trace$step_type))

cat("\nComparison:\n")
cat("  OLD (no reset):     Converged=FALSE, Iters=100, f=0.7657, ||g||=0.299\n")
cat("  NEW (with reset):   Converged=", result$converged,
    ", Iters=", result$iterations,
    ", f=", round(result$value, 4),
    ", ||g||=", sprintf("%.3f", sqrt(sum(result$gradient^2))), "\n")

if (result$converged && result$value < 0.7) {
  cat("\n✓ SUCCESS! Sigma reset prevented stagnation!\n")
} else {
  cat("\n✗ Still having issues\n")
}
