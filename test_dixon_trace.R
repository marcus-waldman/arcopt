devtools::load_all()

# Load exact benchmark functions
source("benchmarks/dixon_price/fn.R")
source("benchmarks/dixon_price/gr.R")
source("benchmarks/dixon_price/hess.R")
source("benchmarks/arcopt_tests/starting_points.R")

# Get starting points using exact benchmark method
starts <- generate_dixon_price_starts(seed = 49)

# Get random 4 (index 6)
start_random_4 <- starts[[6]]
x0 <- start_random_4$x0

cat("Dixon-Price Random 4 Trace Analysis\n")
cat("=====================================\n\n")
cat("Starting point:", paste(round(x0, 3), collapse=", "), "\n")
cat("Initial f:", dixon_price_fn(x0), "\n\n")

# Run with trace=2 to see sigma, rho, step_type evolution
result <- arcopt(x0, dixon_price_fn, dixon_price_gr, dixon_price_hess,
  control = list(maxit = 100, use_momentum = FALSE, trace = 2))

cat("Final results:\n")
cat("Converged:", result$converged, "\n")
cat("Iterations:", result$iterations, "\n")
cat("Final f:", result$value, "\n")
cat("Gradient norm:", sqrt(sum(result$gradient^2)), "\n")
cat("Final sigma:", result$sigma, "\n\n")

cat("Sigma evolution (first 20 and last 20):\n")
n_iter <- length(result$trace$sigma)
if (n_iter <= 40) {
  cat(paste(round(result$trace$sigma, 2), collapse=", "), "\n")
} else {
  cat("First 20:", paste(round(result$trace$sigma[1:20], 2), collapse=", "), "\n")
  cat("Last 20:", paste(round(result$trace$sigma[(n_iter-19):n_iter], 2), collapse=", "), "\n")
}

cat("\nStep type distribution:\n")
print(table(result$trace$step_type))

cat("\nRho values (acceptance ratio) - first 20:\n")
cat(paste(round(result$trace$rho[1:min(20, n_iter)], 4), collapse=", "), "\n")

cat("\nReciprocal condition number (first 20 and last 20):\n")
if (n_iter <= 40) {
  cat(paste(format(result$trace$rcond, scientific=TRUE, digits=2), collapse=", "), "\n")
} else {
  cat("First 20:", paste(format(result$trace$rcond[1:20], scientific=TRUE, digits=2), collapse=", "), "\n")
  cat("Last 20:", paste(format(result$trace$rcond[(n_iter-19):n_iter], scientific=TRUE, digits=2), collapse=", "), "\n")
}

# Find when sigma becomes very large
large_sigma_idx <- which(result$trace$sigma > 100)
if (length(large_sigma_idx) > 0) {
  cat("\n\nIterations where sigma > 100:\n")
  cat("Iterations:", paste(large_sigma_idx - 1, collapse=", "), "\n")  # -1 because iter 0 is first
  cat("Sigma values:", paste(round(result$trace$sigma[large_sigma_idx], 2), collapse=", "), "\n")
}

cat("\nTrace analysis complete!\n")
