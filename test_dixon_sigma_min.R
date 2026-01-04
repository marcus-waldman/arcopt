devtools::load_all()

# Load exact benchmark functions
source("benchmarks/dixon_price/fn.R")
source("benchmarks/dixon_price/gr.R")
source("benchmarks/dixon_price/hess.R")
source("benchmarks/arcopt_tests/starting_points.R")

# Get starting points using exact benchmark method
starts <- generate_dixon_price_starts(seed = 49)
start_random_4 <- starts[[6]]
x0 <- start_random_4$x0

cat("Dixon-Price Random 4 with NEW sigma_min = 1e-6\n")
cat("===============================================\n\n")
cat("Starting point:", paste(round(x0, 3), collapse=", "), "\n")
cat("Initial f:", dixon_price_fn(x0), "\n\n")

# Run with new sigma_min (default is now 1e-6)
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
  cat(paste(format(result$trace$sigma, scientific=TRUE, digits=3), collapse=", "), "\n")
} else {
  cat("First 20:", paste(format(result$trace$sigma[1:20], scientific=TRUE, digits=3), collapse=", "), "\n")
  cat("Last 20:", paste(format(result$trace$sigma[(n_iter-19):n_iter], scientific=TRUE, digits=3), collapse=", "), "\n")
}

cat("\nDid sigma hit the floor (sigma_min = 1e-6)?\n")
min_sigma <- min(result$trace$sigma)
cat("  Minimum sigma:", format(min_sigma, scientific=TRUE, digits=3), "\n")
if (abs(min_sigma - 1e-6) < 1e-10) {
  cat("  YES - sigma hit sigma_min and stayed there\n")
} else {
  cat("  NO - sigma remained above sigma_min\n")
}

cat("\nStep type distribution:\n")
print(table(result$trace$step_type))

cat("\nComparison with OLD sigma_min = 1e-16:\n")
cat("  OLD: Converged=FALSE, Iterations=100, Final f=0.7657, Final sigma=1e-16\n")
cat("  NEW: Converged=", result$converged, ", Iterations=", result$iterations,
    ", Final f=", round(result$value, 4), ", Final sigma=", format(result$sigma, scientific=TRUE, digits=3), "\n")

if (result$converged) {
  cat("\n SUCCESS! New sigma_min prevents stagnation!\n")
} else {
  cat("\n Still stagnating - may need different approach\n")
}
