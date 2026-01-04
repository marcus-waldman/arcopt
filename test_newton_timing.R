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

# Run with trace=2
result <- arcopt(x0, dixon_price_fn, dixon_price_gr, dixon_price_hess,
  control = list(maxit = 100, use_momentum = FALSE, trace = 2))

# Find Newton step iterations
newton_iters <- which(result$trace$step_type == "newton")

cat("Newton Step Analysis for dixon_price random 4\n")
cat("==============================================\n\n")

cat("Total iterations:", result$iterations, "\n")
cat("Newton steps taken:", length(newton_iters), "\n")
cat("Newton iterations (0-indexed):", paste(newton_iters - 1, collapse=", "), "\n\n")

if (length(newton_iters) > 0) {
  cat("Newton step details:\n")
  cat("Iter | Sigma      | rcond      | rho      | Accepted?\n")
  cat("-----|------------|------------|----------|----------\n")
  for (idx in newton_iters) {
    iter <- idx - 1  # 0-indexed
    sigma_val <- result$trace$sigma[idx]
    rcond_val <- result$trace$rcond[idx]
    rho_val <- result$trace$rho[idx]
    accepted <- if (is.na(rho_val)) "init" else if (rho_val >= 0.1) "YES" else "NO"

    cat(sprintf("%4d | %10.2e | %10.2e | %8.4f | %s\n",
                iter, sigma_val, rcond_val, rho_val, accepted))
  }
}

cat("\n\nProblem: Are we attempting Newton when sigma is too small?\n")
cat("Sigma at Newton attempts:\n")
cat("  Min:", min(result$trace$sigma[newton_iters]), "\n")
cat("  Max:", max(result$trace$sigma[newton_iters]), "\n")
cat("  Median:", median(result$trace$sigma[newton_iters]), "\n")

cat("\nrcond at Newton attempts:\n")
cat("  Min:", min(result$trace$rcond[newton_iters]), "\n")
cat("  Max:", max(result$trace$rcond[newton_iters]), "\n")
cat("  Mean:", mean(result$trace$rcond[newton_iters]), "\n")

cat("\nConclusion:\n")
if (any(result$trace$sigma[newton_iters] < 1e-6)) {
  cat("WARNING: Newton attempted when sigma < 1e-6 (essentially zero)\n")
  cat("This may indicate premature Newton steps before adequate regularization.\n")
}
if (any(result$trace$rcond[newton_iters] == 0)) {
  cat("WARNING: Newton attempted with indefinite Hessian (rcond=0)\n")
  cat("Newton step should fail Cholesky factorization on indefinite Hessian.\n")
}
