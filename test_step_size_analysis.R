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

cat("Step Size Overshoot Analysis for dixon_price random 4\n")
cat("======================================================\n\n")

# Compute step norms from lambda = sigma * ||s||
step_norms <- result$trace$lambda / result$trace$sigma
step_norms[1] <- NA  # First iteration is init

cat("Step size statistics (first 30 iterations):\n")
valid_steps <- step_norms[2:min(31, length(step_norms))]
cat("  Min ||s||:", sprintf("%.4f", min(valid_steps, na.rm=TRUE)), "\n")
cat("  Max ||s||:", sprintf("%.4f", max(valid_steps, na.rm=TRUE)), "\n")
cat("  Mean ||s||:", sprintf("%.4f", mean(valid_steps, na.rm=TRUE)), "\n")
cat("  Median ||s||:", sprintf("%.4f", median(valid_steps, na.rm=TRUE)), "\n\n")

# Compute function value changes
f_changes <- diff(result$trace$f)

cat("Function value changes (first 20):\n")
cat("Iter | f_k        | Δf         | ||s||    | rho    | sigma    | Overshoot?\n")
cat("-----|------------|------------|----------|--------|----------|------------\n")
for (i in 2:min(21, length(result$trace$f))) {
  iter <- i - 1
  f_val <- result$trace$f[i]
  delta_f <- f_changes[i-1]
  s_norm <- step_norms[i]
  rho_val <- result$trace$rho[i]
  sigma_val <- result$trace$sigma[i]

  # Check for overshoot indicators:
  # 1. rho >> 1 (much better than predicted)
  # 2. Large step size
  # 3. Function decreasing rapidly
  overshoot_flag <- ""
  if (!is.na(rho_val) && rho_val > 1.5) {
    overshoot_flag <- "rho>>1"
  }
  if (!is.na(s_norm) && s_norm > 5) {
    overshoot_flag <- paste(overshoot_flag, "large_s")
  }

  cat(sprintf("%4d | %10.4f | %10.4f | %8.4f | %6.3f | %8.2e | %s\n",
              iter, f_val, delta_f, s_norm, rho_val, sigma_val, overshoot_flag))
}

cat("\n\nOscillation detection (looking for f going up/down):\n")
sign_changes <- sign(f_changes[-1]) != sign(f_changes[-length(f_changes)])
oscillations <- sum(sign_changes, na.rm = TRUE)
cat("  Number of f direction changes:", oscillations, "\n")
if (oscillations > length(f_changes) * 0.3) {
  cat("  WARNING: High oscillation rate (>30%)\n")
}

cat("\n\nProgress after iteration 8 (when Hessian became indefinite):\n")
if (length(result$trace$f) > 20) {
  f_later <- result$trace$f[20:min(30, length(result$trace$f))]
  cat("  f values (iter 19-29):", paste(round(f_later, 4), collapse=", "), "\n")
  cat("  Range:", sprintf("%.4f to %.4f", min(f_later), max(f_later)), "\n")
  cat("  Variance:", sprintf("%.6f", var(f_later)), "\n")

  if (var(f_later) < 0.001) {
    cat("  STAGNATION DETECTED: f barely changing after iter 8\n")
  }
}

cat("\n\nConclusion:\n")
if (any(step_norms > 10, na.rm = TRUE)) {
  cat("- Large steps detected (||s|| > 10)\n")
}
if (any(result$trace$rho > 2, na.rm = TRUE)) {
  cat("- Very optimistic rho values (> 2) indicate model underestimates reduction\n")
}
if (all(f_changes[1:20] < 0)) {
  cat("- Function monotonically decreasing (no overshooting in f values)\n")
} else {
  cat("- Function has some increases (potential overshooting)\n")
}
