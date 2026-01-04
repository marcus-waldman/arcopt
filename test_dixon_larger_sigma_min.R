devtools::load_all()

# Load exact benchmark functions
source("benchmarks/dixon_price/fn.R")
source("benchmarks/dixon_price/gr.R")
source("benchmarks/dixon_price/hess.R")
source("benchmarks/arcopt_tests/starting_points.R")

# Get starting points
starts <- generate_dixon_price_starts(seed = 49)
x0 <- starts[[6]]$x0

cat("Testing different sigma_min values\n")
cat("==================================\n\n")

for (sigma_min_val in c(1e-4, 1e-3, 0.01, 0.1)) {
  result <- arcopt(x0, dixon_price_fn, dixon_price_gr, dixon_price_hess,
    control = list(maxit = 100, use_momentum = FALSE, trace = 0,
                   sigma_min = sigma_min_val))

  cat(sprintf("sigma_min=%.0e: iters=%3d, f=%.4f, ||g||=%.4f, sigma_final=%.2e\n",
              sigma_min_val, result$iterations, result$value,
              sqrt(sum(result$gradient^2)), result$sigma))
}

cat("\nConclusion:\n")
cat("If larger sigma_min values don't help, this is a fundamental issue\n")
cat("with the problem/algorithm interaction, not just the sigma_min floor.\n")
