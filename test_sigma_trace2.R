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

cat("Starting point:", paste(round(x0, 3), collapse=", "), "\n")
cat("Initial f:", dixon_price_fn(x0), "\n\n")

# Run with trace=1
result <- arcopt(x0, dixon_price_fn, dixon_price_gr, dixon_price_hess,
  control = list(maxit = 100, use_momentum = FALSE, trace = 1))

cat("\nFinal results:\n")
cat("Converged:", result$converged, "\n")
cat("Iterations:", result$iterations, "\n")
cat("Final f:", result$value, "\n")
cat("Gradient norm:", sqrt(sum(result$gradient^2)), "\n")
cat("Final sigma:", result$sigma, "\n")
