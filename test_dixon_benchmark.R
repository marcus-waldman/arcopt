devtools::load_all()

# Load benchmark infrastructure
source("benchmarks/arcopt_tests/starting_points.R")

# Load dixon_price functions
source("benchmarks/dixon_price/fn.R")
source("benchmarks/dixon_price/gr.R")
source("benchmarks/dixon_price/hess.R")

# Get starting points
starts <- generate_dixon_price_starts(seed = 49)

# Get random 4 (index 6: near optimum, at optimum, random 1-6)
start_random_4 <- starts[[6]]
x0 <- start_random_4$x0

cat("Testing dixon_price random 4 via benchmark infrastructure\n")
cat("Starting point:", paste(round(x0, 3), collapse=", "), "\n")
cat("Initial f:", dixon_price_fn(x0), "\n\n")

# Run with same control as benchmark
result <- arcopt(
  x0 = x0,
  fn = dixon_price_fn,
  gr = dixon_price_gr,
  hess = dixon_price_hess,
  control = list(
    trace = FALSE,
    maxit = 1000,
    gtol_abs = 1e-6,
    use_momentum = FALSE
  )
)

cat("Results:\n")
cat("Converged:", result$converged, "\n")
cat("Iterations:", result$iterations, "\n")
cat("Final f:", result$value, "\n")
cat("Gradient norm:", sqrt(sum(result$gradient^2)), "\n")
cat("Max |gradient|:", max(abs(result$gradient)), "\n")
cat("Message:", result$message, "\n")
