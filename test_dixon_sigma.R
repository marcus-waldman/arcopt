devtools::load_all()

# Dixon-Price function
dixon_price <- function(x) {
  n <- length(x)
  (x[1] - 1)^2 + sum(seq(2, n) * (2*x[2:n]^2 - x[1:(n-1)])^2)
}
dixon_price_gr <- function(x) {
  n <- length(x)
  g <- numeric(n)
  g[1] <- 2*(x[1] - 1) - 2*2*(2*x[2]^2 - x[1])
  for (i in 2:(n-1)) {
    g[i] <- 2*i*4*x[i]*(2*x[i]^2 - x[i-1]) - 2*(i+1)*(2*x[i+1]^2 - x[i])
  }
  g[n] <- 2*n*4*x[n]*(2*x[n]^2 - x[n-1])
  g
}
dixon_price_hess <- function(x) {
  n <- length(x)
  H <- matrix(0, n, n)
  for (i in 1:n) {
    if (i == 1) {
      H[1, 1] <- 2 + 4
      H[1, 2] <- -2*2*4*x[2]
    } else if (i == n) {
      H[n, n] <- 2*n*(12*x[n]^2 - 4*x[n-1])
      H[n, n-1] <- -2*n*4*x[n]
    } else {
      H[i, i] <- 2*i*(12*x[i]^2 - 4*x[i-1]) + 2*(i+1)
      H[i, i-1] <- -2*i*4*x[i]
      H[i, i+1] <- -2*(i+1)*4*x[i+1]
    }
  }
  (H + t(H)) / 2
}

# Random 4 starting point (exact same as benchmark)
set.seed(49)  # From generate_dixon_price_starts()
d <- 5
low <- -10
high <- 10
# Generate points 1-3 to advance RNG, then get point 4
runif(d, low, high)  # random 1
runif(d, low, high)  # random 2
runif(d, low, high)  # random 3
x0 <- runif(d, low, high)  # random 4

cat("Starting point:", paste(round(x0, 3), collapse=", "), "\n")
cat("Initial f:", dixon_price(x0), "\n\n")

result <- arcopt(x0, dixon_price, dixon_price_gr, dixon_price_hess,
  control = list(maxit = 1000, use_momentum = TRUE))

cat("\nFinal results:\n")
cat("Converged:", result$converged, "\n")
cat("Iterations:", result$iterations, "\n")
cat("Final f:", result$value, "\n")
cat("Gradient norm:", sqrt(sum(result$gradient^2)), "\n")
cat("Max |gradient|:", max(abs(result$gradient)), "\n")
cat("Final sigma:", result$sigma, "\n")
