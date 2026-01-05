# QN-ARC Comparison Benchmarks
# =============================
#
# Compares standard ARC (with exact Hessian) vs QN-ARC variants
# to evaluate Hessian evaluation savings.

library(arcopt)

# =============================================================================
# Test Functions with Analytic Hessians
# =============================================================================

# Sphere (convex quadratic)
sphere <- list(
  name = "sphere",
  fn = function(x) sum(x^2),
  gr = function(x) 2 * x,
  hess = function(x) 2 * diag(length(x)),
  x0 = function(n) rep(5, n),
  x_opt = function(n) rep(0, n)
)

# Rosenbrock (ill-conditioned valley)
rosenbrock <- list(
  name = "rosenbrock",
  fn = function(x) {
    n <- length(x)
    sum(100 * (x[2:n] - x[1:(n - 1)]^2)^2 + (1 - x[1:(n - 1)])^2)
  },
  gr = function(x) {
    n <- length(x)
    g <- numeric(n)
    for (i in 1:(n - 1)) {
      g[i] <- g[i] - 400 * x[i] * (x[i + 1] - x[i]^2) - 2 * (1 - x[i])
      g[i + 1] <- g[i + 1] + 200 * (x[i + 1] - x[i]^2)
    }
    g
  },
  hess = function(x) {
    n <- length(x)
    h <- matrix(0, n, n)
    for (i in 1:(n - 1)) {
      h[i, i] <- h[i, i] + 1200 * x[i]^2 - 400 * x[i + 1] + 2
      h[i, i + 1] <- h[i, i + 1] - 400 * x[i]
      h[i + 1, i] <- h[i + 1, i] - 400 * x[i]
      h[i + 1, i + 1] <- h[i + 1, i + 1] + 200
    }
    h
  },
  x0 = function(n) {
    x <- rep(-1.2, n)
    x[seq(2, n, by = 2)] <- 1.0
    x
  },
  x_opt = function(n) rep(1, n)
)

# Quadratic (ill-conditioned)
quadratic <- list(
  name = "quadratic",
  fn = function(x) {
    n <- length(x)
    # Condition number ~ n
    0.5 * sum(seq_len(n) * x^2)
  },
  gr = function(x) {
    n <- length(x)
    seq_len(n) * x
  },
  hess = function(x) {
    n <- length(x)
    diag(seq_len(n))
  },
  x0 = function(n) rep(10, n),
  x_opt = function(n) rep(0, n)
)

# Beale's function (n=2, nonconvex)
beale <- list(
  name = "beale",
  fn = function(x) {
    (1.5 - x[1] * (1 - x[2]))^2 +
      (2.25 - x[1] * (1 - x[2]^2))^2 +
      (2.625 - x[1] * (1 - x[2]^3))^2
  },
  gr = function(x) {
    g <- numeric(2)
    y1 <- 1 - x[2]
    y2 <- 1 - x[2]^2
    y3 <- 1 - x[2]^3
    r1 <- 1.5 - x[1] * y1
    r2 <- 2.25 - x[1] * y2
    r3 <- 2.625 - x[1] * y3
    g[1] <- -2 * (r1 * y1 + r2 * y2 + r3 * y3)
    g[2] <- 2 * x[1] * (r1 + 2 * r2 * x[2] + 3 * r3 * x[2]^2)
    g
  },
  hess = function(x) {
    y1 <- 1 - x[2]
    y2 <- 1 - x[2]^2
    y3 <- 1 - x[2]^3
    r1 <- 1.5 - x[1] * y1
    r2 <- 2.25 - x[1] * y2
    r3 <- 2.625 - x[1] * y3

    h11 <- 2 * (y1^2 + y2^2 + y3^2)
    h12 <- 2 * (y1 + 2 * y2 * x[2] + 3 * y3 * x[2]^2 -
      x[1] * (1 + 2 * y2 + 6 * x[2] * y3))
    h22 <- 2 * x[1]^2 * (1 + 4 * x[2]^2 + 9 * x[2]^4) +
      2 * x[1] * (2 * r2 + 6 * r3 * x[2])

    matrix(c(h11, h12, h12, h22), 2, 2)
  },
  x0 = function(n) c(0, 0),
  x_opt = function(n) c(3, 0.5)
)

# Test problems list
test_problems <- list(sphere, rosenbrock, quadratic, beale)


# =============================================================================
# Benchmark Runner
# =============================================================================

run_single_benchmark <- function(problem, method, n, maxit = 1000) {
  x0 <- problem$x0(n)

  # Configure method
  if (method == "standard") {
    control <- list(trace = 0, maxit = maxit)
    use_qn <- FALSE
    hess <- problem$hess
  } else {
    control <- list(trace = 0, maxit = maxit, use_qn = TRUE, qn_method = method)
    use_qn <- TRUE
    hess <- NULL  # Pure QN - no Hessian
  }

  # Timing
  start_time <- Sys.time()

  result <- tryCatch({
    arcopt(x0, problem$fn, problem$gr, hess,
           control = control)
  }, error = function(e) {
    list(converged = FALSE, iterations = NA, evaluations = list(hess = NA),
         message = e$message)
  })

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Extract metrics
  data.frame(
    problem = problem$name,
    n = n,
    method = method,
    converged = result$converged,
    iterations = result$iterations,
    fn_evals = result$evaluations$fn,
    gr_evals = result$evaluations$gr,
    hess_evals = if (is.null(result$evaluations$hess)) 0
                 else result$evaluations$hess,
    time_sec = elapsed,
    stringsAsFactors = FALSE
  )
}


# =============================================================================
# Main Benchmark Loop
# =============================================================================

run_qn_benchmarks <- function(dimensions = c(2, 5, 10),
                               methods = c("standard", "sr1", "bfgs"),
                               n_reps = 3) {
  results <- list()
  i <- 1

  for (problem in test_problems) {
    dims_to_test <- if (problem$name == "beale") 2 else dimensions

    for (n in dims_to_test) {
      for (method in methods) {
        cat(sprintf("Testing %s (n=%d) with %s...\n",
                    problem$name, n, method))

        for (rep in seq_len(n_reps)) {
          result <- run_single_benchmark(problem, method, n)
          result$rep <- rep
          results[[i]] <- result
          i <- i + 1
        }
      }
    }
  }

  do.call(rbind, results)
}


# =============================================================================
# Summarize Results
# =============================================================================

summarize_results <- function(results) {
  # Aggregate by problem, n, method
  agg <- aggregate(
    cbind(converged, iterations, fn_evals, gr_evals, hess_evals, time_sec) ~
      problem + n + method,
    data = results,
    FUN = mean
  )

  # Compute Hessian savings relative to standard
  for (prob in unique(agg$problem)) {
    for (dim in unique(agg$n[agg$problem == prob])) {
      std_hess <- agg$hess_evals[agg$problem == prob &
                                   agg$n == dim &
                                   agg$method == "standard"]
      if (length(std_hess) > 0 && !is.na(std_hess)) {
        mask <- agg$problem == prob & agg$n == dim
        agg$hess_savings[mask] <- 1 - agg$hess_evals[mask] / std_hess
      }
    }
  }

  agg
}


# =============================================================================
# Run Benchmarks if Executed Directly
# =============================================================================

if (interactive() || !exists("skip_benchmark_run")) {
  cat("Running QN-ARC comparison benchmarks...\n\n")

  results <- run_qn_benchmarks(
    dimensions = c(2, 5, 10),
    methods = c("standard", "sr1", "bfgs"),
    n_reps = 1
  )

  summary <- summarize_results(results)

  cat("\n=== BENCHMARK SUMMARY ===\n\n")
  print(summary, digits = 3)

  cat("\n=== KEY FINDINGS ===\n")
  cat("Hessian evaluations - Standard ARC vs QN-ARC:\n")

  # Compare Hessian evaluations
  for (m in c("sr1", "bfgs")) {
    qn_hess <- sum(summary$hess_evals[summary$method == m], na.rm = TRUE)
    std_hess <- sum(summary$hess_evals[summary$method == "standard"], na.rm = TRUE)
    savings <- (1 - qn_hess / std_hess) * 100
    cat(sprintf("  %s: %.1f%% fewer Hessian evaluations\n", m, savings))
  }
}
