# Investigate Rosenbrock saddle point trapping
# =============================================

# Setup
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) > 0) {
  script_dir <- dirname(normalizePath(script_path))
} else {
  script_dir <- getwd()
}
benchmarks_dir <- dirname(script_dir)
pkg_dir <- dirname(benchmarks_dir)

devtools::load_all(pkg_dir)

# Source the benchmark functions
source(file.path(benchmarks_dir, "rosenbrock", "fn.R"))
source(file.path(benchmarks_dir, "rosenbrock", "gr.R"))
source(file.path(benchmarks_dir, "rosenbrock", "hess.R"))

# The classic starting point that fails
x0_classic <- c(-1.2, 1, 1, 1, 1)
x_opt <- rep(1, 5)

cat("=== Rosenbrock Saddle Point Analysis ===\n\n")
cat("Starting point:", x0_classic, "\n")
cat("True optimum:", x_opt, "\n\n")

# Run without momentum
cat("Running WITHOUT momentum...\n")
result_no_mom <- arcopt(
  x0 = x0_classic,
  fn = rosenbrock_fn,
  gr = rosenbrock_gr,
  hess = rosenbrock_hess,
  control = list(trace = FALSE, maxit = 200, gtol_abs = 1e-8)
)

cat("\nWITHOUT MOMENTUM:\n")
cat("  Converged:", result_no_mom$converged, "\n")
cat("  Iterations:", result_no_mom$iterations, "\n")
cat("  Final x:", round(result_no_mom$par, 4), "\n")
cat("  Final f:", result_no_mom$value, "\n")
cat("  Gradient norm:", max(abs(result_no_mom$gradient)), "\n")
cat("  Distance to optimum:", sqrt(sum((result_no_mom$par - x_opt)^2)), "\n\n")

# Check what point it converged to
cat("Analysis of converged point:\n")
cat("  x[1] =", result_no_mom$par[1], "(should be 1)\n")
cat("  This appears to be near (-1, 1, 1, 1, 1) - a local minimum!\n\n")

# Run WITH momentum
cat("Running WITH momentum...\n")
result_mom <- arcopt(
  x0 = x0_classic,
  fn = rosenbrock_fn,
  gr = rosenbrock_gr,
  hess = rosenbrock_hess,
  control = list(
    trace = FALSE,
    maxit = 200,
    gtol_abs = 1e-8,
    use_momentum = TRUE,
    momentum_max = 0.9,
    momentum_c1 = 0.1,
    momentum_c2 = 0.1
  )
)

cat("\nWITH MOMENTUM:\n")
cat("  Converged:", result_mom$converged, "\n")
cat("  Iterations:", result_mom$iterations, "\n")
cat("  Final x:", round(result_mom$par, 4), "\n")
cat("  Final f:", result_mom$value, "\n")
cat("  Gradient norm:", max(abs(result_mom$gradient)), "\n")
cat("  Distance to optimum:", sqrt(sum((result_mom$par - x_opt)^2)), "\n\n")

# Compare
cat("=== COMPARISON ===\n")
cat("                   No Momentum    With Momentum\n")
cat(sprintf("Distance to opt:   %.6f       %.6f\n",
            sqrt(sum((result_no_mom$par - x_opt)^2)),
            sqrt(sum((result_mom$par - x_opt)^2))))
cat(sprintf("Iterations:        %d              %d\n",
            result_no_mom$iterations, result_mom$iterations))
cat(sprintf("Final f:           %.6f       %.6f\n",
            result_no_mom$value, result_mom$value))

# Check Hessian eigenvalues at the converged point
cat("\n=== Hessian Analysis at No-Momentum Solution ===\n")
H_at_sol <- rosenbrock_hess(result_no_mom$par)
eigs <- eigen(H_at_sol, only.values = TRUE)$values
cat("Eigenvalues of Hessian:\n")
print(round(eigs, 4))
if (any(eigs < 0)) {
  cat("NEGATIVE EIGENVALUE DETECTED - this is a SADDLE POINT!\n")
} else {
  cat("All eigenvalues positive - this is a LOCAL MINIMUM\n")
}

cat("\n=== Hessian Analysis at True Optimum ===\n")
H_at_opt <- rosenbrock_hess(x_opt)
eigs_opt <- eigen(H_at_opt, only.values = TRUE)$values
cat("Eigenvalues of Hessian at (1,1,1,1,1):\n")
print(round(eigs_opt, 4))
