devtools::load_all()

# Simple 2D Rosenbrock
rosenbrock_fn <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
rosenbrock_gr <- function(x) {
  c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
    200 * (x[2] - x[1]^2))
}
rosenbrock_hess <- function(x) {
  matrix(c(
    1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
    -400 * x[1], 200
  ), 2, 2)
}

cat("Testing trace level 0 (no trace):\n")
result0 <- arcopt(
  c(-1.2, 1),
  rosenbrock_fn,
  rosenbrock_gr,
  rosenbrock_hess,
  control = list(trace = 0, maxit = 10)
)
cat("  sigma in result:", !is.null(result0$sigma), "\n")
cat("  trace in result:", !is.null(result0$trace), "\n\n")

cat("Testing trace level 1:\n")
result1 <- arcopt(
  c(-1.2, 1),
  rosenbrock_fn,
  rosenbrock_gr,
  rosenbrock_hess,
  control = list(trace = 1, maxit = 10)
)
cat("  sigma:", result1$sigma, "\n")
cat("  iterations:", result1$iterations, "\n")
cat("  trace fields:", paste(names(result1$trace), collapse=", "), "\n")
cat("  trace$f length:", length(result1$trace$f), "\n")
cat("  trace$f values:", paste(round(result1$trace$f, 3), collapse=", "), "\n\n")

cat("Testing trace level 2:\n")
result2 <- arcopt(
  c(-1.2, 1),
  rosenbrock_fn,
  rosenbrock_gr,
  rosenbrock_hess,
  control = list(trace = 2, maxit = 10)
)
cat("  trace fields:", paste(names(result2$trace), collapse=", "), "\n")
cat("  trace$sigma:", paste(round(result2$trace$sigma, 4), collapse=", "), "\n")
cat("  trace$rho:", paste(round(result2$trace$rho, 4), collapse=", "), "\n")
cat("  trace$step_type:", paste(result2$trace$step_type, collapse=", "), "\n")
cat("  trace$solver_used:", paste(result2$trace$solver_used, collapse=", "), "\n")
cat("  trace$rcond:", paste(round(result2$trace$rcond, 4), collapse=", "), "\n")
cat("  trace$lambda:", paste(round(result2$trace$lambda, 4), collapse=", "), "\n\n")

cat("Testing trace level 3:\n")
result3 <- arcopt(
  c(-1.2, 1),
  rosenbrock_fn,
  rosenbrock_gr,
  rosenbrock_hess,
  control = list(trace = 3, maxit = 5)
)
cat("  trace fields:", paste(names(result3$trace), collapse=", "), "\n")
cat("  trace$x dimensions:", paste(dim(result3$trace$x), collapse="x"), "\n")
cat("  trace$s dimensions:", paste(dim(result3$trace$s), collapse="x"), "\n")
cat("  trace$H dimensions:", paste(dim(result3$trace$H), collapse="x"), "\n")
cat("  converge_criteria length:", length(result3$trace$converge_criteria), "\n\n")

cat("All trace tests completed!\n")
