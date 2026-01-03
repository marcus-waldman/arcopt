# Starting Point Generators for Benchmark Functions
# =================================================
#
# Generates 10 starting points per function with varying difficulty:
# - 6 random domain samples
# - 2 classic/literature starting points
# - 2 edge cases (near optimum, far/boundary)

# Helper to create a starting point entry
make_start <- function(x0, difficulty, description) {
  list(
    x0 = x0,
    difficulty = difficulty,
    description = description
  )
}

# -----------------------------------------------------------------------------
# Sphere Function (5D)
# Domain: [-5.12, 5.12]^n, Optimum: (0, ..., 0)
# -----------------------------------------------------------------------------
generate_sphere_starts <- function(seed = 42) {
  set.seed(seed)
  d <- 5
  low <- -5.12
  high <- 5.12

  list(
    # Near optimum (easy)
    make_start(rep(0.01, d), "easy", "near optimum"),
    # Center of domain
    make_start(rep(0, d), "easy", "at optimum (sanity check)"),
    # Random points
    make_start(runif(d, low, high), "medium", "random 1"),
    make_start(runif(d, low, high), "medium", "random 2"),
    make_start(runif(d, low, high), "medium", "random 3"),
    make_start(runif(d, low, high), "medium", "random 4"),
    make_start(runif(d, low, high), "medium", "random 5"),
    make_start(runif(d, low, high), "medium", "random 6"),
    # Corner of domain (hard)
    make_start(rep(high, d), "hard", "corner of domain"),
    # Opposite corner
    make_start(rep(low, d), "hard", "opposite corner")
  )
}

# -----------------------------------------------------------------------------
# Sum of Squares Function (5D)
# Domain: [-10, 10]^n, Optimum: (0, ..., 0)
# -----------------------------------------------------------------------------
generate_sum_of_squares_starts <- function(seed = 43) {
  set.seed(seed)
  d <- 5
  low <- -10
  high <- 10

  list(
    make_start(rep(0.01, d), "easy", "near optimum"),
    make_start(rep(0, d), "easy", "at optimum (sanity check)"),
    make_start(runif(d, low, high), "medium", "random 1"),
    make_start(runif(d, low, high), "medium", "random 2"),
    make_start(runif(d, low, high), "medium", "random 3"),
    make_start(runif(d, low, high), "medium", "random 4"),
    make_start(runif(d, low, high), "medium", "random 5"),
    make_start(runif(d, low, high), "medium", "random 6"),
    make_start(rep(high, d), "hard", "corner of domain"),
    make_start(rep(low, d), "hard", "opposite corner")
  )
}

# -----------------------------------------------------------------------------
# Bohachevsky Function (2D only)
# Domain: [-100, 100]^2, Optimum: (0, 0)
# -----------------------------------------------------------------------------
generate_bohachevsky_starts <- function(seed = 44) {
  set.seed(seed)
  d <- 2
  low <- -100
  high <- 100

  list(
    make_start(c(0.01, 0.01), "easy", "near optimum"),
    make_start(c(0, 0), "easy", "at optimum (sanity check)"),
    make_start(runif(d, low, high), "medium", "random 1"),
    make_start(runif(d, low, high), "medium", "random 2"),
    make_start(runif(d, low, high), "medium", "random 3"),
    make_start(runif(d, low, high), "medium", "random 4"),
    make_start(runif(d, low, high), "medium", "random 5"),
    make_start(runif(d, low, high), "medium", "random 6"),
    make_start(c(high, high), "hard", "corner of domain"),
    make_start(c(low, low), "hard", "opposite corner")
  )
}

# -----------------------------------------------------------------------------
# Rotated Hyper-Ellipsoid Function (5D)
# Domain: [-65.536, 65.536]^n, Optimum: (0, ..., 0)
# -----------------------------------------------------------------------------
generate_rotated_hyper_ellipsoid_starts <- function(seed = 45) {
  set.seed(seed)
  d <- 5
  low <- -65.536
  high <- 65.536

  list(
    make_start(rep(0.01, d), "easy", "near optimum"),
    make_start(rep(0, d), "easy", "at optimum (sanity check)"),
    make_start(runif(d, low, high), "medium", "random 1"),
    make_start(runif(d, low, high), "medium", "random 2"),
    make_start(runif(d, low, high), "medium", "random 3"),
    make_start(runif(d, low, high), "medium", "random 4"),
    make_start(runif(d, low, high), "medium", "random 5"),
    make_start(runif(d, low, high), "medium", "random 6"),
    make_start(rep(high, d), "hard", "corner of domain"),
    make_start(rep(low, d), "hard", "opposite corner")
  )
}

# -----------------------------------------------------------------------------
# Trid Function (5D)
# Domain: [-d^2, d^2]^d = [-25, 25]^5, Optimum: x_i = i*(d+1-i)
# For d=5: x* = (5, 8, 9, 8, 5), f* = -50
# Note: Origin is a saddle point
# -----------------------------------------------------------------------------
generate_trid_starts <- function(seed = 46) {

  set.seed(seed)
  d <- 5
  low <- -25
  high <- 25

  # Compute the true optimum for d=5
  x_opt <- sapply(1:d, function(i) i * (d + 1 - i))

  list(
    # Near optimum
    make_start(x_opt + rnorm(d, 0, 0.1), "easy", "near optimum"),
    # At optimum (sanity check)
    make_start(x_opt, "easy", "at optimum (sanity check)"),
    # Random points
    make_start(runif(d, low, high), "medium", "random 1"),
    make_start(runif(d, low, high), "medium", "random 2"),
    make_start(runif(d, low, high), "medium", "random 3"),
    make_start(runif(d, low, high), "medium", "random 4"),
    make_start(runif(d, low, high), "medium", "random 5"),
    make_start(runif(d, low, high), "medium", "random 6"),
    # Origin (saddle point - interesting test)
    make_start(rep(0, d), "hard", "origin (saddle point)"),
    # Corner
    make_start(rep(high, d), "hard", "corner of domain")
  )
}

# -----------------------------------------------------------------------------
# Powell Singular Function (4D only)
# Domain: [-4, 5]^4, Optimum: (0, 0, 0, 0)
# Classic start: (3, -1, 0, 1)
# Note: Hessian is singular at optimum
# -----------------------------------------------------------------------------
generate_powell_singular_starts <- function(seed = 47) {
  set.seed(seed)
  d <- 4
  low <- -4
  high <- 5

  list(
    # Near optimum
    make_start(rep(0.01, d), "easy", "near optimum"),
    # At optimum (sanity check)
    make_start(rep(0, d), "easy", "at optimum (sanity check)"),
    # Classic MGH starting point (hard)
    make_start(c(3, -1, 0, 1), "hard", "classic MGH start"),
    # Random points
    make_start(runif(d, low, high), "medium", "random 1"),
    make_start(runif(d, low, high), "medium", "random 2"),
    make_start(runif(d, low, high), "medium", "random 3"),
    make_start(runif(d, low, high), "medium", "random 4"),
    make_start(runif(d, low, high), "medium", "random 5"),
    make_start(runif(d, low, high), "medium", "random 6"),
    # Corner
    make_start(rep(high, d), "hard", "corner of domain")
  )
}

# -----------------------------------------------------------------------------
# Rosenbrock Function (5D)
# Domain: [-5, 10]^n, Optimum: (1, ..., 1)
# Classic start: (-1.2, 1, 1, ..., 1)
# -----------------------------------------------------------------------------
generate_rosenbrock_starts <- function(seed = 48) {
  set.seed(seed)
  d <- 5
  low <- -5
  high <- 10

  list(
    # Near optimum
    make_start(rep(1, d) + rnorm(d, 0, 0.01), "easy", "near optimum"),
    # At optimum (sanity check)
    make_start(rep(1, d), "easy", "at optimum (sanity check)"),
    # Classic starting point (challenging)
    make_start(c(-1.2, rep(1, d - 1)), "hard", "classic (-1.2, 1, ...)"),
    # All -1.2 (very hard)
    make_start(rep(-1.2, d), "hard", "all -1.2"),
    # Random points
    make_start(runif(d, low, high), "medium", "random 1"),
    make_start(runif(d, low, high), "medium", "random 2"),
    make_start(runif(d, low, high), "medium", "random 3"),
    make_start(runif(d, low, high), "medium", "random 4"),
    make_start(runif(d, low, high), "medium", "random 5"),
    # Far corner
    make_start(rep(high, d), "hard", "corner of domain")
  )
}

# -----------------------------------------------------------------------------
# Dixon-Price Function (5D)
# Domain: [-10, 10]^n, Optimum: x_i = 2^(-(2^i - 2) / 2^i)
# -----------------------------------------------------------------------------
generate_dixon_price_starts <- function(seed = 49) {
  set.seed(seed)
  d <- 5
  low <- -10
  high <- 10

  # Compute the true optimum
  x_opt <- sapply(1:d, function(i) 2^(-(2^i - 2) / 2^i))

  list(
    # Near optimum
    make_start(x_opt + rnorm(d, 0, 0.01), "easy", "near optimum"),
    # At optimum (sanity check)
    make_start(x_opt, "easy", "at optimum (sanity check)"),
    # Random points
    make_start(runif(d, low, high), "medium", "random 1"),
    make_start(runif(d, low, high), "medium", "random 2"),
    make_start(runif(d, low, high), "medium", "random 3"),
    make_start(runif(d, low, high), "medium", "random 4"),
    make_start(runif(d, low, high), "medium", "random 5"),
    make_start(runif(d, low, high), "medium", "random 6"),
    # Corners
    make_start(rep(high, d), "hard", "corner of domain"),
    make_start(rep(low, d), "hard", "opposite corner")
  )
}

# -----------------------------------------------------------------------------
# Main generator function
# -----------------------------------------------------------------------------
generate_starting_points <- function(func_name) {
  switch(func_name,
    "sphere" = generate_sphere_starts(),
    "sum_of_squares" = generate_sum_of_squares_starts(),
    "bohachevsky" = generate_bohachevsky_starts(),
    "rotated_hyper_ellipsoid" = generate_rotated_hyper_ellipsoid_starts(),
    "trid" = generate_trid_starts(),
    "powell_singular" = generate_powell_singular_starts(),
    "rosenbrock" = generate_rosenbrock_starts(),
    "dixon_price" = generate_dixon_price_starts(),
    stop(paste("Unknown function:", func_name))
  )
}

# -----------------------------------------------------------------------------
# Get optimum for each function (for computing solution error)
# -----------------------------------------------------------------------------
get_optimum <- function(func_name) {
  switch(func_name,
    "sphere" = list(x = rep(0, 5), f = 0),
    "sum_of_squares" = list(x = rep(0, 5), f = 0),
    "bohachevsky" = list(x = c(0, 0), f = 0),
    "rotated_hyper_ellipsoid" = list(x = rep(0, 5), f = 0),
    "trid" = {
      d <- 5
      x <- sapply(1:d, function(i) i * (d + 1 - i))
      f <- -d * (d + 4) * (d - 1) / 6  # = -50 for d=5
      list(x = x, f = f)
    },
    "powell_singular" = list(x = rep(0, 4), f = 0),
    "rosenbrock" = list(x = rep(1, 5), f = 0),
    "dixon_price" = {
      d <- 5
      x <- sapply(1:d, function(i) 2^(-(2^i - 2) / 2^i))
      list(x = x, f = 0)
    },
    stop(paste("Unknown function:", func_name))
  )
}
