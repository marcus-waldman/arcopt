# NIST StRD higher-difficulty problems (those shipped in NISTnls)
# ===============================================================
#
# 10 of the 12 official higher-difficulty NLS problems are available:
# Bennett5, Eckerle4, Hahn1, Lanczos1, Lanczos2, Lanczos3, MGH09, MGH10,
# MGH17, Thurber. (BoxBOD, Rat42, Rat43 are not in NISTnls.)
#
# Each problem entry exposes:
#   mu(par, x)        predicted mean as a function of parameters and predictor
#   start_far         NIST Start I (further from optimum)
#   start_near        NIST Start II (closer to optimum)
#   data_name         NISTnls dataset name (data has columns x, y)
#   par_names         names() of the parameter vector
#
# Objective is least-squares SSR: f(par) = sum_i (y_i - mu(par, x_i))^2.
# Gradient and Hessian are formed from numDeriv on f, since the problems
# are small (n_par <= 7, n_obs <= 154) and we want one-line definitions.

suppressWarnings(library(NISTnls))

nist_problems <- list()

# ---- Bennett5: y = b1 * (b2 + x)^(-1 / b3) ----------------------------
nist_problems$Bennett5 <- list(
  mu = function(par, x) par[[1L]] * (par[[2L]] + x)^(-1 / par[[3L]]),
  start_far  = c(b1 = -2000, b2 = 50,  b3 = 0.8),
  start_near = c(b1 = -1500, b2 = 45,  b3 = 0.85),
  data_name = "Bennett5",
  par_names = c("b1", "b2", "b3")
)

# ---- Eckerle4: y = (b1 / b2) * exp(-0.5 * ((x - b3) / b2)^2) ----------
nist_problems$Eckerle4 <- list(
  mu = function(par, x) {
    (par[[1L]] / par[[2L]]) *
      exp(-0.5 * ((x - par[[3L]]) / par[[2L]])^2)
  },
  start_far  = c(b1 = 1, b2 = 10, b3 = 500),
  start_near = c(b1 = 1.5, b2 = 5, b3 = 450),
  data_name = "Eckerle4",
  par_names = c("b1", "b2", "b3")
)

# ---- Hahn1: y = (b1+b2*x+b3*x^2+b4*x^3)/(1+b5*x+b6*x^2+b7*x^3) -------
nist_problems$Hahn1 <- list(
  mu = function(par, x) {
    (par[[1L]] + par[[2L]] * x + par[[3L]] * x^2 + par[[4L]] * x^3) /
      (1 + par[[5L]] * x + par[[6L]] * x^2 + par[[7L]] * x^3)
  },
  start_far  = c(b1 = 10, b2 = -1, b3 = 0.05, b4 = -1e-5,
                 b5 = -0.05, b6 = 0.001, b7 = -1e-6),
  start_near = c(b1 = 1, b2 = -0.1, b3 = 0.005, b4 = -1e-6,
                 b5 = -0.005, b6 = 0.0001, b7 = -1e-7),
  data_name = "Hahn1",
  par_names = c("b1", "b2", "b3", "b4", "b5", "b6", "b7")
)

# ---- Lanczos1/2/3: y = b1*exp(-b2*x)+b3*exp(-b4*x)+b5*exp(-b6*x) ------
make_lanczos <- function(data_name) {
  list(
    mu = function(par, x) {
      par[[1L]] * exp(-par[[2L]] * x) +
        par[[3L]] * exp(-par[[4L]] * x) +
        par[[5L]] * exp(-par[[6L]] * x)
    },
    start_far  = c(b1 = 1.2, b2 = 0.3, b3 = 5.6, b4 = 5.5,
                   b5 = 6.5, b6 = 7.6),
    start_near = c(b1 = 0.5, b2 = 0.7, b3 = 3.6, b4 = 4.2,
                   b5 = 4,   b6 = 6.3),
    data_name = data_name,
    par_names = c("b1", "b2", "b3", "b4", "b5", "b6")
  )
}
nist_problems$Lanczos1 <- make_lanczos("Lanczos1")
nist_problems$Lanczos2 <- make_lanczos("Lanczos2")
nist_problems$Lanczos3 <- make_lanczos("Lanczos3")

# ---- MGH09: y = b1*(x^2 + x*b2)/(x^2 + x*b3 + b4) ---------------------
nist_problems$MGH09 <- list(
  mu = function(par, x) {
    par[[1L]] * (x^2 + x * par[[2L]]) /
      (x^2 + x * par[[3L]] + par[[4L]])
  },
  start_far  = c(b1 = 25, b2 = 39, b3 = 41.5, b4 = 39),
  start_near = c(b1 = 0.25, b2 = 0.39, b3 = 0.415, b4 = 0.39),
  data_name = "MGH09",
  par_names = c("b1", "b2", "b3", "b4")
)

# ---- MGH10: y = b1 * exp(b2 / (x + b3)) -------------------------------
nist_problems$MGH10 <- list(
  mu = function(par, x) par[[1L]] * exp(par[[2L]] / (x + par[[3L]])),
  start_far  = c(b1 = 2,    b2 = 400000, b3 = 25000),
  start_near = c(b1 = 0.02, b2 = 4000,   b3 = 250),
  data_name = "MGH10",
  par_names = c("b1", "b2", "b3")
)

# ---- MGH17: y = b1 + b2*exp(-x*b4) + b3*exp(-x*b5) --------------------
nist_problems$MGH17 <- list(
  mu = function(par, x) {
    par[[1L]] +
      par[[2L]] * exp(-x * par[[4L]]) +
      par[[3L]] * exp(-x * par[[5L]])
  },
  start_far  = c(b1 = 50, b2 = 150, b3 = -100, b4 = 1,   b5 = 2),
  start_near = c(b1 = 0.5, b2 = 1.5, b3 = -1,   b4 = 0.01, b5 = 0.02),
  data_name = "MGH17",
  par_names = c("b1", "b2", "b3", "b4", "b5")
)

# ---- Thurber: y = (b1+b2*x+b3*x^2+b4*x^3)/(1+b5*x+b6*x^2+b7*x^3) ------
nist_problems$Thurber <- list(
  mu = function(par, x) {
    (par[[1L]] + par[[2L]] * x +
       par[[3L]] * x^2 + par[[4L]] * x^3) /
      (1 + par[[5L]] * x + par[[6L]] * x^2 + par[[7L]] * x^3)
  },
  start_far  = c(b1 = 1000, b2 = 1000, b3 = 400, b4 = 40,
                 b5 = 0.7,  b6 = 0.3,  b7 = 0.03),
  start_near = c(b1 = 1300, b2 = 1500, b3 = 500, b4 = 75,
                 b5 = 1,    b6 = 0.4,  b7 = 0.05),
  data_name = "Thurber",
  par_names = c("b1", "b2", "b3", "b4", "b5", "b6", "b7")
)

# ---- helpers ----------------------------------------------------------

#' Load the (x, y) numeric vectors for a NIST problem.
load_nist_xy <- function(name) {
  e <- new.env()
  do.call(data, list(name, package = "NISTnls", envir = e))
  d <- get(name, envir = e)
  list(x = d$x, y = d$y)
}

#' Build (fn, gr, hess) for arcopt from a problem entry.
make_nist_fns <- function(problem) {
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("package 'numDeriv' is required")
  }
  xy <- load_nist_xy(problem$data_name)
  x <- xy$x
  y <- xy$y
  par_names <- problem$par_names

  fn <- function(par) {
    pp <- stats::setNames(par, par_names)
    res <- y - problem$mu(pp, x)
    sum(res * res)
  }
  gr <- function(par) {
    numDeriv::grad(fn, par, method = "Richardson")
  }
  hess <- function(par) {
    h_mat <- numDeriv::jacobian(gr, par, method = "Richardson")
    (h_mat + t(h_mat)) / 2
  }
  list(fn = fn, gr = gr, hess = hess, x = x, y = y)
}
