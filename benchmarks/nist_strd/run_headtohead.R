# Head-to-head NIST StRD comparison
# =================================
#
# For each problem, run four methods at TWO starts (Start I = far,
# Start II = near):
#   - arcopt
#   - nls(port)            (Levenberg-Marquardt-style globalization)
#   - nls() default        (Gauss-Newton)
#
# Reference SSR = minimum SSR achieved by any (method, start) combo.
# Pass = method's SSR within 1% of the best-known SSR (or absolute
# 1e-6 if SSR is near zero).
#
# This avoids hardcoding certified values: if any of the three methods
# from any starting point lands in the global basin, that defines the
# reference.

source("benchmarks/nist_strd/problems.R")
devtools::load_all(".", quiet = TRUE)

nls_formula <- function(problem) {
  body_text <- switch(problem$data_name,
    Bennett5 = "y ~ b1 * (b2 + x)^(-1/b3)",
    Eckerle4 = "y ~ (b1/b2) * exp(-0.5 * ((x - b3)/b2)^2)",
    Hahn1    = "y ~ (b1 + b2*x + b3*x^2 + b4*x^3) / (1 + b5*x + b6*x^2 + b7*x^3)",
    Lanczos1 = ,
    Lanczos2 = ,
    Lanczos3 = "y ~ b1*exp(-b2*x) + b3*exp(-b4*x) + b5*exp(-b6*x)",
    MGH09    = "y ~ b1 * (x^2 + x*b2) / (x^2 + x*b3 + b4)",
    MGH10    = "y ~ b1 * exp(b2 / (x + b3))",
    MGH17    = "y ~ b1 + b2*exp(-x*b4) + b3*exp(-x*b5)",
    Thurber  = "y ~ (b1 + b2*x + b3*x^2 + b4*x^3) / (1 + b5*x + b6*x^2 + b7*x^3)",
    stop("unknown: ", problem$data_name)
  )
  stats::as.formula(body_text)
}

run_nls <- function(problem, start_vec, alg = c("port", "default")) {
  alg <- match.arg(alg)
  d <- as.data.frame(load_nist_xy(problem$data_name))
  start_list <- as.list(start_vec)
  form <- nls_formula(problem)
  ctrl <- stats::nls.control(maxiter = 500L, warnOnly = TRUE,
                             minFactor = 1e-12)
  args <- list(formula = form, data = d, start = start_list, control = ctrl)
  if (alg == "port") args$algorithm <- "port"

  fit <- tryCatch(do.call(stats::nls, args),
                  error = function(e) NULL,
                  warning = function(w) NULL)
  if (is.null(fit)) {
    return(list(par = rep(NA_real_, length(start_vec)), ssr = NA_real_,
                converged = FALSE, message = "ERROR/WARN"))
  }
  list(par = stats::coef(fit),
       ssr = sum(stats::residuals(fit)^2),
       converged = fit$convInfo$isConv,
       message = fit$convInfo$stopMessage)
}

run_arcopt <- function(problem, start_vec, use_qn = FALSE) {
  fns <- make_nist_fns(problem)
  x0 <- unname(start_vec)
  ctrl <- list(maxit = 500L, gtol_abs = 1e-8)
  hess_arg <- fns$hess
  if (use_qn) {
    ctrl$use_qn <- TRUE
    hess_arg <- NULL                     # QN seeds B_0 from FD on gr
  }
  res <- tryCatch(
    arcopt(x0, fn = fns$fn, gr = fns$gr, hess = hess_arg, control = ctrl),
    error = function(e) NULL)
  if (is.null(res)) {
    return(list(par = rep(NA_real_, length(x0)), ssr = NA_real_,
                converged = FALSE, message = "ERROR",
                iter = NA_integer_, hess = NA_integer_))
  }
  list(par = stats::setNames(res$par, problem$par_names),
       ssr = res$value,
       converged = isTRUE(res$converged) && res$message != "max_iter",
       message = res$message,
       iter = res$iterations,
       hess = res$evaluations$hess)
}

# ---- main loop ----

prob_names <- names(nist_problems)
runs <- list()

for (nm in prob_names) {
  p <- nist_problems[[nm]]
  for (start_label in c("far", "near")) {
    start_vec <- if (start_label == "far") p$start_far else p$start_near
    r_arc    <- run_arcopt(p, start_vec, use_qn = FALSE)
    r_arc_qn <- run_arcopt(p, start_vec, use_qn = TRUE)
    r_port   <- run_nls(p, start_vec, "port")
    r_gn     <- run_nls(p, start_vec, "default")

    runs[[length(runs) + 1L]] <- data.frame(
      problem = nm, start = start_label, method = "arcopt",
      converged = r_arc$converged, ssr = r_arc$ssr,
      iter = r_arc$iter, hess = r_arc$hess,
      message = r_arc$message, stringsAsFactors = FALSE)
    runs[[length(runs) + 1L]] <- data.frame(
      problem = nm, start = start_label, method = "arcopt-qn",
      converged = r_arc_qn$converged, ssr = r_arc_qn$ssr,
      iter = r_arc_qn$iter, hess = r_arc_qn$hess,
      message = r_arc_qn$message, stringsAsFactors = FALSE)
    runs[[length(runs) + 1L]] <- data.frame(
      problem = nm, start = start_label, method = "nls-port",
      converged = r_port$converged, ssr = r_port$ssr,
      iter = NA_integer_, hess = NA_integer_,
      message = r_port$message, stringsAsFactors = FALSE)
    runs[[length(runs) + 1L]] <- data.frame(
      problem = nm, start = start_label, method = "nls-gn",
      converged = r_gn$converged, ssr = r_gn$ssr,
      iter = NA_integer_, hess = NA_integer_,
      message = r_gn$message, stringsAsFactors = FALSE)
  }
}
runs <- do.call(rbind, runs)

# ---- compute per-problem best SSR + pass flag ----

best_by_problem <- tapply(runs$ssr, runs$problem,
                          function(z) min(z, na.rm = TRUE))
runs$best_ssr <- best_by_problem[runs$problem]

# Pass: SSR within 1% of best (relative), or 1e-6 absolute when best is tiny.
runs$rel_excess <- (runs$ssr - runs$best_ssr) / pmax(abs(runs$best_ssr), 1e-8)
runs$pass <- !is.na(runs$ssr) & runs$rel_excess < 0.01

# ---- summary table ----

dir.create("benchmarks/nist_strd/results", showWarnings = FALSE,
           recursive = TRUE)
utils::write.csv(runs,
                 "benchmarks/nist_strd/results/headtohead.csv",
                 row.names = FALSE)

cat("\n========== Pass rates by (method, start) ==========\n")
agg <- aggregate(pass ~ method + start, data = runs,
                 FUN = function(z) sum(z, na.rm = TRUE))
agg$pct <- round(100 * agg$pass / 10, 1)
print(agg, row.names = FALSE)

cat("\n========== Per-problem detail ==========\n")
# Wide table: row=problem, columns are (start,method) pass flags
runs$key <- paste(runs$start, runs$method, sep = ".")
wide <- reshape(runs[, c("problem", "key", "pass")],
                idvar = "problem", timevar = "key", direction = "wide")
names(wide) <- gsub("pass.", "", names(wide))
print(wide, row.names = FALSE)

cat("\n========== best SSR per problem ==========\n")
ssr_table <- aggregate(ssr ~ problem, data = runs[!is.na(runs$ssr), ],
                       FUN = min)
print(ssr_table, row.names = FALSE)

# go/no-go for arcopt and arcopt-qn
for (m in c("arcopt", "arcopt-qn")) {
  for (s in c("far", "near")) {
    n <- sum(runs$pass[runs$method == m & runs$start == s])
    cat(sprintf("%-9s  %s-start pass: %d / 10\n", m, s, n))
  }
}
