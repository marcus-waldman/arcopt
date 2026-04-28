# §4.5 NIST StRD higher-difficulty go/no-go
# =========================================
#
# For each NIST problem:
#   1. Build a robust reference solution by running nls(algorithm="port")
#      from Start II (the closer NIST starting point). If port fails,
#      fall back to Gauss-Newton nls() from Start II. If both fail, we
#      flag the reference as missing and skip the comparison.
#   2. Run arcopt from Start I (the further NIST starting point).
#   3. Pass criterion: arcopt converges AND
#         max |b_arc - b_ref| / max(|b_ref|, 1e-8) < 1e-3
#      AND |SSR_arc - SSR_ref| / max(|SSR_ref|, 1e-8) < 1e-4.
#
# Go/no-go: arcopt matches the reference on >= 80% of problems (8/10).

source("benchmarks/nist_strd/problems.R")
devtools::load_all(".", quiet = TRUE)

#' Build an nls() formula expression from the problem definition.
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

#' Reference fit via nls() from Start II.
fit_reference <- function(problem) {
  d <- as.data.frame(load_nist_xy(problem$data_name))
  start_list <- as.list(problem$start_near)
  form <- nls_formula(problem)

  port <- tryCatch(stats::nls(form, data = d, start = start_list,
                              algorithm = "port",
                              control = stats::nls.control(maxiter = 500L)),
                   error = function(e) NULL,
                   warning = function(w) NULL)
  if (!is.null(port)) {
    return(list(par = stats::coef(port), ssr = sum(stats::residuals(port)^2),
                method = "nls(port)"))
  }
  gn <- tryCatch(stats::nls(form, data = d, start = start_list,
                            control = stats::nls.control(maxiter = 500L,
                                                         warnOnly = TRUE)),
                 error = function(e) NULL)
  if (!is.null(gn)) {
    return(list(par = stats::coef(gn), ssr = sum(stats::residuals(gn)^2),
                method = "nls(gauss-newton)"))
  }
  list(par = NULL, ssr = NA_real_, method = "FAILED")
}

#' arcopt run from Start I.
fit_arcopt <- function(problem, gtol = 1e-8) {
  fns <- make_nist_fns(problem)
  x0 <- unname(problem$start_far)

  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(
    arcopt(x0, fn = fns$fn, gr = fns$gr, hess = fns$hess,
           control = list(maxit = 1000L, gtol_abs = gtol)),
    error = function(e) {
      list(par = rep(NA_real_, length(x0)), value = NA_real_,
           gradient = NA_real_, converged = FALSE,
           iterations = NA_integer_,
           evaluations = list(fn = NA, gr = NA, hess = NA),
           diagnostics = list(solver_mode_final = "ERROR",
                              ridge_switches = NA),
           message = conditionMessage(e))
    })
  t1 <- proc.time()[["elapsed"]]

  list(par = res$par, ssr = res$value, g_inf = max(abs(res$gradient)),
       iter = res$iterations, hess_evals = res$evaluations$hess,
       converged = isTRUE(res$converged),
       mode = res$diagnostics$solver_mode_final,
       message = res$message, time_s = t1 - t0)
}

# ---- run all problems ----

problems <- nist_problems
prob_names <- names(problems)
cat(sprintf("Running arcopt on %d NIST higher-difficulty problems...\n\n",
            length(prob_names)))

rows <- list()
for (nm in prob_names) {
  p <- problems[[nm]]
  ref <- fit_reference(p)
  arc <- fit_arcopt(p)

  if (is.null(ref$par) || any(is.na(arc$par))) {
    par_err <- NA_real_
    ssr_err <- NA_real_
    pass_par <- FALSE
    pass_ssr <- FALSE
  } else {
    # Reorder ref$par to match arc$par order
    ref_par <- ref$par[p$par_names]
    par_err <- max(abs(arc$par - ref_par) /
                     pmax(abs(ref_par), 1e-8))
    ssr_err <- abs(arc$ssr - ref$ssr) / max(abs(ref$ssr), 1e-8)
    pass_par <- par_err < 1e-3
    pass_ssr <- ssr_err < 1e-4
  }

  pass <- arc$converged && pass_par && pass_ssr

  rows[[nm]] <- data.frame(
    problem = nm,
    n_par = length(p$par_names),
    n_obs = length(arc$par_at_arc <- NULL),  # placeholder; n_obs filled below
    converged = arc$converged,
    iter = arc$iter,
    hess = arc$hess_evals,
    g_inf = arc$g_inf,
    rel_par_err = par_err,
    rel_ssr_err = ssr_err,
    ssr_arc = arc$ssr,
    ssr_ref = ref$ssr,
    ref_method = ref$method,
    pass = pass,
    time_s = arc$time_s,
    stringsAsFactors = FALSE
  )
  rows[[nm]]$n_obs <- length(load_nist_xy(p$data_name)$x)

  cat(sprintf("  %-10s n=%3d  np=%d  iter=%4s  hess=%4s  conv=%s  par_err=%.2e  ssr_err=%.2e  PASS=%s\n",
              nm, rows[[nm]]$n_obs, rows[[nm]]$n_par,
              format(arc$iter, width = 4),
              format(arc$hess_evals, width = 4),
              arc$converged,
              ifelse(is.na(par_err), NA_real_, par_err),
              ifelse(is.na(ssr_err), NA_real_, ssr_err),
              pass))
}

results <- do.call(rbind, rows)
rownames(results) <- NULL

# ---- save & summarize ----

dir.create("benchmarks/nist_strd/results", showWarnings = FALSE,
           recursive = TRUE)
utils::write.csv(results,
                 "benchmarks/nist_strd/results/strd_furthest_start.csv",
                 row.names = FALSE)

n_pass <- sum(results$pass, na.rm = TRUE)
n_total <- nrow(results)

cat("\n====================================================\n")
cat(sprintf("PASS RATE: %d / %d  (%.1f%%)\n",
            n_pass, n_total, 100 * n_pass / n_total))
cat(sprintf("Go/no-go (>= 80%% pass): %s\n",
            if (n_pass / n_total >= 0.80) "PASS" else "FAIL"))
cat("====================================================\n")
