# TR-fallback hybrid validation on DPM Heckman
# ============================================
#
# Purpose: verify that the new cubic -> TR hybrid (default ON) converges to
# a true local minimum on the DPM Heckman MAP problem that previously
# stalled arcopt-AD at a saddle with ||g||_inf ~ 1.6.
#
# Compares three configurations:
#   - arcopt-AD-hybrid    (tr_fallback_enabled = TRUE,  new default)
#   - arcopt-AD-nohybrid  (tr_fallback_enabled = FALSE, prior behavior)
#   - trust               (reference: matches briefing's winning method)
#
# Usage (from repo root):
#   Rscript benchmarks/dpm_heckman/validate_tr_hybrid.R [seed] [maxit]

suppressPackageStartupMessages({
  # Use the development source tree rather than an installed version of
  # the package, so that un-installed TR-fallback changes are picked up.
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(arcopt)
  }
})

source(file.path("benchmarks", "dpm_heckman", "setup.R"))

args <- commandArgs(trailingOnly = TRUE)
seed <- if (length(args) >= 1L) as.integer(args[[1]]) else 101L
maxit <- if (length(args) >= 2L) as.integer(args[[2]]) else 200L

n_obs <- 500L
k_trunc <- 20L
k_true <- 2L
p_out <- 5L
p_sel <- 7L
n_excl <- 2L
gtol <- 1e-5

log_path <- file.path("benchmarks", "dpm_heckman",
                      sprintf("validate_tr_hybrid_seed%d.log", seed))
try(file.remove(log_path), silent = TRUE)
log_line <- function(...) {
  msg <- paste0(sprintf(...), "\n")
  cat(msg)
  cat(msg, file = log_path, append = TRUE)
}

log_line("Compiling Stan model...")
model_so <- compile_dpm_heckman_model()

log_line("Simulating seed=%d  n_obs=%d  K_trunc=%d  K_true=%d",
         seed, n_obs, k_trunc, k_true)
data <- simulate_dpm_heckman(n_obs = n_obs, k_true = k_true,
                             p_out = p_out, p_sel = p_sel,
                             n_excl = n_excl, seed = seed)
fns <- make_dpm_heckman_fns(data, model_so, k_trunc = k_trunc)
x0 <- fns$x0
f0 <- fns$fn(x0)
log_line("  n_params=%d  sel_rate=%.2f  f(x0)=%.2f  maxit=%d",
         fns$n_params, data$selection_rate, f0, maxit)

run_arcopt <- function(label, hybrid_on) {
  ctrl <- list(maxit = maxit, gtol_abs = gtol, trace = 0L,
               tr_fallback_enabled = hybrid_on)
  log_line("  [running ] %s", label)
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    arcopt::arcopt(x0 = x0, fn = fns$fn, gr = fns$gr,
                   hess = fns$hess_ad, control = ctrl),
    error = function(e) list(error = conditionMessage(e))
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  if (!is.null(fit$error)) {
    log_line("  %s  ERROR: %s", label, fit$error)
    return(list(label = label, error = fit$error))
  }
  g <- fit$gradient
  eig_min <- tryCatch(
    min(eigen(fns$hess_ad(fit$par), symmetric = TRUE,
              only.values = TRUE)$values),
    error = function(e) NA_real_
  )
  status <- if (is.na(eig_min)) "?" else
    if (eig_min < -1e-6) "saddle" else "local min"
  log_line(paste("  %-22s  f=%10.3f  |g|_inf=%7.2e  iter=%3d  time=%7.2fs",
                 "eig_min=%+9.2e  mode=%s  switches=%d  conv=%s  %s"),
           label, fit$value, max(abs(g)), fit$iterations, elapsed,
           eig_min, fit$solver_mode_final, fit$ridge_switches,
           isTRUE(fit$converged) && max(abs(g)) < gtol, status)
  list(label = label, fit = fit, elapsed = elapsed,
       eig_min = eig_min, status = status)
}

run_trust <- function() {
  if (!requireNamespace("trust", quietly = TRUE)) {
    log_line("  trust package not installed; skipping reference")
    return(NULL)
  }
  log_line("  [running ] trust")
  t0 <- proc.time()[["elapsed"]]
  objfun <- function(p) list(value = fns$fn(p),
                             gradient = fns$gr(p),
                             hessian = fns$hess_ad(p))
  fit <- tryCatch(
    trust::trust(objfun, x0, rinit = 1, rmax = 100, iterlim = maxit),
    error = function(e) list(error = conditionMessage(e))
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  if (!is.null(fit$error)) {
    log_line("  trust  ERROR: %s", fit$error)
    return(NULL)
  }
  g <- fns$gr(fit$argument)
  eig_min <- tryCatch(
    min(eigen(fns$hess_ad(fit$argument), symmetric = TRUE,
              only.values = TRUE)$values),
    error = function(e) NA_real_
  )
  status <- if (is.na(eig_min)) "?" else
    if (eig_min < -1e-6) "saddle" else "local min"
  log_line("  %-22s  f=%10.3f  |g|_inf=%7.2e  iter=%3d  time=%7.2fs  eig_min=%+9.2e  conv=%s  %s",
           "trust", fit$value, max(abs(g)), fit$iterations, elapsed,
           eig_min, fit$converged && max(abs(g)) < gtol, status)
  list(label = "trust", fit = fit, elapsed = elapsed,
       eig_min = eig_min, status = status)
}

# Run all three
res_trust <- run_trust()
res_nohybrid <- run_arcopt("arcopt-AD-nohybrid", hybrid_on = FALSE)
res_hybrid <- run_arcopt("arcopt-AD-hybrid", hybrid_on = TRUE)

# Verdict
log_line("")
log_line("===== VERDICT =====")
if (!is.null(res_trust) && !is.null(res_trust$fit)) {
  log_line("trust:    converged to eig_min=%+.3e at f=%.3f in %d iters",
           res_trust$eig_min, res_trust$fit$value, res_trust$fit$iterations)
}
if (!is.null(res_nohybrid$fit)) {
  log_line("nohybrid: eig_min=%+.3e  f=%.3f  iters=%d  |g|_inf=%.2e  (%s)",
           res_nohybrid$eig_min, res_nohybrid$fit$value,
           res_nohybrid$fit$iterations,
           max(abs(res_nohybrid$fit$gradient)),
           res_nohybrid$status)
}
if (!is.null(res_hybrid$fit)) {
  log_line("hybrid:   eig_min=%+.3e  f=%.3f  iters=%d  |g|_inf=%.2e  (%s)  mode=%s  switches=%d",
           res_hybrid$eig_min, res_hybrid$fit$value,
           res_hybrid$fit$iterations,
           max(abs(res_hybrid$fit$gradient)),
           res_hybrid$status,
           res_hybrid$fit$solver_mode_final,
           res_hybrid$fit$ridge_switches)
}

# Success criteria (briefing §8 point 4):
#   hybrid reaches local min (eig_min > 0) in <= 2x trust iter count
if (!is.null(res_trust) && !is.null(res_trust$fit) &&
    !is.null(res_hybrid$fit)) {
  reached_local_min <- isTRUE(res_hybrid$eig_min > 0)
  within_budget <- res_hybrid$fit$iterations <= 2L * res_trust$fit$iterations
  log_line("")
  log_line("Success criterion: reached_local_min=%s  within_2x_iters=%s",
           reached_local_min, within_budget)
  if (reached_local_min && within_budget) {
    log_line("  >> PASS <<")
  } else {
    log_line("  >> FAIL <<")
  }
}

invisible(NULL)
