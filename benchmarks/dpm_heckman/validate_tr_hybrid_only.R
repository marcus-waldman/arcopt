# Rerun only the hybrid configuration with devtools::load_all so the
# un-installed TR-fallback code is picked up. Trust and nohybrid results
# are already in validate_tr_hybrid_seed101.log from the prior run.

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
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
                      sprintf("validate_tr_hybrid_only_seed%d.log", seed))
try(file.remove(log_path), silent = TRUE)
log_line <- function(msg) {
  cat(msg, "\n", file = log_path, append = TRUE)
  cat(msg, "\n")
}

log_line("Compiling Stan model...")
model_so <- compile_dpm_heckman_model()

log_line(sprintf("Simulating seed=%d  n_obs=%d  K_trunc=%d  K_true=%d",
                 seed, n_obs, k_trunc, k_true))
data <- simulate_dpm_heckman(n_obs = n_obs, k_true = k_true,
                             p_out = p_out, p_sel = p_sel,
                             n_excl = n_excl, seed = seed)
fns <- make_dpm_heckman_fns(data, model_so, k_trunc = k_trunc)
log_line(sprintf("  n_params=%d  sel_rate=%.2f  f(x0)=%.2f  maxit=%d",
                 fns$n_params, data$selection_rate, fns$fn(fns$x0), maxit))

log_line("--- HYBRID (devtools::load_all) ---")
t0 <- proc.time()[["elapsed"]]
fit <- tryCatch(
  arcopt(fns$x0, fns$fn, fns$gr, fns$hess_ad,
         control = list(maxit = maxit, gtol_abs = gtol, trace = 0,
                        tr_fallback_enabled = TRUE)),
  error = function(e) list(error = conditionMessage(e))
)
elapsed <- proc.time()[["elapsed"]] - t0

if (!is.null(fit$error)) {
  log_line(sprintf("ERROR after %.1fs: %s", elapsed, fit$error))
} else {
  g <- fit$gradient
  eig_min <- tryCatch(
    min(eigen(fns$hess_ad(fit$par), symmetric = TRUE,
              only.values = TRUE)$values),
    error = function(e) NA_real_
  )
  status <- if (is.na(eig_min)) "?" else
    if (eig_min < -1e-6) "saddle" else "local min"
  log_line(sprintf("elapsed:           %.2f s", elapsed))
  log_line(sprintf("f:                 %.3f", fit$value))
  log_line(sprintf("|g|_inf:           %.2e", max(abs(g))))
  log_line(sprintf("iterations:        %d", fit$iterations))
  log_line(sprintf("eig_min:           %+.3e", eig_min))
  log_line(sprintf("status:            %s", status))
  log_line(sprintf("solver_mode_final: %s", fit$diagnostics$solver_mode_final))
  log_line(sprintf("ridge_switches:    %d", fit$diagnostics$ridge_switches))
  log_line(sprintf("radius_final:      %s",
                   if (is.na(fit$diagnostics$radius_final)) "NA"
                   else sprintf("%.3e", fit$diagnostics$radius_final)))
  log_line(sprintf("converged:         %s", fit$converged))
  log_line(sprintf("message:           %s", fit$message))
}
