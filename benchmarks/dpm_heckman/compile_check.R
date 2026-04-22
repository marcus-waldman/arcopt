library(bridgestan)
stan_file <- normalizePath("benchmarks/dpm_heckman/dpm_heckman.stan")
cat("Compiling", stan_file, "...\n")
model_so <- compile_model(stan_file)
cat("Compiled to:", model_so, "\n")

# Minimal dummy data to test instantiation
set.seed(1)
N <- 20L; P_out <- 3L; P_sel <- 4L; K <- 5L
data_list <- list(
  N = N, P_out = P_out, P_sel = P_sel, K = K,
  X = matrix(rnorm(N * P_out), N, P_out),
  W = matrix(rnorm(N * P_sel), N, P_sel),
  Z = as.integer(rbinom(N, 1, 0.5)),
  Y = rnorm(N),
  alpha = 1.0,
  sigma_beta = 5.0,
  sigma_gamma = 5.0
)
data_json <- jsonlite::toJSON(data_list, auto_unbox = TRUE, pretty = FALSE)
writeLines(data_json, "benchmarks/dpm_heckman/.dummy_data.json")

model <- StanModel$new(lib = model_so,
                       data = "benchmarks/dpm_heckman/.dummy_data.json",
                       seed = 1L)
cat("Model instantiated. Unconstrained dim:", model$param_unc_num(), "\n")
theta_unc <- rnorm(model$param_unc_num(), sd = 0.1)
lp_grad <- model$log_density_gradient(theta_unc)
cat("log density at random init:", lp_grad$val, "\n")
cat("gradient norm:", sqrt(sum(lp_grad$gradient^2)), "\n")
lp_hess <- model$log_density_hessian(theta_unc)
cat("hessian symmetric:",
    isSymmetric(lp_hess$hessian, tol = 1e-8), "\n")
cat("hessian dim:", dim(lp_hess$hessian), "\n")
cat("OK\n")
