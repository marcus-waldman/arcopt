// Truncated Dirichlet process mixture Heckman selection model
// ===========================================================
//
// N observations; for each i:
//   Z_i in {0, 1}   — selection indicator (1 = selected/observed)
//   Y_i in R        — outcome (meaningful when Z_i = 1)
//   X_i in R^P_out  — covariates for outcome equation
//   W_i in R^P_sel  — covariates for selection equation
//
// Mixture of K truncated stick-breaking components. Within component k:
//   outcome:   Y_i | comp = k, selected = mu^o_{i,k} + eps,   eps ~ N(0, sigma_k^2)
//   selection: Z_i = 1 iff mu^s_{i,k} + u > 0,    u ~ N(0, 1)
//   correlation: corr(eps/sigma_k, u) = rho_k in (-1, 1)
// with mu^o_{i,k} = X_i . beta_k and mu^s_{i,k} = W_i . gamma_k.
//
// For the mixture likelihood we use the standard selection-model form
// (e.g., Van de Ven & Van Pragg 1981), with correlation handled through
// the conditional selection probit.
//
// Weak priors on component parameters regularize empty components so
// the MAP objective has no flat ridges; label switching between active
// components remains as the arcopt-shaped symmetric saddle.

data {
  int<lower=1> N;
  int<lower=1> P_out;
  int<lower=1> P_sel;
  int<lower=1> K;
  matrix[N, P_out] X;
  matrix[N, P_sel] W;
  array[N] int<lower=0, upper=1> Z;
  vector[N] Y;                       // pass 0 for Z_i = 0 rows; unused
  real<lower=0> alpha;               // stick-breaking concentration
  real<lower=0> sigma_beta;          // weak prior SD on beta
  real<lower=0> sigma_gamma;         // weak prior SD on gamma
}

parameters {
  matrix[K, P_out] beta;
  matrix[K, P_sel] gamma;
  vector<lower=0>[K] sigma;
  vector<lower=-1, upper=1>[K] rho;
  vector<lower=0, upper=1>[K - 1] V;
}

transformed parameters {
  vector[K] log_pi;
  {
    log_pi[1] = log(V[1]);
    real log_rem = log1m(V[1]);
    for (k in 2 : (K - 1)) {
      log_pi[k] = log_rem + log(V[k]);
      log_rem += log1m(V[k]);
    }
    log_pi[K] = log_rem;
  }
}

model {
  // Weak priors (regularize empty components, keep the MAP well-posed).
  to_vector(beta) ~ normal(0, sigma_beta);
  to_vector(gamma) ~ normal(0, sigma_gamma);
  sigma ~ lognormal(0, 1);
  // rho implicitly uniform on (-1, 1) via its declared bounds (the
  // Jacobian from the unconstrained parametrization handles this).
  V ~ beta(1, alpha);

  // Mixture likelihood.
  for (i in 1 : N) {
    vector[K] lp;
    for (k in 1 : K) {
      real mu_sel = dot_product(W[i], gamma[k]);
      if (Z[i] == 0) {
        lp[k] = log_pi[k] + std_normal_lcdf(-mu_sel);
      } else {
        real mu_out = dot_product(X[i], beta[k]);
        real eps_std = (Y[i] - mu_out) / sigma[k];
        real z_cond = (mu_sel + rho[k] * eps_std)
                      / sqrt(1 - square(rho[k]));
        lp[k] = log_pi[k]
                + normal_lpdf(Y[i] | mu_out, sigma[k])
                + std_normal_lcdf(z_cond);
      }
    }
    target += log_sum_exp(lp);
  }
}
