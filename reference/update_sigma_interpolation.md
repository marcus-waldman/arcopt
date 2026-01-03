# Update Regularization Parameter (Interpolation-Enhanced CGT)

Adaptively adjusts sigma using interpolation to estimate local Lipschitz
constant from model-function discrepancy. May provide faster adaptation
on large-scale problems with varying curvature (Gould, Porcelli, Toint,
2012).

## Usage

``` r
update_sigma_interpolation(
  sigma_current,
  rho,
  s,
  f_trial,
  m_trial,
  eta1 = 0.1,
  eta2 = 0.9,
  gamma1 = 0.5,
  gamma2 = 2,
  sigma_min = 1e-16,
  sigma_max = 1e+16
)
```

## Arguments

- sigma_current:

  Current regularization parameter (positive scalar)

- rho:

  Ratio of actual to predicted reduction (scalar)

- s:

  Current step vector

- f_trial:

  Function value at trial point f(x + s)

- m_trial:

  Model value at step m(s)

- eta1:

  Acceptance threshold (default: 0.1)

- eta2:

  Very successful threshold (default: 0.9)

- gamma1:

  Sigma decrease factor for very successful steps (default: 0.5)

- gamma2:

  Sigma increase factor for unsuccessful steps (default: 2.0)

- sigma_min:

  Minimum allowed sigma (default: 1e-16)

- sigma_max:

  Maximum allowed sigma (default: 1e16)

## Value

Updated sigma value (scalar)

## Details

Estimates local Lipschitz constant from discrepancy between function and
model: L_hat = 2 \* \|f(x+s) - m(s)\| / \|\|s\|\|^3

The update strategy incorporates this estimate:

- **Very successful** (rho \>= eta2): sigma = max(L_hat/2, sigma/gamma1,
  sigma_min)

- **Moderately successful** (eta1 \<= rho \< eta2): Keep sigma unchanged

- **Unsuccessful** (rho \< eta1): sigma = min(max(L_hat/2,
  gamma2\*sigma), sigma_max)

This may provide better adaptation than Algorithm 2a (classical CGT)
when curvature varies significantly across the domain.
