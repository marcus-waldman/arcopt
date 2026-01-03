# Update Regularization Parameter (Classical CGT)

Adaptively adjusts the regularization parameter sigma based on how well
the cubic model predicted the actual function decrease. Follows the
classical Cartis-Gould-Toint (2011) update strategy.

## Usage

``` r
update_sigma_cgt(
  sigma_current,
  rho,
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

- eta1:

  Acceptance threshold (default: 0.1). Steps with rho \>= eta1 are
  accepted.

- eta2:

  Very successful threshold (default: 0.9). Steps with rho \>= eta2
  trigger sigma decrease.

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

The update strategy is:

- **Very successful** (rho \>= eta2): Decrease sigma by gamma1 factor
  (trust model more)

- **Moderately successful** (eta1 \<= rho \< eta2): Keep sigma unchanged

- **Unsuccessful** (rho \< eta1): Increase sigma by gamma2 factor (trust
  model less)

The ratio rho measures model quality: rho = (f(x) - f(x+s)) /
pred_reduction. Values near 1 indicate excellent model accuracy; values
near 0 or negative indicate poor model fit.
