#' Update Regularization Parameter (Classical CGT)
#'
#' Adaptively adjusts the regularization parameter sigma based on how well the
#' cubic model predicted the actual function decrease. Follows the classical
#' Cartis-Gould-Toint (2011) update strategy.
#'
#' @param sigma_current Current regularization parameter (positive scalar)
#' @param rho Ratio of actual to predicted reduction (scalar)
#' @param eta1 Acceptance threshold (default: 0.1). Steps with rho >= eta1 are accepted.
#' @param eta2 Very successful threshold (default: 0.9). Steps with rho >= eta2 trigger sigma decrease.
#' @param gamma1 Sigma decrease factor for very successful steps (default: 0.5)
#' @param gamma2 Sigma increase factor for unsuccessful steps (default: 2.0)
#' @param sigma_min Minimum allowed sigma (default: 1e-16)
#' @param sigma_max Maximum allowed sigma (default: 1e16)
#'
#' @return Updated sigma value (scalar)
#'
#' @details
#' The update strategy is:
#' \itemize{
#'   \item \strong{Very successful} (rho >= eta2): Decrease sigma by gamma1 factor (trust model more)
#'   \item \strong{Moderately successful} (eta1 <= rho < eta2): Keep sigma unchanged
#'   \item \strong{Unsuccessful} (rho < eta1): Increase sigma by gamma2 factor (trust model less)
#' }
#'
#' The ratio rho measures model quality: rho = (f(x) - f(x+s)) / pred_reduction.
#' Values near 1 indicate excellent model accuracy; values near 0 or negative
#' indicate poor model fit.
#'
#' @keywords internal
update_sigma_cgt <- function(sigma_current, rho,
                             eta1 = 0.1,
                             eta2 = 0.9,
                             gamma1 = 0.5,
                             gamma2 = 2.0,
                             sigma_min = 1e-16,
                             sigma_max = 1e16) {
  if (rho >= eta2) {
    # Very successful: decrease sigma (trust model more)
    sigma_new <- max(sigma_min, gamma1 * sigma_current)
  } else if (rho >= eta1) {
    # Moderately successful: keep sigma
    sigma_new <- sigma_current
  } else {
    # Unsuccessful: increase sigma (trust model less)
    sigma_new <- min(sigma_max, gamma2 * sigma_current)
  }

  return(sigma_new)
}


#' Update Regularization Parameter (Interpolation-Enhanced CGT)
#'
#' Adaptively adjusts sigma using interpolation to estimate local Lipschitz constant
#' from model-function discrepancy. May provide faster adaptation on large-scale
#' problems with varying curvature (Gould, Porcelli, Toint, 2012).
#'
#' @param sigma_current Current regularization parameter (positive scalar)
#' @param rho Ratio of actual to predicted reduction (scalar)
#' @param s Current step vector
#' @param f_trial Function value at trial point f(x + s)
#' @param m_trial Model value at step m(s)
#' @param eta1 Acceptance threshold (default: 0.1)
#' @param eta2 Very successful threshold (default: 0.9)
#' @param gamma1 Sigma decrease factor for very successful steps (default: 0.5)
#' @param gamma2 Sigma increase factor for unsuccessful steps (default: 2.0)
#' @param sigma_min Minimum allowed sigma (default: 1e-16)
#' @param sigma_max Maximum allowed sigma (default: 1e16)
#'
#' @return Updated sigma value (scalar)
#'
#' @details
#' Estimates local Lipschitz constant from discrepancy between function and model:
#' L_hat = 2 * |f(x+s) - m(s)| / ||s||^3
#'
#' The update strategy incorporates this estimate:
#' \itemize{
#'   \item \strong{Very successful} (rho >= eta2): sigma = max(L_hat/2, sigma/gamma1, sigma_min)
#'   \item \strong{Moderately successful} (eta1 <= rho < eta2): Keep sigma unchanged
#'   \item \strong{Unsuccessful} (rho < eta1): sigma = min(max(L_hat/2, gamma2*sigma), sigma_max)
#' }
#'
#' This may provide better adaptation than Algorithm 2a (classical CGT) when
#' curvature varies significantly across the domain.
#'
#' @keywords internal
update_sigma_interpolation <- function(sigma_current, rho, s, f_trial, m_trial,
                                       eta1 = 0.1,
                                       eta2 = 0.9,
                                       gamma1 = 0.5,
                                       gamma2 = 2.0,
                                       sigma_min = 1e-16,
                                       sigma_max = 1e16) {
  # STEP 1: Estimate local Lipschitz constant
  step_norm <- sqrt(sum(s^2))
  if (step_norm < 1e-14) {
    # Avoid division by zero for very small steps
    l_hat <- sigma_current
  } else {
    l_hat <- 2 * abs(f_trial - m_trial) / (step_norm^3)
  }

  # STEP 2: Update with interpolation
  if (rho >= eta2) {
    # Very successful: use interpolation with CGT decrease
    sigma_new <- max(l_hat / 2, gamma1 * sigma_current, sigma_min)
  } else if (rho >= eta1) {
    # Moderately successful: keep sigma
    sigma_new <- sigma_current
  } else {
    # Unsuccessful: use interpolation with CGT increase
    sigma_new <- max(l_hat / 2, gamma2 * sigma_current)
    sigma_new <- min(sigma_new, sigma_max)
  }

  return(sigma_new)
}
