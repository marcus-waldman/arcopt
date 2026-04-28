# Flat-Ridge Detector for Cubic-to-Trust-Region Fallback

Detects when cubic regularization has entered a flat-ridge stagnation
regime: gradient stalls in a near-singular positive-definite region
where the cubic penalty has no bite. Triggers the trust-region fallback
in the main arcopt iteration when fired.

## Details

The detector maintains a sliding window of runtime diagnostics and fires
when **all four** of the following signals hold throughout the window:

1.  **Sigma pinned at floor**: `sigma_k <= 10 * sigma_min` for every
    iteration in the window.

2.  **Near-perfect model**: `|rho_k - 1| < rho_tol` for every iteration
    in the window (model reduction matches actual reduction).

3.  **Gradient stagnant**: The ratio of the most recent `||g||_inf` to
    the oldest in the window exceeds `grad_decrease_max` (less than 10
    percent decrease over the window by default).

4.  **Hessian ill-conditioned**: The smallest Hessian eigenvalue at the
    most recent iteration satisfies `lambda_min < tol_ridge`. This
    includes both the classical "flat-ridge" regime (small positive
    `lambda_min`) and the "stuck-at-indefinite-saddle" regime (any
    negative `lambda_min`), because cubic regularization loses its grip
    in both.

An absolute floor `g_inf_floor` is also enforced: the trigger will not
fire if `||g||_inf` is already below the typical convergence tolerance,
preventing spurious triggers at true local minima.

The detector is an empirical proxy for local-error-bound (EB) violation
(Yue, Zhou, So 2018). Under EB, cubic regularization attains Q-quadratic
convergence even at degenerate minima; the detector fires precisely in
the regime where EB is practically vacuous and no theoretical
convergence guarantee applies.
