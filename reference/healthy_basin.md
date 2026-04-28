# Healthy-Basin Detector for Cubic-to-QN-Polish Transition

Detects when the iterate has entered the quadratic attraction basin of a
strict local minimum, at which point arcopt can switch from cubic
regularization to line-search BFGS (qn_polish mode) and skip further
Hessian evaluations for the remainder of the run.

## Details

The detector maintains a sliding window of runtime diagnostics and fires
when **all five** of the following signals hold throughout the window:

1.  **Newton step accepted**: `step_type == "newton"` at every iteration
    in the window. Indicates the cubic subproblem has not been invoked
    (Cholesky succeeded and the Newton step passed the ratio test).

2.  **High ratio**: `rho_k >= rho_polish` (default `0.9`) throughout the
    window. The quadratic model is accurately predicting actual
    reduction.

3.  **Hessian well-conditioned PD**:
    `lambda_min(H_k) >= lambda_min_polish` (default `1e-3`) throughout
    the window. Not just "Cholesky succeeded" but strictly bounded away
    from zero.

4.  **Superlinear gradient decay**:
    `||g_k||_inf / ||g_{k-1}||_inf <= g_decay_polish` (default `0.5`)
    for every consecutive pair in the window. This is the strongest
    signal of the quadratic attraction basin.

5.  **Above convergence floor**: `||g_0||_inf > g_inf_floor_polish` at
    window start (default `1e-8`). Prevents firing at convergence where
    Newton-accept and gradient-decay are trivially satisfied.

Contrast with the flat-ridge detector
([`flat_ridge`](https://marcus-waldman.github.io/arcopt/reference/flat_ridge.md)):
that detector fires on *stagnation* (gradient not decreasing); this one
fires on *healthy convergence* (gradient decreasing super- linearly).
Both are expected to be mutually exclusive in practice.
