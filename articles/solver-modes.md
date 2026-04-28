# Solver modes: cubic, trust-region fallback, and quasi-Newton polish

``` r
library(arcopt)
```

`arcopt` solves its subproblem in one of three modes, chosen adaptively
at runtime:

| Mode        | Subproblem                                        | Targeted regime                                                |
|-------------|---------------------------------------------------|----------------------------------------------------------------|
| `cubic`     | Eigendecomposition cubic regularization (default) | Indefinite Hessians, saddle points                             |
| `tr`        | Trust-region with adaptive radius                 | Flat ridges where the cubic penalty has decayed to its floor   |
| `qn_polish` | Wolfe-safeguarded BFGS line search (opt-in)       | Strongly convex healthy basins where `hess()` calls are wasted |

The default configuration runs in `cubic` mode and switches *only* if
its detector signals trigger. The trust-region fallback is enabled by
default but dormant on well-behaved problems. Polish mode is opt-in via
`control = list(qn_polish_enabled = TRUE)`. Detector thresholds are
documented in
[`?arcopt_advanced_controls`](https://marcus-waldman.github.io/arcopt/reference/arcopt_advanced_controls.md).

The `result$diagnostics` sublist reports which mode the run finished in
and how many transitions occurred:

- `solver_mode_final` – `"cubic"`, `"tr"`, or `"qn_polish"`
- `ridge_switches` – count of cubic-to-TR fallback firings
- `qn_polish_switches`, `qn_polish_reverts` – bidirectional polish
  transition counts
- `hess_evals_at_polish_switch` – Hessian-evaluation count at the first
  polish transition

## Mode 1: cubic regularization (default)

The cubic subproblem eigendecomposes the true Hessian and writes its
step in the eigenbasis. A single negative eigenvalue is enough to
produce a step in the negative-curvature direction, so saddle escape
with the true Hessian is automatic. As a small example, consider a
two-component Gaussian mixture likelihood started from a symmetric
near-saddle point:

``` r
set.seed(42)
y <- c(rnorm(100, -2, 1), rnorm(100, 2, 1))

mixture_nll <- function(theta) {
  -sum(log(0.5 * dnorm(y, theta[1], 1) + 0.5 * dnorm(y, theta[2], 1)))
}
mixture_gr <- function(theta) {
  p1 <- 0.5 * dnorm(y, theta[1], 1)
  p2 <- 0.5 * dnorm(y, theta[2], 1)
  w1 <- p1 / (p1 + p2)
  w2 <- p2 / (p1 + p2)
  -c(sum(w1 * (y - theta[1])), sum(w2 * (y - theta[2])))
}
mixture_hess <- function(theta, h = 1e-5) {
  d <- 2L
  h_mat <- matrix(0, d, d)
  for (i in seq_len(d)) {
    e_i <- numeric(d)
    e_i[i] <- h
    h_mat[, i] <- (mixture_gr(theta + e_i) - mixture_gr(theta - e_i)) /
      (2 * h)
  }
  0.5 * (h_mat + t(h_mat))
}

res <- arcopt(c(0.01, 0.01), mixture_nll, mixture_gr, mixture_hess,
              control = list(trace = 0))
res$par               # near (-2, 2) modulo label permutation
#> [1]  1.866497 -1.980722
res$diagnostics$solver_mode_final
#> [1] "cubic"
min(eigen(res$hessian, only.values = TRUE)$values)  # > 0: local min
#> [1] 73.04552
```

Started from the symmetric point near the origin, BFGS-class methods
land at the saddle and report convergence; arcopt’s cubic mode steps in
the negative-curvature direction and recovers the class means.

## Mode 2: trust-region fallback (default-on, dormant)

Cubic regularization can stall in *flat-ridge* regimes – iterations
where the regularization parameter `sigma_k` has been driven to its
floor, the step-acceptance ratio `rho_k` is near 1, and yet the gradient
is bounded away from zero because the Hessian is positive definite but
near-singular in some weakly identified subspace. arcopt detects this
with a sliding-window check on four signals and, on a hold of ten
consecutive iterations, swaps the cubic subproblem for a trust-region
subproblem. The transition is one-way and default-on.

For well-conditioned problems the detector simply does not fire, and the
run terminates in cubic mode:

``` r
rosen_fn <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
rosen_gr <- function(x) {
  c(-2 * (1 - x[1]) - 400 * x[1] * (x[2] - x[1]^2),
    200 * (x[2] - x[1]^2))
}
rosen_hess <- function(x) {
  matrix(c(1200 * x[1]^2 - 400 * x[2] + 2, -400 * x[1],
           -400 * x[1],                     200), 2, 2)
}

res <- arcopt(c(-1.2, 1), rosen_fn, rosen_gr, rosen_hess,
              control = list(trace = 0))
res$diagnostics$solver_mode_final
#> [1] "cubic"
res$diagnostics$ridge_switches
#> [1] 0
```

The fallback’s value shows up on problems with weakly identified
parameter subspaces – e.g., Dirichlet-process-style mixtures with unused
stick-breaking segments, where the cubic-only configuration exhausts its
iteration budget at a strongly indefinite saddle but the trust-region
fallback fires once and reaches a true local minimum. That regime is too
data- and dependency-heavy for an inline example; see the JSS manuscript
referenced from the package README for a worked case.

If you want to disable the fallback explicitly:

``` r
arcopt(x0, fn, gr, hess,
       control = list(tr_fallback_enabled = FALSE))
```

## Mode 3: quasi-Newton polish (opt-in)

In a strongly convex healthy basin where the cubic regularization
parameter has dropped to its floor and the Hessian is essentially
constant across iterates, every `hess()` call recomputes second
derivatives that contribute negligibly to step quality. Polish mode
recognizes this regime via a five-signal detector and replaces the cubic
subproblem with a Wolfe line search along a BFGS-approximated Newton
direction, suspending further `hess()` calls until either convergence or
a revert.

Polish is *opt-in* (`qn_polish_enabled = FALSE` by default). Enable it
for long-running smooth problems with expensive Hessians (analytic AD
via Stan, finite differences, or any pipeline where each second
derivative takes appreciable wall-clock time):

``` r
fn   <- function(x) sum(0.5 * x^2 + x^4 / 12)
gr   <- function(x) x + x^3 / 3
hess <- function(x) diag(1 + x^2)
x0   <- rep(10, 5)

res_default <- arcopt(x0, fn, gr, hess,
                      control = list(maxit = 50, trace = 0))
res_polish  <- arcopt(x0, fn, gr, hess,
                      control = list(maxit = 50, trace = 0,
                                     qn_polish_enabled = TRUE))

c(default_hess = res_default$evaluations$hess,
  polish_hess  = res_polish$evaluations$hess)
#> default_hess  polish_hess 
#>            9            6

res_polish$diagnostics$solver_mode_final     # "qn_polish"
#> [1] "qn_polish"
res_polish$diagnostics$qn_polish_switches    # 1
#> [1] 1
```

On this five-dimensional smooth-convex problem the polish detector fires
once after the initial cubic-Newton approach steps, and the remaining
iterations consume zero `hess()` calls – typically a saving of 30-50% of
the Hessian budget at no loss of accuracy.

## Choosing modes

The defaults are designed so most users do not need to touch the mode
switches. The cubic mode handles indefiniteness and saddle escape; the
trust-region fallback handles flat-ridge stagnation; and polish is
reserved for problems where `hess()` is genuinely expensive and the
basin is large enough for the detector window to fill. A practical
heuristic:

- If the iterate enters a saddle and stalls – you have the right
  default. Make sure `hess()` returns the *true* Hessian, not a
  positive-definite proxy.
- If `arcopt` exhausts its iteration budget on a problem that should be
  locally well-conditioned – inspect
  `result$diagnostics$ridge_switches`. If it is 0 and the gradient is
  bounded away from zero, you may be in a flat-ridge regime; the
  TR-fallback should already be on by default.
- If a Hessian-supplied run is fast in clock time but the per-call
  Hessian is expensive (Stan AD, FD), enable `qn_polish_enabled = TRUE`
  and check `result$diagnostics$qn_polish_switches`.

For any of these cases the diagnostic fields tell you which mode
actually carried the run, and the per-mode controls in
[`?arcopt_advanced_controls`](https://marcus-waldman.github.io/arcopt/reference/arcopt_advanced_controls.md)
provide the tuning surface.
