# Apply Box Constraints to Step

Truncates a proposed step to ensure the new iterate stays within box
constraints, then projects for numerical safety.

## Usage

``` r
apply_box_constraints(x_current, s, lower, upper)
```

## Arguments

- x_current:

  Current iterate

- s:

  Proposed unconstrained step

- lower:

  Lower bounds (use -Inf for unbounded)

- upper:

  Upper bounds (use Inf for unbounded)

## Value

List with components:

- s_bounded:

  Truncated step that respects bounds

- x_new:

  New iterate projected to feasible region

## Details

This implements Algorithm 7a from the design specification. The
approach:

1.  Compute maximum feasible step length alpha_max

2.  Truncate step: s_bounded = alpha_max \* s

3.  Project x_k + s_bounded to box for numerical safety
