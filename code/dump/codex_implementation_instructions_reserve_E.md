# Codex implementation brief: reserve E in the Case5 workflow

## Task

Extend the existing Case5.ipynb analysis to include firm battery reserve E in:

1. Optimal bid selection.
2. Allocation of a specified bid.

Use reserve_E_bid_and_allocation.md as the mathematical specification. Implement using the same code organization, marginal models, dependence conventions, numerical methods, plotting style, and result formats as Case5.ipynb wherever applicable.

**Availability note:** Case5.ipynb was not provided when this brief was written. No notebook cells, helper names, solver choices, or existing numerical results have been inspected. The steps below are requirements to map onto the actual notebook, not claims about its contents.

## 1. Inspect the project before implementing

- Read applicable AGENTS.md instructions and the supplied discussion file.
- Locate and read Case5.ipynb and any modules it imports from the project. If it is unavailable, request that notebook; do not create a substitute and claim it follows the original structure.
- Identify the existing parameter/configuration cells, units, marginal distributions, dependence models, sampling or integration routines, quantile conventions, optimization routines, allocation objective, and plotting helpers.
- Establish whether the baseline is known dependence, a comparison across copulas, or a certified robust calculation under a specified ambiguity set.
- Identify whether baseline delivery uses fixed individual targets or permits asset-to-asset surplus redispatch. The main reserve model in this brief uses fixed targets with battery recourse. If the notebook uses pooled delivery, implement the pooled model as its direct extension and label the fixed-target model as an additional comparison.
- Record this mapping in a short opening Markdown cell. Do not invent economic objectives or dependence assumptions.

Create Case5_with_reserve.ipynb by default, preserving Case5.ipynb. Follow project instructions if they explicitly require an in-place extension. Avoid unnecessary refactoring or new dependencies.

## 2. Modeling contract

Use nonnegative energies Xi, targets ci, reserve E, and one common delivery interval. Default alpha=0.9. E is firm usable energy with sufficient discharge capability for the interval.

For fixed assigned targets, define

```text
shortfall_i = max(c_i - X_i, 0)
total_shortfall = sum_i shortfall_i
success = total_shortfall <= E
reserve_used = min(total_shortfall, E)
uncovered_shortfall = max(total_shortfall - E, 0)
```

The reserve is shared across assets and allocated after shortfalls are observed. Do not give every asset an independent reserve E. Use a documented numerical tolerance consistently and show that it does not materially alter the selected result.

The targets include energy that may be supplied by the battery. Therefore b=sum(c); do not add E to that bid again.

If there are additional physical or target constraints, preserve and document them. Clarify whether backed targets can exceed an asset's standalone capacity. The additive benchmark bound in the discussion assumes unrestricted nonnegative backed targets; check that assumption before asserting it in the notebook.

## 3. Reusable calculations

Reuse existing helper patterns and naming. The following conceptual interface describes needed behavior; adapt its signatures to the notebook rather than duplicating existing infrastructure:

```python
assigned_shortfall(samples, c)           # samples shape (n, d), output (n,)
reserve_success_probability(c, E, ...)  # uses the existing probability engine
required_reserve(c, alpha, ...)         # lower alpha-quantile of shortfall
optimize_bid_with_reserve(E, alpha, ...)
allocate_bid_with_reserve(b, E, alpha, ...)
```

Where samples are used, `assigned_shortfall` is conceptually

```python
np.maximum(c[None, :] - samples, 0.0).sum(axis=1)
```

Preserve existing integration or analytic calculations when they are appropriate. Do not replace an existing robust probability calculation with a single simulated joint distribution.

For required reserve, implement

```text
inf {e >= 0 : P(total_shortfall <= e) >= alpha}
```

With empirical samples, choose a quantile/order-statistic convention that actually satisfies the empirical inequality; interpolated quantiles can give the wrong coverage with atoms or small samples. Treat E=0 using the same non-strict event as the baseline.

## 4. (a) Bid optimization

Implement

```text
maximize    sum(c)
over        c >= 0
subject to  P(sum_i max(c_i-X_i, 0) <= E) >= alpha
            existing applicable constraints
```

Use the notebook's optimizer or contour-search method where appropriate. The empirical chance constraint can be discontinuous; do not assume a smooth convex problem or report a local optimizer as globally optimal without justification. Include feasible starting points, search bounds, convergence status, and final probability checks. For an existing two-asset contour-grid approach, extend it directly and report resolution limitations.

Use a reserve sweep containing E=0 and several positive values scaled to the existing asset energies. Choose scales from the notebook rather than arbitrary units. Reuse common samples across reserve levels where sampling is appropriate.

Compare:

- Baseline AND result at E=0.
- Conservative additive benchmark b_AND(0)+E, under its stated assumptions.
- Optimized assigned-deficit bid b*(E).
- Optional pooled benchmark, distinctly labeled.

Where the assumptions hold, b*(E) should be nondecreasing and no smaller than the additive feasible benchmark. Numerical violations require investigation, not smoothing or relabeling.

## 5. (b) Allocation for a specified bid

Implement feasibility for

```text
c >= 0
sum(c) = b
P(sum_i max(c_i-X_i, 0) <= E) >= alpha
existing applicable constraints
```

Use the baseline allocation objective if it has one. If the notebook only finds feasible allocations, keep that terminology and report the feasible set or a documented feasible representative. Do not call the first feasible point optimal.

As a reserve-oriented allocation analysis, implement where consistent with scope:

```text
E_min(b) = min over c >= 0, sum(c)=b of required_reserve(c, alpha)
```

Then compare E_min(b) with the available E. Label this as minimizing required reserve, not as an unstated economic allocation objective.

Choose b values that demonstrate baseline feasibility, feasibility enabled by reserve, and infeasibility for the selected reserve. Determine these from computed results. For d=2, parameterize c2=b-c1 and reuse existing one-dimensional allocation-search/plotting patterns.

## 6. Dependence and robustness

Retain the Case5 dependence setting. If the project studies unknown dependence, use its stated ambiguity set and implement the constraint

```text
inf over admissible joint laws P of P(total_shortfall <= E) >= alpha
```

A finite set of copulas gives a comparison or a finite-family robust result only. Clearly identify that scope. If a certified ambiguity-set solver is unavailable, provide the implementable comparison and explicitly leave the full robust extension unresolved; do not label it solved.

If Case5 uses a rearrangement algorithm for aggregate quantiles, it cannot be reused unchanged on sum(Xi) for the assigned-deficit model. The relevant random variables are Zi=(ci-Xi)+, and their distributions change with c. Derive the appropriate transformed marginals and solver orientation before adapting the method. Preserve distinctions between heuristic numerical bounds and certified guarantees.

## 7. Figures and result tables

Follow existing notebook styling. Produce:

1. Bid against E, including baseline/additive/optimized curves and the pooled benchmark if implemented.
2. For d=2, contours of P(total_shortfall<=E)=alpha in allocation space for selected E, with maximizing targets and sum(c)=b lines. If exact equality is absent on a discrete grid, show the feasible boundary and say so.
3. A fixed-bid allocation figure along sum(c)=b, showing reliability or required reserve and the feasible allocations.
4. A realized-capacity-space example showing the event total_shortfall<=E for one fixed c. Label it separately from allocation-space contours.
5. A concise table with dependence setting, E, bid, c, optimization probability, validation probability where applicable, expected reserve used, uncovered-shortfall probability, and solver status.

Do not draw a Kendall contour and label its probability level as battery energy. If a Kendall comparison is included, label it an optional geometric comparison rather than the reserve model.

## 8. Meaningful verification

Use lightweight checks tied to modeling risks:

- E=0 recovers the original AND event and the baseline numerical solution within the established resolution/tolerance.
- For fixed c, success probability does not decrease with E. For fixed E, increasing any component of c cannot improve success probability.
- One-asset unconstrained case: the optimal reliable target increases by exactly E relative to its E=0 value.
- Deterministic example: c=(6,4), X=(5,3) gives total shortfall 2; E=1 fails and E=2 succeeds. This checks shared-budget accounting.
- Distinguish pooling: c=(6,4), X=(8,2), E=0 fails the assigned-target condition but satisfies aggregate delivery of 10 under free pooling.
- If using the illustrative two-state example from the discussion, verify that targets (2,2) with E=2 succeed for both (10,0) and (0,10).
- Confirm every returned allocation sums to its bid and satisfies the declared constraints.
- Where an additive benchmark is applicable, use its constructed allocation as a feasible candidate or numerical lower bound.

For sample-based optimization, evaluate selected decisions on an independent validation sample. Report a binomial uncertainty interval; an estimated success rate near 0.9 does not certify a true 90% guarantee. If the study requires a statistical certification rule, specify it before using that label. A nominal 90% solution need not have a validation lower confidence bound above 0.9.

Do not add a large testing framework solely for this extension. Run the notebook in order and inspect the key outputs. If source data, dependencies, or the original notebook are missing, state the exact blocked step and do not fabricate execution results.

## 9. Deliverables

- Case5_with_reserve.ipynb, following Case5.ipynb's organization.
- Any minimal helper changes required by the existing project architecture.
- Figures and result tables in the project's established output locations.
- A brief findings note documenting the operational model, dependence scope, numerical method, bid gains, allocation behavior, validation, and remaining limitations.

The implementation is complete when both bid selection and fixed-bid allocation use E consistently, the zero-reserve baseline is recovered, and the reported reliability is supported at the scope claimed.
