# Reserve E study: computed findings

Notebook: `Case5_with_reserve.ipynb`. Mathematical specification: `../reserve_E_bid_and_allocation.md`. The original Case5 notebook and marginal generator were not modified.

## Setup and methods

The executed example uses two ordinary assets, taking the first two latent parameter pairs from the current Case5 50-asset generation (seed 0), alpha=0.9 and Gaussian-copula rho=0.5. CapacityMarginal is imported from library/generate_marginals.py. Targets are fixed, surplus cannot be transferred, and a single firm reserve E covers the sum of realized positive deficits. The bid already equals sum(c); E is never added twice.

The seven-section organization follows Case5. At E=0 the robust baseline reuses its marginal-quantile MILP and active-set polishing. For E>0 the unknown-dependence calculation uses a marginal transport LP on lower/upper capacity grids plus valid deterministic-reserve-split certificates. The optimized robust curve is a feasible lower estimate, not the exact continuous robust optimum.

The Gaussian comparison uses conditional quadrature in the default two-asset example, including latent intervals associated with zero atoms. Higher-dimensional configurations explicitly fall back to empirical Gaussian probabilities. A finite direction search and bisection maximize sum(c) subject to the selected survival constraint; this is best-found optimization, not a global-optimality certificate. The separate fixed-bid study minimizes sufficient/required reserve on the same direction grid.

## Computed comparison

| Dependence calculation | E=0 bid | Largest E | Bid at largest E | Gain over own additive baseline |
|---|---:|---:|---:|---:|
| robust | 11.065823 | 6.386664 | 17.452487 | 0.000000 |
| gaussian | 11.371020 | 6.386664 | 19.923718 | 2.166033 |

The robust search found only the additive gain for this example. This is a result for the configured conservative model and numerical search, not a proof that arbitrary portfolios can never obtain additional reserve benefits. The Gaussian comparison exhibits a gain beyond its own additive baseline and approaches the separate pooled reference.

## Checks and limitations

All notebook cells executed without errors. Checks cover shared-budget accounting, the E=0 AND event, strict marginal failure at zero, monotonicity, additive feasibility, exact empirical quantile coverage where used, fixed-bid equality, the one-asset reserve identity, and returned allocation feasibility. Additional checks confirmed mixed zero-atom lower/upper bounds and Gaussian quadrature against the independent AND identity and Monte Carlo.

The first empirical Gaussian optimization showed out-of-sample undercoverage. It was replaced by deterministic conditional quadrature for the executed example. Final optimization probabilities are at least 0.9; independent 50,000-draw Gaussian validation intervals include 0.9. Those intervals are uncertainty summaries, not separate statistical certifications. Validation reserve-use statistics refer to the Gaussian law even for robustly chosen allocations.

The robust probability bound at a selected middle-E allocation was 0.900000 to 0.914286; a different, finer but nonnested grid gave 0.900000 to 0.915789. Thus the feasible lower bound is stable in this check, while the upper probability bound remains approximate and grid-dependent. This is not a global bid error bound.

Transport state count grows exponentially; d=2 is deliberate for an executable full contour study, and the notebook guards excessive state counts. It does not silently claim to solve the full 50-asset robust reserve optimization.

## Saved artifacts

- ../figures/Case5_with_reserve_results.csv: all reserve levels, allocations, probability bounds, Gaussian validation intervals and reserve-use summaries.
- ../figures/Case5_with_reserve_bid_vs_E.png: bid, additive gain and validation comparison.
- ../figures/Case5_with_reserve_contours.png: allocation-space survival contours.
- ../figures/Case5_with_reserve_fixed_bid.png: fixed-bid survival comparison.
- All five figures and output tables remain embedded in the executed notebook.
