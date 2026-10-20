# Confidence-adjusted aggregation of two Gaussian flexibility envelopes under Energinet P90: effect of the coefficient vector a, the haircut construction, and the unknown-dependence worst case

The coefficient vector a is the wrong place to put asset confidence: scaling X_i by a_i < 1 shrinks its mean and standard deviation together. The construction you describe (fix x1 at its P90, push x2 down) does not lower α. It raises the aggregate nominal reliability above 0.9, which is the intended buffer. It also has to start from the unknown-dependence worst-case bid, not from the comonotonic bid VaR^M. In the worked example (X1~N(10,2²), X2~N(15,4²), ε=0.1) the VaR^M bid of 17.31 is only 81.4% reliable under the worst-case copula. The worst-case bid is 15.335. Confidence is best expressed as an ambiguity band on each asset's lower-tail quantile function, calibrated from its delivery history. For d=2 this works out exactly: a Kolmogorov–Smirnov (KS) band of width η_2 on X2 is the same as running the Makarov worst case at ε − η_2.

## TL;DR
- **The direction of a, not its size, carries the information, and derating is not "low confidence".** VaR_α(a·X) is positively homogeneous. Its gradient equals the Euler allocation and the first-order reliability method (FORM) design point. a_i<1 scales both μ_i and σ_i, so the coefficient of variation stays the same. It can match a mean shift or variance inflation at one quantile level only. Away from that level, and under any non-comonotonic copula, the three give different bids.
- **Your haircut raises α; it does not lower it, and only the bid level matters under netting.** With a=1 and comonotonic dependence, α′ = Φ(|z_ε| + h/σ_Z). A haircut h = 2.18 on X2 gives α′ = 0.95, the same as derating a_2 to 0.779. Under pooled delivery, where along the line a·x=c you put the point does not affect aggregate reliability. It only matters later, for disaggregation (joint survival probabilities, design points).
- **Put confidence into the worst case as a per-asset lower-tail ambiguity set.** Worst-case lower VaR is monotone in the marginal quantile functions on (0,ε). A KS band of width η_i therefore takes η_i out of the risk budget: with η_2=0.05, the bid falls from 15.335 to 13.41. That is still ≥90% reliable for any copula and any true F_2 within the band. η_i can be calibrated from probability integral transform (PIT) exceedance history, with credibility shrinkage and online adaptive conformal updates. Banking, actuarial science, forecast verification and structural reliability each supply part of this machinery.

## Conventions and notation
VaR_α(Z) = sup{z : P(Z ≥ z) ≥ α} = F_Z^{-1}(1−α), with ε = 1−α = 0.1 and z_ε = Φ^{-1}(0.1) = −1.2816. A higher VaR is better. φ(z_ε) = 0.1755, so φ(z_ε)/ε = 1.755. Worst case (written underline-VaR) is the infimum over the Fréchet class; best case is the supremum. Example marginals are X1~N(10,2²) (high confidence) and X2~N(15,4²) (low confidence). Their 10% quantiles are q1 = 7.437 and q2 = 9.874, so VaR^M = 17.311. Expected shortfall (ES) is the conditional tail mean; the lower-tail version is E[X | X ≤ F^{-1}(ε)]. A market time unit (MTU) is 15 minutes here. First-order stochastic dominance is written FSD.

## Part A — Rigorous analysis

### A1. Effect of a on Z and on VaR_α(Z) under the Gaussian model
For X ~ N(μ, Σ), Z = a·X ~ N(a·μ, aᵀΣa), so

VaR_α(a) = a·μ + z_ε σ_Z(a), σ_Z(a) = √(aᵀΣa).

**Gradient = Euler allocation = FORM design point.** ∂VaR/∂a_i = μ_i + z_ε (Σa)_i/σ_Z. For jointly Gaussian (X, Z), E[X_i | Z = c] = μ_i + Cov(X_i,Z)(c − μ_Z)/σ_Z². Setting c = VaR gives μ_i + z_ε(Σa)_i/σ_Z, so the gradient equals E[X_i | Z = VaR_α]. The most probable point (Mahalanobis-closest to μ) on the hyperplane {a·x = VaR} is x* = μ + z_ε Σa/σ_Z, the same vector. Its reliability index is β = |z_ε|. With independent components, x*_i = μ_i − β α_i σ_i, where α_i = a_iσ_i/σ_Z. This is exactly the design-value format of EN 1990 Annex C, discussed in Part B. By Euler's theorem Σ a_i x*_i = VaR_α(a).
- Example, a = 1. With ρ = 1, x* = (7.437, 9.874) = (q1, q2): your anchor point is the comonotonic design point. With ρ = 0, x* = (8.854, 10.415), which sums to 19.269.

**Homogeneity and true degrees of freedom.** VaR_α(ta) = t·VaR_α(a) for t > 0. Writing a = ‖a‖u, the direction u (in d=2, one angle) fixes the risk per unit; ‖a‖ only rescales the product. Your "three degrees of freedom (a1, a2, α)" therefore reduce as follows:
- α is fixed at the regulatory floor of 0.9, or chosen above it as a buffer.
- The reliability of a derated product with respect to its own random variable a·X does not depend on ‖a‖.
- The only physically meaningful output is the bid level c. The pair (a1, a2) is just one way of producing a c.

**What 0 < a_i < 1 means.** a_i X_i ~ N(a_iμ_i, a_i²σ_i²). It is a derating: the asset offers a fraction a_i of itself. You are right that a_i > 1 has no physical meaning in this sense. However, a_i < 1 shrinks the spread, which is the opposite of what low confidence (possible distribution shift or misquantification) should do.

**A scaling alone cannot represent a mean shift or misquantification.** At a single level ε:
- Derating gives the quantile a_i q_i.
- A mean shift δ_i gives q_i − δ_i.
- Variance inflation κ_i gives μ_i + z_ε κ_i σ_i.

They coincide at that level iff δ_i = (1−a_i)q_i, respectively κ_i = (μ_i − a_iq_i)/(|z_ε|σ_i). At other levels they differ. In particular, derating cannot reproduce a pure shift N(μ_i−δ, σ_i²) or a pure inflation N(μ_i, κ²σ_i²) for every u in (0,1), because a_i μ_i = μ_i − δ and a_i σ_i = σ_i cannot both hold with a_i ≠ 1.

### A2. Your geometric construction, made rigorous (and corrected)

**Geometry.** {x : a·x = c} is a level line of Z with normal a. Varying c moves it parallel to itself, and each c corresponds to one α under a given joint law. "Restricting to the a-vector through the origin" therefore adds nothing: every line with normal a is parallel, and c alone sets the level.

The level line is also not a survival contour. The joint survival set {x : P(X1 ≥ x1, X2 ≥ x2) = α} is a curve. At (q1,q2) its value is 0.9 under ρ = 1 but 0.81 under independence.

**Haircut formula.** Anchor at (q1, q2) on a·x = VaR^M = a1q1 + a2q2. Keep x1 = q1 and lower x2 by h. The new bid is c′ = VaR^M − a2h. Under the Gaussian model with correlation ρ,

α′(h) = P(Z ≥ c′) = Φ( [ |z_ε|(a1σ1 + a2σ2) + a2h ] / σ_Z(ρ) ), σ_Z(ρ)² = a1²σ1² + a2²σ2² + 2ρa1a2σ1σ2.

At ρ = 1, σ_Z = a1σ1 + a2σ2 and this reduces to your expected form α′ = Φ(|z_ε| + a2h/σ_Z). The haircut needed for a target α′ is

h = [ z_{α′} σ_Z(ρ) − |z_ε|(a1σ1 + a2σ2) ] / a2, where z_{α′} = Φ^{-1}(α′).

Example (a = 1, ρ = 1, σ_Z = 6):
- α′ = 0.95 needs h = 6(1.6449 − 1.2816) = 2.180. The x2 coordinate becomes 7.694, which is X2's 3.39% quantile.
- α′ = 0.99 needs h = 6.269.
- At ρ = 0 the formula gives h = −0.333 for α′ = 0.95. A negative value means VaR^M is already more than 95% reliable when the assets are independent.

**Correction 1: you do not "bring the α threshold down".** Lowering the bid raises α′, as your own phrase "alpha mass above the contour becomes larger" says. What goes down is the quantile level read off the low-confidence asset (10% → 3.39% for X2). Aggregate α goes up. These are the same operation seen from two sides, and the aggregate requirement α ≥ 0.9 never relaxes.

**Correction 2: under netting, where the point sits on the line is irrelevant.** P(a·X ≥ c′) depends only on c′. Lowering x2 by h at fixed x1 and lowering x1 by (a2/a1)h at fixed x2 give the same aggregate reliability. "We trust X1, so fix x1" has no aggregate meaning. It matters only for disaggregation, i.e. allocating the bid to assets through joint quadrant probabilities, design points or Euler shares.

**"Moving down x2" vs "shrinking a2".** Choosing a2′ = 1 − h/q2 makes q1 + a2′q2 = c′, so the two bids coincide at the calibration level. Example: h = 2.180 gives a2′ = 0.7792. The equivalence breaks in three ways:
1. At other levels ε: q1(ε) + a2′q2(ε) ≠ q1(ε) + q2(ε) − h.
2. Under any ρ < 1, because σ_Z changes. At ρ = 0, VaR_0.9(X1 + 0.7792X2) = 16.942, whereas VaR_0.9(X1 + X2) − h = 17.089, which is the mean-shift result.
3. Derating changes variance; a haircut is a pure location change.

### A3. Pointwise dominance
If X_i ≥ 0 almost surely and 0 ≤ a_i ≤ 1, then a·X ≤ 1·X pointwise. Hence P(1·X ≥ c) ≥ P(a·X ≥ c) for every c and every joint law, so VaR_α(a·X) ≤ VaR_α(1·X) under any copula, and the same holds for the worst and best cases.

A more general version needs only FSD. Suppose G_i ≤_st F_i, meaning G_i^{-1} ≤ F_i^{-1}. Couple X_i′ = G_i^{-1}(U_i) and X_i = F_i^{-1}(U_i) through the same uniform vector U. Then Σ X_i′ ≤ Σ X_i for every copula, so underline-VaR(G) ≤ underline-VaR(F).

The Gaussian model violates positivity, since P(X2 < 0) = Φ(−3.75) ≈ 8.8×10⁻⁵. The statement is then exact only conditional on positivity, or wherever the relevant quantiles are positive (q1, q2 > 0 here). For the real, strictly positive envelopes it holds exactly.

For d = 2, the worst case below uses quantiles only at levels in (0, ε]. Lower-tail order of the quantile functions on (0, ε] is therefore enough. That is why variance inflation, which is not FSD-ordered, is still conservative.

### A4. Confidence adjustment inside the unknown-dependence worst case (core result)
**Baseline worst case.** underline-VaR = sup_{a+b=ε}{F1^{-1}(a) + F2^{-1}(b)}. The first-order condition is σ2φ(z1) = σ1φ(z2), and Gaussian quantile functions are concave on u < 0.5, so the stationary point is the maximum.
- Result: 15.335 at (6.224, 9.112), with (a*, b*) = (0.0295, 0.0705).
- The VaR^M bid of 17.311 has worst-case shortfall probability 0.186, i.e. reliability 0.814. I obtained this by inverting the sharp d=2 bound: z1 = −1.620, z2 = −1.113.
- So h = 17.311 − 15.335 = 1.976 is needed before any confidence adjustment just to reach 0.9 for every copula. The comonotonic bid is not conservative.

**How each adjustment enters.** Replace F2^{-1} by the adjusted quantile function G2^{-1} inside the sup:
- Mean shift δ: G2^{-1} = F2^{-1} − δ. The argmax does not move and the bid falls by exactly δ.
- Derating a2: G2^{-1} = a2F2^{-1}. The first-order condition becomes σ1φ(z2) = a2σ2φ(z1), and the worst-case point shifts.
- Variance inflation κ: G2^{-1} = μ2 + κσ2z. The first-order condition becomes σ1φ(z2) = κσ2φ(z1).

Each adjustment below is calibrated so that X2's adjusted 10% quantile equals 7.694 (the 3.39% quantile, the h = 2.18 case): a2 = 0.7792, δ = 2.180, κ = 1.425.

| Adjustment of X2 | Bid, comonotonic | Bid, Gaussian ρ=0 | Bid, Fréchet worst case | Worst-case point (x1, x2); (a*, b*) |
|---|---|---|---|---|
| None | 17.311 | 19.269 | 15.335 | (6.224, 9.112); (0.0295, 0.0705) |
| Derate a2=0.7792 | 15.131 | 16.942 | 13.346 | (6.414, 6.932); (0.0365, 0.0635) |
| Mean shift δ=2.18 | 15.131 | 17.089 | 13.156 | (6.224, 6.932); (0.0295, 0.0705) |
| Variance ×κ=1.425 | 15.131 | 17.257 | 12.886 | (5.952, 6.935); (0.0215, 0.0785) |

Nominal reliability of each worst-case bid, evaluated with the reported (trusted) marginals:

| Worst-case bid | Comonotonic | Gaussian ρ=0 | Fréchet worst case |
|---|---|---|---|
| 15.335 (none) | 0.946 | 0.985 | 0.900 |
| 13.346 (derate) | 0.974 | 0.995 | 0.951 |
| 13.156 (mean shift) | 0.976 | 0.996 | 0.955 |
| 12.886 (variance) | 0.978 | 0.997 | 0.959 |

These are my calculations. Normal cdf values are to 3–4 decimals, and the bids should be read as accurate to about ±0.01.

Interpretation:
1. **Comonotonic additivity hides the difference.** All three adjustments give 15.131 under ρ = 1 because only the 10% quantiles enter.
2. **The ranking reverses between Gaussian and worst-case dependence.** At ρ = 0, derating is most conservative and variance inflation least. Under the worst-case copula it is the other way round: the Makarov sup reads the quantile functions at levels a*, b* < ε, where variance inflation has the fattest lower tail.
3. **The high-confidence asset's worst-case point moves.** x1* = 6.224, 6.414 and 5.952 across the three cases. "Fixing x1" is therefore not invariant once dependence is unknown.
4. **The buffer.** Adjusting X2 to its 3.39% quantile buys about 0.95–0.96 worst-case nominal reliability with respect to the reported marginals, and 0.97–0.98 under comonotonicity.

**Clean special case: a KS band subtracts from the risk budget.** Suppose asset i's true cdf may lie anywhere in {G : sup|G − F_i| ≤ η_i}, with physical lower support 0. The FSD-minimal member is G_i^L = min(1, F_i + η_i), whose quantile function is F_i^{-1}(u − η_i) for u > η_i and 0 otherwise. The worst case over the band and all copulas is then the worst case at the marginals G_i^L. This is exact, by the coupling argument in A3, because G_i^L belongs to the band. For d = 2, provided the optimum is interior,

underline-VaR_ε(band) = sup_{a′+b′ = ε−η1−η2}{F1^{-1}(a′) + F2^{-1}(b′)} = underline-VaR_{ε−η1−η2}(F1, F2).

Example: η1 = 0 and η2 = 0.05 give a bid of 13.412 at z1 = −2.164, z2 = −1.815, with (a*, b*) = (0.0153, 0.0347) in reported levels. This bid is
- 95% worst-case reliable under the reported marginals, and
- ≥90% reliable for every copula and every true F2 within η2 of the reported F2.

This turns "confidence" into a probability budget that can be checked against history. For d ≥ 3 the identity "worst case over the band = worst case at G_i^L" still holds, but the result is no longer a simple level shift. It has to be computed with the rearrangement algorithm (RA) of Embrechts–Puccetti–Rüschendorf (2013).\[1\]\[2\] For large d, the ES-type asymptote applies with the lower-tail ES of each G_i^L in place of the Gaussian term μ_i − σ_iφ(z_ε)/ε.

**Other ambiguity sets.**
- **Moments only.** With mean and variance known but nothing else, the worst lower ε-quantile is μ − √((1−ε)/ε)σ, the Cantelli/El Ghaoui–Oks–Oustry factor, which is 3 at ε = 0.1.\[3\]\[4\] With unknown dependence this gives 25 − 3·6 = 7. This is far too conservative, because moment information discards the shape.
- **Wasserstein-1 ball of radius θ.** The worst lower ε-quantile y*(θ) solves ∫_{F(y)}^{ε}(F^{-1}(u) − y)du = θ. This is my derivation from the partial-mass-transport structure; check it before relying on it.
- **Pointwise envelope.** For any set B_i, using the pointwise lower quantile envelope L_i(u) = inf_{G∈B_i}G^{-1}(u) inside the Makarov sup gives a valid but possibly non-sharp bound, since sup inf ≤ inf sup. It is sharp when some G in B_i attains L_i, as with KS bands or mean-shift families.

### A5. Recommended framework
The bid should be

c* = inf { VaR_α(Σ a_i X_i) : X_i ~ G_i ∈ B_i(η_i), copula arbitrary }, with α = 0.9.

Here B_i is a lower-tail ambiguity set around the asset's self-quantified envelope F_i, and η_i is that asset's confidence metric, estimated from delivery history.

**Established:**
- Sharpness and the formula for d = 2 (Makarov 1981; Rüschendorf 1982).\[5\]\[6\]
- FSD monotonicity, which gives the reduction to G_i^L.
- The RA for d ≥ 3.\[1\]

**My derivations, to be verified before publication:**
- The KS budget-subtraction identity.
- The Wasserstein quantile characterisation.

**Open questions:**
- Choosing η_i with serially correlated 15-minute data.
- Joint calibration across MTUs, given the three-month frequency check (Part B).
- Strategic misreporting: assets can game the band.
- Whether the d ≥ 3 large-portfolio asymptote under bands keeps the ES form.

**Disaggregation hand-off.** Euler shares or design points apply under a model. Under unknown dependence, the Makarov point (a*, b*) gives a natural per-asset reservation: x1* = 6.224 and x2* = 9.112 in the unadjusted case.

## Part B — Confidence/trust metrics from historical performance

**Energinet's P90 rule (primary text).** Section 10.2.2 of Energinet's "Prequalification of units and aggregated portfolios" (Doc. 13/80940-106, v2.1.2, valid from 5 Aug 2024) sets the following:\[7\]
- The forecast corresponds to "the 10th percentile of a probabilistic forecast".\[7\]\[8\]
- "forecast resolution must be at least same time resolution as MTU".\[7\]
- "three months' worth of forecast data must be provided".\[7\]\[8\]
- "the forecast value must be lower or equal to the actual for at least 90% of the time".\[7\]\[8\]\[9\]
- If the forecast is finer than the MTU, "the lowest value within each MTU must be used as the highest allowable capacity".\[7\]
- After prequalification, "Energinet will continuously monitor if the sold capacity is available, and if deemed necessary an audit of the forecast may be conducted."\[7\]

For aggregated portfolios, Section 10.1 says "The reserve providing entity must always provide a response that meets the technical requirements, while the individual resources in the group on their own do not necessarily have to".\[7\]\[8\]\[10\] Pooled compliance is therefore what counts, which matches the half-plane argument in A2.

The often-quoted ex-post three-month check with exclusion is not stated in v2.1.2. It appears in secondary DTU work: Gade, Bindner & Kazempour 2024 write "Energinet may check ex-post … e.g., three months", and Lunde, Damm, Gade & Kazempour ("Aggregator of Electric Vehicles Bidding in Nordic FCR-D Markets: A Chance-Constrained Program", IEEE Trans. Power Syst. 2025) write "Energinet typically conducts these ex-post checks in 3-month periods." The older v2.0 (valid from 1 Sept 2023), §9.2.1, says that "there is a 90% chance that the sold capacity or more is available. This is when the prognosis is assumed to be correct." It also states: "If a market participant repeatedly, in good faith, does not deliver the sold reserve-capacity, then the participant will be excluded from participating in the market, until an approved prognosis can be approved by Energinet", and "A minimum of three months of data is required." Energinet's Danish page lists a retired February 2025 version, so a newer document exists; I have not verified its wording.\[11\]

P90 is therefore a frequency-coverage test over a sequence of MTUs. For three months of 15-minute MTUs (n ≈ 8,832), the binomial standard deviation of the exceedance rate at p = 0.1 is about 0.0032, before correcting for serial correlation.

| Field / metric | How computed from history | How it modifies quantile / reliability | Fit for self-quantified, private envelopes, unknown dependence |
|---|---|---|---|
| Basel traffic light (BCBS Basel Framework MAR99, Table 2) | Count VaR exceptions in 250 days at 99% | Green 0–4 → plus 0.00; yellow 5/6/7/8/9 → 0.40/0.50/0.65/0.75/0.85; red ≥10 → 1.00, added to multiplier 3 (version in force until 31 Dec 2022).\[12\]\[13\]\[14\] Current MAR99.17 uses multipliers 1.50; 1.70–1.92; 2.00\[15\] | Good template: exception count → discrete conservatism multiplier; per asset or per pool; no dependence needed |
| Kupiec POF / Christoffersen conditional coverage (Kupiec 1995 J. Derivatives 3:73–84; Christoffersen 1998 Int. Econ. Rev. 39:841–862) | Likelihood-ratio tests: exception rate = ε (χ²(1)); independence of exceptions (Markov); combined LR_CC = LR_POF + LR_CCI | Test, not adjustment: decides whether η̂_i is significant; flags clustering (important for 15-min data) | Direct; uses only exceedance indicators 1{x_t < q_t}\[16\]\[17\]\[18\] |
| Bühlmann / Bühlmann–Straub credibility | Z_i = n_i/(n_i+k), k = expected process variance / variance of hypothetical means (Bühlmann & Straub 1970) | Shrink an asset's own exceedance rate toward the pool: p̃_i = Z_i p̂_i + (1−Z_i)p̄; set η_i from p̃_i | Excellent for new assets with short history; needs a pool of similar assets\[19\]\[20\] |
| Bayesian Beta–Binomial posterior | Prior Beta(a0,b0) on exceedance rate; s exceedances in n | Posterior mean = Z p̂ + (1−Z)·prior mean with Z = n/(n+a0+b0), the exact credibility form; an upper posterior quantile gives a one-sided η_i | Same data needs; gives an uncertainty statement |
| PIT calibration (Gneiting, Balabdaoui & Raftery 2007, JRSS-B 69:243–268) | u_t = F_{i,t}(x_{i,t}); probabilistic calibration ⇔ uniform PIT | Recalibrate: read reported envelope at level Ĥ_i^{-1}(ε), Ĥ_i = empirical PIT cdf; tail deficit η̂_i = max(0, Ĥ_i(ε)−ε) | Privacy-friendly: asset reports only a scalar u_t per MTU; one-sided lower tail is what matters\[21\]\[22\] |
| KS/DKW band (my proposal; DKW–Massart inequality) | Band half-width = observed tail deficit + sampling margin √(ln(2/β)/(2n)) (≈0.0145 for n=8,832, β=0.05, ignoring autocorrelation) | Exactly subtracts η_i from risk budget in d=2 worst case (A4) | Matches unknown-dependence framework exactly |
| Adaptive conformal inference (ACI; Gibbs & Candès 2021, arXiv 2106.00170; DtACI, JMLR 2024) | ε_{t+1} = ε_t + γ(ε − err_t), err_t = 1{miss} | Online quantile-level correction per asset; long-run coverage guarantee under arbitrary distribution shift; γ ~ size of shift | Ideal for nonstationary assets; per-asset ε_{i,t} feeds the Makarov sup\[23\]\[24\]\[25\] |
| Wasserstein distributionally robust optimization (DRO) (Mohajerin Esfahani & Kuhn 2018, Thm 3.5; Gao & Kleywegt 2023; Gade et al. 2024 for P90 aggregators) | Radius θ from finite-sample concentration (shrinks with N; curse of dimensionality) or cross-validation | Worst case over ball; in Gade et al., θ is the aggregator's "conservativeness" in a joint chance constraint | Works on joint data, which the user lacks; apply per marginal and combine with Makarov\[9\]\[26\]\[27\]\[28\]\[29\] |
| Moment / variance-constrained bounds (El Ghaoui, Oks & Oustry 2003, Oper. Res. 51:543–556; Bernard, Rüschendorf & Vanduffel 2017, JRI 84:923–959) | Bounds on means/covariances or on Var(sum) | √((1−ε)/ε) factor; a variance cap on the sum tightens dependence-free bounds (extended RA) | Moments alone are too loose (bid 7 in example); a pool-variance cap estimated from aggregate history is a useful add-on\[3\]\[30\]\[31\] |
| Robust distortion risk measures (Bernard, Pesenti & Vanduffel 2024, Math. Finance 34:774–818) | Wasserstein ball plus known mean/variance around a reference | Quasi-explicit worst-case quantile functions (VaR included, Sect. 4.1) | Per-marginal robust quantile function to insert into Makarov\[32\]\[33\]\[34\] |
| FORM / partial factors (EN 1990 Annex C; JRC 2024 report) | Design value x_d = μ − α_X β σ; fixed α_R=0.8, α_E=−0.7, β=3.8 (50 yr), non-dominant 0.4× | Fixed sensitivity factors let components be designed independently and conservatively | Direct analogue of per-asset haircuts without dependence knowledge; your x* formula is this format\[35\] |
| Capacity-market de-rating (GB) | Technology-class availability from historical data (e.g., Short-Term Operating Reserve (STOR) data for demand-side response (DSR)) | Literal a_i: the National Energy System Operator's (NESO) 2024 Electricity Capacity Report states a 79% DSR de-rating factor. Per the UK government call for evidence "Capacity Market: consumer-led flexibility", it rose "over 7% to 79% since the 2022 ECR, largely due to STOR moving from seasonal fixed contract procurement to a day-ahead auction" | Shows a_i is used in practice, but per class, not per asset history |
| Demand-response reliability index (Venizelou et al. 2023, Energies 16:2815) | Per-customer reliability and fairness indices from past events | Prioritises reliable customers in dispatch | Heuristic weighting; no quantile semantics\[36\] |
| Conformal for a Danish aggregator (Pipada Sunil Kumar, Pourmousavi, Liisberg & Lemos-Vinasco, arXiv 2601.14663, Jan 2026) | Monte Carlo dropout + conformal prediction | Calibrated intervals for an mFRR bidding model; their MCP variant kept 71/69/68% of perfect-information profit (unadjusted/hourly/daily penalty, Table 3), "while satisfying regulatory reliability constraints" | Evidence that conformal calibration works under Energinet rules |

## Recommendations
1. **Bid from the worst case, not from VaR^M.** Use the d=2 Makarov formula, or the RA for d ≥ 3. Drop a as a confidence device and keep a_i ∈ {0,1}, or fractional offers, for the separate decision of which assets to include.
2. **Define the confidence metric as a lower-tail coverage deficit.** For each asset, η_i = max(0, Ĥ_i(ε) − ε) plus a sampling margin, computed from PIT values over a rolling three-month window, which matches Energinet's prequalification window.\[8\]
3. **Shrink η_i for new assets** with Bühlmann credibility (or the Beta–Binomial equivalent) toward the pool, and run ACI online to track drift.
4. **Use the result directly.** In d=2 the bid is underline-VaR at ε − Σ η_i. Report the excess nominal reliability with respect to the reported marginals, which is 0.95 for η_2 = 0.05, and backtest the pool with Kupiec and Christoffersen tests.
5. **For disaggregation, start from the worst-case point (a*, b*)** of the adjusted problem. Do not start from (q1, q2): it is optimal only under comonotonicity, and it moves as η changes.

## Caveats
- All numeric results are my own calculations. Normal cdf values are to 3–4 decimals.
- The KS budget-subtraction identity and the Wasserstein quantile formula are derived here, not quoted from a source. Verify them before publication.
- Serial correlation in 15-minute data shrinks the effective sample size, so the DKW and binomial margins above are optimistic.
- The Energinet ex-post three-month exclusion rule comes from secondary sources. A post-February-2025 version of the prequalification document exists that I have not checked.\[9\]\[11\]
- The GB 79% DSR de-rating figure is taken from the UK government call for evidence "Capacity Market: consumer-led flexibility", which attributes it to NESO's 2024 Electricity Capacity Report; I have not read the report itself.
- Gaussian marginals violate positivity, so the pointwise-dominance result is exact only for the real, positive envelopes.

## Sources

1. [(PDF) Model Uncertainty and VaR Aggregation](https://www.researchgate.net/publication/257211953_Model_Uncertainty_and_VaR_Aggregation)
2. [Paul Embrechts The Rearrangement Algorithm](https://people.math.ethz.ch/~embrecht/ftp/Rearrangement_Algorithm.pdf)
3. [Worst-case reinsurance strategy with likelihood ratio uncertainty](https://www.cambridge.org/core/journals/astin-bulletin-journal-of-the-iaa/article/worstcase-reinsurance-strategy-with-likelihood-ratio-uncertainty/E20BD22100D9C1EC853947080C61BA1B)
4. [Worst-case VaR and robust portfolio optimization with interval random uncertainty set - ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0957417410005142)
5. [The beautiful art of rearranging matrices](https://link.springer.com/article/10.1007/s10479-025-06799-y)
6. [Model-free bounds on Value-at-Risk using extreme value information and statistical distances](https://arxiv.org/pdf/1610.09734)
7. <https://www.energinet.dk/media/ox0gqmvw/gaeldende-prequalification-of-units-and-aggregated-portfolios.pdf>
8. <https://en.energinet.dk/media/3obnf0q4/prequalification-of-units-and-aggregated-portfolios.pdf>
9. <https://arxiv.org/html/2404.12807>
10. <https://en.energinet.dk/media/r4rbirq2/prequalfication-of-fcr-n-and-fcr-d.pdf>
11. [Prækvalifikation og test](https://www.energinet.dk/el/balancering-og-systemydelser/adgang-til-systemydelsesmarkederne/praekvalifikation-og-test/)
12. [Application guidance | Bank for International Settlements](https://www.bis.org/committees/bcbs/basel-framework/standard/mar/99/inforce/2019-12-15/published/2019-12-15)
13. [Backtesting VaR](https://analystprep.com/study-notes/frm/part-2/market-risk-measurement-and-management/backtesting-var/)
14. [Backtesting VaR: How to Validate Your Risk Model](https://ryanoconnellfinance.com/backtesting-var/)
15. <https://www.bis.org/committees/bcbs/basel-framework/81307/chapter.pdf>
16. [Backtesting Lambda Value at Risk](https://arxiv.org/pdf/1602.07599)
17. [cc - Conditional coverage mixed test for value-at-risk (VaR) backtesting - MATLAB](https://www.mathworks.com/help/risk/varbacktest.cc.html)
18. [Overview of VaR Backtesting - MATLAB & Simulink](https://www.mathworks.com/help/risk/overview-of-var-backtesting.html)
19. [Interval Estimation of the Credibility Factor](https://variancejournal.org/article/141987-interval-estimation-of-the-credibility-factor)
20. [Scaling Laws, Tabular Data and Actuarial Ratemaking Models](https://arxiv.org/pdf/2609.03106)
21. [Probabilistic forecasts, calibration and sharpness](https://www.semanticscholar.org/paper/Probabilistic-forecasts,-calibration-and-sharpness-Gneiting-Balabdaoui/bfaf3af966e2fdd758dfd82b9732bc3831f664a3)
22. [Area-covering postprocessing of ensemble precipitation forecasts using topographical and seasonal conditions](https://arxiv.org/pdf/1912.11827)
23. [Conformal Inference for Online Prediction with Arbitrary Distribution Shifts](https://arxiv.org/pdf/2208.08401)
24. [Online conformal inference with retrospective adjustment for faster adaptation to distribution shift](https://arxiv.org/pdf/2511.04275)
25. [Adaptive Conformal InferenceUnder Distribution Shift](https://arxiv.org/html/2106.00170v2)
26. [Finite-Sample Guarantees for Wasserstein Distributionally Robust Optimization: Breaking the Curse of Dimensionality](https://arxiv.org/pdf/2009.04382)
27. [\[1505.05116\] Data-Driven Distributionally Robust Optimization Using the Wasserstein Metric: Performance Guarantees and Tractable Reformulations](https://ar5iv.labs.arxiv.org/html/1505.05116)
28. [Aggregator of Electric Vehicles Bidding in Nordic FCR-D Markets: A Chance-Constrained Program](https://www.researchgate.net/publication/390682328_Aggregator_of_Electric_Vehicles_Bidding_in_Nordic_FCR-D_Markets_A_Chance-Constrained_Program)
29. [Wasserstein Distributionally Robust Optimization and Its Tractable Regularization Formulation](https://link.springer.com/article/10.1007/s10957-025-02896-x)
30. [Value‐at‐Risk Bounds With Variance Constraints - Bernard - 2017 - Journal of Risk and Insurance - Wiley Online Library](https://onlinelibrary.wiley.com/doi/10.1111/jori.12108)
31. [Value-at-Risk Bounds with Variance Constraints by Carole Bernard, Ludger Rüschendorf, Steven Vanduffel :: SSRN](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2342068)
32. [Robust distortion risk measures - Bernard - 2024 - Mathematical Finance - Wiley Online Library](https://onlinelibrary.wiley.com/doi/10.1111/mafi.12414)
33. [Robust distortion risk measures - Bernard - 2024 - Mathematical Finance - Wiley Online Library](https://doi.org/10.1111/mafi.12414)
34. [Worst-case distortion risk measures of transformed losses with uncertain distributions lying in Wasserstein balls](https://resolve.cambridge.org/core/journals/astin-bulletin-journal-of-the-iaa/article/worstcase-distortion-risk-measures-of-transformed-losses-with-uncertain-distributions-lying-in-wasserstein-balls/38372AA02C3384976985EC50AB05CF19)
35. [Reliability background of the Eurocodes](https://eurocodes.jrc.ec.europa.eu/sites/default/files/2024-11/JRC_Reliability_report_final_23Oct2024_with%20identifiers.pdf)
36. [DSO-Aggregator Demand Response Cooperation Framework towards Reliable, Fair and Secure Flexibility Dispatch](https://doi.org/10.3390/en16062815)
