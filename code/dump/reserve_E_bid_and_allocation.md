# Adding firm reserve energy E to bidding and allocation

## 1. Scope and assumptions

An aggregator controls nonnegative random available energies
\(\mathbf X=(X_1,\ldots,X_d)\) over one delivery interval. It assigns deterministic energy targets \(\mathbf c\ge0\) and has firm, usable backup energy \(E\ge0\) from a battery. Reliability is \(\alpha=0.9\).

All quantities use the same energy units and delivery interval. The reserve is available after asset shortfalls are observed and can be directed to any asset's deficit. Battery power, efficiency, state of charge, and timing restrictions must either already be included in the definition of usable E or be modeled separately.

The main model covers each assigned asset deficit with the battery. It does not automatically transfer surplus from another asset to that deficit. Section 7 gives the alternative when such pooling is allowed.

## 2. Baseline: the AND survival constraint

Without backup, simultaneous delivery of all assigned targets is the event

\[
A_{\mathbf c}=\bigcap_{i=1}^d\{X_i\ge c_i\}.
\]

The feasible allocation region and its 90% contour are

\[
\mathcal C_0=\{\mathbf c\ge0:\Pr(A_{\mathbf c})\ge\alpha\},
\qquad
S_0=\{\mathbf c\ge0:\Pr(A_{\mathbf c})=\alpha\}.
\]

The inequality defines the optimization domain. For continuous distributions, an attained maximum will commonly lie on the probability contour. With atoms or discretized probabilities, equality need not be attainable, so retain the inequality in the optimization.

For a chosen allocation, \(\sum_i c_i\) is deterministic. The probability concerns whether the random capacities meet those targets.

## 3. Relaxing AND with a shared reserve

Define the total assigned shortfall

\[
D_{\mathbf c}(\mathbf X)=\sum_{i=1}^d(c_i-X_i)_+,
\qquad (z)_+=\max(z,0).
\]

The battery restores every assigned target if and only if

\[
D_{\mathbf c}(\mathbf X)\le E.
\]

Equivalently, delivered energy credited to the assigned targets is

\[
\sum_i\min(X_i,c_i)+E\ge\sum_i c_i.
\]

This expression treats E as available backup; unused reserve need not be discharged.

The relaxed feasible region and contour are

\[
\boxed{\mathcal C_E=\{\mathbf c\ge0:
\Pr(D_{\mathbf c}(\mathbf X)\le E)\ge\alpha\}},
\]

\[
\boxed{S_E=\{\mathbf c\ge0:
\Pr(D_{\mathbf c}(\mathbf X)\le E)=\alpha\}}.
\]

At E=0 this recovers the original AND condition exactly. Increasing E enlarges the feasible region:

\[
E_2\ge E_1\quad\Longrightarrow\quad
\mathcal C_{E_1}\subseteq\mathcal C_{E_2}.
\]

For example, targets (6,4) and capacities (5,3) give a shortfall of 2. They fail the original AND event but deliver the full target of 10 with E=2.

## 4. (a) Optimal bid selection

For targets that include the energy the battery may need to supply, select

\[
\boxed{
b^*(E)=\max_{\mathbf c\ge0}\left\{
\sum_i c_i:
\Pr\!\left(\sum_i(c_i-X_i)_+\le E\right)\ge0.9
\right\}.
}
\]

The maximizing vector \(\mathbf c^*(E)\) determines both the bid and its allocation. The reserve enters the probability constraint, and the bid is \(\sum_i c_i^*(E)\). Do not add E to this sum again.

Geometrically, maximize the linear function \(\sum_i c_i\) over \(\mathcal C_E\). In two dimensions, move a line \(c_1+c_2=b\) outward until it reaches the feasible boundary. Do not assume this region is convex or that a local stationary point is the global maximum.

### Conservative additive benchmark

Let \(\mathbf c^0\) be a baseline feasible allocation. Keeping its asset targets and scheduling the firm battery separately gives the reliable benchmark

\[
b_{\mathrm{add}}(E)=\sum_i c_i^0+E.
\]

If the baseline solution is optimal, this is \(b^*(0)+E\). With unrestricted nonnegative backed targets, it is also feasible in the relaxed model: choose a deterministic vector \(\mathbf r\ge0\) with \(\sum_i r_i=E\) and set \(\mathbf c=\mathbf c^0+\mathbf r\). Whenever the baseline AND event holds, the new shortfall is at most E. Therefore,

\[
b^*(E)\ge b^*(0)+E
\]

under these assumptions. Additional target caps or assignment restrictions must be checked separately before using this bound.

Joint optimization can improve on the additive benchmark because the battery can cover different deficits in different realizations. For illustration, let capacities be (10,0) or (0,10), each with probability 0.5. The 90%-reliable AND bid without reserve is 0. With E=2, targets (2,2) succeed in both states, supporting a bid of 4 rather than the additive benchmark of 2. This is a discrete illustration, not a numerical result from Case5.ipynb.

## 5. (b) Allocation for a given bid

For a specified bid b and reserve E, find

\[
\boxed{
\mathbf c\ge0,\qquad \sum_i c_i=b,\qquad
\Pr\!\left(\sum_i(c_i-X_i)_+\le E\right)\ge0.9.
}
\]

Equality assigns exactly the bid. Using \(\sum_i c_i\ge b\) is a sufficient formulation but may assign unnecessarily high targets. Under this model, targets can be reduced to sum to b without worsening reliability.

These conditions define a feasible allocation. An optimal allocation requires a separate objective. Reuse the allocation objective in the existing case study if it has one; examples include minimizing a documented allocation cost or minimizing reserve required for a fixed bid.

One useful additional quantity is the least reserve required for an allocation:

\[
E_{\mathrm{req}}(\mathbf c)
=\inf\{e\ge0:\Pr(D_{\mathbf c}(\mathbf X)\le e)\ge0.9\}.
\]

Then the reserve needed for a bid is

\[
E_{\min}(b)=\min_{\mathbf c\ge0,\,\sum_i c_i=b}
E_{\mathrm{req}}(\mathbf c).
\]

If the minimum is attained, a bid b is feasible with backup E exactly when \(E_{\min}(b)\le E\). This is an allocation objective with a clear operational meaning.

## 6. Geometry and connection to the hazard-scenario discussion

For fixed targets c, the success region in realized-capacity space is

\[
G_{\mathbf c,E}
=\{\mathbf x\ge0:D_{\mathbf c}(\mathbf x)\le E\}.
\]

Its complementary risk region is

\[
H_{\mathbf c,E}
=\{\mathbf x\ge0:D_{\mathbf c}(\mathbf x)>E\},
\qquad \Pr(\mathbf X\in H_{\mathbf c,E})\le0.1.
\]

Do not confuse this boundary in realized-capacity space with \(S_E\), the 90% contour in allocation space.

The success region has a union-of-AND representation:

\[
G_{\mathbf c,E}
=\bigcup_{\substack{\mathbf e\ge0\\\sum_i e_i\le E}}
\{\mathbf x\ge0:\mathbf x\ge\mathbf c-\mathbf e\}.
\]

The reserve allocation e may differ between realizations. The budget is shared: assigning E independently to every asset would incorrectly allow up to dE total reserve.

This geometry resembles the Survival Kendall union of AND regions discussed by Salvadori et al. However, here the candidate thresholds are determined by a physical energy budget. A Kendall probability level is not an amount of energy, and no identification between E and that level is assumed. The reserve construction does not require a Kendall model.

Reference: Salvadori, Durante, De Michele, Bernardi, and Petrella (2016), A multivariate copula-based framework for dealing with hazard scenarios and failure probabilities, Water Resources Research, DOI: 10.1002/2015WR017225, especially Sections 2-3.1. The reserve model above is our proposed application, not a claim that the paper derives it.

## 7. Alternative when asset surplus can be pooled freely

If any asset's surplus can offset another asset's deficit, and all available energy can be dispatched, the success condition instead is

\[
\sum_iX_i+E\ge b.
\]

The optimal pooled bid is

\[
b_{\mathrm{pool}}^*(E)
=\sup\{b:\Pr(\sum_iX_i+E\ge b)\ge0.9\}
=b_{\mathrm{pool}}^*(0)+E.
\]

For targets summing to b, assigned-deficit success implies pooled success, but the converse need not hold. Keep the two operational models separate in calculations and plots. The assigned-deficit model is relevant when original targets remain in force and the battery supplies their deficits; the pooled model assumes sufficient information and redispatch to use surplus.

## 8. Unknown dependence

The previous probabilities are evaluated under a specified joint law P. If only marginal distributions and partial dependence information are known, define the admissible set \(\mathcal P\) explicitly and replace the constraint by

\[
\boxed{
\inf_{P\in\mathcal P}
P(D_{\mathbf c}(\mathbf X)\le E)\ge0.9.
}
\]

Use the same robust constraint in both bid selection and allocation. A calculation for one copula or a finite collection of sampled copulas is not a guarantee over all admissible dependence structures. The reserve E is taken as firm; uncertain reserve requires an additional random variable and its dependence assumptions.

## 9. Case-study outputs

Compare the baseline AND solution, the conservative additive benchmark, the optimized assigned-deficit reserve solution, and a clearly labeled pooled benchmark where appropriate. Report bid, allocation, E, success probability, reserve use, and shortfall beyond reserve. Show \(b^*(E)\), the allocation-space feasible regions, and a fixed-bid allocation study. Distinguish illustrative dependence assumptions from certified robust results.
