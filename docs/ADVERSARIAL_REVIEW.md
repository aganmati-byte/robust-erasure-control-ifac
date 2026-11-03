# Independent adversarial editorial review — consolidated release

**Recommendation: NO-GO at present; reserve a distinct final scientific priority check, rather than initiating Elsevier submission.**

## Strengths, supported by actual results
- Tight, deliberately **controller-class-restricted** mean-square threshold formula, with a written Hardy-space converse and continuum constructive proof.
- Faithful treatment of impossible stationary intervals (algebraic cutoff≥1 is not a feasible probability).
- Exact physical-time expected quadratic state/input costs, not misidentified with finite-horizon Monte Carlo.
- A cost-aware stationary design on the same plant/noise assumptions, showing a measured 0.339402→0.280345 reduction of the *1001-grid maximum of exact per-plant cost* at q=.96; this is not a global cost-optimality claim.
- Separate exact-rational continuum certificates, 19 unit/regression tests, documented seeds, independently verified two endpoint rational Lyapunov solutions.
- More scientifically fair benchmark by separating oracle known-a, nominal deadbeat, stationary, periodic and limited past-measurement baselines.

## Highest risk: classical prior-art subsumption
Silva and Pulgar (2011), Automatica 47:1729–1736, additionally give a general necessary-and-sufficient mean-square stability characterization and an SNR-constrained equivalent control problem for scalar erasure channels. Elia 2005 solves mean-square stabilization of fixed LTI plants over fading/erasure architectures with robust synthesis and necessary/sufficient conditions. Tannenbaum 1980 and Khargonekar–Poolla–Tannenbaum 1985 predate the present robust margin/interpolation and periodic phenomena. Okano–Ishii 2012/2013/2017 studied uncertain scalar quantized/lossy systems. A credible paper must prove that its exact interval-specific restricted-class formula is not a previously derived corollary, **not merely note differences of title/assumptions**. Accessible evidence neither proves identity nor clears priority. Require a specialist priority certificate and subscription-level snowballing.

## Remaining technical limitations
1. The cost-aware SLSQP optimization is a finite-grid **local** solve; positivity on a continuum is rigorously certified but its worst-cost global interval maximum is sampled, and its optimality across L1 is not claimed.
2. Competing optimized periodic polynomial is designed on finite grid and only the selected numerical polynomial has a full interval certification. Long transients and actuator peaks remain important.
3. The two-sample memory-three certificate concerns **two poles, not an interval**; extrapolation is forbidden.
4. No independent human proof referee has attested the original and new arguments. No explicit hardware roundoff robustness; treating machine-stored gains as ideal dyadic rationals changes the mathematical object being certified.
5. An independent Claude manuscript is not available for overlap comparison. Manuscript self-declaration and author identity remain human responsibility.
6. The live detailed *journal-specific* Guide for Authors was not fully retrievable, though the Elsevier elsarticle manuscript compiles.
7. Some source references are only available as preprint abstracts or previews; a full theorem-level comparison for every cited recent paper is not claimed.

## Negative-result reporting
- The feasibility-optimal stationary policy is not cost-optimal at a fixed feasible q. The locally optimized policy redistributes costs: for a=1.2 the cost becomes higher, while the worst tested values decrease.
- Exact threshold for ordinary Chebyshev is **not** the global periodic minimum and generally exceeds known-a limit.
- Fixed-gain memory-two benchmark in v3 was MS unstable for challenging examples even if horizon-limited costs looked finite; don't cite such traces as stable.
- Artificial physical benchmark was deliberately omitted because a closely matched published plant/quantizer/reference policy was not source-verified.

## Concrete gates for reconsidering submission
Obtain (i) human theorem-level comparison against Elia (2005), Okano–Ishii and interpolation/simultaneous-stabilization corollaries; (ii) independently checked algebraic proof dossier; (iii) complete external time-varying-paper overlap check; (iv) author, disclosure and journal-format approvals. If prior art subsumes stationary result, shift to a genuinely nontrivial certified cost or policy-family theorem and provide an unmistakable new theorem before submission.
