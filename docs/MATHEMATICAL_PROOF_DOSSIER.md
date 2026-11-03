# Mathematical verification dossier — fixed uncertain pole / Bernoulli erasures

Date 20 October 2026. This dossier is a checkable independent derivation and regression summary, **not a signed independent human proof certificate**. Main text retains full proofs. Numerical certificate statements are algorithmically exact **for rational encodings of the stored gains**, not error bounds for physical floating-point implementations.

## Definitions and quantifiers
Fixed `a∈[a1,a2]⊂(1,∞)` for an entire trajectory; iid `γ_k~Ber(q)` independent of initial state and noise, iid zero-mean `w_k` with variance `σ²>0`. When `γ_k=1`, the controller sees **exact state** `x_k` before choosing `u_k`. No observation of `a`. Information class `L1`: `u_k=-h_age x_latest` (with zero input before first receipt), all deterministic age gains `h_0,h_1,…` constant across receptions, even if explicitly selected for q. Policies cannot learn a from a time history in L1. Reception-indexed periodic gains and past-sample feedback are separate, strictly larger/different classes. `MS stable`: physical-time `sup_k E x_k²<∞`. Mean square *asymptotic zero* is not claimed in the presence of `σ²>0`.

## Complete key proof checks and boundary conditions

**Lemma necessity.** Conditioning on information available since the last received observation shows variance includes unrecoverable `σ² Σ_{j=0}^{k-1}[(1-q)a²]^j`. Therefore necessary `(1-q)a²<1`; equality also fails (linear divergence). The result applies to arbitrary causal controllers with the same delayed information constraint, not an oracle able to measure unreceived noise.

**Exact renewal condition, last reception.** Let the geometric gap `G>=1` have `Pr(G=g)=q(1-q)^{g-1}` and homogeneous multiplier `m_g(a)=a^g−Σ_{j=0}^{g−1}a^{g−1−j}h_j`. Accumulated noise variance is `n_g=σ² Σ_{j=0}^{g−1}a^{2j}`. The arrival-level recursion is `V_{n+1}=ρ V_n+ν`, with `ρ=Σ Pr(G=g)m_g²` and `ν=σ²/[1−(1−q)a²]`, provided latter finite. Physical-time renewal of moments is explicitly proved in LaTeX, and `sup E x_k²<∞` iff `ra²<1` and `ρ<1`. Importantly use **strict** inequality, and do not substitute a sampled finite horizon for this criterion.

**Hardy representation.** Set `r=(1−q)a²`, `rho_q=√(1−q)`, `P(w)=Σ_{j≥0}h_j w^j`, `p(z)=P(rho_q z)`, `σ_i=rho_q a_i`. The identity `ρ(a,q,h)=q||(p−a) k_{sigma}||²_{H²}` is obtained coefficient by coefficient from `Σ_{g≥1} e_g t^g=(t−Q(zt))/(1−t)`, with `Q(w)=wP(w)`, `e_g=a^{-g}m_g`, `z=1/a`. It remains a formal-series identity in `[0,+∞]`, avoiding silently assuming p already belongs to H².

**Stationary two-endpoint converse.** From `u_i=(p−a_i)k_{σ_i}` and each `||u_i||<1/√q`, evaluate `(1−σ_1 z)u_1−(1−σ_2 z)u_2=a_2−a_1` at `z=s=(σ_1+σ_2)/2`. The evaluation Cauchy-Schwarz inequality gives `d< 2√(1−s²)/√q`, hence `d² q < 4[1−(1−q)(a1+a2)²/4]`. Algebra is equivalent to `q>[(a1+a2)²−4]/(4a1a2)`. This is a lower bound **on class L1 only**; interval necessity follows from endpoints.

**Stationary constructive converse.** Let `bar=(a1+a2)/2`, `X=d²/4`, `h0=bar`, `hj=X*bar^{j−1}` for `j>=1`. For each a, compute `ρ(a,q,h)` as norm of the rational Hardy transfer. The detailed elementary comparison in the manuscript proves `ρ(a,q,h)≤qX/[1−(1−q)bar²]`, with equality at endpoints. Together with `(1−q)a2²<1`, this proves whole-interval sufficiency. No grid assumption occurs. This schedule is **feasibility-rate optimal in L1**, not quadratic-cost optimal at fixed q. If algebraic cutoff≥1 (`d≥2`), rate feasible set in `(0,1)` empty; infimum in extended reals is `+∞`, **not an attainable probability >1**.

**At q=1, d=2 or greater.** Single-sample feedback reduces to robust static scalar feedback, possible iff `d<2`, strict. The q=1 proposition for L2 relates to a *different* policy class; cannot project its width-4 condition back onto L1. Formula for Chebyshev impulse schedule at q=1 is interpreted as deterministic arrival-period cycles; for d≥4 the ordinary schedule cannot stabilize the entire continuum for any arrival q.

**Chebyshev periodic (restricted family).** Gains are `c_j=bar + (d/2)cos[(2j−1)π/(2m)]` at receptions, u=0 between arrivals. For fixed a, reception multiplier per phase `η_j=q(a−c_j)²/[1−(1−q)a²]`. Product contraction iff `Πη_j<1`. The ordinary Chebyshev root polynomial satisfies `sup_{a∈I}|Π(a−c_j)|^{2/m}=(d²/16)2^{2/m}=B_m`, with equality at a2, where `K(a,q)=q/[1−(1−q)a²]` is also maximal for a>1. Thus `q_m=(a2²−1)/(a2²−B_m)` if `B_m<1`, strict inequality q>q_m. If `B_m≥1`, even q=1 fails uniform robust contraction. This does **not** extend to arbitrary periodic gains or learning controllers.

**Two-measurement memory-three.** Rational Lyapunov endpoint conditions `P ≻0` and `P−E[C_G^T P C_G]≻0` at a=6/5 and a=3, q=9/10, reconstructed in `src/verify_memory3.py` with exact SymPy rationals. These establish simultaneous stabilization of the **set of two plants** only. Finite-grid checks do not establish an interval-wide claim. No blanket theorem for optimal two-sample memory follows.

**Exact periodic interval verifier.** For stored period-six gains, `src/verify_weighted_exact.py` computes the rational polynomial homogeneous period multiplier and proves `Π_j η_j(a,q)<1` for every a across each interval by rational polynomial inequalities, with the full continuum and strict inequality. Seeds and gains included. These are certificates of selected designs; numerical LP search alone does not certify feasibility beyond its nodes.

**Stationary exact asymptotic costs.** For `hj` from `h0=L`, `hj=B*c^{j−1}` for j≥1, `m_g=U a^{g−1}+V c^{g−1}` at a≠c, where `U=a−L−B/(a−c)` and `V=B/(a−c)`. Exact geometric summation gives
`ρ=q[U²/D_a + 2 U V/D_ac + V²/D_c]` where `D_a=1−ra²`, `D_ac=1−rac`, `D_c=1−rc²` and `r=1−q`. The removable center has exact limit
`ρ(c)=q[(c−L)²/D_a − 2(c−L)Brc/D_a² + B² r(1+rc²)/D_a³]`.
With ν=σ²/D_a and arrival variance `Z=ν/(1−ρ)`, exact physical-time limits:
`Ex²=(q+rρ)Z+rν`, `Eu²=qZ(L²+rB²/D_c)`.
Assumptions: strictly `ra²<1, rc²<1, ρ<1` and independent zero-mean finite-variance disturbances. Same noise and weights for both stationary policies. The 0.280345 and 0.339402 *worst-cost values* are finite **1001-point grid maxima of exact per-plant costs**, not continuum-certified maxima. The 3500×500 paired Monte Carlo supplies an independent **finite-time** sanity check, not proof of stability or asymptotic cost equality at that horizon.

**New exact continuum polynomial proof.** For rational a, L, B, c, q let `N=(a−L)(a−c)−B`. Then define
`H(a)=((a−c)² D_a D_ac D_c − q[N² D_ac D_c +2 N B D_a D_c + B² D_a D_ac])/(a−c)²`.
The numerator and first derivative vanish at a=c; exact symbolic polynomial division yields polynomial H (degree ≤3) and identity `H=(1−ρ)D_a D_ac D_c` on all a by continuity. If `r*a2²<1`, each D factor is positive on a∈[a1,a2] (with c∈ interval). On each of 64 rational subintervals, express H in Bernstein coefficients computed using Python `fractions.Fraction`. Exact positivity of **every** coefficient proves `H>0` for the **whole continuous interval**, not merely nodes. The full rational coefficient witness is exported at `results/cost_aware_bernstein_full_rational_certificate.json` with the exact degree-three polynomial H, every rational subinterval endpoint, and **all 256 exact rational Bernstein coefficients**, each explicitly positive. Therefore readers can check signs by integer arithmetic without running the optimizer. The exported decimal minimum coefficient is informative; the underlying exact positivity witness uses exact dyadic binary64 gain values. Different physical finite-precision gains require a separate perturbation margin.

**Formal independence and residual risk.** No external human mathematician or formal proof assistant has independently endorsed the paper. The archived auxiliary rate-aware/universality corollary cannot be treated as verified. A priority literature check and independent audit of the endpoint Hardy inequality, closed-form rational cost limits, and periodic renewal equivalence remain mandatory before the GO decision.

## Scripts and evidence
- `python -m pytest -q tests/` — regression and inconsistency guards.
- `python src/verify_memory3.py` — exact two endpoint Lyapunov matrices.
- `python src/verify_weighted_exact.py results` — exact period-six polynomial certificates.
- `python code/run_consolidation.py --profile standard` — local optimizer, rational Bernstein certificate, deterministic cost curves and paired MC.
- `results/cost_aware_stationary.json` — stored hexadecimal gain encodings, rational interval and certificate metadata.
- `results/weighted_periodic_exact_certificates.json` — period-six rational certificates.
- `manuscript/*.tex` — actual full proofs and mathematical statements.
