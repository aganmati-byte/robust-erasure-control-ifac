# Claim-to-evidence traceability — manuscript consolidation

| Claim in manuscript | Mathematical / machine-readable support | Evidence tier |
|---|---|---|
| Missing measurements imply lower bound `(1−q)a²<1` | Lemma `lem:any` proof in manuscript; elementary conditional noise-variance bound | Analytic theorem; human check pending |
| Stationary renewal MS iff `ra²<1` and `ρ<1` | `lem:ms` and `eq:renewal` with full proof; `src/robust_erasure.py` consistency tests | Analytic proof |
| Hardy H² identity | `lem:circle` coefficient-by-coefficient derivation | Analytic proof |
| Stationary critical `max(known_a2, interval)` for class L1 | `thm:main` endpoint converse and constructive continuum proof; `tests/test_robust_erasure.py` regression, original `thresholds.csv` | Analytic class-restricted theorem; **priority unknown** |
| Width ≥ 2 precludes stationary robust MS | `thm:main` algebraic inequality and empty-set convention | Analytic corollary |
| Memory-two works for d<4 when q near 1 | `prop:near1` compactness of second-moment operator near q=1 | Analytic existence proof, no explicit q0 |
| Memory-three two endpoint design a=1.2 and a=3 at q=.9 | `prop:cert` and `src/verify_memory3.py`, `results/memory3_exact_endpoint_certificates.json` | Exact rational **endpoint only** certificate |
| Ordinary Chebyshev exact family threshold | `prop:cheb-periodic` and Chebyshev polynomial identity; `results/thresholds.csv` | Analytic restricted-schedule threshold |
| Optimized period-six gains stabilize three continua | `src/verify_weighted_exact.py`, `results/weighted_periodic_exact_certificates.json`, `results/weighted_periodic_designs.csv` | Exact rational continuum for **specific stored idealized gains** |
| All v3 empirical comparison simulations and bootstrap intervals | `run_benchmark.py --profile standard`; `results/controller_performance.csv`, `results/paired_bootstrap_finite_horizon.csv`, `results/run_metadata.json` | Finite horizon empirical, not asymptotic stability proof |
| Exact stationary physical-time state+input moments | `prop:cost`, `prop:geomcost`, `code/cost_aware_stationary.py` and 20 unit tests | Analytic identity, fixed plant |
| Cost-aware gain design with max 1001-grid J ≈.280345 vs .339402 | `code/run_consolidation.py`; `results/cost_aware_stationary.json`, `results/cost_aware_bernstein_full_rational_certificate.json`, `results/cost_aware_exact_curves.csv`, `fig05_cost_aware_stationary.png` | **Exact cost per grid plant**, finite-grid extrema; local SLSQP |
| Whole-interval MS stability of cost-aware stationary gain | `prop:geomcost`, `bernstein_certificate` in `code/cost_aware_stationary.py`, exact positive Bernstein dyadic rational coefficients, 64 cells | Exact rational strict continuum stability; no finite-precision hardware guarantee |
| Independent finite-horizon cross-check of new policy | `results/cost_aware_paired_simulation.csv`, seed 20261064, 3500 paths horizon 500 for standard run | Finite-horizon Monte Carlo |
| Foundational and 2023–2026 comparisons | `docs/NOVELTY_MATRIX_2026.md`, `docs/VERIFIED_BIBLIOGRAPHY.md` | Public source metadata, selected accessible theorem text; incomplete subscription priority review |
| Separate Claude work non-overlap | `docs/OVERLAP_ASSESSMENT.md` | **Not assessed directly:** external manuscript unavailable |
| Scientific submission readiness | `docs/SUBMISSION_CHECKLIST.md`, `docs/ADVERSARIAL_REVIEW.md` | **NO-GO**, priority remains decisive |

**Never present** a numerical dense grid as a continuous-interval certificate, or finite horizon paired Monte Carlo as an asymptotic mean-square stability proof, or an endpoint certificate as robustness on the whole interval.
