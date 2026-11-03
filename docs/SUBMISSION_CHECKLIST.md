# Submission readiness checklist — 20 October 2026

Statuses: **COMPLETE** = executable/artifact verified within stated scope, **EXTERNAL** = needs human, institutional, publisher or subscription verification, **BLOCKING** = prevents an ethical GO now.

| Requirement | Status | Evidence / next action |
|---|---|---|
| Preserve original v3 unchanged | COMPLETE | `source_v3_unchanged/`, checksum manifest on release |
| Fix periodic zero-log convention and stationary gain-sequence bug | COMPLETE | v3 `src/robust_erasure.py` and tests |
| Whole-interval weighted periodic stability | COMPLETE | exact rational stored-gain verification, three intervals; no global optimality |
| Two-observation memory-three endpoint Lyapunov | COMPLETE | rational certificates for two endpoint plants; not all a |
| Cost-aware stationary finite-dimensional comparison | COMPLETE | local 31-node SLSQP, exact cost identities, interval-wide exact Bernstein certificate (including all positive rational coefficients), separate 1001-grid cost estimates |
| Reproducible tests, figures, CSV and Colab | COMPLETE | package, standard profile saved, quick-mode notebook smoke execution |
| Written mathematical theorem proofs | COMPLETE WITH LIMITATIONS | complete LaTeX and dossier, code-independent arguments for core theorems; some computer-assisted certificates require reading rational outputs |
| Independent human math review | EXTERNAL / BLOCKING | have a qualified control theorist check converse, renewal assumptions, strict boundaries and evidence |
| Exact bibliography priority of stationary L1 theorem | BLOCKING | could be implied by Elia (2005), Tannenbaum (1980), classical simultaneous/Hardy interpolation; no definitive priority review |
| Scopus / Web of Science comprehensive citation search | EXTERNAL | subscriptions unavailable; do not claim searched |
| Journal-specific Guide for Authors verification | EXTERNAL | official journal webpage not fully accessible (403); Elsevier elsarticle used as best-supported draft |
| Professional LaTeX/PDF | COMPLETE AS PREPARATION | compiled 22-page elsarticle preprint; PDF and figures checked |
| Figure / data reproducibility and defined evidence tiers | COMPLETE | original v3 benchmark + new stationary extension; distinct seeds and classes |
| Fair broad SOTA comparison incl adaptive/MPC | DISCLOSED SCIENTIFIC LIMITATION | fundamentally different information and inference assumptions; no false comparison |
| Quantified finite-precision actuator robustness | SCIENTIFIC LIMITATION | exact dyadic interpreted as ideal real gain; perturbation-to-hardware margin not established |
| Reference journal-format details, highlights, metadata | EXTERNAL | publisher submission portal checklist needed |
| All author approvals, contributor identities, conflicts, funding | BLOCKING / EXTERNAL | user and any genuine collaborators must attest; do not invent |
| Honest AI assistance disclosure and publisher policy | EXTERNAL | draft in manuscript and template; verify required current form |
| Text/proof/experiment overlap with separate Claude paper | BLOCKING / EXTERNAL | external draft not available; human side-by-side check needed |
| Cover letter | COMPLETE DRAFT | must be edited and signed only after blocking issues solved |

## Decision: **NO-GO FOR SUBMISSION TODAY**

The mathematical and executable package has materially improved, but **scientific priority of the headline exact stationarity threshold has not been established**, and the paper's distinctive result may be an explicit specialization of classic robust stochastic control synthesis. This is a decisive scientific issue, not merely an inaccessible administrative checkbox. No clear counterexample to the scoped stationary derivation was found during numerical testing; it is not independently peer-certified. This NO-GO is not a recommendation to discard the work: continue by an expert theorem-to-theorem mapping and strengthen/reframe the contribution if the principal result is already covered, before submitting to IFAC JSC.
