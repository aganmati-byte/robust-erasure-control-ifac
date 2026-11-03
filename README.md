# Robust erasure control — IFAC collaboration (private pre-publication)

**Status: research draft / not an approved submission.** Shared Overleaf manuscript is authoritative. The GitHub repository is for reproducible algorithms, certificates, experiments and coauthor collaboration; it is intentionally separate from the time-varying/irregular-observation project.

**Scientific scope.** Fixed but unknown scalar pole on an interval, independent Bernoulli *sensor-to-controller* packet erasures, reliable control actuation, age-indexed stationary feedback and periodic alternatives. Finite-horizon Monte Carlo costs, exact moments and continuum-stability certificates are different evidence types.

## Reproduce on CPU

Python 3.10+ recommended. From repository root:

```bash
python -m pip install -r requirements.txt
python -m pytest -q tests/ figure4/test_figure4_archived.py
python src/verify_memory3.py
python src/verify_weighted_exact.py results
```

The above validates the archived numerical and rational evidence **without** re-running the expensive standard benchmark.

Run the complete standard analysis (can overwrite existing `results/`):

```bash
python run_benchmark.py --profile standard
python code/run_consolidation.py --profile standard
python figure4/recompute_figure4.py
```

Run `python figure4/recompute_figure4.py` **after** the main benchmark, since it requires `results/controller_performance.csv` and stored optimized gains. Figure 4 uses 1,200 trajectories of length 300, with 400 bootstrap resamples at each plant/controller combination. The figure shows 95% **pointwise** bootstrap intervals on *finite-horizon mean costs*; it is not a stability certificate. Unstable cases may have unreliable sampling intervals.

Output: `figure4/fig03_finite_horizon_mc_fivepoint_readable.png`; the current manuscript asset is `figures/fig03_finite_horizon_mc.png`.

### Notebook

- `notebooks/IFAC_core_preNov02_CPU_COLAB.ipynb`: earlier self-contained original CPU notebook for the core experiment, **pre-dating the November 2 Figure 4 update**.
- `notebooks/IFAC_Figure4_replot_saved_CIs.ipynb`: quick Colab plot regeneration from the embedded already-computed Figure 4 CSV.
- The full Figure 4 Monte Carlo can be recomputed from `figure4/recompute_figure4.py` (notebook visualization and full simulation are different tasks).

## Folder guide

- `src/`, `code/`, `tests/`: analytical methods, routines and automated checks.
- `results/`: stored original standard-profile numerical tables and exact certificates.
- `figure4/`: **November 2** five-pole raw data, bootstrap errors and regeneration script.
- `figures/`: current display figures, including corrected Figure 4.
- `manuscript/`: dated LaTeX snapshot only; use shared Overleaf for latest text.
- `docs/`: mathematical proof dossier, traceability, novelty review and submission gate checklist.
- `zenodo/`: metadata and instructions for an **unpublished Zenodo draft**; not a published record.

## Scientific / editorial safeguards

The stationary critical rate is exact **only for the stationary last-received-state class**. Weighted-periodic certificates apply to specified gains stored as exact dyadic rational numbers; they are not proofs of global optimality. The 2-observation memory-three Lyapunov evidence is *endpoint-only*. Publication priority relative to classical robust-synthesis literature remains subject to scientific assessment. See `docs/RELEASE_REVIEW_GATES.md`.

**Coauthor review:** Proposed authors are Abdelilah Ganmati and Alexander Vasilyev (corresponding author for the article, subject to final approval). No public release, Zenodo publication, or journal submission has been performed.
