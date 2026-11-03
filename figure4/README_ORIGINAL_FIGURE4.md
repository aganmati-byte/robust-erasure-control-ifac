# Figure 4 revision for coauthor review — 2 November 2026

Five equally spaced plant poles now appear for **all** q in {0.92, 0.96, 0.98}. For q=0.98, the stationary baseline is omitted at the endpoint poles only because its fixed-plant MS test fails there (rho ≈ 2.689). The two-received mean cost at a=4.5 is 68021.1479 and is flagged with an off-scale marker; the full-scale figure and complete CSV retain its value.

## Files
- `fig03_finite_horizon_mc.png`: ready to upload to Overleaf (same figure filename).
- `fig03_finite_horizon_mc_fivepoint_fullscale.png`: full-range figure with all results.
- `figure4_fivepoint_SE_CI95.csv`: 73 controller-plant combinations with the sample mean, sample standard deviation, standard error and pointwise bootstrap percentile 95% CI.
- `figure4_trajectory_costs.npz` and `figure4_trajectory_cost_index.csv`: 1,200 trajectory costs per combination, enabling alternate figures, custom confidence intervals, and paired policy comparisons without rerunning simulation.
- `figure4_paired_differences_optimized_minus_cheb.csv`: paired trajectory mean difference and uncertainty, for exploratory analysis.
- `figure4_all5_with_uncertainty.py`: self-contained reproduction script; requires the included `bench/...` source and stored design gains. Run `python figure4_all5_with_uncertainty.py` after installing dependencies in `requirements_figure4.txt`.

## Protocol and limitations
- N=1200 independent trajectories, T=300 steps, process noise N(0,0.01), x0~N(0,1), quadratic cost weight 0.05. Bernoulli arrivals with rate q.
- Identical packet sequences, process noise and x0 across all controllers at each pole; original seed offsets retained for j=0,2,4. The 43 previously reported values reproduce to floating-point precision.
- Per-policy standard error = sample SD of the 1,200 trajectory costs / sqrt(1,200). Shaded bands are pointwise bootstrap percentile 95% confidence intervals, B=400 resamples per combination, and are **not** ±SE bands.
- For unstable policies, heavy tails/rare events can make Monte Carlo standard errors and bootstrap intervals unreliable. Do not infer MS stability or asymptotic cost from these finite-horizon estimates.
- The plot does not write the word 'mean' on its axes, but its point estimates are still arithmetic sample means. Do not conceal this in the caption or Methods.
