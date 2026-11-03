"""Figure 4 follow-up to coauthor Alexander Vasilyev, 2026-11-02.
Reproduce all pre-existing 43 original 3-grid points, and 2 new quarter poles for
q=.92/.96, under identical CPU simulation. Save per-architecture SE and
pointwise bootstrap 95% confidence intervals. Does NOT certify stability.

Run from repository root using: python figure4/recompute_figure4.py
"""
import sys, json, pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'src'))
from robust_erasure import sim_one, cheb_nodes, stationary_rho_exact
SEED=20261020
N=1200
T=300
B=400
SCENARIOS=[(1.2,3.0,0.92),(1.2,3.0,0.96),(1.5,4.5,0.98)]
POLICIES=['oracle','stationary','cheb','optimized','two_received']
LABEL={'oracle':'Known-plant oracle','stationary':'Stationary','cheb':'Chebyshev','optimized':'Weighted periodic','two_received':'Two received states'}
designs=pd.read_csv(HERE.parent/'results/weighted_periodic_designs.csv')
# Old 3-point experiments used j=0/2/4 with seed offsets 0/1/2.
# Preserve these exact random streams. New j=1/3 have offsets 3/4.
seed_offset={0:0,1:3,2:1,3:4,4:2}
records=[]
raw_costs={}
for lo,hi,q in SCENARIOS:
    row=designs[(np.isclose(designs.a1,lo))&(np.isclose(designs.a2,hi))&(np.isclose(designs.q,q))].iloc[0]
    optimized=np.asarray(json.loads(row.nodes),float)
    cheb=cheb_nodes(lo,hi,6)
    for j in range(5):
        a=float(lo+(hi-lo)*j/4)
        seed=SEED+int(100*lo)+int(1000*q)+seed_offset[j]
        rng=np.random.default_rng(seed)
        packets=rng.random((N,T))<q
        noise=rng.normal(0,.1,(N,T))
        x0=rng.normal(0,1,N)
        for policy in POLICIES:
            if policy=='stationary' and q==0.98 and j in (0,4):continue
            nodes=optimized if policy=='optimized' else cheb if policy=='cheb' else None
            values=np.asarray(sim_one(a,q,lo,hi,packets,noise,x0,policy,nodes=nodes)['J'],float)
            if not np.all(np.isfinite(values)):raise ArithmeticError((lo,hi,q,a,policy))
            raw_costs[f"traj_{len(records):03d}"]=values
            # Seed distinct by scenario, pole and policy; reproducible trajectory bootstrap
            bid=POLICIES.index(policy)
            bs=np.random.default_rng(seed+200_000+bid)
            # Use bounded peak memory via 4 groups of 100 bootstraps
            boot=[]
            for batch in range(4):
                inds=bs.integers(0,N,size=(100,N))
                boot.extend(values[inds].mean(axis=1).tolist())
            qlo,qhi=np.quantile(boot,[.025,.975])
            records.append({'a1':lo,'a2':hi,'q':q,'a':a,'grid_index':j,'policy':policy,
                'seed':seed,'N':N,'T':T,'mean_cost':values.mean(), 'median_cost':np.median(values),
                'sample_SD':np.std(values,ddof=1),'sample_SE':np.std(values,ddof=1)/np.sqrt(N),
                'bootstrap_CI95_low':qlo,'bootstrap_CI95_high':qhi,
                'CI_method':'percentile_bootstrap_of_trajectory_means_B400_pointwise',
                'interpretation':'finite_horizon_MC_only'})
        print(f'completed q={q:.2f} a={a:.2f}',flush=True)
df=pd.DataFrame(records)
out=HERE/'figure4_fivepoint_SE_CI95.csv'
df.to_csv(out,index=False,float_format='%.17g')
np.savez_compressed(HERE/'figure4_trajectory_costs.npz',**raw_costs)
df[['a1','a2','q','a','policy','seed']].assign(npz_key=[f"traj_{i:03d}" for i in range(len(df))]).to_csv(HERE/'figure4_trajectory_cost_index.csv',index=False)
# Cross-check existing q=.92,.96,.98 original values, preserved random streams.
original=pd.read_csv(HERE.parent/'results/controller_performance.csv')
merged=original.merge(df,how='inner',left_on=['a1','a2','q','a','controller'],right_on=['a1','a2','q','a','policy'])
assert len(merged)==43,(len(merged), merged[['q','a','policy']].head())
np.testing.assert_allclose(merged.mean_J_MC.to_numpy(),merged.mean_cost.to_numpy(),rtol=5e-12,atol=1e-12)
print('Regression PASS: all 43 original benchmark costs exactly reproduced to rounding.')

# Two displays: readable journal figure (with explicit offscale outlier)
# and full-scale comparison showing all values.
def make_plot(mode):
    fig,axs=plt.subplots(1,3,figsize=(16.5,4.7),layout='constrained')
    for ax,(lo,hi,q) in zip(axs,SCENARIOS):
        sub=df[(np.isclose(df.a1,lo)) & (np.isclose(df.a2,hi))]
        sub=sub[np.isclose(sub.q,q)]
        for policy in POLICIES:
            ss=sub[sub.policy==policy].sort_values('a')
            if ss.empty:continue
            if mode=='readable' and policy=='two_received' and np.isclose(q,.98):
                ss=ss[ss.a<hi-1e-8]   # state the excluded outlier explicitly
            line,=ax.plot(ss.a,ss.mean_cost,marker='o',markersize=3.5,lw=1.3,label=LABEL[policy])
            ax.fill_between(ss.a,ss.bootstrap_CI95_low,ss.bootstrap_CI95_high,
                        color=line.get_color(),alpha=.15,lw=0)
        ax.set(title=rf'$a\in[{lo:g},{hi:g}],\ q={q:.2f}$', xlabel=r'Plant pole $a$')
        ax.set_xticks(np.linspace(lo,hi,5));ax.set_yscale('log');ax.grid(alpha=.18)
        ax.set_ylabel(r'$J_T$ (finite horizon, $T=300$)')
        ax.tick_params(axis='x',labelsize=8)
        if mode=='readable' and np.isclose(q,.98):
            ax.set_ylim(.008,120)
            extraordinary=sub[(sub.policy=='two_received')&(np.isclose(sub.a,hi))].iloc[0]
            ax.scatter([hi],[100],marker='^',s=36,color='tab:purple',clip_on=False)
            # The extreme point is labelled in the figure caption and data table;
            # do not put an illegible annotation over the legend.
            # The caption gives the exact excluded value; the marker flags it visually.
        ax.legend(fontsize=6.9,loc='upper center',ncol=1,framealpha=.85)
    fig.suptitle('Paired finite-horizon costs; shaded pointwise 95% bootstrap CIs',fontsize=11)
    name=f'fig03_finite_horizon_mc_fivepoint_{mode}'
    fig.savefig(HERE/f'{name}.png',dpi=220,bbox_inches='tight')
    fig.savefig(HERE/f'{name}.pdf',bbox_inches='tight')
    plt.close(fig)
make_plot('readable');make_plot('fullscale')
metadata={'seed':SEED,'N':N,'T':T,'B':B,'process_noise_variance':.01,'initial_state':'N(0,1)',
'control_weight':.05,'q':[.92,.96,.98], 'plant_grid':'a1+(a2-a1)*j/4 for j=0,1,2,3,4',
'stationary_q098_endpoints':'omitted from display; fixed-plant MS unstable (rho>1)',
'stationary_endpoints_rho_q098':[stationary_rho_exact(a,.98,1.5,4.5) for a in [1.5,4.5]],
'standard_error':'sample_SD/sqrt(1200) across independent trajectory costs J_T',
'figure_band':'pointwise percentile 95% bootstrap CI of mean cost, 400 trajectory resamples per controller and plant',
'paired_within_plant':True,'seed_offset_by_index':seed_offset,
'preexisting_three_pole_values_regression':'PASS 43 comparisons',
'warning':'Mean/SE and bootstrap intervals are finite-horizon sampling uncertainty, NOT mean-square-stability certificates; non-stable controllers can have severe rare-event heavy tails; intervals may be unreliable there.'}
(HERE/'figure4_fivepoint_metadata.json').write_text(json.dumps(metadata,indent=2))
print('Total rows:',len(df),'Outputs saved:',out)
