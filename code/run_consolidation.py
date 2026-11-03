"""CPU deterministic supplementary experiment: threshold-optimal vs cost-aware stationary.
All identical noise, packet and initial-state arrays; exact expected costs and finite horizon MC separate.
"""
from pathlib import Path
import sys,json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cost_aware_stationary import stationary_cost_geom, optimize_cost_aware, bernstein_certificate, rho_geom
from export_full_bernstein_witness import export_witness

ROOT=Path(__file__).resolve().parents[1]

def simulator(a,q,packet,noise,xinit,L,B,c):
    n,T=packet.shape;xx=xinit.copy();last=np.zeros(n);age=np.zeros(n,dtype=int);received=np.zeros(n,dtype=bool)
    ss=np.zeros(n);uu=np.zeros(n);maxx=np.abs(xx);maxu=np.zeros(n)
    for k in range(T):
        got=packet[:,k];last=np.where(got,xx,last);received|=got;age=np.where(got,0,age+1)
        gains=np.where(age==0,L, B*np.power(c,np.minimum(np.maximum(age-1,0),500)))
        u=np.where(received,-gains*last,0.)
        ss+=xx**2;uu+=u**2;maxx=np.maximum(maxx,abs(xx));maxu=np.maximum(maxu,abs(u))
        xx=a*xx+u+noise[:,k]
    return {'J_samples':(ss+0.05*uu)/T,'Ex2_samples':ss/T,'Eu2_samples':uu/T,'peakx_samples':maxx,'peaku_samples':maxu}

def main(profile='standard'):
    out=ROOT/'results';out.mkdir(exist_ok=True)
    opt=optimize_cost_aware()
    (out/'cost_aware_stationary.json').write_text(json.dumps(opt,indent=2))
    export_witness(opt['lo'],opt['hi'],opt['q'],opt['L'],opt['B'],opt['c'],output=out/'cost_aware_bernstein_full_rational_certificate.json')
    lo,hi,q=opt['lo'],opt['hi'],opt['q'];c=(lo+hi)/2;X=(hi-lo)**2/4
    xx=np.linspace(lo,hi,101)
    rows=[]
    for a in xx:
        for label,L,B in [('threshold_optimal',c,X),('cost_aware_locally_optimized',opt['L'],opt['B'])]:
            R=stationary_cost_geom(a,q,L,B,c)
            rows.append({'a':float(a),'q':q,'controller':label,**R})
    with (out/'cost_aware_exact_curves.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    fig,ax=plt.subplots(figsize=(7.7,4.8))
    for label in ['threshold_optimal','cost_aware_locally_optimized']:
        r=[d for d in rows if d['controller']==label]
        ax.plot([i['a'] for i in r],[i['J'] for i in r],label=label.replace('_',' '))
    ax.set(xlabel='Uncertain fixed pole a',ylabel='Exact limiting physical-time quadratic cost',title='Feasibility-optimal vs cost-aware stationary policy')
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'fig05_cost_aware_stationary.png',dpi=170);plt.close(fig)
    n,T={'quick':(500,250),'standard':(3500,500)}[profile]
    rng=np.random.default_rng(20261020+44)
    packet=rng.random((n,T))<q;noise=rng.normal(0,.1,(n,T));xinit=rng.normal(0,1,n)
    simulation=[]
    for a in [lo,c,hi]:
        for label,L,B in [('threshold_optimal',c,X),('cost_aware_locally_optimized',opt['L'],opt['B'])]:
            s=simulator(a,q,packet,noise,xinit,L,B,c);J=s['J_samples']
            analytic=stationary_cost_geom(a,q,L,B,c)
            simulation.append({'a':a,'q':q,'controller':label,'N':n,'T':T,'finite_horizon_mean_J':float(np.mean(J)),
               'finite_horizon_median_J':float(np.median(J)), 'finite_horizon_peak_x_mean':float(np.mean(s['peakx_samples'])),
               'finite_horizon_peak_u_mean':float(np.mean(s['peaku_samples'])),
               'asymptotic_exact_J':analytic['J'],'exact_asymptotic_Ex2':analytic['Ex2'],'exact_asymptotic_Eu2':analytic['Eu2'],
               'interpretation':'finite horizon MC mean is NOT estimator of steady-state in early/transient observation windows'})
    with (out/'cost_aware_paired_simulation.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=simulation[0].keys());writer.writeheader();writer.writerows(simulation)
    stats={'profile':profile,'seed':20261064,'n_trajectories':n,'horizon':T,'noise_std':.1,
           'initial_std':1.,'q':q,'source_v3':'untouched','stationary_gain_convention':'h0=L; h_j=B*c^(j-1) j>=1',
           'author_note':'SLSQP local finite-grid cost optimization, exact Bernstein continuum stability, exact cost formulas; not global cost optimality'}
    (out/'consolidation_metadata.json').write_text(json.dumps(stats,indent=2))
    print('Certified stable:',opt['certificate']['certified'],'Worst J costaware',opt['max_dense_grid_J'], 'original',opt['benchmark_J_hstar'])
    print('Simulation rows:',simulation)
    return opt
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--profile',choices=['quick','standard'],default='standard');args=p.parse_args();main(args.profile)
