"""Run: python run_benchmark.py --profile standard. CPU, deterministic, exported results."""
import argparse,sys,os,json,platform,hashlib,time
from pathlib import Path
import numpy as np,pandas as pd,scipy,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from robust_erasure import *

def main(profile='standard'):
 start=time.time();root=Path(__file__).resolve().parent;out=root/'results';out.mkdir(exist_ok=True)
 cfg={'quick':(200,100,301,10001),'standard':(1200,300,1201,40001),'extended':(5000,600,2001,120001)}[profile]
 N,T,gridN,certN=cfg;seed=20261020
 scenarios=[(1.2,3.,.92,[1.2,2.1,3.]),(1.2,3.,.96,[1.2,2.1,3.]),(1.5,4.5,.98,[1.5,3.0,4.5])]
 policies=['oracle','stationary','nominal','cheb','optimized','two_received','open']
 summary=[];stab=[];poly=[];regr=[];optmap={};paired_cis=[]
 for lo,hi,q,poles in scenarios:
  mnodes=cheb_nodes(lo,hi,6)
  opt,coeff,lp=weighted_periodic_lp(lo,hi,q,m=6,n_grid=gridN)
  cert=certificate_weighted_poly(lo,hi,q,coeff,n_grid=certN)
  optmap[(lo,hi,q)]=opt
  poly.append({'a1':lo,'a2':hi,'q':q,'period':6,'lp_grid_peak':lp,'continuum_bound_float':cert['bound'],'continuum_bound_positive_margin':1-cert['bound'],'continuum_bound_status':cert['status'],'continuum_contraction':cert['certifies_strict_contraction'],'nodes':json.dumps(opt.tolist()),'coeff_chebyshev':json.dumps(coeff.tolist())})
  agrid=np.linspace(lo,hi,801)
  for label,nodes in [('cheb',mnodes),('optimized',opt)]:
   logs=np.array([periodic_log_growth(a,q,nodes) for a in agrid]);finite=logs[np.isfinite(logs)]
   stab.append({'a1':lo,'a2':hi,'q':q,'controller':label,'continuum_status':('analytic_exact_cheb' if label=='cheb' else 'derivative_bound_float') if (q>q_cheb(lo,hi,6) if label=='cheb' else cert['certifies_strict_contraction']) else 'no_continuum_stability_certificate','max_grid_loggrowth':max(finite),'continuum_stable':(q>q_cheb(lo,hi,6) if label=='cheb' else cert['certifies_strict_contraction'])})
  stab.append({'a1':lo,'a2':hi,'q':q,'controller':'stationary','continuum_status':'analytic_stationary_theorem' if q>q_stationary(lo,hi) else 'analytically_not_robust','max_grid_loggrowth':np.nan,'continuum_stable':q>q_stationary(lo,hi)})
  for idx,a in enumerate(poles):
   rg=np.random.default_rng(seed+idx+int(lo*100)+int(q*1000));packets=rg.random((N,T))<q;noise=rg.normal(0,.1,(N,T));x0=rg.normal(0,1,N)
   by_pol={}
   for pol in policies:
    nodes=opt if pol=='optimized' else mnodes if pol=='cheb' else None
    outp=sim_one(a,q,lo,hi,packets,noise,x0,pol,nodes=nodes)
    by_pol[pol]=outp['J']
    if pol=='stationary':r=stationary_rho_exact(a,q,lo,hi);e=exact_cost_stationary(a,q,lo,hi)
    elif pol=='two_received':r=two_received_spectral_radius(a,q,lo,hi);e={'status':'no_exact_cost_implemented','Ex2':np.nan,'Eu2':np.nan,'J':np.nan}
    else:
     ns=nodes if nodes is not None else [a] if pol=='oracle' else [(lo+hi)/2] if pol=='nominal' else [0]
     lg=periodic_log_growth(a,q,ns)
     r=float(np.exp(lg)) if lg>-np.inf and np.isfinite(lg) else (0. if lg==-np.inf else np.inf)
     e=exact_cost_impulse(a,q,ns) if pol!='open' else {'status':'unstable_or_boundary','Ex2':np.inf,'Eu2':0.,'J':np.inf}
    if pol=='two_received':status='exact_spectral_fixed_plant_stable' if r<1-1e-10 else 'exact_spectral_fixed_plant_unstable'
    else:status='exact_fixed_plant_stable' if r<1-1e-10 else 'exact_fixed_plant_unstable'
    vals=outp['J'];summary.append({'a1':lo,'a2':hi,'a':a,'q':q,'controller':pol,'per_plant_status':status,'second_moment_multiplier':r,'finite_horizon_T':T,'N':N,'mean_J_MC':np.mean(vals),'median_J_MC':np.median(vals),'p95_J_MC':np.quantile(vals,.95) if np.all(np.isfinite(vals)) else np.nan,'MC_finite_trajectory_fraction':np.isfinite(vals).mean(),'mean_x2_MC':np.mean(outp['Jx']),'mean_u2_MC':np.mean(outp['Ju']),'mean_peak_x_MC':np.mean(outp['maxx']),'mean_peak_u_MC':np.mean(outp['maxu']),'exact_E_x2':e.get('Ex2',np.nan),'exact_E_u2':e.get('Eu2',np.nan),'exact_J':e.get('J',np.nan),'exact_cost_status':e['status']})
   # Paired bootstrap of *finite-horizon* differences: descriptive only.
   for A,B in [('optimized','cheb'),('optimized','stationary'),('stationary','cheb')]:
    x=by_pol[A];y=by_pol[B]
    if np.isfinite(x).all() and np.isfinite(y).all():
     D=x-y;rb=np.random.default_rng(seed+404+idx)
     means=np.array([D[rb.integers(0,N,N)].mean() for _ in range(300)])
     paired_cis.append({'a1':lo,'a2':hi,'a':a,'q':q,'A':A,'B':B,'N':N,'horizon':T,'mean_JA_minus_JB':D.mean(),'CI95_percentile_low':np.quantile(means,.025),'CI95_percentile_high':np.quantile(means,.975),'interpretation':'trajectory_bootstrap_finite_horizon_only_not_asymptotic'})
 result=pd.DataFrame(summary);result.to_csv(out/'controller_performance.csv',index=False)
 pd.DataFrame(paired_cis).to_csv(out/'paired_bootstrap_finite_horizon.csv',index=False)
 pd.DataFrame(stab).to_csv(out/'continuum_stability.csv',index=False)
 pd.DataFrame(poly).to_csv(out/'weighted_periodic_designs.csv',index=False)
 # Exact-rational proof of continuum bound, independently from any floating test.
 # IEEE-754 stored gains are interpreted as exact dyadic rationals.
 from verify_weighted_exact import from_saved_csv
 rational_certs=from_saved_csv(out)
 assert len(rational_certs)==len(poly) and all(x['bound_less_than_one_exact'] for x in rational_certs)
 (out/'weighted_periodic_exact_certificates.json').write_text(json.dumps(rational_certs,indent=2))
 for row,cert in zip(poly,rational_certs):
  row['continuum_bound_exact_rational']=float(cert['continuum_bound_decimal'])
  row['continuum_certified_exact']=cert['bound_less_than_one_exact']
 pd.DataFrame(poly).to_csv(out/'weighted_periodic_designs.csv',index=False)
 # Exact stationary versus finite-horizon MC: diagnostics never used as same estimand
 thresholds=[]
 for lo,hi,q,_ in scenarios:
  thresholds.append({'a1':lo,'a2':hi,'q':q,'q_known_a2':q_known(hi),'q_stationary':q_stationary(lo,hi),'q_cheb_m6':q_cheb(lo,hi,6),'q_cheb_m32':q_cheb(lo,hi,32),'stationary_robust':q>q_stationary(lo,hi),'cheb_m6_robust':q>q_cheb(lo,hi,6)})
 pd.DataFrame(thresholds).to_csv(out/'thresholds.csv',index=False)
 # figure 1 stability / phase diagram
 fig,axs=plt.subplots(1,2,figsize=(12,4.3))
 for ax,(lo,hi,q,_) in zip(axs,scenarios[::2]):
  aa=np.linspace(lo,hi,1201)
  nc=cheb_nodes(lo,hi,6);nw=optmap[(lo,hi,q)]
  for name,nodes in [('Chebyshev m=6',nc),('Weighted LP m=6',nw)]:
   xx=np.array([periodic_log_growth(float(a),q,nodes) for a in aa]);xx=np.where(np.isneginf(xx),np.nan,xx)
   ax.plot(aa,xx,label=name)
  ax.axhline(0,ls='--',color='black',lw=1);ax.grid(alpha=.2);ax.set(title=f'[{lo}, {hi}] q={q}; fixed-policy exact second moments',xlabel='a',ylabel='log geometric MS multiplier / arrival');ax.legend(fontsize=8)
 fig.tight_layout();fig.savefig(out/'fig01_periodic_stability.png',dpi=180);plt.close(fig)
 # figure 2 exact stationary and periodic cost on stable fixed plants; mask infinity
 fig,axs=plt.subplots(1,3,figsize=(13,4.2),sharey=False)
 for ax,(lo,hi,q,_) in zip(axs,scenarios):
  a=np.linspace(lo,hi,121)
  for label in ['oracle','stationary','cheb','optimized']:
   nodes=optmap[(lo,hi,q)] if label=='optimized' else cheb_nodes(lo,hi,6) if label=='cheb' else None
   js=[]
   for aa in a:
    e=exact_cost_stationary(aa,q,lo,hi) if label=='stationary' else exact_cost_impulse(aa,q,[aa] if label=='oracle' else nodes)
    js.append(e['J'] if np.isfinite(e['J']) and e['J']<1e9 else np.nan)
   ax.plot(a,js,label=label)
  ax.grid(alpha=.2);ax.set(title=f'[{lo},{hi}] q={q}',xlabel='a',ylabel='exact asymptotic $E[x^2+.05u^2]$');ax.set_yscale('log');ax.legend(fontsize=8)
 fig.tight_layout();fig.savefig(out/'fig02_exact_asymptotic_costs.png',dpi=180);plt.close(fig)
 # fig 3 paired mc finite horizon
 fig,axs=plt.subplots(1,3,figsize=(13,4.2))
 for ax,(lo,hi,q,_) in zip(axs,scenarios):
  for p in ['oracle','stationary','cheb','optimized','two_received']:
   ss=result[(result.a1==lo)&(result.a2==hi)&(result.q==q)&(result.controller==p)].sort_values('a')
   ax.plot(ss.a,ss.mean_J_MC,'-o',label=p)
  ax.set(yscale='log',xlabel='a',ylabel='finite-horizon sampled cost',title=f'[{lo},{hi}] q={q}');ax.grid(alpha=.2);ax.legend(fontsize=7)
 fig.tight_layout();fig.savefig(out/'fig03_finite_horizon_mc.png',dpi=180);plt.close(fig)
 # fig 4 thresholds
 fig,ax=plt.subplots(figsize=(7,4.5))
 hi=3.;aa=np.linspace(1.01,hi-.01,140)
 ax.plot(hi-aa,[q_stationary(lo,hi) for lo in aa],label='stationary exact')
 for m in [4,6,32]:ax.plot(hi-aa,[q_cheb(lo,hi,m) for lo in aa],label=f'ordinary Cheb m={m}')
 ax.axhline(q_known(hi),ls=':',c='black',label='known-plant necessary')
 ax.axhline(1.,ls='--',c='grey');ax.set(xlabel='interval width',ylabel='critical arrival probability',ylim=(.85,1.15));ax.legend(fontsize=8);ax.grid(alpha=.2)
 fig.tight_layout();fig.savefig(out/'fig04_analytic_thresholds.png',dpi=180);plt.close(fig)
 # fixed-plant spectral memory stability table
 rows=[]
 for lo,hi,q,poles in scenarios:
  for a in np.linspace(lo,hi,501):rows.append({'a1':lo,'a2':hi,'q':q,'a':a,'two_received_spectral_radius':two_received_spectral_radius(a,q,lo,hi)})
 pd.DataFrame(rows).to_csv(out/'two_received_finite_grid.csv',index=False)
 manifest={'profile':profile,'seed':seed,'n_trajectories':N,'horizon':T,'optimization_grid':gridN,'certificate_grid':certN,'noise_std':.1,'state_init_std':1.,'quadratic_weight_u':.05,'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'duration_sec':round(time.time()-start,2),'notes':'No claim that float derivative bound is formally verified interval arithmetic. No SOTA status.'}
 (out/'run_metadata.json').write_text(json.dumps(manifest,indent=2))
 print('COMPLETE',manifest['duration_sec'],'sec',result.shape,'tests rerun separately with pytest')
 print(pd.DataFrame(poly)[['a1','a2','q','lp_grid_peak','continuum_bound_float','continuum_contraction']].to_string(index=False))
 print(result[['a1','a2','a','q','controller','per_plant_status','exact_J','mean_J_MC']].to_string(index=False,max_rows=23))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--profile',choices=['quick','standard','extended'],default='standard');main(p.parse_args().profile)
