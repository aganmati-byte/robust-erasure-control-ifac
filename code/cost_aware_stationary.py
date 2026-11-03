"""Reproducible cost-aware stationary last-reception control and exact polynomial certificates.
Controller h_0=L and h_j=B*c**(j-1) (j>=1), c fixed at interval midpoint.
Stability certificates use rational coefficients and exact Bernstein positivity on interval cells.
"""
import json
from fractions import Fraction as F
import numpy as np
from scipy.optimize import minimize
import sympy as s

def rho_geom(a,q,L,B,c):
    r=1-q
    if r*a*a>=1 or r*c*c>=1:return float('inf')
    if abs(a-c)<1e-4:
        # The infinite-series expression has severe cancellation near the pole coincidence;
        # sum its positive terms with an explicit convergence criterion instead.
        alpha=1.;cum=0.;p=q
        for g in range(1,2000):
            gain=L if g==1 else B*c**(g-2)
            alpha=a*alpha-gain
            v=p*alpha*alpha;cum+=v
            if g>80 and v<1e-24:break
            p*=r
        return float(cum)
    v=B/(a-c);u=(a-L)-v
    return float(q*(u*u/(1-r*a*a)+2*u*v/(1-r*a*c)+v*v/(1-r*c*c)))

def stationary_cost_geom(a,q,L,B,c,sigma2=.01,weight=.05):
    r=1-q
    rho=rho_geom(a,q,L,B,c)
    if rho>=1 or not np.isfinite(rho):return {'rho':rho,'Ex2':float('inf'),'Eu2':float('inf'),'J':float('inf')}
    V=sigma2/(1-r*a*a)/(1-rho)
    Ex=(q+r*rho)*V+r*sigma2/(1-r*a*a)
    Eu=q*V*(L*L+r*B*B/(1-r*c*c))
    return {'rho':rho,'Ex2':float(Ex),'Eu2':float(Eu),'J':float(Ex+weight*Eu)}

def _rat(v):return s.Rational(F.from_float(v).numerator,F.from_float(v).denominator) if isinstance(v,(float,np.floating)) else s.Rational(str(v))

def rho_numerator_rational(lo,hi,q,L,B,c):
    """Returns H(a)=(1-rho(a))*D_aa*D_ac*D_cc, exact polynomial degree <=3.
    All denominators are positive when r*hi**2<1 and c in [lo,hi].
    """
    a=s.Symbol('a')
    q,L,B,c=map(_rat,(q,L,B,c));r=1-q
    A=1-r*a*a;C=1-r*a*c;D=1-r*c*c
    N=(a-L)*(a-c)-B
    raw=s.Poly(s.expand((a-c)**2*A*C*D-q*(N*N*C*D+2*N*B*A*D+B*B*A*C)),a)
    div,rem=s.div(raw,s.Poly((a-c)**2,a))
    assert rem.is_zero and div.degree()<=3
    return div

def bernstein_certificate(lo,hi,q,L,B,c,segments=64):
    """Exact rational interval certificate for max rho<1, via Bernstein positivity.
    Certifies *strict* whole-continuum inequality if every Bernstein coeff >0.
    """
    lo,hi=map(lambda x:s.Rational(str(x)),(lo,hi))
    ratq=s.Rational(str(q));r=1-ratq
    if not (r*hi*hi<1):return {'certified':False,'reason':'dropout necessary condition'}
    poly=rho_numerator_rational(lo,hi,q,L,B,c);deg=3
    a=s.Symbol('a');t=s.Symbol('t');cells=[];mini=None
    for j in range(segments):
        left=lo+(hi-lo)*j/segments;right=lo+(hi-lo)*(j+1)/segments
        f=s.Poly(poly.as_expr().subs(a,left+(right-left)*t),t)
        b=[f.nth(k) for k in range(deg+1)]
        bern=[s.factor(sum(b[k]*s.binomial(i,k)/s.binomial(deg,k) for k in range(i+1))) for i in range(deg+1)]
        loc_min=min(bern);mini=loc_min if mini is None or loc_min<mini else mini
        if loc_min<=0:cells.append(j)
    return {'certified':not bool(cells),'segments':segments,'bad_segments':cells,
        'min_Bernstein_coefficient_decimal':str(s.N(mini,15)),
        'stability_criterion':'H(a)>0 on each cell, D_aa D_ac D_cc>0, H=(1-rho)*D_aa D_ac D_cc',
        'exact_rational':True, 'gain_hex':{'L':float(L).hex(),'B':float(B).hex(),'c':float(c).hex()},
        'a1':str(lo),'a2':str(hi),'q':str(ratq)}

def optimize_cost_aware(lo=1.2,hi=3.,q=.96,n_grid=31,weight=.05):
    agrid=np.linspace(lo,hi,n_grid);c=(lo+hi)/2;X=(hi-lo)**2/4
    z0=np.array([c,X,0.7])
    # Epigraph formulation: max over grid of exact asymptotic J, with per-plant strict rho constraints.
    # Interpret as grid-based design; certify continuum separately with Bernstein method.
    def costs(z):return np.array([stationary_cost_geom(a,q,z[0],z[1],c,weight=weight)['J'] for a in agrid])
    def constraint(z):
        ro=np.array([rho_geom(a,q,z[0],z[1],c) for a in agrid]);return np.r_[0.999-ro,z[2]-costs(z)]
    res=minimize(lambda z:z[2],z0,method='SLSQP',constraints=[{'type':'ineq','fun':constraint}],
          bounds=[(c-.15,c+.15),(X-.5,X+.5),(0.00001,100.)],options={'maxiter':260,'ftol':2e-10})
    candidates=[('optimal_slsqp_grid',res.x)] if res.success else []
    candidates.append(('stationary_threshold_optimal',z0))
    for label,z in candidates:
        cert=bernstein_certificate(lo,hi,q,float(z[0]),float(z[1]),c)
        if cert['certified']:
            xx=np.linspace(lo,hi,1001)
            ys=[stationary_cost_geom(a,q,float(z[0]),float(z[1]),c,weight=weight)['J'] for a in xx]
            return {'design':label,'q':q,'lo':lo,'hi':hi,'c':c,'L':float(z[0]),'B':float(z[1]),
                'optimization_grid_points':n_grid,'objective_grid_bound':float(z[2]),
                'max_dense_grid_J':float(np.max(ys)), 'max_dense_grid_a':float(xx[np.argmax(ys)]),
                'max_grid_rho':float(max(rho_geom(a,q,float(z[0]),float(z[1]),c) for a in xx)),
                'benchmark_J_hstar':float(max(stationary_cost_geom(a,q,c,X,c,weight=weight)['J'] for a in xx)),
                'certificate':cert,'statement':'SLSQP grid local result, not globally cost-optimal; continuum strict stability separately certified.'}
    raise RuntimeError('No Bernstein-certified candidate found, despite analytically feasible hstar')

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--output',default='cost_aware_stationary.json');args=p.parse_args()
    result=optimize_cost_aware()
    from pathlib import Path
    Path(args.output).write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='certificate'},indent=2));print(result['certificate'])
    assert result['certificate']['certified']
