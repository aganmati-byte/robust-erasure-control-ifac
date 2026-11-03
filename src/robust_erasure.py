"""Scalar erasure-channel benchmark: analytic renewal moments, validated numeric tests.
All 'certified' continuum margins below employ an explicit global derivative bound,
not merely a sampled maximum. Theorem statements require separate mathematical proofs.
"""
import numpy as np
from scipy.optimize import linprog
from numpy.polynomial import chebyshev as C

def q_known(a): return 1.0-1.0/float(a)**2
def q_stationary(lo,hi): return max(q_known(hi),((lo+hi)**2-4)/(4*lo*hi))
def cheb_nodes(lo,hi,m):
    j=np.arange(m);return (lo+hi)/2+(hi-lo)/2*np.cos((2*j+1)*np.pi/(2*m))
def B_cheb(lo,hi,m): return (hi-lo)**2/16*2**(2/m)
def q_cheb(lo,hi,m):
    B=B_cheb(lo,hi,m)
    return (hi**2-1)/(hi**2-B) if B<1 else float('nan')
def stationary_gains(lo,hi,count):
    assert count>=1
    z=(lo+hi)/2;X=(hi-lo)**2/4
    return np.r_[z,X*z**np.arange(count-1)]
def m_of_h(a,h):
    out=[];v=1.
    for gain in h: v=a*v-gain;out.append(v)
    return np.array(out)
def stationary_rho_exact(a,q,lo,hi):
    r=1-q;mid=(lo+hi)/2;X=(hi-lo)**2/4;t=a-mid
    if r*a*a>=1 or r*mid*mid>=1: return np.inf
    if abs(t)<1e-6*max(mid,1):
        # Avoid cancellation near midpoint via a quadratically convergent expansion or a long, bounded tail.
        if abs(t)<1e-13:
            b=r*mid*mid
            return q*r*X**2*(1+b)/(1-b)**3
        # recurrence summation: geometric decay; ensure tail controlled
        m=1.;v=0.;rho=0.;p=q
        for g in range(1,15000):
            gain=mid if g==1 else X*mid**(g-2)
            m=a*m-gain
            term=p*m*m
            rho+=term
            if g>50 and term<1e-15*max(1,rho) and (r*a*a)<.999:break
            p*=r
        return rho
    A=t-X/t;B=X/t
    return q*(A*A/(1-r*a*a)+2*A*B/(1-r*a*mid)+B*B/(1-r*mid*mid))
def stationary_rho_series(a,q,lo,hi,n=600):
    r=1-q
    if r*a*a>=1 or r*((lo+hi)/2)**2>=1:return np.inf
    h=stationary_gains(lo,hi,n);m=m_of_h(a,h)
    return q*float(np.sum(np.power(r,np.arange(n))*m*m))
def impulse_rho(a,q,c):
    den=1-(1-q)*a*a
    return q*(a-c)**2/den if den>0 else np.inf
def periodic_log_growth(a,q,nodes):
    den=1-(1-q)*a*a
    if den<=0:return np.inf
    dif=np.abs(a-np.asarray(nodes,dtype=float))
    if np.any(dif==0):return -np.inf
    return float(np.log(q/den)+2*np.log(dif).mean())
def periodic_reception_cycle(a,q,nodes,noise_var=.01):
    """Return second moments at reception before applying phase j, or None if unstable."""
    nodes=np.asarray(nodes);den=1-(1-q)*a*a
    if den<=0:return None
    rho=q*(a-nodes)**2/den;b=noise_var/den
    # v_m = product(rho)*v_0 + b*(1+rho[-1]+rho[-1]*rho[-2]+...)
    z=0.;mult=1.
    for t in range(len(nodes)): z=rho[t]*z+b;mult*=rho[t]
    if not mult<1-1e-12:return None
    # z = mult*v_0 + additive ; calculate additive z starting at v_0=0
    v=np.empty(len(nodes));v[0]=z/(1-mult)
    for t in range(1,len(nodes)):v[t]=rho[t-1]*v[t-1]+b
    return v
def exact_cost_impulse(a,q,nodes,noise_var=.01,weight=.05):
    den=1-(1-q)*a*a
    v=periodic_reception_cycle(a,q,nodes,noise_var)
    if v is None:return {'status':'unstable_or_boundary','Ex2':np.inf,'Eu2':np.inf,'J':np.inf,'peak_phase_v':np.inf}
    nodes=np.asarray(nodes);r=1-q
    ex=float(np.mean((q+q*r*(a-nodes)**2/den)*v)+r*noise_var/den)
    eu=float(q*np.mean(nodes**2*v))
    return {'status':'exact_fixed_plant_MS','Ex2':ex,'Eu2':eu,'J':ex+weight*eu,'peak_phase_v':float(v.max())}
def exact_cost_stationary(a,q,lo,hi,noise_var=.01,weight=.05):
    r=1-q;mid=(lo+hi)/2;X=(hi-lo)**2/4;den=1-r*a*a;dn=1-r*mid**2
    rho=stationary_rho_exact(a,q,lo,hi)
    if den<=0 or dn<=0 or not rho<1-1e-12:return {'status':'unstable_or_boundary','Ex2':np.inf,'Eu2':np.inf,'J':np.inf}
    b=noise_var/den;v=b/(1-rho)
    ex=(q+r*rho)*v+r*b
    eu=q*v*(mid**2+X**2*r/dn)
    return {'status':'exact_fixed_plant_MS','Ex2':ex,'Eu2':eu,'J':ex+weight*eu,'arrival_v':v,'rho':rho}
def two_received_matrix(a,q,mid,beta):
    r=1-q
    if r*a*a>=1:return None
    K1=q/(1-r*a);K2=q/(1-r*a*a);D=a-mid
    return np.array([[K2*D**2,-2*K2*D*beta,K2*beta**2],
                     [K1*D,-K1*beta,0.],[1.,0.,0.]])
def two_received_spectral_radius(a,q,lo,hi):
    M=two_received_matrix(a,q,(lo+hi)/2,(hi-lo)**2/16)
    return float(max(abs(np.linalg.eigvals(M)))) if M is not None else np.inf

def weighted_periodic_lp(lo,hi,q,m=6,n_grid=1201):
    """Grid-optimized weighted monic polynomial, no claim of continuum optimality."""
    a=np.linspace(lo,hi,n_grid);x=(2*a-lo-hi)/(hi-lo)
    den=1-(1-q)*a*a
    if np.any(den<=0):raise ValueError('q below known-plant necessary threshold')
    # weighted factor >=0; homogeneous return per-cycle MS factor <=1 iff max weight*|P| <1
    weight=(q/den)**(m/2)
    leading=(hi-lo)**m/2**(2*m-1)
    bas=C.chebvander(x,m)
    A=np.column_stack([weight[:,None]*bas[:,:m],-np.ones(len(a))])
    bb=-leading*weight*bas[:,m]
    Aub=np.vstack([A, np.column_stack([-weight[:,None]*bas[:,:m],-np.ones(len(a))])]);bub=np.r_[bb,-bb]
    result=linprog(np.r_[np.zeros(m),1.],A_ub=Aub,b_ub=bub,bounds=[(None,None)]*m+[(0,None)],method='highs',options={'dual_feasibility_tolerance':1e-8,'primal_feasibility_tolerance':1e-8})
    if not result.success:raise RuntimeError(result.message)
    coeff=np.r_[result.x[:m],leading];roots=C.chebroots(coeff)
    if max(abs(np.imag(roots)))>1e-7:raise RuntimeError('nonreal roots in optimizer')
    nodes=np.sort((lo+hi)/2+(hi-lo)/2*np.real(roots))
    if np.min(nodes)<lo-1e-7 or np.max(nodes)>hi+1e-7:raise RuntimeError('root outside interval')
    return nodes,coeff,float(result.fun)

def certificate_weighted_poly(lo,hi,q,coeff,n_grid=100001):
    """Mathematically conservative continuum upper bound via ||T'_k||<=k^2.
    Uses IEEE float outward slack: not formal interval arithmetic/verified rounding.
    Thus labelled ANALYTIC_DERIVATIVE_BOUND+FLOAT, not computer-assisted proof.
    """
    m=len(coeff)-1;r=1-q
    if r*hi*hi>=1:return {'bound':np.inf,'status':'necessary_dropout_failed'}
    a=np.linspace(lo,hi,n_grid);x=(2*a-lo-hi)/(hi-lo)
    P=C.chebval(x,coeff);den=1-r*a*a;W=(q/den)**(m/2)
    Y=np.abs(W*P); max_grid=float(Y.max())
    Pb=sum(abs(v) for v in coeff);Pd=sum(abs(coeff[k])*k*k for k in range(1,m+1))*2/(hi-lo)
    Wmax=(q/(1-r*hi*hi))**(m/2)
    # |W'| <= Wmax*m*r*hi/(1-r*hi^2).
    Lips=Wmax*(Pd+Pb*m*r*hi/(1-r*hi*hi))
    bound=max_grid+Lips*(hi-lo)/(2*(n_grid-1))
    return {'bound':float(bound),'grid_max':max_grid,'derivative_bound':float(Lips),'status':'derivative_continuum_bound_float','certifies_strict_contraction':bool(bound<1-1e-10)}

def sim_one(a,q,lo,hi,packets,noise,x0,policy,nodes=None,cap=None):
    N,T=packets.shape;x=np.asarray(x0,dtype=float).copy();last=np.zeros(N);prev=np.zeros(N)
    age=np.zeros(N,dtype=int);received=np.zeros(N,dtype=bool);count=np.zeros(N,dtype=int)
    mid=(lo+hi)/2;X=(hi-lo)**2/4;beta=X/4
    jx=np.zeros(N);ju=np.zeros(N);maxx=np.abs(x);maxu=np.zeros(N)
    for k in range(T):
        arrive=packets[:,k];prev=np.where(arrive&received,last,prev);last=np.where(arrive,x,last)
        received |=arrive;age=np.where(arrive,0,age+1);count+=arrive
        if policy=='oracle':u=np.where(arrive,-a*last,0.)
        elif policy=='nominal':u=np.where(arrive,-mid*last,0.)
        elif policy=='open':u=np.zeros(N)
        elif policy=='stationary':
            # Avoid np.power overflow if not received yet or extraordinarily long gap
            gain=np.where(age==0,mid,X*np.power(mid,np.minimum(np.maximum(age-1,0),600)))
            u=np.where(received,-gain*last,0.)
        elif policy in ['cheb','optimized']:
            ix=np.maximum(count-1,0)%len(nodes)
            u=np.where(arrive,-nodes[ix]*last,0.)
        elif policy=='two_received':u=np.where(arrive,-mid*last-beta*prev,0.)
        else:raise ValueError(policy)
        if cap is not None:u=np.clip(u,-cap,cap)
        jx+=x*x;ju+=u*u;maxx=np.maximum(maxx,np.abs(x));maxu=np.maximum(maxu,np.abs(u))
        x=a*x+u+noise[:,k]
    return {'J':(jx+.05*ju)/T,'Jx':jx/T,'Ju':ju/T,'maxx':maxx,'maxu':maxu,'end':x*x}
