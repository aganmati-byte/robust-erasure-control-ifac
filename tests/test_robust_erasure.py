import sys, pathlib
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from robust_erasure import *

def test_known_threshold():
 assert abs(q_known(3)-8/9)<1e-14
 assert abs(q_stationary(1.2,3)-.9472222222222)<1e-11
 assert abs(q_cheb(1.5,4.5,128)-.9780820811)<1e-9

def test_gain_sequence_dimensions_and_recursion():
 h=stationary_gains(1.2,3,5)
 assert h.shape==(5,)
 assert np.allclose(h[:3],[2.1,.81,1.701])
 assert np.allclose(m_of_h(3,h)[:2],[.9,1.89])

def test_logzero_is_negative_infinity():
 assert periodic_log_growth(2.1,.96,[2.1,1.2])==-np.inf

def test_exact_stationary_rho_vs_series():
 for a in [1.2,1.5,2.1,2.1000001,2.5,3.]:
  for q in [.96,.99]:
   exact=stationary_rho_exact(a,q,1.2,3.)
   approx=stationary_rho_series(a,q,1.2,3.,n=400)
   assert abs(exact-approx)<max(1,abs(exact))*1e-7,(a,q,exact,approx)

def test_exact_periodic_vs_renewal():
 for a in [1.5,2.4,3.]:
  nodes=cheb_nodes(1.2,3,4);v=periodic_reception_cycle(a,.96,nodes)
  if v is not None:
   den=1-(1-.96)*a*a; rho=.96*(a-nodes)**2/den;b=.01/den
   assert np.allclose(np.roll(v,-1),rho*v+b,rtol=1e-10,atol=1e-10)

def test_two_received_moments_regression():
 assert abs(two_received_spectral_radius(3.,.92,1.2,3.)-2.1864445227)<1e-7
 assert abs(two_received_spectral_radius(3.,.96,1.2,3.)-.7390078497)<1e-7

def test_weighted_poly_real_roots_and_actual_improvement():
 lo,hi,q=1.2,3.,.94
 nodes,c,_=weighted_periodic_lp(lo,hi,q,m=6,n_grid=801)
 assert len(nodes)==6
 assert np.all((nodes>lo)&(nodes<hi))
 a=np.linspace(lo,hi,10001)
 opt=max(periodic_log_growth(z,q,nodes) for z in a)
 old=max(periodic_log_growth(z,q,cheb_nodes(lo,hi,6)) for z in a)
 assert opt<=old+1e-4,(opt,old)
 result=certificate_weighted_poly(lo,hi,q,c,n_grid=30001)
 assert result['grid_max']<=result['bound']

def test_cost_matches_known_deadbeat_identity():
 a=3.;q=.96;sigma=.01;den=1-(1-q)*a*a
 cost=exact_cost_impulse(a,q,[a],noise_var=sigma)
 assert np.isclose(cost['Ex2'],sigma/den)
 assert np.isclose(cost['Eu2'],q*a*a*sigma/den)

def test_no_finite_cost_when_unstable():
 assert exact_cost_impulse(3,.92,[2.1])['status']=='unstable_or_boundary'
 assert exact_cost_stationary(3,.92,1.2,3.)['status']=='unstable_or_boundary'

def test_paired_identical_noise():
 rg=np.random.default_rng(1);packets=rg.random((8,20))<.95;noise=rg.normal(0,.1,(8,20));x0=rg.normal(size=8)
 v=sim_one(3,.95,1.2,3,packets,noise,x0,'oracle');w=sim_one(3,.95,1.2,3,packets,noise,x0,'oracle')
 assert np.array_equal(v['J'],w['J'])

def test_exact_rational_endpoint_lyapunov_certificate():
 import sympy as S
 from verify_memory3 import lyapunov_certificate
 h=list(map(S.Rational,['1977/1000','143/125','741/250']))
 k=list(map(S.Rational,['21/100','179/1000','9/200']))
 for pole in ['6/5','3']:
  proof=lyapunov_certificate(S.Rational(pole),S.Rational('9/10'),h,k)
  assert proof['P11_positive'] and proof['detP_positive']


def test_exact_rational_continuum_certificate():
 from verify_weighted_exact import exact_continuum_certificate
 for lo,hi,q in [(1.2,3.,.92),(1.5,4.5,.98)]:
  nodes,_,_=weighted_periodic_lp(lo,hi,q,m=6,n_grid=1201)
  report=exact_continuum_certificate(lo,hi,q,nodes,n_grid=1201)
  assert report['exact_rational'] and report['bound_less_than_one_exact']
  numerator,denominator=map(int,report['bound_fraction'].split('/'))
  assert numerator < denominator
