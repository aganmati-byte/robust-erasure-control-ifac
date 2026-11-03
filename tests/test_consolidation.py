import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import numpy as np
import sympy as s
from cost_aware_stationary import rho_geom,stationary_cost_geom,bernstein_certificate,optimize_cost_aware,rho_numerator_rational
from robust_erasure import stationary_rho_exact,exact_cost_stationary,periodic_log_growth

def test_general_matches_v3_multiplier():
 for a in [1.2,1.5,2.,2.1,2.2,2.7,3.]:
  lo,hi,q=1.2,3.,.96;c=(lo+hi)/2;X=(hi-lo)**2/4
  assert np.isclose(rho_geom(a,q,c,X,c),stationary_rho_exact(a,q,lo,hi),atol=1e-10)

def test_general_matches_v3_cost():
 for a in [1.2,2.1,2.8,3.]:
  lo,hi,q=1.2,3.,.96;c=(lo+hi)/2;X=(hi-lo)**2/4
  one=stationary_cost_geom(a,q,c,X,c);two=exact_cost_stationary(a,q,lo,hi)
  for k in ['Ex2','Eu2','J']:assert np.isclose(one[k],two[k],atol=1e-10)

def test_renewal_vs_direct_series():
 a,q,L,B,c=2.8,.96,2.103,.815,2.1
 v=1.;s1=0.;p=q
 for g in range(1,300):
  v=a*v-(L if g==1 else B*c**(g-2))
  s1+=p*v*v;p*=1-q
 assert np.isclose(rho_geom(a,q,L,B,c),s1,atol=1e-11)

def test_exact_rational_identity():
 a=s.Symbol('a');q,L,B,c=.96,2.104692202603161,.8143494389928447,2.1
 H=rho_numerator_rational(1.2,3.,q,L,B,c).as_expr()
 for aa in [s.Rational(6,5),s.Rational(21,10),s.Integer(3)]:
  rr=rho_geom(float(aa),q,L,B,c)
  den=(1-s.Rational(4,100)*aa**2)*(1-s.Rational(4,100)*aa*s.Rational(c))*(1-s.Rational(4,100)*s.Rational(c)**2)
  assert abs(float(H.subs(a,aa)/den)-(1-rr))<1e-8

def test_bernstein_certificate():
 opt=optimize_cost_aware()
 assert opt['certificate']['certified']
 assert opt['max_dense_grid_J']<opt['benchmark_J_hstar']

def test_uncertified_large_wrong_gain():
 cert=bernstein_certificate(1.2,3.,.96,5.,0.,2.1,segments=16)
 assert cert['certified'] is False

def test_zero_periodic_log():
 assert periodic_log_growth(3.,.96,[1.2,3.])==-np.inf

def test_full_rational_bernstein_witness():
    from export_full_bernstein_witness import export_witness
    import sympy as s
    witness=export_witness(1.2,3.0,.96,2.104692202603161,.8143494389928447,2.1,segments=64)
    assert witness['certified_all_strictly_positive_rational_coefficients']
    assert len(witness['cells'])==64
    assert all(s.Rational(x)>0 for cell in witness['cells'] for x in cell['bernstein_coefficients_exact'])
