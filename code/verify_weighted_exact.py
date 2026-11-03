"""Exact-rational continuous-interval stability certificate for even-period impulse gains.

The actual IEEE-754 controller gains are interpreted as exact dyadic rational
numbers (Fraction.from_float); every inequality is evaluated using Python
arbitrary-precision rational arithmetic. Valid for even m where W(a) is rational.
No interval package, floating-point epsilon or sampled-only claim is required.
This proves a whole-interval stability assertion for the exact real values
represented by the stored binary64 nodes; a numerical hardware execution is
still subject to finite precision and model assumptions.
"""
from fractions import Fraction as F
import json
from pathlib import Path

def _mul(p,q):
 out=[F(0)]*(len(p)+len(q)-1)
 for i,a in enumerate(p):
  for j,b in enumerate(q):out[i+j]+=a*b
 return out

def _cheb_coeff_from_power(p):
 # t^k represented in first-kind Chebyshev basis, exactly
 out=[F(0)]*len(p);power=[F(1)]
 for k,pk in enumerate(p):
  for j,v in enumerate(power):out[j]+=pk*v
  if k+1<len(p):
   nxt=[F(0)]*(len(power)+1)
   for j,c in enumerate(power):
    if j==0:nxt[1]+=c
    else:nxt[j+1]+=c/2;nxt[j-1]+=c/2
   power=nxt
 return out

def _dec(x,n=12):return float(x.numerator/x.denominator) if n is None else format(float(x),f'.{n}g')

def exact_continuum_certificate(lo,hi,q,nodes,n_grid=1201):
 """All comparisons and the resulting bound are exact rational operations."""
 lo,hi,q=map(lambda x:F(str(x)),(lo,hi,q))
 nodes=[F.from_float(float(x)) for x in nodes]
 m=len(nodes)
 if not m or m%2:raise ValueError('Even polynomial degree is required for rational W')
 if not lo<hi or not F(0)<q<=1:raise ValueError('Invalid plant/rate')
 r=1-q
 if r*hi**2>=1:raise ValueError('Necessary dropout bound fails')
 mid=(lo+hi)/2;half=(hi-lo)/2
 p=[F(1)]
 for c in nodes:p=_mul(p,[mid-c,half])
 coeff=_cheb_coeff_from_power(p)
 S0=sum(map(abs,coeff));S1=sum(abs(coeff[k])*k*k for k in range(1,m+1))/half
 Wmax=(q/(1-r*hi*hi))**(m//2)
 L=Wmax*(S1+S0*m*r*hi/(1-r*hi*hi))
 peak=F(0);peak_index=None
 for i in range(n_grid):
  t=F(-1)+F(2*i,n_grid-1)
  a=mid+half*t
  v=F(0)
  for c in reversed(p):v=v*t+c
  y=abs(v)*(q/(1-r*a*a))**(m//2)
  if y>peak:peak=y;peak_index=i
 bound=peak+L*(hi-lo)/(2*(n_grid-1))
 return {
  'a1':str(lo),'a2':str(hi),'q':str(q),'period':m,
  'grid_points':n_grid,'exact_rational':True,'nodes_float_hex':[float(x).hex() for x in nodes],
  'grid_peak_decimal':_dec(peak),'lipschitz_decimal':_dec(L),
  'continuum_bound_decimal':_dec(bound),
  'continuum_margin_decimal':_dec(1-bound),
  'grid_peak_index':peak_index,
  'bound_less_than_one_exact':bool(bound<1),
  'bound_fraction':str(bound.numerator)+'/'+str(bound.denominator),
  'exact_comparison':'bound_fraction < 1 is checked using integer comparison with arbitrary precision',
  'scope':'whole interval, for exact dyadic rationals representing stored binary64 gains; assumes ideal real-arithmetic controller'
 }

def from_saved_csv(directory):
 import csv
 data=[]
 with (Path(directory)/'weighted_periodic_designs.csv').open() as f:
  for row in csv.DictReader(f):
   data.append(exact_continuum_certificate(row['a1'],row['a2'],row['q'],json.loads(row['nodes'])))
 return data

if __name__=='__main__':
 import sys
 target=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parents[1]/'results'
 certs=from_saved_csv(target)
 (target/'weighted_periodic_exact_certificates.json').write_text(json.dumps(certs,indent=2))
 print('EXACT RATIONAL CONTINUUM CERTIFICATES')
 for c in certs:print(c['a1'],c['a2'],c['q'],c['continuum_bound_decimal'],c['bound_less_than_one_exact'])
 assert all(c['bound_less_than_one_exact'] for c in certs)
