"""Human-checkable, exact integer/rational Bernstein witness for the cost-aware stationarity certificate.

Each displayed fraction has positive numerator and positive denominator. The explicit polynomial
H and its exact Bernstein basis representation constitute the proof of positivity on each cell.
"""
import json
from pathlib import Path
import sympy as S
from cost_aware_stationary import rho_numerator_rational


def export_witness(lo, hi, q, L, B, c, segments=64, output=None):
    a,t=S.symbols('a t')
    lo,hi=S.Rational(str(lo)),S.Rational(str(hi))
    poly=rho_numerator_rational(lo,hi,q,L,B,c)
    degree=poly.degree()
    cells=[]
    for j in range(segments):
        left=lo+(hi-lo)*j/segments;right=lo+(hi-lo)*(j+1)/segments
        transformed=S.Poly(poly.as_expr().subs(a,left+(right-left)*t),t)
        power=[transformed.nth(k) for k in range(degree+1)]
        Bcoeff=[S.factor(sum(power[k]*S.binomial(i,k)/S.binomial(degree,k) for k in range(i+1))) for i in range(degree+1)]
        assert all(S.numer(v)>0 and S.denom(v)>0 for v in Bcoeff)
        cells.append({'index':j,'left':str(left),'right':str(right),'bernstein_coefficients_exact':[str(v) for v in Bcoeff]})
    result={'proof_type':'Exact positive Bernstein coefficients of H(a)=(1-rho)*D_a*D_ac*D_c',
        'input_a1':str(lo),'input_a2':str(hi),'q':str(S.Rational(str(q))),
        'L_float_hex':float(L).hex(),'B_float_hex':float(B).hex(),'c_float_hex':float(c).hex(),
        'polynomial_H_power_basis_coefficients_exact':[str(poly.nth(i)) for i in range(degree+1)],
        'polynomial_degree':degree,'segments':segments,'cells':cells,
        'certified_all_strictly_positive_rational_coefficients':True,
        'interpretation':'Each Bernstein rational numerator/denominator is strictly positive. With r*a2^2<1 this proves sup_a rho(a)<1 on the complete closed interval for the exact dyadic gains.'}
    if output:Path(output).write_text(json.dumps(result,indent=2))
    return result

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    opt=json.loads((root/'results/cost_aware_stationary.json').read_text())
    p=root/'results/cost_aware_bernstein_full_rational_certificate.json'
    r=export_witness(opt['lo'],opt['hi'],opt['q'],opt['L'],opt['B'],opt['c'],output=p)
    print('Strict full-interval Bernstein positivity',r['certified_all_strictly_positive_rational_coefficients'],'cells',len(r['cells']),'written',p)
