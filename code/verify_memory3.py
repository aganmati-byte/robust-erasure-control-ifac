"""Independent exact rational Lyapunov verification for the three-age two-measurement policy.
Exact endpoint certificates only: NOT a continuum proof.
"""
import sympy as S

def moment_coeff(a,q,h,k):
 a=S.Rational(a);q=S.Rational(q);r=1-q
 alpha=S.Integer(1);beta=S.Integer(0);alphas=[];betas=[]
 for j in range(3):
  alpha=a*alpha-h[j];beta=a*beta-k[j];alphas.append(S.factor(alpha));betas.append(S.factor(beta))
 if not r*a*a<1:raise ValueError('dropout noise moment does not exist')
 def moment(p):
  return S.factor(q*sum(r**j*p(alphas[j],betas[j]) for j in range(2))
                   +q*r*r*p(alphas[2],betas[2])/(1-r*a*a))
 # alpha and beta first moments have different geometric tails
 ealpha=S.factor(q*sum(r**j*alphas[j] for j in range(2))+q*r*r*alphas[2]/(1-r*a))
 ebeta=S.factor(q*sum(r**j*betas[j] for j in range(2))+q*r*r*betas[2]/(1-r*a))
 return moment(lambda x,y:x*x), moment(lambda x,y:x*y),moment(lambda x,y:y*y),ealpha,ebeta

def lyapunov_certificate(a,q,h,k):
 A,B,C,D,E=moment_coeff(a,q,h,k)
 # L(P) = [[pA+2sD+t, pB+sE],[pB+sE, pC]]
 p,s,t=S.symbols('p s t')
 sol=S.solve([p-(p*A+2*s*D+t)-1,s-(p*B+s*E),t-p*C-1],[p,s,t],dict=True)[0]
 P=S.Matrix([[sol[p],sol[s]],[sol[s],sol[t]]]); Q=S.eye(2)
 assert S.simplify(P-S.Matrix([[sol[p]*A+2*sol[s]*D+sol[t],sol[p]*B+sol[s]*E],[sol[p]*B+sol[s]*E,sol[p]*C]]))==Q
 checks={'a':str(a),'P11':str(P[0,0]),'P12':str(P[0,1]),'P22':str(P[1,1]),'P11_positive':bool(P[0,0]>0),'detP_positive':bool(P.det()>0),'exact_lyapunov_difference_identity':True}
 return checks

if __name__=='__main__':
 import json
 H=list(map(S.Rational,['1977/1000','143/125','741/250']))
 K=list(map(S.Rational,['21/100','179/1000','9/200']))
 out=[lyapunov_certificate(S.Rational(a),S.Rational('9/10'),H,K) for a in ['6/5','3']]
 print(json.dumps(out,indent=2))
 assert all(z['P11_positive'] and z['detP_positive'] for z in out)
