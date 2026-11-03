from pathlib import Path
r=Path(__file__).resolve().parent
src=(r/'manuscript/v3_reference.tex').read_text()
src=src.replace(r'\documentclass[11pt]{article}',r'\documentclass[preprint,12pt,numbers]{elsarticle}')
src=src.replace(r'\usepackage[margin=0.95in]{geometry}'+'\n','')
src=src.replace(r'\usepackage[numbers,sort&compress]{natbib}'+'\n','')
old=src[src.index(r'\title{'):src.index(r'\begin{abstract}')]+r'\begin{abstract}'
new=r'''\journal{IFAC Journal of Systems and Control}
\begin{document}
\begin{frontmatter}
\title{Robust Mean-Square Stabilization of an Uncertain Scalar Pole under Bernoulli Erasures: Controller-Class Limits, Certified Periodic Policies, and Quadratic Performance}
\author[uiz]{Abdelilah Ganmati}
\ead{a.ganmati@uiz.ac.ma}
\address[uiz]{ENSA Agadir, Ibn Zohr University, Agadir, Morocco}
\begin{abstract}'''
assert old in src
src=src.replace(old,new)
abstractend=r'''\end{abstract}
\begin{keyword}
Mean-square stability \sep Packet erasures \sep Simultaneous stabilization \sep Periodic control \sep Exact rational certification \sep Quadratic performance
\end{keyword}
\end{frontmatter}'''
src=src.replace(r'\end{abstract}',abstractend,1)
src=src.replace('We additionally derive closed-form asymptotic state and input second moments for the stationary gain table and exact periodic-cycle moments for impulse schedules, separating them from paired finite-horizon simulations.', 'We additionally derive exact physical-time quadratic costs for a parametric family of stationary last-reception controllers, and exhibit an independently certified cost-aware stationary baseline that improves worst-case asymptotic cost while preserving continuum stability. Exact reception-cycle moments and finite-horizon simulations are kept distinct.')
src=src.replace('The analysis reveals trade-offs among stability, transient amplification, and actuation effort, but it does not establish optimality over all periodic, adaptive, or full-history controllers.','The results expose stability--performance trade-offs, but neither the cost-aware stationary design nor the periodic designs is claimed globally optimal over general controller classes.')
src=src.replace('periodic feedback compensation and its potential advantages are also classical \\citep{khargonekar1985periodic}.', 'periodic feedback compensation and its potential advantages are also classical \\citep{khargonekar1985periodic}.')
src=src.replace('These settings differ materially from the fixed scalar interval with fully observed received samples and unbounded geometric erasure gaps studied here.', 'These settings differ materially from our fixed scalar interval with fully observed received samples and unbounded geometric erasure gaps. Classical analytic simultaneous stabilization, including modern interpolation formulations \\citep{cui2024}, is important prior work; no general novelty is claimed for interpolation itself. A 2026 study of packet-loss-adaptive triggering likewise considers a different finite-bit event-triggered information structure \\citep{ding2026}.')
src=src.replace('(ii) we give an exact closed-form physical-time stationary second-moment cost identity;', '(ii) we give exact physical-time moments for a tractable family of stationary controllers and a separately certified cost-aware comparison;')
src=src.replace('Our originality claim is therefore narrower than a fundamental limit for all causal controllers. A full priority check of the stationary theorem remains a prerequisite to submission.', 'Our originality claim is restricted to the specific fixed-parameter, last-reception controller class; establishing bibliographic priority against all broader simultaneous-stabilization results remains a prerequisite to submission.')
# Condense definition to remove imprecise q-dependent δ corollary unless independently vetted; preserve mathematical depth in ancillary proof appendix instead.
start=src.index('\\begin{corollary}[per-plant gap with a known rate]')
end=src.index('\\subsection{A different information structure:',start)
archival=src[start:end]
(r/'docs/archived_rate_aware_corollary_for_external_check.tex').write_text(archival)
src=src[:start]+r'''\begin{remark}[Scope of rate awareness]
The explicit threshold-attaining sequence in Theorem~\ref{thm:main} does not depend on $q$. This fact concerns feasibility at a \emph{fixed} $q$ and does not establish a universal performance optimality theorem or settle gap-to-plant-specific-bound questions under every possible cross-rate quantifier. The earlier draft's auxiliary rate-gap corollary is retained in the research archive for separate proof review and is not required for the present article.
\end{remark}

'''+src[end:]
newblock=r'''
\subsection{Cost-aware stationary design and exact interval certification}\label{sec:costaware}
For a fixed design pole $c\in[a_1,a_2]$ consider a tractable age-indexed stationary family
\begin{equation}\label{eq:geomfamily}
h_0=L,\qquad h_j=Bc^{j-1}\quad(j\ge1),
\end{equation}
where $L,B$ are independent of the unknown $a$. This family contains the critical-rate design in Theorem~\ref{thm:main} upon taking $c=L=\bar a$ and $B=d^2/4$.

\begin{proposition}[Exact stationary cost and polynomial robustness test]\label{prop:geomcost}
Let $r=1-q$, $D_a=1-ra^2$, $D_{ac}=1-rac$, $D_c=1-rc^2$, and assume $ra^2<1$, $rc^2<1$. For $a\ne c$, put $U=a-L-B/(a-c)$ and $V=B/(a-c)$. The homogeneous reception multiplier of \eqref{eq:geomfamily} equals
\begin{equation}\label{eq:geomrho}
\rho_{L,B,c}(a,q)=q\left(\frac{U^2}{D_a}+\frac{2UV}{D_{ac}}+\frac{V^2}{D_c}\right).
\end{equation}
At $a=c$, its continuous value is
\begin{equation}\label{eq:geomcenter}
q\left[\frac{(a-L)^2}{D_a}-\frac{2(a-L)Bra}{D_a^2}+\frac{B^2r(1+ra^2)}{D_a^3}\right].
\end{equation}
Whenever $\rho_{L,B,c}(a,q)<1$, set $\nu=\sigma^2/D_a$ and $Z=\nu/(1-\rho_{L,B,c})$. The exact limiting physical-time moments are
\begin{equation}\label{eq:geomcost}
\mathbb Ex^2=(q+r\rho_{L,B,c})Z+r\nu,\qquad
\mathbb Eu^2=qZ\left(L^2+\frac{rB^2}{D_c}\right).
\end{equation}
Define $N(a)=(a-L)(a-c)-B$ and the polynomial
\begin{align}\label{eq:Hpoly}
H(a)=\frac{1}{(a-c)^2}\bigl[&(a-c)^2D_aD_{ac}D_c\\[-2mm]
&-q\bigl(N(a)^2D_{ac}D_c+2N(a)BD_aD_c+B^2D_aD_{ac}\bigr)\bigr].\nonumber
\end{align}
The numerator is divisible by $(a-c)^2$ identically; $H$ is a real polynomial of degree at most three and
\begin{equation}\label{eq:Hidentity}H(a)=(1-\rho_{L,B,c}(a,q))D_aD_{ac}D_c\end{equation}
for all $a$ including $a=c$ by continuity. Consequently, if $ra_2^2<1$ and $H(a)>0$ on $[a_1,a_2]$, then \eqref{eq:geomfamily} robustly MS-stabilizes the complete interval. Strict positivity may be certified by expressing $H$ in Bernstein form on a finite partition and checking that every coefficient is strictly positive.
\end{proposition}
\begin{proof}
On a reception gap of length $g$, the homogeneous state coefficient is
$M_g=(a-L)a^{g-1}-B\sum_{j=1}^{g-1}a^{g-1-j}c^{j-1}$.
For $a\ne c$ the geometric divided difference equals $(a^{g-1}-c^{g-1})/(a-c)$, yielding $M_g=Ua^{g-1}+Vc^{g-1}$, and summing $q\sum_{g\ge1}r^{g-1}M_g^2$ gives \eqref{eq:geomrho}. When $a=c$, $M_g=(a-L)a^{g-1}-B(g-1)a^{g-2}$ and the standard identities $\sum_{n\ge0}nt^n=t/(1-t)^2$ and $\sum_{n\ge0}n^2t^n=t(1+t)/(1-t)^3$ give \eqref{eq:geomcenter}. The reception variance recursion is $Z_{n+1}=\rho Z_n+\nu$; the geometric-age physical-time decomposition used in Proposition~\ref{prop:cost} gives \eqref{eq:geomcost}, including the exact $j\ge1$ input-gain series. Combining the terms of \eqref{eq:geomrho} over $D_aD_{ac}D_c$ and clearing the apparent $(a-c)^{-2}$ pole gives \eqref{eq:Hpoly}. Substitution $a=c$ cancels both the numerator and its first derivative, so polynomial division is exact; the identity \eqref{eq:Hidentity} then extends by continuity. Under $ra_2^2<1$ and $c\in[a_1,a_2]$, each factor $D_a,D_{ac},D_c$ is positive throughout the interval. Finally, the convex-hull property of a Bernstein representation of a polynomial on any subinterval places its value between the minimum and maximum Bernstein coefficients; strictly positive coefficients suffice for strict positivity on the closed subinterval, and a finite covering proves the claim.
\end{proof}

We optimize the \emph{worst-case exact asymptotic} $J=\mathbb Ex^2+0.05\mathbb Eu^2$ within the two-parameter family $(L,B)$ with $c=2.1$ and $(a_1,a_2,q)=(1.2,3.0,0.96)$. This is an epigraph-based \emph{local} nonlinear solve with 31 interval nodes; no global optimality of the cost solution is asserted. The stored binary64 gains are then treated as exact dyadic rational numbers in an independent symbolic polynomial-and-Bernstein verifier using 64 subintervals. Its nonnegative-integer coefficient comparisons are exact and certify $H>0$ \emph{over the continuum}, independently of the numerical optimizer. This is a certificate for the exact dyadic controller values, not an assertion of floating-point hardware robustness.

\begin{table}[t]\centering\small
\caption{Matched stationary policies on $[1.2,3.0]$, $q=0.96$, $\sigma^2=0.01$, weight $0.05$. The displayed worst-interval exact cost values were independently checked on a 1001-point grid (not globally certified cost maxima). Both stability certificates are interval-wide; the lower cost is not a theorem of cost optimality.}\label{tab:costaware}
\begin{tabular}{@{}lrrrc@{}}\toprule
Policy & $L$ & $B$ & max-grid $J_\infty$ & Continuum MS stability \\\midrule
Critical-rate $h^*$ & 2.100000 & 0.810000 & 0.339402 & analytic \\
Cost-aware local & 2.104692 & 0.814349 & 0.280345 & exact Bernstein \\
\bottomrule\end{tabular}
\end{table}
\begin{figure}[t]\centering
\includegraphics[width=.88\linewidth]{fig05_cost_aware_stationary.png}
\caption{Exact limiting quadratic costs for the two stationary controllers under identical noise and channels. Values are deterministic evaluations of Proposition~\ref{prop:geomcost}, not finite-horizon sample means. The optimized policy is a \emph{locally optimized grid design}, with independently certified continuum stability.}\label{fig:costaware}
\end{figure}
'''
anchor=r'\subsection{Rate-aware weighted polynomial periodic baseline}'
assert anchor in src
src=src.replace(anchor,newblock+'\n'+anchor)
src=src.replace('with seed 20261020, $\\sigma^2=0.01$', 'with seed 20261020, $\\sigma^2=0.01$')
src=src.replace('All results in this section were regenerated from the v3 source on a CPU', 'The original v3 comparison results were regenerated from the preserved v3 source on a CPU')
src=src.replace('For $\\Theta=[1.5,4.5]$ at $q=0.98$', 'The cost-aware stationary experiment uses a separate seed 20261064, 3500 paired trajectories and horizon 500; its time-domain paired simulation is exported separately from its deterministic asymptotic moment identities. At the upper plant endpoint $a=3$, the exact asymptotic cost changes from $0.339402$ to $0.280345$; this reduction is offset by a cost increase at $a=1.2$ (from $0.230494$ to $0.280345$). The policy is designed against the maximum over the interval, not a selected favorable endpoint.\n\nFor $\\Theta=[1.5,4.5]$ at $q=0.98$')
src=src.replace('We do not report imaginary reproductions of them.', 'We do not report reproductions of these incommensurate baselines.')
src=src.replace('The strongest conclusion supported by this work is a sharp controller-class-specific robust feasibility threshold, plus exact interpretable moment and performance comparisons.', 'The strongest conclusion supported by this work is a sharply delimited controller-class-specific robust feasibility threshold, together with certified periodic alternatives and exact interpretable moment and performance comparisons. Within the stationary family, the new cost-aware design quantifies the difference between optimizing feasibility at a critical rate and optimizing a quadratic performance objective at a fixed feasible rate.')
src=src.replace('An independent mathematical priority review, and research-priority verification remain required before submission;', 'Independent proof and mathematical-priority reviews remain required before submission;')
src=src.replace('For the class of stationary last-observation linear controllers, a tight robust mean-square critical reception probability is analytically tractable.', 'For the restricted class of stationary last-observation linear controllers, a tight robust mean-square critical reception probability is analytically tractable; the feasibility-optimal gain schedule need not minimize quadratic cost. A locally cost-aware schedule improves tested worst-interval quadratic performance while preserving a separately certified continuum contraction condition.')
src=src.replace('\\paragraph{Generative AI disclosure}\nThe research draft and computational implementation were developed with AI assistance, under author responsibility for mathematical verification, data integrity and final manuscript approval. The exact form of any publisher-required disclosure will be confirmed using the accessible submission instructions.', '\\paragraph{Declaration of generative AI and AI-assisted technologies}\nAI assistance was used in code development, manuscript drafting and computational checking; all claims, calculations, authorship and final editorial responsibility remain with the human author. The final wording and any required submission-form disclosure are subject to the author\'s verification against the publisher\'s current policy.')
src=src.replace('\\bibitem[Okano and Ishii(2015)]{okano2015} K.~Okano, H.~Ishii. Stabilization of uncertain systems using quantized and lossy observations and uncertain control inputs. arXiv:1506.03518, 2015.', '\\bibitem[Okano and Ishii(2017)]{okano2015} K.~Okano, H.~Ishii. Stabilization of uncertain systems using quantized and lossy observations and uncertain control inputs. \\emph{Automatica}, 81:261--269, 2017. doi:10.1016/j.automatica.2017.03.036.')
src=src.replace('\\bibitem[Tannenbaum(1980)]{tannenbaum1980}', '\\bibitem[Okano and Ishii(2013)]{okano2013} K.~Okano, H.~Ishii. Networked control of uncertain systems over data rate limited and lossy channels. \\emph{IEICE Trans. Fundamentals}, E96-A(5):853--860, 2013. doi:10.1587/transfun.E96.A.853.\n\\bibitem[Cui and Lindquist(2024)]{cui2024} Y.~Cui, A.~Lindquist. A general solution to the simultaneous stabilization problem by analytic interpolation. arXiv:2405.13440, 2024.\n\\bibitem[Ding et al.(2026)]{ding2026} X.~Ding, Y.~Zhang, J.~Xue, J.~Zhou. Mean-square exponential ultimate boundedness of networked control systems with packet-loss-adaptive event triggering and scale-normalized encoding. \\emph{Electronic Research Archive}, 34(11):8689--8718, 2026. doi:10.3934/era.2026366.\n\\bibitem[Tannenbaum(1980)]{tannenbaum1980}')
# Elsevier class handles natbib; ensure no figure parse mistakes
# Remove duplicated Cui--Lindquist citation introduced by v3 bibliographic carryover.
item=r'\bibitem[Cui and Lindquist(2024)]{cui2024} Y.~Cui, A.~Lindquist. A general solution to the simultaneous stabilization problem by analytic interpolation. arXiv:2405.13440, 2024.'
assert src.count(item)==2,src.count(item)
src=src.replace(item+'\n','',1)
# Align section labels with journal's prescribed major-section names.
src=src.replace(r'\section{Methods: plant, information, and controller classes}',r'\section{Methods}'+'\n'+r'\subsection{Plant, information, and controller classes}')
src=src.replace(r'\section{Results: exact stationary-class bounds}',r'\section{Results}'+'\n'+r'\subsection{Exact stationary-class bounds}')
# Remove unused, not-independently-vetted cross-rate performance-gap definitions.
qstart=src.index(r'\begin{align*}'+'\n'+r'\qrob(\Theta)&=')
qend=src.index(r'\end{align*}',qstart)+len(r'\end{align*}')
qreplacement=r'\[\qrob(\Theta)=\inf\{q\in(0,1):\exists h\in\Lone\text{ MS-stabilizing every }a\in\Theta\text{ at rate }q\}.\]'
src=src[:qstart]+qreplacement+src[qend:]
oldq='Thus $\\qrob$ is the feasible class critical reception rate, not a real-valued algebraic cutoff if no admissible rate exists; $\\delta^*_{\\rm aware}$ and $\\delta^*_{\\rm univ}$ measure how close a single design comes to every plant\'s own critical rate, with and without knowledge of $q$.'
assert oldq in src
src=src.replace(oldq,'Thus $\\qrob$ is the feasible class critical reception rate, not a real-valued algebraic cutoff if no admissible rate exists.')
# Include two especially close classical dropout-channel synthesis papers.
oldintro='Necessary channel loss conditions are classical for known plants'
src=src.replace(oldintro,'Robust stochastic feedback synthesis over erasure channels was addressed by Elia and by Silva and Pulgar \\citep{elia2005,silvapulgar2011}; separate necessary channel loss conditions are classical for known plants')
oldbib=r'\bibitem[Elia(2005)]{elia2005} N.~Elia. Remote stabilization over fading channels. \emph{Systems \& Control Letters}, 54(3):237--249, 2005.'
assert oldbib in src
src=src.replace(oldbib,oldbib[:-1]+', doi:10.1016/j.sysconle.2004.08.009.')
ref=r'\bibitem[Silva and Pulgar(2011)]{silvapulgar2011} E.~I.~Silva, S.~A.~Pulgar. Control of LTI plants over erasure channels. \emph{Automatica}, 47(8):1729--1736, 2011. doi:10.1016/j.automatica.2011.02.050.'
where=r'\bibitem[Elia and Eisenbeis(2011)]{eliaeisenbeis2011}'
assert where in src
src=src.replace(where,ref+'\n'+where,1)
(r/'manuscript/IFAC_consolidated_submission_draft.tex').write_text(src)
print('new document chars',len(src),'lines',src.count('\n'))
