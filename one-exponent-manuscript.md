# A single semiclassical exponent governs the exponentially small scales at bosonic configuration crossings

**George Gallagher**

*Manuscript draft prepared for submission to SciPost Physics — June 12, 2026*

---

## Abstract

We study avoided crossings between two competing bosonic configurations — distinct interacting sectors exchanged by a pair-transfer coupling — and show that all exponentially small scales of the problem are governed by a single semiclassical exponent: the minimal gap, the imaginary part of the exceptional point in the complexified drive plane, the peak fidelity susceptibility, and the convergence radius of the Schrieffer–Wolff expansion. The exponent is the instanton action of the lower sheet of the $2\times 2$ matrix symbol at tuned double degeneracy, defined chart-free by the characteristic-variety condition $(h_A-E_c)(h_B-E_c)=|w(x,\chi)|^2$, and evaluated in closed form. We prove the law for the direct coupling matrix element and certify the spectral law numerically under a pre-registered, bracketed protocol across six landscapes and two coupling algebras, including a blind prediction at never-tested parameters accurate to $10^{-5}$, and through a three-layer campaign of engineered well degeneracies establishing an instanton minimum rule. An exact Feshbach reduction identifies the certified quadrature as the literal symbol of an effective scalar Jacobi operator, yielding a concrete route to a full proof. The result is the discrete, configuration-coexistence counterpart of the rigorous theory of molecular predissociation widths, realized in a setting it has not addressed, with metrological observables it does not treat.

---

## 1. Introduction

Many-body quantum systems frequently organize themselves into competing *configurations*: collective structures of distinct character — spherical and deformed shapes in nuclei [1,2], normal and intruder shell-model spaces, distinct mean-field branches in interacting bosonic ensembles — which coexist over a range of parameters and exchange the role of ground state at a crossing. The rigorous understanding of such organization remains among the deep open problems of mathematical physics; Simon's challenge to derive the shell model from first principles [3] is its emblem. Near a configuration crossing, the physically interesting quantities are uniformly *exponentially small* in the system size: the avoided-crossing gap, the parametric sensitivity of the spectrum, the analyticity radius of perturbation theory. Exponentially small quantities are precisely where uncontrolled numerics mislead and where semiclassical structure, when it exists, is most powerful.

Two mature bodies of theory border this problem without entering it. Within a *single* collective sector, exponentially small tunneling splittings between symmetry-broken mean-field wells are classical material: in the Lipkin–Meshkov–Glick (LMG) model [1] the $e^{-cN}$ splitting follows from degenerate perturbation theory in coherent states [4], and instanton evaluations — including the discreteness correction $s\to s+\tfrac12$ to the prefactor [5] — date to the 1990s. In the continuum, the Helffer–Sjöstrand theory of multiple wells [6] and its descendant, the rigorous theory of molecular predissociation [7–10], establish that resonance *widths* of $2\times2$ systems of one-dimensional Schrödinger operators are exponentially small with exponents given by Agmon distances [11] of the degenerate metric $\max(0,\min(V_1,V_2))\,dx^2$ — the Agmon geometry of the lower symbol sheet — for off-diagonal couplings as general as first-order pseudodifferential operators [9]. Separately, an active debate concerns the metrological meaning of exceptional points (EPs) and the divergence of quantum Fisher information near them [12–16].

This paper addresses the configuration-crossing problem these literatures bracket: **two confining sectors of different particle content, coupled by pair transfer, tuned to double degeneracy** — coexistence rather than decay, a gap rather than a width, on a discrete bosonic Fock chain rather than in continuum position space. Our contributions:

1. **A one-exponent law.** The four scales above share the single exponent $\mathcal{S}[H]$, the instanton action of the matrix symbol's lower sheet at degeneracy, with a closed-form chart evaluation (Sec. 2–3). The familiar local correction action is demoted to a coordinate chart of this invariant; an earlier constant-coupling chart is falsified by measurement, and the invariant *retrodicts* that falsification as necessary.
2. **One proven pillar.** For the direct coupling matrix element we prove the Agmon law with an elementary discrete argument (Theorem 3), error $\mathcal{O}(M^{2/3})$ in the logarithm.
3. **A certification standard.** The spectral law is verified under pre-registration: closed-form predictions locked before measurement, two independent extrapolation models required to bracket them, extended-precision gaps to $N=280$ ($\mathrm{gap}\sim10^{-38}$), and a strict cross-check blocker rule (Sec. 4). A blind test at never-touched parameters lands within $10^{-5}$.
4. **An instanton minimum rule.** Through engineered degenerate landscapes we establish $\Gamma=\min(S_1,S_2)$ and identify a three-layer protection mechanism — coupling-induced repulsion detuning, zero-point detuning, and symmetric hybridization with sign-definite couplings (Sec. 5).
5. **A proof skeleton.** An exact Feshbach reduction maps the two-sector problem to a scalar Jacobi double well whose symbol *is* the certified quadrature, opening a route to the full theorem by predissociation-school methods (Sec. 6).

## 2. Model and the invariant

Sector $B$ contains $M=N+2$ bosons in two modes $(s,t)$ with
$$H_B=\varepsilon\, n_t+\frac{g_2}{M}n_t^2+\frac{g_3}{M^2}n_t^3+\frac{g_4}{M^3}n_t^4-\frac{\kappa_B}{M}Q^2,\qquad Q=s^\dagger t+t^\dagger s,$$
sector $A$ contains $N$ bosons with $H_A=\varepsilon_A n_t-\frac{\kappa_A}{N}Q_A^2$, and the sectors are coupled by pair transfer $W=\omega\,(s^\dagger s^\dagger+t^\dagger t^\dagger)$ in the even-$n_t$ parity sector. The full family is $H(s)=\big(H_A\oplus(H_B+s)\big)+W$, with $s$ the drive tuned through the crossing. In interleaved ordering the matrix is pentadiagonal, enabling controlled inertia-bisection spectroscopy at arbitrary precision.

With $x=n_t/M$ and conjugate phase $\chi$, the principal matrix symbol is
$$h(x,\chi)=\begin{pmatrix} h_A & w\\ \bar w & h_B+\bar s\end{pmatrix},\quad
\begin{aligned}
h_{A,B}&=V_{A,B}(x)-2\kappa_{A,B}\,x(1-x)\,(1+\cos 2\chi),\\
w(x,\chi)&=\omega\big[(1-x)+x\,e^{2i\chi}\big],
\end{aligned}$$
with $V_B=\varepsilon x+g_2x^2+g_3x^3+g_4x^4$, $V_A=\varepsilon_A x$. The $\chi$-dependence of $w$ — forced by the operator content of $t^\dagger t^\dagger$, which moves $n_t$ — is the structural point: under the barrier ($\chi=i\lambda$),
$$|w(x,i\lambda)|^2=\omega^2\big[P(x)+R(x)\cosh 2\lambda\big],\qquad P=(1-x)^2+x^2,\;R=2x(1-x).$$

**Definition (invariant).** $\mathcal{S}[H]:=\lim_{M\to\infty}-\tfrac1M\log \mathrm{gap}_{\min}$, with the offset $\bar s$ and crossing energy $E_c$ fixed by double degeneracy of the dressed lower sheet. The candidate identity, the subject of this paper, is that $\mathcal{S}[H]$ equals $\int_0^{x_*}\lambda(x)\,dx$ along the path defined chart-free by
$$\big(h_A-E_c\big)\big(h_B+\bar s-E_c\big)=|w(x,\chi)|^2. \tag{1}$$

**Closed-form chart.** Substituting $c=\cosh 2\lambda$, Eq. (1) is *linear* in $c$ for $\kappa_A=0$,
$$c(x)=\frac{a\,B_1-\omega^2 P}{a\,K+\omega^2 R},\qquad a=\varepsilon_A x-E_c,\;K=2\kappa_B x(1-x),\;B_1=V_B+\bar s-E_c-K,$$
and *quadratic* for $\kappa_A>0$: $K_AK_B\,c^2-(A_1K_B+B_1K_A+\omega^2R)\,c+(A_1B_1-\omega^2P)=0$. The physical branch is **derived, not chosen**: the lower-sheet condition $a+b>0$ forces the smaller root (at the value $c_*$ where $a+b=0$ one has $ab=-a^2<|w|^2$, so the quadratic is negative there and $c_*$ separates the roots). Boundary data are exact: $h_A(0)=0$ gives $E_{\rm left}(\bar s)=\tfrac{\bar s}{2}-\sqrt{\tfrac{\bar s^2}{4}+\omega^2}$ in closed form.

**Self-protection lemma.** Since $E_c$ is *defined* as the global minimum of the dressed lower sheet, the real-section discriminant of (1) is nonnegative everywhere: open classically allowed pockets cannot occur at the gap-min energy. They arise only for excited crossings, delimiting the law's stated scope.

**Chart versus invariant.** Evaluating the coupling at its diagonal value ($|w|^2\to\omega^2$) defines the naive chart $S_{\rm const}$. Measurement falsifies it; the invariant formulation *requires* that falsification, since $\omega^2$ is not the gauge-invariant coupling. The surviving chart $S_{\rm corr}$, Eq. (1) with the full $\cosh$ structure, is what every certification below tests.

## 3. Proven results

**Theorem 1 (EP–fidelity lock).** In the $2\times2$ reduction at the crossing, the complexified degeneracies sit at $z_\pm=s_*\pm 2i w/\delta$ and the peak fidelity susceptibility obeys $4\chi_F^{\max}\,(\mathrm{Im}\,z_{\rm EP})^2=1$ exactly; deviations in the full model are controlled by a Feshbach transfer bound and certified to $10^{-6}$. This is the Hermitian-family cousin of relations long discussed for non-Hermitian EP sensing [12–16]; we use it only as the lock binding the observable class $\{\mathrm{gap},\,\mathrm{Im}\,z_{\rm EP},\,(4\chi_F^{\max})^{-1/2},\,r_{\rm SW}\}$ to a common exponent, not as a contribution to that debate.

**Theorem 2 (exact anchor).** In the solvable case ($\varepsilon=0,\kappa_B=1,g_i=0$) the direct coupling admits the closed form $w_{\rm direct}=2\sqrt2\,\omega\sqrt{(N{+}1)(N{+}2)}\;2^{-(N+2)/2}$; certified to relative error $1.3\times10^{-15}$.

**Theorem 3 (Agmon law for the matrix element).** Under (H1) $\kappa_B>0$, (H2) a unique nondegenerate interior minimum $x_*$ of $f=V_B-4\kappa_Bx(1-x)$, (H3) $f>E_0$ elsewhere,
$$-\tfrac1M\log\!\big(w_{\rm direct}/\omega M\big)\;\longrightarrow\;S_{\rm Agmon}=\tfrac12\int_0^{x_*}\!\operatorname{arccosh}\frac{f(x)-E_0}{2\kappa_Bx(1-x)}\,dx,$$
with error $\mathcal{O}(M^{2/3})$ in the logarithm. The proof is elementary and fully discrete: Perron–Frobenius positivity of the gauged chain, Gershgorin energy bounds, a concentration estimate, a two-sided Riccati sandwich $e^{\theta-\mu}\le R_j\le e^{\theta+\mu}$ with a contracting error recursion, and an exact boundary-layer treatment. $S_{\rm Agmon}=\mathcal{S}|_{\omega\to0}$: the matrix element and the spectral scales form two classes joined at the weak-coupling limit.

## 4. Certification protocol and results for the spectral law

*Certified* means, throughout: the closed-form prediction is computed and locked **before** measurement; gaps are obtained by banded $LDL^\top$ inertia bisection in adaptive multiprecision (working precision scaled to $\mathcal{S}N$), cross-validated against double precision wherever both apply, with any single cross-check failure blocking the series; local slopes $\Gamma_{\rm loc}(N)$ on $N=64\ldots280$ are extrapolated by two independent models, $G+b/N$ and $G+b/N+c\ln N/N$, whose pair must **bracket** the prediction.

| Configuration | Prediction $\mathcal{S}$ | Result |
|---|---|---|
| C1 ($\varepsilon{=}0,\kappa_B{=}1$), $\omega{=}0.2$ | 0.25179 | bracketed; $|\Delta|\le4\times10^{-4}$ class |
| C3 ($\varepsilon{=}-0.2,g_2{=}2$), $\omega{=}0.1$ | 0.20305 | bracketed; same class |
| **C4 blind** ($\varepsilon{=}-0.5,g_2{=}3,\kappa_B{=}0.8$), $\omega{=}0.13/0.05$ | 0.15013 / 0.19505 | $|\Delta|=1\times10^{-5}$ / $2.2\times10^{-4}$ |
| Dynamical $A$ ($\kappa_A{=}0.2$; quadratic branch rule) | 0.09308 | bracket $[0.09205,0.09383]$, $|\Delta|=1.4\times10^{-4}$ |
| Quartic dip landscape, $\omega{=}0.08$ | 0.45514 | $|\Delta|=1.8\times10^{-4}$; the sub-degenerate dip is a spectator |

Configuration C4 was constructed adversarially (first $\kappa_B\neq1$, new landscape, new couplings) after the formula was frozen; its agreement at $10^{-5}$ is the strongest single datum. A five-attack falsification battery — estimator closure against the independently measured quantum crossing slope ($|\Delta|\le2.2\times10^{-4}$), prefactor-chart covariance ($\Gamma_\infty$ spread $\le8\times10^{-5}$ with $b$ shifting by exactly one unit per power of $N$), window-split extrapolation, pipeline-branch cross-validation, and the blind test — was executed in full; all five attacks failed to break the law.

The finite-size flow obeys $d\Gamma/d\ln N=-(\Gamma-\mathcal{S})$: a line of fixed points with correction-to-scaling eigenvalue $-1$ and a non-universal amplitude $b\in[0.24,0.44]$ flowing with $\omega$ — the descendant of the LMG discreteness correction [5], whose first-principles derivation we leave open.

## 5. The instanton minimum rule

Degenerate landscapes are where naive semiclassics fails first; we engineered them deliberately, in a quartic family $f'(x)\propto(x-r_1)(x-r_b)(x-r_2)$ with exactly degenerate wells, and pre-registered each layer.

**Layer 1 — generic degeneracy.** The coupling itself lifts classical degeneracy: level repulsion $\sim\omega^2/\Delta h$ differs between wells, detuning the dressed minima *extensively*. The blind construction selects the favored well ($S_{\rm sel}=0.27997$ vs alternative $0.41920$); direct spectroscopy confirms the predicted architecture — a single ground anticrossing carrying the selected well's character ($\langle x\rangle_B=0.266$), the alternative exiled to the excited gap ($\langle x\rangle_B=0.627$) on the predicted side; measured exponent $0.27680$.

**Layer 2 — classically tuned degeneracy.** Cancelling the repulsion detuning by a calibrated tilt leaves an $\mathcal{O}(1)$ zero-point detuning, *measured* to be $0.41$ and independent of $M$ across three sizes. The gap-min remains a clean two-level anticrossing: measured $0.27369$ on $S_1=0.28068$.

**Layer 3 — quantum-tuned degeneracy.** Tuning per $N$ realizes the genuine merged three-level point. The crossing doublet hybridizes symmetrically ($\langle x\rangle_B=0.462$, the well midpoint); Perron–Frobenius positivity makes the two well-couplings *add*, forbidding destructive interference; on the co-moving resonance manifold the local exponent rises cleanly into the predicted band, $0.27975\to0.28496$ against $S_1\in[0.28424,0.28564]$. The merger floor itself is a measurement: the anticrossings cannot approach closer than twice the inter-well tunneling amplitude, and the floor's decay rate $0.115$ matches the independently estimated in-barrier action $0.106$.

**Result.** $\Gamma=\min(S_1,S_2)$ at every layer of fine-tuning; only the spectral architecture reorganizes. A methodological corollary, learned from a transient $+0.037$ exponent anomaly later dissected as an artifact: resonance manifolds that drift with $N$ must be tracked with co-moving tuning, as fixed-parameter windows sweep through resonance and contaminate fits.

## 6. Toward a proof: exact reduction to a scalar double well

For $\kappa_A=0$ the $A$-sector amplitude is exactly slaved, $\psi_A(n)=[w_1\psi_B(n)+w_2\psi_B(n{+}2)]/(d_A(n)-E)$, with denominator uniformly positive near the crossing. Feshbach elimination is therefore exact and local, collapsing the problem to a **scalar Jacobi chain** for $\psi_B$; since $(H-E)$ in the $\mathbb{Z}_2$-gauged basis is an M-matrix below the crossing and Schur complements of M-matrices are M-matrices, the effective chain inherits sign-definite hops and Perron–Frobenius structure automatically. Two exact identities anchor the reduction. First, the effective chain's forbidden-region decay condition reproduces $\omega^2[P+R\cosh2\lambda]/a$ — **the certified quadrature is the literal symbol of the effective scalar operator.** Second, the dressed boundary value $f_{\rm eff}(0)=h_B(0)+\bar s-\omega^2/(h_A(0)-E_c)=E_c$ exactly, by the degeneracy relation $E_c^2-\bar sE_c-\omega^2=0$: the eliminated $A$-branch reappears as a boundary half-well at precisely the crossing energy. The two-sector gap problem *is* a scalar double-well tunneling problem.

The remaining program is then standard-adjacent: (i) the Riccati sandwich of Theorem 3 transfers to the effective chain, with the energy dependence of its coefficients pinned to an exponentially small self-consistency window; (ii) the genuinely new estimate is a threshold boundary-well lemma at $x=0$; (iii) the final gap formula follows from a Helffer–Sjöstrand interaction matrix on the scalar double well — exactly the step where the predissociation-school machinery [7–10] applies, transplanted from width-at-crossing to gap-at-degeneracy geometry. For $\kappa_A>0$ the Schur complement remains well defined with Combes–Thomas [17] exponential off-diagonal decay, and the same architecture applies with a long-range-tolerant sandwich. We state this as a skeleton with named gaps, not a theorem.

## 7. Discussion

**Relation to prior art.** Within one sector, exponentially small splittings are LMG-classical [4,5]; in the continuum, lower-sheet Agmon exponents for $2\times2$ systems are proven theorems of predissociation theory [7–10], including pseudodifferential couplings. Our setting — coexistence of two confining configurations at tuned degeneracy, on a discrete Fock chain, with the certified collapse of spectral, exceptional-point, and metrological scales into one number, and an established selection rule under engineered degeneracy — is addressed by neither, while inheriting tools from both. We regard the present law as the configuration-coexistence realization of that circle of ideas, and Sec. 6 as the bridge by which its proof should be completed.

**Experimental relevance.** The model class is directly realizable in two-mode bosonic platforms (double-well or two-component condensates, circuit-QED Kerr resonator pairs) where both the gap and parametric susceptibility are measurable; the one-exponent collapse predicts that independently measured sensitivity and splitting must agree at exponential accuracy.

**Limitations, stated plainly.** The spectral law is certified, not proven; Theorem 3 is the sole proven pillar of its mechanism. The model class is minimal — one algebra family, one coupling operator, even-parity sector — and the distance to the many-body problems that motivate it [3] remains large. One reproducibility flag from the verification campaign (a single scratch evaluation discrepant with two agreeing production runs) remains open and is disclosed in the archive. Open problems: completion of the proof along Sec. 6; first-principles derivation of $b(\omega)$, for which the measured zero-point offset ($0.41$) and inter-well action ($0.115$) are calibration data; the barrier-washout boundary $\omega\to\omega_c$; and structural extension beyond $\mathfrak{su}(2)$ sectors.

---

**Data and code availability.** All certified numbers are reproducible from a single canonical module (`pipeline.py`) in the program archive, together with seven research documents, twelve drivers, and the complete falsification record, including the failures.

**Acknowledgments and AI disclosure.** This research program was conducted in an adversarial human–AI collaboration with Claude Fable 5 (Anthropic), which contributed derivations, proofs, all numerical pipelines, falsification design, and manuscript drafting under the author's direction; in accordance with publisher policies on authorship and accountability, it is not listed as an author. Every analytical and numerical claim was verified within the program's pre-registration protocol; the author takes full responsibility for the content.

## References

[1] H. J. Lipkin, N. Meshkov, A. J. Glick, *Nucl. Phys.* **62**, 188 (1965).
[2] K. Heyde, J. L. Wood, *Rev. Mod. Phys.* **83**, 1467 (2011).
[3] B. Simon, *J. Math. Phys.* **41**, 3523 (2000).
[4] Finite-size LMG symmetry breaking and the exponentially small gap, arXiv:1709.07657, Appendix A, and references therein.
[5] Instanton analysis of spin tunneling in the Lipkin model, arXiv:cond-mat/0102073, and the 1997 instanton literature cited therein.
[6] B. Helffer, J. Sjöstrand, *Comm. Part. Diff. Eq.* **9**, 337 (1984).
[7] M. Klein, *Ann. Phys.* **178**, 48 (1987).
[8] A. Grigis, A. Martinez, *Anal. PDE* **7**, 1027 (2014); arXiv:1205.5196.
[9] S. Fujiié, A. Martinez, T. Watanabe, *J. Diff. Eq.* **260**, 4051 (2016); **262**, 5880 (2017).
[10] S. Ashida, *Asymptot. Anal.* (2018); arXiv:1707.07963.
[11] S. Agmon, *Lectures on Exponential Decay of Solutions of Second-Order Elliptic Equations* (Princeton, 1982).
[12] W. D. Heiss, *J. Phys. A* **45**, 444016 (2012).
[13] J. Wiersig, *Phys. Rev. Lett.* **112**, 203901 (2014).
[14] D. Anderson, M. Shah, L. Fan, *Phys. Rev. Applied* **19**, 034059 (2023).
[15] Scaling of quantum Fisher information for quantum exceptional-point sensors, arXiv:2404.03803.
[16] Fundamental limits of non-Hermitian sensing from quantum Fisher information, arXiv:2603.10614 (2026).
[17] J. M. Combes, L. Thomas, *Comm. Math. Phys.* **34**, 251 (1973).

*Draft note: bracketed-author references [4,5,15,16] to be completed with full author lists at final verification; one further literature pass (Hagedorn–Joye adiabatic asymptotics; IBM configuration-mixing numerics) is scheduled before submission.*
