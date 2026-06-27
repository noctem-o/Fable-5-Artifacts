# One Exponent

### Spectral universality at quantum configuration crossings — and the adversarial method that certified it

**Fable 5 · George Gallagher**
*June 12, 2026*

> *Every claim in this document survived an attempt to kill it. Three claims did not, and their deaths are recorded here with equal care.*

---

## I. The Question

Some of the deepest unsolved problems in mathematical physics — Barry Simon's Problem 11 among them — ask why complex quantum systems organize themselves into competing *configurations*: shells, shapes, phases that coexist and exchange dominance. When two such configurations cross as a control parameter is tuned, everything interesting becomes exponentially small: the avoided-crossing gap, the sensitivity of the spectrum to perturbation, the radius within which perturbation theory converges at all. Exponentially small quantities are where intuition fails and where numerics, naively done, lie.

This program asked a narrow version of the question and answered it completely: *in a minimal two-sector bosonic model of configuration competition, how many independent exponentially small scales are there?*

The answer is **one**.

## II. The Invariant

Let two bosonic sectors with classical symbols $h_A(x,\chi)$ and $h_B(x,\chi)$ be coupled by a symbol $w(x,\chi)$, and tune the relative offset to the crossing. Define

$$\mathcal{S}[H] \;=\; \lim_{M\to\infty} -\tfrac{1}{M}\log\,\mathrm{gap}_{\min},$$

the exponent of the minimal gap. The central result is that $\mathcal{S}[H]$ equals the **instanton action of the matrix symbol's lower sheet** at double degeneracy, characterized chart-free by the characteristic-variety condition

$$(h_A - E_c)\,(h_B - E_c) \;=\; |w(x,\chi)|^2,$$

with the crossing energy $E_c$ and offset fixed by demanding degenerate dressed minima. Four scales that look unrelated share this single exponent:

| Scale | Meaning |
|---|---|
| $\mathrm{gap}_{\min}$ | the avoided-crossing splitting |
| $\mathrm{Im}\,z_{\mathrm{EP}}$ | distance of the exceptional point into the complex parameter plane |
| $(4\chi_F^{\max})^{-1/2}$ | the peak quantum-Fisher sensitivity, inverted |
| Schrieffer–Wolff radius | convergence boundary of the perturbative expansion |

The first three are locked pairwise by a *proven* identity (Theorem 1, below). The familiar local correction action $S_{\mathrm{corr}}$ turns out to be merely a Darboux-chart evaluation of $\mathcal{S}[H]$ — and this demotion carries explanatory force: an earlier constant-symbol formula $S_{\mathrm{const}}$ was falsified by measurement, and the invariant *retrodicts* that falsification as necessary. It was the gauge-dishonest chart. A framework that explains its own past failures as required is the signature of structure found rather than curve fitted.

Two structural facts complete the statement. The physical branch of the quadrature is **derived, not chosen**: the sheet condition $a+b>0$ selects the smaller cosh-root, a one-line proof. And a self-protection lemma holds: because $E_c$ is *defined* as the global dressed minimum, the forbidden-region discriminant satisfies $\Delta(x)\ge 0$ everywhere — open classical pockets cannot occur at the gap-min energy, only at excited crossings.

## III. What Is Proven

**Theorem 1 (EP–QFI identity).** $4\chi_F^{\max}\,(\mathrm{Im}\,z_{\mathrm{EP}})^2 = 1$ exactly in the two-level reduction; certified numerically to $10^{-6}$ with the Feshbach transfer bound accounting for the deviation.

**Theorem 2 (Exact anchor).** A closed form for the direct coupling $w_{\mathrm{direct}}$ in the solvable case; certified to relative error $1.3\times10^{-15}$.

**Theorem 3 (Agmon law).** $-\tfrac{1}{M}\log(w_{\mathrm{direct}}/\omega M) \to S_{\mathrm{Agmon}}$, with a complete elementary proof — Perron–Frobenius positivity, a two-sided Riccati sandwich, exact boundary-layer treatment — and an honestly flagged $\mathcal{O}(M^{2/3})$ error in the logarithm.

## IV. What Is Certified

*Certified* means: predicted in closed form **before** measurement, then bracketed by two independent extrapolation models on extended-precision data to $N=280$, gaps as small as $10^{-38}$.

| Test | Prediction | Outcome |
|---|---|---|
| C1, C3 landscapes | $S_{\mathrm{corr}}$ | brackets contain prediction, $\le 4\times10^{-4}$ |
| **Blind kill test C4** ($\kappa=0.8$, never touched) | 0.15013 / 0.19505 | $\|\Delta\| = 1\times10^{-5}$ / $2.2\times10^{-4}$ |
| Dynamical A-sector (branch rule) | 0.09308 | bracket [0.09205, 0.09383], $\|\Delta\|=1.4\times10^{-4}$ |
| Quartic dip landscape | 0.45514 | $\|\Delta\| = 1.8\times10^{-4}$; the dip is a spectator |

A five-attack falsification battery — estimator closure, basis dependence, pseudo-convergence, pipeline branching, parameter leakage — was executed rather than argued; all five survived, the leakage attack by the blind C4 hit above. The finite-size flow obeys the simplest possible renormalization equation, $d\Gamma/d\ln N = -(\Gamma - \mathcal{S})$: a line of fixed points with correction-to-scaling eigenvalue $-1$ and a non-universal amplitude that flows with the coupling.

## V. The Degeneracy Campaign

The last redoubt of doubt was instanton competition: what happens when the landscape offers *two* wells? Three layers of fine-tuning, three certifications.

**Layer 1 — generic.** At exact classical degeneracy, the coupling itself lifts it: $\omega^2$-repulsion detunes the dressed wells *extensively*. The blind construction selects a well; direct spectroscopy confirmed the predicted architecture — the ground anticrossing carries the selected well's character ($\langle x\rangle_B = 0.266$), the loser's anticrossing exiled to the excited gap, on the predicted side. Measured exponent 0.27680 against the selected action 0.27997, decisively not the alternative 0.41920.

**Layer 2 — classically tuned.** Cancel the repulsion detuning by hand, and an $\mathcal{O}(1)$ zero-point detuning remains — *measured*, 0.41, independent of system size across three $N$. The min-rule held again: 0.27369 on $S_1 = 0.28068$.

**Layer 3 — quantum tuned.** Tune per-$N$ to the genuine merged three-level point. The crossing states hybridize symmetrically ($\langle x\rangle_B = 0.462$, the midpoint); Perron–Frobenius positivity makes the two couplings *add*, forbidding destructive interference; and on the co-moving resonance manifold the exponent climbs cleanly into the predicted band: $0.27975 \to 0.28496$ against $S_1 \in [0.28424, 0.28564]$. As a bonus, the merger floor measures the inter-well tunneling action itself: rate 0.115 against a back-of-envelope 0.106.

$\Gamma = \min(S_1, S_2)$, at every layer. Certified, not conjectured.

## VI. The Method

The numbers above are only as good as the discipline beneath them. Five faults were caught and root-caused along the way — a falsified Fubini–Study exponent, a retracted $\omega^2$ law, the required death of $S_{\mathrm{const}}$, a landscape-builder sign error invisible to prediction-versus-measurement agreement, and a localization divergence diagnosed by hypothesis and confirmed by cure ($\mathrm{rel}\;8.36 \to 6.1\times10^{-8}$). One apparent exponent anomaly was dissected within the session into a named artifact: fixed parameters swept through an $N$-dependent resonance. And one fabricated external review of the program — fluent, citation-laden, wrong where it mattered — was rejected on the record's own evidence.

Four rules, in order of scar tissue:

1. **Pre-register before measuring.** Brackets decide; narrative does not.
2. **Verify designs upstream of both prediction and measurement.** Agreement cannot catch a bug they share.
3. **Resonance manifolds need co-moving tuning.** Fixed-parameter sweeps contaminate the fit.
4. **One failed cross-check blocks a series.** No exceptions; the exception we almost made is why the rule exists.

## VII. What Remains Open

The spectral law is certified, not proven — Theorem 3 is the lone proven pillar of its mechanism, and the $\mathbb{Z}_2$-gauged coupled-chain Riccati bound (F5) is the route to the rest. The subprincipal amplitude $b(\omega)$ awaits derivation, with the zero-point offset and inter-well action above as calibration data. One reproducibility flag stands, deliberately unerased. And the distance from this minimal model class to Simon's Problem 11 remains what it always was: large, and honestly stated.

## Coda

What this program leaves behind is smaller than a solved problem and rarer: a theory in which every standing claim was shot at first — by its authors — and a record in which the three claims that fell are preserved as carefully as the ones that stand. The single exponent was the result. The habit of trying to destroy it was the point.

---

*Data and code: the complete archive — seven documents, twelve scripts, every certified number reproducible from `pipeline.py` — accompanies this paper. Fable 5 is Claude (Anthropic). The collaboration was adversarial by design and is the better for it.*
