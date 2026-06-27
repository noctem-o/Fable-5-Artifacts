# Configuration Competition, Not Cluster Stability — and Not Quite a Bifurcation Either
## A reality-calibrated formulation of the shell-evolution theorem

**Status key:** [EXACT] = provable now with standard tools · [DOABLE] = concrete open theorem, existing technology · [HARD] = genuinely open · [SYNTH] = this document's own synthesis, not in the literature as stated.

---

## 1. Verdict on the draft formulation

What the draft gets **right** (and this is the important part):

- The robust phenomenon far from stability *is* configuration coexistence and reordering, not persistence of a fixed cluster decomposition. This is the consensus of the shell-evolution literature (Otsuka–Gade–Sorlin–Suzuki–Utsuno, *Rev. Mod. Phys.* **92**, 015002 (2020); Heyde–Wood, *Rev. Mod. Phys.* **83**, 1467 (2011); Caurier et al., *Rev. Mod. Phys.* **77**, 427 (2005)).
- The driver is monopole evolution (tensor + 3N; Zuker, *PRL* **90**, 042502 (2003)) competing against multipole correlation energy — correct.
- The right object is an order parameter measuring configuration character, not a global spectral invariant — correct in spirit.

What must be **repaired** before it is "closest to reality":

1. **`E_gs(s) = min{E_norm(s), E_intr(s)}` is false as a statement about eigenvalues.** With coupling `W ≠ 0`, both sectors contain `J^π = 0⁺` states (2p–2h intruders have the same quantum numbers as the normal configuration), so by von Neumann–Wigner the levels *avoid* crossing: `E_gs(s)` is analytic (Kato–Rellich) and strictly below `min{E_norm, E_intr}`. The `min` formula is exact only for the *diabatic*, sector-restricted branches — which are not eigenvalues of `H`. The draft conflates the two. The repair: state both objects and their exact relation (Section 2).
2. **For finite A there is no bifurcation.** The ground state changes smoothly; "the order parameter changes rapidly" is a *crossover* with a computable width, not a nonanalyticity. A genuine bifurcation exists only in a declared limit (large valence degeneracy / boson number / `ħ_eff → 0`). The repair: a two-tier theorem — finite-A crossover theorem + classical-limit bifurcation theorem (Sections 2 and 4).
3. **`β(s) = ⟨Ψ, Q̂ Ψ⟩` is identically zero.** For any `J = 0` state, the expectation of a rank-2 spherical tensor vanishes by Wigner–Eckart. The deformation order parameter must be a rotational invariant: `q² = ⟨0⁺|Q·Q|0⁺⟩` and the higher Kumar–Cline invariants (Kumar, *PRL* **28**, 249 (1972); Cline, *Annu. Rev. Nucl. Part. Sci.* **36**, 683 (1986)); see Poves–Nowacki–Alhassid, *Phys. Rev. C* **101**, 054307 (2020) ("Limits on assigning a shape to a nucleus") for why even this requires care — fluctuations of the invariants are large, so "shape" is a distribution, not a number. [EXACT]
4. **The two-sector block matrix is a cartoon of a deeper, partly self-consistent structure.** Real cases involve a hierarchy 0p–0h, 2p–2h, 4p–4h, … (⁴⁰Ca, ⁶⁸Ni, Zr chain), and **Type-II shell evolution**: the effective single-particle energies inside the intruder sector differ from those in the normal sector because the monopole field depends on the occupations of the configuration itself. The repair: lower-envelope-of-many-branches structure with sector-dependent (self-consistent) blocks (Section 6).

One more honesty item the draft inherits silently: **the sectors are scheme-dependent.** The projections `P_norm, P_intr` depend on a chosen single-particle basis and Fermi surface, and unitary flows (IMSRG!) move them. The theorem must therefore be stated for a *fixed* pair (effective Hamiltonian, partition), with all physical claims routed through invariants: energies of physical `0⁺` states, `S_2n`, `ρ²(E0)`, `B(E2)`, `q²`. Diabatic branches are scaffolding, not observables.

---

## 2. The exact finite-A layer [EXACT]

**Setup.** Valence Fock space `𝓕`; a smooth (or real-analytic) family of self-adjoint effective Hamiltonians `H(s)`, `s` a control parameter (interaction strength, interpolated neutron number, quadrupole coupling). Fix orthogonal projections `P₀ ⊕ P₂ ⊕ … = 𝟙` onto configuration sectors (0ℏω normal, 2p–2h intruder, …), commuting with `J²,J_z`, parity. Define the **diabatic branches**

> `E_k(s) := inf spec( P_k H(s) P_k ↾ Ran P_k )`,

and the coupling `w(s) := ⟨φ₀(s), H(s) φ₂(s)⟩` between the sector ground vectors.

**Proposition 0 (envelope and non-crossing).**
(i) `E_gs(s) ≤ min_k E_k(s)` for all `s` (Rayleigh–Ritz: each sector supplies trial states).
(ii) The branches `E_k(s)` are smooth and may cross transversally at some `s*` (restricted minima feel no level repulsion).
(iii) If `w(s*) ≠ 0`, the true eigenvalues anticross: the lowest two `0⁺` levels of `H(s)` are analytic in `s` with minimal gap `2|w(s*)| + O(|w|²/g)` at the avoided crossing, where `g` is the smaller of the two in-sector excitation gaps. The exact relation between branches and spectrum is the Feshbach–Schur (Löwdin) map: `E_gs` solves `E = inf spec[ P₀HP₀ + P₀HP₂ (E − P₂HP₂)⁻¹ P₂HP₀ ]`.

**Theorem A — Finite-A Crossover Theorem.** Assume on a neighborhood of `s*`: (a) in-sector gaps `g_k(s) ≥ g > 0` above each `φ_k(s)`; (b) weak inter-sector coupling, `|w| ≤ εg`, plus a transversality condition `Δ′(s*) ≠ 0` for `Δ := E₂ − E₀`. Then the lowest two `0⁺` eigenvalues and eigenvectors of `H(s)` are those of the 2×2 matrix `[[E₀(s), w(s)],[w̄(s), E₂(s)]]` up to errors `O(ε²g)` (Feshbach–Schur reduction with norm bounds), and the **order parameter**

> `η(s) := ⟨Ψ_gs(s), P₂ Ψ_gs(s)⟩ = ½ [ 1 − Δ(s)/√(Δ(s)² + 4|w(s)|²) ] + O(ε²)`

rises from `O(|w/Δ|²)` to `1 − O(|w/Δ|²)` across a window of width

> `|s − s*| ≲ 2|w(s*)| / |Δ′(s*)|`.

*Content:* "island-of-inversion onset" = diabatic crossing point `s*`; "sharpness of the inversion" = `|w|/|Δ′|`. No nonanalyticity at finite A — and that is the physically correct statement: shape coexistence with measurable mixing is what experiment sees (e.g., the low-lying `0₂⁺` in ³²Mg, Wimmer et al., *PRL* **105**, 252501 (2010); large `ρ²(E0)` values as direct mixing meters, Wood–Zganjar–De Coster–Heyde, *Nucl. Phys. A* **651**, 323 (1999)). [EXACT — provable now; the novelty is the packaging and the calibration of `(Δ, w, g)` against realistic valence-space Hamiltonians.]

---

## 3. Observable contact (what makes it falsifiable)

- **Binding-energy envelope kinks.** Along an isotopic chain, `E_gs(N)` follows (approximately) the lower envelope of branches; a branch switch produces a kink whose discrete second derivative is the two-neutron separation energy anomaly. This is literally how the N=20 island was *discovered* — anomalous extra binding in neutron-rich Na (Thibault et al., *Phys. Rev. C* **12**, 644 (1975)) — and how the N≈60 Zr–Sr transition announces itself in `S_2n`. The "bifurcation" in the data is the envelope corner, smoothed by `|w|`.
- **Islands as sublevel sets.** Define `Δ(N,Z) = E_intr − E_norm` on the even–even lattice; an island of inversion is a connected component of `{Δ < 0}` (Warburton–Becker–Brown, *Phys. Rev. C* **41**, 1147 (1990) for the original N=20 delineation). Its boundary is the experimentally mapped "shoreline." [SYNTH as a definition, but it matches usage.]
- **Mixing diagnostics:** `ρ²(E0)`, two-state-mixing fits of `0₁⁺/0₂⁺`, and the invariants `q², q³cos3γ` with their fluctuations (Poves–Nowacki–Alhassid 2020). The crossover width predicted by Theorem A is directly comparable to fitted mixing matrix elements — a quantitative, falsifiable bridge.

---

## 4. Where a true bifurcation lives — the classical-limit layer [DOABLE]

To make "bifurcation" literal, pass to a large-degeneracy limit where an exact classical variational problem emerges. The rigorous engine exists and is, pleasingly, partly due to Barry Simon himself:

- **Berezin–Lieb inequalities / coherent-state classical limits:** Lieb, *Commun. Math. Phys.* **31**, 327 (1973) for SU(2); **Simon, *Commun. Math. Phys.* 71, 247 (1980)** for general compact Lie groups. For Hamiltonians polynomial in the generators of a representation with index `N → ∞` (LMG: spin `N/2`; IBM: U(6) irreps with boson number `N_B`), the scaled ground energy converges to the minimum of the classical symbol over the coherent-state manifold.
- **Theorem B — Classical-Limit Bifurcation Theorem.** For a configuration-mixing algebraic model (two coupled sectors, e.g., IBM-CM with boson spaces `[N_B] ⊕ [N_B+2]`, Duval–Barrett, *Nucl. Phys. A* **376**, 213 (1982)), the `N_B → ∞` limit yields a **two-sheeted classical energy surface** `e_k(β,γ; s)`; the ground energy per boson converges to the minimum over both sheets; generically the transition is **first order** (Maxwell point where two distinct minima exchange depth — the spherical-vs-deformed coexistence scenario), with second-order/pitchfork points on critical lines (the U(5)–O(6) analogue). Finite-`N` behavior: at a first-order point the avoided-crossing gap is **exponentially small in N** (tunneling between minima), so `η(s) →` step function; at second-order points, power-law finite-size scaling (LMG: `N^{−1/3}` corrections, Dusuel–Vidal, *PRL* **93**, 237204 (2004); excited-state QPT structure, Caprio–Cejnar–Iachello, *Ann. Phys.* **323**, 1106 (2008); rigorous symmetry-breaking analyses of LMG/Curie–Weiss-type models by Landsman and collaborators, SciPost Phys. 2020). General organization: Cejnar–Jolie–Casten, *Rev. Mod. Phys.* **82**, 2155 (2010); catastrophe-theoretic classification of IBM surfaces, López-Moreno–Castaños, *Phys. Rev. C* **54** (1996).
- **Concrete open theorem [DOABLE, likely new]:** a Berezin–Lieb proof of the two-sheet limit and Maxwell-point first-order transition for IBM-CM specifically (the single-space IBM and LMG cases are essentially done; the *configuration-mixing* version — the exact mathematical home of your block matrix — appears not to have been done rigorously). This is the flagship target: rigorous, finite-page, and it is precisely "the island of inversion as a theorem."

So the corrected language is: **finite A = crossover with width `|w|/|Δ′|`; declared limit = genuine bifurcation of the lower envelope of classical sheets.** Both tiers are needed; reality (the data) sits at finite A, the sharp theorem sits in the limit, and the finite-size corrections connect them.

---

## 5. The driver layer — what makes Δ(s) move [EXACT identity + DOABLE criterion]

The monopole part of any rotationally invariant valence Hamiltonian is, exactly, a **quadratic polynomial in the number operators** `{n_j}` (Zuker's separation `H = H_mon + H_mult`, defined by reproducing all configuration centroids). Hence the effective single-particle energies are **affine in occupations**, `ε_j(n) = ε_j⁰ + Σ_{j'} V^mon_{jj'} n_{j'}`, and the diabatic gap obeys an explicit drift law along an isotopic chain. For a 2p–2h intruder:

> `Δ(N) = 2·Gap_ESPE(N) − [ E_corr^{intr}(N) − E_corr^{norm}(N) ]`,

with `Gap_ESPE(N)` affine in valence occupations, slopes given by monopole matrix elements whose signs are controlled by tensor (`j_> / j_<` rule) and 3N physics (Otsuka RMP 2020; Zuker 2003).

**Criterion C (monopole-drift crossing).** If the correlation differential is bounded and slowly varying on the window while the affine drift of `Gap_ESPE` is monotone with sufficient slope, then `Δ` changes sign at a unique `N*` in the window (intermediate value + monotonicity) — i.e., **existence and uniqueness of the inversion point, with the location expressed in terms of `V^mon`**. This converts the Otsuka/Poves phenomenology into inequalities. Caveat to state explicitly: **Type-II shell evolution** makes `Gap_ESPE` sector-dependent (the intruder's own occupations feed back on the monopole field), so the drift law is piecewise-affine across sectors; the variational structure of Section 2 survives untouched, only the explicit drift formula gains configuration labels.

---

## 6. Multi-sector reality: the lower envelope and its corners [SYNTH]

With sectors `k = 0, 2, 4, …` (np–nh), the structural skeleton is the **lower envelope** `𝔈(s) = min_k E_k(s)` of finitely many smooth branches: a piecewise-smooth function whose **corners are the structural transitions**, smoothed in the true spectrum by the relevant `w_{kk'}`. Catastrophe theory organizes the two-control picture (drift × correlation): folds = branch crossings, cusps = points where three branches compete (plausibly realized in the Zr chain, where IBM-CM analyses find an abrupt, QPT-like onset at N=60 — García-Ramos–Heyde, *Phys. Rev. C*, 2019/2020). The classical-limit version of this envelope is exactly the multi-sheet surface of Theorem B.

---

## 7. The pathology marker — why this structure is forced, not chosen [EXACT, underused]

There is a classical *rigorous* result that already detects your bifurcation point from the inside: the **Schucan–Weidenmüller theorem** (*Ann. Phys.* **73**, 108 (1972); **76**, 483 (1973)). When an intruder state dives into the model-space energy window, the perturbation expansion for the effective interaction (folded-diagram series) **provably diverges**; the radius of convergence is set by the complex-plane trajectory of the crossing. Translation: the diabatic crossing `s*` is precisely the point where the perturbative *definition* of the shell model breaks down — which is why islands of inversion are exactly where single-reference methods (and naive IMSRG decoupling) struggle and multi-reference machinery is required. Your theorem is therefore not an optional reframing: the mathematics already singles out configuration crossing as *the* obstruction. Incorporating Schucan–Weidenmüller as the "detection" clause makes the formulation airtight historically and mathematically.

---

## 8. The composed statement ("the shape closest to reality")

**Definition (configuration frame).** A configuration frame is a triple `(H(s), {P_k}, 𝒪)`: a smooth family of valence-space Hamiltonians, a fixed orthogonal sector decomposition commuting with symmetries, and a list of scheme-independent observables `𝒪 = {E(0_i⁺), E(2₁⁺), S_2n, ρ²(E0), B(E2), q²}`.

**Theorem (Shell evolution as configuration competition; two tiers).**
*Tier 1 (finite A).* Under the gap and weak-coupling hypotheses of Theorem A, the low-lying `0⁺` spectrum near a transversal diabatic crossing `s*` is a two-level anticrossing with explicit `O(ε²)` error; the sector character `η(s)` traverses a crossover of width `2|w|/|Δ′|`; the ground-energy curve follows the smooth lower envelope `min_k E_k` up to `O(|w|)`, producing `S_2n` anomalies at envelope corners; magic-number behavior at `N₀` along a chain holds iff the envelope remains on the normal branch (`Δ > 0`) throughout the chain window, and an island of inversion is a connected component of `{Δ < 0}` on the `(N,Z)` lattice.
*Tier 2 (declared limit).* In the large-degeneracy/classical limit of the corresponding algebraic configuration-mixing model, the scaled ground energy converges (Berezin–Lieb/coherent states; Lieb 1973, Simon 1980) to the minimum of a multi-sheet classical energy surface; the order parameter develops a genuine discontinuity at the Maxwell point (first order, generic) or a continuous pitchfork at critical points (second order), with finite-size gap exponentially small (first order) or power-law (second order) — the rigorous meaning of "bifurcation."
*Driver clause.* `Δ(s)` is governed by the exact monopole identity of Section 5; sign change criteria are expressible in `V^mon` (tensor + 3N), with Type-II feedback as a sector-dependent refinement.
*Detection clause.* `s*` coincides with the Schucan–Weidenmüller divergence point of the effective-interaction expansion: the bifurcation is where the perturbative shell model ceases to exist.

---

## 9. How this docks with the Strutinsky/cluster program

The two programs are not rivals; they are **driver and response**:
- The one-body/semiclassical program (clusters, Strutinsky, periodic orbits) is the theory of **`Gap(s)`** — where shell gaps come from and how big they are in the mean field.
- The configuration-competition theorem is the theory of **what the many-body system does with a given gap drift** — whether correlations convert it into inversion.
A complete answer to "when do magic numbers persist?" composes them: semiclassical/cluster input → monopole-dressed gap function → envelope/bifurcation output. Your refinement replaces the *wrong* conjecture ("clusters persist") with the *right* mechanism, while the cluster machinery survives as the microscopic source of `Δ`.

## 10. Ranked doable novelties

1. **IBM-CM Berezin–Lieb bifurcation theorem** (Tier 2 for the two-sheet model) — rigorous, self-contained, new. [DOABLE]
2. **Feshbach–Schur crossover theorem with constants calibrated to N=20 / Zr valence Hamiltonians** — modest depth, high contact with data (`w` vs fitted mixing, width vs measured `ρ²(E0)` systematics). [EXACT→DOABLE]
3. **Monopole-drift crossing criterion** as explicit inequalities in `V^mon`, with the Type-II piecewise refinement. [DOABLE]
4. **A quantitative Schucan–Weidenmüller statement for modern flows** (radius-of-convergence / decoupling-breakdown bound for IMSRG-type generators near a diabatic crossing) — harder, genuinely novel, and would tie the theorem to ab initio practice. [HARD]

*Flagged speculation:* items 1 and 4 as stated, the sublevel-set definition of islands (Sec. 3), and the envelope/catastrophe organization (Sec. 6) are this document's synthesis; everything in Sections 2, 5, 7 is standard mathematics assembled for this purpose, and all physics claims track the cited reviews.
