# The IBM-CM Berezin–Lieb Bifurcation Theorem
## A complete research dossier: statement, proof architecture, verified toolkit, and reality calibration

**Status key:** [VERIFIED] = primary source checked this session · [STANDARD] = canonical, high-confidence · [DOABLE] = open theorem, existing technology suffices · [FRONTIER] = open, at the current edge of the field · [SYNTH] = this dossier's own construction.

---

## 0. The claim, and the novelty audit

**Claim.** The Interacting Boson Model with Configuration Mixing (IBM-CM) — the block-structured Hamiltonian on `ℋ_[N] ⊕ ℋ_[N+2]` that is the operational mathematical home of "normal vs. intruder competition" — admits a **rigorous** large-N classical limit in which the ground state undergoes a **genuine first-order bifurcation** (Maxwell point of a two-sheeted classical energy surface over `ℂP⁵`), with finite-N rounding controlled by **explicitly computable, exponentially small** coherent-state overlaps. Proving this would make "the island of inversion" — and by extension the *failure boundary of magic numbers* — a theorem.

**Novelty audit [VERIFIED this session].** The classical limit of the IBM is, in the entire physics literature, *heuristic*: coherent-state expectation values with `O(1/N)` contraction terms discarded (Dieperink–Scholten–Iachello, *PRL* **44**, 1747 (1980); Hatch–Levit, *PRC* **24**, 684 (1981); van Roosmalen's thesis; reviewed with the explicit "neglecting the O(N⁻¹) terms" step in Cejnar–Jolie–Casten, *RMP* **82**, 2155 (2010), §2.1). The configuration-mixing classical limit is likewise heuristic: the **matrix coherent-state method** (Frank–Castaños–Van Isacker–Padilla, AIP Conf. Proc. 638 (2002); Frank–Van Isacker–Vargas, *PRC* **69**, 034323 (2004); Frank–Van Isacker–Iachello, *PRC* **73**, 061302 (2006) — "Type II" QPTs). **No Berezin–Lieb-grade or Toeplitz-grade rigorous treatment of the IBM exists, let alone of IBM-CM.** Meanwhile, every analytic ingredient needed has matured independently in mathematical physics (Sections 2–3). The theorem sits in a genuine gap between two communities that have not met.

**The Simon irony, doubled.** This is an attack on Barry Simon's Problem 11, and its two tiers rest on two theorems of Barry Simon: the compact-group classical limit (*Commun. Math. Phys.* **71**, 247–276 (1980) [VERIFIED], extending Lieb, *CMP* **31**, 327 (1973)), which powers the bifurcation tier; and the "flea on the elephant" analysis of low-lying double-well eigenvalues (*J. Funct. Anal.* **63**, 123 (1985), with Jona-Lasinio–Martinelli–Scoppola, *CMP* **80**, 223 (1981)), which is precisely the mathematics of *which configuration a finite nucleus near the Maxwell point actually adopts* — the modern rigorous form being the Landsman-school SSB program (van de Ven–Groenenboom–Reuvers–Landsman, *SciPost Phys.* **8**, 022 (2020) [VERIFIED]).

---

## 1. Evidence hierarchy used throughout (as requested: prestige-ranked analysis)

All claims below are weighted by the following hierarchy, applied in order; institutional provenance is recorded but is subordinate to venue and replication.

| Tier | Source class | Instances in this dossier |
|---|---|---|
| **E1** | Evaluated databases & flagship refereed venues: ENSDF/AME (BNL-NNDC), *Annals of Math.*, *CMP*, *PRL*, *RMP*, Princeton UP monographs | Lieb '73, Simon '80/'85, Boutet de Monvel–Guillemin '81, Bordemann–Meinrenken–Schlichenmaier '94, DSI '80, Togashi et al. '16, the three RMPs ('05, '10, '11, '20), Zr evaluated level data |
| **E2** | Strong field journals & major-school monographs: *PRC*, *NPA*, *J. Funct. Anal.*, *J. Geom. Anal.*, *SciPost*, Cambridge UP, Springer/Birkhäuser | Iachello–Arima '87, Duval–Barrett '81/'82, Frank–Van Isacker–Vargas '04, Gavrielov–Leviatan–Iachello '22, Charles '03, Deleporte '19–'21, Ma–Marinescu '07, Le Floch '18, Landsman school '20 |
| **E3** | Preprints/proceedings from established groups | Deleporte–Le Floch (2025), recent exponentially-small-Toeplitz-eigenvalue work (2025–26), Leviatan–Gavrielov Bose-Fermi extension (2025) |
| Provenance map | Princeton (Lieb, Simon, BdM–Guillemin) · Caltech (Simon, later) · Yale–Tokyo (Iachello–Arima; Gavrielov–Leviatan–Iachello with Hebrew U.) · Tokyo/RIKEN/CNS (Otsuka school) · GANIL–UNAM (Van Isacker–Frank) · Ghent–Georgia Tech (Heyde–Wood) · Prague–Cologne–Yale (Cejnar–Jolie–Casten) · Paris-Saclay/Strasbourg/Sorbonne (Deleporte, Charles, Le Floch) · Radboud/Trento (Landsman, van de Ven) · UCLA (Biskup–Chayes–Starr) | — |

---

## 2. The model, precisely

**Hilbert space.** `ℋ_N^tot = ℋ_[N] ⊕ ℋ_[N+2]`, where `ℋ_[N]` is the symmetric irrep `[N]` of `U(6)` (one `s` and five `d` bosons; `dim ℋ_[N] = \binom{N+5}{5} ~ N⁵/120`). The semiclassical parameter is `ħ_eff = 1/N`. Geometrically `ℋ_[N] ≅ H⁰(ℂP⁵, 𝒪(N))` — holomorphic sections of the N-th power of the hyperplane bundle over the coadjoint orbit `ℂP⁵ = U(6)/(U(1)×U(5))`. This identification is what unlocks the modern toolkit.

**Hamiltonian (standard ECQF-CM form, as fitted to data by Gavrielov–Leviatan–Iachello, *PRC* **105**, 014305 (2022) [VERIFIED]):**

> `Ĥ = P̂_A [ ε_A n̂_d − κ_A Q̂^{χ_A}·Q̂^{χ_A} ] P̂_A  +  P̂_B [ ε_B n̂_d − κ_B Q̂^{χ_B}·Q̂^{χ_B} + κ'_B L̂·L̂ + Δ ] P̂_B  +  Ŵ`,
> `Ŵ = ω [ s†s† + (d†·d†)^{(0)}-type pair ] + h.c.`, mapping `[N] → [N+2]`,

with `Q̂^χ_μ = (s†d̃ + d†s)^{(2)}_μ + χ (d†d̃)^{(2)}_μ`. Sector A ("normal") lives in `[N]`; sector B ("intruder," microscopically 2p–2h proton excitation across Z=40 rebosonized as two extra bosons) lives in `[N+2]` (Duval–Barrett, *PLB* **100**, 223 (1981); *NPA* **376**, 213 (1982)).

**Coherent states and symbols.** `|N; β,γ⟩ = (N!)^{-1/2} (b_c†)^N |0⟩`, `b_c† = (1+β²)^{-1/2}[ s† + β cosγ d₀† + 2^{-1/2} β sinγ (d₂† + d₋₂†) ]`. Lower (covariant) symbol of the in-sector Hamiltonians, per boson, in the limit:

> `e(β,γ; ε,κ,χ) = ε β²/(1+β²) − κ [ 4β² − 4√(2/7) χ β³ cos3γ + (2/7) χ² β⁴ ] / (1+β²)²` (+ subleading one-body pieces),

the canonical IBM energy surface (Iachello–Arima, *The Interacting Boson Model*, Cambridge UP 1987; conventions as in Cejnar–Jolie–Casten). Upper and lower symbols of any fixed-degree polynomial in `u(6)` generators differ at relative `O(1/N)` — the Berezin–Lieb mechanism.

**The mixing symbol — exact and instructive [SYNTH, elementary].** For `Ŵ = ω_s s†s† + ω_d (d†·d̃†-pair) + h.c.`,

> `⟨N+2; β,γ | Ŵ | N; β,γ⟩ = √((N+1)(N+2)) · [ ω_s + ω_d β² ] / (1+β²)`,

γ-independent for the SO(6)-scalar pair; for `ω_s = ω_d = ω` it is exactly `ω√((N+1)(N+2))` — *shape-independent extensive mixing*. This reproduces (and slightly sharpens) the `Ω(β)` entries of the matrix coherent-state method, whose 2×2 (and 3×3, for [N]⊕[N+2]⊕[N+4]) potential-energy matrices

> `𝔈(β,γ) = [[ E_A(β,γ), Ω(β) ], [ Ω(β), E_B(β,γ)+Δ ]]`

are standard in applications (Frank–Van Isacker–Vargas '04; used with Gogny-mapped surfaces by Nomura–Rodríguez-Guzmán–Robledo and successors). **The theorem's job is to prove that the spectrum of `Ĥ` is governed by `𝔈` with controlled errors.**

**The two scaling regimes (the key structural decision):**
- **(R-I) Extensive mixing** — `ω` fixed, `Δ = Nδ`: the off-diagonal is `O(N)`, same order as the sheets. Classical limit = `min_{ℂP⁵} λ₋(h(β,γ))`, the lower eigenvalue of the intensive 2×2 symbol. The bifurcation lives in the geometry of `λ₋`.
- **(R-II) Physical scaling** — `ω, Δ` fixed in MeV (as in all fits), control parameter `s` (neutron number / interaction drift) moving the *intensive* sheet minima `μ_A(s), μ_B(s)` through a transversal crossing at `s*`: ground energy follows `N·min{μ_A, μ_B} + O(1)`; the crossover window is `O(1/N)`; the avoided-crossing gap is exponentially small (Theorem C below). Real Zr nuclei (N ≈ 3–8 bosons) sit here.

---

## 3. The verified toolkit (what exists, exactly)

1. **Berezin–Lieb inequalities & compact-group classical limit [E1, VERIFIED].** Lieb, *CMP* **31**, 327 (1973) (SU(2)); **Simon, *CMP* 71, 247 (1980)**: extension to *general compact Lie groups*, via coherent states on maximal-weight orbits, including the proof that every bounded operator has a diagonal (upper-symbol) representation — precisely the `U(6)` case needed. Refinement to coherent-state *matrix elements* of `e^{-βH}`: Biskup–Chayes–Starr (*CMP*, 2007), built to upgrade classical phase transitions to quantum ones at large spin — the same upgrade pattern needed here.
2. **Berezin–Toeplitz spectral theory on compact Kähler manifolds [E1/E2, VERIFIED].** Boutet de Monvel–Guillemin (*Ann. Math. Studies* 99, Princeton UP 1981); Bordemann–Meinrenken–Schlichenmaier (*CMP* **165**, 281 (1994)): `‖T_N(f)‖ → ‖f‖_∞`, symbol calculus with full `1/N` expansion; Charles (*CMP* **239**, 1 (2003)); Ma–Marinescu (Birkhäuser 2007); Le Floch (Springer 2018, the pedagogical entry point).
3. **Low-energy localization, wells, and exponential precision [E2/E3, VERIFIED].** Deleporte: *Low-energy spectrum of Toeplitz operators* (thesis, Strasbourg 2019) — Morse **and Morse–Bott** well asymptotics, "quantum selection" when subprincipal terms discriminate degenerate classical minima; miniwell paper (*CMP*); *Toeplitz operators with analytic symbols* (*J. Geom. Anal.* **31**, 3915 (2021)) and Charles (arXiv:1912.06819): Bergman kernel and symbol calculus with `O(e^{-cN})` errors on real-analytic Kähler manifolds (ℂP⁵ qualifies); Agmon-type *forbidden-region decay* for Toeplitz eigenfunctions (Deleporte, arXiv:2001.07921). Active 2025–26 frontier: distribution of exponentially small Toeplitz eigenvalues on projective manifolds; Deleporte–Le Floch non-self-adjoint Bohr–Sommerfeld (arXiv:2504.00965); Deleporte–Vũ Ngọc uniform well asymptotics. **The community that proves exactly these lemmas is active right now.**
4. **Strict deformation quantization & rigorous SSB in mean-field models [E2, VERIFIED].** van de Ven–Groenenboom–Reuvers–Landsman, *SciPost Phys.* **8**, 022 (2020): finite-N Curie–Weiss ↔ discretized double-well Schrödinger operator; "SSB is tied to a limit but must occur, approximately, before the limit"; van de Ven, *J. Math. Phys.* **61**, 121901 (2020) (classical limit of mean-field theories); Landsman–Moretti–van de Ven, *Rev. Math. Phys.* **32** (2020); continuing CMP/LMP papers (Drago–van de Ven 2024–25). This is the C*-algebraic chassis for Theorem B's lower bound and for the finite-N "which-configuration" selection (the flea mechanism).
5. **Finite-size scaling at algebraic-model QPTs [STANDARD].** Dusuel–Vidal, *PRL* **93**, 237204 (2004): `N^{1/3}` scaling at the LMG critical point via continuous unitary transformations; excited-state QPT structure: Caprio–Cejnar–Iachello, *Ann. Phys.* **323**, 1106 (2008). These are the exponents the rigorous Toeplitz analysis must (and, by Deleporte's degenerate-well results, does) reproduce.
6. **The phenomenological target, quantified [E1/E2, VERIFIED].** Togashi–Tsunoda–Otsuka–Shimizu, *PRL* **117**, 172502 (2016): MCSM identifies the Zr N=60 shape change as a QPT driven by Type-II shell evolution (self-reinforcing `νg_{7/2} ↔ πg_{9/2}` monopole feedback). Gavrielov–Leviatan–Iachello, *PRC* **105**, 014305 (2022): full IBM-CM fit of ⁹²⁻¹¹⁰Zr against levels, S₂n, E2, E0, isotope shifts, magnetic moments; "intertwined QPTs" = abrupt Type-II configuration crossing + gradual Type-I shape evolution *within* the intruder sheet (U(5)→SU(3)→O(6)). Experimental sharpness: `E(2₁⁺)` falls ~1.2 MeV (⁹⁸Zr) → ~0.2 MeV (¹⁰⁰Zr); `B(E2)` rises ~×10 across two neutrons; kinks in S₂n and charge radii (evaluated data, NNDC; lifetime program ⁹⁸⁻¹⁰⁴Zr 2024).

---

## 4. The theorem suite [SYNTH — the dossier's construction]

**Theorem A (single-sector classical limit; the missing "rigorization of 1980").** *Let `Ĥ_N` be a fixed self-adjoint polynomial of bounded degree in the `u(6)` generators on `[N]`, intensively normalized. Then `E_0(N)/N → min_{ℂP⁵} e(β,γ)` and `−(βN)^{-1} log Tr e^{-βĤ_N} → min`-free-energy analogue, with rate `O(1/N)`; moreover the low-lying spectrum near a nondegenerate minimum is `e_min + N^{-1}(harmonic frequencies) + O(N^{-2})`.*
*Proof inputs:* Lieb–Simon upper/lower symbols (existence of diagonal representation: Simon '80, Appendix 2); BMS '94 for the expansion; Deleporte wells for the low-lying refinement. **Status: [DOABLE], essentially an assembly paper — and, astonishingly, absent from the literature as a stated theorem.**

**Theorem B (IBM-CM matrix-symbol bifurcation; regime R-I).** *With `Δ_N = Nδ`, `ω` fixed: `E_0(N)/N → min_{ℂP⁵} λ₋(h(β,γ;s))` where `h` is the 2×2 intensive symbol matrix; the minimizer data `(β*, γ*, θ*)` (with `θ*` the configuration mixing angle of the minimizing eigenvector) is the classical order parameter. Along a generic one-parameter family `s`, the ground-state phase diagram has fold lines and first-order (Maxwell) lines where two global minima of `λ₋` exchange depth; across a Maxwell line the limit order parameters `(β*, θ*)` jump.*
*Proof inputs (the genuinely new mathematics):* (i) **inter-irrep Toeplitz calculus**: `Ŵ : H⁰(ℂP⁵,𝒪(N)) → H⁰(ℂP⁵,𝒪(N+2))` is a Toeplitz-type operator whose symbol is a section of `𝒪(2)`; the needed composition/norm asymptotics live inside Ma–Marinescu's bundle-valued framework — assembly required, no conceptual obstruction; (ii) **lower bound**: either a matrix-valued upper-symbol decomposition (Simon-'80-style diagonal representation, extended to the 2-sector algebra) or the strict-deformation-quantization bundle route (extend van de Ven JMP '20 to the continuous field of C*-algebras over `1/N` with fiber `C(ℂP⁵) ⊗ M₂(ℂ)` at `0`); (iii) upper bound: matrix coherent-state trial states — already the physics method. **Status: [DOABLE→FRONTIER boundary]; this is the flagship.**

**Theorem C (physical regime R-II: envelope, `O(1/N)` rounding, and the exact tunneling exponent).** *With `ω, Δ` fixed and intensive minima `μ_A(s), μ_B(s)+0` crossing transversally at `s*` with distinct shapes `β_A* ≠ β_B*`: (a) `E_0(N,s) = N min{μ_A(s), μ_B(s)} + O(1)` uniformly near `s*`; (b) the configuration weight `η_N(s) = ⟨Ψ_0, P_B Ψ_0⟩` converges to the step function, with transition window of width `O(1/N)`; (c) the minimal avoided-crossing gap obeys*

> `log Δ_gap(N) = log(ωN) − N·Γ + O(1)`,  `Γ = ½ log[ (1+β_A*²)(1+β_B*²) / (1+β_A*β_B* cos(γ_A*−γ_B*)·…)² ]`,

*i.e., the tunneling exponent is the explicit logarithm of the coherent-state overlap per boson; for a spherical-to-deformed crossing (`β_A*=0`): `Γ = ½ log(1+β_B*²)`.*
*Why this is easier than standard double-well tunneling:* the two "wells" live in **different sectors coupled at first order by `Ŵ`** — the splitting is a *perturbative matrix element with an exponentially small overlap factor*, not an instanton through a barrier. Leading log-asymptotics need only: sector ground states exponentially close to coherent+Bogoliubov states (Deleporte's analytic calculus), in-sector gap bounds (Theorem A refinement), and Feshbach–Schur. **Status: (a),(b) [DOABLE]; (c) at log-scale [DOABLE], full prefactor [FRONTIER] — and it is the directly falsifiable output: it predicts the ⁰₂⁺–0₁⁺ mixing scale from fitted `(ω, N, β*)` alone.**

**Theorem D (critical scaling at intra-sheet Type-I points).** *At a quartic (spherical-to-deformed critical) point of a sheet, the low-lying gap scales as `N^{-1/3}` with explicit eigenfunction concentration; along O(6)/γ-unstable directions the minima are Morse–Bott circles and the Deleporte Morse–Bott well analysis applies verbatim.* Matches the Dusuel–Vidal exponent; upgrades "intertwined QPTs" to rigorous intra-sheet critical asymptotics. **Status: [DOABLE by adaptation].**

**Corollary E (observables; the reality interface).** Lower-envelope kink ⇒ S₂n discontinuity of slope `|μ_A' − μ_B'|·N`; `E(0₂⁺)(s)` has a minimum at `s*` equal to `Δ_gap(N)`; `⟨Q·Q⟩` (the rotational invariant — *not* `⟨Q⟩`, which vanishes by Wigner–Eckart) jumps across the window; `ρ²(E0)` peaks `∝` mixing² inside the `O(1/N)` window. Each maps onto a measured Zr signature (Section 3, item 6).

---

## 5. Proof architecture — the lemma chain

1. **L1 (Geometry/quantization dictionary).** `[N] ≅ H⁰(ℂP⁵, 𝒪(N))`; polynomial `u(6)` operators = Berezin–Toeplitz operators with polynomial symbols; upper/lower symbol difference `O(1/N)` with constants from BMS '94. [STANDARD]
2. **L2 (Single-sheet wells).** Nondegenerate minima of `e(β,γ)` ⇒ low-lying clusters at `N e_min + harmonic ladder + O(1/N)`; Morse–Bott version for γ-unstable. (Deleporte.) [DOABLE]
3. **L3 (In-sector gap).** Away from Type-I critical points, the sector-restricted gap above the sector ground state is bounded below by `c/N`·(harmonic frequency) — uniform in `s` on compacts avoiding criticality. [DOABLE]
4. **L4 (Inter-irrep symbol calculus).** Norm and composition asymptotics for Toeplitz maps `𝒪(N) → 𝒪(N+2)` with section-valued symbols; `‖Ŵ‖ = ωN(1+o(1))·max-symbol`; exponential off-diagonal kernel decay in the analytic category. [DOABLE — assembly within Ma–Marinescu/Charles/Deleporte technology; the one lemma with no exact off-the-shelf statement]
5. **L5 (Coherent overlap estimates).** `⟨N;β',γ'|N;β,γ⟩ = [ (1 + β'β cos Θ)/√((1+β'²)(1+β²)) ]^N` exactly (single-boson overlap to the N-th power); sector ground states are `O(e^{-cN})`-close to squeezed coherent states ⇒ the matrix element `⟨gs_B|Ŵ|gs_A⟩ = C ωN e^{-NΓ}(1+o(1))` at log-scale. [DOABLE]
6. **L6 (Feshbach–Schur reduction).** With L3 + L5: two-level reduction with `O(‖W_eff‖²/gap)` errors ⇒ Theorem C(b),(c). [STANDARD]
7. **L7 (Matrix Berezin–Lieb / SDQ lower bound).** Either extend Simon-'80's diagonal-representation theorem to the direct-sum module with `𝒪(2)`-twisted off-diagonal symbols, or build the continuous C*-bundle with fiber `C(ℂP⁵)⊗M₂` and prove the ground-state-sequence convergence à la van de Ven. ⇒ Theorem B. [FRONTIER-adjacent; the heart of the new contribution]
8. **L8 (Maxwell-point selection — the flea).** Near degeneracy of the two sheets, exponentially small terms (the `Ŵ` overlap, or symmetry-breaking perturbations) select the finite-N ground state; rigorous template: Jona-Lasinio–Martinelli–Scoppola '81, Simon '85, modernized in SciPost 8, 022 (2020). Explains the empirically delicate, "in-between" character of ⁹⁸Zr. [DOABLE by transfer]

**Failure modes to watch:** (i) L4 constants could degrade near `β → ∞` chart boundaries — work in the compact ℂP⁵ picture throughout; (ii) at the *intertwined* point where the Type-II crossing coincides with a Type-I critical point (the Zr empirical situation!), L3 fails — the uniform analysis there is genuinely [FRONTIER] and genuinely interesting (a "critical Maxwell point"); (iii) angular-momentum projection: all statements above are for energies/order parameters at `O(1)`–`O(N)` resolution; rotational band fine structure is `O(1/N)` and needs equivariant (SO(3)-reduced) refinements — defer.

---

## 6. Reality calibration — the Zr chain as the laboratory

The mapping from theorem to nature, with verified anchors:

- **Control parameter:** neutron number 52→70 ⇒ drift of `(ε, κ, χ, Δ)` (the GLI-fitted family). The Type-II crossing sits between N=58 and N=60 — exactly where `E(2₁⁺)` collapses ~1.2→~0.2 MeV and `B(E2)` jumps ~×10 [VERIFIED].
- **Boson numbers:** `N = 3–8` along the chain — small! This is a *feature*: the theorem's finite-N corrections (`O(1/N)` window, `e^{-NΓ}` gap) are not academic decorations; they are the difference between the idealized step and the measured rounding, and small `N` makes them measurable. The Curie–Weiss↔double-well correspondence (SciPost '20) shows rigorous statements remain meaningful and checkable at `N ~ 60`, let alone `N ~ 8` with exact diagonalization as ground truth.
- **Falsifiable output #1 (gap formula):** `Δ_gap ≈ C ωN(1+β_B*²)^{-N/2}` vs. the fitted mixing/`0₂⁺` systematics of the GLI Zr fit and the measured `ρ²(E0)` values. A log-scale match across the chain would be the first *rigorous-formula-level* contact between Kähler quantization and nuclear data. [SYNTH]
- **Falsifiable output #2 (abruptness gradient):** the theorem says sharpness is controlled by `Γ(β_A*, β_B*)` and the window by `1/N|μ'_A − μ'_B|`. Empirically the N=60 transition is sharpest at Z=40 and softens toward Mo/Kr/Ru [VERIFIED]. Conjecture: the softening tracks the convergence of the two sheets' shapes (`Γ → small`) and slope difference — extractable from existing IBM-CM fits of the neighboring chains. A clean cross-chain test. [SYNTH]
- **Microscopic bridge (mitigating "it's only an effective model"):** (a) MCSM/Type-II shell evolution gives the configuration-resolved microscopy the IBM-CM coarse-grains (Togashi '16); (b) EDF-mapped IBM-CM (Gogny surfaces → IBM-CM matrix, Nomura et al. line of work) replaces fitted parameters with mean-field-derived ones, tightening the chain `chiral/EDF → IBM-CM → theorem`. The theorem rigorizes the last arrow; the first arrow remains physics — state this honestly, always.

---

## 7. Attack plan (ordered, with deliverables)

1. **Paper 0 — "The classical limit of the interacting boson model, rigorously."** Theorem A + L1–L3. Self-contained; closes a 45-year-old gap between DSI '80 and Simon '80, who published the two halves *in the same year* without ever being joined. Low risk, immediate citation surface in two communities.
2. **Paper 1 (flagship) — "Configuration mixing as a two-sheeted classical limit: envelope, rounding, and exponentially small gaps."** Theorem C + L4–L6 + Corollary E, calibrated on the GLI Zr parameters; includes the explicit `Γ` formula and the cross-chain abruptness conjecture. This is "the island of inversion as a theorem," regime R-II.
3. **Paper 2 — "Matrix Berezin–Lieb and the bifurcation theorem."** Theorem B via L7 (SDQ-bundle route), with the fold/Maxwell classification of the `λ₋` landscape (connecting to the catastrophe-theoretic IBM analyses). Natural collaboration surface: Landsman school (SDQ), Deleporte/Charles/Le Floch (Toeplitz), Leviatan–Gavrielov (model & data).
4. **Paper 3 — "Critical Maxwell points."** Theorem D + the intertwined case where Type-I criticality meets the Type-II crossing — the genuinely new analysis problem this program donates back to pure semiclassics.
5. **Stretch — prefactor-level tunneling** for the inter-sector gap (full Helffer–Sjöstrand-grade asymptotics in the Kähler category): the one item parked at today's [FRONTIER], with the 2025–26 exponentially-small-Toeplitz-eigenvalue literature converging on it from the mathematics side.

---

## 8. Why this is the highest-probability route to "mapping reality"

Every alternative rigorous attack on shell evolution either lives too far from data (abstract many-body limits with no fitted counterpart) or too far from theorems (ab initio computation with uncertainty bands but no proofs). The IBM-CM bifurcation theorem is unique in that **the identical mathematical object** — the 2×2 matrix surface over `ℂP⁵` — is simultaneously (i) a Toeplitz-quantization problem the current mathematical community can finish, (ii) the working tool already fitted, parameter by parameter, to the best-measured shape-transition chain in the nuclear chart, and (iii) the exact formalization of the configuration-competition structure identified as "closest to reality" in the previous stage of this program. The theorem does not approximate the phenomenology; it *is* the phenomenology, promoted to mathematics.

**Honesty box.** [SYNTH] items: Theorems A–D as stated, the `Γ` formula and both falsifiable outputs, the "critical Maxwell point" problem, and the two-communities-gap thesis. Everything in Sections 2–3 is verified or canonical as labeled. The bridge from the nucleon Hamiltonian to IBM-CM remains outside the theorem's scope — this program rigorizes the effective theory's bifurcation, which is the correct, honest target for Simon Problem 11's "make mathematical sense of" clause at the configuration-competition layer.
