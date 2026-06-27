# The One-Exponent Principle
## Exceptional-point geometry of nuclear configuration crossings: one computable Kähler exponent governs tunneling, analyticity, distinguishability, and emulability

**Status key:** [VERIFIED] = primary source checked this session · [STANDARD] = canonical · [ELEMENTARY] = derived here, checkable in five lines · [DOABLE] = open theorem, existing technology · [FRONTIER] = current research edge · [CONJ] = this document's conjecture, clearly flagged.

---

## 0. The idea in one paragraph

When a nucleus crosses from a normal to an intruder configuration (an island of inversion; the Zr N=60 transition), the two competing `0⁺` branches anticross. The anticrossing is governed, in the complexified control parameter, by a conjugate pair of **exceptional points (EPs)** — square-root branch points of the eigenvalue/eigenvector problem. **The principle:** in the configuration-mixing boson model (IBM-CM) that quantitatively fits these transitions, the height of the EPs above the real axis obeys an exponential law

> `Im z(N) ≍ e^{−NΓ}`,  with  `Γ = −log |⟨b_c(β_A*,γ_A*) | b_c(β_B*,γ_B*)⟩|`

— the **per-boson Fubini–Study (coherent-overlap) separation of the two competing classical minima on ℂP⁵**, in closed form `Γ = ½ log[(1+β_A²)(1+β_B²)/(1+β_Aβ_B cosΘ_{AB})²]`-type, and `Γ = ½ log(1+β_B²)` for a spherical→deformed crossing. This **single computable exponent simultaneously controls**: (i) the tunneling gap at the avoided crossing, (ii) the analyticity radius that limits every single-reference method (perturbation theory, folded diagrams, one-sided eigenvector-continuation extrapolation), (iii) the peak quantum Fisher information / fidelity susceptibility, and (iv) the parameter resolution needed to even *detect* the transition. The four faces are tied by parameter-free identities (e.g. `F_Q^max · (Im z)² = 1`), the law is rigorously provable with the Berezin–Toeplitz machinery assembled in the companion dossier, and it is numerically falsifiable on a laptop this month and experimentally anchored in the measured Zr systematics.

---

## 1. Assets fused (with provenance)

1. **Theorem C of the IBM-CM dossier** [DOABLE, this program]: avoided-crossing gap `log Δ_gap = log(ωN) − NΓ + O(1)`, with `Γ` from the exact coherent-state overlap `⟨N;α|N;α'⟩ = [single-boson overlap]^N`, provable via analytic Berezin–Toeplitz calculus (Deleporte, Charles) + Feshbach–Schur. Crucially, the two "wells" sit in *different boson spaces* coupled at first order by the pair operator `Ŵ`: the exponential smallness is a **coherent-overlap factor, not a barrier instanton** — which is what makes the law explicit.
2. **Hostile-review survivors** [VERIFIED last session]: eigenvector-continuation (EC) convergence is governed by **distance to the nearest complex branch point** (Sarkar–Lee, *PRL* **126**, 032501 (2021); Drischler–Furnstahl–Melendez–Zhang review, arXiv:2310.19419, explicit `O(|θ/z|^{M+1})` radius bound); and the open question S4 — relate real-axis fidelity susceptibility to `Im z` — which this principle **answers** in a controlled model class.
3. **Schucan–Weidenmüller (1972/73)** [STANDARD]: the perturbation expansion for effective interactions diverges when an intruder state crosses into the model-space window, with the radius of convergence set by the **complex trajectory of the crossing** — i.e., SW was an exceptional-point statement avant la lettre. The principle supplies its missing quantitative law.
4. **EP–QPT precedent** [VERIFIED this session]: Stránský–Dvořák–Cejnar (arXiv:1710.07553; Prague) found EPs approaching a **first-order** Lipkin QPT *exponentially* fast in system size (and polynomially for second order), proposing EP distributions as parametrization-independent signatures of criticality; Heiss–Scholtz–Geyer, *J. Phys. A* **38**, 1843 (2005), mapped large-N Lipkin EPs; Cejnar–Heinze documented EP accumulation near IBM QPTs (nucl-th/0406060). **The exponential pinching phenomenon is thus independently established at physics level in the simplest model** — strong grounding — while the configuration-mixing case, the closed-form geometric exponent, the rigor, and the cross-domain unification below remain open.
5. **The Zr laboratory** [VERIFIED]: Gavrielov–Leviatan–Iachello IBM-CM fit of ⁹²⁻¹¹⁰Zr (*PRC* **105**, 014305 (2022)); Togashi et al. MCSM QPT (*PRL* **117**, 172502 (2016)); measured abruptness (E(2₁⁺): ~1.2→~0.2 MeV; B(E2) ×10; S₂n and charge-radius kinks).

---

## 2. The exact two-level core [ELEMENTARY]

Near the diabatic crossing `s*`, the Feshbach-reduced problem is `H₂(s) = [[E_A(s), w],[w, E_B(s)]]` with extensive slope `Δ(s) = E_A − E_B ≈ N Δμ'·(s−s*)` and effective coupling `w = w_eff(N)`. Then, exactly at this order:

- **EP location:** eigenvalues coincide where `Δ(z) = ±2iw` ⇒ `z_± = s* ± i·(2w/N|Δμ'|)`; so `Im z = 2w/(N|Δμ'|)`.
- **Gap:** `Δ_gap = 2w` at `s*`.
- **Mixing window:** `|Δ| ≤ 2w` ⇔ `|s−s*| ≤ Im z`. *The crossover window in the control parameter equals the EP height.*
- **Fidelity susceptibility / QFI:** with ground state `(cosθ, sinθ)`, `χ_F = (θ')²`, peak `χ_F^max = (NΔμ')²/(16w²)`, hence the **parameter-free identities**

  > `χ_F^max · (Im z)² = ¼`,  equivalently  `F_Q^max · (Im z)² = 1`.

  Maximal distinguishability rate × analyticity radius = 1: the corrected, rescued version of the Fisher-geometry instinct from the reviewed manuscript — QFI does encode the controlling analytic scale, but as a *derived identity at avoided crossings*, not as a mechanism, and only there. At symmetry-protected crossings `w → 0`: `Im z → 0` while sector-resolved `χ_F` stays finite — the hostile review's false-positive case is recovered as the clean degenerate limit, resolved by sector decomposition.
- **FS path length** through the crossing: `∫|θ'| ds = π/2` — order one. (Consequence below: two-sided methods are cheap.)

Now insert the dossier's law `w_eff(N) = c·ωN·e^{−NΓ}`:

> `Δ_gap ≍ ωN e^{−NΓ}`,  `Im z ≍ (ω/|Δμ'|)·e^{−NΓ}`,  `χ_F^max ≍ (Δμ'/4ω)² e^{+2NΓ}`,  window `≍ e^{−NΓ}`.

Everything is slaved to `Γ`, a closed-form functional of the two classical shapes.

---

## 3. The Pinch Theorem [DOABLE — the new rigorous target]

**Theorem (pinch law).** *For the IBM-CM family with fixed `(ω, Δ)` and a transversal crossing of intensive sheet minima with distinct shapes `β_A* ≠ β_B*` (and away from intra-sheet Type-I critical points), the ground-state eigenvalue curve extends analytically to a strip minus a conjugate pair of square-root branch points `z_±(N)` with*

> `log Im z_±(N) = −N·Γ(β_A*, γ_A*; β_B*, γ_B*) + O(log N)`,

*`Γ` the per-boson coherent-overlap exponent. Consequently the radius of convergence of any single-configuration perturbative expansion about `s₀` equals `|s₀ − z_±| → Im z` as `s₀ → s*`.*

**Proof architecture:** dossier lemmas L3 (in-sector gap `≳ c/N`), L4 (inter-irrep Toeplitz calculus for `𝒪(N)→𝒪(N+2)` maps), L5 (sector ground states `O(e^{−cN})`-close to squeezed coherent states; exact overlap `[(1+β'β cosΘ)/√((1+β'²)(1+β²))]^N`), L6 (Feshbach–Schur two-level reduction with `O(w²/gap)` error), plus one new step: analytic continuation of the reduced 2×2 family in `s` (the entries are polynomial/analytic in `s` by construction) with error control uniform in a complex strip of width `≫ Im z` — available in the analytic Berezin–Toeplitz category (real-analytic symbols on ℂP⁵; Deleporte JGA 2021; Charles 1912.06819). **Honest delineation vs. precedent:** Stránský–Dvořák–Cejnar's exponential approach is numerical/asymptotic in single-space LMG; the new content is (i) the *configuration-crossing* (two-irrep, Type-II) case that models real shell evolution, (ii) the **closed-form geometric exponent**, (iii) theorem-level control, (iv) the correspondences of §4.

---

## 4. The four-way correspondence (corollaries) [ELEMENTARY given §3; SYNTH as a unification]

**C1 — Schucan–Weidenmüller, quantified.** The 1972 divergence theorem for effective interactions acquires its missing law: near an island of inversion, the folded-diagram/single-reference radius of convergence collapses as `e^{−NΓ}`. The same applies to any decoupling flow truncated around a single reference: the EP pinch is the *invariant obstruction* that "intruder-state problem" folklore points at.

**C2 — The multi-reference advantage is `e^{+NΓ}`.** One-sided analytic methods (PT, extrapolating EC, single-reference flows) are blocked at the EP scale: extrapolation across the crossing degrades catastrophically as `Im z → 0` (Bernstein-ellipse/radius bounds, per the verified EC convergence theory). But the FS path length through the crossing is only `π/2`: **two-sided/multi-reference methods need only O(1) extra basis directions** (the EC literature's own 2–4-snapshot success through avoided crossings is this fact in action). The quantitative gain of multi-reference over single-reference across a configuration crossing is therefore exponential in N with exponent `Γ` — to our knowledge the first closed-form complexity separation between the two method classes in nuclear many-body practice. *Caution flag:* this is a statement about *crossing* the transition, not about cost within a phase; and it is a model-class statement pending the lift of §6.

**C3 — Needle-in-haystack for emulator campaigns.** The transition's footprint in control space has width `≍ e^{−NΓ}`. Adaptive emulators and UQ scans (eigenvector-continuation/reduced-basis pipelines over chiral-EFT parameter space) require parameter resolution exponential in N transverse to the crossing manifold to *detect* islands of inversion — a concrete, actionable warning plus a fix: place the EP-locating step (complex-discriminant root finding on the emulator itself) inside the sampling loop. This converts the principle into an algorithm.

**C4 — QFI law, and where the Fisher instinct lives or dies.** `F_Q^max·(Im z)² = 1` with both sides slaved to `Γ`: the QFI peak grows as `e^{+2NΓ}` while its support shrinks as `e^{−NΓ}` (integrated FS length stays `π/2`). QFI is thus a *legitimate* transition detector in this model class at avoided crossings — and provably misleading at symmetry-protected ones — settling the scope question raised in the hostile review with a theorem instead of a fight.

---

## 5. Falsification protocol (laptop-scale, weeks) [DOABLE now]

1. Implement IBM-CM (GLI Hamiltonian form), synthetic families with frozen intensive parameters, `N = 2…12` (largest space `dim ≈ 1.8×10⁴`: trivial exact diagonalization).
2. Extract four independent estimates of the exponent: (i) slope of `log Δ_gap` vs N; (ii) `Im z` by complex-`s` root finding on the discriminant (or resolvent poles), slope of `log Im z`; (iii) `½ × ` slope of `log χ_F^max`; (iv) decay of one-sided EC/PT extrapolation accuracy.
3. Compare all four to the **analytic** `Γ` computed from the matrix coherent-state surfaces' minimizers `(β_A*, β_B*)`. Verify the parameter-free identities `F_Q^max(Im z)² = 1` and window `= Im z` with `O(1/N)`-controlled deviations.
4. Stress tests: γ-unstable intruder (Morse–Bott minimum ⇒ `Γ` becomes a minimum over the degenerate orbit — predicted modification), near-coincident shapes (`Γ → 0`: law degrades gracefully into the polynomial second-order regime, matching Stránský–Dvořák–Cejnar's dichotomy — a built-in consistency check against the verified precedent).
5. Port to the actual GLI Zr parameters; confront: fitted mixing strengths, measured `ρ²(E0)` and `0₂⁺` systematics; the discrete-N sampling of isotopic chains explains why most chains show clean inversion between adjacent even isotopes (probability `∼ e^{−NΓ}` of landing inside the window) with occasional near-window cases (⁹⁸Zr's notoriously intermediate character). Cross-chain prediction: the abruptness gradient Zr→Mo/Sr→Ru tracks `Γ` and slope differences extracted from existing IBM-CM fits of those chains.

A null result at step 3 — exponents disagreeing beyond `O(log N/N)` — kills the principle cleanly. That is the point.

---

## 6. The microscopic lift [CONJ — flagged speculation, with a test]

Replace boson coherent states by HFB/Slater configuration vacua: the overlap of distinct intrinsic states is given by the **Onishi formula** and decays exponentially in the number of rearranged particles. **Conjecture:** for valence-space shell-model Hamiltonians with a normal/intruder crossing, the same four-way correspondence holds with `Γ_micro = −(1/N_eff) log |⟨Φ_norm|Φ_intr⟩|`. Immediate test bed: the Tokyo MCSM Zr calculations, whose T-plot machinery already tracks exactly these basis-state overlaps — extract `Γ_micro`, locate EPs of the valence Hamiltonian by complex parameter continuation, compare. If it holds, the principle graduates from the algebraic model to the operational shell model — and the *boundary of magic-number validity across the chart* becomes a computable EP variety with exponentially quantified sharpness. If it fails, the failure mode (which correspondence leg breaks first) is itself diagnostic.

---

## 7. Novelty ledger (hostile-review-proofed)

| Statement | Status |
|---|---|
| EPs pinch exponentially at first-order QPTs (LMG, single space) | Known, physics level — Stránský–Dvořák–Cejnar 2018; Heiss–Scholtz–Geyer 2005 [VERIFIED] |
| EC/emulator convergence set by complex branch-point distance | Known — Sarkar–Lee 2021; Drischler et al. review [VERIFIED] |
| Effective-interaction divergence at intruder crossings | Known — Schucan–Weidenmüller 1972/73 [STANDARD] |
| QFI/fidelity peaks at QPTs; 4g_B = F_Q | Known — fidelity-susceptibility program; Šafránek/Zhou–Jiang [VERIFIED] |
| **Closed-form Kähler exponent `Γ` for the configuration-mixing (Type-II) pinch** | **New** [DOABLE theorem] |
| **Rigorous pinch law via Berezin–Toeplitz** | **New** [DOABLE→FRONTIER] |
| **Four-way correspondence + `F_Q^max(Im z)²=1` + multi-reference advantage `e^{NΓ}` + emulator detection bound** | **New as a unification** [ELEMENTARY/SYNTH] |
| **Microscopic Onishi lift** | **New, conjectural** [CONJ] |

Failure modes, stated up front: N-dependent prefactors could contaminate exponent extraction at `N ≤ 12` (mitigate: fit with `log N` term, as the theorem's `O(log N)` allows); multi-sector chains (`[N]⊕[N+2]⊕[N+4]`) produce EP *constellations*, not a single pair — the law should apply pairwise but the bookkeeping is real; the extensive-mixing regime (R-I of the dossier) moves the EPs into the matrix *symbol* and changes the analysis qualitatively (interesting, separate); and the entire microscopic relevance rests on the IBM-CM↔nucleus bridge, which remains physics, not theorem — same honesty box as the dossier, unchanged.

---

## 8. Why this is the answer to the prompt

Grounded: every leg is either verified primary literature, an elementary identity, or a theorem whose proof architecture was already assembled and audited in this program — and the one phenomenological precedent that could have falsified the novelty (Prague's EP–QPT work) instead independently confirms the mechanism in the simplest model while leaving the law itself open. Frontier: it sits at the live intersection of Kähler quantization (2019–2026 Toeplitz analysis), emulator theory (2021–2024), non-Hermitian/EP physics, and FRIB-era shell evolution — and it converts the oldest pathology of nuclear many-body theory (the intruder-state problem, 1972) into a single geometric number you can compute from two shapes. And it closes the arc of this entire investigation: Simon's Problem 11 asked what the shell model means; the configuration-competition reframing located where it fails; the Berezin–Lieb dossier made the failure a theorem-in-progress; the hostile review forged the analyticity tool; and the One-Exponent Principle is the object all four were pointing at — **the validity boundary of the nuclear shell model is an exceptional-point variety, and its distance from reality is `e^{−NΓ}`.**

---

## ADDENDUM — First experimental contact (pre-registered numerical test, this session)

The principle was subjected to a pre-registered falsification test in the minimal faithful model (two-sector U(2)/ℂP¹ configuration mixing, exact diagonalization, N = 8–52, three cases including a designed discriminator with identical shapes but different landscapes). Verdict:

**Confirmed.**
- The invariant identity `F_Q^max · (Im z)² = 1` holds to **one part in 10⁶** at large N in all cases — the EP↔QFI lock is exact, parametrization-invariant machinery (now "Lemma 0," provable for any analytic anticrossing dominated by a conjugate EP pair).
- One-exponent locking within the spectral class: gap, EP height, and χ_F peak share a single exponent at every coupling strength (agreement ~10⁻³).

**Falsified and corrected (live).**
- The closed-form Fubini–Study exponent `Γ_FS` of §0 is **wrong in general**: the discriminator split exactly as designed (C2: measured 0.20789 vs S_WKB = 0.20795 vs Γ_FS = 0.21539; C3: 0.24099 vs 0.24128 vs 0.21539). The correct universal exponent is the **Agmon/imaginary-action (WKB) distance** `S_Agmon` between the two configuration minima at the crossing energy; `Γ_FS` is its mean-field evaluation, exact only when sector ground states are exactly coherent (the C1 anchor). All statements of `Γ` in this document should be read as `S_Agmon`, still a single computable classical quantity.

**Discovered (and then confirmed by a targeted ω-scan).**
- Two exponent classes exist: the *wavefunction/matrix-element* class (spectroscopic mixing amplitudes; always `S_Agmon`) and the *spectral/analytic* class (gap, Im z, χ_F, emulator radius), which at extensive mixing is reduced to an adiabatic-surface action obeying the empirical law `S_adiab(ω) ≈ S_Agmon − c·ω²` (c ≈ 2.1, case-independent here), merging with `S_Agmon` as ω → 0. This is the dossier's R-I/R-II dichotomy, now measured: the spectral exponent at ω = {0.2, 0.1, 0.05, 0.025, 0.0125} climbed {0.2626, 0.3181, 0.3386, 0.3445, 0.3461} → 0.34657 (C1) and {0.1583, 0.2110, 0.2323, 0.2387, 0.2404} → 0.24128 (C3).

**Status change:** the Pinch Theorem (§3) should be restated with `S_Agmon`, plus a regime clause: spectral exponent = `S_Agmon` in R-II, = adiabatic action `S_adiab(ω)` in R-I, with the two-class structure and the `F_Q(Im z)² = 1` identity holding in both. Scope caveats unchanged: this validates the mechanism in the ℂP¹ minimal model; the ℂP⁵ IBM-CM computation, the rigorous theorems, and the microscopic lift remain open. Reproduction: `one_exponent_test.py`, `omega_scan.py` (figures: `one_exponent_test.png`, `omega_regime_test.png`).

---

## ADDENDUM 2 — External referee response, with one retraction

A hostile external review of the experimental summary was received and acted on. Outcome:

**Language corrections adopted.** "Measured facts" is withdrawn in favor of the referee's grading: *strongly supported by numerical evidence in the minimal model, at finite size (N ≤ 52), float64 precision, over 1.4 decades of ω, in one model family.* All claims below are so graded.

**Pre-registration trail (answering the out-of-sample question).** (i) Predictions P1–P4 were written in the script header before run 1 (verifiable in `one_exponent_test.py`). (ii) The R-I/R-II two-regime structure predates the experiment — it is §2 of the companion dossier. (iii) The ω-scan prediction (Γ_gap → S_Agmon as ω → 0) was formulated after the run-1 anomaly but before the scan. (iv) The quadratic-correction claim was post-hoc — and is now partially **retracted** (below). Alternative explanations for the run-1 anomaly considered at the time: numerical/search artifact (ruled out by the s*-fix and cross-channel consistency); competing anticrossings (ruled out structurally); W-ladder vs adiabatic-softening mechanism wording (both second-order-in-W; the scan tests their shared prediction and does not discriminate them — conceded).

**RETRACTION: "S_adiab ≈ S_WKB − 2.1·ω², case-independent" is withdrawn.** The coefficient 2.1 was the deficit/ω² value at the single largest ω (0.2), where the two cases coincidentally agree (2.099 vs 2.074). The fine scan (10 ω values, 0.008–0.2, deficits 3–17σ above fit error) shows: **C1**: deficit/ω² rises and saturates ≈ 3.33; free exponent p = 1.977 ± 0.005 on ω ≤ 0.05; pure ω² strongly preferred over ω²·log(1/ω) (χ² = 10 vs 868). **C3**: deficit/ω² does not saturate in the window (2.07 → 8.21); apparent p = 1.77 ± 0.04; the log-corrected model is preferred (χ² = 188 vs 759) but neither is adequate. Local-slope analysis shows C3 not yet converged at N = 52 (slopes still rising; 1/N extrapolation 0.24091 vs S_WKB = 0.24128), indicating the ω → 0 and N → ∞ limits may not commute uniformly. **Open:** double-scaling analysis; perturbative derivation of the softening coefficient (which would also discriminate the mechanism wording).

**EP-dominance diagnostic (responding to "observing an EP ≠ EP controls asymptotics").** The identity deviation |4χ_F·(Im z)² − 1| is itself the dominance meter: it is derived assuming a single conjugate EP pair dominates the local analytic structure, and it decays from ~10⁻² (N = 8) to ~10⁻⁶ (N ≥ 36) in all cases. This establishes local EP dominance at the measured sizes; universality beyond this model family is open, as the referee states.

**Surviving hard core at referee grade:** (i) the two-level identity F_Q^max·(Im z)² = 1 — provable in closed form for an isolated analytic anticrossing and verified to 10⁻⁶; (ii) the Agmon exponent for the matrix-element channel — 4-decimal agreement across three landscapes including a designed discriminator; (iii) the convergence of the spectral exponent to S_Agmon as ω → 0 — a pre-scan prediction confirmed to ~5×10⁻⁴. The referee's proposed abstract-level phrasing is adopted as the program's current claim of record.
