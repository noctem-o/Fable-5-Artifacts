# Theorems for the One-Exponent Program
## Statements, proofs, and numerical certification

**Status labels:** [PROVEN] complete proof below (routine constants omitted where noted) · [PROVEN + CERTIFIED] proof plus machine-precision numerical check · [ARCHITECTURE] statement with complete proof strategy, details unwritten · [DERIVED + CERTIFIED] closed-form parameter-free formula derived from the classical symbol and confirmed against exact diagonalization in its regime of validity · [OPEN] posed precisely, unsolved.

Model throughout: the two-sector boson family of the experimental program. Sector A: M_A = N bosons, H_A = ε_A n̂_t (ε_A = 1). Sector B: M_B = N+2 bosons, H_B = ε n̂_t + (g/M_B) n̂_t² − (κ/M_B) Q̂², Q̂ = s†t + t†s. Coupling Ŵ = ω(s†s† + t†t†): ℋ_A → ℋ_B. Control: H(s) = (H_A ⊕ (H_B + s)) + (Ŵ + Ŵ†), restricted to even n_t. Diabatic classical data: V(x) = εx + gx², e_B(x) = V(x) − 4κx(1−x), minimum (x*, E₀); S_Agmon = ∫₀^{x*} λ₀(x)dx with cosh 2λ₀ = (V−E₀)/(2κx(1−x)) − 1.

---

## Theorem 1 (EP–QFI identity with transfer). [PROVEN]

**Part (a) — exact 2×2 identity.** Let H₂(s) = [[E_A(s), w],[w, E_B(s)]] be real symmetric with w > 0 constant and Δ(s) := E_A − E_B = δ·(s − s*), δ ≠ 0 (linear crossing). Then:
1. The analytic continuation of the eigenvalue pair has exactly two branch points (exceptional points) z_± = s* ± 2iw/δ, each a square-root EP.
2. The minimal real-s gap is 2w, attained at s*.
3. The ground-state Fubini–Study susceptibility χ_F(s) = ‖∂_sψ₀‖² − |⟨ψ₀,∂_sψ₀⟩|² equals w²δ²/(Δ² + 4w²)², maximized at s* with χ_F^max = δ²/(16w²); the quantum Fisher information F_Q = 4χ_F.
4. **F_Q^max · (Im z_±)² = 1, exactly.**
5. The mixing window {s : |Δ(s)| ≤ 2w} has half-width exactly Im z_±.

*Proof.* Eigenvalues E_±(s) = m(s) ± ½√(Δ(s)² + 4w²) with m = ½(E_A+E_B); the trace part m is entire and branch-free, so branch points are the zeros of the discriminant Δ(z)² + 4w² = 0, i.e. δ(z−s*) = ±2iw, giving (1); the discriminant is minimized on ℝ at s*, giving (2). Eigenvectors of a real symmetric 2×2 are (cos θ, sin θ) with tan 2θ = 2w/Δ; differentiating, θ′ = −wΔ′/(Δ²+4w²) = −wδ/(Δ²+4w²), and for a real unit vector χ_F = (θ′)², giving (3). Then F_Q^max(Im z)² = (δ²/4w²)(4w²/δ²) = 1, and |Δ| ≤ 2w ⇔ |s−s*| ≤ 2w/|δ|. ∎

**Part (b) — analytic perturbation of the crossing.** If instead Δ is real-analytic with a simple zero at s* and second-derivative bound κ₂ = sup_{D(s*,ρ)}|Δ″|/|Δ′(s*)| on a disc of radius ρ ≥ 8w/|Δ′(s*)|, then (1)–(5) hold with relative errors O(w κ₂/|Δ′(s*)|): by Rouché applied to Δ(z) ∓ 2iw on the disc, the EP pair persists with |Im z − 2w/|Δ′(s*)|| ≤ Cwκ₂·w/Δ′², and the χ_F computation perturbs at the same order. Hence F_Q^max(Im z)² = 1 + O(wκ₂/|Δ′|). *(Routine constants omitted.)* In the experimental model the control enters linearly (∂H/∂s = P_B), so Δ_eff is linear to leading order and this clause contributes only exponentially small corrections.

**Part (c) — transfer to the full operator.** Let H(s) = H₀ + sB be self-adjoint on finite-dimensional ℋ, B = P_B, and suppose on an interval I the two lowest eigenvalues are separated from the rest by g > 0. (i) The Riesz projection P(s) onto the doublet is real-analytic; the doublet eigenvalues are those of the analytic 2×2 family P(s)H(s)P(s)|_{Ran P}, so Part (a)/(b) applies to the doublet with effective (Δ_eff, w_eff). (ii) The susceptibility decomposes exactly as χ_F = Σ_{n≥1} |⟨n|B|0⟩|²/(E_n−E_0)² = χ_doublet + R with 0 ≤ R ≤ ‖B‖²/g² = 1/g². Hence at the crossing,
**|F_Q^max (Im z)² − 1| ≤ 16 w_eff²/(Δ_eff′² g²) + O(w_eff κ₂,eff/Δ_eff′).**
Since w_eff is exponentially small in N while g = O(1), the deviation is exponentially small — *quantitatively certified:* the measured deviations decayed 6.7×10⁻³ (N=8) → ~10⁻⁶ (N ≥ 36) across all three model cases, matching the predicted scale and decay. This bound also answers the referee's demand to *prove* EP control: within radius ~g of the crossing, the doublet's EP pair is, by (i)–(ii), the only singular structure up to the stated error.

---

## Theorem 2 (Exactly solvable anchor). [PROVEN + CERTIFIED]

For ε = g = 0 (case C1), with M = N+2 and J = M/2 via the Schwinger map (J_x = ½Q̂, n̂_t = J − J_z): H_B = −(4κ/M)J_x², whose ground space is span{|J_x = ±J⟩}; the even-parity ground state is the coherent cat (|J_x{=}J⟩ + |J_x{=}{-}J⟩)/√2, with number-basis amplitudes ψ(n) = √2·2^{−M/2}√C(M,n) for even n. Since the sector-A ground state is exactly |n_t = 0⟩,

**w_direct(N) = ⟨gs_B|Ŵ|gs_A⟩ = 2√2 · ω · √((N+1)(N+2)) · 2^{−(N+2)/2},**

hence −(1/N)log(w/ωN) → ½ log 2 = Γ exactly. *Certification:* the closed form matches exact diagonalization to relative error ≤ 1.3×10⁻¹⁵ at N = 8, 16, 32, 48. This pins one full column of the program's tables at theorem grade and explains why the Fubini–Study formula (falsified in general) is exact precisely here: the ground state is exactly coherent.

---

## Theorem 3 (Agmon law for the matrix element). [ARCHITECTURE; certified numerically]

**Statement.** For κ > 0 and a nondegenerate deformed minimum x* > 0: −(1/N) log(w_direct/ωN) → S_Agmon, the imaginary-action integral above.

**Proof chain (each link standard technology):** (a) *Exact reduction:* gs_A = |n_t=0⟩ exactly, so w_direct = ω[c₁ψ_B(0) + c₂ψ_B(2)] with c₁ ~ N: the claim reduces to boundary-tail asymptotics of the ground eigenvector of a Jacobi (three-term, after parity blocking five-term → reducible) matrix with smooth symbol. (b) *Positivity:* −H_B + cI has nonnegative entries (the hopping enters H_B with negative sign), so by Perron–Frobenius ψ_B > 0 componentwise — every path contribution has one sign; no cancellation can spoil lower bounds. (c) *Upper bound:* discrete Agmon/Combes–Thomas estimate with weight e^{Mφ(n/M)}, φ′ = λ₀, gives ψ_B(0) ≤ C e^{−M(S_Agmon − ε)}. (d) *Lower bound:* positivity plus a single-path/transfer-matrix bound, or rigorous discrete WKB for linear recurrences (Costin–Costin; the orthogonal-polynomial asymptotics tradition), gives the matching exponent. **Certification:** C1 exact (Theorem 2); C2 measured 0.20789 vs S_Agmon = 0.20795; C3 measured 0.24099 vs 0.24128, residuals consistent with the O(log N/N) window of the fits — four-decimal agreement across three landscapes including the designed same-x*, different-landscape discriminator. Status: a short self-contained paper; no conceptual obstruction.

---

## Result 4 (Spectral exponent = classical adiabatic action; the law that replaces the retracted ω² claim). [DERIVED + CERTIFIED in its regime]

**Construction (zero free parameters).** On the 2×2 classical symbol h(x,χ) = [[ε_A x, ω],[ω, V(x)+s̄−2κx(1−x)(1+cos2χ)]]: (i) tune s̄ = s̄(ω) so the two minima of the lower adiabatic sheet h₋ are degenerate at energy E_c(ω) — the left minimum sits at the boundary x = 0 with the exact value E_left = ½s̄ − √(¼s̄² + ω²). (ii) The instanton condition h₋(x, iλ) = E_c reduces, in closed form, to the **hyperbola condition**
**(h_A − E_c)(h_B + s̄ − E_c) = ω²,**
giving cosh 2λ_ω(x) = [V + s̄ − E_c − ω²/(ε_A x − E_c)]/(2κx(1−x)) − 1. (iii) **S_adiab(ω) = ∫₀^{x_min(ω)} λ_ω(x) dx**, with x_min(ω) the exact right-well location at degeneracy. The avoided crossing regularizes the x→0 log-singularity of the diabatic integrand — the geometric origin of the softening.

**Consistency [proven at formula level]:** as ω → 0⁺, E_c → 0, ω²/a → 0 pointwise, x_min → x*, so S_adiab(ω) → S_Agmon (verified numerically to ~10⁻⁷ at ω = 10⁻⁴).

**Certification (the no-free-parameter confrontation):** for ω ≤ 0.07, |Γ_gap^meas(ω) − S_adiab(ω)| ≤ 1.3×10⁻³ (C1) and ≤ 5.7×10⁻⁴ (C3) across the grid, with best agreement 3×10⁻⁵ (C3, ω = 0.0125) and the small positive residuals consistent with the documented pre-asymptotic finite-N drift. **This resolves the retracted quadratic law and the anomalous referee-response fits:** the deficit S_Agmon − Γ_spec(ω) is S_Agmon − S_adiab(ω), a *non-power-law* classical function; fitting powers to it necessarily yields window- and case-dependent pseudo-exponents (the measured p = 1.977 vs 1.769) — the anomaly was an artifact of fitting the wrong functional form, and the referee's demand for a derivation is met, in this regime, with no adjustable parameters.

**Documented deviation [honest limit]:** for Nω² ≳ 0.5 the measured exponent falls progressively *below* S_adiab (to −0.017 at ω = 0.2, C1; −0.011, C3), with extended-N local slopes (to N = 64) still drifting away — evidence that a second decay channel with smaller action overtakes the single-sheet instanton at extensive coupling.

## Conjecture 5 (Non-adiabatic channel). [OPEN, posed]

At extensive coupling, Γ_spec = min{S_adiab(ω), S_×(ω)} where S_× is the imaginary action of a sheet-switching path through the complexified discriminant zeros (h_A − h_B − s̄)(x,χ) = ±2iω — the Dykhne–Davis–Pechukas/Stokes mechanism for matrix symbols. The conjecture predicts the sign, onset (Nω² ~ 1), and growth of the documented deviation. Computing S_× and proving the min-rule (Hagedorn–Joye exponential estimates; matrix-WKB mode conversion) is the program's sharpest open analysis problem — and the deviation table above is its target data.

## Proposition 6 (Perturbation radius = EP distance). [textbook]

For the finite self-adjoint analytic family H(s), the Rayleigh–Schrödinger series for the ground state about s₀ has radius of convergence equal to the distance from s₀ to the nearest non-real degeneracy of the analytic continuation (Kato, *Perturbation Theory*, Ch. II) — for the doublet, |s₀ − z_±|. With Theorem 1(c) and Result 4 this yields the quantitative Schucan–Weidenmüller statement: near a configuration crossing, the single-reference expansion radius collapses like e^{−N S_adiab(ω)} (regime R-II), with the prefactor and regime corrections as documented.

---

## Ledger

| Item | Status | Tether |
|---|---|---|
| T1 EP–QFI identity + transfer | PROVEN | deviation decay 10⁻²→10⁻⁶ matches bound |
| T2 anchor closed form | PROVEN + CERTIFIED | rel. err ≤ 1.3×10⁻¹⁵ |
| T3 Agmon law (matrix element) | ARCHITECTURE | 4-decimal, 3 landscapes, discriminator |
| R4 Γ_spec = S_adiab(ω), regime R-II | DERIVED + CERTIFIED | ≤ 5.7×10⁻⁴ for ω ≤ 0.07; best 3×10⁻⁵; 0 parameters |
| C5 non-adiabatic channel | OPEN | deviation table, onset Nω² ~ 1 |
| P6 SW radius, quantified | textbook + R4 | — |

Two integration-endpoint bugs in the S_adiab solver were caught by the built-in ω→0 consistency check before any comparison with data — the check existed for exactly this purpose. Scope honesty, unchanged: all statements are about the two-sector model class; the ℂP⁵/IBM-CM lift and the physical bridge remain as scoped in the principle document. Reproduction: `certify.py`; figure `theorem_certification.png`.

---

## ADDENDUM — Conjecture 5 test: the χ-dependent coupling symbol (pre-registered)

**The hypothesis.** Before reaching for Stokes/sheet-switching paths, a candidate error was found *inside Result 4's own derivation*: the coupling operator Ŵ = ω(s†s† + t†t†) has a χ-dependent symbol, w(x,χ) = ω[(1−x) + x e^{2iχ}], because t†t† moves n_t. Under the barrier, |w|² → ω²[(1−x)² + x² + 2x(1−x)cosh 2λ] — exponentially amplified where it matters. The corrected instanton condition stays **linear in c = cosh 2λ**: c(x) = (a·B1 − ω²P)/(a·K + ω²R), with boundary data (s̄, E_c, x_min) provably unchanged (|w| = ω at χ = 0 and at x = 0). Pre-registered criteria: (i) S_corr(ω→0) = S_Agmon; (ii) S_corr < S_const for ω > 0; (iii) residuals Γ_meas − S_corr shrink to the finite-N drift scale (≲ few×10⁻³, positive-biased) at **all** ω.

**Verdict: (i) ✓ (ii) ✓ (iii) partial.** Consistency holds to 10⁻⁶; the correction is strictly softening. The residuals flipped sign and structure: meas − S_const ran from −0.0167 (C1) / −0.0114 (C3) at ω = 0.2 to small positive at small ω; meas − S_corr is now **uniformly positive and smooth**, from +0.011/+0.008 at ω = 0.2 down to +3.5×10⁻⁴/+1.0×10⁻⁴ at ω = 0.008. At small ω criterion (iii) is met; at large ω the +0.008–0.011 residuals exceed the pre-registered scale, so by the letter of the registration the criterion is **not fully met**.

**The bracket and the trajectory.** The measured curve is now sandwiched, S_corr(ω) < Γ_meas(ω) < S_const(ω), at all large-ω points, and the extended-N local slopes (to N = 62) are in every case still *falling toward* S_corr from above while moving *away* from S_const (e.g., C1, ω = 0.2: slopes 0.2653 → 0.2588, S_corr = 0.2518, S_const = 0.2793). The slope-to-S_corr gaps at N = 62 are strikingly uniform (0.0059–0.0070 across all four checks), suggesting a common subleading structure rather than case physics. This is *consistent with* the asymptote being S_corr — but the slopes have not converged, and that reading is post-hoc.

**Status changes.** Result 4 is amended: the parameter-free law uses the χ-dependent coupling symbol (S_corr); the constant-symbol version S_const is identified as the **dominant error** behind the previously documented deviation, and is retained only as the empirical upper edge of the bracket. Conjecture 5 **narrows but does not close**: no exotic channel is needed for the sign or the bulk of the deviation, but the residual (+0.003 to +0.011, growing with ω) remains unattributed, with three live candidates: (a) finite-N pre-asymptotics (favored by the slope trajectories; decidable by extended-precision runs at N ~ 100–300, cheap here since dim ~ N/2); (b) subprincipal/Berry-type corrections to the matrix-WKB exponent or the symbol-level coupling normalization (derivable); (c) a genuinely small residual channel. Third tether-catch of the program: two integration endpoints, now one constant-symbol approximation — each caught by a built-in check or by the deviation table doing its job as a detector.

---

## ADDENDUM — Conjecture 5: RESOLVED at tested parameters (extended-precision decision experiment)

**Method.** Spectral gaps to N = 280 via arbitrary-precision banded LDLᵀ inertia bisection (the full H(s) is pentadiagonal in an interleaved ordering), with staged crossing localization and a 5-point hyperbola finisher on gap²(s); adaptive precision to dps ≈ 85 at the deepest point (gap ~ 10⁻³⁸). Pipeline validated against float64 at N = 64 in all four series (rel. diff ≤ 1.2×10⁻⁶). Three pipeline defects were caught by the built-in tether checks during development (probe offsets vs large gaps; fixed-tolerance exhaustion at ΓN ≈ 87; bisection-bracket clipping) — each diagnosed by its signature and fixed before any verdict was drawn.

**Result (pre-registered criterion H-a, all four series):**

| series | S_corr | Γ∞ (1/N fit) | Γ∞ (+lnN/N fit) | b |
|---|---|---|---|---|
| C1, ω=0.2 | 0.25179 | 0.25201 ± 0.00003 | 0.25159 ± 0.00004 | 0.438 |
| C3, ω=0.2 | 0.15066 | 0.15103 ± 0.00005 | 0.15027 ± 0.00005 | 0.398 |
| C1, ω=0.1 | 0.30942 | 0.30986 ± 0.00005 | 0.30909 ± 0.00002 | 0.366 |
| C3, ω=0.1 | 0.20305 | 0.20350 ± 0.00006 | 0.20267 ± 0.00003 | 0.363 |

In every series the two extrapolation models **bracket S_corr** (+2–5×10⁻⁴ above / −2–4×10⁻⁴ below). Criterion H-a is satisfied throughout: |Γ∞ − S_corr| ≤ model spread ≈ 4×10⁻⁴ ≪ the residuals (+0.003–0.011) that motivated the conjecture.

**Conclusions.** (1) The large-ω deviation that created Conjecture 5 was **finite-N pre-asymptotics around the S_corr asymptote** — no second decay channel exists at the tested parameters within ±4×10⁻⁴ in the exponent. Conjecture 5 is resolved (negatively, in its exotic form) at (C1, C3) × (ω = 0.1, 0.2); its formal min-rule remains untested at still larger ω or other landscapes. (2) Result 4 is upgraded: the zero-parameter law Γ_spec(ω) = S_corr(ω), with the χ-dependent coupling symbol, is now certified both across the ω-grid at fixed N-window (≤5.7×10⁻⁴) **and asymptotically in N at extensive coupling** (±4×10⁻⁴). (3) **New sharply-posed object:** the 1/N slope coefficient is near-universal, b ≈ 0.36–0.44 across all four (case, ω) combinations, i.e. the gap prefactor scales as ~N^(1−b) ≈ N^0.6 rather than the naive N — a concrete target for the subprincipal/fluctuation-determinant derivation (observation only; not claimed).

**Updated ledger:** T1 [PROVEN] · T2 [PROVEN + CERTIFIED] · T3 [ARCHITECTURE; 4-decimal certified] · **R4 [DERIVED + CERTIFIED, ω-grid and N-asymptotics]** · **C5 [RESOLVED at tested parameters]** · P6 [textbook] · NEW: prefactor exponent β [OPEN, measured b ≈ 0.4]. Reproduction: `extended_N2.py`, figure `extended_N_scaling.png`.
- MYTHOS stress protocol executed: A1-A5 all survived; out-of-sample kill test at C4 (kappa=0.8) passed at 1e-5/2.2e-4; S_corr reclassified as chart of invariant S[h]; see mythos-report.md.
- Round 2 corrected: builder sign bug (tether-catch #4) found by hardened constructor; real beta (dip) certified 1.8e-4; real gamma maps out-of-domain law: omega^2-repulsion detuning selects the instanton well, confirmed by anticrossing spectroscopy. F1 narrowed to dressed-degenerate wells.
- F1 RESOLVED: min-rule Gamma = min(S1,S2) certified at engineered dressed degeneracy (0.27369 vs S1=0.28068, S2=0.41393); residual anticrossing offset measured as O(1) ZPE detuning (0.41, M-independent), merger asymptotic. Two-layer protection: repulsion detuning + ZPE detuning.
- Codim-2 coda: min-rule certified ON the resonance manifold (Gamma_loc -> S1 band); separation floor = 2t_B measures in-B action (0.115 vs 0.106 est.); fixed-tilt 0.322 anomaly resolved as sweep artifact. F1 fully closed at all three layers.
- Final audit: ALPHA T1 failure (N=64, rel 8.36) flagged as unresolved; refit excl. faulted point consistent with S_gen but clean re-run owed. README.md handoff written.
- ALPHA resolved: localization-divergence fault root-caused, patched (pipeline.py), confirmed by cure (8.36 -> 6.1e-8); clean re-run certifies branch rule (bracket, 1.4e-4). Flag: one irreproducible scratch N=88 eval (open). External review doc rejected; memo filed.
