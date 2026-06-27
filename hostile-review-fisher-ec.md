# Hostile Review: "Fisher-Geometry Control of Eigenvector Continuation"
## Adversarial audit, and what survives in the ashes

**Reviewer stance:** maximally hostile, sources verified against primary literature this session. The brief is not to be fair; it is to find every way the section is wrong, imprecise, or non-novel, and then to state — honestly — the residue that withstands the assault.

**One-line verdict:** There is a correct and worthwhile theorem in here, but it is *not* the theorem stated, the proof of the headline subtlety is **backwards** relative to the established convergence theory it cites, and the framing oversells a relabeling as a mechanism. About 40% survives; the surviving 40% is genuinely useful and, suitably narrowed, publishable as a *geometric reformulation* — never as a new convergence criterion.

---

## 1. What the section claims (steelmanned)

Stripped to load-bearing assertions:
- **(C1)** EC is a reduced-basis/snapshot method, hence bounded below by the Kolmogorov m-width of the state manifold `ℳ`.
- **(C2)** The relevant geometry on constant-rank patches is Fisher/Bures; locally `d_B² = ¼ J_{ij} dλⁱdλʲ + O(‖dλ‖³)`.
- **(C3)** Therefore EC accuracy is "controlled by" the Bures-width compressibility of `ℳ`, and EC "degrades whenever the Fisher geometry broadens or becomes singular."
- **(C4)** At rank-changing / critical points the QFI is discontinuous while Bures is the continuous completion; basis size "generically increases."
- **(C5)** (Corollary) QFI enhancement/singularity ⇒ Bures width broadens ⇒ EC error worsens at fixed m; so the Fisher spectrum is "a practical diagnostic for EC breakdown."

I will grant C1 and the bare metric identity in C2, then dismantle the inferential spine C3–C5.

---

## 2. The kill shots

### Kill shot A — The metric used is the wrong metric for EC, and the established theory says so explicitly. [decisive]

EC's *proven* convergence theory is **not** governed by Bures/Fisher geometry. Sarkar–Lee (*PRL* **126**, 032501 (2021) [VERIFIED]) and the Drischler–Furnstahl–Melendez–Zhang emulator review (arXiv:2310.19419 [VERIFIED]) establish that EC convergence is governed by **Kolmogorov n-width in the Hilbert-space (ℓ²/vector) norm**, controlled by **analyticity of `|ψ(λ)⟩` and the distance to the nearest branch point in the complexified parameter plane**. The review states the mechanism in plain terms: the sharp bend in the Hubbard ground-state energy near `U/t = −3.8` is caused by an avoided crossing whose **branch points lie near the real axis**, and *that* — proximity of complex-plane singularities — is what sets EC error, with explicit bound `O(|θ/z|^{M+1})`, `z` the nearest non-analytic point.

The Bures metric and the vector-norm width are **inequivalent** as analytic control objects:
- Bures distance saturates: `d_B ≤ √2` always (fidelity ∈ [0,1]). It is **insensitive to the global phase and the embedding scale** that the EC Gram/norm-matrix analysis actually exploits.
- The EC subspace is spanned by *vectors* `|ψ(λ_a)⟩`, and its quality is a statement about the **secant/derivative structure of the vector-valued analytic curve**, not about the projective-state distinguishability metric. Two analytic curves with identical Bures geometry can have radically different vector-norm n-widths (e.g., differ by an analytic-in-λ phase `e^{iα(λ)}` that is invisible to `d_B` but changes the span and its Gram matrix). So C3 mis-identifies the controlling functional.

**The section commits a metric substitution that the field's own convergence proofs forbid.** This is not a quibble: the entire inferential chain C3–C5 is built on the claim that *Bures* width controls EC, and the proven theory says *vector-norm* width (driven by complex analyticity) does. Bures width is at best a lower bound proxy on constant-rank patches, and a loose one.

### Kill shot B — The criticality claim is provably backwards. [decisive]

C4–C5 assert that EC degrades at critical/level-crossing points because the QFI is enhanced/singular there. The established EC literature asserts — and proves by construction — the **opposite** in the most important case:

> Drischler–Furnstahl review [VERIFIED]: "the eigenvectors can still be defined as **analytic functions in the neighborhood of these exact level crossings**." And in the chemistry-EC literature (arXiv:2305.00060 [VERIFIED]): *"unless there is a **symmetry-protected** level crossing, the eigenvectors and eigenvalues along a given Hamiltonian path are continuous,"* and a **mixed ground/excited training set** makes EC accurate through avoided crossings and dissociation degeneracies precisely where fidelity-type diagnostics collapse.

Meanwhile, the fidelity/QFI quantity the section wants to use as a "breakdown diagnostic" has a known, published pathology at exactly the targeted points: **global-state fidelity drops to zero at every continuous level crossing** (Kwok–Ho–Gu, *PRA* **78**, 062302 (2008) [VERIFIED]), which is why the fidelity community had to invent *partial-state* fidelity. So:
- QFI/Bures is **maximally singular** at continuous level crossings, yet
- EC remains **analytic and well-behaved** there (with appropriate, standard training-set choices).

Therefore C5's "diagnostic" fires loudest exactly where EC does *not* break. The proposed criterion has false positives at the canonical critical points. A diagnostic that is maximally triggered where the method works is worse than no diagnostic.

The section's *one* defensible criticality instinct — that EC needs more basis vectors near **avoided** crossings with branch points near the real axis — is real, but it is (i) already in Sarkar–Lee/Drischler et al., and (ii) explained by **complex-plane analyticity**, not by Fisher singularity on the real axis. The section reaches a partially-right conclusion through demonstrably wrong reasoning, and generalizes it to the wrong (symmetry-protected/continuous-crossing) case.

### Kill shot C — The factor is written inverted; the convention is not innocuous here. [serious]

The section writes `d_B² = ¼ J_{ij} dλ dλ`. The verified standard relation is `4 g_B = F_Q` i.e. `d_B²= g_B dλdλ = ¼ F_Q dλdλ` (Šafránek 2017; Zhou–Jiang, arXiv:1910.08473 [VERIFIED]; "omitting a factor of four" caveats throughout the metrology literature). If the section's `J` denotes the **QFI** `F_Q` (as its name "quantum Fisher information tensor" states), then `d_B² = ¼ F_Q dλdλ` is **correct**. But Lemma 1 then conflates `J` with "the quantum Fisher metric" (= `g_B`, which is `¼ F_Q`), so the text uses `J` to mean both `F_Q` and `g_B` within three lines. With `J = g_B`, the formula `d_B² = ¼ J dλdλ` is **wrong by a factor of 4**. This is a genuine internal inconsistency, not pedantry: the whole section is *about* the Fisher/Bures identity, so getting the identity's bookkeeping ambiguous is a structural flaw at the exact point of maximal scrutiny.

### Kill shot D — "Reduced-basis ⇒ width lower bound" is true but vacuous as stated. [framing]

Lemma 2 (C1) is correct and is the one rigorous statement — but it is a **tautology imported from approximation theory**: *any* linear subspace method is bounded below by the Kolmogorov width of the target set in the relevant norm (Pinkus, *n-Widths in Approximation Theory*, 1985; Kolmogorov 1936). The section presents this as if Bures is the natural norm; it is not (Kill shot A). In the **correct** norm (vector ℓ²), this lemma is exactly Sarkar–Lee's starting point and adds nothing new. In the **Bures** norm, the lemma is true but disconnected from EC's actual error (which is measured by the generalized-eigenvalue/Rayleigh-quotient residual in vector norm, not by `inf_φ d_B`). The quantity `e_m(λ) = inf_{φ∈V_m} d_B(ρ(λ), |φ⟩⟨φ|)` is **not** the EC energy error and **not** the EC eigenvector residual; it is a projective-distance surrogate the section never connects to either EC output. The objects are defined, but the bridge to what EC reports is missing.

### Kill shot E — Constant-rank assumption excludes the entire phenomenon of interest. [scope]

The whole apparatus is restricted to `Λ₀` where `ρ(λ)` has **constant rank** and `J` is positive-definite. For ground-state EC the object is a **pure state** (`rank 1`) along the path — fine — but then the QFI/Bures machinery degenerates to the **Fubini–Study metric** on `ℂP^{d-1}` (verified: Bures = FS on pure states), and the elaborate density-operator framing is overkill: nothing mixed ever appears. Conversely, the *interesting* EC-stress events (true crossings, where the targeted eigenstate's character changes) are **exactly** the rank-change / degeneracy points the assumptions **exclude**. So the theorem is stated precisely on the set where nothing dramatic happens, and goes silent on the set it advertises (C4). The "Bures provides the continuous completion" gesture does not rescue this: continuity of the metric completion does not imply the EC subspace dimension behaves continuously, and no such implication is proved (the proof "sketch" asserts "typically requires a larger basis" with no argument).

---

## 3. Novelty audit

| Claim | Status | Prior art (verified) |
|---|---|---|
| EC = reduced-basis, width-bounded | **Not novel** | Sarkar–Lee 2021; Drischler et al. 2023 (explicitly "Kolmogorov N-width") |
| EC convergence is geometric/analytic | **Not novel** | Sarkar–Lee 2021 (analyticity, branch-point distance, "differential folding") |
| `4 g_B = F_Q`, discontinuity at rank change | **Not novel** (and is textbook) | Šafránek 2017; Zhou–Jiang 2019; Bures/Helstrom |
| QFI as criticality diagnostic | **Not novel**, and misapplied | Gu fidelity-susceptibility program 2008–2010; Kwok–Ho–Gu 2008 |
| *Bures*-width (not vector-width) controls EC | **Novel but false** (Kill shot A,B) | contradicted by Sarkar–Lee |
| RG as semigroup flow of distinguishability geometry | **Speculative, not new as metaphor** | information-geometric RG (Beny–Osborne 2015); not established as predictive |

Net: every *true* statement is prior art; the *novel* statement is the false one. That is the worst possible novelty profile.

---

## 4. What survives in the ashes

After the demolition, a non-empty and actually-useful residue remains. Stated honestly and narrowly:

**S1 — The width lower bound, in the correct norm.** EC is a reduced-basis method and is bounded below by the Kolmogorov m-width of `{|ψ(λ)⟩ : λ∈Λ₀}` **in Hilbert-space norm**. This is correct, standard, and the right skeleton. *Keep it; drop Bures from it.* [SURVIVES — but is Sarkar–Lee, not new.]

**S2 — A legitimate geometric reformulation on pure-state patches.** On the (generic, away-from-crossing) set where the targeted eigenstate is non-degenerate, the state curve lives on `ℂP^{d-1}` with the **Fubini–Study** metric, and EC's local secant geometry can be expressed in FS terms. This is a *coordinate-free restatement* of the analytic theory, and it has modest value: it makes the training-point optimization problem ("place snapshots to minimize sup FS-secant distance") geometrically natural, connecting to known greedy/width-optimal sampling. **This is the real, defensible contribution** — but it must be sold as *reformulation/geometric insight*, not as a new convergence-controlling mechanism, and it must use **Fubini–Study (pure-state)**, not the mixed-state Bures/QFI apparatus, which never engages. [SURVIVES, narrowed.]

**S3 — A corrected criticality statement.** The honest, defensible version of C4–C5: *EC subspace dimension must grow when `|ψ(λ)⟩` has branch points in the complex λ-plane lying close to the real training interval* — i.e., near **avoided** crossings with small gaps, where the analytic continuation radius shrinks. The controlling quantity is **distance to the nearest complex branch point**, not real-axis QFI singularity. (At **symmetry-protected/continuous** crossings, by contrast, EC stays accurate with a mixed training set, and QFI-type diagnostics *mislead*.) [SURVIVES only after inverting the mechanism and restricting to avoided crossings.]

**S4 — A genuinely open, correctly-posed question.** Is there a quantitative relation between (a) the real-axis behavior of the **fidelity susceptibility** `χ_F(λ)` (the leading QFI coefficient) and (b) the **imaginary part of the nearest branch point** `Im z(λ)` that actually controls EC? Heuristically both are governed by the gap to the nearest excited state, so a relation should exist — `χ_F ~ Σ |⟨n|∂_λ H|0⟩|²/(E_n−E_0)²` diverges as the real-axis gap closes, and the branch point pinches the real axis as the same gap closes. **Pinning this down rigorously** (when does fidelity-susceptibility growth predict EC basis-size growth, and when does it give false positives?) is a real, narrow, publishable theorem — and it would *correct* the present section rather than vindicate it. This is the version worth pursuing. [SURVIVES as a research target, reframed.]

**What does NOT survive:** the headline that *Bures/Fisher* geometry controls EC (S-A); the criterion that QFI singularity diagnoses EC breakdown (false-positive at continuous crossings, S-B); the mixed-state density-operator framing (never engages, S-E); the factor bookkeeping as written (S-C); and any suggestion that the result is new (S-novelty). The "RG as distinguishability-geometry semigroup" discussion is unobjectionable as motivation but carries zero proved content and should be labeled as speculation, not as a "mathematically natural object" that the analysis supports.

---

## 5. Minimal repair to make it correct and honest

1. **Replace Bures with Fubini–Study throughout**; delete the constant-rank density-operator scaffolding (it collapses to pure states anyway). Roughly halves the section and removes the factor ambiguity and the rank-change hand-wave in one move.
2. **Re-cite the actual theory**: Sarkar–Lee (*PRL* 126, 032501) and Drischler et al. (arXiv:2310.19419) own EC convergence; state your result as a *geometric reformulation* of theirs, not a replacement.
3. **Invert the criticality claim**: control object is `Im z(λ)` (nearest complex branch point), which shrinks at **avoided** (small-gap) crossings; explicitly note QFI gives **false positives** at symmetry-protected/continuous crossings, where EC is fine with mixed training.
4. **Demote the corollary** to a conjecture (S4) about `χ_F` vs `Im z`, clearly flagged as open, and either prove the easy direction (gap-closing ⇒ both grow) or drop the "practical diagnostic" language.
5. **Strike** the implied novelty; the honest contribution is "a Fubini–Study restatement of EC convergence that makes optimal-snapshot placement a width problem," which is a fine, modest paragraph — not a theorem environment with five numbered results.

---

## 6. Bottom line

The section dresses a **tautology** (reduced-basis ⇒ width bound) and a **textbook identity** (`4g_B = F_Q`) as a new theorem, attaches a **novel-but-false** claim (Bures width controls EC), and tops it with a **backwards** criticality criterion (QFI singularity ⇒ breakdown) that fires hardest exactly where the established theory proves EC keeps working. What remains after the fire is real but small: EC lives on a projective state-manifold, its convergence is a Fubini–Study width problem, and the honest open question is the precise link between real-axis fidelity-susceptibility and the complex branch-point distance that genuinely governs the method. Publish *that* — three honest sentences and one open conjecture — and retract the rest. The deliberately narrow disclaimer at the end ("not a gravity program, a spectral-geometry program") is wise; the section should extend the same discipline to its own central claim and concede that the controlling geometry is **analytic (complex-plane), realized metrically as Fubini–Study**, not **Fisher/Bures distinguishability** on the real axis.

---

### Verified source ledger
- S. L. Braunstein, C. M. Caves, *PRL* **72**, 3439 (1994) — QFI/Bures, Cramér–Rao [standard].
- D. Šafránek, *Phys. Rev. A* **95**, 052320 (2017); S. Zhou, L. Jiang, arXiv:1910.08473 (2019) — exact `4g_B=F_Q`, discontinuity at rank change [VERIFIED].
- A. Sarkar, D. Lee, *PRL* **126**, 032501 (2021) — EC convergence via analyticity / vector-norm width / differential folding [VERIFIED].
- C. Drischler, R. J. Furnstahl, J. A. Melendez, X. Zhang et al., "Eigenvector Continuation and Projection-Based Emulators," arXiv:2310.19419 (Rev. Mod. Phys.-track) — EC as reduced-basis, explicit Kolmogorov N-width, branch-point control, analyticity through crossings [VERIFIED].
- H.-M. Kwok, C.-S. Ho, S.-J. Gu, *PRA* **78**, 062302 (2008) — global fidelity collapses at every continuous level crossing; partial-state fidelity needed [VERIFIED].
- Chemistry EC, arXiv:2305.00060 (2023) — eigenvectors analytic unless symmetry-protected crossing; mixed training through avoided crossings/dissociation [VERIFIED].
- Pinkus, *n-Widths in Approximation Theory* (1985); Kolmogorov (1936) — width lower bounds for linear methods [standard].
