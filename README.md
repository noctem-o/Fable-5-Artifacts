# One-Exponent Program — Index & Status Ledger (final handoff)
Central result: for two-sector bosonic configuration crossings, ALL exponentially small spectral scales
(gap, Im z_EP, (4 chi_F^max)^{-1/2}, SW radius) share ONE exponent S[H] = the instanton action of the
2x2 matrix symbol's lower sheet at double degeneracy, condition (h_A - E_c)(h_B - E_c) = |w(x,chi)|^2.
S_corr is a Darboux-chart evaluation of this invariant, not fundamental.

https://github.com/noctem-o/One-Exponent-Program_Fable-5/blob/main/one-exponent-manuscript.md

## Status ledger
| Claim | Status | Key evidence | Where |
|---|---|---|---|
| EP-QFI identity 4chi_F(Im z)^2 = 1 | PROVEN + certified 1e-6 | Theorem 1 | theorems-and-certification.md |
| C1 anchor w_direct closed form | PROVEN + certified 1.3e-15 | Theorem 2 | same |
| Agmon law for w_direct (err O(M^2/3)) | PROVEN | Theorem 3, elementary proof | theorem3-proof.md |
| Spectral law Gamma_inf = S[h] | DERIVED + CERTIFIED | C1/C3 <=4e-4; BLIND C4 (kappa=0.8): 1e-5, 2.2e-4; real BETA (dip): 1.8e-4 | mythos-report.md |
| Branch rule (kapA>0): a+b>0 = smaller root | DERIVED + supported* | ALPHA refit excl. faulted pt: brackets S_gen=0.09308 (+4.6e-4/-?) | this file, *caveat below |
| Min-rule Gamma = min(S1,S2): 3 layers | CERTIFIED | repulsion detuning; ZPE (0.41 measured); hybridized resonance manifold Gamma_loc -> S1 | mythos-report.md F1 sections |
| RG: dGamma/dlnN = -(Gamma - S), y=-1 | MEASURED | fixed line; b non-universal, flows with omega | mythos-report.md |
| Bonus measurables | MEASURED | ZPE offset 0.41 (M-indep); in-B action S_B: floor rate 0.115 vs 0.106 est. | codim2 coda |
*ALPHA caveat: T1 cross-check FAILED at N=64 (mp/f64 rel 8.36) -> outlier slope; cause unresolved.
Refit on remaining 8 slopes is consistent with S_gen, but a clean re-run with the fault diagnosed is OWED.

## Reproducibility
pipeline.py = canonical corrected machinery (builders, mp inertia-bisection gap, construct_invariant).
Drivers: mythos.py (round 1 + C4 kill), mythos2.py (round 2, annotated), harden2.py (corrected beta/gamma),
f1_kill.py (dressed degeneracy), codim2.py (quantum degeneracy). Figures alongside.

## Methodology rules (earned, in order of scar tissue)
1. Pre-register before measuring; brackets decide. 2. Design-verification asserts UPSTREAM of both
prediction and measurement (pred-vs-meas agreement cannot catch shared bugs). 3. Resonance manifolds
need co-moving tuning; fixed-parameter sweeps contaminate window fits. 4. Every T1 failure blocks
certification of its series -- no exceptions (see ALPHA).

## Open frontiers (each wants a fresh session)
F5 proof: spectral lower bound via Z2-gauged coupled-chain Riccati comparison. | b(omega) subprincipal
derivation (calibration: ZPE 0.41, S_B 0.115). | ALPHA T1 fault diagnosis + clean re-run. | Structural
escape: CP^2 / other coupling operators. | F2 boundary: omega -> omega_c washout mapping.

## Honest scope
"Certified" = numerical, pre-registered, bracketed -- not proof (Theorem 3 is the only proven pillar
of the spectral story's mechanism). Minimal model class; distance to Simon's Problem 11 remains large.

## UPDATE (final session): ALPHA caveat resolved
Fault root-caused and fixed: for large gaps the mp localization loop diverged (eps_s ~ delta^2 can
never meet <= 0.7*gap_hi; delta grows x12/iter; finisher samples the linear branch). Patch: skip mp
refinement when gapf > 1e-6 (f64 seed refined to 5e-13). Confirmed by cure: N=64 rel 8.36 -> 6.1e-8.
Clean full re-run: slopes smooth, reproduce original N>=88 chain to 5 decimals, fits 0.09383/0.09205
BRACKET S_gen = 0.09308, |mean-S| = 1.4e-4. ALPHA branch rule: DERIVED + CERTIFIED.
Residual flag (open): one isolated N=88 scratch evaluation disagreed with both full runs by 12.7%
(cause undetermined; does not touch the certified series). External review doc received and REJECTED
(fabrications incl. an invented fault diagnosis) -- see review-rejection-memo.md.
