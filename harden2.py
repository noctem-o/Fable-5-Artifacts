"""CORRECTED ROUND 2: real BETA (dip landscape, as designed) and real GAMMA (triple-degenerate).
Design-verification asserts added so a builder bug can never again hide upstream of both
prediction and measurement."""
import numpy as np, time
from pipeline import (Vfun, quartic_from_roots, construct_invariant, build_dense,
                      f64_seed, mp_gap, run_series, fit_two, NGRID)
from scipy.optimize import brentq

def verify_design(p, roots, tol=1e-9):
    V = Vfun(p); hB0 = lambda x: V(x) - 4 * p['kapB'] * x * (1 - x)
    d = 1e-6
    for r in roots:
        der = (hB0(r + d) - hB0(r - d)) / (2 * d)
        assert abs(der) < 1e-4, f'design FAIL: hB0\'({r}) = {der}'
    print(f'  design verified: hB0\'(roots) = 0 at {roots}; '
          f'depths f(x_m)={hB0(roots[0]):.4f}, f(barrier)={hB0(roots[1]):.4f}, f(x*)={hB0(roots[2]):.4f}')

print('=' * 78)
print('REAL BETA -- quartic dip landscape, sign bug fixed; design verified before anything else')
bpars, _ = quartic_from_roots(0.22, 0.40, 0.62, 0.40, kapB=1.0)
BETA = dict(bpars); BETA.update(epsA=1.0, kapA=0.0, om=0.08)
verify_design(BETA, [0.22, 0.40, 0.62])
res = construct_invariant(BETA, verbose=True)
Sb = res[0]['S_sel']
print(f'  PREDICTION (pre-registered): S = {Sb:.5f}  (dip is a spectator; single instanton through it)')
s0, gf, _ = f64_seed(64, BETA)
import mpmath
gm = mp_gap(64, BETA, Sb)
print(f'  T1 sanity N=64: f64 = {gf:.6e}  mp = {mpmath.mp.nstr(gm, 8)}  rel = {abs(float(gm)-gf)/gf:.2e}', flush=True)
resB = run_series(BETA, Sb, 'REAL BETA om=0.08')

print('\n' + '=' * 78)
print('REAL GAMMA -- triple-degenerate landscape (H3 violated), the genuine out-of-domain probe')
r1g, r2g = 0.25, 0.65
def depths(rb):
    pars, _ = quartic_from_roots(r1g, rb, r2g, 0.25, kapB=1.0)
    V = Vfun(pars); hB0 = lambda x: V(x) - 4 * x * (1 - x)
    return hB0(r1g) - hB0(r2g)
rbst = brentq(depths, r1g + 0.03, r2g - 0.03, xtol=1e-12)
gpars, _ = quartic_from_roots(r1g, rbst, r2g, 0.25, kapB=1.0)
GAMMA = dict(gpars); GAMMA.update(epsA=1.0, kapA=0.0, om=0.08)
verify_design(GAMMA, [r1g, rbst, r2g])
resG = construct_invariant(GAMMA, verbose=True)
sel = resG[0]
print(f"  BLIND CONSTRUCTION OUTPUT: selects x = {sel['x_sel']:.4f}, S_sel = {sel['S_sel']:.5f}, "
      f"F1 flag = {sel['F1_flag']}")
acts = {f"{a[0]:.3f}": a[2] for a in sel['actions']}
print(f'  candidate actions: {acts}')
print('  PRE-REGISTERED: gap-min exponent = S_sel (repulsion-detuned well selection); '
      'second anticrossing expected in E2-E1, displaced by ~M x (dressed-min splitting).')

N = 40; M = N + 2
H0, dA, dB = build_dense(N, GAMMA)
EA0 = np.linalg.eigvalsh(H0[:dA, :dA])[0]; EB0 = np.linalg.eigvalsh(H0[dA:, dA:])[0]
sc = EA0 - EB0; D = 3.0 + 0.06 * N
nB_int = np.array([2 * j for j in range(dB)], float) / M
def spect(s):
    H = H0.copy(); H[dA:, dA:] += s * np.eye(dB)
    E, U = np.linalg.eigh(H)
    info = []
    for k in (0, 1, 2):
        wA = float(np.sum(U[:dA, k] ** 2))
        xB = float(np.sum(U[dA:, k] ** 2 * nB_int) / max(np.sum(U[dA:, k] ** 2), 1e-300))
        info.append((wA, xB))
    return E[1] - E[0], E[2] - E[1], info
grid = np.linspace(sc - D, sc + D, 3001)
g01 = np.zeros(len(grid)); g12 = np.zeros(len(grid))
for i, s in enumerate(grid):
    g01[i], g12[i], _ = spect(s)
for tag, gv, idx in [('E1-E0', g01, 0), ('E2-E1', g12, 1)]:
    locm = [i for i in range(1, len(grid) - 1) if gv[i] < gv[i - 1] and gv[i] <= gv[i + 1]]
    print(f'  {tag}: {len(locm)} anticrossing(s) in window')
    for i in locm:
        s_best = grid[i]; step = grid[1] - grid[0]
        for _ in range(14):
            g3l = np.linspace(s_best - 2.2 * step, s_best + 2.2 * step, 61)
            v3 = np.array([spect(s)[idx] for s in g3l])
            s_best = g3l[np.argmin(v3)]; step = g3l[1] - g3l[0]
            if step < 1e-12: break
        a, b, info = spect(s_best)
        gmin = a if idx == 0 else b
        ks = (0, 1) if idx == 0 else (1, 2)
        print(f'    s* = {s_best:.5f} (s*/M = {s_best/M:.5f}): gap = {gmin:.3e}; '
              f'state{ks[0]} (A {info[ks[0]][0]:.2f}, <x>_B {info[ks[0]][1]:.3f}); '
              f'state{ks[1]} (A {info[ks[1]][0]:.2f}, <x>_B {info[ks[1]][1]:.3f})', flush=True)

Ns_g = np.arange(20, 53, 4).astype(float)
gaps = []
for Nv in Ns_g:
    _, gv, _ = f64_seed(int(Nv), GAMMA)
    gaps.append(gv)
gaps = np.array(gaps); ok = gaps > 1e-12
L = np.log(gaps[ok]) - np.log(np.sqrt((Ns_g[ok] + 1) * (Ns_g[ok] + 2)))
A = np.vstack([np.ones_like(Ns_g[ok]), Ns_g[ok]]).T
cg, *_ = np.linalg.lstsq(A, L, rcond=None)
print(f"\n  MEASURED (f64 window N=20-52, {ok.sum()} pts): Gamma = {-cg[1]:.5f}"
      f"   vs blind S_sel = {sel['S_sel']:.5f}"
      f"   vs alternatives {acts}", flush=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-whitegrid')
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
Nmid, loc, c1f, c2f = resB
axes[0].plot(Nmid, loc, 'o', ms=5, label='measured slopes (mp)')
xs2 = np.linspace(Nmid[0], 3000, 300)
axes[0].plot(xs2, c1f[0] + c1f[1] / xs2, '--', lw=1.2, label=f'$\\Gamma_\\infty$ = {c1f[0]:.4f}')
axes[0].axhline(Sb, ls=':', lw=1.6, color='k', label=f'pre-registered = {Sb:.4f}')
axes[0].set_xscale('log'); axes[0].set_xlabel('N'); axes[0].set_ylabel(r'$\Gamma_{\rm loc}(N)$')
axes[0].set_title('REAL BETA: dip landscape (bug fixed, design verified)', fontsize=10, fontweight='bold')
axes[0].legend(fontsize=8)
axes[1].semilogy(grid / M, g01, lw=1.0, label='$E_1-E_0$')
axes[1].semilogy(grid / M, g12, lw=1.0, label='$E_2-E_1$')
axes[1].set_xlabel('$s/M$'); axes[1].set_ylabel('gap')
axes[1].set_title('REAL GAMMA: spectroscopy at N=40 (triple-degenerate)', fontsize=10, fontweight='bold')
axes[1].legend(fontsize=8)
plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/mythos2_corrected.png', dpi=150, bbox_inches='tight')
print('figure saved.')
