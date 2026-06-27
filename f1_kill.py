"""TRUE F1 KILL TEST: dressed-degenerate wells (codimension-1, engineered by tuning a linear
tilt delta so the two dressed sheet minima are EXACTLY equal -- repulsion detuning cancelled).
Pre-registered: (P1) the two anticrossings (A<->x_m in E1-E0, A<->x* in E2-E1, separated by
0.024 in s/M in the detuned gamma) MERGE to |Ds*/M| < 1e-3; (P2) gap-min exponent = min(S1,S2)
= S1. Falsification signatures: Gamma << S1 (cooperative enhancement) or near S2 / intermediate."""
import numpy as np
from pipeline import Vfun, quartic_from_roots, construct_invariant, build_dense, f64_seed
from scipy.optimize import brentq

r1g, r2g = 0.25, 0.65
def base(rb):
    pars, _ = quartic_from_roots(r1g, rb, r2g, 0.25, kapB=1.0)
    V = Vfun(pars); h = lambda x: V(x) - 4 * x * (1 - x)
    return pars, h
rbst = brentq(lambda rb: base(rb)[1](r1g) - base(rb)[1](r2g), r1g + 0.03, r2g - 0.03, xtol=1e-12)

def make(delta):
    pars, _ = base(rbst)
    p = dict(pars); p['eps'] = pars['eps'] + delta
    p.update(epsA=1.0, kapA=0.0, om=0.08)
    return p
def split(delta):
    res = construct_invariant(make(delta))[0]
    ms = sorted(res['mins'], key=lambda t: t[0])     # by position: [x_m-side, x*-side]
    return ms[0][1] - ms[1][1], res
dst = brentq(lambda d: split(d)[0], -0.08, 0.0, xtol=1e-12)
_, res = split(dst)
ms = sorted(res['mins'], key=lambda t: t[0])
acts = sorted(res['actions'], key=lambda t: t[0])
S1, S2 = acts[0][2], acts[1][2]
P = make(dst)
print(f'tuned tilt delta = {dst:.6f}: dressed minima {ms[0][1]:.7f} (x={ms[0][0]:.3f}) vs '
      f'{ms[1][1]:.7f} (x={ms[1][0]:.3f})  -> degenerate; F1 flag = {res["F1_flag"]}')
print(f'candidate actions: S1(0->{acts[0][0]:.3f}) = {S1:.5f},  S2(0->{acts[1][0]:.3f}) = {S2:.5f}')
print(f'PRE-REGISTERED: P1 anticrossings merge (|Ds*/M| < 1e-3, vs 0.024 detuned); '
      f'P2 Gamma = min = {S1:.5f}', flush=True)

N = 40; M = N + 2
H0, dA, dB = build_dense(N, P)
EA0 = np.linalg.eigvalsh(H0[:dA, :dA])[0]; EB0 = np.linalg.eigvalsh(H0[dA:, dA:])[0]
sc = EA0 - EB0; D = 3.0 + 0.06 * N
nB = np.array([2 * j for j in range(dB)], float) / M
def spect(s):
    H = H0.copy(); H[dA:, dA:] += s * np.eye(dB)
    E, U = np.linalg.eigh(H)
    info = [(float(np.sum(U[:dA, k] ** 2)),
             float(np.sum(U[dA:, k] ** 2 * nB) / max(np.sum(U[dA:, k] ** 2), 1e-300))) for k in (0, 1, 2)]
    return E[1] - E[0], E[2] - E[1], info
grid = np.linspace(sc - D, sc + D, 3001)
g01 = np.zeros(len(grid)); g12 = np.zeros(len(grid))
for i, s in enumerate(grid):
    g01[i], g12[i], _ = spect(s)
stars = {}
for tag, gv, idx in [('E1-E0', g01, 0), ('E2-E1', g12, 1)]:
    locm = [i for i in range(1, len(grid) - 1) if gv[i] < gv[i - 1] and gv[i] <= gv[i + 1]]
    best = min(locm, key=lambda i: gv[i])
    s_best = grid[best]; step = grid[1] - grid[0]
    for _ in range(14):
        g3 = np.linspace(s_best - 2.2 * step, s_best + 2.2 * step, 61)
        v3 = np.array([spect(s)[idx] for s in g3])
        s_best = g3[np.argmin(v3)]; step = g3[1] - g3[0]
        if step < 1e-12: break
    a, b, info = spect(s_best)
    gm = a if idx == 0 else b
    ks = (0, 1) if idx == 0 else (1, 2)
    stars[tag] = s_best / M
    print(f'  {tag}: deepest anticrossing s*/M = {s_best/M:.5f}, gap = {gm:.3e}; '
          f'states <x>_B = {info[ks[0]][1]:.3f}, {info[ks[1]][1]:.3f} (A-weights {info[ks[0]][0]:.2f}, {info[ks[1]][0]:.2f})')
print(f'  P1 verdict: |Ds*/M| = {abs(stars["E1-E0"]-stars["E2-E1"]):.5f}  (detuned gamma had 0.024)', flush=True)

Ns = np.arange(20, 53, 4).astype(float)
gaps = np.array([f64_seed(int(n), P)[1] for n in Ns])
ok = gaps > 1e-12
L = np.log(gaps[ok]) - np.log(np.sqrt((Ns[ok] + 1) * (Ns[ok] + 2)))
A = np.vstack([np.ones_like(Ns[ok]), Ns[ok]]).T
cg, *_ = np.linalg.lstsq(A, L, rcond=None)
print(f'  P2 verdict: MEASURED Gamma (f64, N=20-52) = {-cg[1]:.5f}  vs S1 = {S1:.5f}, S2 = {S2:.5f}')

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(6.2, 4.2))
ax.semilogy(grid / M, g01, lw=1.0, label='$E_1-E_0$')
ax.semilogy(grid / M, g12, lw=1.0, label='$E_2-E_1$')
ax.set_xlabel('$s/M$'); ax.set_ylabel('gap')
ax.set_title('Dressed-degenerate wells: merged 3-level anticrossing (N=40)', fontsize=10, fontweight='bold')
ax.legend(fontsize=8); plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/f1_kill_test.png', dpi=150, bbox_inches='tight')
print('figure saved.')
