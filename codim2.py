"""CODIMENSION-2 PROBE: tune tilt so the two anticrossings MERGE at N=40 (quantum degeneracy:
classical + ZPE detuning both cancelled). Pre-registered: even at the genuine 3-level point,
Gamma(gap-min) = S1 = min over candidate actions (Perron-Frobenius positivity of both couplings
forbids destructive interference in the ground sector). Falsification: Gamma < S1 (cooperation)
or > S1 (suppression). Note the tuned point is N-dependent (offset ~ 0.41/M), so the window
N=20-52 sits within O(0.2 extensive) of merger throughout -- stated, not hidden."""
import numpy as np
from pipeline import Vfun, quartic_from_roots, construct_invariant, build_dense, f64_seed
from scipy.optimize import brentq

r1g, r2g = 0.25, 0.65
def base(rb):
    pars, _ = quartic_from_roots(r1g, rb, r2g, 0.25, kapB=1.0); V = Vfun(pars)
    return pars, (lambda x: V(x) - 4 * x * (1 - x))
rbst = brentq(lambda rb: base(rb)[1](r1g) - base(rb)[1](r2g), 0.28, 0.62, xtol=1e-12)
def make(t):
    p = dict(base(rbst)[0]); p['eps'] += t; p.update(epsA=1.0, kapA=0.0, om=0.08)
    return p

N = 40; M = N + 2
def stars(t):
    P = make(t)
    H0, dA, dB = build_dense(N, P)
    sc = np.linalg.eigvalsh(H0[:dA, :dA])[0] - np.linalg.eigvalsh(H0[dA:, dA:])[0]
    def sp(s, idx):
        H = H0.copy(); H[dA:, dA:] += s * np.eye(dB)
        E = np.linalg.eigvalsh(H); return E[idx + 1] - E[idx]
    g = np.linspace(sc - 5.4, sc + 5.4, 1501); out = []
    for idx in (0, 1):
        v = np.array([sp(s, idx) for s in g])
        lm = [i for i in range(1, len(g) - 1) if v[i] < v[i - 1] and v[i] <= v[i + 1]]
        sb = g[min(lm, key=lambda i: v[i])]; step = g[1] - g[0]
        for _ in range(12):
            g3 = np.linspace(sb - 2.2 * step, sb + 2.2 * step, 41)
            v3 = np.array([sp(s, idx) for s in g3])
            sb = g3[np.argmin(v3)]; step = g3[1] - g3[0]
        out.append((sb, sp(sb, idx)))
    return out
d = lambda t: stars(t)[0][0] - stars(t)[1][0]
scan = np.linspace(-0.105, -0.045, 7)
dv = [abs(d(t)) for t in scan]
t0 = scan[int(np.argmin(dv))]
fine = np.linspace(t0 - 0.006, t0 + 0.006, 7)
dv2 = [abs(d(t)) for t in fine]
tst = fine[int(np.argmin(dv2))]
print(f'scan |Ds*|: coarse min {min(dv):.4f} at {t0:.4f}; refined min {min(dv2):.4f} at {tst:.5f}')
st = stars(tst)
print(f'quantum-degenerate tilt = {tst:.6f}: merged at s*/M = {st[0][0]/M:.5f} / {st[1][0]/M:.5f} '
      f'(Ds* = {abs(st[0][0]-st[1][0]):.2e})')
print(f'  3-level point gaps: E1-E0 = {st[0][1]:.3e},  E2-E1 = {st[1][1]:.3e}')
P = make(tst)
H0, dA, dB = build_dense(N, P)
H = H0.copy(); H[dA:, dA:] += st[0][0] * np.eye(dB)
E, U = np.linalg.eigh(H)
nB = np.array([2 * j for j in range(dB)], float) / M
for k in (0, 1, 2):
    wA = float(np.sum(U[:dA, k] ** 2))
    xB = float(np.sum(U[dA:, k] ** 2 * nB) / max(np.sum(U[dA:, k] ** 2), 1e-300))
    print(f'  state {k}: A-weight {wA:.2f}, <x>_B = {xB:.3f}')
res = construct_invariant(P)[0]
acts = sorted(res['actions'], key=lambda a: a[0])
S1, S2 = acts[0][2], acts[1][2]
print(f'  constructor candidates: S1(0->{acts[0][0]:.3f}) = {S1:.5f}, S2(0->{acts[1][0]:.3f}) = {S2:.5f}')
print(f'  PRE-REGISTERED: Gamma = min = {S1:.5f}', flush=True)
Ns = np.arange(20, 53, 4).astype(float)
gaps = np.array([f64_seed(int(n), P)[1] for n in Ns])
ok = gaps > 1e-12
L = np.log(gaps[ok]) - np.log(np.sqrt((Ns[ok] + 1) * (Ns[ok] + 2)))
A = np.vstack([np.ones_like(Ns[ok]), Ns[ok]]).T
cg, *_ = np.linalg.lstsq(A, L, rcond=None)
print(f'  MEASURED Gamma (f64, N=20-52) = {-cg[1]:.5f}  vs S1 = {S1:.5f}, S2 = {S2:.5f}')
