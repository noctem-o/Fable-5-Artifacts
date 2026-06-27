"""EXTENDED-PRECISION N-SCALING v2 (hardened pipeline; same pre-registered criteria as v1).
Fixes from v1 tether-catches:
  - adaptive stage offsets: delta scaled to predicted gap (v1 failed when gap > 1e-8 at small N);
  - adaptive precision: dps = 35 + 0.46*S_corr*N + 10 (v1 fixed tolerances died at Gamma*N ~ 87);
  - hyperbola finisher: 5-point parabola fit to gap^2(s), exact for the local 2x2, robust when the
    crossing center is known only to ~gap accuracy; quartic model error ~ 27*kappa2*gap (negligible).
Criteria unchanged: H-a closes Conjecture 5 if |Gamma_inf - S_corr| <= 3*(stderr + model spread)."""
import numpy as np
from mpmath import mp, mpf
import time, math

EPS_A = 1.0

def build_dense(N, eps, kappa, g, omega):
    M = N + 2
    nf = np.arange(M + 1)
    Q = np.zeros((M + 1, M + 1))
    for k in range(M):
        Q[k + 1, k] = np.sqrt((k + 1) * (M - k)); Q[k, k + 1] = Q[k + 1, k]
    HBf = eps * np.diag(nf.astype(float)) + (g / M) * np.diag(nf.astype(float) ** 2) - (kappa / M) * (Q @ Q)
    iB = np.where(nf % 2 == 0)[0]
    HB = HBf[np.ix_(iB, iB)]; nB = nf[iB]
    nA = np.arange(0, N + 1, 2)
    HA = np.diag(EPS_A * nA.astype(float))
    dA, dB = len(nA), len(nB)
    posB = {v: i for i, v in enumerate(nB)}
    W = np.zeros((dB, dA))
    for j, m in enumerate(nA):
        ns = N - m
        W[posB[m], j] += omega * np.sqrt((ns + 1) * (ns + 2))
        if m + 2 in posB:
            W[posB[m + 2], j] += omega * np.sqrt((m + 1) * (m + 2))
    H0 = np.zeros((dA + dB, dA + dB))
    H0[:dA, :dA] = HA; H0[dA:, dA:] = HB
    H0[dA:, :dA] = W;  H0[:dA, dA:] = W.T
    return H0, dA, dB

def f64_seed(N, eps, kappa, g, omega):
    H0, dA, dB = build_dense(N, eps, kappa, g, omega)
    EB0 = np.linalg.eigvalsh(H0[dA:, dA:])[0]
    s_center = -EB0
    def gap(s):
        H = H0.copy(); H[dA:, dA:] += s * np.eye(dB)
        E = np.linalg.eigvalsh(H)
        return E[1] - E[0]
    D = 3.0 + 0.06 * N
    grid = np.linspace(s_center - D, s_center + D, 401)
    vals = np.array([gap(s) for s in grid])
    s_best = grid[np.argmin(vals)]; step = grid[1] - grid[0]
    for _ in range(16):
        if step < 5e-13: break
        grid = np.linspace(s_best - 2.5 * step, s_best + 2.5 * step, 81)
        vals = np.array([gap(s) for s in grid])
        s_best = grid[np.argmin(vals)]; step = grid[1] - grid[0]
    def E0f(s):
        H = H0.copy(); H[dA:, dA:] += s * np.eye(dB)
        return np.linalg.eigvalsh(H)[0]
    return s_best, float(vals.min()), E0f

def mp_bands(N, eps, kappa, g, omega):
    M = N + 2
    nsz = N + 3
    epsm, gm, km, om = mpf(eps), mpf(g), mpf(kappa), mpf(omega)
    Mm = mpf(M)
    d0 = [mpf(0)] * nsz; d1 = [mpf(0)] * nsz; d2 = [mpf(0)] * nsz
    isB = [False] * nsz
    def posB(n): return n if n <= N else N + 2
    def posA(n): return n + 1
    for n in range(0, M + 1, 2):
        i = posB(n)
        q2 = mpf((n + 1) * (M - n) + n * (M - n + 1))
        d0[i] = epsm * n + gm * mpf(n * n) / Mm - km * q2 / Mm
        isB[i] = True
    for n in range(0, N + 1, 2):
        d0[posA(n)] = mpf(EPS_A) * n
    for n in range(0, M - 1, 2):
        d2[posB(n)] = -km * mp.sqrt(mpf((n + 1) * (n + 2) * (M - n) * (M - n - 1))) / Mm
    for n in range(0, N + 1, 2):
        ns = N - n
        d1[posB(n)] = om * mp.sqrt(mpf((ns + 1) * (ns + 2)))
        d1[posA(n)] = om * mp.sqrt(mpf((n + 1) * (n + 2)))
    return d0, d1, d2, isB

def count_below(bands, s, sigma):
    d0, d1, d2, isB = bands
    n = len(d0); neg = 0
    Dp1 = mpf(0); Dp2 = mpf(0); l1p = mpf(0); l2p1 = mpf(0); l2p2 = mpf(0)
    guard = mpf(10) ** (-(mp.dps + 5))
    for i in range(n):
        di = d0[i] - sigma
        if isB[i]: di += s
        di -= Dp1 * l1p * l1p + Dp2 * l2p2 * l2p2
        if di == 0: di = guard
        l1c = (d1[i] - Dp1 * l1p * l2p1) / di if i + 1 < n else mpf(0)
        l2c = d2[i] / di if i + 2 < n else mpf(0)
        if di < 0: neg += 1
        Dp2 = Dp1; Dp1 = di
        l2p2 = l2p1; l2p1 = l2c; l1p = l1c
    return neg

def eig_k(bands, s, k, lo, hi, tol):
    lo, hi = mpf(lo), mpf(hi)
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if count_below(bands, s, mid) > k: hi = mid
        else: lo = mid
    return (lo + hi) / 2

def mp_gap(N, eps, kappa, g, omega, Sc):
    dps_need = 35 + int(0.46 * Sc * N) + 10
    mp.dps = max(45, dps_need)
    s0f, gapf, E0f = f64_seed(N, eps, kappa, g, omega)
    bands = mp_bands(N, eps, kappa, g, omega)
    pref = math.sqrt((N + 1) * (N + 2))
    if gapf > 1e-10:                                   # float64 gap is reliable: tight band
        gap_lo = mpf(gapf) * mpf('0.8'); gap_hi = mpf(gapf) * mpf('1.3')
    else:                                              # deep regime: S_corr-based band
        gap_lo = mpf(2 * omega) * mpf(pref) * mp.e ** (-mpf(Sc) * N)
        gap_hi = gap_lo * 25
    def E_k(s, k, tol, pad=1e-5):
        c = E0f(float(s))
        return eig_k(bands, s, k, c - 1e-5, c + pad, tol)
    # staged localization of s* (diabatic-branch line intersection), adaptive offsets
    s = mpf(float(s0f))
    delta = mpf(max(1e-6, 50 * float(gap_hi)))
    eps_s = delta
    for _ in range(5):
        tau = max(mpf('0.06') * delta ** 2, gap_lo * mpf('1e-7'), mpf(10) ** (-(mp.dps - 6)))
        EL = [E_k(s - 2 * delta, 0, tau), E_k(s - delta, 0, tau)]
        ER = [E_k(s + delta, 0, tau), E_k(s + 2 * delta, 0, tau)]
        mB = (EL[1] - EL[0]) / delta; bB = EL[1] - mB * (s - delta)
        mA = (ER[1] - ER[0]) / delta; bA = ER[0] - mA * (s + delta)
        s = (bA - bB) / (mB - mA)
        eps_s = mpf('0.12') * delta ** 2 + 3 * tau
        if eps_s <= mpf('0.7') * gap_hi: break
        delta = max(50 * gap_hi, 12 * eps_s)
    # hyperbola finisher: parabola fit to gap^2 over 5 points spanning the crossing
    h = mpf('1.2') * gap_hi
    tolE = gap_lo * mpf('1e-6')
    padE = max(1e-5, 6.0 * float(h))
    us, g2s = [], []
    for j in (-2, -1, 0, 1, 2):
        sp = s + j * h
        e0 = E_k(sp, 0, tolE, padE); e1 = E_k(sp, 1, tolE, padE)
        us.append(mpf(j)); g2s.append(((e1 - e0) / h) ** 2)   # work in units of h for conditioning
    # least-squares quadratic g2 = A u^2 + B u + C  (normal equations, exact small sums)
    S0 = mpf(5); S1 = sum(us); S2 = sum(u * u for u in us); S3 = sum(u ** 3 for u in us); S4 = sum(u ** 4 for u in us)
    T0 = sum(g2s); T1 = sum(u * y for u, y in zip(us, g2s)); T2 = sum(u * u * y for u, y in zip(us, g2s))
    Mm = mp.matrix([[S4, S3, S2], [S3, S2, S1], [S2, S1, S0]])
    rhs = mp.matrix([T2, T1, T0])
    sol = mp.lu_solve(Mm, rhs)
    A, B, C = sol[0], sol[1], sol[2]
    g2min = (C - B * B / (4 * A)) * h * h
    if g2min <= 0:
        return None
    return mp.sqrt(g2min)

from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import quad
def S_corr(eps, kappa, g, omega):
    V = lambda x: eps * x + g * x * x
    hB0 = lambda x: V(x) - 4 * kappa * x * (1 - x)
    r = minimize_scalar(hB0, bounds=(1e-9, 1 - 1e-9), method='bounded')
    sbar0 = -r.fun
    hA = lambda x: EPS_A * x
    E_left = lambda sb: 0.5 * sb - np.sqrt(0.25 * sb * sb + omega * omega)
    def E_right(sb):
        f = lambda x: 0.5 * (hA(x) + hB0(x) + sb) - np.sqrt(0.25 * (hA(x) - hB0(x) - sb) ** 2 + omega ** 2)
        rr = minimize_scalar(f, bounds=(1e-9, 1 - 1e-9), method='bounded')
        return rr.fun, rr.x
    sb = brentq(lambda s: E_left(s) - E_right(s)[0], sbar0 - 0.5, sbar0 + 1.0, xtol=1e-13)
    Ec = E_left(sb); xmin = E_right(sb)[1]
    def lam(x):
        a = EPS_A * x - Ec
        K = 2 * kappa * x * (1 - x)
        B1 = V(x) + sb - Ec - K
        P = (1 - x) ** 2 + x ** 2
        R = 2 * x * (1 - x)
        c = (a * B1 - omega * omega * P) / (a * K + omega * omega * R)
        return 0.5 * np.arccosh(c) if c > 1.0 else 0.0
    return quad(lam, 1e-8, xmin, limit=600, points=[xmin * 0.5, xmin * 0.9])[0]

CASES = [('C1', 0.0, 1.0, 0.0), ('C3', -0.2, 1.0, 2.0)]
NGRID = [64, 88, 112, 136, 160, 184, 208, 232, 256, 280]
OMS = [0.2, 0.1]

print('=== T1 sanity: hardened mp pipeline vs float64 at N=64 (both omegas, both cases) ===', flush=True)
for om in OMS:
    for name, eps, kappa, g in CASES:
        Sc = S_corr(eps, kappa, g, om)
        _, gf, _ = f64_seed(64, eps, kappa, g, om)
        gm = mp_gap(64, eps, kappa, g, om, Sc)
        rel = abs(float(gm) - gf) / gf
        print(f'  {name} om={om}: f64 = {gf:.6e}   mp = {mp.nstr(gm, 8)}   rel.diff = {rel:.2e}', flush=True)

results = {}
for om in OMS:
    for name, eps, kappa, g in CASES:
        Sc = S_corr(eps, kappa, g, om)
        t0 = time.time()
        lgs = []
        for N in NGRID:
            gm = mp_gap(N, eps, kappa, g, om, Sc)
            lg = float(mp.log(gm)) - float(np.log(np.sqrt((N + 1) * (N + 2))))
            lgs.append(lg)
        Ns = np.array(NGRID, float); L = np.array(lgs)
        loc = (L[:-1] - L[1:]) / np.diff(Ns)
        Nmid = 0.5 * (Ns[:-1] + Ns[1:])
        A1 = np.vstack([np.ones_like(Nmid), 1.0 / Nmid]).T
        c1, *_ = np.linalg.lstsq(A1, loc, rcond=None)
        r1 = loc - A1 @ c1; s1e = np.sqrt((r1 @ r1 / max(len(Nmid) - 2, 1)) * np.linalg.inv(A1.T @ A1)[0, 0])
        A2 = np.vstack([np.ones_like(Nmid), 1.0 / Nmid, np.log(Nmid) / Nmid]).T
        c2, *_ = np.linalg.lstsq(A2, loc, rcond=None)
        r2 = loc - A2 @ c2; s2e = np.sqrt((r2 @ r2 / max(len(Nmid) - 3, 1)) * np.linalg.inv(A2.T @ A2)[0, 0])
        results[(name, om)] = (Nmid, loc, Sc, c1, s1e, c2, s2e)
        print(f'\n{name} om={om}  [{time.time()-t0:.0f}s]: slopes {np.array2string(loc, precision=5)}')
        print(f'  S_corr = {Sc:.5f}')
        print(f'  fit G+b/N         : G_inf = {c1[0]:.5f} +/- {s1e:.5f}  (G-S_corr = {c1[0]-Sc:+.5f}, b = {c1[1]:.3f})')
        print(f'  fit G+b/N+c lnN/N : G_inf = {c2[0]:.5f} +/- {s2e:.5f}  (G-S_corr = {c2[0]-Sc:+.5f})', flush=True)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-whitegrid')
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
for ax, om in zip(axes, OMS):
    for name, col in zip(['C1', 'C3'], ['C0', 'C1']):
        Nmid, loc, Sc, c1, s1e, c2, s2e = results[(name, om)]
        ax.plot(Nmid, loc, 'o', color=col, ms=5, label=f'{name}: slopes (mp, to N=280)')
        xs = np.linspace(Nmid[0], 3000, 500)
        ax.plot(xs, c1[0] + c1[1] / xs, '--', color=col, lw=1.2,
                label=f'{name}: $\\Gamma_\\infty$={c1[0]:.4f} (1/N fit)')
        ax.axhline(Sc, color=col, ls=':', lw=1.6)
    ax.set_xscale('log')
    ax.set_xlabel('N'); ax.set_title(f'$\\omega$ = {om}   (dotted: derived $S_{{corr}}$)', fontsize=10.5, fontweight='bold')
    ax.legend(fontsize=7.5)
axes[0].set_ylabel(r'local slope $\Gamma_{\rm loc}(N)$')
fig.suptitle('Decision experiment (hardened pipeline): spectral exponents vs derived $S_{corr}(\\omega)$',
             fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/extended_N_scaling.png', dpi=150, bbox_inches='tight')
print('\nfigure saved.', flush=True)
