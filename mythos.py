"""MYTHOS PHASE-1 FALSIFICATION BATTERY (executed, not argued).
A1 closure   : (i) uniqueness/monotonicity of the degeneracy root; (ii) classical sbar vs the
               MEASURED quantum crossing s*(N)/M (independent observable, never compared before).
A2 basis     : Gamma_inf refit under prefactor choices {1, sqrt((N+1)(N+2)), (N+1)(N+2)} -- a
               representation artifact would shift Gamma_inf; exact theory says only b shifts.
A3 finite-N  : window-split extrapolation (low half vs high half of N-grid); pseudo-convergence
               would show window-dependent Gamma_inf drift exceeding errors.
A4 branching : compute one mid-N point under BOTH pipeline band branches; methodology imprint
               would show branch-dependent gap.
A5 leakage   : S_corr signature is (eps,kappa,g,omega) only -- no N, no measured data. KILL TEST:
               pre-registered prediction at NEVER-TESTED point C4 = (eps=-0.5, kappa=0.8, g=3.0),
               omega in {0.13, 0.05}; then blind measurement of the N->inf asymptote to N=280.
               PASS iff the two-model extrapolation bracket contains S_corr or |mean-S_corr|<=6e-4.
"""
import numpy as np
from mpmath import mp, mpf
import math, time
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import quad

EPS_A = 1.0

# ---------- machinery (verbatim from the validated extended_N2 pipeline) ----------
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
    M = N + 2; nsz = N + 3
    epsm, gm, km, om = mpf(eps), mpf(g), mpf(kappa), mpf(omega); Mm = mpf(M)
    d0 = [mpf(0)] * nsz; d1 = [mpf(0)] * nsz; d2 = [mpf(0)] * nsz; isB = [False] * nsz
    def posB(n): return n if n <= N else N + 2
    def posA(n): return n + 1
    for n in range(0, M + 1, 2):
        i = posB(n)
        d0[i] = epsm * n + gm * mpf(n * n) / Mm - km * mpf((n + 1) * (M - n) + n * (M - n + 1)) / Mm
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

def mp_gap(N, eps, kappa, g, omega, Sc, force_band=None):
    dps_need = 35 + int(0.46 * Sc * N) + 10
    mp.dps = max(45, dps_need)
    s0f, gapf, E0f = f64_seed(N, eps, kappa, g, omega)
    bands = mp_bands(N, eps, kappa, g, omega)
    pref = math.sqrt((N + 1) * (N + 2))
    use_f64 = (gapf > 1e-10) if force_band is None else (force_band == 'f64')
    if use_f64:
        gap_lo = mpf(gapf) * mpf('0.8'); gap_hi = mpf(gapf) * mpf('1.3')
    else:
        gap_lo = mpf(2 * omega) * mpf(pref) * mp.e ** (-mpf(Sc) * N)
        gap_hi = gap_lo * 25
    def E_k(s, k, tol, pad=1e-5):
        c = E0f(float(s))
        return eig_k(bands, s, k, c - 1e-5, c + pad, tol)
    s = mpf(float(s0f))
    delta = mpf(max(1e-6, 50 * float(gap_hi)))
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
    h = mpf('1.2') * gap_hi
    tolE = gap_lo * mpf('1e-6')
    padE = max(1e-5, 6.0 * float(h))
    us, g2s = [], []
    for j in (-2, -1, 0, 1, 2):
        sp = s + j * h
        e0 = E_k(sp, 0, tolE, padE); e1 = E_k(sp, 1, tolE, padE)
        us.append(mpf(j)); g2s.append(((e1 - e0) / h) ** 2)
    S0 = mpf(5); S1 = sum(us); S2 = sum(u * u for u in us); S3 = sum(u ** 3 for u in us); S4 = sum(u ** 4 for u in us)
    T0 = sum(g2s); T1 = sum(u * y for u, y in zip(us, g2s)); T2 = sum(u * u * y for u, y in zip(us, g2s))
    sol = mp.lu_solve(mp.matrix([[S4, S3, S2], [S3, S2, S1], [S2, S1, S0]]), mp.matrix([T2, T1, T0]))
    A, B, C = sol[0], sol[1], sol[2]
    g2min = (C - B * B / (4 * A)) * h * h
    return (mp.sqrt(g2min) if g2min > 0 else None), float(s)

def S_corr_full(eps, kappa, g, omega):
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
    F = lambda s: E_left(s) - E_right(s)[0]
    sb = brentq(F, sbar0 - 0.5, sbar0 + 1.0, xtol=1e-13)
    Ec = E_left(sb); xmin = E_right(sb)[1]
    def lam(x):
        a = EPS_A * x - Ec
        K = 2 * kappa * x * (1 - x)
        B1 = V(x) + sb - Ec - K
        P = (1 - x) ** 2 + x ** 2
        R = 2 * x * (1 - x)
        c = (a * B1 - omega * omega * P) / (a * K + omega * omega * R)
        return 0.5 * np.arccosh(c) if c > 1.0 else 0.0
    S = quad(lam, 1e-8, xmin, limit=600, points=[xmin * 0.5, xmin * 0.9])[0]
    return S, sb, F

def fit_two(Nmid, loc):
    A1 = np.vstack([np.ones_like(Nmid), 1.0 / Nmid]).T
    c1, *_ = np.linalg.lstsq(A1, loc, rcond=None)
    r1 = loc - A1 @ c1; s1 = np.sqrt((r1 @ r1 / max(len(Nmid) - 2, 1)) * np.linalg.inv(A1.T @ A1)[0, 0])
    A2 = np.vstack([np.ones_like(Nmid), 1.0 / Nmid, np.log(Nmid) / Nmid]).T
    c2, *_ = np.linalg.lstsq(A2, loc, rcond=None)
    r2 = loc - A2 @ c2; s2 = np.sqrt((r2 @ r2 / max(len(Nmid) - 3, 1)) * np.linalg.inv(A2.T @ A2)[0, 0])
    return c1, s1, c2, s2

# ================= A5: PRE-REGISTERED OUT-OF-SAMPLE KILL TEST =================
C4 = dict(eps=-0.5, kappa=0.8, g=3.0)          # never-tested landscape; kappa != 1 for the first time
OMS4 = [0.13, 0.05]                             # never-tested couplings
print('=' * 78)
print('A5 KILL TEST -- pre-registered predictions at NEVER-TESTED point C4 (eps=-0.5, kappa=0.8, g=3.0)')
preds = {}
for om in OMS4:
    S, sb, _ = S_corr_full(C4['eps'], C4['kappa'], C4['g'], om)
    preds[om] = (S, sb)
    print(f'  PREDICTION: S_corr(C4, omega={om}) = {S:.5f}   (classical sbar = {sb:.6f})')
print('  PASS criterion: two-model bracket contains S_corr, or |mean - S_corr| <= 6e-4.')
NGRID = [64, 88, 112, 136, 160, 184, 208, 232, 256, 280]
A5 = {}
for om in OMS4:
    S = preds[om][0]; t0 = time.time(); lgs = []
    for N in NGRID:
        gm, _ = mp_gap(N, C4['eps'], C4['kappa'], C4['g'], om, S)
        lgs.append(float(mp.log(gm)) - float(np.log(np.sqrt((N + 1) * (N + 2)))))
    Ns = np.array(NGRID, float); L = np.array(lgs)
    loc = (L[:-1] - L[1:]) / np.diff(Ns); Nmid = 0.5 * (Ns[:-1] + Ns[1:])
    c1, s1, c2, s2 = fit_two(Nmid, loc)
    A5[om] = (Nmid, loc, c1, s1, c2, s2, L)
    brk = (min(c1[0], c2[0]) <= S <= max(c1[0], c2[0]))
    print(f'\n  MEASURED C4 omega={om} [{time.time()-t0:.0f}s]: slopes {np.array2string(loc, precision=5)}')
    print(f'    G+b/N      : {c1[0]:.5f} +/- {s1:.5f}  (diff {c1[0]-S:+.5f}, b = {c1[1]:.3f})')
    print(f'    G+b/N+ln/N : {c2[0]:.5f} +/- {s2:.5f}  (diff {c2[0]-S:+.5f})')
    print(f'    bracket contains S_corr: {brk}   |mean-S| = {abs(0.5*(c1[0]+c2[0])-S):.5f}')

# ================= A1: CLOSURE =================
print('\n' + '=' * 78)
print('A1 CLOSURE -- (i) degeneracy-root uniqueness; (ii) classical sbar vs measured quantum s*(N)/M')
for tag, eps, kappa, g, om in [('C1', 0.0, 1.0, 0.0, 0.2), ('C4', C4['eps'], C4['kappa'], C4['g'], 0.13)]:
    S, sb, F = S_corr_full(eps, kappa, g, om)
    grid = np.linspace(sb - 0.4, sb + 0.6, 2001)
    vals = np.array([F(s) for s in grid])
    signchg = int(np.sum(np.diff(np.sign(vals)) != 0))
    mono = bool(np.all(np.diff(vals) < 0))
    Ns = np.array([24, 32, 40, 48, 56, 64])
    sstars = np.array([f64_seed(int(N), eps, kappa, g, om)[0] for N in Ns])
    Ms = Ns + 2.0
    Afit = np.vstack([Ms, np.ones_like(Ms), 1.0 / Ms]).T
    cf, *_ = np.linalg.lstsq(Afit, sstars, rcond=None)
    print(f'  {tag} om={om}: root sign-changes = {signchg} (want 1), F monotone-decreasing = {mono}')
    print(f'    sbar_classical = {sb:.6f}   slope of measured s*(N) vs M = {cf[0]:.6f}   diff = {cf[0]-sb:+.2e}')

# ================= A2: BASIS / PREFACTOR =================
print('\n' + '=' * 78)
print('A2 BASIS -- Gamma_inf under prefactor choices (C4, omega=0.13 series reused)')
Nmid, loc, c1, s1, c2, s2, L = A5[0.13]
Ns = np.array(NGRID, float)
for pname, lp in [('pref=1', 0.0 * Ns), ('pref=sqrt((N+1)(N+2))', np.log(np.sqrt((Ns + 1) * (Ns + 2)))),
                  ('pref=(N+1)(N+2)', np.log((Ns + 1) * (Ns + 2)))]:
    Lp = (L + np.log(np.sqrt((Ns + 1) * (Ns + 2)))) - lp     # rebuild log(gap) then strip chosen pref
    locp = (Lp[:-1] - Lp[1:]) / np.diff(Ns)
    cp, sp, cp2, sp2 = fit_two(Nmid, locp)
    print(f'  {pname:>24s}: G+b/N -> {cp[0]:.5f} (b={cp[1]:+.3f});  +ln/N -> {cp2[0]:.5f}')

# ================= A3: WINDOW SPLIT =================
print('\n' + '=' * 78)
print('A3 FINITE-N -- window-split extrapolation (C4, omega=0.13)')
mask_lo = Nmid <= 160; mask_hi = Nmid > 160
for tag, m in [('low  N<=160', mask_lo), ('high N>160 ', mask_hi)]:
    cw, sw, cw2, sw2 = fit_two(Nmid[m], loc[m])
    print(f'  {tag}: G+b/N -> {cw[0]:.5f} +/- {sw:.5f};  +ln/N -> {cw2[0]:.5f} +/- {sw2:.5f}')

# ================= A4: BRANCH BOUNDARY =================
print('\n' + '=' * 78)
print('A4 BRANCHING -- same point computed under both pipeline band branches (C4, omega=0.13)')
S = preds[0.13][0]
for N in [88]:
    g_f, _ = mp_gap(N, C4['eps'], C4['kappa'], C4['g'], 0.13, S, force_band='f64')
    g_p, _ = mp_gap(N, C4['eps'], C4['kappa'], C4['g'], 0.13, S, force_band='pred')
    rel = abs(float(g_f) - float(g_p)) / float(g_f)
    print(f'  N={N}: f64-band gap = {mp.nstr(g_f, 9)}   pred-band gap = {mp.nstr(g_p, 9)}   rel.diff = {rel:.2e}')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(7, 4.6))
for om, col in zip(OMS4, ['C0', 'C1']):
    Nmid, loc, c1, s1, c2, s2, L = A5[om]
    S = preds[om][0]
    ax.plot(Nmid, loc, 'o', color=col, ms=5, label=f'C4 $\\omega$={om}: measured slopes')
    xs = np.linspace(Nmid[0], 3000, 400)
    ax.plot(xs, c1[0] + c1[1] / xs, '--', color=col, lw=1.2)
    ax.axhline(S, color=col, ls=':', lw=1.6, label=f'pre-registered $S_{{corr}}$ = {S:.4f}')
ax.set_xscale('log'); ax.set_xlabel('N'); ax.set_ylabel(r'$\Gamma_{\rm loc}(N)$')
ax.set_title('A5 out-of-sample kill test: never-tested parameters (C4), blind asymptote',
             fontsize=10, fontweight='bold')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/mythos_kill_test.png', dpi=150, bbox_inches='tight')
print('\nfigure saved.')
