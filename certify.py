"""THEOREM CERTIFICATION SUITE
(A) Theorem 2 (anchor): closed-form w_direct for C1 vs numerics -> machine precision expected.
(B) Claim 4 derivation: classical adiabatic action S_adiab(omega) on the lower sheet of the
    2x2 symbol  h(x,chi) = [[eps_A*x, omega],[omega, V(x)+sbar-2k x(1-x)(1+cos2chi)]].
    Instanton condition derived in closed form:  (h_A - E_c)(h_B + sbar - E_c) = omega^2
    along the path, with sbar tuned so the two adiabatic minima are degenerate at E_c.
    Consistency: S_adiab(omega->0) -> S_WKB.
(C) Confrontation: measured Gamma_gap(omega) (N-fits, 20<=N<=52) vs S_adiab(omega), both cases;
    PLUS extended-N (up to 64) local slopes at omega in {0.2, 0.1} -- the adiabatic-regime test:
    prediction is that local slopes drift toward S_adiab(omega) as N grows (LZ parameter N*omega^2 grows).
"""
import numpy as np
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import quad

EPS_A = 1.0

# ---------------- shared machinery (as certified in earlier runs) ----------------
def build_full(M, eps, kappa, g):
    n = np.arange(M + 1)
    Q = np.zeros((M + 1, M + 1))
    for k in range(M):
        Q[k + 1, k] = np.sqrt((k + 1) * (M - k)); Q[k, k + 1] = Q[k + 1, k]
    H = eps * np.diag(n.astype(float)) + (g / M) * np.diag(n.astype(float) ** 2) - (kappa / M) * (Q @ Q)
    return H, n

def sector(M, eps, kappa, g):
    Hf, n = build_full(M, eps, kappa, g)
    idx = np.where(n % 2 == 0)[0]
    He, ne = Hf[np.ix_(idx, idx)], n[idx]
    E, V = np.linalg.eigh(He)
    return He, ne, E, V

def build_W(N, nA, nB, omega):
    posB = {v: i for i, v in enumerate(nB)}
    W = np.zeros((len(nB), len(nA)))
    for j, m in enumerate(nA):
        ns = N - m
        W[posB[m], j] += omega * np.sqrt((ns + 1) * (ns + 2))
        if m + 2 in posB:
            W[posB[m + 2], j] += omega * np.sqrt((m + 1) * (m + 2))
    return W

def measure_gap(N, eps, kappa, g, omega):
    MB = N + 2
    HA, nA, EA, VA = sector(N, EPS_A, 0.0, 0.0)
    HB, nB, EB, VB = sector(MB, eps, kappa, g)
    W = build_W(N, nA, nB, omega)
    dA, dB = len(nA), len(nB)
    w_direct = abs(VB[:, 0] @ W @ VA[:, 0])
    s_center = EA[0] - EB[0]
    H0 = np.zeros((dA + dB, dA + dB))
    H0[:dA, :dA] = HA; H0[dA:, dA:] = HB
    H0[dA:, :dA] = W;  H0[:dA, dA:] = W.T
    def gap(s):
        H = H0.copy(); H[dA:, dA:] += s * np.eye(dB)
        E = np.linalg.eigvalsh(H)
        return E[1] - E[0]
    D = 3.0 + 0.06 * N
    grid = np.linspace(s_center - D, s_center + D, 401)
    vals = np.array([gap(s) for s in grid])
    s_best = grid[np.argmin(vals)]; step = grid[1] - grid[0]
    for _ in range(14):
        if step < max(5e-4 * vals.min(), 5e-14): break
        grid = np.linspace(s_best - 2.5 * step, s_best + 2.5 * step, 81)
        vals = np.array([gap(s) for s in grid])
        s_best = grid[np.argmin(vals)]; step = grid[1] - grid[0]
    return w_direct, float(vals.min())

def fit_gamma(Ns, ys):
    Ns = np.asarray(Ns, float); L = np.log(ys) - np.log(np.sqrt((Ns + 1) * (Ns + 2)))
    A = np.vstack([np.ones_like(Ns), Ns]).T
    c, *_ = np.linalg.lstsq(A, L, rcond=None)
    r = L - A @ c; cov = (r @ r / max(len(Ns) - 2, 1)) * np.linalg.inv(A.T @ A)
    return -c[1], np.sqrt(cov[1, 1])

def S_wkb(eps, kappa, g):
    f = lambda x: eps * x + g * x * x - 4 * kappa * x * (1 - x)
    r = minimize_scalar(f, bounds=(1e-9, 1 - 1e-9), method='bounded')
    xs, E0 = r.x, r.fun
    lam = lambda x: 0.5 * np.arccosh(max((eps * x + g * x * x - E0) / (2 * kappa * x * (1 - x)) - 1.0, 1.0))
    return quad(lam, 0, xs, limit=400, points=[xs * 0.5, xs * 0.9])[0]

# ---------------- (A) Theorem 2: closed-form anchor ----------------
print('=== (A) Theorem 2 certification: C1 closed form  w = 2*sqrt(2)*omega*sqrt((N+1)(N+2))*2^{-(N+2)/2}')
om0 = 0.2
for N in [8, 16, 32, 48]:
    w_meas, _ = (lambda r: r)(measure_gap(N, 0.0, 1.0, 0.0, om0))[0], None
    w_meas = measure_gap(N, 0.0, 1.0, 0.0, om0)[0]
    w_pred = 2 * np.sqrt(2) * om0 * np.sqrt((N + 1) * (N + 2)) * 2 ** (-(N + 2) / 2)
    print(f'  N={N:>3}:  w_meas={w_meas:.12e}   w_pred={w_pred:.12e}   rel.err={abs(w_meas-w_pred)/w_pred:.2e}')

# ---------------- (B) classical adiabatic action ----------------
def S_adiab(eps, kappa, g, omega):
    V = lambda x: eps * x + g * x * x
    hB0 = lambda x: V(x) - 4 * kappa * x * (1 - x)          # chi = 0 diabatic curve
    r = minimize_scalar(hB0, bounds=(1e-9, 1 - 1e-9), method='bounded')
    x0, E0d = r.x, r.fun
    sbar0 = -E0d
    hA = lambda x: EPS_A * x
    def E_left(sb):                                          # exact 2x2 at x = 0
        return 0.5 * sb - np.sqrt(0.25 * sb * sb + omega * omega)
    def E_right(sb):
        f = lambda x: 0.5 * (hA(x) + hB0(x) + sb) - np.sqrt(0.25 * (hA(x) - hB0(x) - sb) ** 2 + omega ** 2)
        rr = minimize_scalar(f, bounds=(1e-9, 1 - 1e-9), method='bounded')
        return rr.fun, rr.x
    F = lambda sb: E_left(sb) - E_right(sb)[0]
    sb = brentq(F, sbar0 - 0.5, sbar0 + 1.0, xtol=1e-13)
    Ec = E_left(sb)
    xmin = E_right(sb)[1]                       # exact right-well location at degeneracy
    def lam(x):
        a = hA(x) - Ec
        arg = (V(x) + sb - Ec - omega * omega / a) / (2 * kappa * x * (1 - x)) - 1.0
        return 0.5 * np.arccosh(arg) if arg > 1.0 else 0.0
    # path: x=0 (boundary well) -> barrier -> right well at xmin(omega); integrate exactly to xmin.
    xR = xmin
    val, _ = quad(lam, 1e-10, xR, limit=600, points=[xR * 0.5, xR * 0.9])
    return val, sb, Ec

print('\n=== (B) S_adiab solver:  omega->0 consistency check ===')
for name, eps, kappa, g in [('C1', 0.0, 1.0, 0.0), ('C3', -0.2, 1.0, 2.0)]:
    S0 = S_wkb(eps, kappa, g)
    Sa, _, _ = S_adiab(eps, kappa, g, 1e-4)
    print(f'  {name}: S_WKB = {S0:.6f}   S_adiab(1e-4) = {Sa:.6f}   diff = {abs(S0-Sa):.2e}')

OMEGAS = np.array([0.2, 0.15, 0.1, 0.07, 0.05, 0.035, 0.025, 0.0175, 0.0125, 0.008])
NLIST = list(range(20, 53, 4))
print('\n=== (C) confrontation: measured Gamma_gap(omega)  vs  derived S_adiab(omega) ===')
conf = {}
for name, eps, kappa, g in [('C1', 0.0, 1.0, 0.0), ('C3', -0.2, 1.0, 2.0)]:
    S0 = S_wkb(eps, kappa, g)
    print(f'\n{name}  (S_WKB = {S0:.5f})')
    print(f'{"omega":>8} {"Gamma_gap (meas)":>18} {"S_adiab (derived)":>18} {"meas - derived":>15} {"N*om^2 range":>14}')
    rows = []
    for om in OMEGAS:
        gaps = [measure_gap(N, eps, kappa, g, om)[1] for N in NLIST]
        Gg, eg = fit_gamma(NLIST, gaps)
        Sa, _, _ = S_adiab(eps, kappa, g, om)
        rows.append((om, Gg, eg, Sa))
        print(f'{om:>8.4f} {Gg:>12.5f}+/-{eg:.5f} {Sa:>18.5f} {Gg-Sa:>15.5f}   {NLIST[0]*om**2:>5.2f}-{NLIST[-1]*om**2:.2f}')
    conf[name] = (S0, rows)

print('\n=== (C2) adiabatic-regime test: extended-N local slopes at omega in {0.2, 0.1} ===')
ext = {}
for name, eps, kappa, g in [('C1', 0.0, 1.0, 0.0), ('C3', -0.2, 1.0, 2.0)]:
    for om in (0.2, 0.1):
        Ns = np.array(list(range(20, 65, 4)), float)
        gaps = np.array([measure_gap(int(N), eps, kappa, g, om)[1] for N in Ns])
        L = np.log(gaps) - np.log(np.sqrt((Ns + 1) * (Ns + 2)))
        loc = (L[:-1] - L[1:]) / np.diff(Ns)
        Sa, _, _ = S_adiab(eps, kappa, g, om)
        ext[(name, om)] = (Ns, loc, Sa)
        print(f'{name} omega={om}: local slopes N=22..62: {np.array2string(loc, precision=5)}')
        print(f'      S_adiab = {Sa:.5f}   last local slope = {loc[-1]:.5f}   (N*om^2 at end = {Ns[-1]*om**2:.1f})')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-whitegrid')
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
ax = axes[0]
for name, col in zip(['C1', 'C3'], ['C0', 'C1']):
    S0, rows = conf[name]
    om = [r[0] for r in rows]; Gg = [r[1] for r in rows]; Sa = [r[3] for r in rows]
    ax.semilogx(om, Gg, 'o', color=col, ms=6, label=f'{name}: measured $\\Gamma_{{gap}}$ (N=20-52 fit)')
    ax.semilogx(om, Sa, '-', color=col, lw=2, label=f'{name}: derived $S_{{adiab}}(\\omega)$')
    ax.axhline(S0, color=col, ls=':', lw=1.2)
ax.set_xlabel(r'$\omega$'); ax.set_ylabel('exponent')
ax.set_title('Derived classical adiabatic action vs measured exponent\n(dotted: $S_{WKB}$ asymptotes)',
             fontsize=10, fontweight='bold')
ax.legend(fontsize=7.5)
ax = axes[1]
for (name, om), (Ns, loc, Sa) in ext.items():
    col = 'C0' if name == 'C1' else 'C1'
    ls = '-' if om == 0.2 else '--'
    ax.plot(0.5 * (Ns[:-1] + Ns[1:]), loc, ls, marker='o', ms=4, color=col, label=f'{name}, $\\omega$={om}')
    ax.axhline(Sa, color=col, ls=':' if om == 0.1 else '-.', lw=1.0)
ax.set_xlabel('N'); ax.set_ylabel(r'local slope $\Gamma_{\rm loc}(N)$')
ax.set_title('Adiabatic-regime test: local slopes vs derived $S_{adiab}$\n(horizontal lines)',
             fontsize=10, fontweight='bold')
ax.legend(fontsize=7.5)
plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/theorem_certification.png', dpi=150, bbox_inches='tight')
print('\nfigure saved.')
