"""CONJECTURE 5 TEST (pre-registered).
Run discovered: measured spectral exponent falls BELOW S_adiab(const-omega) for N*omega^2 >~ 0.5.
Candidate resolution found analytically BEFORE this run: the coupling operator
W = omega(s+s+ + t+t+) has a chi-DEPENDENT symbol, w(x,chi) = omega[(1-x) + x e^{2i chi}],
because t+t+ moves n_t. Under the barrier (chi = i lambda):
   |w|^2 -> omega^2 [ (1-x)^2 + x^2 + 2x(1-x) cosh 2lambda ]   (exponentially amplified),
and the instanton condition  a*b = |w|^2  stays LINEAR in c = cosh 2lambda:
   c(x) = ( a*B1 - omega^2 P ) / ( a*K + omega^2 R ),
   a = eps_A x - E_c,  B1 = V + sbar - E_c - K,  K = 2 kappa x(1-x),
   P = (1-x)^2 + x^2,  R = 2x(1-x).
The boundary data (E_left, E_right, sbar, E_c, x_min) are UNCHANGED since |w| = omega at chi=0 and x=0.

PRE-REGISTERED HYPOTHESIS H-coupling:
  (i)  S_corr(omega->0) = S_Agmon  (consistency);
  (ii) S_corr(omega) < S_adiab(omega) for omega > 0;
  (iii) the corrected residuals  Gamma_meas - S_corr  shrink to the finite-N drift scale
        (|.| <~ few x 1e-3, positive-biased) ACROSS ALL omega, removing the systematic
        -0.017 / -0.011 deviations.  If a systematic residual survives at large omega,
        a genuine second (Stokes) channel remains and Conjecture 5 stands for the remainder.
"""
import numpy as np
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import quad

EPS_A = 1.0

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
    return float(vals.min())

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

def boundary_data(eps, kappa, g, omega):
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
    return sb, E_left(sb), E_right(sb)[1], V

def S_const(eps, kappa, g, omega):                # previous formula (constant coupling symbol)
    sb, Ec, xmin, V = boundary_data(eps, kappa, g, omega)
    def lam(x):
        a = EPS_A * x - Ec
        arg = (V(x) + sb - Ec - omega * omega / a) / (2 * kappa * x * (1 - x)) - 1.0
        return 0.5 * np.arccosh(arg) if arg > 1.0 else 0.0
    return quad(lam, 1e-8, xmin, limit=600, points=[xmin * 0.5, xmin * 0.9])[0]

def S_corr(eps, kappa, g, omega):                 # chi-dependent coupling symbol (closed form, linear in c)
    sb, Ec, xmin, V = boundary_data(eps, kappa, g, omega)
    def lam(x):
        a  = EPS_A * x - Ec
        K  = 2 * kappa * x * (1 - x)
        B1 = V(x) + sb - Ec - K
        P  = (1 - x) ** 2 + x ** 2
        R  = 2 * x * (1 - x)
        c  = (a * B1 - omega * omega * P) / (a * K + omega * omega * R)
        return 0.5 * np.arccosh(c) if c > 1.0 else 0.0
    return quad(lam, 1e-8, xmin, limit=600, points=[xmin * 0.5, xmin * 0.9])[0]

CASES = [('C1', 0.0, 1.0, 0.0), ('C3', -0.2, 1.0, 2.0)]
OMEGAS = np.array([0.2, 0.15, 0.1, 0.07, 0.05, 0.035, 0.025, 0.0175, 0.0125, 0.008])
NLIST = list(range(20, 53, 4))

print('=== consistency (i): omega -> 0 ===')
for name, eps, kappa, g in CASES:
    print(f'  {name}: S_WKB = {S_wkb(eps,kappa,g):.6f}   S_corr(1e-4) = {S_corr(eps,kappa,g,1e-4):.6f}')

print('\n=== confrontation with chi-dependent coupling symbol (zero free parameters) ===')
store = {}
for name, eps, kappa, g in CASES:
    print(f'\n{name}')
    print(f'{"omega":>8} {"Gamma_meas":>12} {"S_const":>10} {"S_corr":>10} {"meas-S_const":>13} {"meas-S_corr":>12}')
    rows = []
    for om in OMEGAS:
        gaps = [measure_gap(N, eps, kappa, g, om) for N in NLIST]
        Gg, _ = fit_gamma(NLIST, gaps)
        Sc = S_const(eps, kappa, g, om)
        Sx = S_corr(eps, kappa, g, om)
        rows.append((om, Gg, Sc, Sx))
        print(f'{om:>8.4f} {Gg:>12.5f} {Sc:>10.5f} {Sx:>10.5f} {Gg-Sc:>13.5f} {Gg-Sx:>12.5f}')
    store[name] = rows

print('\n=== extended-N endpoint check (last local slopes from previous run) vs S_corr ===')
prev = {('C1', 0.2): 0.25881, ('C1', 0.1): 0.31540, ('C3', 0.2): 0.15703, ('C3', 0.1): 0.20897}
for (name, om), slope in prev.items():
    eps, kappa, g = [(c[1], c[2], c[3]) for c in CASES if c[0] == name][0]
    print(f'  {name} omega={om}: last local slope (N=62) = {slope:.5f}   S_corr = {S_corr(eps,kappa,g,om):.5f}'
          f'   S_const = {S_const(eps,kappa,g,om):.5f}')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(7.5, 5))
for name, col in zip(['C1', 'C3'], ['C0', 'C1']):
    rows = store[name]
    om = [r[0] for r in rows]
    ax.semilogx(om, [r[1] for r in rows], 'o', color=col, ms=6, label=f'{name}: measured')
    ax.semilogx(om, [r[2] for r in rows], '--', color=col, lw=1.5, label=f'{name}: $S_{{const}}$ (old)')
    ax.semilogx(om, [r[3] for r in rows], '-', color=col, lw=2.2, label=f'{name}: $S_{{corr}}$ ($\\chi$-dep. coupling)')
ax.set_xlabel(r'$\omega$'); ax.set_ylabel('spectral exponent')
ax.set_title('Conjecture-5 test: corrected coupling symbol vs measured exponents\n(zero free parameters)',
             fontsize=10.5, fontweight='bold')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/conjecture5_test.png', dpi=150, bbox_inches='tight')
print('\nfigure saved.')
