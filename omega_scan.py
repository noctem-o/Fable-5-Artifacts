"""Follow-up: regime test. Hypothesis from the surprise in run 1:
at extensive mixing (omega fixed) the SPECTRAL exponent (gap/EP/chi) is the
adiabatic-surface action S_adiab(omega) < S_WKB, while the bare matrix element
w_direct always carries S_WKB. Prediction: Gamma_gap(omega) -> S_WKB as omega -> 0,
while Gamma_w stays pinned at S_WKB for all omega."""
import numpy as np
from scipy.optimize import minimize_scalar
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

def S_wkb(eps, kappa, g):
    f = lambda x: eps * x + g * x * x - 4 * kappa * x * (1 - x)
    r = minimize_scalar(f, bounds=(1e-9, 1 - 1e-9), method='bounded')
    xs, E0 = r.x, r.fun
    lam = lambda x: 0.5 * np.arccosh(max((eps * x + g * x * x - E0) / (2 * kappa * x * (1 - x)) - 1.0, 1.0))
    return quad(lam, 0, xs, limit=400, points=[xs * 0.5, xs * 0.9])[0]

def measure_lite(N, eps, kappa, g, omega):
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

CASES = [('C1', 0.0, 1.0, 0.0), ('C3', -0.2, 1.0, 2.0)]
OMEGAS = [0.2, 0.1, 0.05, 0.025, 0.0125]
NLIST = list(range(20, 53, 4))

results = {}
print(f'{"case":>5} {"omega":>8} {"Gamma_gap":>20} {"Gamma_w":>20} {"S_WKB":>9}')
for name, eps, kappa, g in CASES:
    S = S_wkb(eps, kappa, g)
    for om in OMEGAS:
        ws, gs = [], []
        for N in NLIST:
            w, gp = measure_lite(N, eps, kappa, g, om)
            ws.append(2 * w); gs.append(gp)
        Gg, eg = fit_gamma(NLIST, gs)
        Gw, ew = fit_gamma(NLIST, ws)
        results[(name, om)] = (Gg, eg, Gw, ew, S)
        print(f'{name:>5} {om:>8.4f} {Gg:>12.5f}+/-{eg:.5f} {Gw:>12.5f}+/-{ew:.5f} {S:>9.5f}')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-whitegrid')
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
for ax, (name, eps, kappa, g) in zip(axes, CASES):
    S = S_wkb(eps, kappa, g)
    oms = OMEGAS
    Gg = [results[(name, o)][0] for o in oms]
    Gw = [results[(name, o)][2] for o in oms]
    ax.semilogx(oms, Gg, 'o-', label=r'$\Gamma_{\rm gap}(\omega)$  (spectral exponent)', ms=7)
    ax.semilogx(oms, Gw, 's--', label=r'$\Gamma_{w}(\omega)$  (matrix-element exponent)', ms=7)
    ax.axhline(S, color='k', lw=1.5, ls=':', label=f'$S_{{WKB}}$ = {S:.4f}')
    ax.set_xlabel(r'mixing strength $\omega$')
    ax.set_title(f'{name}: does $\\Gamma_{{\\rm gap}}\\to S_{{WKB}}$ as $\\omega\\to 0$?',
                 fontsize=10.5, fontweight='bold')
    ax.legend(fontsize=9)
axes[0].set_ylabel(r'measured exponent $\Gamma$')
fig.suptitle('Regime test: do the two exponent classes merge in the weak-mixing (R-II) limit?',
             fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/omega_regime_test.png', dpi=150, bbox_inches='tight')
print('figure saved.')
