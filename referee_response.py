"""Referee-response analysis.
Demands addressed:
 (R1) Is the correction exponent exactly 2? -> free-p fit of  S_WKB - Gamma_gap(omega) = c * omega^p
      over ~1.4 decades, with stderr; fitted on full range AND on the small-omega half (asymptotic p).
 (R2) Pure omega^2 vs omega^2*log(1/omega) -> residual comparison of both one-parameter models.
 (R3) Finite-size scaling for the exponent claims -> local slopes Gamma_loc(N) per consecutive-N pair,
      with 1/N extrapolation, instead of a single global fit.
"""
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

CASES = [('C1', 0.0, 1.0, 0.0), ('C3', -0.2, 1.0, 2.0)]
OMEGAS = np.array([0.2, 0.15, 0.1, 0.07, 0.05, 0.035, 0.025, 0.0175, 0.0125, 0.008])
NLIST = list(range(20, 53, 4))

store = {}
print('=== fine omega scan ===')
print(f'{"case":>5} {"omega":>8} {"Gamma_gap":>22} {"deficit":>12} {"deficit/om^2":>13}')
for name, eps, kappa, g in CASES:
    S = S_wkb(eps, kappa, g)
    rows = []
    for om in OMEGAS:
        data = [measure_gap(N, eps, kappa, g, om) for N in NLIST]
        gaps = [d[1] for d in data]
        Gg, eg = fit_gamma(NLIST, gaps)
        deficit = S - Gg
        rows.append((om, Gg, eg, deficit, gaps))
        print(f'{name:>5} {om:>8.4f} {Gg:>14.5f} +/- {eg:.5f} {deficit:>12.6f} {deficit/om**2:>13.4f}')
    store[name] = (S, rows)

print()
print('=== (R1) free-exponent fit:  deficit = c * omega^p  (weighted by fit stderr) ===')
for name in [c[0] for c in CASES]:
    S, rows = store[name]
    om = np.array([r[0] for r in rows]); dfc = np.array([r[3] for r in rows]); sig = np.array([r[2] for r in rows])
    ok = dfc > 3 * sig
    def freefit(mask):
        x, y, w = np.log(om[mask]), np.log(dfc[mask]), (dfc[mask] / sig[mask]) ** 2
        A = np.vstack([np.ones_like(x), x]).T
        Wm = np.diag(w)
        cov = np.linalg.inv(A.T @ Wm @ A)
        c = cov @ A.T @ Wm @ y
        r = y - A @ c
        s2 = (r @ Wm @ r) / max(mask.sum() - 2, 1)
        return c[1], np.sqrt(s2 * cov[1, 1]), np.exp(c[0])
    p_all, ep_all, c_all = freefit(ok)
    small = ok & (om <= 0.05)
    p_sm, ep_sm, c_sm = freefit(small)
    print(f'{name}:  p(full 1.4 decades) = {p_all:.3f} +/- {ep_all:.3f}   '
          f'p(omega<=0.05, asymptotic) = {p_sm:.3f} +/- {ep_sm:.3f}   c = {c_sm:.3f}')
    # (R2) one-parameter model comparison on the asymptotic window
    def chi2(model):
        m = small
        c = np.sum(dfc[m] * model(om[m]) / sig[m] ** 2) / np.sum(model(om[m]) ** 2 / sig[m] ** 2)
        return np.sum(((dfc[m] - c * model(om[m])) / sig[m]) ** 2), c
    chi_q, c_q = chi2(lambda o: o ** 2)
    chi_l, c_l = chi2(lambda o: o ** 2 * np.log(1.0 / o))
    print(f'      model comparison (omega<=0.05): chi2[c*om^2] = {chi_q:.2f} (c={c_q:.3f})   '
          f'chi2[c*om^2*log(1/om)] = {chi_l:.2f} (c={c_l:.3f})')

print()
print('=== (R3) finite-size convergence: local slopes Gamma_loc(N), gap channel ===')
for name in [c[0] for c in CASES]:
    S, rows = store[name]
    for om_pick in (0.05, 0.0125):
        row = [r for r in rows if abs(r[0] - om_pick) < 1e-12][0]
        gaps = np.array(row[4]); Ns = np.array(NLIST, float)
        L = np.log(gaps) - np.log(np.sqrt((Ns + 1) * (Ns + 2)))
        loc = (L[:-1] - L[1:]) / np.diff(Ns)
        Nmid = 0.5 * (Ns[:-1] + Ns[1:])
        A = np.vstack([np.ones_like(Nmid), 1.0 / Nmid]).T
        c, *_ = np.linalg.lstsq(A, loc, rcond=None)
        print(f'{name} omega={om_pick}: local slopes {np.array2string(loc, precision=5)}')
        print(f'      1/N extrapolation -> Gamma_inf = {c[0]:.5f}   (S_WKB = {S:.5f}, '
              f'S_WKB - deficit_fit = {S - row[3]:.5f})')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-whitegrid')
fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
ax = axes[0]
for name, mk in zip(['C1', 'C3'], ['o', 's']):
    S, rows = store[name]
    om = np.array([r[0] for r in rows]); dfc = np.array([r[3] for r in rows]); sig = np.array([r[2] for r in rows])
    ax.errorbar(om, dfc, yerr=sig, fmt=mk, ms=6, label=f'{name}: deficit $S_{{WKB}}-\\Gamma_{{gap}}(\\omega)$', capsize=3)
oref = np.linspace(0.008, 0.2, 100)
ax.plot(oref, 2.1 * oref ** 2, 'k--', lw=1.5, label=r'$2.1\,\omega^2$ reference')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r'$\omega$'); ax.set_ylabel('exponent deficit')
ax.set_title('Correction law on log-log axes (slope = p)', fontsize=10.5, fontweight='bold')
ax.legend(fontsize=8)
ax = axes[1]
for name, col in zip(['C1', 'C3'], ['C0', 'C1']):
    S, rows = store[name]
    for om_pick, ls in [(0.05, '--'), (0.0125, '-')]:
        row = [r for r in rows if abs(r[0] - om_pick) < 1e-12][0]
        gaps = np.array(row[4]); Ns = np.array(NLIST, float)
        L = np.log(gaps) - np.log(np.sqrt((Ns + 1) * (Ns + 2)))
        loc = (L[:-1] - L[1:]) / np.diff(Ns)
        Nmid = 0.5 * (Ns[:-1] + Ns[1:])
        ax.plot(Nmid, loc, ls, marker='o', ms=4, color=col, label=f'{name}, $\\omega$={om_pick}')
    ax.axhline(S, color=col, ls=':', lw=1.5)
ax.set_xlabel('N (midpoint of slope pair)'); ax.set_ylabel(r'local slope $\Gamma_{\rm loc}(N)$')
ax.set_title('Finite-size convergence of the gap exponent\n(dotted lines: $S_{WKB}$ per case)',
             fontsize=10.5, fontweight='bold')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/referee_response_analysis.png', dpi=150, bbox_inches='tight')
print('\nfigure saved.')
