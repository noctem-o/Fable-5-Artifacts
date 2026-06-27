"""
PRE-REGISTERED NUMERICAL TEST OF THE ONE-EXPONENT PRINCIPLE
===========================================================
Minimal faithful model: two-sector U(2) boson model ("IBM-CM on CP^1").

  Sector A ("normal"):  M_A = N bosons (s,t),  H_A = eps_A * n_t          -> spherical gs = |n_t=0>, E=0 exactly
  Sector B ("intruder"): M_B = N+2 bosons,     H_B = eps_B*n_t + (g/M)*n_t^2 - (kappa/M)*Q^2,  Q = s^t t + t^t s
  Coupling W = omega * (s^t s^t + t^t t^t):  maps sector A -> sector B (adds 2 bosons), preserves n_t parity.
  Control:  H(s) = [H_A  W^T; W  H_B + s*1].   Everything restricted to even n_t.

PRE-REGISTERED PREDICTIONS (locked before running):
  P1 (structure): the minimal gap, twice the direct matrix element <gsB|W|gsA>, and the
      exceptional-point height Im z (located on the FULL matrix, complex s) all coincide:
      gap = 2*w_direct = Im z   (slope of diabatic difference is exactly 1 here).
  P2 (identity): F_Q^max * (Im z)^2 = 1, i.e. 4 * chi_F(s*) * (Im z)^2 = 1.
  P3 (the law): all of {gap, w_direct, Im z} decay as  N * exp(-N*Gamma), and chi_F(s*)
      grows as exp(+2*N*Gamma)/N^2, with ONE exponent Gamma per model.
  P4 (the formula -- two rival hypotheses):
      H_FS    : Gamma = Gamma_FS  = (1/2) log(1/(1-x*))     [naive coherent/Fubini-Study overlap]
      H_AGMON : Gamma = S_WKB     = int_0^{x*} (1/2) arccosh[(V(x)-E0)/(2 kappa x(1-x)) - 1] dx
      Case C1 (eps_B=0): ground state is an EXACT coherent cat  => both hypotheses coincide (calibration).
      Case C2 (eps_B=1.2, g=0):  x*=0.35.  Squeezed gs; hypotheses may or may not split.
      Case C3 (eps_B=-0.2, g=2): SAME x*=0.35 (same Gamma_FS as C2!) but different landscape
               => S_WKB differs from C2.  C3 is the designed discriminator between H_FS and H_AGMON.
"""

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.integrate import quad

rng = np.random.default_rng(0)
EPS_A = 1.0
OMEGA = 0.2

# ---------------------------------------------------------------- sectors
def build_full(M, eps, kappa, g):
    n = np.arange(M + 1)
    Q = np.zeros((M + 1, M + 1))
    for k in range(M):
        Q[k + 1, k] = np.sqrt((k + 1) * (M - k))
        Q[k, k + 1] = Q[k + 1, k]
    H = eps * np.diag(n.astype(float)) + (g / M) * np.diag(n.astype(float) ** 2) \
        - (kappa / M) * (Q @ Q)
    return H, n

def even_restrict(H, n):
    idx = np.where(n % 2 == 0)[0]
    return H[np.ix_(idx, idx)], n[idx]

def sector(M, eps, kappa, g):
    Hf, n = build_full(M, eps, kappa, g)
    He, ne = even_restrict(Hf, n)
    E, V = np.linalg.eigh(He)
    return He, ne, E, V

# ---------------------------------------------------------------- coupling
def build_W(N, nA, nB):
    """W = OMEGA*(s+s+ + t+t+): sector A (M=N) even basis -> sector B (M=N+2) even basis."""
    dA, dB = len(nA), len(nB)
    posB = {v: i for i, v in enumerate(nB)}
    W = np.zeros((dB, dA))
    for j, m in enumerate(nA):                      # m = n_t in A
        ns = N - m
        W[posB[m], j]     += OMEGA * np.sqrt((ns + 1) * (ns + 2))   # s+s+
        if m + 2 in posB:
            W[posB[m + 2], j] += OMEGA * np.sqrt((m + 1) * (m + 2)) # t+t+
    return W

# ---------------------------------------------------------------- analytics
def classical(eps, kappa, g):
    f = lambda x: eps * x + g * x * x - 4 * kappa * x * (1 - x)
    r = minimize_scalar(f, bounds=(1e-9, 1 - 1e-9), method='bounded')
    return r.x, r.fun

def gamma_FS(xs):
    return 0.5 * np.log(1.0 / (1.0 - xs))

def S_wkb(eps, kappa, g):
    xs, E0 = classical(eps, kappa, g)
    def lam(x):
        G = (eps * x + g * x * x - E0) / (2 * kappa * x * (1 - x)) - 1.0
        return 0.5 * np.arccosh(max(G, 1.0))
    val, _ = quad(lam, 0, xs, limit=400, points=[xs * 0.5, xs * 0.9])
    return val, xs

# ---------------------------------------------------------------- per-N measurements
def measure(N, eps, kappa, g):
    MB = N + 2
    HA, nA, EA, VA = sector(N, EPS_A, 0.0, 0.0)
    HB, nB, EB, VB = sector(MB, eps, kappa, g)
    W = build_W(N, nA, nB)
    dA, dB = len(nA), len(nB)

    w_direct = abs(VB[:, 0] @ W @ VA[:, 0])
    s_center = EA[0] - EB[0]                       # diabatic crossing estimate (exact to O(w^2))

    H0 = np.zeros((dA + dB, dA + dB))
    H0[:dA, :dA] = HA
    H0[dA:, dA:] = HB
    H0[dA:, :dA] = W
    H0[:dA, dA:] = W.T
    PB = np.zeros(dA + dB); PB[dA:] = 1.0

    def full_eigh(s):
        H = H0.copy()
        H[dA:, dA:] += s * np.eye(dB)
        return np.linalg.eigh(H)

    def gap(s):
        E, _ = full_eigh(s)
        return E[1] - E[0]

    # locate s*: 2nd-order shifts displace the crossing from s_center by O(omega^2 N),
    # so do a wide coarse scan, then iterative zoom down to the anticrossing floor.
    D = 3.0 + 0.06 * N
    grid = np.linspace(s_center - D, s_center + D, 401)
    vals = np.array([gap(s) for s in grid])
    s_best = grid[np.argmin(vals)]
    step = grid[1] - grid[0]
    for _ in range(14):
        if step < max(5e-4 * vals.min(), 5e-14):
            break
        grid = np.linspace(s_best - 2.5 * step, s_best + 2.5 * step, 81)
        vals = np.array([gap(s) for s in grid])
        s_best = grid[np.argmin(vals)]
        step = grid[1] - grid[0]
    s_star, gap_min = s_best, float(vals.min())

    # chi_F at s*  via sum over states of |<n|P_B|0>|^2/(E_n-E_0)^2
    E, V = full_eigh(s_star)
    pb0 = V[:, 0] * PB
    amp = V.T @ pb0
    chi = float(np.sum((amp[1:] ** 2) / (E[1:] - E[0]) ** 2))

    # exceptional point: Newton on analytic gfun(z) = (dlam)^2, full complex-symmetric matrix
    def gfun(z):
        H = H0.astype(complex).copy()
        H[dA:, dA:] += z * np.eye(dB)
        ev = np.linalg.eigvals(H)
        ev = ev[np.argsort(ev.real)][:4]
        best, pair = np.inf, None
        for i in range(len(ev)):
            for j in range(i + 1, len(ev)):
                d = abs(ev[i] - ev[j])
                if d < best:
                    best, pair = d, (ev[i], ev[j])
        return (pair[0] - pair[1]) ** 2

    imz = np.nan
    if gap_min > 3e-6:                              # float64 EP-locating window
        z = s_star + 1.3j * gap_min
        ok = True
        for _ in range(60):
            h = max(abs(z.imag) * 1e-3, 1e-10)
            gp = (gfun(z + h) - gfun(z - h)) / (2 * h)
            if gp == 0:
                ok = False; break
            dz = gfun(z) / gp
            z -= dz
            if abs(z - s_star) > 100 * gap_min or z.imag < 0:
                ok = False; break
            if abs(dz) < 1e-5 * max(abs(z.imag), gap_min):
                break
        if ok:
            imz = abs(z.imag)

    return dict(N=N, w=w_direct, gap=gap_min, chi=chi, imz=imz, s_star=s_star)

# ---------------------------------------------------------------- fits
def fit_gamma(Ns, ys, mode):
    Ns = np.asarray(Ns, float); ys = np.asarray(ys, float)
    pref = np.sqrt((Ns + 1) * (Ns + 2))
    if mode == 'decay':       # y ~ pref * exp(-Gamma N)
        L = np.log(ys) - np.log(pref)
        sgn = -1.0
    else:                     # chi ~ exp(+2 Gamma N)/pref^2
        L = 0.5 * (np.log(ys) + 2 * np.log(pref))
        sgn = +1.0
    A = np.vstack([np.ones_like(Ns), Ns]).T
    coef, *_ = np.linalg.lstsq(A, L, rcond=None)
    res = L - A @ coef
    dof = max(len(Ns) - 2, 1)
    cov = (res @ res / dof) * np.linalg.inv(A.T @ A)
    return sgn * coef[1], np.sqrt(cov[1, 1])

# ---------------------------------------------------------------- run
CASES = [('C1 (exact-coherent anchor)', 0.0, 1.0, 0.0),
         ('C2 (squeezed, x*=0.35)',     1.2, 1.0, 0.0),
         ('C3 (engineered, same x*)',  -0.2, 1.0, 2.0)]
NLIST = list(range(8, 53, 4))
FITMIN = 20

print('=' * 78)
print('PRE-REGISTERED PREDICTIONS')
print('=' * 78)
preds = {}
for name, eps, kappa, g in CASES:
    S, xs = S_wkb(eps, kappa, g)
    gFS = gamma_FS(xs)
    preds[name] = (gFS, S, xs)
    print(f'{name:30s}  x*={xs:.4f}   Gamma_FS={gFS:.5f}   S_WKB={S:.5f}')
print()

results = {}
for name, eps, kappa, g in CASES:
    rows = [measure(N, eps, kappa, g) for N in NLIST]
    results[name] = rows
    print('-' * 78)
    print(name)
    print(f'{"N":>4} {"2*w_direct":>12} {"gap_min":>12} {"Im z":>12} {"Imz/gap":>9} {"4*chi*Imz^2":>12}')
    for r in rows:
        ident = 4 * r['chi'] * r['imz'] ** 2 if np.isfinite(r['imz']) else np.nan
        ratio = r['imz'] / r['gap'] if np.isfinite(r['imz']) else np.nan
        print(f"{r['N']:>4} {2*r['w']:>12.4e} {r['gap']:>12.4e} "
              f"{(r['imz'] if np.isfinite(r['imz']) else float('nan')):>12.4e} "
              f"{ratio:>9.4f} {ident:>12.6f}")

print()
print('=' * 78)
print('EXPONENT EXTRACTION (fit window N >= %d; quoted as Gamma +/- stderr)' % FITMIN)
print('=' * 78)
summary = {}
for name, eps, kappa, g in CASES:
    rows = [r for r in results[name] if r['N'] >= FITMIN and r['gap'] > 1e-11]
    Ns   = [r['N'] for r in rows]
    g_w,  e_w  = fit_gamma(Ns, [2 * r['w'] for r in rows], 'decay')
    g_g,  e_g  = fit_gamma(Ns, [r['gap'] for r in rows], 'decay')
    g_c,  e_c  = fit_gamma(Ns, [r['chi'] for r in rows], 'grow')
    ep_rows = [r for r in rows if np.isfinite(r['imz'])]
    if len(ep_rows) >= 3:
        g_e, e_e = fit_gamma([r['N'] for r in ep_rows], [r['imz'] for r in ep_rows], 'decay')
    else:
        g_e, e_e = np.nan, np.nan
    gFS, S, xs = preds[name]
    summary[name] = (g_w, g_g, g_e, g_c, gFS, S)
    print(f'\n{name}:   PREDICTIONS  Gamma_FS = {gFS:.5f}   S_WKB = {S:.5f}')
    print(f'  Gamma from  w_direct : {g_w:.5f} +/- {e_w:.5f}')
    print(f'  Gamma from  gap      : {g_g:.5f} +/- {e_g:.5f}')
    print(f'  Gamma from  Im z(EP) : {g_e:.5f} +/- {e_e:.5f}')
    print(f'  Gamma from  chi_F    : {g_c:.5f} +/- {e_c:.5f}')

# ---------------------------------------------------------------- figure
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-whitegrid')

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=False)
for ax, (name, eps, kappa, g) in zip(axes, CASES):
    rows = results[name]
    Ns = np.array([r['N'] for r in rows], float)
    pref = np.sqrt((Ns + 1) * (Ns + 2))
    ax.semilogy(Ns, [2 * r['w'] for r in rows], 'o', label=r'$2\,w_{\rm direct}$', ms=6)
    ax.semilogy(Ns, [r['gap'] for r in rows], 's', mfc='none', label='gap$_{\\min}$', ms=8)
    epN = [r['N'] for r in rows if np.isfinite(r['imz'])]
    epV = [r['imz'] for r in rows if np.isfinite(r['imz'])]
    ax.semilogy(epN, epV, 'x', label=r'Im$\,z$ (EP)', ms=8, color='crimson')
    gFS, S, xs = preds[name]
    for G, lab, stl in [(gFS, r'$\Gamma_{FS}$', '--'), (S, r'$S_{WKB}$', ':')]:
        ax.semilogy(Ns, 2 * OMEGA * pref * np.exp(-G * Ns) *
                    (rows[3]['gap'] / (2 * OMEGA * pref[3] * np.exp(-G * Ns[3]))),
                    stl, lw=2, label=f'slope {lab} = {G:.4f}')
    ax.set_title(name, fontsize=11, fontweight='bold')
    ax.set_xlabel('N (boson number)')
    ax.legend(fontsize=8, loc='lower left')
axes[0].set_ylabel('energy scale (log)')
fig.suptitle('One-Exponent test: three independent measurements vs two rival exponent formulas',
             fontsize=12.5, fontweight='bold')
plt.tight_layout()
plt.savefig('/mnt/user-data/outputs/one_exponent_test.png', dpi=150, bbox_inches='tight')
print('\nfigure saved.')
