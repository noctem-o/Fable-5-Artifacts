"""MYTHOS ROUND 2 -- STRUCTURAL ATTACKS (pre-registered).
ALPHA: dynamical A-sector (kappa_A > 0, both sheets bend). Generalized characteristic-variety
  condition (h_A - E_c)(h_B + sbar - E_c) = |w|^2 becomes QUADRATIC in c = cosh(2 lambda):
     KA*KB*c^2 - (A1*KB + B1*KA + om^2 R)*c + (A1*B1 - om^2 P) = 0,
  A1 = epsA x - E_c - KA, B1 = V + sbar - E_c - KB, KA = 2 kapA x(1-x), KB = 2 kapB x(1-x).
  BRANCH RULE (derived, not fitted): physical root has a+b>0  <=>  c < (A1+B1)/(KA+KB)  <=> smaller root.
  P-ALPHA: Gamma_inf = S_gen, |Delta| <= 6e-4. Consistency: S_gen(kapA->0) = S_corr (old).
BETA: quartic landscape with sub-degenerate intermediate dip (f' roots 0.22, 0.40, 0.62; barrier 0.4,
  dip 0.17 above E0). In-domain (H3 holds). P-BETA: Gamma_inf = S_corr(quartic), dip is a spectator.
GAMMA (out-of-domain probe): degenerate triple well f(x_m) = E0 = f(x*). (H3) violated.
  P-GAMMA (low confidence): global gap-min exponent = S1 + S2 (virtual intermediate, O(1) detuning);
  document phenomenology either way."""
import numpy as np
from mpmath import mp, mpf
import math, time
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import quad

# ---------------- generalized builders: V quartic, kappa_A >= 0 ----------------
def Vfun(p):
    return lambda x: p['eps'] * x + p['g2'] * x * x + p['g3'] * x ** 3 + p['g4'] * x ** 4

def build_sector_full(Mtot, eps_lin, kap, g2, g3, g4):
    nf = np.arange(Mtot + 1)
    Q = np.zeros((Mtot + 1, Mtot + 1))
    for k in range(Mtot):
        Q[k + 1, k] = np.sqrt((k + 1) * (Mtot - k)); Q[k, k + 1] = Q[k + 1, k]
    nn = nf.astype(float)
    H = eps_lin * np.diag(nn) + (g2 / Mtot) * np.diag(nn ** 2) + (g3 / Mtot ** 2) * np.diag(nn ** 3) \
        + (g4 / Mtot ** 3) * np.diag(nn ** 4) - (kap / Mtot) * (Q @ Q)
    iE = np.where(nf % 2 == 0)[0]
    return H[np.ix_(iE, iE)]

def build_dense(N, p):
    M = N + 2
    HB = build_sector_full(M, p['eps'], p['kapB'], p['g2'], p['g3'], p['g4'])
    HA = build_sector_full(N, p['epsA'], p['kapA'], 0.0, 0.0, 0.0)
    dA, dB = HA.shape[0], HB.shape[0]
    om = p['om']
    W = np.zeros((dB, dA))
    for j in range(dA):                       # A even value n = 2j -> B even index
        m = 2 * j; ns = N - m
        W[j, j] += om * np.sqrt((ns + 1) * (ns + 2))          # s+s+: B(n) at index n/2 = j
        if j + 1 < dB:
            W[j + 1, j] += om * np.sqrt((m + 1) * (m + 2))    # t+t+: B(n+2) at index j+1
    H0 = np.zeros((dA + dB, dA + dB))
    H0[:dA, :dA] = HA; H0[dA:, dA:] = HB
    H0[dA:, :dA] = W;  H0[:dA, dA:] = W.T
    return H0, dA, dB

def f64_seed(N, p):
    H0, dA, dB = build_dense(N, p)
    EA0 = np.linalg.eigvalsh(H0[:dA, :dA])[0]
    EB0 = np.linalg.eigvalsh(H0[dA:, dA:])[0]
    s_center = EA0 - EB0
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

def mp_bands(N, p):
    M = N + 2; nsz = N + 3
    om = mpf(p['om']); kB = mpf(p['kapB']); kA = mpf(p['kapA'])
    e1, g2, g3, g4 = mpf(p['eps']), mpf(p['g2']), mpf(p['g3']), mpf(p['g4'])
    eA = mpf(p['epsA']); Mm = mpf(M); Nm = mpf(N)
    d0 = [mpf(0)] * nsz; d1 = [mpf(0)] * nsz; d2 = [mpf(0)] * nsz; isB = [False] * nsz
    def pB(n): return n if n <= N else N + 2
    def pA(n): return n + 1
    for n in range(0, M + 1, 2):
        i = pB(n); nm = mpf(n)
        d0[i] = e1 * nm + g2 * nm ** 2 / Mm + g3 * nm ** 3 / Mm ** 2 + g4 * nm ** 4 / Mm ** 3 \
                - kB * mpf((n + 1) * (M - n) + n * (M - n + 1)) / Mm
        isB[i] = True
    for n in range(0, N + 1, 2):
        d0[pA(n)] = eA * n - kA * mpf((n + 1) * (N - n) + n * (N - n + 1)) / Nm
    for n in range(0, M - 1, 2):
        d2[pB(n)] = -kB * mp.sqrt(mpf((n + 1) * (n + 2) * (M - n) * (M - n - 1))) / Mm
    for n in range(0, N - 1, 2):
        d2[pA(n)] = -kA * mp.sqrt(mpf((n + 1) * (n + 2) * (N - n) * (N - n - 1))) / Nm
    for n in range(0, N + 1, 2):
        ns = N - n
        d1[pB(n)] = om * mp.sqrt(mpf((ns + 1) * (ns + 2)))
        d1[pA(n)] = om * mp.sqrt(mpf((n + 1) * (n + 2)))
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

def mp_gap(N, p, Sc):
    mp.dps = max(45, 35 + int(0.46 * Sc * N) + 10)
    s0f, gapf, E0f = f64_seed(N, p)
    bands = mp_bands(N, p)
    pref = math.sqrt((N + 1) * (N + 2))
    if gapf > 1e-10:
        gap_lo = mpf(gapf) * mpf('0.8'); gap_hi = mpf(gapf) * mpf('1.3')
    else:
        gap_lo = mpf(2 * p['om']) * mpf(pref) * mp.e ** (-mpf(Sc) * N)
        gap_hi = gap_lo * 25
    def E_k(s, k, tol, pad=1e-5):
        c = E0f(float(s))
        return eig_k(bands, s, k, c - 1e-5, c + pad, tol)
    s = mpf(float(s0f))
    # FIX (ALPHA N=64 fault): for large gaps the localization loop diverges (eps_s ~ delta^2
    # can never satisfy <= 0.7*gap_hi, delta grows x12+/iter, s exits the crossing; the
    # finisher then samples the linear branch -> gap overestimated ~x10). The f64 seed is
    # already refined to 5e-13 there, so skip mp refinement when gapf is float64-trustworthy.
    delta = mpf(max(1e-6, 50 * float(gap_hi)))
    for _ in range(0 if gapf > 1e-6 else 5):
        tau = max(mpf('0.06') * delta ** 2, gap_lo * mpf('1e-7'), mpf(10) ** (-(mp.dps - 6)))
        EL = [E_k(s - 2 * delta, 0, tau), E_k(s - delta, 0, tau)]
        ER = [E_k(s + delta, 0, tau), E_k(s + 2 * delta, 0, tau)]
        mB = (EL[1] - EL[0]) / delta; bB = EL[1] - mB * (s - delta)
        mA = (ER[1] - ER[0]) / delta; bA = ER[0] - mA * (s + delta)
        s = (bA - bB) / (mB - mA)
        if mpf('0.12') * delta ** 2 + 3 * tau <= mpf('0.7') * gap_hi: break
        delta = max(50 * gap_hi, 12 * (mpf('0.12') * delta ** 2 + 3 * tau))
    h = mpf('1.2') * gap_hi
    tolE = gap_lo * mpf('1e-6'); padE = max(1e-5, 6.0 * float(h))
    us, g2s = [], []
    for j in (-2, -1, 0, 1, 2):
        sp = s + j * h
        e0 = E_k(sp, 0, tolE, padE); e1 = E_k(sp, 1, tolE, padE)
        us.append(mpf(j)); g2s.append(((e1 - e0) / h) ** 2)
    S1 = sum(us); S2 = sum(u * u for u in us); S3 = sum(u ** 3 for u in us); S4 = sum(u ** 4 for u in us)
    T0 = sum(g2s); T1 = sum(u * y for u, y in zip(us, g2s)); T2 = sum(u * u * y for u, y in zip(us, g2s))
    S0c = mpf(5)
    sol = mp.lu_solve(mp.matrix([[S4, S3, S2], [S3, S2, S1], [S2, S1, S0c]]), mp.matrix([T2, T1, T0]))
    Aq, Bq, Cq = sol[0], sol[1], sol[2]
    g2min = (Cq - Bq * Bq / (4 * Aq)) * h * h
    return mp.sqrt(g2min) if g2min > 0 else None

# ---------------- generalized quadrature: quadratic in c, branch a+b>0 (= smaller root) ----------------
def S_gen(p):
    om = p['om']; V = Vfun(p)
    epsA, kapA, kapB = p['epsA'], p['kapA'], p['kapB']
    hA0 = lambda x: epsA * x - 4 * kapA * x * (1 - x)
    hB0 = lambda x: V(x) - 4 * kapB * x * (1 - x)
    xs = np.linspace(1e-6, 1 - 1e-6, 4001)
    fb = np.array([hB0(x) for x in xs])
    x0 = xs[np.argmin(fb)]
    r = minimize_scalar(hB0, bounds=(max(1e-9, x0 - 0.02), min(1 - 1e-9, x0 + 0.02)), method='bounded')
    sbar0 = -r.fun
    E_left = lambda sb: 0.5 * sb - np.sqrt(0.25 * sb * sb + om * om)
    def E_right(sb):
        f = lambda x: 0.5 * (hA0(x) + hB0(x) + sb) - np.sqrt(0.25 * (hA0(x) - hB0(x) - sb) ** 2 + om ** 2)
        vals = np.array([f(x) for x in xs])
        xm = xs[np.argmin(vals)]
        rr = minimize_scalar(f, bounds=(max(1e-9, xm - 0.02), min(1 - 1e-9, xm + 0.02)), method='bounded')
        return rr.fun, rr.x
    sb = brentq(lambda s: E_left(s) - E_right(s)[0], sbar0 - 0.6, sbar0 + 1.2, xtol=1e-13)
    Ec = E_left(sb); xmin = E_right(sb)[1]
    def lam(x):
        KA = 2 * kapA * x * (1 - x); KB = 2 * kapB * x * (1 - x)
        A1 = epsA * x - Ec - KA
        B1 = V(x) + sb - Ec - KB
        P = (1 - x) ** 2 + x ** 2; R = 2 * x * (1 - x)
        if kapA == 0.0:
            den = A1 * KB + om * om * R
            c = (A1 * B1 - om * om * P) / den if den > 0 else 1.0
        else:
            aq = KA * KB
            bq = -(A1 * KB + B1 * KA + om * om * R)
            cq = A1 * B1 - om * om * P
            disc = bq * bq - 4 * aq * cq
            if disc < 0: return 0.0
            c = (-bq - np.sqrt(disc)) / (2 * aq)     # smaller root: F(c*)<0 at a+b=0 proves it lies below c*
        return 0.5 * np.arccosh(c) if c > 1.0 else 0.0
    S, _ = quad(lam, 1e-8, xmin, limit=800, points=[xmin * 0.3, xmin * 0.6, xmin * 0.9])
    return S, sb, Ec, xmin, lam

# ---------------- quartic landscape designers ----------------
def quartic_from_roots(r1, rb, r2, barrier, kapB=1.0):
    c3 = 1.0; c2 = -(r1 + rb + r2); c1 = r1 * rb + r1 * r2 + rb * r2; c0 = -r1 * rb * r2
    F = lambda x: c3 * x ** 4 / 4 + c2 * x ** 3 / 3 + c1 * x ** 2 / 2 + c0 * x
    A = barrier / (F(rb) - F(r2))
    pars = dict(eps=A * c0 + 4 * kapB, g2=A * c1 / 2 - 4 * kapB, g3=A * c2 / 3, g4=A * c3 / 4, kapB=kapB)
    fvals = dict(fm=A * F(r1), fb=A * F(rb), fs=A * F(r2))
    return pars, fvals

def fit_two(Nmid, loc):
    A1 = np.vstack([np.ones_like(Nmid), 1.0 / Nmid]).T
    c1, *_ = np.linalg.lstsq(A1, loc, rcond=None)
    r1 = loc - A1 @ c1; s1 = np.sqrt((r1 @ r1 / max(len(Nmid) - 2, 1)) * np.linalg.inv(A1.T @ A1)[0, 0])
    A2 = np.vstack([np.ones_like(Nmid), 1.0 / Nmid, np.log(Nmid) / Nmid]).T
    c2, *_ = np.linalg.lstsq(A2, loc, rcond=None)
    r2 = loc - A2 @ c2; s2 = np.sqrt((r2 @ r2 / max(len(Nmid) - 3, 1)) * np.linalg.inv(A2.T @ A2)[0, 0])
    return c1, s1, c2, s2

NGRID = [64, 88, 112, 136, 160, 184, 208, 232, 256, 280]
def run_series(p, Sc, tag):
    t0 = time.time(); lgs = []
    for N in NGRID:
        gm = mp_gap(N, p, Sc)
        lgs.append(float(mp.log(gm)) - float(np.log(np.sqrt((N + 1) * (N + 2)))))
    Ns = np.array(NGRID, float); L = np.array(lgs)
    loc = (L[:-1] - L[1:]) / np.diff(Ns); Nmid = 0.5 * (Ns[:-1] + Ns[1:])
    c1, s1, c2, s2 = fit_two(Nmid, loc)
    brk = min(c1[0], c2[0]) <= Sc <= max(c1[0], c2[0])
    print(f'\n  MEASURED {tag} [{time.time()-t0:.0f}s]: slopes {np.array2string(loc, precision=5)}')
    print(f'    G+b/N      : {c1[0]:.5f} +/- {s1:.5f}  (diff {c1[0]-Sc:+.5f}, b = {c1[1]:.3f})')
    print(f'    G+b/N+ln/N : {c2[0]:.5f} +/- {s2:.5f}  (diff {c2[0]-Sc:+.5f})')
    print(f'    bracket contains prediction: {brk}   |mean-S| = {abs(0.5*(c1[0]+c2[0])-Sc):.5f}', flush=True)
    return Nmid, loc, c1, c2

def construct_invariant(p, verbose=False):
    om = p['om']; V = Vfun(p)
    epsA, kapA, kapB = p['epsA'], p['kapA'], p['kapB']
    hA0 = lambda x: epsA * x - 4 * kapA * x * (1 - x)
    hB0 = lambda x: V(x) - 4 * kapB * x * (1 - x)
    sheet = lambda x, sb: 0.5 * (hA0(x) + hB0(x) + sb) - np.sqrt(0.25 * (hA0(x) - hB0(x) - sb) ** 2 + om ** 2)
    xs = np.linspace(1e-6, 1 - 1e-6, 6001)
    fb = np.array([hB0(x) for x in xs]); sbar0 = -fb.min()
    E_left = lambda sb: 0.5 * sb - np.sqrt(0.25 * sb * sb + om * om)
    def sheet_global(sb):
        v = sheet(xs, sb); i = np.argmin(v)
        r = minimize_scalar(lambda x: sheet(x, sb),
                            bounds=(max(1e-9, xs[i] - 0.01), min(1 - 1e-9, xs[i] + 0.01)), method='bounded')
        return r.fun, r.x
    F = lambda s: E_left(s) - sheet_global(s)[0]
    sgrid = np.linspace(sbar0 - 0.8, sbar0 + 1.6, 481)
    Fv = np.array([F(s) for s in sgrid])
    roots = [brentq(F, sgrid[i], sgrid[i + 1], xtol=1e-13)
             for i in range(len(sgrid) - 1) if Fv[i] * Fv[i + 1] < 0]
    if not roots:
        raise RuntimeError('no degeneracy root')
    out = []
    for sb in roots:
        Ec = E_left(sb)
        sv = sheet(xs, sb)
        mins = []
        for i in range(1, len(xs) - 1):
            if sv[i] < sv[i - 1] and sv[i] <= sv[i + 1]:
                r = minimize_scalar(lambda x: sheet(x, sb),
                                    bounds=(xs[i - 1], xs[i + 1]), method='bounded')
                mins.append((r.x, r.fun))
        mins.sort(key=lambda t: t[1])
        Delta = lambda x: (hA0(x) - Ec) * (hB0(x) + sb - Ec) - om * om   # w(x,0) = om for this symbol
        Dv = np.array([Delta(x) for x in xs])
        proots = [brentq(Delta, xs[i], xs[i + 1], xtol=1e-13)
                  for i in range(len(xs) - 1) if Dv[i] * Dv[i + 1] < 0]
        def lam(x):
            KA = 2 * kapA * x * (1 - x); KB = 2 * kapB * x * (1 - x)
            A1 = epsA * x - Ec - KA; B1 = V(x) + sb - Ec - KB
            P = (1 - x) ** 2 + x ** 2; R = 2 * x * (1 - x)
            if kapA == 0.0:
                den = A1 * KB + om * om * R
                c = (A1 * B1 - om * om * P) / den if den > 0 else 1.0
            else:
                aq = KA * KB; bq = -(A1 * KB + B1 * KA + om * om * R); cq = A1 * B1 - om * om * P
                disc = bq * bq - 4 * aq * cq
                if disc < 0: return 0.0
                c = (-bq - np.sqrt(disc)) / (2 * aq)
            return 0.5 * np.arccosh(c) if c > 1.0 else 0.0
        actions = []
        for (xw, vw) in mins:
            pts = sorted([r for r in proots if 1e-8 < r < xw]) + [xw * 0.5]
            S = quad(lam, 1e-8, xw, limit=900, points=sorted(set(pts)))[0]
            actions.append((xw, vw, S))
        x_sel, v_sel = mins[0]
        F1 = len(mins) > 1 and (mins[1][1] - mins[0][1]) < 5e-3
        out.append(dict(sb=sb, Ec=Ec, mins=mins, patch_roots=proots, actions=actions,
                        x_sel=x_sel, S_sel=[a[2] for a in actions if a[0] == x_sel][0], F1_flag=F1))
        if verbose:
            print(f'  degeneracy root sbar = {sb:.6f}, Ec = {Ec:.6f}')
            for (xw, vw, S) in actions:
                print(f'    sheet min at x = {xw:.4f}: dressed value = {vw:.6f} (vs Ec {Ec:.6f}), S(0->x) = {S:.5f}')
            print(f'    Delta sign changes at: {[f"{r:.4f}" for r in proots]}   F1 flag: {F1}')
    return out

