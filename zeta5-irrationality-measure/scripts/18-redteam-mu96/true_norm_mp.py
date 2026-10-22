"""18-redteam-mu96: ||B_K|| via a multiprecision (mpmath, unbounded exponent range) Stieltjes/Lanczos procedure,
for K where sqrt(weights) span more than float64's range (K >~ 450). Validated against the float64 version at K=320.
Weights: float64 log-weights (relative weight errors ~1e-13 change ||B_K|| by at most that relative amount, since
(1-eta)G <= G' <= (1+eta)G implies |log||B'|| - log||B||| <= eta).
"""
import sys, time, math
import numpy as np
import mpmath as mp
sys.path.insert(0, './scripts/18-redteam-mu96')
from true_norm_lanczos import nodes_weights, log_weight, norm_from_recurrence


def stieltjes_mp(tau, logomega, h, dps=40):
    mp.mp.dps = dps
    n = len(tau)
    T = np.array([mp.mpf(float(t)) for t in tau], dtype=object)
    L = float(logomega.max())
    s = np.array([mp.exp(mp.mpf(float(l - L)) / 2) for l in logomega], dtype=object)
    nrm = mp.sqrt(mp.fsum(s * s))
    logpi0 = -0.5 * L - float(mp.log(nrm))
    v_prev = np.array([mp.mpf(0)] * n, dtype=object)
    v = s / nrm
    al, be = [], [mp.mpf(0)]
    V_keep = {0: v}
    for k in range(h):
        w = T * v
        a = mp.fsum(v * w)
        al.append(a)
        w = w - a * v - be[k] * v_prev
        if k + 1 < h:
            b = mp.sqrt(mp.fsum(w * w))
            be.append(b)
            v_prev, v = v, w / b
            if k + 1 in (h // 2, h - 1):
                V_keep[k + 1] = v
    # orthogonality diagnostics
    diag = {}
    for (i, j) in [(0, h - 1), (h // 2, h - 1), (0, h // 2)]:
        diag[(i, j)] = float(abs(mp.fsum(V_keep[i] * V_keep[j])))
    return np.array([float(a) for a in al]), np.array([float(b) for b in be]), logpi0, diag, al, be


def run(K, n=64, fine=1.0, ymax_keep=1.3, dps=40):
    t0 = time.time()
    h = 37 * K // 40
    y, wq = nodes_weights(K, n, fine)
    keep = y <= ymax_keep
    y, wq = y[keep], wq[keep]
    lo = log_weight(y, K) + np.log(wq)
    tau = y * y
    al_f, be_f, logpi0, diag, al_m, be_m = stieltjes_mp(tau, lo, h, dps)
    t1 = time.time()
    # recurrence with full-precision coefficients
    lam, it = norm_from_recurrence_mp(K, al_m, be_m, logpi0, dps=dps + 10)
    return dict(K=K, nodes=len(y), diag=diag, lognorm=float(mp.log(abs(lam))), it=it, t1=t1 - t0,
                t=time.time() - t0)


def norm_from_recurrence_mp(K, alm, bem, logpi0, dps=50, iters=60, seed=0):
    from fractions import Fraction
    from true_norm_lanczos import r_coeffs
    mp.mp.dps = dps
    N, h = 3 * K // 40, 37 * K // 40
    r = r_coeffs(K)
    js = list(range(N + 1, K + 1))
    p0 = mp.exp(mp.mpf(logpi0))
    Pi = []
    for j in js:
        x = -mp.mpf(j) ** 2 / K ** 2
        row = [p0]
        prev, cur = mp.mpf(0), p0
        for k in range(h - 1):
            nxt = ((x - alm[k]) * cur - (bem[k] * prev if k > 0 else 0)) / bem[k + 1]
            prev, cur = cur, nxt
            row.append(cur)
        Pi.append(row)
    rr = [mp.mpf(r[j].numerator) / r[j].denominator for j in js]
    K5 = mp.mpf(K) ** 5
    rng = np.random.default_rng(seed)
    v = [mp.mpf(float(z)) for z in rng.standard_normal(h)]
    lam_old = None
    for it in range(iters):
        Pv = [mp.fsum(Pi[i][a] * v[a] for a in range(h)) for i in range(len(js))]
        DPv = [rr[i] * Pv[i] for i in range(len(js))]
        Av = [K5 * mp.fsum(Pi[i][a] * DPv[i] for i in range(len(js))) for a in range(h)]
        lam = mp.fsum(Av[a] * v[a] for a in range(h)) / mp.fsum(z * z for z in v)
        nv = mp.sqrt(mp.fsum(z * z for z in Av))
        v = [z / nv for z in Av]
        if lam_old is not None and abs(lam - lam_old) < abs(lam) * mp.mpf(10) ** (-25):
            break
        lam_old = lam
    return lam, it


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        K = int(arg)
        res = run(K)
        print("K=%d nodes=%d (mp Stieltjes) : log||B_K|| = %.10f ; rate = %.8f ; orth diag %s ; [%.0fs/%.0fs]" % (
            K, res['nodes'], res['lognorm'], res['lognorm'] / K, {str(k): '%.1e' % v for k, v in res['diag'].items()},
            res['t1'], res['t']), flush=True)
