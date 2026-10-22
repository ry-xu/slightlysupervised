"""18-redteam-mu96: ||B_K|| for larger K via orthonormal polynomials (own implementation).

||B_K|| = K^5 * || Pi^T diag(r) Pi ||_2,  Pi_{ja} = pi_a(x_j), x_j = -(j/K)^2 (N<j<=K), pi_a orthonormal in L^2(W_K dtau)
(16's Lemma 1(c); identical to the definition B_K = G^-1/2 M_K G^-1/2 by Lemma 1(a),(b), which we verified exactly).
  int p(tau)^2 W_K(tau) dtau = int_0^oo p(y^2)^2 Omega(y^2) w(K y) dy.
Discretisation: composite Gauss-Legendre in y (graded panels near 0); float64 Lanczos with full reorthogonalisation
on the vector sqrt(omega_i) (all quantities O(1)); recurrence coefficients -> pi_a(x_j) by forward recurrence in mpmath;
spectral norm by power iteration in mpmath.
r_j = u_j^4 prod_{l<=N}(u_l^2-u_j^2)^5 / prod_{N<k<=K,k!=j}(u_k^2-u_j^2)   (Fauzan (C.3)).
"""
import sys, time, math
import numpy as np
import mpmath as mp
from fractions import Fraction


def panels(K, fine=1.0):
    edges = [0.0]
    y = 1.0 / (4 * K)
    edges.append(y)
    while y < 0.05:
        y *= 1.0 + 0.2 / fine
        edges.append(min(y, 0.05))
    y = edges[-1]
    while y < 1.2:
        y = min(y + 0.01 / fine, 1.2)
        edges.append(y)
    while y < 4.0:
        y = min(y + 0.05 / fine, 4.0)
        edges.append(y)
    return np.array(edges)


def nodes_weights(K, n=64, fine=1.0):
    x, w = np.polynomial.legendre.leggauss(n)
    E = panels(K, fine)
    ys, ws = [], []
    for a, b in zip(E[:-1], E[1:]):
        ys.append(0.5 * (b - a) * x + 0.5 * (b + a))
        ws.append(0.5 * (b - a) * w)
    return np.concatenate(ys), np.concatenate(ws)


def log_weight(y, K):
    N = 3 * K // 40
    u2 = (np.arange(1, K + 1) / K) ** 2
    t = y * y
    lo = 6 * np.sum(np.log(t[:, None] + u2[None, :N]), axis=1) - np.sum(np.log(t[:, None] + u2[None, :]), axis=1)
    Y = K * y
    x = np.exp(-2 * np.pi * Y)
    lw = math.log((2 * math.pi) ** 4 / 12) + 5 * np.log(Y) - 2 * np.pi * Y + np.log(1 + 11 * x + 11 * x * x + x ** 3) \
        - 5 * np.log(-np.expm1(-2 * np.pi * Y))
    return lo + lw


LOGS = 480 * math.log(2)  # Lanczos vectors are kept at norm S = 2^480 (extends float64 range; see notes)


def lanczos(tau, logomega, h):
    """float64 Lanczos with full reorthogonalisation on T = diag(tau), start sqrt(omega); vectors scaled to norm S
    so that entries down to S*e^-1490 (~ e^-1077) stay representable. Returns J coefficients and log pi_0."""
    L = logomega.max()
    lw = 0.5 * (logomega - L)                       # log sqrt(omega_i / max)
    lnrm = 0.5 * math.log(np.sum(np.exp(2 * lw)))    # log || sqrt(omega/max) ||
    logpi0 = -0.5 * L - lnrm
    S = math.exp(LOGS)
    v0 = np.exp(lw - lnrm + LOGS)                    # = S * sqrt(omega)/||sqrt(omega)||
    n_under = int(np.sum((lw - lnrm + LOGS) < -744))
    V = np.zeros((len(tau), h))
    V[:, 0] = v0
    al = np.zeros(h)
    be = np.zeros(h)  # be[k] couples k-1,k
    S2 = S * S
    for k in range(h):
        w = tau * V[:, k]
        al[k] = (V[:, k] @ w) / S2
        w -= al[k] * V[:, k]
        if k > 0:
            w -= be[k] * V[:, k - 1]
        for _ in range(2):
            w -= V[:, :k + 1] @ ((V[:, :k + 1].T @ w) / S2)
        if k + 1 < h:
            be[k + 1] = np.linalg.norm(w) / S
            V[:, k + 1] = w / be[k + 1]
    lanczos.n_under = n_under
    return al, be, logpi0, V / S


def r_coeffs(K):
    N = 3 * K // 40
    out = {}
    for j in range(N + 1, K + 1):
        num = Fraction(j ** 4, K ** 4)
        for l in range(1, N + 1):
            num *= Fraction(l * l - j * j, K * K) ** 5
        den = Fraction(1)
        for k in range(N + 1, K + 1):
            if k != j:
                den *= Fraction(k * k - j * j, K * K)
        out[j] = num / den
    return out


def norm_from_recurrence(K, al, be, logpi0, dps=40, iters=60, seed=0):
    mp.mp.dps = dps
    N, h = 3 * K // 40, 37 * K // 40
    r = r_coeffs(K)
    js = list(range(N + 1, K + 1))
    Pi = []
    alm = [mp.mpf(float(a)) for a in al]
    bem = [mp.mpf(float(b)) for b in be]
    p0 = mp.exp(mp.mpf(logpi0))
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
        lam = mp.fsum(Av[a] * v[a] for a in range(h)) / mp.fsum(v[a] ** 2 for a in range(h))
        nv = mp.sqrt(mp.fsum(z * z for z in Av))
        v = [z / nv for z in Av]
        if lam_old is not None and abs(lam - lam_old) < abs(lam) * mp.mpf(10) ** (-25):
            break
        lam_old = lam
    # Frobenius norm as a cross-check (upper bound)
    return lam, it


def run(K, n=64, fine=1.0, dps=40):
    t0 = time.time()
    N, h = 3 * K // 40, 37 * K // 40
    y, wq = nodes_weights(K, n, fine)
    lo = log_weight(y, K) + np.log(wq)
    tau = y * y
    al, be, logpi0, V = lanczos(tau, lo, h)
    orth = np.abs(V.T @ V - np.eye(h)).max()
    print("   [K=%d] nodes underflowing in start vector: %d ; log-weight range over nodes %.1f" % (
        K, lanczos.n_under, float(lo.max() - lo.min())), flush=True)
    t1 = time.time()
    lam, it = norm_from_recurrence(K, al, be, logpi0, dps=dps)
    return dict(K=K, nodes=len(y), orth=orth, lognorm=float(mp.log(abs(lam))), sign=int(mp.sign(lam)), iters=it,
                t_lanczos=t1 - t0, t_total=time.time() - t0)


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        K = int(arg)
        for (n, fine) in [(64, 1.0), (96, 1.5)]:
            res = run(K, n, fine)
            print("K=%d nodes=%d orth=%.1e : log||B_K|| = %.10f ; rate = %.8f ; sign %+d ; it %d ; [%.1fs/%.1fs]" % (
                K, res['nodes'], res['orth'], res['lognorm'], res['lognorm'] / K, res['sign'], res['iters'],
                res['t_lanczos'], res['t_total']), flush=True)
