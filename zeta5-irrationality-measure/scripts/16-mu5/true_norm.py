"""16-mu5: numerical value of ||B_K|| (Fauzan App. C) for larger K, via orthonormal polynomials.  [NUMERICAL]

||B_K|| = K^5 sup_p |sum_j r_j p(x_j)^2| / int_0^oo p^2 W_K,   x_j = -(j/K)^2, N < j <= K   (notes/16-mu5.md, Lemma 1)
In an orthonormal basis pi_0..pi_{h-1} of L^2(W_K) this is the spectral norm of  B = K^5 Pi^T diag(r) Pi,
Pi_{ja} = pi_a(x_j).  W_K(t) dt = Omega(y^2) w(K y) dy with y = sqrt t, Omega(t) = prod_{l<=N}(t+u_l^2)^6 / prod_{j<=K}(t+u_j^2).
Discretization: composite Gauss-Legendre in y; Stieltjes procedure in mpmath; pi_a(x_j) by the three-term recurrence.
Output: (1/K) log ||B_K|| (eigenvalues via mpmath for moderate h) and (1/K) log ||B||_F.
"""
import sys, math, time, json
import mpmath as mp


def gl_nodes(n, dps):
    with mp.workdps(dps + 10):
        xs, ws = [], []
        # Gauss-Legendre on [-1,1] via mpmath
        import mpmath.calculus.quadrature as q
        gl = q.GaussLegendre(mp.mp)
        # degree m gives 3*2^(m-1) nodes; pick m with enough nodes
        m = 1
        while 3 * 2 ** (m - 1) < n:
            m += 1
        nodes = gl.calc_nodes(m, mp.mp.prec)
        return [(mp.mpf(x), mp.mpf(w)) for x, w in nodes]


def run(K, dps=40, panels=None, nper=48, eig=True):
    mp.mp.dps = dps
    N = 3 * K // 40
    h = K - N
    u = [mp.mpf(j) / K for j in range(K + 1)]
    if panels is None:
        panels = [mp.mpf(i) / 400 for i in range(0, 41)] + [mp.mpf(1) / 10 + mp.mpf(i) / 50 for i in range(1, 71)]
    base = gl_nodes(nper, dps)
    ys, ws = [], []
    for lo, hi in zip(panels[:-1], panels[1:]):
        c, d = (lo + hi) / 2, (hi - lo) / 2
        for x, w in base:
            ys.append(c + d * x)
            ws.append(d * w)
    # log weights
    lw = []
    for y, w in zip(ys, ws):
        t = y * y
        lo = 6 * mp.fsum(mp.log(t + u[l] ** 2) for l in range(1, N + 1)) - mp.fsum(mp.log(t + u[j] ** 2) for j in range(1, K + 1))
        Y = K * y
        q = mp.exp(-2 * mp.pi * Y)
        ww = (2 * mp.pi) ** 4 * Y ** 5 / 12 * q * (1 + 11 * q + 11 * q * q + q ** 3) / (1 - q) ** 5
        lw.append(lo + mp.log(ww) + mp.log(w))
    lmax = max(lw)
    Wn = [mp.exp(x - lmax) for x in lw]          # normalized discrete weights (overall factor e^lmax)
    ts = [y * y for y in ys]
    # Stieltjes
    a, b = [], []
    p_prev = [mp.mpf(0)] * len(ts)
    nrm0 = mp.sqrt(mp.fsum(Wn))
    p = [1 / nrm0] * len(ts)
    bprev = mp.mpf(0)
    for k in range(h):
        ak = mp.fsum(w * t * x * x for w, t, x in zip(Wn, ts, p))
        a.append(ak)
        if k == h - 1:
            break
        q = [(t - ak) * x - bprev * xp for t, x, xp in zip(ts, p, p_prev)]
        bk = mp.sqrt(mp.fsum(w * x * x for w, x in zip(Wn, q)))
        b.append(bk)
        p_prev, p = p, [x / bk for x in q]
        bprev = bk
    # orthonormal w.r.t. sum Wn, i.e. w.r.t. W_K up to the factor e^{lmax}: int p^2 W_K = e^{lmax} sum p^2 Wn
    # evaluate at x_j
    Pi = []
    for j in range(N + 1, K + 1):
        x = -u[j] ** 2
        vals = [1 / nrm0]
        prev = mp.mpf(0)
        for k in range(h - 1):
            nxt = ((x - a[k]) * vals[-1] - (b[k - 1] if k > 0 else 0) * prev) / b[k]
            prev = vals[-1]
            vals.append(nxt)
        Pi.append(vals)
    # r_j closed form
    r = []
    for j in range(N + 1, K + 1):
        sgn = (-1) ** (N + (j - N - 1))
        lr = mp.log(2) + (2 * K - 12 * N - 6) * mp.log(K) + 6 * mp.log(mp.factorial(j + N)) - 6 * mp.log(mp.factorial(j - N - 1)) \
            - mp.log(mp.factorial(K - j)) - mp.log(mp.factorial(K + j))
        r.append(sgn * mp.exp(lr))
    # B = K^5 e^{-lmax} Pi^T diag(r) Pi
    scale = mp.mpf(K) ** 5 * mp.exp(-lmax)
    Bm = mp.matrix(h, h)
    for aa in range(h):
        for bb in range(aa, h):
            s = mp.fsum(r[i] * Pi[i][aa] * Pi[i][bb] for i in range(h))
            Bm[aa, bb] = Bm[bb, aa] = scale * s
    fro = mp.sqrt(mp.fsum(Bm[i, j] ** 2 for i in range(h) for j in range(h)))
    out = dict(K=K, h=h, logFro=float(mp.log(fro)), rateFro=float(mp.log(fro)) / K)
    if eig:
        ev = mp.eigsy(Bm, eigvals_only=True)
        mx = max(abs(e) for e in ev)
        out.update(logNorm=float(mp.log(mx)), rate=float(mp.log(mx)) / K,
                   npos=sum(1 for e in ev if e > 0), nneg=sum(1 for e in ev if e < 0))
    return out


if __name__ == "__main__":
    Ks = [int(x) for x in sys.argv[1].split(',')]
    eig = (sys.argv[2] != 'noeig') if len(sys.argv) > 2 else True
    dps = int(sys.argv[3]) if len(sys.argv) > 3 else 40
    for K in Ks:
        t0 = time.time()
        res = run(K, dps=dps, eig=eig)
        res['time'] = round(time.time() - t0, 1)
        print(json.dumps(res), flush=True)
