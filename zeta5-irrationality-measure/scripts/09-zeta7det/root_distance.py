"""Distance from zeta(k) to the nearest root of the exact Hankel polynomial Delta_K(X):
||B_K|| = 1/dist  (Fauzan App. C: Q_K(X) = Q_K(xi) prod (1+(X-xi) beta_j)), so rho_true = log(1/dist)/K."""
import sys, json, math
sys.path.insert(0, './scripts/09-zeta7det')
import mpmath as mp
from zgen import Family, run_gen, fauzan_like

k = int(sys.argv[1]); Km, Nm, r, hm = float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
ns = [int(x) for x in sys.argv[6].split(',')]
fam = Family('zeta', k)
for n in ns:
    K, N, h = int(round(Km * n)), int(round(Nm * n)), int(round(hm * n))
    res, D, P = run_gen(fam, fauzan_like(fam, K, N, r), h)
    cs = [int(c.p) for c in P.coeffs()]
    mx = max(abs(c).bit_length() for c in cs if c)
    mp.mp.dps = int(mx * 0.30103) + 2000
    xi = mp.zeta(k)
    poly = lambda x: mp.polyval(list(reversed(cs)), x)
    # roots are real (B_K symmetric); find nearest sign change of P around xi on a log grid, then bisect
    P0 = poly(xi)
    def side(sgn):
        lo = None
        for m in range(1, 4000):
            dlt = mp.mpf(10) ** (-m / 4.0)
            pass
        # scan from tiny to large delta
        prev = mp.mpf(10) ** (-(mp.mp.dps // 3))
        if mp.sign(poly(xi + sgn * prev)) != mp.sign(P0):
            return prev
        e = prev
        while e < 10:
            e2 = e * 2
            if mp.sign(poly(xi + sgn * e2)) != mp.sign(P0):
                a, b = e, e2
                for _ in range(60):
                    m = (a + b) / 2
                    if mp.sign(poly(xi + sgn * m)) != mp.sign(P0):
                        b = m
                    else:
                        a = m
                return (a + b) / 2
            e = e2
        return mp.inf
    d = min(side(1), side(-1))
    nreal = None
    out = dict(k=k, K=K, h=h, logP=res['logP'], log_dist=float(mp.log(d)), rho=float(-mp.log(d)) / K,
               eps=-res['logP'] / K ** 2, n_real_roots=nreal, deg=len(cs) - 1)
    print(json.dumps(out), flush=True)
