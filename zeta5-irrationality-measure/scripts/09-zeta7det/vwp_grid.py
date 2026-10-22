"""Grid over Zudilin very-well-poised (r=5, q=9) Hankel determinants (07's vwp_hankel.forms):
L(P) = (1/4!) sum_{t>=1-h1} (R P(u))^{(4)}(t), u = t(t+h0); values in Q + Q zeta(7).
Record log P_n(zeta(7)) (primitive normalization) for n=1..nmax, fit c2 ~ log P / n^2 trend."""
import sys, json, time, itertools
sys.path.insert(0, './scripts/07-zeta')
import mpmath as mp
from vwp_hankel import forms, hankel_poly, content

out = sys.argv[1]
nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 3
r = 5
cases = []
for a in (1, 2):
    for b in (1, 2, 3):
        S = 5 * a + 4 * b
        for extra in (1, 2, 4, 6):
            eta0 = S // 2 + extra
            dmax = eta0 - S / 2
            for frac in (0.5, 1.0):
                kap = max(1, int(dmax * frac))
                cases.append((eta0, [a] * 5 + [b] * 4, kap))
for eta0, etas, kap in cases:
    rows = []
    for n in range(1, nmax + 1):
        d = kap * n
        t0 = time.time()
        try:
            A, B = forms(eta0, etas, r, n, 2 * d - 2)
        except AssertionError as e:
            rows.append(None); continue
        J = hankel_poly(A, B, d)
        if J == 0:
            rows.append('zero'); continue
        c = content(J)
        P = J / c
        cs = P.coeffs()
        mx = max(float(mp.log(abs(mp.mpf(int(x.p))))) for x in cs if x != 0)
        mp.mp.dps = int(mx / 2.3) + 200
        z = mp.zeta(7)
        v = mp.mpf(0)
        for x in reversed(cs):
            v = v * z + mp.mpf(int(x.p))
        lv = float(mp.log(abs(v))) if v != 0 else None
        rows.append(dict(n=n, d=d, logP=lv, log_cont=float(mp.log(mp.mpf(int(c.p))) - mp.log(int(c.q))), t=time.time() - t0))
    rec = dict(eta0=eta0, etas=etas, kappa=kap, res=rows)
    print(json.dumps(rec), flush=True)
    with open(out, 'a') as fh:
        fh.write(json.dumps(rec) + "\n")
