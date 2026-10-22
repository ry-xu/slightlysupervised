"""Quadrature check: mu_{zeta(k)}(W t^m) == int_0^oo W(y^2) y^{2m} w_k(y) dy for small weights."""
import sys
sys.path.insert(0, './scripts/09-zeta7det')
import mpmath as mp
from zdet import entries, fauzan_exps

mp.mp.dps = 30


def wk(k, y):
    s = mp.polylog(1 - k, mp.e ** (-2 * mp.pi * y))
    return 2 * (2 * mp.pi) ** (k - 1) * y ** k * s / mp.factorial(k - 1)


for k in (3, 5, 7, 9):
    for (K, N, r) in ((3, 1, 2), (4, 1, 3), (5, 2, 2)):
        e = fauzan_exps(K, N, r)
        a, b = entries(k, e, 3)
        W = lambda t: mp.fprod([(t + l * l) ** v for l, v in e.items()])
        z = mp.zeta(k)
        errs = []
        for m in range(5):
            I = mp.quad(lambda y: W(y * y) * y ** (2 * m) * wk(k, y), [0, 0.5, 2, 6, mp.inf])
            pred = mp.mpf(int(a[m].p)) / int(a[m].q) + mp.mpf(int(b[m].p)) / int(b[m].q) * z
            errs.append(float(abs(I - pred) / abs(I)))
        print(f"k={k} K={K} N={N} r={r}: max rel err {max(errs):.2e}")
