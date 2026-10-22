"""
Same probe as sheet_probe.py with arb matrices (faster, rigorous balls).   [NUMERICAL]
Configurations: quantiles of the kappa=2 equilibrium measure (design (20,3,1,2,1)).
Reports log|R_n| - log(V(1/x) prod|x|^g) per n^2.
"""
import sys, json
import numpy as np
from flint import arb, arb_mat, ctx
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from sheet_probe import equilibrium_quantiles, log_far


def logR_arb(xs, g):
    n = len(xs)
    ctx.prec = int(n * n * 1.2) + 200
    rows = []
    X = [arb(float(x)) for x in xs]
    M = arb_mat(n, n)
    for i, x in enumerate(X):
        if xs[i] < 0:
            xf, xn = arb(1) / 2, arb(1) / 2
        else:
            xf, xn = arb(-1) / 4, arb(5) / 4
        p = x ** (-g)          # x^{-g}
        inv = 1 / x
        pinv = inv ** (-g)     # x^{g}
        for r in range(n):
            M[r, i] = xf * pinv + xn * p
            p = p * x
            pinv = pinv * inv
    d = M.det()
    return float(abs(d).log().mid()), d.rad()

if __name__ == "__main__":
    ns = [int(v) for v in sys.argv[1:]]
    for n in ns:
        g = n // 10
        xs, val = equilibrium_quantiles(n)
        lr, rad = logR_arb(xs, g)
        lf = log_far(xs, g)
        print(json.dumps(dict(n=n, g=g, logR=lr, logFar=lf, diff=lr - lf, diff_n2=(lr - lf) / n ** 2,
                              npos=int((xs > 0).sum()))), flush=True)
