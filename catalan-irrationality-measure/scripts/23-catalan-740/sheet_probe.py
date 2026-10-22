"""
Probe how loose OpenAI's kappa=2 sheet bound |R_n(x)| <= e^{O(n log n)} V(1/x) prod|x|^g is
at configurations distributed like the kappa=2 equilibrium measure.            [NUMERICAL]

R_n(x) = det[ xi_f(x_i) x_i^{g-r} + xi_n(x_i) x_i^{r-g} ]_{r<n, i<=n}
 (xi_f,xi_n) = (1/2,1/2) for x<0 and (-1/4,5/4) for x>0   (OpenAI (65)).
usage: sheet_probe.py n0list    (design (20,3,1,2,1); g = n/10)
"""
import sys, json
import numpy as np
import mpmath as mp
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from laplace_explore import build


def equilibrium_quantiles(n, Mx=160, Ms=100):
    P, (ex, exp_, es) = build(1.0, Mx, Ms)
    val, z = P.solve(iters=3000, step=0.5, polish=True)
    w = z[P.blk(0)]; t = z[P.blk(1)]
    # theta=1: sigma+tau is the kappa=2 measure; merge sigma (on (0,1) grid) into the x grid
    xm = (ex[:-1] + ex[1:]) / 2
    # x grid positive part coincides with exp_ grid
    npos = len(exp_) - 1
    mu = t.copy()
    mu[-npos:] += w
    cdf = np.concatenate([[0], np.cumsum(mu)])
    cdf /= cdf[-1]
    q = (np.arange(n) + 0.5) / n
    xs = np.interp(q, cdf, ex)
    return xs, val


def logR(xs, g, dps):
    mp.mp.dps = dps
    n = len(xs)
    M = mp.matrix(n, n)
    for i, x in enumerate(xs):
        x = mp.mpf(x)
        xf, xn = (mp.mpf(1) / 2, mp.mpf(1) / 2) if x < 0 else (mp.mpf(-1) / 4, mp.mpf(5) / 4)
        for r in range(n):
            M[r, i] = xf * x ** (g - r) + xn * x ** (r - g)
    d = mp.det(M)
    return float(mp.log(abs(d)))


def log_far(xs, g):
    xs = np.asarray(xs, float)
    n = len(xs)
    s = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            s += np.log(abs(1 / xs[i] - 1 / xs[j]))
    return s + g * np.log(np.abs(xs)).sum()


if __name__ == "__main__":
    ns = [int(v) for v in sys.argv[1:]] or [20, 40, 60, 80]
    for n in ns:
        g = n // 10
        xs, val = equilibrium_quantiles(n)
        dps = int(n * n * 0.5) + 50
        lr = logR(xs, g, dps)
        lf = log_far(xs, g)
        print(json.dumps(dict(n=n, g=g, logR=lr, logFar=lf, diff=lr - lf, diff_n2=(lr - lf) / n ** 2,
                              npos=int((xs > 0).sum()))), flush=True)
