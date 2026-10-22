"""Test hypotheses for the n^2-size of R_n at equilibrium-quantile configurations. [NUMERICAL]"""
import sys, json
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from sheet_probe import equilibrium_quantiles, log_far
from sheet_probe2 import logR_arb

def pairsum(xs, f):
    s = 0.0
    n = len(xs)
    for i in range(n):
        for j in range(i + 1, n):
            s += f(xs[i], xs[j])
    return s

for n in [int(v) for v in sys.argv[1:]]:
    xs, val = equilibrium_quantiles(n)
    for gfrac in [0.0, 0.1]:
        g = int(round(gfrac * n))
        lr, rad = logR_arb(xs, g)
        lf = log_far(xs, g)
        neg = xs[xs < 0]; pos = xs[xs > 0]
        Jnn = pairsum(neg, lambda a, b: np.log(1 - a * b))
        Jpp = pairsum(pos, lambda a, b: np.log(1 - a * b))
        Jall = pairsum(xs, lambda a, b: np.log(1 - a * b))
        Jnp = Jall - Jnn - Jpp
        print(json.dumps(dict(n=n, g=g, d=(lr - lf) / n**2, Jnn=Jnn / n**2, Jpp=Jpp / n**2, Jnp=Jnp / n**2,
                              resid=(lr - lf - Jnn) / n**2)), flush=True)
