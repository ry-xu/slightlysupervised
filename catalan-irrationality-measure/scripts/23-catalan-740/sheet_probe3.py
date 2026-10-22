"""Sheet-bound looseness vs g at a fixed configuration (kappa=2 equilibrium quantiles of (20,3,1,2,1)).
[NUMERICAL]"""
import sys, json
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from sheet_probe import equilibrium_quantiles, log_far
from sheet_probe2 import logR_arb
n = int(sys.argv[1])
xs, val = equilibrium_quantiles(n)
for gfrac in [0.0, 0.05, 0.1, 0.15, 0.2]:
    g = int(round(gfrac * n))
    lr, rad = logR_arb(xs, g)
    lf = log_far(xs, g)
    print(json.dumps(dict(n=n, g=g, diff=lr - lf, diff_n2=(lr - lf) / n ** 2)), flush=True)
