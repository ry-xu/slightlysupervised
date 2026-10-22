"""Calibration of the sheet-looseness on simple configurations (g=0 and g=n/10). [NUMERICAL]
 A: all negative, uniform quantiles of (-0.9,-0.1)
 B: all positive, uniform quantiles of (0.1,0.9)
 C: half/half: negatives uniform (-0.9,-0.1), positives uniform (0.1,0.9)
Prints d = (log|R_n| - log far)/n^2 and Jall = sum_{i<j} log(1-x_i x_j)/n^2, Jnn, Jpp."""
import sys, json
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from sheet_probe import log_far
from sheet_probe2 import logR_arb

def pairsum(xs):
    xs = np.asarray(xs); M = np.log(1 - np.outer(xs, xs)); 
    return (M.sum() - np.trace(M)) / 2

def quant(a, b, m):
    return a + (b - a) * (np.arange(m) + 0.5) / m

for n in [int(v) for v in sys.argv[1:]]:
    confs = {'A': quant(-0.9, -0.1, n), 'B': quant(0.1, 0.9, n),
             'C': np.concatenate([quant(-0.9, -0.1, n // 2), quant(0.1, 0.9, n - n // 2)])}
    for name, xs in confs.items():
        for g in [0, n // 10]:
            if name == 'A' and g > 0:
                continue
            lr, rad = logR_arb(xs, g)
            lf = log_far(xs, g)
            neg = xs[xs < 0]; pos = xs[xs > 0]
            print(json.dumps(dict(conf=name, n=n, g=g, d=(lr - lf) / n**2, Jall=pairsum(xs) / n**2,
                                  Jnn=pairsum(neg) / n**2 if len(neg) > 1 else 0, Jpp=pairsum(pos) / n**2 if len(pos) > 1 else 0)), flush=True)
