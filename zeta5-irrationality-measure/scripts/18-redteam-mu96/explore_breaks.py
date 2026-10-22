"""Float: f = 2U - V at all breakpoints and gap-wise sampled extrema (exploration only)."""
import sys, time
import numpy as np
sys.path.insert(0, '.')
from common import *
A, B, C = load_measure()
fm = FloatMeasure(A, B, C)
f = lambda t: 2 * fm.U(t) - V_float(t)
bps = np.sort(np.array([float(x) for x in A] + [float(x) for x in B]))
fb = np.array([f(t) for t in bps])
print("at breakpoints: min %.12f at %.9f ; max %.12f at %.9f" % (fb.min(), bps[fb.argmin()], fb.max(), bps[fb.argmax()]))
# per gap, sample 41 points
mins, maxs = [], []
t0 = time.time()
for k in range(len(bps) - 1):
    l, r = bps[k], bps[k + 1]
    ts = np.linspace(l, r, 41)[1:-1]
    v = np.array([f(t) for t in ts])
    mins.append((v.min(), ts[v.argmin()])); maxs.append((v.max(), ts[v.argmax()]))
mn = min(mins); mx = max(maxs)
print("gap samples: min %.12f at %.9f ; max %.12f at %.9f  (%.1fs)" % (mn[0], mn[1], mx[0], mx[1], time.time() - t0))
# distribution of where the minimum of f in S occurs
print("inf_S estimate: %.12f ; slack to m0=-6.660814: %.3e" % (min(fb.min(), mn[0]), min(fb.min(), mn[0]) + 6.660814))
print("sup estimate (in S): %.12f ; slack to M0=-6.6607665: %.3e" % (max(fb.max(), mx[0]), -6.6607665 - max(fb.max(), mx[0])))
