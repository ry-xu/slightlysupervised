"""Float exploration (not a certificate): landscape of f = 2U^rho - V on [0,2] and H = 2U^rho - Vt on Gamma."""
import sys, math, time
import numpy as np
sys.path.insert(0, '.')
from common import *

A, B, C = load_measure()
adm = admissibility(A, B, C)
print("admissibility:", adm)
fm = FloatMeasure(A, B, C)
lo, hi = float(A[-1]), float(B[-1])
f = lambda t: 2 * fm.U(t) - V_float(t)
# grid: dense inside S near edges, plus outside
ts = np.unique(np.concatenate([np.linspace(0, 2, 20001), np.geomspace(1e-9, 0.01, 20001),
                               np.linspace(0.70, 0.80, 100001), np.linspace(lo, hi, 200001)]))
t0 = time.time()
vals = np.array([f(t) for t in ts])
print("grid eval %.1fs, n=%d" % (time.time() - t0, len(ts)))
i = np.argmax(vals)
print("sup f on grid: %.12f at t=%.9f" % (vals[i], ts[i]))
inS = (ts >= lo) & (ts <= hi)
j = np.argmin(np.where(inS, vals, np.inf))
print("inf_S f on grid: %.12f at t=%.9f" % (vals[j], ts[j]))
print("max_S f on grid: %.12f" % np.max(vals[inS]))
# local refinement of sup
from scipy.optimize import minimize_scalar
r = minimize_scalar(lambda t: -f(t), bounds=(max(0, ts[i] - 2e-3), min(2, ts[i] + 2e-3)), method='bounded', options={'xatol': 1e-13})
print("refined sup: %.13f at t=%.12f" % (-r.fun, r.x))
print("claimed M0 = -6.6607665  slack = %.3e" % (-6.6607665 - (-r.fun)))
print("claimed m0 = -6.660814   slack = %.3e" % (vals[j] - (-6.660814)))
# contour
x1 = -1.069911
Hf = lambda z: 2 * fm.Uc(complex(z)) - Vt_float(z)
print("H(x1) = %.15f" % Hf(x1))
# real-axis minimum left of -1
xs = np.linspace(-1.2, -1.0001, 4001)
hv = np.array([Hf(x) for x in xs])
k = np.argmin(hv)
print("min_{x<-1} H on grid: %.12f at %.7f" % (hv[k], xs[k]))
r2 = minimize_scalar(lambda x: Hf(x), bounds=(xs[k] - 1e-3, xs[k] + 1e-3), method='bounded', options={'xatol': 1e-12})
print("refined saddle: x=%.10f H=%.15f" % (r2.x, r2.fun))
# along Gamma+
ys = np.linspace(0, 0.2, 2001)
left = np.array([Hf(complex(x1, y)) for y in ys])
print("left edge: max %.12f at y=%.5f; monotone decreasing: %s" % (left.max(), ys[np.argmax(left)], bool(np.all(np.diff(left) < 0))))
xs2 = np.linspace(x1, 1.0, 4001)
top = np.array([Hf(complex(x, 0.2)) for x in xs2])
print("top edge: max %.12f at x=%.5f" % (top.max(), xs2[np.argmax(top)]))
right = np.array([Hf(complex(1.0, y)) for y in ys])
print("right edge: max %.12f at y=%.5f" % (right.max(), ys[np.argmax(right)]))
