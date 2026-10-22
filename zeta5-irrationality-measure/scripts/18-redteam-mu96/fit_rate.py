"""18-redteam-mu96: asymptotic rate of log||B_K|| from own computations (exact K=40,80; Lanczos K<=320; mp K=480,640)."""
import numpy as np
data = {40: 292.5613171408397, 80: 591.0372993642756, 160: 1188.2728511251, 240: 1785.6413312423,
        320: 2383.0664002111, 480: 3577.9992183777, 640: 4772.9897927749, 960: 7163.0547545318,
        1280: 9553.1780391228}
Ks = np.array(sorted(data)); L = np.array([data[k] for k in Ks])
print("K, log||B_K||, rate:")
for k, l in zip(Ks, L):
    print("  %4d  %.10f  %.8f" % (k, l, l / k))
def fit(sel, fix_c=None):
    k = Ks[sel].astype(float); l = L[sel]
    if fix_c is None:
        A = np.vstack([k, -np.log(k), -np.ones_like(k)]).T
        sol, *_ = np.linalg.lstsq(A, l, rcond=None)
        res = l - A @ sol
        return sol, np.abs(res).max()
    A = np.vstack([k, -np.ones_like(k)]).T
    sol, *_ = np.linalg.lstsq(A, l + fix_c * np.log(k), rcond=None)
    res = l + fix_c * np.log(k) - A @ sol
    return sol, np.abs(res).max()
for name, sel in [("K>=160", Ks >= 160), ("K>=240", Ks >= 240), ("K>=320", Ks >= 320), ("K>=480", Ks >= 480), ("K>=640", Ks >= 640)]:
    sol, r = fit(sel)
    print("fit R K - c log K - d  [%s]: R = %.6f, c = %.3f, d = %.3f, max resid %.2e" % (name, sol[0], sol[1], sol[2], r))
for c in [0.5]:
    for name, sel in [("K>=320", Ks >= 320), ("K>=640", Ks >= 640), ("K>=960", Ks >= 960)]:
        sol, r = fit(sel, fix_c=c)
        print("fit with c=%.1f fixed [%s]: R = %.6f, d = %.3f, max resid %.2e" % (c, name, sol[0], sol[1], r))
# successive differences: (L(K2)-L(K1))/(K2-K1) + c*(log K2 - log K1)/(K2-K1) with c=1/2
print("difference quotients (c=1/2 corrected):")
for i in range(len(Ks) - 1):
    k1, k2 = Ks[i], Ks[i + 1]
    dq = (L[i + 1] - L[i]) / (k2 - k1) + 0.5 * (np.log(k2) - np.log(k1)) / (k2 - k1)
    print("  [%d,%d]: %.6f" % (k1, k2, dq))
rho2_cert = 466851 / 62500
print("certified Theorem-2 rate rho_2 = %.6f ; deficit (rho_2 K - log||B||) by K:" % rho2_cert,
      ["%d:%.3f" % (k, rho2_cert * k - l) for k, l in zip(Ks, L)])
