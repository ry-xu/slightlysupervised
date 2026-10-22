"""15-redteam-mu5: consistency sandwich for the analytic constant (float, not a certificate).
For any positive rho of mass lam:  E(rho) := I(rho) - int V drho  <=  max_rho E = U* - C*  <=  J(rho) := lam*sup(2U-V) - I(rho).
So E(rho)+C* <= U* <= (certified U).  A violation would indicate a bug."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from measure_common import load_09, load_fauzan, V_f
lam = 37 / 40
Cs = -2 * lam + 12 * (3/40) * lam * (1 - math.log(3/40)) + 3 * lam ** 2 - 2 * lam ** 2 * math.log(2 * lam)
def E(a, b, c, n=20000):
    k = np.arange(1, n + 1); cosk = np.cos((2 * k - 1) * math.pi / (2 * n))
    S = 0.0; I = 0.0; IV = 0.0
    for aj, bj, cj in zip(a, b, c):
        aj, bj, cj = float(aj), float(bj), float(cj)
        S1 = S + cj; I += (S1 * S1 - S * S) * math.log((bj - aj) / 4); S = S1
        pts = (aj + bj) / 2 + (bj - aj) / 2 * cosk
        IV += cj * np.mean([V_f(t) for t in pts[::10]]) if n > 5000 else cj * np.mean([V_f(t) for t in pts])
    return I - IV
for name, m, Ucert in [('fauzan', load_fauzan(), -1.3669955645), ('09', load_09('3p'), -1.3843989249)]:
    e = E(*m)
    print(f"{name}: E(rho)+C* = {e + Cs:.7f}  <=  U*  <=  certified {Ucert:.7f}   consistent: {e + Cs <= Ucert}")
