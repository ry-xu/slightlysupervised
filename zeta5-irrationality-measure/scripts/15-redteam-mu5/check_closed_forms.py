"""15-redteam-mu5: check the closed forms used in the certificate against their defining integrals (mpmath quad).
(A.1)  U^{omega[a,b]}(t) = int log|t-u| / (pi sqrt((u-a)(b-u))) du
(A.5)  V(t) = 2 pi sqrt t + int_0^1 log(t+u^2) du - 6 int_0^al log(t+u^2) du"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from measure_common import U_arc_f, V_f
mp.mp.dps = 30
a, b = mp.mpf('0.003906748086'), mp.mpf('0.008992695531')
worst = 0
for t in [0, 0.001, 0.0039, 0.005, 0.0089927, 0.02, 0.3, 1.7]:
    t = mp.mpf(t)
    # substitute u = m + r cos(theta): d omega = d theta / pi
    m, r = (a + b) / 2, (b - a) / 2
    f = lambda th: mp.log(abs(t - m - r * mp.cos(th))) / mp.pi
    pts = [0, mp.pi]
    if a < t < b:
        pts = [0, mp.acos((t - m) / r), mp.pi]
    val = mp.quad(f, pts)
    worst = max(worst, abs(val - U_arc_f(float(t), float(a), float(b))))
print("max |(A.1) closed form - quadrature| =", mp.nstr(worst, 3))
al = mp.mpf(3) / 40
worst = 0
for t in [0, 1e-6, 0.00592, 0.1, 0.75, 2]:
    t = mp.mpf(t)
    val = 2 * mp.pi * mp.sqrt(t) + mp.quad(lambda u: mp.log(t + u * u), [0, al, 1]) \
        - 6 * mp.quad(lambda u: mp.log(t + u * u), [0, al])
    worst = max(worst, abs(val - V_f(float(t))))
print("max |(A.5) closed form - quadrature| =", mp.nstr(worst, 3))
