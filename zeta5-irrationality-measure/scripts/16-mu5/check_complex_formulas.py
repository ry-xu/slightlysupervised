"""16-mu5: independent check (mpmath quadrature) of the complex closed forms used in contour_certify.py:
   Vt(tau) = int_0^1 log|tau+u^2| du - 6 int_0^al log|tau+u^2| du,
   U^{omega[a,b]}(tau) = int log|tau-u| d omega[a,b](u)  (arcsine),
   F'(tau) components: int_0^a du/(tau+u^2) and int d omega(u)/(tau-u)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
from fractions import Fraction as Fr
from flint import acb, arb
from certlib import A, alA
from contour_certify import Vt, dVt, C

mp.mp.dps = 30
al = mp.mpf(3) / 40
pts = [(-1.069911, 0.0), (-1.069911, 0.07), (-0.5, 0.2), (-0.98, 0.05), (0.3, 0.2), (1.0, 0.1), (-0.003, 0.2)]
worst = 0
for x, y in pts:
    tau = mp.mpc(x, y)
    q = lambda a: mp.quad(lambda u: mp.log(abs(tau + u * u)), mp.linspace(0, a, 9))
    Vq = q(1) - 6 * q(al)
    Va = Vt(C(Fr(x).limit_denominator(10**9), Fr(y).limit_denominator(10**9)))
    d1 = abs(Vq - mp.mpf(float(Va.mid())))
    # derivative check
    dq = mp.quad(lambda u: 1 / (tau + u * u), [0, 1]) - 6 * mp.quad(lambda u: 1 / (tau + u * u), [0, al])
    da = dVt(C(Fr(x).limit_denominator(10**9), Fr(y).limit_denominator(10**9)))
    d2 = abs(dq - mp.mpc(float(da.real.mid()), float(da.imag.mid())))
    worst = max(worst, d1)
    print(f"tau={x}+{y}i: Vt quad = {mp.nstr(Vq, 15)}, closed = {Va.mid().str(15)}, |diff| = {mp.nstr(d1, 3)};  dVt |diff| = {mp.nstr(d2, 3)}")
# arcsine potential and Stieltjes transform at a complex point
a_, b_ = mp.mpf('0.1'), mp.mpf('0.6')
m, r = (a_ + b_) / 2, (b_ - a_) / 2
for (x, y) in [(-1.0699, 0.0), (-0.5, 0.2), (0.3, 0.2), (1.0, 0.1)]:
    tau = mp.mpc(x, y)
    dens = lambda th: 1 / mp.pi  # u = m + r cos th, d omega = dth/pi
    Uq = mp.quad(lambda th: mp.log(abs(tau - (m + r * mp.cos(th)))) / mp.pi, [0, mp.pi])
    z = (tau - m) / r
    w = z + mp.sqrt(z - 1) * mp.sqrt(z + 1)
    Uc = mp.log(r / 2) + mp.log(abs(w))
    Sq = mp.quad(lambda th: 1 / (tau - (m + r * mp.cos(th))) / mp.pi, [0, mp.pi])
    Sc = 1 / (mp.sqrt(tau - a_) * mp.sqrt(tau - b_))
    print(f"arcsine at {x}+{y}i: U quad-closed = {mp.nstr(abs(Uq - Uc), 3)};  Stieltjes quad-closed = {mp.nstr(abs(Sq - Sc), 3)}")
