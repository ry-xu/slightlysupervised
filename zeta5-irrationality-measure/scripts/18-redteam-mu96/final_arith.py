"""18-redteam-mu96: arithmetic constant A_M, decay eps, and the final mu threshold (exact rationals + arb).

Inputs (each re-derived/re-certified in this directory unless stated):
  I_out = 127751/96000                       (5.10)  [Fauzan; 15 recomputed exactly]
  int_3^20 R/x^3 = (5.18)                    [Fauzan App. B; 15 recomputed exactly]
  int_20^4000 R/x^3 (exact rational)         [09 and 15, two independent exact computations; file from 15]
  int_4000^oo R/x^3 <= -lam/4000 + (3 lam/(2 al) + 13/16)/4000^2   [re-derived here, see note]
  int_M^oo R/x^3 >= -lam/M + (2923/240 - 1/4)/M^2 - 32/M^3  for 40 | M     (5.17)
  (5.11): limsup K^-2 log m_{K,M} <= I_out + 6 lam/M + int_3^M R/x^3        (Prop 5.2, any fixed M >= 40)
  U <= lam*M0 - I(rho) + C*  with M0 = -6.6607665 (certified in run_real_certs.py)
  rho_2 = M1c - m0 = 0.808802 + 6.660814 (certified in run_contour_cert.py / run_real_certs.py)
"""
import sys
from flint import fmpq, arb, ctx
sys.path.insert(0, './scripts/18-redteam-mu96')
from common import load_measure, energy_A2, Cstar, ALPHA, LAM, ROOT

ctx.prec = 200
al, lam = ALPHA, LAM
Iout = fmpq(127751, 96000)
I318 = fmpq(322437603634266857629, 7535670527041937280000)
num, den = open(ROOT + "/data/15-redteam-mu5/int_20_4000.txt").read().split()
I20_4000 = fmpq(int(num), int(den))
tail4000 = -lam / 4000 + (3 * lam / (2 * al) + fmpq(13, 16)) / 4000 ** 2
Astar = Iout + I318 + I20_4000 + tail4000


def A_M(M):
    assert M % 40 == 0 and M > 4000
    return Astar + 7 * lam / M - (fmpq(2923, 240) - fmpq(1, 4)) / M ** 2 + fmpq(32) / M ** 3


print("int_20^4000 = %.16f ; tail bound at 4000 = %.12e" % (float(I20_4000), float(tail4000)))
print("A*          = %.16f" % float(Astar))
A1e5 = A_M(100000)
A4e8 = A_M(400000000)
print("A_{1e5}     = %.16f   (15: 1.31351917198018)" % float(A1e5))
print("A_{4e8}     = %.16f   (16: 1.313454439360566)" % float(A4e8))
print("gain A_1e5 - A_4e8 = %.6e" % float(A1e5 - A4e8))

A, B, C = load_measure()
M0 = fmpq(-66607665, 10 ** 7)
I = energy_A2(A, B, C)
U = arb(lam) * arb(M0) - I + Cstar()
print("U <= %s" % U.upper())
eps_max = -(arb(A4e8) + U)
print("eps_max = -(A_M + U) = %s" % eps_max)
eps_09 = -(arb(A1e5) + U)
print("(same U, M=1e5): eps = %s" % eps_09)
rho2 = fmpq(808802, 10 ** 6) - fmpq(-6660814, 10 ** 6)
print("rho_2 = %s = %s" % (rho2, float(rho2)))
thr = arb(rho2) * arb(lam) / eps_max
print("threshold rho_2*lam/eps_max = %s" % thr)
# 16's explicit choices
eps = fmpq(9005205399, 125 * 10 ** 9)
c = fmpq(513596269, 40000000)
mu = fmpq(9591, 100)
print("eps = %s < eps_max : %s" % (eps, bool(arb(eps) < eps_max)))
print("eps*c - lam = %s (>0: %s)" % (float(eps * c - lam), eps * c - lam > 0))
print("mu - rho2*c = %s (>0: %s)" % (float(mu - rho2 * c), mu - rho2 * c > 0))
# admissible-c interval for mu = 95.91: lam/eps < c < mu/rho2
print("c-interval for mu=95.91: (%.9f, %.9f)" % (float(lam / eps), float(mu / rho2)))
# smallest 2-decimal mu
for m in [fmpq(9590, 100), fmpq(95909, 1000), fmpq(9591, 100)]:
    ok = arb(m) > thr
    print("mu = %s admissible (mu > threshold): %s" % (float(m), bool(ok)))
# sensitivity
d_eps = fmpq(1, 10 ** 6)
print("d threshold / d eps ~ %.3f per 1e-6 in eps" % float(arb(rho2) * arb(lam) * (1 / (eps_max - arb(d_eps)) - 1 / eps_max)).__float__() if False else "")
print("margin of 95.91 over threshold: %s" % (arb(mu) - thr))
