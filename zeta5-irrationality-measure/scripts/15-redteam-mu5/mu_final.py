"""15-redteam-mu5: final arithmetic for mu(zeta(5)) from the independently certified constants.

Framework (Fauzan App. C): pick eps < -(A_M + U) [M = 100000], rho_ > kappa (C.6), and c with
     eps*c > lam   and   rho_*c < mu.     Then |zeta(5) - a/b| > b^-mu for all large b.
So any mu > rho_*lam/eps works.  All comparisons below are done with exact rationals or certain arb comparisons.
Inputs (upper bounds, from my own runs):
   U        <= -1.38439892493401715650888837988593764 + 5.3e-36     (u_certify.py, M0 = -6.6600532)
   A_100000 <= 1.31351917198018 + 1.5e-15                            (final_A.py)
   A_200    <= 1.34553519400640 + 6e-16
"""
from fractions import Fraction as Fr
from flint import arb, fmpq, ctx

ctx.prec = 200
A = lambda q: arb(fmpq(Fr(q).numerator, Fr(q).denominator))
lam, al = A(Fr(37, 40)), A(Fr(3, 40))
pi = arb.pi()

# rigorous upper bounds as rationals (rounded UP)
U_ub = Fr('-1.3843989249340')        # certified value -1.38439892493401715... < this
A1e5_ub = Fr('1.3135191719802')      # certified 1.31351917198018 +/- 1.5e-15 < this
A200_ub = Fr('1.3455351940065')      # certified 1.34553519400640 +/- 6e-16 < this
assert A(U_ub) > arb('-1.38439892493401715650888837988593764 +/- 6e-36')
assert A(A1e5_ub) > arb('1.31351917198018 +/- 1.5e-15') and A(A200_ub) > arb('1.34553519400640 +/- 6e-16')

# kappa of (C.6) with a = 1/18, b = 1/3
a, b = A(Fr(1, 18)), A(Fr(1, 3))
F = 2 * pi * b.sqrt() + (1 + b).log() - 6 * al * a.log()
Y = (2 + a + b) / (b - a)
Rr = Y + (Y * Y - 1).sqrt()
kappa = F + lam * (1 - (al * lam).log()) + 2 * lam * Rr.log()
print(f"kappa(1/18,1/3) = {kappa.str(15)};  Y = {Y.str(10)} (paper 43/5);  kappa < 6933/500: {kappa < A(Fr(6933,500))}")
print(f"  paper (A.11): 2773171335773/200000000000 < kappa < 6932928339433/500000000000: "
      f"{A(Fr(2773171335773,200000000000)) < kappa < A(Fr(6932928339433,500000000000))}")
rho_ = Fr(6933, 500)

eps_lb = -(A1e5_ub + U_ub)       # eps can be anything < this
print(f"\n-(A_100000 + U) >= {float(eps_lb):.10f}")
thr = rho_ * Fr(37, 40) / eps_lb
print(f"mu threshold rho*lam/eps  = {float(thr):.6f}  (with rho = 6933/500);  with rho -> kappa: "
      f"{(kappa * lam / A(eps_lb)).str(10)}")
# 09's concrete choice
eps09 = Fr(7087, 100000); c09 = Fr(5221, 400)
print(f"09 choice: eps = 7087/100000 < -(A+U): {eps09 < eps_lb};  eps*c - lam = {eps09*c09 - Fr(37,40)} > 0: "
      f"{eps09*c09 > Fr(37,40)};  rho*c = {rho_*c09} = {float(rho_*c09):.6f} < 181: {rho_*c09 < 181}")
# Fauzan's original
epsF, cF = Fr(79, 1600), Fr(75, 4)
print(f"Fauzan: eps*c - lam = {epsF*cF - Fr(37,40)} (paper 1/24000), rho*c - 260 = {rho_*cF - 260} (paper -1/80);"
      f" his threshold rho*lam/eps = {float(rho_*Fr(37,40)/epsF):.4f}")
# Theorem 2.1-type decay at M = 200 with the new constants
print(f"\nirrationality decay: -1600(A_200 + U) >= {float(-1600*(A200_ub + U_ub)):.4f}  (Fauzan 139/5 = 27.8)")
# what integer mu is the best this route gives, and slack
for mu in (180, 181):
    # need c in (lam/eps, mu/rho)
    lo, hi = Fr(37, 40) / eps_lb, Fr(mu) / rho_
    print(f"mu = {mu}: admissible c interval ({float(lo):.6f}, {float(hi):.6f}) nonempty: {lo < hi}")
