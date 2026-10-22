"""16-mu5: final certified assembly of the irrationality-measure bound for zeta(5).

  mu(zeta(5)) <= rho_K * lambda / eps   for any rational rho_K with ||B_K|| <= e^{(rho_K + o(1)) K}
                                         and any eps < -(A_M + U)
  (Fauzan App. C.2, with Prop. C.1 replaced by Theorem 1 / Theorem 2 of notes/16-mu5.md).

Everything below is recomputed here, in exact rationals or arb/acb ball arithmetic:
  * the test measure rho (1600 nested arcsine components, exact rationals, data/16-mu5/rho1600.exact.json):
      - admissibility for Fauzan's Lemma 6.1 (positivity, mass, nesting, support, regularization constant <= 60),
      - sup_{t>=0}(2U^rho - V) < M0  and  U = lam*M0 - I(rho) + C*                      [Lemma 6.1 / Prop. 6.3]
      - inf_{supp rho}(2U^rho - V) > m0                                                    [Theorems 1, 2]
      - sup_{s in [al,1]}(2U^rho(-s^2) + phi(s)) < M1                                      [Theorem 1]
      - max_{tau in Gamma} (2U^rho(tau) - Vt(tau)) < M1c  on the polygonal contour Gamma  [Theorem 2]
  * A_M for M = 4*10^8 from Fauzan's (5.11), (5.17), (5.18), I_out = 127751/96000 (5.10), the exact
    int_20^4000 R/x^3 of 15-redteam-mu5 (data/15-redteam-mu5/int_20_4000.txt, = 09's value), and the tail bound
    int_4000^oo R/x^3 <= -lam/4000 + (3lam/(2al) + 13/16)/4000^2 (09/15; re-derived in notes/16-mu5.md);
  * the final inequalities eps*c > lam and rho_K*c < mu, in exact rationals.
"""
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import arb, fmpq, ctx
from certlib import Measure, A, LAM, AL, lamA, alA, CSTAR, PI
from certify import certify_sup, certify_inf, certify_neg
from contour_certify import CMeasure, certify_path

ctx.prec = 200
ROOT = '.'
P = dict(fn=f'{ROOT}/data/16-mu5/rho1600.exact.json', M0='-6.6607665', m0='-6.660814', M1='0.8691543',
         x1='-1.069911', Y='1/5', x2='1', M1c='0.808802', Mcut=4 * 10 ** 8, mu1='96.69', mu2='95.91')
for arg in sys.argv[1:]:
    k, v = arg.split('=')
    P[k] = v
M0, m0, M1, M1c = Fr(P['M0']), Fr(P['m0']), Fr(P['M1']), Fr(P['M1c'])
x1, Y, x2 = Fr(P['x1']), Fr(P['Y']), Fr(P['x2'])

t0 = time.time()
print(f"== 16-mu5 final assembly; measure {P['fn']}")
meas = Measure.from_exact_json(P['fn'])
ok = meas.check()
assert ok['positive'] and ok['mass_exact'] and ok['nested'] and ok['support_in_(0,2)'] and ok['reg_const'] <= 60

# ---------------- analytic constant (Lemma 6.1 / Prop. 6.3 with the new measure) ----------------
certify_sup(meas, M0)
I = meas.energy()
U = lamA * A(M0) - I + CSTAR
print(f"   I(rho) = {I.str(25)}\n   C*     = {CSTAR.str(25)}\n   U      = lam*M0 - I + C* = {U.str(25)}")

# ---------------- norm-rate constants ----------------
certify_inf(meas, m0)
certify_neg(meas, M1)
certify_path(CMeasure(meas), x1, Y, x2, M1c)
rho1, rho2 = M1 - m0, M1c - m0
print(f"   Theorem 1 (residue moduli):  rho_K = M1  - m0 = {rho1} = {float(rho1):.10f}")
print(f"   Theorem 2 (contour):         rho_K = M1c - m0 = {rho2} = {float(rho2):.10f}   (Fauzan: 6933/500 = 13.866)")

# ---------------- arithmetic constant A_M ----------------
I_out = Fr(127751, 96000)                                                   # (5.10)
I_3_20 = Fr(322437603634266857629, 7535670527041937280000)                    # (5.18)
p, q = open(f'{ROOT}/data/15-redteam-mu5/int_20_4000.txt').read().split()
I_20_4000 = Fr(int(p), int(q))                                                # exact (15-redteam-mu5, exact_R.py)
X = 4000
tail = -LAM / X + (3 * LAM / (2 * AL) + Fr(13, 16)) / X ** 2                  # upper bound for int_4000^oo R/x^3
Astar = I_out + I_3_20 + I_20_4000 + tail
Mcut = int(P['Mcut'])
assert Mcut % 40 == 0
A_M = Astar + 7 * LAM / Mcut - (Fr(2923, 240) - Fr(1, 4)) / Mcut ** 2 + Fr(32) / Mcut ** 3   # (5.11) + (5.17)
A1e5 = Astar + 7 * LAM / 100000 - (Fr(2923, 240) - Fr(1, 4)) / 100000 ** 2 + Fr(32) / 100000 ** 3
print(f"   A* (exact rational) = {float(Astar):.15f};  A_M (M = {Mcut}) = {float(A_M):.15f}")
print(f"   [check] int_20^4000 = {float(I_20_4000):.13f} (09/15: -0.0598432940999);  A_100000 = {float(A1e5):.13f} (15: 1.3135191719802)")

# ---------------- final inequalities ----------------
eps_max = -(A(A_M) + U)
print(f"   eps_max = -(A_M + U) = {eps_max.str(20)}")
eps = Fr(int(float(eps_max.lower().mid()) * 10 ** 12) - 1, 10 ** 12)
assert A(eps) < eps_max
for name, rr, mu in (('Theorem 1', rho1, Fr(P['mu1'])), ('Theorem 2', rho2, Fr(P['mu2']))):
    thr = A(rr) * lamA / eps_max
    c_lo, c_hi = LAM / eps, mu / rr
    assert c_lo < c_hi, (name, float(c_lo), float(c_hi))
    c = Fr(round((c_lo + c_hi) / 2 * 10 ** 9), 10 ** 9)
    assert eps * c > LAM and rr * c < mu
    print(f"   [{name}] threshold rho_K*lam/eps_max = {thr.str(12)};  eps = {eps}, c = {c}:"
          f" eps*c - lam = {float(eps*c - LAM):.3e} > 0, mu - rho_K*c = {float(mu - rr*c):.3e} > 0  ==> mu(zeta(5)) <= {float(mu)}")
print(f"==> mu(zeta(5)) <= {P['mu2']}  [{time.time()-t0:.0f}s]")
for label, rr, ee in [('Fauzan (rho=6933/500, eps=79/1600)', Fr(6933, 500), Fr(79, 1600)),
                      ('09 (rho=6933/500, eps=0.0708797530)', Fr(6933, 500), Fr('0.0708797530')),
                      ('Thm 1 rho, 09 eps', rho1, Fr('0.0708797530')),
                      ('Thm 1 rho, new eps', rho1, eps),
                      ('Thm 2 rho, 09 eps', rho2, Fr('0.0708797530')),
                      ('Thm 2 rho, new eps', rho2, eps)]:
    print(f"   {label:40s}: rho*lam/eps = {float(rr * LAM / ee):.4f}")
