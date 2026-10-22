"""Final arithmetic for the irrationality-measure bound of zeta(5) via Fauzan's App. C with sharpened constants.
Inputs (certified elsewhere):  U <= U_cert (arcsine_certify.py level 2),  A_100000 <= A_M (exact_A.py, X=4000).
Checks Fauzan (5.14) -1/2 <= Q <= 13/8 for Q = R - x F on exact sample points of [20, 4000]."""
import sys
from fractions import Fraction as Fr
sys.path.insert(0, './scripts/09-zeta7det')
from exact_A import R, lam, al, fl, breakpoints
import random

def F(x):
    return 4 * lam + 2 * lam * (x - fl(x)) - 12 * lam * (al * x - fl(al * x))

pts = breakpoints(20, 4000)
random.seed(3)
qmin, qmax = Fr(10), Fr(-10)
for l, r in zip(pts, pts[1:]):
    for t in (Fr(1, 3), Fr(2, 3)):
        x = l + t * (r - l)
        Q = R(x) - x * F(x)
        qmin = min(qmin, Q); qmax = max(qmax, Q)
print("Q range on exact interior points of [20,4000]:", float(qmin), float(qmax), " (Fauzan (5.14): [-1/2, 13/8])")
assert qmin >= Fr(-1, 2) and qmax <= Fr(13, 8)

U_cert = Fr(-13843989092, 10 ** 10)    # upper bound from arcsine_certify.py level 3p (J=121), rounded up
A_M = Fr(13135191720, 10 ** 10)        # upper bound from exact_A.py (rounded up)
eps = -(A_M + U_cert)
print("eps >=", float(eps))
lamF = Fr(37, 40); rho = Fr(6933, 500)  # Fauzan (C.1)
eps_t = Fr(7087, 100000); c = Fr(5221, 400); mu = 181
assert eps_t < eps
print("eps_t*c - lam =", eps_t * c - lamF, ">0:", eps_t * c > lamF)
print("rho*c =", rho * c, float(rho * c), "< mu =", mu, ":", rho * c < mu)
print("best possible with these constants: mu > rho*lam/eps =", float(rho * lamF / eps))
