"""15-redteam-mu5: exact proof-check of Fauzan's identity (5.12)-(5.13), R(x) = x F(x) + Q(x), on [3,60],
then for all x >= 20 by 40-periodicity of both sides (Q of (5.13) depends only on {x},{al x},{2Hx},{2x},{2lam x};
R - xF is 40-periodic on [20,oo) by exact_R.py (3a),(3b)).  On each piece both sides are polynomials of
degree <= 2 in x, so agreement at 3 distinct points of a piece proves agreement on the piece.
Motivation: C. W. Lai (ai.vixra.org 2609.0080, erratum (8)) claims this identity 'fails by a term of order x'."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from exact_R import R, F, Q_fauzan, breakpoints
bps = breakpoints(3, 60)
bad = 0
for l, r in zip(bps[:-1], bps[1:]):
    for w in (Fr(1, 5), Fr(1, 2), Fr(4, 5)):
        x = l + (r - l) * w
        if R(x) != x * F(x) + Q_fauzan(x):
            bad += 1
print(f"pieces on [3,60]: {len(bps)-1}; identity R = xF + Q(5.13) violated at {bad} of {3*(len(bps)-1)} test points")
# far out: spot checks at large x (exact)
import random
random.seed(1)
far = [Fr(random.randint(10**8, 10**9), random.choice([7, 11, 13, 1000003])) for _ in range(300)]
print("spot checks at x ~ 1e5..1e9:", all(R(x) == x * F(x) + Q_fauzan(x) for x in far))
