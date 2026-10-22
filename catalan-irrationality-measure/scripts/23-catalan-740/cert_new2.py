"""
Rigorous certificate for one field (X or Y) with the second-order test (certify4) plus the
MVT/SING tests of certify3; anchored complex powers.            [rigorous: arb/acb, 256-bit]
usage: cert_new2.py trial.json FIELD CAP [lo hi]   (optional angle sub-range for parallel runs)
"""
import sys, json, time
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from certify3 import *
import certify4 as c4

d = json.load(open(sys.argv[1]))
field = sys.argv[2]
cap = float(Fr(sys.argv[3]))
p = Tr(d['Pc'], d['P']); v = Tr(d['Vc'], d['V'])
c4.add_L2(p); c4.add_L2(v)
t0 = time.time()
a0, a1, a2 = arb(0), arb(157079) / 100000, arb(157080) / 100000
res = []
if field == 'X':
    Fpd = lambda t: c4.X_point_d1(p, v, t); F2 = lambda P: c4.X_d2ball(p, v, P)
    Fd = lambda P: X_dball(p, v, P); Fs = lambda P: X_sing(p, v, P)
    rng = [(a0, a1), (a1, a2), (a2, arb(PI.upper()))]
elif field == 'Y':
    Fpd = lambda t: c4.Y_point_d1(v, t); F2 = lambda P: c4.Y_d2ball(v, P)
    Fd = lambda P: Y_dball(v, P); Fs = lambda P: Y_sing(v, P)
    rng = [(a0, a1), (a1, arb(HALFPI.upper()))]
if len(sys.argv) > 5 and sys.argv[4] == 'upper':        # Y on [lo, pi/2+]
    lo = arb(Fr(sys.argv[5]).numerator) / Fr(sys.argv[5]).denominator
    rng = [(lo, a1), (a1, arb(HALFPI.upper()))]
elif len(sys.argv) > 5:
    lo = arb(Fr(sys.argv[4]).numerator) / Fr(sys.argv[4]).denominator
    hi = arb(Fr(sys.argv[5]).numerator) / Fr(sys.argv[5]).denominator
    rng = [(lo, hi)]
for (a, b) in rng:
    r = c4.certify_field2(Fpd, F2, Fd, Fs, a, b, cap)
    res.append((f'{field} [{float(a.mid()):.6f}, {float(b.mid()):.6f}]',) + r)
    print(res[-1], flush=True)
ok = all(r[1] for r in res)
print(json.dumps(dict(field=field, cap=sys.argv[3], range=sys.argv[4:6], ok=ok,
                      worst=max(r[2] for r in res if r[2] is not None) if any(r[2] is not None for r in res) else None,
                      pieces=sum(r[3] for r in res), seconds=round(time.time() - t0, 1))))
