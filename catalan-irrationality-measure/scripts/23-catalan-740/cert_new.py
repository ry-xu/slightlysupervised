"""
Rigorous certificate for one field of an exported trial (Litsey format).    [rigorous: arb/acb]
usage: cert_new.py trial.json FIELD CAP          FIELD in {X, Y, M, NORM}
  X: sup over x in [-1,0)u(0,1) of the row field <= CAP
  Y: sup over s in (0,1) of the column field <= CAP
  M: sup over 401/2500 <= |x| <= 1 of the marked residual <= CAP
  NORM: prints enclosures of ||p||^2, ||v||^2, 2||p||^2+||v||^2/2 and p_1 (CAP ignored)
CAP is a decimal string, converted exactly to a rational.
"""
import sys, json, time
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from certify3 import *

d = json.load(open(sys.argv[1]))
field = sys.argv[2]
cap = float(Fr(sys.argv[3])) if len(sys.argv) > 3 else None
p = Tr(d['Pc'], d['P']); v = Tr(d['Vc'], d['V'])
t0 = time.time()
a0, a1, a2 = arb(0), arb(157079) / 100000, arb(157080) / 100000
if field == 'NORM':
    NP = p.norm2(); NV = v.norm2()
    K = 2 * NP + NV / 2
    print(json.dumps(dict(NP=str(NP), NV=str(NV), K=str(K), K_upper=upper(K), p1=str(p.first()),
                          p1_lower=float(p.first().lower()))))
    sys.exit(0)
res = []
if field == 'X':
    Fp = lambda t: X_point(p, v, t); Fd = lambda P: X_dball(p, v, P); Fs = lambda P: X_sing(p, v, P)
    res.append(('X [0, 1.57079]',) + certify_field(Fp, Fd, Fs, a0, a1, cap))
    res.append(('X [1.57079, 1.57080]',) + certify_field(Fp, Fd, Fs, a1, a2, cap))
    res.append(('X [1.57080, pi]',) + certify_field(Fp, Fd, Fs, a2, arb(PI.upper()), cap))
elif field == 'Y':
    Fp = lambda t: Y_point(v, t); Fd = lambda P: Y_dball(v, P); Fs = lambda P: Y_sing(v, P)
    res.append(('Y [0, 1.57079]',) + certify_field(Fp, Fd, Fs, a0, a1, cap))
    res.append(('Y [1.57079, pi/2+]',) + certify_field(Fp, Fd, Fs, a1, arb(HALFPI.upper()), cap))
elif field == 'M':
    tg = Q(Fr(401, 2500)).acos()
    Fp = lambda t: M_point(p, t); Fd = lambda P: M_dball(p, P); Fs = lambda P: M_sing(p, P)
    res.append(('M [0, tg]',) + certify_field(Fp, Fd, Fs, a0, arb(tg.upper()), cap))
    res.append(('M [pi-tg, pi]',) + certify_field(Fp, Fd, Fs, arb((PI - tg).lower()), arb(PI.upper()), cap))
ok = all(r[1] for r in res)
for r in res:
    print(r, flush=True)
print(json.dumps(dict(field=field, cap=sys.argv[3], ok=ok, worst=max(r[2] for r in res if r[2] is not None),
                      pieces=sum(r[3] for r in res), seconds=round(time.time() - t0, 1))))
