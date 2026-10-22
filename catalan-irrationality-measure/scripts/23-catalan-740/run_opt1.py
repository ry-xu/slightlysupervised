"""Re-optimize Litsey's second-trial family (same bases), exchange method.  [NUMERICAL]"""
import sys, json, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from trial import Trial
from optdual import Param, solve
from fractions import Fraction as Fr
d = json.load(open(sys.argv[1]))
p = Trial.from_litsey(d['Pc'], d['P']); v = Trial.from_litsey(d['Vc'], d['V'])
Kp = int(sys.argv[2]) if len(sys.argv) > 2 else 160
Kv = int(sys.argv[3]) if len(sys.argv) > 3 else 160
par = Param(Kp, [z for z, r in p.tails], Kv, [z for z, r in v.tails])
th0 = par.from_trials(p, v)
t0 = time.time()
theta, info = solve(par, th0, iters=int(sys.argv[4]) if len(sys.argv) > 4 else 8, p1min=0.181, Rcap=-0.0524)
info.pop('hist')
info['U'] = -info['dual']; info['time'] = time.time() - t0
print(json.dumps(info))
np.save(sys.argv[5] if len(sys.argv) > 5 else '/tmp/theta.npy', theta)
