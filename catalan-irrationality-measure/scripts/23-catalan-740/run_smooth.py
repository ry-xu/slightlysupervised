"""Driver for smoothed optimisation.  [NUMERICAL]
usage: run_smooth.py Kp Kv mode out.npz [init.npz]
 mode 'litsey' : Litsey's bases only
 mode 'rich'   : Litsey's bases + log-spaced real bases near -1 and +1 for p and v
"""
import sys, json, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from trial import Trial
from optdual import Param
from optsmooth import Smooth
SP = '/private/tmp/claude-501/-Users-ryanxu-Documents-dev/313d84c6-3fcb-46b6-9f6f-56b316f7c448/scratchpad/'
d = json.load(open(SP + 'litsey/strongtrial.json'))
p = Trial.from_litsey(d['Pc'], d['P']); v = Trial.from_litsey(d['Vc'], d['V'])
Kp, Kv, mode, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
tp = [z for z, r in p.tails]; tv = [z for z, r in v.tails]
if mode == 'rich':
    deltas = np.geomspace(3e-7, 0.2, 16)
    for dlt in deltas:
        for sgn in (-1, 1):
            z = complex(round(sgn * (1 - dlt), 9))
            if all(abs(z - w) > 1e-12 for w in tp): tp.append(z)
            if all(abs(z - w) > 1e-12 for w in tv): tv.append(z)
par = Param(Kp, tp, Kv, tv)
th0 = par.from_trials(p, v)
if len(sys.argv) > 5 and sys.argv[5] != 'none':
    z = np.load(sys.argv[5]); th0 = z['theta']
S = Smooth(par, mx=30001, ms=20001)
h, mx, my = S.hard(th0)
print('start U =', -h, 'dim', par.dim, flush=True)
t0 = time.time()
taus = tuple(float(t) for t in (sys.argv[6].split(',') if len(sys.argv) > 6 else ['3e-6','1e-6','3e-7','1e-7','3e-8']))
th = S.run(th0, taus=taus, maxiter=4000)
h, mx, my = S.hard(th)
print(json.dumps(dict(U=-h, MX=mx, MY=my, norm=float(th @ S.H @ th), time=time.time() - t0)))
np.savez(out, theta=th, Kp=Kp, Kv=Kv, tp=np.array(tp), tv=np.array(tv))
