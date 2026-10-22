"""Driver: start from Litsey's second trial, enlarge the family, optimise.  [NUMERICAL]
usage: run_opt2.py Kp Kv iters extra_p extra_v out.npz
 extra_*: comma-separated extra real tail bases (e.g. -0.99999,0.999999) or 'none'"""
import sys, json, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from trial import Trial
from optdual import Param
from optdual2 import optimize
SP = '/private/tmp/claude-501/-Users-ryanxu-Documents-dev/313d84c6-3fcb-46b6-9f6f-56b316f7c448/scratchpad/'
d = json.load(open(SP + 'litsey/strongtrial.json'))
p = Trial.from_litsey(d['Pc'], d['P']); v = Trial.from_litsey(d['Vc'], d['V'])
Kp, Kv, iters = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
ep = [] if sys.argv[4] == 'none' else [complex(float(t)) for t in sys.argv[4].split(',')]
ev = [] if sys.argv[5] == 'none' else [complex(float(t)) for t in sys.argv[5].split(',')]
tp = [z for z, r in p.tails] + [z for z in ep if z not in [w for w, _ in p.tails]]
tv = [z for z, r in v.tails] + [z for z in ev if z not in [w for w, _ in v.tails]]
par = Param(Kp, tp, Kv, tv)
th0 = par.from_trials(p, v)
t0 = time.time()
theta, info = optimize(par, th0, iters=iters, p1min=0.181, Rcap=-0.0524)
info['time'] = time.time() - t0
print(json.dumps(info))
np.savez(sys.argv[6], theta=theta, Kp=Kp, Kv=Kv, tp=np.array(tp), tv=np.array(tv))
