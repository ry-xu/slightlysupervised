"""Initial trial = primal moments (k<=K) + Litsey's tails (k>K behaviour).   [NUMERICAL]
usage: init_from_primal.py primal.npz Kp Kv out.npz"""
import sys, json
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from trial import Trial, dual_value, xgrid, sgrid
from optdual import Param
SP = '/private/tmp/claude-501/-Users-ryanxu-Documents-dev/313d84c6-3fcb-46b6-9f6f-56b316f7c448/scratchpad/'
z = np.load(sys.argv[1]); Kp, Kv = int(sys.argv[2]), int(sys.argv[3])
d = json.load(open(SP + 'litsey/strongtrial.json'))
pL = Trial.from_litsey(d['Pc'], d['P']); vL = Trial.from_litsey(d['Vc'], d['V'])
tp = [zz for zz, r in pL.tails]; tv = [zz for zz, r in vL.tails]
par = Param(Kp, tp, Kv, tv)
# tails-only trials
pt = Trial(np.zeros(0), pL.tails); vt = Trial(np.zeros(0), vL.tails)
cp = z['pstar'][:Kp] - pt.seq(Kp)
cv = z['vstar'][:Kv] - vt.seq(Kv)
p = Trial(cp, pL.tails); v = Trial(cv, vL.tails)
th = par.from_trials(p, v)
out = dual_value(p, v, xgrid(400001), sgrid(400001))
print(json.dumps({k: float(val) for k, val in out.items()}), 'U=', -out['dual'])
np.savez(sys.argv[4], theta=th, Kp=Kp, Kv=Kv, tp=np.array(tp), tv=np.array(tv))
