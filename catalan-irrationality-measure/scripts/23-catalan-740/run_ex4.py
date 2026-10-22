"""Driver for Exchange4.  usage: run_ex4.py Kp Kv iters out [init.npy|litsey] [extra tail bases for v,p]  [NUMERICAL]"""
import sys, json, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from trial import Trial
from optdual import Param
from optdual4 import Exchange4
SP = '/private/tmp/claude-501/-Users-ryanxu-Documents-dev/313d84c6-3fcb-46b6-9f6f-56b316f7c448/scratchpad/'
d = json.load(open(SP + 'litsey/strongtrial.json'))
pL = Trial.from_litsey(d['Pc'], d['P']); vL = Trial.from_litsey(d['Vc'], d['V'])
Kp, Kv, iters, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
init = sys.argv[5] if len(sys.argv) > 5 else 'litsey'
tp = [z for z, r in pL.tails]; tv = [z for z, r in vL.tails]
par = Param(Kp, tp, Kv, tv)
th0 = par.from_trials(pL, vL) if init == 'litsey' else np.load(init)
t0 = time.time()
E = Exchange4(par, mx=40001, ms=40001)
print('setup', round(time.time() - t0, 1), 'start U', -E.obj(th0)[0], 'dim', par.dim, flush=True)
th, U = E.run(th0, iters=iters, save=out + '.npy')
np.savez(out + '.npz', theta=th, Kp=Kp, Kv=Kv, tp=np.array(tp), tv=np.array(tv))
print(json.dumps(dict(U=U, time=time.time() - t0)))
