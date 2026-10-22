"""Exchange optimisation, family: Kp modes, Kv modes, Litsey's tail bases (+ optional extra).  [NUMERICAL]
usage: run_ex3.py Kp Kv iters out.npz [init.npz|litsey]"""
import sys, json, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from trial import Trial
from optdual import Param
from optdual3 import Exchange
SP = '/private/tmp/claude-501/-Users-ryanxu-Documents-dev/313d84c6-3fcb-46b6-9f6f-56b316f7c448/scratchpad/'
d = json.load(open(SP + 'litsey/strongtrial.json'))
pL = Trial.from_litsey(d['Pc'], d['P']); vL = Trial.from_litsey(d['Vc'], d['V'])
Kp, Kv, iters, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
init = sys.argv[5] if len(sys.argv) > 5 else 'litsey'
tp = [z for z, r in pL.tails]; tv = [z for z, r in vL.tails]
par = Param(Kp, tp, Kv, tv)
if init == 'litsey':
    th0 = par.from_trials(pL, vL)
else:
    th0 = np.load(init)['theta'] if init.endswith('.npz') else np.load(init)
t0 = time.time()
E = Exchange(par, mx=40001, ms=40001)
print('setup', time.time() - t0, 'start U', -E.obj(th0)[0], flush=True)
th, U = E.run(th0, iters=iters, save=out + '.cur.npy')
np.savez(out, theta=th, Kp=Kp, Kv=Kv, tp=np.array(tp), tv=np.array(tv))
print(json.dumps(dict(U=U, time=time.time() - t0)))
