"""Exchange4 on a richer family, warm-started from a previous iterate.  [NUMERICAL]
usage: run_ex4_rich.py prev.npz Kp Kv iters out"""
import sys, json, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from optdual import Param
from optdual4 import Exchange4
prev = np.load(sys.argv[1]); Kp, Kv, iters, out = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
par0 = Param(int(prev['Kp']), [complex(z) for z in prev['tp']], int(prev['Kv']), [complex(z) for z in prev['tv']])
p0, v0 = par0.to_trials(prev['theta'])
extra = []
for w in (1.558, 1.5633):
    for rho in (0.99, 0.995, 0.998):
        extra.append(complex(round(rho * np.cos(w), 6), round(rho * np.sin(w), 6)))
for w in (0.0141,):
    for rho in (0.99, 0.997, 0.999):
        extra.append(complex(round(rho * np.cos(w), 6), round(rho * np.sin(w), 6)))
extra += [complex(-0.99999), complex(-0.999997)]
tp = list(par0.tp) + [z for z in extra if z not in par0.tp]
tv = list(par0.tv) + [z for z in extra if z not in par0.tv]
par = Param(Kp, tp, Kv, tv)
th0 = par.from_trials(p0, v0)
t0 = time.time()
E = Exchange4(par, mx=60001, ms=60001)
print('setup', round(time.time() - t0, 1), 'start U', -E.obj(th0)[0], 'dim', par.dim, flush=True)
th, U = E.run(th0, iters=iters, save=out + '.npy')
np.savez(out + '.npz', theta=th, Kp=Kp, Kv=Kv, tp=np.array(tp), tv=np.array(tv))
print(json.dumps(dict(U=U, time=time.time() - t0)))
