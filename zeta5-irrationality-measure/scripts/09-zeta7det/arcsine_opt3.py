"""Level-3 refinement (J=121) of the nested arcsine measure, warm-started from level 2."""
import sys, json
import numpy as np
sys.path.insert(0, './scripts/09-zeta7det')
from arcsine_opt import refine, optimize, lam
d = json.load(open('./data/09-zeta7det/arcsine_opt.json'))
a2, b2, c2 = np.array(d['2']['a']), np.array(d['2']['b']), np.array(d['2']['c'])
a, b = refine(a2, b2, 1)
cinit = np.zeros(len(a)); cinit[0::2] = c2 / 2; cinit[1::2] = (c2[:-1] + c2[1:]) / 4
cinit = np.maximum(cinit, 1e-6); cinit *= lam / cinit.sum()
best = optimize(a, b, cinit, iters=40, trust=0.002, verbose=False)
print("level 3: J=%d U_bound=%.7f sup=%.7f I=%.7f min c=%.2e" % (len(a), best[0], best[2], best[3], best[1].min()))
d['3p'] = dict(a=a.tolist(), b=b.tolist(), c=best[1].tolist(), U_bound=best[0], M0=best[2], I=best[3])
json.dump(d, open('./data/09-zeta7det/arcsine_opt.json', 'w'))
