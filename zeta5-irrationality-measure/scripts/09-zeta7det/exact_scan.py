"""Exact Fauzan-type determinants for families (K, N, r, h) = (Km n, Nm n, r, hm n): log P(zeta(k)) and content.
usage: python exact_scan.py ks Km Nm r hm nmax out.jsonl [B]"""
import sys, json, time
sys.path.insert(0, './scripts/09-zeta7det')
from zdet import run, fauzan_exps

ks = [int(x) for x in sys.argv[1].split(',')]
Km, Nm, r, hm, nmax = map(float, sys.argv[2:7])
r = int(r); nmax = int(nmax)
out = sys.argv[7]
for n in range(1, nmax + 1):
    K, N, h = int(round(Km * n)), int(round(Nm * n)), int(round(hm * n))
    for k in ks:
        t0 = time.time()
        res, D, P = run(k, fauzan_exps(K, N, r), h, B=4 * K)
        res.update(n=n, K=K, N=N, r=r, h=h, wall=time.time() - t0)
        res['vp'] = {str(p): v for p, v in res.get('vp', {}).items()}
        print(json.dumps({x: res[x] for x in res if x != 'vp'}), flush=True)
        with open(out, 'a') as fh:
            fh.write(json.dumps(res) + "\n")
