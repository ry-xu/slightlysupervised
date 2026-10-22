"""Exact test of the combined functional a*mu3(-tR)+b*mu5(R) (values in Q + Q(a zeta3 + b zeta5)).
usage: python run_comb35.py a b Km Nm r hm nmax out"""
import sys, json
sys.path.insert(0, './scripts/09-zeta7det')
from zgen import Comb35, run_gen, fauzan_like
a, b = int(sys.argv[1]), int(sys.argv[2])
Km, Nm, r, hm, nmax = float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7])
out = sys.argv[8]
fam = Comb35(a, b)
for n in range(1, nmax + 1):
    K, N, h = int(round(Km * n)), int(round(Nm * n)), int(round(hm * n))
    res, D, P = run_gen(fam, fauzan_like(fam, K, N, r), h, B=4 * K + 10)
    res.update(n=n, K=K, N=N, r=r, a=a, b=b, logP_K2=res['logP'] / K ** 2)
    print(json.dumps({x: res[x] for x in res if x != 'vp'}), flush=True)
    with open(out, 'a') as fh:
        fh.write(json.dumps(res) + "\n")
