"""Scan the generalized Fauzan bound model: margin(k,r,alpha) = A_inf(k,r,alpha,lam=1-alpha) + U*(r,alpha,lam)
(U* = true equilibrium value, A_inf = Fauzan-type local bounds with M -> oo, exact tail)."""
import sys, json
sys.path.insert(0, './scripts/09-zeta7det')
from bookkeeping import A_inf
from equilibrium import U_fauzan

out = sys.argv[1] if len(sys.argv) > 1 else None
rows = []
for r in (6, 8, 10, 12):
    for al in (0.05, 0.075, 0.1, 0.125, 0.15):
        lam = 1 - al
        try:
            U, info = U_fauzan(r, al, lam, n=4096, guess=(1e-4, 0.8))
        except Exception as e:
            print("U fail", r, al, e); continue
        rec = dict(r=r, al=al, lam=lam, U=U, a=info['a'], b=info['b'], minw=info['minw'])
        for k in (5, 7):
            A, inf = A_inf(k, r, al, lam, X=400.0, nout=6000)
            rec[f'A{k}'] = A
            rec[f'margin{k}'] = A + U
            rec[f'Iout{k}'] = inf['Iout']
            rec[f'Iin{k}'] = inf['Iin'] + inf['tail']
        print(json.dumps(rec), flush=True)
        rows.append(rec)
        if out:
            with open(out, 'a') as fh:
                fh.write(json.dumps(rec) + "\n")
