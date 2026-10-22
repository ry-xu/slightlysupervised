"""
[NUMERICAL] sanity check, independent of the proofs: for q = 2, 3 there is no relation
    a*u*e + b*u + c*e + d = 0,   a,b,c,d in F_p[theta], deg <= dmax  (not all zero).
Unknowns: 4(dmax+1) coefficients.  Equations: coefficients of x^n, -dmax <= n <= L, of
a*(ue) + b*u + c*e + d  (a = sum a_i theta^i = sum a_i x^{-i}).
Full column rank  <=>  no nonzero relation of degree <= dmax.
We also report the smallest L for which the rank is full ("defect"), compared with the
generic value 3*dmax+3 (number of unknowns minus the dmax+1 coefficients of d which
just absorb the polar part):  L_min ~ 3 dmax + O(1) means (1,u,e,ue) behave generically.
"""
import sys, os, time
import numpy as np
import flint
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ff import *

def build(q, p, dmax, L, u, e, ue):
    rows = []
    for n in range(-dmax, L + 1):
        r = []
        for ser in (ue, u, e):
            for i in range(dmax + 1):
                m = n + i
                r.append(int(ser[m]) if 0 <= m < len(ser) else 0)
        for i in range(dmax + 1):          # d contributes x^{-i}: coefficient at n = -i
            r.append(1 if n == -i else 0)
        rows.append(r)
    return flint.nmod_mat(rows, p)

for q, p, dlist in [(2, 2, [20, 60, 150, 300]), (3, 3, [20, 60, 150])]:
    Nser = 8 * max(dlist) + 400
    u = u_series(q, Nser, p)
    e = e_series(q, Nser, p)
    ue = mul(u, e, p)
    for dmax in dlist:
        t0 = time.time()
        ncols = 4 * (dmax + 1)
        L = 4 * dmax + 60
        Mx = build(q, p, dmax, L, u, e, ue)
        rk = Mx.rank()
        # find minimal L with full rank (binary search on row-truncation)
        lo, hi = 0, L
        while lo < hi:
            mid = (lo + hi) // 2
            if build(q, p, dmax, mid, u, e, ue).rank() == ncols:
                hi = mid
            else:
                lo = mid + 1
        print(f'q={q} dmax={dmax}: rank {rk}/{ncols} (full={rk==ncols}); minimal L for full rank = {lo} '
              f'(generic ~ 3*dmax+3 = {3*dmax+3}); {time.time()-t0:.1f}s', flush=True)
