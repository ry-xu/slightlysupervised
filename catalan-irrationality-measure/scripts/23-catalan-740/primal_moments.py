"""Primal kappa=2 equilibrium (piecewise-constant densities) and its dual-optimal moment sequences
p*_k = int T_k dmu,  v*_k = int x^k dmu - 2 int T_k(s) dnu   (exact cell integrals).   [NUMERICAL]
usage: primal_moments.py Mx Ms kmax out.npz"""
import sys, json, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from laplace_explore import build

def cheb_cell_avgs(edges, kmax):
    """average of T_k over each cell, k=1..kmax, via antiderivative (T_{k+1}/(k+1) - T_{k-1}/(k-1))/2."""
    th_a = np.arccos(np.clip(edges[:-1], -1, 1)); th_b = np.arccos(np.clip(edges[1:], -1, 1))
    L = edges[1:] - edges[:-1]
    out = np.zeros((kmax, len(L)))
    for k in range(1, kmax + 1):
        def F(th):
            if k == 1:
                return np.cos(2 * th) / 4
            return 0.5 * (np.cos((k + 1) * th) / (k + 1) - np.cos((k - 1) * th) / (k - 1))
        out[k - 1] = (F(th_b) - F(th_a)) / L
    return out

def pow_cell_avgs(edges, kmax):
    a, b = edges[:-1], edges[1:]
    L = b - a
    out = np.zeros((kmax, len(L)))
    pa = a.copy(); pb = b.copy()
    for k in range(1, kmax + 1):
        pa = pa * a; pb = pb * b       # a^{k+1}
        out[k - 1] = (pb - pa) / ((k + 1) * L)
    return out

Mx, Ms, kmax, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
P, (ex, exp_, es) = build(1.0, Mx, Ms)
val, z = P.solve(iters=3000, step=0.5, polish=True)
sig = z[P.blk(0)]; tau = z[P.blk(1)]; nu = z[P.blk(2)]
npos = len(exp_) - 1
mu = tau.copy(); mu[-npos:] += sig
A = cheb_cell_avgs(ex, kmax) @ mu
B = pow_cell_avgs(ex, kmax) @ mu
Cc = cheb_cell_avgs(es, kmax) @ nu
pstar = A; vstar = B - 2 * Cc
print(json.dumps(dict(Phi=val, U=-val)))
np.savez(out, pstar=pstar, vstar=vstar, ex=ex, es=es, mu=mu, nu=nu, Phi=val)
