"""
Smoothed (log-sum-exp) L-BFGS optimizer for the kappa=2 dual, Litsey trial format.  [NUMERICAL]
   f_tau(theta) = theta^T H theta + tau*LSE(X/tau) + tau*LSE(Y/tau)      (convex, smooth)
with tau-continuation.  Grids are fixed fine Chebyshev-type grids plus refinement near 0, +-1.
"""
import sys, json, time
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from trial import LOG2
from optdual import Param


def fine_xgrid(m):
    th = np.linspace(0, np.pi, m)
    x = -np.cos(th)
    extra = [np.linspace(-0.25, 0.25, m // 4),
             -1 + np.geomspace(1e-14, 1e-2, 3000), 1 - np.geomspace(1e-14, 1e-2, 3000)]
    x = np.concatenate([x] + extra)
    x = x[(x != 0) & (x < 1) & (x >= -1)]
    return np.unique(x)


def fine_sgrid(m):
    th = np.linspace(0, np.pi, m)
    s = (1 - np.cos(th)) / 2
    extra = [np.geomspace(1e-14, 1e-2, 3000), 1 - np.geomspace(1e-15, 1e-2, 4000)]
    s = np.concatenate([s] + extra)
    s = s[(s > 0) & (s < 1)]
    return np.unique(s)


class Smooth:
    def __init__(self, par, mx=60001, ms=40001, ridge=1e-13):
        self.par = par
        self.H = par.H(); self.H = 0.5 * (self.H + self.H.T) + ridge * np.eye(par.dim)
        self.xg = fine_xgrid(mx); self.sg = fine_sgrid(ms)
        self.WX, self.AX = par.Xrows(self.xg)
        self.WY, self.AY = par.Yrows(self.sg)

    def fields(self, th):
        return self.WX + self.AX @ th, self.WY + self.AY @ th

    def hard(self, th):
        fX, fY = self.fields(th)
        return th @ self.H @ th + fX.max() + fY.max() - 0.7 * LOG2, fX.max(), fY.max()

    def soft(self, th, tau):
        fX, fY = self.fields(th)
        lx = logsumexp(fX / tau); ly = logsumexp(fY / tau)
        wx = np.exp(fX / tau - lx); wy = np.exp(fY / tau - ly)
        val = th @ self.H @ th + tau * lx + tau * ly
        grad = 2 * self.H @ th + self.AX.T @ wx + self.AY.T @ wy
        return val, grad

    def run(self, th0, taus=(1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7), maxiter=3000, verbose=True):
        th = th0.copy()
        for tau in taus:
            res = minimize(lambda t: self.soft(t, tau), th, jac=True, method='L-BFGS-B',
                           options=dict(maxiter=maxiter, maxcor=50, ftol=1e-16, gtol=1e-12))
            th = res.x
            if verbose:
                h, mx, my = self.hard(th)
                print(f"tau={tau:.0e}: U_hard={-h:.10f} MX={mx:.9f} MY={my:.9f} nit={res.nit}", flush=True)
        return th
