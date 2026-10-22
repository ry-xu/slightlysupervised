"""
Proximal exchange method with an interior-point solver for the small dual QP.   [NUMERICAL]

Problem:  min_theta  F(theta) = theta^T H theta + max_x X(x;theta) + max_s Y(s;theta)
Iteration: candidate sets I (x-points), J (s-points) = local maxima near the top (accumulated);
  subproblem  min_theta theta^T H theta + beta|theta-theta_t|^2 + max_I(a_i theta+w_i) + max_J(b_j theta+w_j)
  dual: max_{lam in D_I x D_J} lam.(w + G c0) - 1/4 lam^T G Hb^{-1} G^T lam ,  c0 = Hb^{-1} beta theta_t
  theta = c0 - 1/2 Hb^{-1} G^T lam.
"""
import sys, time, json
import numpy as np
from scipy.linalg import cho_factor, cho_solve
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from trial import LOG2
from optdual import Param, local_maxima
from optsmooth import fine_xgrid, fine_sgrid


def simplex_qp_ip(Q, c, groups, iters=200, tol=1e-13):
    """maximize lam.c - 1/4 lam^T Q lam  s.t. lam >= 0, sum over each group = 1.
    Primal-dual interior point on  min 1/4 l^T Q l - c.l ."""
    n = len(c)
    G = len(groups)
    E = np.zeros((G, n))
    for gi, idx in enumerate(groups):
        E[gi, idx] = 1.0
    lam = np.zeros(n)
    for idx in groups:
        lam[idx] = 1.0 / len(idx)
    z = np.ones(n) * 1e-1          # dual slack for lam >= 0
    y = np.zeros(G)
    Qh = 0.5 * Q                    # gradient of 1/4 l^T Q l is 1/2 Q l
    for it in range(iters):
        rd = Qh @ lam - c - E.T @ y - z          # dual residual
        rp = E @ lam - 1                          # primal residual
        mu = lam @ z / n
        if max(np.abs(rd).max(), np.abs(rp).max()) < tol and mu < tol:
            break
        sigma = 0.1
        # Newton system:  [Qh + Z/L, -E^T; E, 0] [dl; dy] = [-rd - (z - sigma mu / lam); -rp] ... (standard)
        D = z / lam
        K = np.zeros((n + G, n + G))
        K[:n, :n] = Qh + np.diag(D)
        K[:n, n:] = -E.T
        K[n:, :n] = E
        rhs = np.concatenate([-rd - z + sigma * mu / lam, -rp])
        try:
            sol = np.linalg.solve(K, rhs)
        except np.linalg.LinAlgError:
            sol = np.linalg.lstsq(K, rhs, rcond=None)[0]
        dl = sol[:n]; dy = sol[n:]
        dz = -z + sigma * mu / lam - D * dl
        # step length
        a = 1.0
        neg = dl < 0
        if neg.any():
            a = min(a, 0.99 * np.min(-lam[neg] / dl[neg]))
        neg = dz < 0
        if neg.any():
            a = min(a, 0.99 * np.min(-z[neg] / dz[neg]))
        lam = lam + a * dl; y = y + a * dy; z = z + a * dz
    return lam


class Exchange:
    def __init__(self, par, mx=40001, ms=40001, ridge=1e-12):
        self.par = par
        H = par.H()
        self.H = 0.5 * (H + H.T) + ridge * np.eye(par.dim)
        self.xg = fine_xgrid(mx); self.sg = fine_sgrid(ms)
        self.WX, self.AX = par.Xrows(self.xg)
        self.WY, self.AY = par.Yrows(self.sg)

    def obj(self, th):
        fX = self.WX + self.AX @ th; fY = self.WY + self.AY @ th
        return th @ self.H @ th + fX.max() + fY.max() - 0.7 * LOG2, fX, fY

    def run(self, th0, iters=40, beta=1e-3, tolact=1e-4, maxact=(800, 800), verbose=True, save=None):
        th = th0.copy()
        actX = np.empty(0); actY = np.empty(0); prevX = np.empty(0); prevY = np.empty(0)
        val, fX, fY = self.obj(th)
        for it in range(iters):
            nX = local_maxima(fX, self.xg, tolact, maxact[0])
            nY = local_maxima(fY, self.sg, tolact, maxact[1])
            actX = np.unique(np.concatenate([prevX, nX])) if it > 0 else nX
            actY = np.unique(np.concatenate([prevY, nY])) if it > 0 else nY
            prevX, prevY = nX, nY
            WXa, AXa = self.par.Xrows(actX)
            WYa, AYa = self.par.Yrows(actY)
            # prune points far below max
            keepX = (WXa + AXa @ th) >= fX.max() - 3 * tolact
            keepY = (WYa + AYa @ th) >= fY.max() - 3 * tolact
            actX, WXa, AXa = actX[keepX], WXa[keepX], AXa[keepX]
            actY, WYa, AYa = actY[keepY], WYa[keepY], AYa[keepY]
            Hb = self.H + beta * np.eye(self.par.dim)
            cf = cho_factor(Hb)
            Gm = np.vstack([AXa, AYa])
            HG = cho_solve(cf, Gm.T)                    # (d, n)
            c0 = cho_solve(cf, beta * th)
            Q = Gm @ HG
            c = np.concatenate([WXa, WYa]) + Gm @ c0
            nI = len(actX); nJ = len(actY)
            lam = simplex_qp_ip(Q, c, [np.arange(nI), nI + np.arange(nJ)])
            cand = c0 - 0.5 * HG @ lam
            # line search (convex along segment)
            best_t, best_v = 0.0, val
            for t in [1.0, 0.5, 0.25, 0.1]:
                v2, _, _ = self.obj(th + t * (cand - th))
                if v2 < best_v - 1e-14:
                    best_t, best_v = t, v2
                    break
            if best_t == 0:
                beta *= 4
            elif best_t == 1.0:
                beta = max(beta / 3, 1e-10)
            th = th + best_t * (cand - th)
            val, fX, fY = self.obj(th)
            if save is not None:
                np.save(save, th)
            if verbose:
                print(f"it {it}: U={-val:.10f} MX={fX.max():.9f} MY={fY.max():.9f} norm={th@self.H@th:.9f} "
                      f"t={best_t} beta={beta:.1e} |I|={nI} |J|={nJ}", flush=True)
        return th, -val
