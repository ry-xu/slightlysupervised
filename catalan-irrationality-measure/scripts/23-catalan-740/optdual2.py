"""
Proximal exchange method for the kappa=2 dual (Litsey trial format).   [NUMERICAL]

Each step solves, for the current candidate sets I (row maxima), J (column maxima), M (marked):
   min_theta  theta^T H theta + beta |theta - theta_t|^2 + max_I (a_i theta + w_i) + max_J (b_j theta + w_j)
   s.t. (optional) marked rows  m_l theta + w_l <= Rcap,   p_1 >= p1min (as a penalty-free linear row)
via its dual: a concave QP over two simplices (+ nonneg multipliers for the extra rows),
solved by SLSQP; theta recovered in closed form.  The true objective (fine grids) is then line-searched.
"""
import sys, time, json
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from trial import xgrid, sgrid, LOG2
from optdual import Param, local_maxima

GAPX = 401 / 2500


def solve_sub(H, beta, th_t, AI, wI, BJ, wJ, Ex=None, ex=None):
    """Ex theta <= ex extra constraints (rows).  Returns theta."""
    d = H.shape[0]
    Hb = H + beta * np.eye(d)
    L = np.linalg.cholesky(Hb)
    def Hinv(g):
        return np.linalg.solve(L.T, np.linalg.solve(L, g))
    nI, nJ = len(wI), len(wJ)
    nE = 0 if Ex is None else len(ex)
    G = np.vstack([AI, BJ] + ([Ex] if nE else []))      # rows: features
    W = np.concatenate([wI, wJ] + ([-ex] if nE else []))
    # theta(lam) = Hb^{-1}(beta th_t - 1/2 G^T lam)
    c0 = Hinv(beta * th_t)
    HG = np.array([Hinv(g) for g in G])                  # (n, d)
    Q = G @ HG.T                                         # (n, n)  G Hb^{-1} G^T
    r = G @ c0
    # dual objective: lam.(W + r) - 1/4 lam^T Q lam   (+ const)  maximize
    n = nI + nJ + nE
    def f(lam):
        val = lam @ (W + r) - 0.25 * lam @ Q @ lam
        grad = (W + r) - 0.5 * Q @ lam
        return -val, -grad
    lam0 = np.concatenate([np.full(nI, 1.0 / nI), np.full(nJ, 1.0 / nJ), np.zeros(nE)])
    cons = [{'type': 'eq', 'fun': lambda l: l[:nI].sum() - 1, 'jac': lambda l: np.concatenate([np.ones(nI), np.zeros(nJ + nE)])},
            {'type': 'eq', 'fun': lambda l: l[nI:nI + nJ].sum() - 1, 'jac': lambda l: np.concatenate([np.zeros(nI), np.ones(nJ), np.zeros(nE)])}]
    res = minimize(f, lam0, jac=True, method='SLSQP', bounds=[(0, None)] * n, constraints=cons,
                   options={'maxiter': 2000, 'ftol': 1e-16})
    lam = np.maximum(res.x, 0)
    theta = c0 - 0.5 * HG.T @ lam
    return theta, res


def optimize(par, theta0, iters=30, beta=1e-4, p1min=None, Rcap=None, verbose=True, ngrid=200001,
             tolact=5e-4, maxact=(600, 300, 200), ridge=1e-12):
    H = par.H()
    H = 0.5 * (H + H.T) + ridge * np.eye(par.dim)
    xg = xgrid(ngrid); sg = sgrid(ngrid)
    xgM = xg[np.abs(xg) >= GAPX]
    WXg, AXg = par.Xrows(xg)
    WYg, AYg = par.Yrows(sg)
    WMg, AMg = par.Mrows(xgM)
    p1 = par.p1row()

    def true_obj(th):
        fX = WXg + AXg @ th; fY = WYg + AYg @ th
        val = th @ H @ th + fX.max() + fY.max() - 0.7 * LOG2
        pen = 0.0
        if Rcap is not None:
            pen += 10 * max(0.0, (WMg + AMg @ th).max() - Rcap)
        if p1min is not None:
            pen += 10 * max(0.0, p1min - p1 @ th)
        return val, pen

    theta = theta0.copy()
    actX = np.empty(0); actY = np.empty(0); actM = np.empty(0)
    for it in range(iters):
        fX = WXg + AXg @ theta; fY = WYg + AYg @ theta; fM = WMg + AMg @ theta
        val, pen = true_obj(theta)
        if verbose:
            print(f"it {it}: U={-val:.10f} pen={pen:.2e} MX={fX.max():.9f} MY={fY.max():.9f} "
                  f"R={fM.max():.6f} p1={p1@theta:.4f} beta={beta:.1e} |I|={len(actX)} |J|={len(actY)}", flush=True)
        actX = np.unique(np.concatenate([actX, local_maxima(fX, xg, tolact, maxact[0])]))
        actY = np.unique(np.concatenate([actY, local_maxima(fY, sg, tolact, maxact[1])]))
        Ex = None; ex = None
        rows = []; rhs = []
        if Rcap is not None:
            actM = np.unique(np.concatenate([actM, local_maxima(fM, xgM, 1e-2, maxact[2])]))
            WM, AM = par.Mrows(actM)
            rows.append(AM); rhs.append(Rcap - WM)
        if p1min is not None:
            rows.append(-p1[None, :]); rhs.append(np.array([-p1min]))
        if rows:
            Ex = np.vstack(rows); ex = np.concatenate(rhs)
        WX, AX = par.Xrows(actX)
        WY, AY = par.Yrows(actY)
        cand, res = solve_sub(H, beta, theta, AX, WX, AY, WY, Ex, ex)
        best_t, best = 0.0, val + pen
        for tt in [1.0, 0.5, 0.25, 0.1, 0.03]:
            v2, p2 = true_obj(theta + tt * (cand - theta))
            if v2 + p2 < best - 1e-13:
                best_t, best = tt, v2 + p2
                break
        if best_t == 0.0:
            beta *= 4
        elif best_t == 1.0:
            beta = max(beta / 2, 1e-9)
        theta = theta + best_t * (cand - theta)
        # prune active sets to points still near the max
        fX = WXg + AXg @ theta; fY = WYg + AYg @ theta
        WX, AX = par.Xrows(actX); WY, AY = par.Yrows(actY)
        actX = actX[(WX + AX @ theta) >= fX.max() - 20 * tolact]
        actY = actY[(WY + AY @ theta) >= fY.max() - 20 * tolact]
    val, pen = true_obj(theta)
    fX = WXg + AXg @ theta; fY = WYg + AYg @ theta; fM = WMg + AMg @ theta
    return theta, dict(U=-val, pen=pen, MX=fX.max(), MY=fY.max(), norm=theta @ H @ theta, R=fM.max(),
                       p1=float(p1 @ theta))
