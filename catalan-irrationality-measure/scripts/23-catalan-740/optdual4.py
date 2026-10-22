"""
Exchange method with an exact primal active-set solver for the dual subproblem.     [NUMERICAL]

Subproblem (candidate sets I, J from the fine grids):
   min_theta theta^T H theta + beta|theta-theta_t|^2 + max_{i in I}(a_i.theta + w_i) + max_{j in J}(b_j.theta + w_j)
Dual: min_lam f(lam) = 1/4 lam^T Q lam - c.lam over lam >= 0, sum_I lam = 1, sum_J lam = 1,
   Q = G Hb^{-1} G^T, c = w + G c0, c0 = beta Hb^{-1} theta_t;  theta = c0 - 1/2 Hb^{-1} G^T lam.
Active-set (Wolfe-type) solver: support S grows by the most negative reduced gradient, equality-
constrained solves on S, ratio test removes components hitting zero.
"""
import sys, time
import numpy as np
from scipy.linalg import cho_factor, cho_solve
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from trial import LOG2
from optdual import Param, local_maxima
from optsmooth import fine_xgrid, fine_sgrid


def simplex_qp_active(G, c, grp, HbF, tol=1e-12, maxit=5000, verbose=False):
    """G: (n,d) rows, c: (n,), grp: (n,) group id 0/1, HbF: cho_factor of Hb.
    returns lam (n,), u = Hb^{-1} G^T lam."""
    n, d = G.shape
    S = []
    for g in (0, 1):
        idx = np.where(grp == g)[0]
        S.append(int(idx[np.argmax(c[idx])]))
    lam = np.zeros(n); lam[S] = 1.0
    HGS = {}                       # cache Hb^{-1} g_i
    def hg(i):
        if i not in HGS:
            HGS[i] = cho_solve(HbF, G[i])
        return HGS[i]
    for it in range(maxit):
        Sarr = np.array(S)
        HG = np.array([hg(i) for i in S]).T          # (d, s)
        u = HG @ lam[Sarr]
        grad = 0.5 * (G @ u) - c                     # (n,)
        # multipliers y_g = mean of grad on S within group (equal at a stationary point)
        y = np.zeros(2)
        for g in (0, 1):
            m = grp[Sarr] == g
            y[g] = grad[Sarr[m]].min() if m.any() else 0.0
        red = grad - y[grp]
        red[Sarr] = 0.0
        i_star = int(np.argmin(red))
        if red[i_star] > -tol:
            break
        S.append(i_star)
        # inner loop: solve equality QP on S, ratio test
        while True:
            Sarr = np.array(S)
            s = len(S)
            HG = np.array([hg(i) for i in S]).T
            Qss = G[Sarr] @ HG
            E = np.zeros((s, 2)); E[np.arange(s), grp[Sarr]] = 1.0
            K = np.zeros((s + 2, s + 2))
            K[:s, :s] = 0.5 * Qss + 1e-14 * np.eye(s)
            K[:s, s:] = -E; K[s:, :s] = E.T
            rhs = np.concatenate([c[Sarr], [1.0, 1.0]])
            try:
                sol = np.linalg.solve(K, rhs)
            except np.linalg.LinAlgError:
                sol = np.linalg.lstsq(K, rhs, rcond=None)[0]
            new = sol[:s]
            if (new >= -1e-15).all():
                lam[:] = 0; lam[Sarr] = np.maximum(new, 0)
                break
            old = lam[Sarr]
            dirn = new - old
            neg = dirn < 0
            ratios = np.where(neg, old / np.where(neg, -dirn, 1), np.inf)
            k = int(np.argmin(ratios)); t = ratios[k]
            lam[Sarr] = old + t * dirn
            lam[Sarr[k]] = 0.0
            S = [i for i in S if lam[i] > 1e-16]
            # keep at least one per group
            for g in (0, 1):
                if not any(grp[i] == g for i in S):
                    idx = np.where(grp == g)[0]
                    j = int(idx[np.argmax(c[idx])]); S.append(j); lam[j] = 1.0
    Sarr = np.array(S)
    HG = np.array([hg(i) for i in S]).T
    u = HG @ lam[Sarr]
    if verbose:
        print(f"   active-set: |S|={len(S)} its={it}", flush=True)
    return lam, u


class Exchange4:
    def __init__(self, par, mx=40001, ms=40001, ridge=1e-12):
        self.par = par
        H = par.H()
        self.H = 0.5 * (H + H.T) + ridge * np.eye(par.dim)
        self.xg = fine_xgrid(mx); self.sg = fine_sgrid(ms)
        self.WX, self.AX = par.Xrows(self.xg)
        self.WY, self.AY = par.Yrows(self.sg)

    @staticmethod
    def cands(f, tol, maxc, nb=3):
        i = np.where((f[1:-1] >= f[:-2]) & (f[1:-1] >= f[2:]))[0] + 1
        i = np.concatenate([i, [0, len(f) - 1]])
        i = i[f[i] >= f.max() - tol]
        if len(i) > maxc // (2 * nb + 1):
            i = i[np.argsort(-f[i])[:maxc // (2 * nb + 1)]]
        idx = (i[:, None] + np.arange(-nb, nb + 1)[None, :]).ravel()
        idx = idx[(idx >= 0) & (idx < len(f))]
        return np.unique(idx)

    def obj(self, th):
        fX = self.WX + self.AX @ th; fY = self.WY + self.AY @ th
        return th @ self.H @ th + fX.max() + fY.max() - 0.7 * LOG2, fX, fY

    def run(self, th0, iters=30, beta=1e-8, tolc=(0.03, 0.03), maxc=10**9, verbose=True, save=None):
        th = th0.copy()
        val, fX, fY = self.obj(th)
        for it in range(iters):
            t0 = time.time()
            iX = np.where(fX >= fX.max() - tolc[0])[0]
            iY = np.where(fY >= fY.max() - tolc[1])[0]
            G = np.vstack([self.AX[iX], self.AY[iY]])
            w = np.concatenate([self.WX[iX], self.WY[iY]])
            grp = np.concatenate([np.zeros(len(iX), int), np.ones(len(iY), int)])
            Hb = self.H + beta * np.eye(self.par.dim)
            HbF = cho_factor(Hb)
            c0 = cho_solve(HbF, beta * th)
            c = w + G @ c0
            lam, u = simplex_qp_active(G, c, grp, HbF)
            cand = c0 - 0.5 * u
            best_t, best_v = 0.0, val
            for t in [1.0, 0.5, 0.25, 0.1, 0.03]:
                v2, _, _ = self.obj(th + t * (cand - th))
                if v2 < best_v - 1e-14:
                    best_t, best_v = t, v2
                    break
            if best_t == 0:
                beta *= 10
            elif best_t == 1.0:
                beta = max(beta / 10, 1e-12)
            th = th + best_t * (cand - th)
            val, fX, fY = self.obj(th)
            if save is not None:
                np.save(save, th)
            if verbose:
                print(f"it {it}: U={-val:.10f} MX={fX.max():.9f} MY={fY.max():.9f} norm={th@self.H@th:.9f} "
                      f"t={best_t} beta={beta:.1e} |I|={len(iX)} |J|={len(iY)} supp={int((lam>0).sum())} "
                      f"{time.time()-t0:.1f}s", flush=True)
        return th, -val
