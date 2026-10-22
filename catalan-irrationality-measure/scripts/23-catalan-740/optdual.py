"""
Exchange-method optimizer for the kappa=2 energy dual of the (20,3,1,2,1) design in Litsey's
trial format (finite modes + geometric tails).                                   [NUMERICAL]

minimize  2||p||^2 + 1/2||v||^2 + max_x X(x) + max_s Y(s)
subject to (optional)  p_1 >= p1min,  marked residual M(x) <= Rcap on |x| >= 401/2500.

Parametrisation: theta = [p-modes (Kp), p-tails (real bases: 1 param; complex: 2),
                          v-modes (Kv), v-tails].
The objective is convex in theta (PSD quadratic + max of affine).
"""
import sys, json, time
import numpy as np
from numpy.polynomial import chebyshev as Ch
from scipy.optimize import minimize, LinearConstraint
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from trial import Trial, Xfield, Yfield, Mfield, xgrid, sgrid, LOG2

GAPX = 401 / 2500


class Param:
    def __init__(self, Kp, tails_p, Kv, tails_v):
        self.Kp, self.Kv = Kp, Kv
        self.tp = [complex(z) for z in tails_p]
        self.tv = [complex(z) for z in tails_v]
        self.np_t = sum(1 if z.imag == 0 else 2 for z in self.tp)
        self.nv_t = sum(1 if z.imag == 0 else 2 for z in self.tv)
        self.dim_p = Kp + self.np_t
        self.dim = self.dim_p + Kv + self.nv_t

    # ---- sequences of basis elements (for Gram) ----
    @staticmethod
    def _tail_cols(zs):
        cols = []
        for z in zs:
            if z.imag == 0:
                cols.append((z, 1.0 + 0j))
            else:
                cols.append((z, 1.0 + 0j))
                cols.append((z, 1j))
        return cols  # each col: sequence Re(r z^k)

    def gram_block(self, K, zs):
        cols = self._tail_cols(zs)
        n = K + len(cols)
        G = np.zeros((n, n))
        k = np.arange(1, K + 1)
        G[:K, :K] = np.diag(1.0 / k)
        for a, (z, r) in enumerate(cols):
            seq = (r * z ** k).real
            G[:K, K + a] = seq / k
            G[K + a, :K] = seq / k
            for b, (w, s) in enumerate(cols):
                # sum_k Re(r z^k) Re(s w^k)/k = 1/2 Re[ r s * sum (zw)^k/k + r conj(s) sum (z conj w)^k/k ]
                val = 0.5 * (-(r * s) * np.log(1 - z * w) - (r * np.conj(s)) * np.log(1 - z * np.conj(w))).real
                G[K + a, K + b] = val
        return G

    def H(self):
        Gp = self.gram_block(self.Kp, self.tp)
        Gv = self.gram_block(self.Kv, self.tv)
        H = np.zeros((self.dim, self.dim))
        H[:self.dim_p, :self.dim_p] = 2 * Gp
        H[self.dim_p:, self.dim_p:] = 0.5 * Gv
        return H

    # ---- feature rows ----
    def Tfeat(self, K, zs, x):
        """rows: d T(u,x)/d theta_u  for x array."""
        x = np.asarray(x, float)
        k = np.arange(1, K + 1)
        Tk = Ch.chebvander(x, K)[:, 1:] / k
        cols = []
        for z in zs:
            L = np.log(1 - 2 * x * z + z * z)
            if z.imag == 0:
                cols.append(-0.5 * L.real)
            else:
                cols.append(-0.5 * L.real)          # r = 1
                cols.append(-0.5 * (1j * L).real)   # r = i
        return np.hstack([Tk, np.array(cols).T]) if cols else Tk

    def Sfeat(self, K, zs, x):
        x = np.asarray(x, float)
        k = np.arange(1, K + 1)
        Sk = np.vander(x, K + 1, increasing=True)[:, 1:] / k
        cols = []
        for z in zs:
            L = np.log(1 - x * z)
            if z.imag == 0:
                cols.append(-L.real)
            else:
                cols.append(-L.real)
                cols.append(-(1j * L).real)
        return np.hstack([Sk, np.array(cols).T]) if cols else Sk

    def Xrows(self, x):
        """X(x) = W(x) + A theta."""
        x = np.asarray(x, float)
        W = 0.4 * np.log(np.abs(x)) + 0.1 * np.log(1 - x) - 1.35 * np.log(1 + x * x)
        A = np.hstack([-4 * self.Tfeat(self.Kp, self.tp, x), -self.Sfeat(self.Kv, self.tv, x)])
        return W, A

    def Yrows(self, s):
        s = np.asarray(s, float)
        W = 0.15 * np.log(s) + 0.05 * np.log(1 - s)
        A = np.hstack([np.zeros((len(s), self.dim_p)), 2 * self.Tfeat(self.Kv, self.tv, s)])
        return W, A

    def Mrows(self, x):
        x = np.asarray(x, float)
        W = -0.9 * np.log(1 + x * x)
        A = np.hstack([-2 * self.Tfeat(self.Kp, self.tp, x), np.zeros((len(x), self.dim - self.dim_p))])
        return W, A

    def p1row(self):
        a = np.zeros(self.dim)
        a[0] = 1.0  # mode k=1 (Kp >= 1)
        off = self.Kp
        for z in self.tp:
            if z.imag == 0:
                a[off] = z.real; off += 1
            else:
                a[off] = z.real; a[off + 1] = (1j * z).real; off += 2
        return a

    def to_trials(self, theta):
        def unpack(th, K, zs):
            c = th[:K]
            tails = []
            off = K
            for z in zs:
                if z.imag == 0:
                    tails.append((z, th[off])); off += 1
                else:
                    tails.append((z, th[off] + 1j * th[off + 1])); off += 2
            return Trial(c, tails)
        return unpack(theta[:self.dim_p], self.Kp, self.tp), unpack(theta[self.dim_p:], self.Kv, self.tv)

    def from_trials(self, p, v):
        """project trials with matching bases into theta."""
        th = np.zeros(self.dim)
        th[:min(self.Kp, len(p.c))] = p.c[:self.Kp]
        off = self.Kp
        dz = {complex(z): r for z, r in p.tails}
        for z in self.tp:
            r = dz.get(z, 0)
            if z.imag == 0:
                th[off] = r.real; off += 1
            else:
                th[off] = r.real; th[off + 1] = r.imag; off += 2
        base = self.dim_p
        th[base:base + min(self.Kv, len(v.c))] = v.c[:self.Kv]
        off = base + self.Kv
        dz = {complex(z): r for z, r in v.tails}
        for z in self.tv:
            r = dz.get(z, 0)
            if z.imag == 0:
                th[off] = r.real; off += 1
            else:
                th[off] = r.real; th[off + 1] = r.imag; off += 2
        return th


def local_maxima(f, grid, tol, maxn):
    i = np.where((f[1:-1] >= f[:-2]) & (f[1:-1] >= f[2:]))[0] + 1
    cand = list(i)
    cand += [0, len(f) - 1]
    cand = np.array(sorted(set(cand)))
    fm = f.max()
    cand = cand[f[cand] >= fm - tol]
    if len(cand) > maxn:
        cand = cand[np.argsort(-f[cand])[:maxn]]
    return grid[cand]


def refine_max(rowfun, theta, pts, half):
    """golden-section-free local refinement on tiny windows (returns refined points)."""
    out = []
    for x0, h in zip(pts, half):
        xs = np.linspace(x0 - h, x0 + h, 41)
        W, A = rowfun(xs)
        f = W + A @ theta
        out.append(xs[np.argmax(f)])
    return np.array(out)


def solve(par, theta0, iters=12, p1min=None, Rcap=None, wM=0.0, verbose=True, xg=None, sg=None):
    H = par.H()
    # ridge for conditioning
    H = H + 1e-14 * np.eye(par.dim)
    xg = xgrid(100001) if xg is None else xg
    sg = sgrid(100001) if sg is None else sg
    xgM = xg[np.abs(xg) >= GAPX]
    WXg, AXg = par.Xrows(xg)
    WYg, AYg = par.Yrows(sg)
    WMg, AMg = par.Mrows(xgM)
    actX = np.empty(0); actY = np.empty(0); actM = np.empty(0)
    theta = theta0.copy()
    hist = []
    for it in range(iters):
        fX = WXg + AXg @ theta; fY = WYg + AYg @ theta; fM = WMg + AMg @ theta
        val = theta @ H @ theta + fX.max() + fY.max()
        hist.append(val)
        if verbose:
            print(f"it {it}: dual={val - 0.7*LOG2:.10f}  MX={fX.max():.9f} MY={fY.max():.9f} "
                  f"norm={theta@H@theta:.9f} R={fM.max():.6f}", flush=True)
        nX = local_maxima(fX, xg, 2e-3, 400)
        nY = local_maxima(fY, sg, 2e-3, 200)
        actX = np.unique(np.concatenate([actX, nX]))
        actY = np.unique(np.concatenate([actY, nY]))
        if Rcap is not None:
            nM = local_maxima(fM, xgM, 2e-3, 100)
            actM = np.unique(np.concatenate([actM, nM]))
        WX, AX = par.Xrows(actX)
        WY, AY = par.Yrows(actY)
        d = par.dim
        # variables z = [theta, MX, MY]
        nz = d + 2
        HH = np.zeros((nz, nz)); HH[:d, :d] = H
        cvec = np.zeros(nz); cvec[d] = 1; cvec[d + 1] = 1
        rowsA = []; ub = []; lb = []
        A1 = np.hstack([AX, -np.ones((len(actX), 1)), np.zeros((len(actX), 1))])
        rowsA.append(A1); ub.append(-WX); lb.append(-np.inf * np.ones(len(actX)))
        A2 = np.hstack([AY, np.zeros((len(actY), 1)), -np.ones((len(actY), 1))])
        rowsA.append(A2); ub.append(-WY); lb.append(-np.inf * np.ones(len(actY)))
        if p1min is not None:
            a = np.concatenate([par.p1row(), [0, 0]])
            rowsA.append(a[None, :]); lb.append(np.array([p1min])); ub.append(np.array([np.inf]))
        if Rcap is not None and len(actM):
            WM, AM = par.Mrows(actM)
            A3 = np.hstack([AM, np.zeros((len(actM), 2))])
            rowsA.append(A3); ub.append(Rcap - WM); lb.append(-np.inf * np.ones(len(actM)))
        Acon = np.vstack(rowsA); lbv = np.concatenate(lb); ubv = np.concatenate(ub)
        z0 = np.concatenate([theta, [fX.max(), fY.max()]])
        obj = lambda z: (z[:d] @ H @ z[:d] + z[d] + z[d + 1], np.concatenate([2 * H @ z[:d], [1, 1]]))
        hess = lambda z: 2 * HH
        res = minimize(obj, z0, jac=True, hess=hess, method='trust-constr',
                       constraints=[LinearConstraint(Acon, lbv, ubv)],
                       options=dict(maxiter=3000, gtol=1e-13, xtol=1e-14, verbose=0))
        cand = res.x[:d]
        # line search on the true (grid) objective, convex along the segment
        def trueobj(th):
            return th @ H @ th + (WXg + AXg @ th).max() + (WYg + AYg @ th).max() + \
                (0 if Rcap is None else 1e3 * max(0.0, (WMg + AMg @ th).max() - Rcap))
        best_t, best_v = 0.0, trueobj(theta)
        for tt in [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125]:
            vv = trueobj(theta + tt * (cand - theta))
            if vv < best_v:
                best_t, best_v = tt, vv
        if verbose:
            print(f"   step t={best_t}  qp_status={res.status}", flush=True)
        theta = theta + best_t * (cand - theta)
    fX = WXg + AXg @ theta; fY = WYg + AYg @ theta; fM = WMg + AMg @ theta
    val = theta @ H @ theta + fX.max() + fY.max() - 0.7 * LOG2
    return theta, dict(dual=val, MX=fX.max(), MY=fY.max(), norm=theta @ H @ theta, R=fM.max(),
                       p1=float(par.p1row() @ theta), hist=hist, nactX=len(actX), nactY=len(actY))
