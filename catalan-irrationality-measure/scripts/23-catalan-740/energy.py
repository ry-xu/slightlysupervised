"""
23-catalan-740: primal (piecewise-constant density) optimum of multi-population
logarithmic energy functionals arising from OpenAI-005-type real-place bounds.

All functionals are concave quadratics in the population measures, maximized over
nonnegative measures of prescribed masses.  The discretized optimum is a LOWER bound
for the true supremum (every admissible certificate is an upper bound).  [NUMERICAL]

Populations used:
  'x'   : nodes on (-1,1)  (cells split at 0)
  'xp'  : nodes on (0,1)
  's'   : column variables on (0,1)

Kernels (cell averages of):
  Klog(x,x')  = log|x-x'|
  KJ(x,x')    = log(1-x x')
  Kxs(x,s)    = log(1-2 x s + x^2)
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize


def phi(u):
    u = np.asarray(u, float)
    out = np.zeros_like(u)
    nz = u != 0
    out[nz] = -(u[nz] ** 2) / 2 * np.log(np.abs(u[nz])) + 0.75 * u[nz] ** 2
    return out


def loglog_cells(ea, eb):
    """K[c,c'] = average of log|x-y| over cell c (grid ea) x cell c' (grid eb)."""
    a = ea[:-1][:, None]; b = ea[1:][:, None]
    c = eb[:-1][None, :]; d = eb[1:][None, :]
    I = phi(b - d) - phi(b - c) - phi(a - d) + phi(a - c)
    return I / ((b - a) * (d - c))


def smooth_cells(ex, ey, f, m=10):
    t, w = leggauss(m)
    xa, xb = ex[:-1], ex[1:]
    ya, yb = ey[:-1], ey[1:]
    X = (xa[:, None] + xb[:, None]) / 2 + (xb - xa)[:, None] / 2 * t[None, :]
    Y = (ya[:, None] + yb[:, None]) / 2 + (yb - ya)[:, None] / 2 * t[None, :]
    F = f(X[:, :, None, None], Y[None, None, :, :])
    return np.einsum('i,j,aibj->ab', w, w, F) / 4.0


def cell_avg_1d(edges, f, m=16):
    t, w = leggauss(m)
    a, b = edges[:-1], edges[1:]
    X = (a[:, None] + b[:, None]) / 2 + (b - a)[:, None] / 2 * t[None, :]
    return (f(X) * w[None, :]).sum(1) / 2.0


def graded(a, b, M):
    th = np.linspace(0, np.pi, M + 1)
    u = (1 - np.cos(th)) / 2
    return a + (b - a) * u


def grids(Mx, Ms):
    exn = graded(-1, 0, Mx // 2)
    exp_ = graded(0, 1, Mx // 2)
    ex = np.concatenate([exn, exp_[1:]])
    es = graded(0, 1, Ms)
    return ex, exp_, es


def J_ker(x, y):
    return np.log(np.maximum(1 - x * y, 1e-300))


def XS_ker(x, s):
    return np.log(np.maximum(1 - 2 * x * s + x * x, 1e-300))


class QP:
    """maximize const + sum_a f_a.w_a + sum_{a,b} w_a^T Q[a][b] w_b  (Q symmetric blockwise:
    we store full blocks, objective uses 1/2 sum_{a,b} w_a^T S_ab w_b with S = Q + Q^T)."""

    def __init__(self, sizes, masses):
        self.sizes = sizes
        self.masses = masses
        self.off = np.concatenate([[0], np.cumsum(sizes)])
        n = self.off[-1]
        self.H = np.zeros((n, n))   # objective = const + f.z + z^T H z  (H symmetric)
        self.f = np.zeros(n)
        self.const = 0.0

    def blk(self, a):
        return slice(self.off[a], self.off[a + 1])

    def add_quad(self, a, b, K, coef):
        """add coef * w_a^T K w_b to the objective."""
        if a == b:
            Ks = 0.5 * (K + K.T)
            self.H[self.blk(a), self.blk(a)] += coef * Ks
        else:
            self.H[self.blk(a), self.blk(b)] += 0.5 * coef * K
            self.H[self.blk(b), self.blk(a)] += 0.5 * coef * K.T

    def add_lin(self, a, v, coef=1.0):
        self.f[self.blk(a)] += coef * v

    def value(self, z):
        return self.const + self.f @ z + z @ self.H @ z

    def grad(self, z):
        return self.f + 2 * self.H @ z

    def solve(self, iters=4000, step=0.5, polish=True, verbose=False, z0=None):
        na = len(self.sizes)
        if z0 is None:
            z = np.concatenate([np.full(s, m / s) for s, m in zip(self.sizes, self.masses)])
        else:
            z = z0.copy()
        for it in range(iters):
            g = self.grad(z)
            for a in range(na):
                sl = self.blk(a)
                if self.masses[a] == 0:
                    continue
                ga = g[sl]
                w = z[sl] * np.exp(step * (ga - ga.max()))
                z[sl] = w / w.sum() * self.masses[a]
            if verbose and it % 1000 == 0:
                print(it, self.value(z), flush=True)
        v_eg = self.value(z)
        if not polish:
            return v_eg, z
        n = self.off[-1]
        cons = []
        for a in range(na):
            sl = self.blk(a)
            vec = np.zeros(n); vec[sl] = 1.0
            cons.append({'type': 'eq', 'fun': (lambda zz, vec=vec, m=self.masses[a]: vec @ zz - m),
                         'jac': (lambda zz, vec=vec: vec)})
        obj = lambda zz: (-self.value(zz), -self.grad(zz))
        res = minimize(obj, z, jac=True, method='SLSQP', bounds=[(0, None)] * n, constraints=cons,
                       options={'maxiter': 1000, 'ftol': 1e-15})
        zz = np.maximum(res.x, 0)
        # renormalise masses exactly
        for a in range(na):
            sl = self.blk(a)
            if self.masses[a] > 0:
                zz[sl] *= self.masses[a] / zz[sl].sum()
        return max(v_eg, self.value(zz)), zz
