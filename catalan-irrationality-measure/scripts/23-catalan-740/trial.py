"""
Trial sequences in Litsey's format (CatalanMeasure/Cert/FieldEval.lean):
    u_k = c_k + Re sum_j r_j z_j^k      (k >= 1; c_k = 0 for k > K)
and the kappa=2 energy-dual fields of the (20,3,1,2,1) design (OpenAI 005 (78), shifted weights
R = <2/5, 1/10, 27/20, 0, 1/5>, C = <3/20, 1/20, 0>):

  X(x) = 2/5 log|x| + 1/10 log(1-x) - 27/20 log(1+x^2) - 4 T(p,x) - S(v,x)   (row field)
  Y(s) = 3/20 log s + 1/20 log(1-s) + 2 T(v,s)                               (column field)
  dual(p,v) = -(7/10) log 2 + 2||p||^2 + 1/2 ||v||^2 + sup X + sup Y   (>= -U_target required)
  marked residual  M(x) = -(9/10) log(1+x^2) - 2 T(p,x)     on 401/2500 <= |x| <= 1

with ||u||^2 = sum u_k^2/k,  T(u,x) = sum u_k T_k(x)/k,  S(u,x) = sum u_k x^k/k.
Float64 evaluation (exploration).  Certified evaluation is in certify.py.
"""
import json
import numpy as np
from numpy.polynomial import chebyshev as C
from fractions import Fraction as Fr

LOG2 = np.log(2.0)


class Trial:
    def __init__(self, c, tails):
        """c: array of c_1..c_K (floats); tails: list of (z complex, r complex)."""
        self.c = np.asarray(c, float)
        self.tails = [(complex(z), complex(r)) for z, r in tails]

    @staticmethod
    def from_litsey(cs, tl):
        c = np.array([v / 1e8 for v in cs])
        tails = []
        for zr, zi, rr, ri in tl:
            tails.append((complex(float(Fr(zr)), float(Fr(zi))), complex(float(Fr(rr)), float(Fr(ri)))))
        return Trial(c, tails)

    def seq(self, kmax):
        k = np.arange(1, kmax + 1)
        u = np.zeros(kmax)
        K = min(len(self.c), kmax)
        u[:K] += self.c[:K]
        for z, r in self.tails:
            u += (r * z ** k).real
        return u

    def norm2(self):
        K = len(self.c)
        k = np.arange(1, K + 1)
        t = np.zeros(K)
        for z, r in self.tails:
            t += (r * z ** k).real
        val = np.sum((self.c ** 2 + 2 * self.c * t) / k)
        # tail-tail: sum_k t_k^2/k = -1/2 Re sum_{j,l} [r_j r_l log(1-z_j z_l) + r_j conj(r_l) log(1-z_j conj(z_l))]
        s = 0.0
        for zj, rj in self.tails:
            for zl, rl in self.tails:
                s += (rj * rl * np.log(1 - zj * zl) + rj * np.conj(rl) * np.log(1 - zj * np.conj(zl))).real
        return val - 0.5 * s

    def T(self, x):
        x = np.asarray(x, float)
        K = len(self.c)
        coef = np.concatenate([[0.0], self.c / np.arange(1, K + 1)])
        val = C.chebval(x, coef)
        for z, r in self.tails:
            val = val - 0.5 * (r * np.log(1 - 2 * x * z + z * z)).real
        return val

    def S(self, x):
        x = np.asarray(x, float)
        K = len(self.c)
        coef = np.concatenate([[0.0], self.c / np.arange(1, K + 1)])
        val = np.polynomial.polynomial.polyval(x, coef)
        for z, r in self.tails:
            val = val - (r * np.log(1 - x * z)).real
        return val


def Xfield(p, v, x):
    x = np.asarray(x, float)
    return 0.4 * np.log(np.abs(x)) + 0.1 * np.log(1 - x) - 1.35 * np.log(1 + x * x) - 4 * p.T(x) - v.S(x)


def Yfield(v, s):
    s = np.asarray(s, float)
    return 0.15 * np.log(s) + 0.05 * np.log(1 - s) + 2 * v.T(s)


def Mfield(p, x):
    x = np.asarray(x, float)
    return -0.9 * np.log(1 + x * x) - 2 * p.T(x)


def xgrid(m=200001):
    th = np.linspace(0, np.pi, m)
    x = -np.cos(th)
    x = x[(x != 0) & (x < 1)]
    # refine near 0 too
    xs = np.concatenate([x, np.linspace(-0.3, 0.3, 60001)])
    xs = xs[(xs != 0) & (xs < 1)]
    return np.unique(xs)


def sgrid(m=200001):
    th = np.linspace(0, np.pi, m)
    s = (1 - np.cos(th)) / 2
    return s[(s > 0) & (s < 1)]


def dual_value(p, v, xg=None, sg=None):
    xg = xgrid() if xg is None else xg
    sg = sgrid() if sg is None else sg
    MX = Xfield(p, v, xg).max()
    MY = Yfield(v, sg).max()
    K = 2 * p.norm2() + 0.5 * v.norm2()
    return dict(norm=K, MX=MX, MY=MY, dual=-0.7 * LOG2 + K + MX + MY, p1=p.seq(1)[0])


if __name__ == "__main__":
    import sys
    d = json.load(open(sys.argv[1]))
    p = Trial.from_litsey(d['Pc'], d['P'])
    v = Trial.from_litsey(d['Vc'], d['V'])
    out = dual_value(p, v)
    xg = xgrid()
    sel = np.abs(xg) >= 401 / 2500
    out['R'] = float(Mfield(p, xg[sel]).max())
    out['NP'] = p.norm2(); out['NV'] = v.norm2()
    out['U'] = -out['dual']
    print(json.dumps(out, indent=1))
