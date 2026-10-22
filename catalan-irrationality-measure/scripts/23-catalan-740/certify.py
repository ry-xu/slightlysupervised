"""
Rigorous (arb ball arithmetic) certificates for kappa=2 trial data in Litsey's exact format
(CatalanMeasure/Cert/FieldEval.lean, Cert/ShiftedFields.lean, Cert/MarkedResidual.lean):

  u_k = cs_{k-1}/10^8 + Re sum_j r_j z_j^k,   cs in Z, z_j, r_j in Q + Q i, |z_j| < 1.

Certified quantities (upper bounds unless noted):
  NP >= ||p||^2, NV >= ||v||^2           (||u||^2 = sum_k u_k^2/k, closed form)
  MX >= sup_{x in [-1,0) u (0,1)} X(x),   X = 2/5 log|x| + 1/10 log(1-x) - 27/20 log(1+x^2) - 4T(p,x) - S(v,x)
  MY >= sup_{s in (0,1)} Y(s),            Y = 3/20 log s + 1/20 log(1-s) + 2T(v,s)
  RM >= sup_{401/2500 <= |x| <= 1} M(x),  M = -(9/10) log(1+x^2) - 2T(p,x)
  p_1 (exact rational).

Method on each closed piece [a,b] (exact rational endpoints):
  sup_{[a,b]} F <= F(m) + h * sup_{[a,b]} |F'|,   m = (a+b)/2, h = (b-a)/2,
  F(m) evaluated as an arb ball at the exact rational m, F' enclosed by evaluating its closed
  form on the ball [a,b] (Chebyshev U-recurrence + rational tail terms).  Pieces failing the
  target are bisected.  Pieces touching the singular endpoints (x=0, x=1, s=0, s=1) are handled
  by monotone majorants (see sing_* below).
"""
import sys, json
from fractions import Fraction as Fr
from flint import arb, acb, ctx

ctx.prec = 160


def Q(fr):
    fr = Fr(fr)
    return arb(fr.numerator) / fr.denominator


class ExactTrial:
    def __init__(self, cs, tails):
        """cs: list of ints (scaled 1e8); tails: list of [zr, zi, rr, ri] as Fraction-compatible strings."""
        self.cs = [int(c) for c in cs]
        self.K = len(self.cs)
        self.c = [arb(c) / 10**8 for c in self.cs]
        self.tails = []
        for zr, zi, rr, ri in tails:
            z = acb(Q(zr), Q(zi)); r = acb(Q(rr), Q(ri))
            assert Fr(zr) ** 2 + Fr(zi) ** 2 < 1
            self.tails.append((z, r))

    def first(self):
        """exact u_1 as Fraction."""
        return None

    def norm2(self):
        K = self.K
        s = arb(0)
        for k in range(1, K + 1):
            t = arb(0)
            for z, r in self.tails:
                t += (r * z ** k).real
            ck = self.c[k - 1]
            s += (ck * ck + 2 * ck * t) / k
        tt = acb(0)
        for zj, rj in self.tails:
            for zl, rl in self.tails:
                tt += rj * rl * (1 - zj * zl).log() + rj * rl.conjugate() * (1 - zj * zl.conjugate()).log()
        return s - tt.real / 2

    # ---- T(u,x) = sum c_k T_k(x)/k - 1/2 Re sum r log(1-2xz+z^2) ; derivative sum c_k U_{k-1}(x) + Re sum r z/(1-2xz+z^2)
    def T(self, x):
        K = self.K
        # Clenshaw for sum_{k=1}^K a_k T_k(x), a_k = c_k/k
        b1 = arb(0); b2 = arb(0)
        for k in range(K, 0, -1):
            b0 = self.c[k - 1] / k + 2 * x * b1 - b2
            b2 = b1; b1 = b0
        val = x * b1 - b2  # sum_{k>=1} a_k T_k = x*b1 - b2  (with a_0 = 0)
        for z, r in self.tails:
            val -= (r * (1 - 2 * x * z + z * z).log()).real / 2
        return val

    def dT(self, x):
        """sum c_k U_{k-1}(x) + Re sum r z/(1-2xz+z^2); x may be a wide ball."""
        K = self.K
        # Clenshaw for sum_{j=0}^{K-1} c_{j+1} U_j(x)
        b1 = arb(0); b2 = arb(0)
        for j in range(K - 1, -1, -1):
            b0 = self.c[j] + 2 * x * b1 - b2
            b2 = b1; b1 = b0
        val = b1  # sum_{j>=0} c_{j+1}U_j = b_0 (U-Clenshaw: S = b0)
        for z, r in self.tails:
            val += (r * z / (1 - 2 * x * z + z * z)).real
        return val

    def S(self, x):
        K = self.K
        acc = arb(0)
        for k in range(K, 0, -1):
            acc = acc * x + self.c[k - 1] / k
        val = acc * x
        for z, r in self.tails:
            val -= (r * (1 - x * z).log()).real
        return val

    def dS(self, x):
        K = self.K
        acc = arb(0)
        for k in range(K, 0, -1):
            acc = acc * x + self.c[k - 1]
        val = acc  # sum c_k x^{k-1}
        for z, r in self.tails:
            val += (r * z / (1 - x * z)).real
        return val


# ------------------------------------------------------------------ fields
def Xval(p, v, x):
    return Q(Fr(2, 5)) * abs(x).log() + Q(Fr(1, 10)) * (1 - x).log() - Q(Fr(27, 20)) * (1 + x * x).log() \
        - 4 * p.T(x) - v.S(x)


def dXval(p, v, x):
    return Q(Fr(2, 5)) / x - Q(Fr(1, 10)) / (1 - x) - Q(Fr(27, 20)) * 2 * x / (1 + x * x) - 4 * p.dT(x) - v.dS(x)


def Yval(v, s):
    return Q(Fr(3, 20)) * s.log() + Q(Fr(1, 20)) * (1 - s).log() + 2 * v.T(s)


def dYval(v, s):
    return Q(Fr(3, 20)) / s - Q(Fr(1, 20)) / (1 - s) + 2 * v.dT(s)


def Mval(p, x):
    return -Q(Fr(9, 10)) * (1 + x * x).log() - 2 * p.T(x)


def dMval(p, x):
    return -Q(Fr(9, 10)) * 2 * x / (1 + x * x) - 2 * p.dT(x)


def ball(a, b):
    """arb ball containing the rational interval [a,b]."""
    m = (Fr(a) + Fr(b)) / 2
    h = (Fr(b) - Fr(a)) / 2
    return Q(m) + Q(h * (1 + Fr(1, 10**12))) * arb(0, 1)   # slightly enlarged radius


def upper(x):
    return float(x.mid()) + float(x.rad())


def cert_sup(F, dF, a, b, target, maxdepth=60, stats=None):
    """prove sup_{[a,b]} F <= target by bisection. returns (ok, worst_bound, npieces)."""
    stack = [(Fr(a), Fr(b), 0)]
    worst = -1e300
    n = 0
    while stack:
        lo, hi, dep = stack.pop()
        m = (lo + hi) / 2
        h = (hi - lo) / 2
        fm = F(Q(m))
        d = dF(ball(lo, hi))
        bnd = fm + Q(h) * abs(d)
        ub = upper(bnd)
        if ub <= target and d.is_finite():
            worst = max(worst, ub); n += 1
            continue
        if dep >= maxdepth:
            return False, ub, n, (float(lo), float(hi))
        stack.append((lo, m, dep + 1)); stack.append((m, hi, dep + 1))
    return True, worst, n, None


def sing_piece(G, c, e, side, d0, target, ratio=Fr(1, 4), dmin=Fr(1, 10**40)):
    """Prove c*log(dist) + G(x) <= target for x between e (excluded) and e+side*d0, c>0.
    dist = |x-e|.  Pieces [d_{j+1}, d_j] with d_{j+1} = ratio*d_j; last piece (0, dmin] uses
    G on the closed ball and log(dmin)."""
    d = Fr(d0)
    worst = -1e300
    n = 0
    while True:
        lo = d * ratio if d > dmin else Fr(0)
        a = e + side * lo; b = e + side * d
        a, b = (a, b) if a <= b else (b, a)
        g = G(ball(a, b))
        bnd = c * Q(d).log() + g
        ub = upper(bnd)
        if not (ub <= target and g.is_finite()):
            return False, ub, n, (float(a), float(b))
        worst = max(worst, ub); n += 1
        if lo == 0:
            break
        d = lo
    return True, worst, n, None


def cert_trial(p, v, MX, MY, RM, d0=Fr(1, 1000), d1=Fr(1, 10**6), verbose=True):
    """Certify sup X <= MX, sup Y <= MY, sup_{|x|>=401/2500} M <= RM. Returns dict."""
    out = {}
    def Xs(x): return Xval(p, v, x)
    def dXs(x): return dXval(p, v, x)
    # singular parts of X: near 0 (c=2/5) and near 1 (c=1/10)
    G0 = lambda x: Q(Fr(1, 10)) * (1 - x).log() - Q(Fr(27, 20)) * (1 + x * x).log() - 4 * p.T(x) - v.S(x)
    G1 = lambda x: Q(Fr(2, 5)) * abs(x).log() - Q(Fr(27, 20)) * (1 + x * x).log() - 4 * p.T(x) - v.S(x)
    r = []
    r.append(('X[-1,-d0]',) + cert_sup(Xs, dXs, -1, -d0, MX))
    r.append(('X(-d0,0)',) + sing_piece(G0, Q(Fr(2, 5)), Fr(0), -1, d0, MX))
    r.append(('X(0,d0)',) + sing_piece(G0, Q(Fr(2, 5)), Fr(0), 1, d0, MX))
    r.append(('X[d0,1-d1]',) + cert_sup(Xs, dXs, d0, 1 - d1, MX))
    r.append(('X(1-d1,1)',) + sing_piece(G1, Q(Fr(1, 10)), Fr(1), -1, d1, MX))
    Ys = lambda s: Yval(v, s); dYs = lambda s: dYval(v, s)
    H0 = lambda s: Q(Fr(1, 20)) * (1 - s).log() + 2 * v.T(s)
    H1 = lambda s: Q(Fr(3, 20)) * s.log() + 2 * v.T(s)
    r.append(('Y(0,d0)',) + sing_piece(H0, Q(Fr(3, 20)), Fr(0), 1, d0, MY))
    r.append(('Y[d0,1-d1]',) + cert_sup(Ys, dYs, d0, 1 - d1, MY))
    r.append(('Y(1-d1,1)',) + sing_piece(H1, Q(Fr(1, 20)), Fr(1), -1, d1, MY))
    if RM is not None:
        Ms = lambda x: Mval(p, x); dMs = lambda x: dMval(p, x)
        g = Fr(401, 2500)
        r.append(('M[-1,-g]',) + cert_sup(Ms, dMs, -1, -g, RM))
        r.append(('M[g,1]',) + cert_sup(Ms, dMs, g, 1, RM))
    ok = all(t[1] for t in r)
    if verbose:
        for t in r:
            print(t, flush=True)
    return ok, r
