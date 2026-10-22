"""
Rigorous sup-certificates for the kappa=2 trial fields (Litsey format), angle variable,
adaptive bisection with two bounds per piece.                           [rigorous: arb/acb]

x = cos(theta), theta in [0, pi];  s = cos(phi), phi in [0, pi/2].
Fields (design (20,3,1,2,1), kappa = 2, lambda = 0; Litsey's shifted fields):
  X(x) = 2/5 log|x| + 1/10 log(1-x) - 27/20 log(1+x^2) - 4T(p,x) - S(v,x)
  Y(s) = 3/20 log s + 1/20 log(1-s) + 2T(v,s)
  M(x) = -(9/10) log(1+x^2) - 2T(p,x)
For a piece [a,b] (angles) we accept if either
  (MVT)  F(m) + h sup_{[a,b]} |F'| <= target   (m = exact midpoint, h = half-length), or
  (SING) sum_i c_i log(sup_{[a,b]} dist_i) + G([a,b]) <= target, where F = sum_i c_i log dist_i + G,
         c_i > 0 and G is enclosed by ball evaluation (valid on every piece).
Pieces are bisected until accepted.  Ball enclosures of x, 1-x, |x|, sin use monotonicity of cos on [0,pi].
"""
import sys, json, math
from fractions import Fraction as Fr
from flint import arb, acb, ctx

ctx.prec = 256
PI = arb.pi()
HALFPI = PI / 2


def Q(fr):
    fr = Fr(fr)
    return arb(fr.numerator) / fr.denominator


def upper(x):
    return float(x.mid()) + float(x.rad())


def U(lo, hi):
    return arb.union(lo, hi)


class Piece:
    """angle piece [a,b] subset [0,pi] (a,b exact arb points)."""
    def __init__(self, a, b):
        self.a, self.b = a, b
        self.t = U(a, b)
        ca, cb = a.cos(), b.cos()                   # cos decreasing on [0,pi]
        lo = cb.lower(); hi = ca.upper()
        if hi > 1 or a.lower() <= 0: hi = arb(1)
        if lo < -1 or b.upper() >= PI.lower(): lo = arb(-1)
        self.x = U(lo, hi)
        # 1 - x = 2 sin^2(t/2), increasing in t
        oa = 2 * (a / 2).sin() ** 2; ob = 2 * (b / 2).sin() ** 2
        lo1 = oa.lower()
        if lo1 < 0: lo1 = arb(0)
        self.omx = U(lo1, ob.upper())
        # 1 + x = 2 cos^2(t/2), decreasing in t
        pa = 2 * (a / 2).cos() ** 2; pb = 2 * (b / 2).cos() ** 2
        lo2 = pb.lower()
        if lo2 < 0: lo2 = arb(0)
        self.opx = U(lo2, pa.upper())
        self.sin = self.t.sin()
        if self.sin.lower() < 0:
            self.sin = U(arb(0), self.sin.upper())
        self.absx_up = abs(self.x).upper()


class Tr:
    def __init__(self, cs, tails):
        self.c = [arb(int(c)) / 10**8 for c in cs]
        self.K = len(self.c)
        self.ck = [self.c[k] / (k + 1) for k in range(self.K)]
        self.L1 = sum((abs(c) for c in self.c), arb(0))
        self.L1k = sum((abs(self.c[k]) * (k + 1) for k in range(self.K)), arb(0))
        self.tails = []
        for zr, zi, rr, ri in tails:
            assert Fr(zr) ** 2 + Fr(zi) ** 2 < 1
            self.tails.append((acb(Q(zr), Q(zi)), acb(Q(rr), Q(ri))))

    def cheb_point(self, t):
        """powers of e^{it} re-anchored every 16 steps (fresh cos/sin of (k+1)t) to stop the
        box-wrapping growth (|cos t|+|sin t|)^k of repeated complex multiplication."""
        w = acb(t.cos(), t.sin()); pw = acb(1)
        A = arb(0); B = arb(0)
        for k in range(self.K):
            if (k + 1) % 16 == 0:
                kt = (k + 1) * t
                pw = acb(kt.cos(), kt.sin())
            else:
                pw = pw * w
            A += self.ck[k] * pw.real
            B += self.c[k] * pw.imag
        return A, B

    def cheb_piece(self, P):
        m = arb(((P.a + P.b) / 2).mid())
        h = arb.union(P.b - m, m - P.a).upper()
        A, B = self.cheb_point(m)
        return A + h * self.L1 * arb(0, 1), B + h * self.L1k * arb(0, 1)

    def Ttail(self, x, omx=None):
        """-1/2 Re sum r log(1-2xz+z^2); for real x in [-1,1].  Uses 1-2xz+z^2 = (1-z)^2 + 2z(1-x)."""
        val = arb(0)
        for z, r in self.tails:
            if omx is None:
                arg = 1 - 2 * x * z + z * z
            else:
                arg = (1 - z) ** 2 + 2 * z * omx
            val -= (r * arg.log()).real / 2
        return val

    def dTtail_dx(self, x, omx=None):
        val = arb(0)
        for z, r in self.tails:
            arg = (1 - 2 * x * z + z * z) if omx is None else ((1 - z) ** 2 + 2 * z * omx)
            val += (r * z / arg).real
        return val

    def S(self, x):
        acc = arb(0)
        for k in range(self.K - 1, -1, -1):
            acc = acc * x + self.ck[k]
        val = acc * x
        for z, r in self.tails:
            val -= (r * (1 - x * z).log()).real
        return val

    def dS(self, x):
        acc = arb(0)
        for k in range(self.K - 1, -1, -1):
            acc = acc * x + self.c[k]
        val = acc
        for z, r in self.tails:
            val += (r * z / (1 - x * z)).real
        return val

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

    def first(self):
        val = self.c[0]
        for z, r in self.tails:
            val += (r * z).real
        return val


C25, C110, C2720, C320, C120, C910 = Q(Fr(2, 5)), Q(Fr(1, 10)), Q(Fr(27, 20)), Q(Fr(3, 20)), Q(Fr(1, 20)), Q(Fr(9, 10))


# ---------------------------------------------------------------- X
def X_point(p, v, t):
    x = t.cos(); omx = 2 * (t / 2).sin() ** 2
    A, B = p.cheb_point(t)
    return C25 * abs(x).log() + C110 * omx.log() - C2720 * (1 + x * x).log() - 4 * (A + p.Ttail(x, omx)) - v.S(x)


def X_dball(p, v, P):
    A, B = p.cheb_piece(P)
    x, sn = P.x, P.sin
    dT = -B - sn * p.dTtail_dx(x, P.omx)
    dS = -sn * v.dS(x)
    dW = -sn * (C25 / x - C110 / P.omx - C2720 * 2 * x / (1 + x * x))
    return -4 * dT - dS + dW


def X_sing(p, v, P):
    A, B = p.cheb_piece(P)
    G = -C2720 * (1 + P.x * P.x).log() - 4 * (A + p.Ttail(P.x, P.omx)) - v.S(P.x)
    return C25 * P.absx_up.log() + C110 * P.omx.upper().log() + G


# ---------------------------------------------------------------- Y  (phi in [0, pi/2])
def Y_point(v, t):
    s = t.cos(); oms = 2 * (t / 2).sin() ** 2
    A, B = v.cheb_point(t)
    return C320 * s.log() + C120 * oms.log() + 2 * (A + v.Ttail(s, oms))


def Y_dball(v, P):
    A, B = v.cheb_piece(P)
    s, sn = P.x, P.sin
    dT = -B - sn * v.dTtail_dx(s, P.omx)
    return -sn * (C320 / s - C120 / P.omx) + 2 * dT


def Y_sing(v, P):
    A, B = v.cheb_piece(P)
    G = 2 * (A + v.Ttail(P.x, P.omx))
    return C320 * P.x.upper().log() + C120 * P.omx.upper().log() + G


# ---------------------------------------------------------------- M
def M_point(p, t):
    x = t.cos(); omx = 2 * (t / 2).sin() ** 2
    A, B = p.cheb_point(t)
    return -C910 * (1 + x * x).log() - 2 * (A + p.Ttail(x, omx))


def M_dball(p, P):
    A, B = p.cheb_piece(P)
    x, sn = P.x, P.sin
    dT = -B - sn * p.dTtail_dx(x, P.omx)
    return C910 * sn * 2 * x / (1 + x * x) - 2 * dT


def M_sing(p, P):
    A, B = p.cheb_piece(P)
    return -C910 * (1 + P.x * P.x).log() - 2 * (A + p.Ttail(P.x, P.omx))


def certify_field(Fp, Fd, Fs, a, b, target, maxdepth=70, maxpieces=3_000_000):
    stack = [(a, b, 0)]
    worst = -1e300; n = 0
    while stack:
        lo, hi, dep = stack.pop()
        P = Piece(lo, hi)
        m = arb(((lo + hi) / 2).mid())
        h = arb.union(hi - m, m - lo).upper()
        ok = False
        try:
            fm = Fp(m); d = Fd(P)
            if fm.is_finite() and d.is_finite():
                ub = upper(fm + h * abs(d).upper())
                if ub <= target:
                    ok = True
        except Exception:
            pass
        if not ok:
            try:
                sb = Fs(P)
                if sb.is_finite():
                    ub = upper(sb)
                    if ub <= target:
                        ok = True
            except Exception:
                pass
        if ok:
            worst = max(worst, ub); n += 1
            if n > maxpieces:
                return False, worst, n, 'too many pieces'
            continue
        if dep >= maxdepth:
            return False, None, n, (float(lo.mid()), float(hi.mid()))
        stack.append((lo, m, dep + 1)); stack.append((m, hi, dep + 1))
    return True, worst, n, None


def cert_all(p, v, MX, MY, RM, verbose=True):
    res = []
    # X on theta in [0, pi]; split at pi/2 (x=0) using exact rational bracketing points around pi/2
    a0, a1 = arb(0), arb(Fr(157079, 100000).numerator) / 100000     # < pi/2
    a2 = arb(157080) / 100000                                          # > pi/2
    Xp = lambda t: X_point(p, v, t); Xd = lambda P: X_dball(p, v, P); Xs = lambda P: X_sing(p, v, P)
    res.append(('X [0, 1.57079]',) + certify_field(Xp, Xd, Xs, a0, a1, MX))
    res.append(('X [1.57079, 1.57080]',) + certify_field(Xp, Xd, Xs, a1, a2, MX))
    res.append(('X [1.57080, pi]',) + certify_field(Xp, Xd, Xs, a2, arb(PI.upper()), MX))
    Yp = lambda t: Y_point(v, t); Yd = lambda P: Y_dball(v, P); Ys = lambda P: Y_sing(v, P)
    res.append(('Y [0, 1.57079]',) + certify_field(Yp, Yd, Ys, a0, a1, MY))
    res.append(('Y [1.57079, pi/2]',) + certify_field(Yp, Yd, Ys, a1, arb(HALFPI.lower()), MY))
    if RM is not None:
        tg = Q(Fr(401, 2500)).acos()
        Mp = lambda t: M_point(p, t); Md = lambda P: M_dball(p, P); Ms = lambda P: M_sing(p, P)
        res.append(('M [0, tg]',) + certify_field(Mp, Md, Ms, a0, arb(tg.upper()), RM))
        res.append(('M [pi-tg, pi]',) + certify_field(Mp, Md, Ms, arb((PI - tg).lower()), arb(PI.upper()), RM))
    ok = all(r[1] for r in res)
    if verbose:
        for r in res:
            print(r, flush=True)
    return ok, res
