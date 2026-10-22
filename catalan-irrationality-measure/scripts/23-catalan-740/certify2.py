"""
Rigorous sup-certificates for the kappa=2 trial fields in the ANGLE variable (no wrapping):
  x = cos(theta), theta in [0, pi];  s = cos(phi), phi in [0, pi/2].
Chebyshev parts via omega = e^{i theta}:  T_k(x) = Re omega^k,  d/dtheta [T_k(cos theta)/k] = -sin(k theta).
Ball radii grow only linearly in k.                                     [rigorous: arb/acb]

Fields (Litsey's shifted fields, design (20,3,1,2,1), kappa=2, lambda=0):
  X(x) = 2/5 log|x| + 1/10 log(1-x) - 27/20 log(1+x^2) - 4T(p,x) - S(v,x)
  Y(s) = 3/20 log s + 1/20 log(1-s) + 2T(v,s)
  M(x) = -(9/10) log(1+x^2) - 2T(p,x)          (marked residual, |x| >= 401/2500)
A piece [t1,t2] (angle) is accepted if  F(t_m) + h * sup_{[t1,t2]} |dF/dt| <= target.
Singular endpoints (x=0, x=1, s=0, s=1) use monotone majorants c*log(dist) + G(ball).
"""
import sys, json, math
from fractions import Fraction as Fr
from flint import arb, acb, ctx

ctx.prec = 256
PI = arb.pi()


def Q(fr):
    fr = Fr(fr)
    return arb(fr.numerator) / fr.denominator


def upper(x):
    return float(x.mid()) + float(x.rad())


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

    # Chebyshev finite parts in the angle variable; t is an arb (ball ok)
    def cheb_parts(self, t):
        """returns enclosures of (A, B) over the ball t, with A = sum c_k/k cos(k t), B = sum c_k sin(k t).
        Evaluated at the exact midpoint by complex powers, plus a Lipschitz term
        |A(t)-A(m)| <= h*sum|c_k|,  |B(t)-B(m)| <= h*sum k|c_k|   (h = radius of t)."""
        m = arb(t.mid())
        h = arb(t.rad())
        w = acb(m.cos(), m.sin())          # e^{im}
        pw = acb(1)
        A = arb(0); B = arb(0)
        for k in range(self.K):
            pw = pw * w
            A += self.ck[k] * pw.real
            B += self.c[k] * pw.imag
        if h != 0:
            A += h * self.L1 * arb(0, 1)
            B += h * self.L1k * arb(0, 1)
        return A, B

    def Ttail(self, x):
        val = arb(0)
        for z, r in self.tails:
            val -= (r * (1 - 2 * x * z + z * z).log()).real / 2
        return val

    def dTtail_dx(self, x):
        val = arb(0)
        for z, r in self.tails:
            val += (r * z / (1 - 2 * x * z + z * z)).real
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
        """u_1 enclosure."""
        val = self.c[0]
        for z, r in self.tails:
            val += (r * z).real
        return val


# ------------------------- X in theta (x = cos theta) -------------------------
def X_theta(p, v, t, parts=('W', 'T', 'S')):
    x = t.cos()
    A, B = p.cheb_parts(t)
    T = A + p.Ttail(x)
    val = -4 * T - v.S(x)
    if 'W' in parts:
        val += Q(Fr(2, 5)) * abs(x).log() + Q(Fr(1, 10)) * (1 - x).log() - Q(Fr(27, 20)) * (1 + x * x).log()
    return val


def dX_theta(p, v, t):
    x = t.cos(); sn = t.sin()
    A, B = p.cheb_parts(t)
    # d/dt[-4T] = -4*(-B) ... : d/dt A = -B ; d/dt Ttail(cos t) = -sin t * dTtail_dx
    dT = -B - sn * p.dTtail_dx(x)
    dS = -sn * v.dS(x)
    dW = -sn * (Q(Fr(2, 5)) / x - Q(Fr(1, 10)) / (1 - x) - Q(Fr(27, 20)) * 2 * x / (1 + x * x))
    return -4 * dT - dS + dW


def M_theta(p, t):
    x = t.cos()
    A, B = p.cheb_parts(t)
    return -Q(Fr(9, 10)) * (1 + x * x).log() - 2 * (A + p.Ttail(x))


def dM_theta(p, t):
    x = t.cos(); sn = t.sin()
    A, B = p.cheb_parts(t)
    dT = -B - sn * p.dTtail_dx(x)
    return -Q(Fr(9, 10)) * (-sn) * 2 * x / (1 + x * x) - 2 * dT


# ------------------------- Y in phi (s = cos phi) -------------------------
def Y_phi(v, t):
    s = t.cos()
    A, B = v.cheb_parts(t)
    return Q(Fr(3, 20)) * s.log() + Q(Fr(1, 20)) * (1 - s).log() + 2 * (A + v.Ttail(s))


def dY_phi(v, t):
    s = t.cos(); sn = t.sin()
    A, B = v.cheb_parts(t)
    dT = -B - sn * v.dTtail_dx(s)
    return -sn * (Q(Fr(3, 20)) / s - Q(Fr(1, 20)) / (1 - s)) + 2 * dT


def tball(a, b):
    """arb ball containing the angle interval [a,b] (a,b arb, b>a)."""
    m = (a + b) / 2
    h = (b - a) / 2
    return m + (h.mid() * (1 + 1e-12) + h.rad()) * arb(0, 1) if False else arb.union(a, b)


def cert_sup_angle(F, dF, a, b, target, maxdepth=64):
    """sup over angle interval [a,b] (arb endpoints, exact-ish) of F <= target."""
    stack = [(a, b, 0)]
    worst = -1e300; n = 0
    while stack:
        lo, hi, dep = stack.pop()
        m = (lo + hi) / 2
        m = arb(m.mid())                       # exact midpoint representative
        h = arb.union(hi - m, m - lo).upper()  # >= max distance to m
        fm = F(m)
        d = dF(arb.union(lo, hi))
        bnd = fm + h * abs(d).upper()
        ub = upper(bnd)
        if ub <= target and d.is_finite() and fm.is_finite():
            worst = max(worst, ub); n += 1
            continue
        if dep >= maxdepth:
            return False, ub, n, (float(lo.mid()), float(hi.mid()))
        stack.append((lo, m, dep + 1)); stack.append((m, hi, dep + 1))
    return True, worst, n, None


def sing_angle(G, c, dist_of_ball, a, b, toward_a, target, ratio=0.25, mindist=1e-40):
    """Prove c*log(dist) + G <= target on the angle interval between the singular endpoint and the
    other end.  dist_of_ball(ball) returns an upper bound (arb) for the distance function on the ball.
    toward_a: True if the singular endpoint is a.  Pieces shrink geometrically toward the singular end.
    Last piece: G on the closed ball, c*log(dist upper)."""
    worst = -1e300; n = 0
    L = b - a
    frac = arb(1)
    while True:
        lo_frac = frac * ratio
        if upper(lo_frac * L) < mindist:
            lo_frac = arb(0)
        if toward_a:
            seg = arb.union(a + lo_frac * L, a + frac * L)
        else:
            seg = arb.union(b - frac * L, b - lo_frac * L)
        g = G(seg)
        dist = dist_of_ball(seg)
        bnd = c * dist.upper().log() + g
        ub = upper(bnd)
        if not (ub <= target and g.is_finite()):
            return False, ub, n, (float(seg.lower()), float(seg.upper()))
        worst = max(worst, ub); n += 1
        if lo_frac == 0:
            break
        frac = lo_frac
    return True, worst, n, None


def cert_all(p, v, MX, MY, RM, verbose=True, eps_th=None):
    """Full certificate.  Angle splitting:
       X: theta in [0, pi]: x=1 at theta=0 (singular, c=1/10, dist=1-x), x=0 at pi/2 (c=2/5, dist=|x|),
          x=-1 at theta=pi (regular).
       Y: phi in (0, pi/2): s=1 at phi=0 (c=1/20), s=0 at pi/2 (c=3/20).
       M: theta in [0, theta_g] U [pi-theta_g, pi], theta_g = arccos(401/2500).
    """
    res = []
    half = PI / 2
    d0 = arb(1) / 200      # angular half-width of the singular window at x=0 / s=0
    d1 = arb(1) / 100      # angular window at x=1 / s=1
    Xf = lambda t: X_theta(p, v, t); dXf = lambda t: dX_theta(p, v, t)
    # singular pieces for X
    G0 = lambda t: X_theta(p, v, t, parts=('T', 'S')) + Q(Fr(1, 10)) * (1 - t.cos()).log() - Q(Fr(27, 20)) * (1 + t.cos() ** 2).log()
    G1 = lambda t: X_theta(p, v, t, parts=('T', 'S')) + Q(Fr(2, 5)) * abs(t.cos()).log() - Q(Fr(27, 20)) * (1 + t.cos() ** 2).log()
    distx0 = lambda t: abs(t.cos())          # |x|
    distx1 = lambda t: 1 - t.cos()           # 1-x
    res.append(('X near x=1',) + sing_angle(G1, Q(Fr(1, 10)), distx1, arb(0), d1, True, MX))
    res.append(('X [d1, pi/2-d0]',) + cert_sup_angle(Xf, dXf, d1, half - d0, MX))
    res.append(('X near x=0 (+)',) + sing_angle(G0, Q(Fr(2, 5)), distx0, half - d0, half, False, MX))
    res.append(('X near x=0 (-)',) + sing_angle(G0, Q(Fr(2, 5)), distx0, half, half + d0, True, MX))
    res.append(('X [pi/2+d0, pi]',) + cert_sup_angle(Xf, dXf, half + d0, PI, MX))
    Yf = lambda t: Y_phi(v, t); dYf = lambda t: dY_phi(v, t)
    H1 = lambda t: Q(Fr(3, 20)) * t.cos().log() + 2 * (v.cheb_parts(t)[0] + v.Ttail(t.cos()))
    H0 = lambda t: Q(Fr(1, 20)) * (1 - t.cos()).log() + 2 * (v.cheb_parts(t)[0] + v.Ttail(t.cos()))
    dists1 = lambda t: 1 - t.cos()
    dists0 = lambda t: t.cos()
    res.append(('Y near s=1',) + sing_angle(H1, Q(Fr(1, 20)), dists1, arb(0), d1, True, MY))
    res.append(('Y [d1, pi/2-d0]',) + cert_sup_angle(Yf, dYf, d1, half - d0, MY))
    res.append(('Y near s=0',) + sing_angle(H0, Q(Fr(3, 20)), dists0, half - d0, half, False, MY))
    if RM is not None:
        tg = Q(Fr(401, 2500)).acos()   # angle of x = g ; |x| >= g  <=> theta in [0,tg] U [pi-tg, pi]
        Mf = lambda t: M_theta(p, t); dMf = lambda t: dM_theta(p, t)
        res.append(('M [0, tg]',) + cert_sup_angle(Mf, dMf, arb(0), arb(tg.upper()), RM))
        res.append(('M [pi-tg, pi]',) + cert_sup_angle(Mf, dMf, arb((PI - tg).lower()), PI, RM))
    ok = all(r[1] for r in res)
    if verbose:
        for r in res:
            print(r, flush=True)
    return ok, res
