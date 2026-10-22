"""
Second-order certifier (adds a Taylor test to certify3).                      [rigorous: arb/acb]

On a piece [a,b] (angle t), m = exact midpoint, h = half-length, accept if one of
  (T2)   F(m) + |F'(m)| h + (1/2) h^2 sup_{[a,b]} |F''|  <= cap
  (MVT)  F(m) + h sup_{[a,b]} |F'|                         <= cap     (certify3)
  (SING) sum c_i log(sup dist_i) + G([a,b])                <= cap     (certify3)
F'(m) is a point enclosure (exact midpoint).  sup|F''| on the piece: the Chebyshev part
d^2/dt^2 [sum (c_k/k) cos kt] = -sum k c_k cos kt is enclosed by its midpoint value plus
h * sum k^2 |c_k|; the elementary and tail parts by ball evaluation on the piece.
Valid by Taylor's theorem with Lagrange remainder (F is C^2 on pieces avoiding the singular points;
pieces where any enclosure is non-finite are not accepted by T2).
"""
import sys
from fractions import Fraction as Fr
from flint import arb, acb, ctx
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import certify3 as c3
from certify3 import Q, upper, Piece, PI, HALFPI, C25, C110, C2720, C320, C120, C910


def add_L2(tr):
    tr.L2k = sum((abs(tr.c[k]) * (k + 1) ** 2 for k in range(tr.K)), arb(0))


def cheb_point3(tr, t):
    """A = sum c_k/k cos kt, B = sum c_k sin kt, Bp = sum k c_k cos kt  at a point t."""
    w = acb(t.cos(), t.sin()); pw = acb(1)
    A = arb(0); B = arb(0); Bp = arb(0)
    for k in range(tr.K):
        if (k + 1) % 16 == 0:
            kt = (k + 1) * t
            pw = acb(kt.cos(), kt.sin())
        else:
            pw = pw * w
        A += tr.ck[k] * pw.real
        B += tr.c[k] * pw.imag
        Bp += tr.c[k] * (k + 1) * pw.real
    return A, B, Bp


def tail_d2(tr, x):
    """d^2/dx^2 of -1/2 Re sum r log(1-2xz+z^2) = Re sum r 2 z^2/(1-2xz+z^2)^2 (x real ball)."""
    val = arb(0)
    for z, r in tr.tails:
        q = 1 - 2 * x * z + z * z
        val += (r * 2 * z * z / (q * q)).real
    return val


def S_d2(tr, x):
    """second derivative of S(v,x) = sum c_k x^k/k - Re sum r log(1-xz): sum (k-1)c_k x^{k-2} + Re sum r z^2/(1-xz)^2."""
    acc = arb(0)
    for k in range(tr.K - 1, 0, -1):          # term index k (0-based) -> power k-1 coefficient k*c_{k+1}... use explicit
        acc = acc * x + tr.c[k] * k
    val = acc
    for z, r in tr.tails:
        q = 1 - x * z
        val += (r * z * z / (q * q)).real
    return val


# ---------------------------------------------------------------- Y(phi)
def Y_point_d1(v, t):
    s = t.cos(); sn = t.sin(); oms = 2 * (t / 2).sin() ** 2
    A, B, Bp = cheb_point3(v, t)
    val = C320 * s.log() + C120 * oms.log() + 2 * (A + v.Ttail(s, oms))
    d1 = -sn * (C320 / s - C120 / oms) + 2 * (-B - sn * v.dTtail_dx(s, oms))
    return val, d1


def Y_d2ball(v, P):
    m = arb(((P.a + P.b) / 2).mid())
    h = arb.union(P.b - m, m - P.a).upper()
    _, _, Bp = cheb_point3(v, m)
    Bp_ball = Bp + h * v.L2k * arb(0, 1)
    s, sn = P.x, P.sin
    cs = s                     # cos(phi) = s
    Wp = C320 / s - C120 / P.omx
    Wpp = -C320 / (s * s) - C120 / (P.omx * P.omx)
    Tp = v.dTtail_dx(s, P.omx)
    Tpp = tail_d2(v, s)
    d2W = -cs * Wp + sn * sn * Wpp
    d2T = -Bp_ball + (-cs * Tp + sn * sn * Tpp)
    return d2W + 2 * d2T


# ---------------------------------------------------------------- X(theta)
def X_point_d1(p, v, t):
    x = t.cos(); sn = t.sin(); omx = 2 * (t / 2).sin() ** 2
    A, B, Bp = cheb_point3(p, t)
    val = C25 * abs(x).log() + C110 * omx.log() - C2720 * (1 + x * x).log() - 4 * (A + p.Ttail(x, omx)) - v.S(x)
    Wp = C25 / x - C110 / omx - C2720 * 2 * x / (1 + x * x)
    d1 = -sn * Wp - 4 * (-B - sn * p.dTtail_dx(x, omx)) - (-sn * v.dS(x))
    return val, d1


def X_d2ball(p, v, P):
    m = arb(((P.a + P.b) / 2).mid())
    h = arb.union(P.b - m, m - P.a).upper()
    _, _, Bp = cheb_point3(p, m)
    Bp_ball = Bp + h * p.L2k * arb(0, 1)
    x, sn = P.x, P.sin
    cs = x
    Wp = C25 / x - C110 / P.omx - C2720 * 2 * x / (1 + x * x)
    Wpp = -C25 / (x * x) - C110 / (P.omx * P.omx) - C2720 * (2 * (1 + x * x) - 4 * x * x) / ((1 + x * x) ** 2)
    Tp = p.dTtail_dx(x, P.omx); Tpp = tail_d2(p, x)
    Sp = v.dS(x); Spp = S_d2(v, x)
    d2W = -cs * Wp + sn * sn * Wpp
    d2T = -Bp_ball + (-cs * Tp + sn * sn * Tpp)
    d2S = -cs * Sp + sn * sn * Spp
    return d2W - 4 * d2T - d2S


def certify_field2(Fpd, F2, Fd, Fs, a, b, target, maxdepth=70, maxpieces=20_000_000):
    stack = [(a, b, 0)]
    worst = -1e300; n = 0; nT2 = 0
    while stack:
        lo, hi, dep = stack.pop()
        P = Piece(lo, hi)
        m = arb(((lo + hi) / 2).mid())
        h = arb.union(hi - m, m - lo).upper()
        ok = False; ub = None
        try:
            fm, d1 = Fpd(m)
            if fm.is_finite() and d1.is_finite():
                d2 = F2(P)
                if d2.is_finite():
                    ub = upper(fm + abs(d1).upper() * h + h * h * abs(d2).upper() / 2)
                    if ub <= target:
                        ok = True; nT2 += 1
        except Exception:
            pass
        if not ok:
            try:
                fm, _ = Fpd(m)
                d = Fd(P)
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
