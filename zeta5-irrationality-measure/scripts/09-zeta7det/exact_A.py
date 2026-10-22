"""Rigorous upper bound for Fauzan's arithmetic constant A_M (zeta(5): K=40n, N=3n, h=37n, r=6), sharpening
his tail bound (5.16).  A_M <= I_out + int_3^X R/x^3 [exact rational, App. B method] + T(X) + 7 lam/M - ...
Tail (X in 40Z): R(x) = x F(x) + Q(x), F = 4lam + 2lam{x} - 12lam{alpha x} (Fauzan (5.12)), -1/2 <= Q <= 13/8 ((5.14)).
   int_X^oo F/x^2 = -lam/X + 2 int_X^oo G/x^3,  G = lam(f^2-f) - (6lam/alpha)(g^2-g) in [-lam/4, 3lam/(2alpha)]
   => int_X^oo R/x^3 <= -lam/X + (3lam/(2alpha) + 13/16)/X^2.
And the M-cutoff terms exactly as in Fauzan (5.17),(5.20):  + 7 lam/M - (2923/240 - 1/4)/M^2 + 32/M^3.
R is affine on each piece between consecutive breakpoints c x in Z, c in {2,2al,2lam,2H,4al,2(1-al),2(1+al)}
(Fauzan App. B); we verify affinity at a third point on every piece."""
import sys
from fractions import Fraction as Fr
from math import floor, ceil
from flint import arb, ctx

ctx.prec = 200
lam = Fr(37, 40); al = Fr(3, 40); H = 1 + 2 * al


def fl(x):
    return Fr(floor(x))


def ell(x, z):
    return fl(x - z) + fl(x + z) + 1


def frac(x):
    return x - fl(x)


def d0(u):
    return min(u, 1 - u)


def Gamma(x):
    T = fl(2 * H * x); s = H * x - T / 2; q = fl(2 * x); npl = (2 * x - q) / 2
    bps = sorted(set([Fr(0), Fr(1, 2), d0(frac(x)), d0(frac(al * x))]))
    I = Fr(0)
    for a, b in zip(bps, bps[1:]):
        if b <= a:
            continue
        z = (a + b) / 2
        l = ell(x, z); bb = 3 * ell(al * x, z)
        I += (T - bb) * (T + bb - l - 5) * (b - a)
    return I + s * (2 * T - q - 5) + max(s - npl, Fr(0))


def J(u):
    m = fl(2 * u)
    return m * u - m * (m + 1) / 4


def Nf(x):
    return 2 * lam * x * fl(x) - 12 * lam * x * fl(al * x) - 2 * J(lam * x)


def R(x):
    return -Gamma(x) - Nf(x)


def breakpoints(x0, x1):
    cs = [Fr(2), 2 * al, 2 * lam, 2 * H, 4 * al, 2 * (1 - al), 2 * (1 + al)]
    pts = {Fr(x0), Fr(x1)}
    for c in cs:
        k0 = ceil(x0 * c); k1 = floor(x1 * c)
        for k in range(k0, k1 + 1):
            v = Fr(k) / c
            if x0 < v < x1:
                pts.add(v)
    return sorted(pts)


def integral(x0, x1):
    pts = breakpoints(x0, x1)
    tot_exact = Fr(0)
    tot = arb(0)
    bad = 0
    for l, r in zip(pts, pts[1:]):
        x1_, x2_, x3_ = (2 * l + r) / 3, (l + 2 * r) / 3, (l + r) / 2
        R1, R2 = R(x1_), R(x2_)
        a = (R2 - R1) / (x2_ - x1_); b = R1 - a * x1_
        if R(x3_) != a * x3_ + b:
            bad += 1
        piece = a * (1 / l - 1 / r) + b / 2 * (1 / l ** 2 - 1 / r ** 2)
        if x1 <= 20:
            tot_exact += piece
        tot += arb(piece.numerator) / arb(piece.denominator)
    return tot, tot_exact, len(pts) - 1, bad


if __name__ == "__main__":
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    assert X % 40 == 0
    I320, ex, npc, bad = integral(3, 20)
    print("int_3^20 R/x^3 =", ex, "pieces", npc, "non-affine", bad,
          " == (5.18)?", ex == Fr(322437603634266857629, 7535670527041937280000))
    I20X, _, npc2, bad2 = integral(20, X)
    print(f"int_20^{X} R/x^3 = {I20X}  pieces {npc2} non-affine {bad2}")
    tail_up = -arb(lam.numerator) / lam.denominator / X + (arb(3) * arb(lam.numerator) / lam.denominator / (2 * arb(al.numerator) / al.denominator) + arb(13) / 16) / (X * X)
    print(f"tail bound int_{X}^oo <= {tail_up}")
    Iout = arb(127751) / 96000
    Astar = Iout + I320 + I20X + tail_up
    print("A*_new (upper bound) =", Astar, "   Fauzan A* =", float(Fr(9928298118277006344769, 7535670527041937280000)))
    for M in (200, 100000):
        AM = Astar + arb(7) * arb(37) / 40 / M - (arb(2923) / 240 - arb(1) / 4) / (M * M) + arb(32) / (M ** 3)
        print(f"A_{M}^new <= {AM}")
