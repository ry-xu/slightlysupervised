"""15-redteam-mu5: independent exact evaluation of Fauzan's arithmetic constant pieces.

R(x) = -Gamma(x) - N(x)   (Fauzan (5.4)-(5.6)), implemented directly from the definitions:
  l(x,z) = floor(x-z) + floor(x+z) + 1,   b(x,z) = 3 l(alpha x, z)
  T = floor(2Hx), s = Hx - T/2, q = floor(2x), n+ = (2x-q)/2
  Gamma(x) = int_0^{1/2} (T-b)(T+b-l-5) dz + s(2T-q-5) + (s-n+)_+
  J(u) = m u - m(m+1)/4, m = floor(2u);  N(x) = 2 lam x floor(x) - 12 lam x floor(alpha x) - 2 J(lam x)
The z-integral is done EXACTLY by splitting (0,1/2) at the jump points of the integrand (not via Fauzan's (B.1)).

Breakpoints in x (derived here): R is affine between consecutive points of
  (1/2)Z u (20/3)Z u (10/23)Z u (20/37)Z u (10/3)Z u (40/43)Z
(floors of x, alpha x, 2Hx, 2x, 2 lam x; z-jump order changes f=g, f=1-g, f=1/2, g=1/2; tau=sigma).
I use Fauzan's (B.2) set {2, 2al, 2lam, 2H, 4al, 2(1-al), 2(1+al)} which contains all of these, and verify
affinity on every piece with a third (midpoint) evaluation.

Tasks:
  (1) int_3^20 R/x^3 exactly; compare (5.18) and Table 3.
  (2) int_20^X R/x^3 exactly for X = 4000 (compare 09: -0.0598432941).
  (3) periodicity: R(40k+y) = R(y) + k (R(y+40)-R(y)) [k^2 coefficient cancels, algebra in notes];
      verify R(y+40) - R(y) = 40 F(y) on all pieces of [20,60)  =>  Q = R - xF is 40-periodic on [20,oo).
      Then compute sup/inf of Q over a period EXACTLY (Q is quadratic on each piece).
  (4) check Fauzan's closed form (5.13) for Q against R - xF at many exact points.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from math import floor

LAM = Fr(37, 40); AL = Fr(3, 40); H = 1 + 2 * AL   # H = 23/20
HALF = Fr(1, 2)


def fl(x):
    return x.numerator // x.denominator


def frac(x):
    return x - fl(x)


def ell(x, z):
    return fl(x - z) + fl(x + z) + 1


def zint(x):
    """exact int_0^{1/2} (T-b)(T+b-l-5) dz for rational x"""
    T = fl(2 * H * x)
    ax = AL * x
    cuts = {Fr(0), HALF}
    for u in (x, ax):
        f = frac(u)
        for z in (f, 1 - f):
            if 0 < z < HALF:
                cuts.add(z)
    cuts = sorted(cuts)
    tot = Fr(0)
    for z0, z1 in zip(cuts[:-1], cuts[1:]):
        zm = (z0 + z1) / 2
        l = ell(x, zm)
        b = 3 * ell(ax, zm)
        tot += (z1 - z0) * (T - b) * (T + b - l - 5)
    return tot


def Gamma(x):
    T = fl(2 * H * x)
    s = H * x - Fr(T, 2)
    q = fl(2 * x)
    npl = (2 * x - q) / 2
    return zint(x) + s * (2 * T - q - 5) + max(s - npl, Fr(0))


def Jf(u):
    m = fl(2 * u)
    return m * u - Fr(m * (m + 1), 4)


def Nf(x):
    return 2 * LAM * x * fl(x) - 12 * LAM * x * fl(AL * x) - 2 * Jf(LAM * x)


def R(x):
    return -Gamma(x) - Nf(x)


def F(x):
    return 4 * LAM + 2 * LAM * frac(x) - 12 * LAM * frac(AL * x)


def Q_fauzan(x):
    """Fauzan (5.13) with (B.1):  Q = [tau(tau-sigma) - (tau-sigma)_+ + eta(1-eta)]/2 + 9 int B^2 - 3 int AB"""
    f, g = frac(x), frac(AL * x)
    tau, sigma, eta = frac(2 * H * x), frac(2 * x), frac(2 * LAM * x)
    d0 = lambda u: min(u, 1 - u)
    e = lambda u: 1 if u <= HALF else -1
    intB2 = d0(g) * (1 - 2 * d0(g))
    intAB = e(f) * e(g) * (min(d0(f), d0(g)) - 2 * d0(f) * d0(g))
    return (tau * (tau - sigma) - max(tau - sigma, Fr(0)) + eta * (1 - eta)) / 2 + 9 * intB2 - 3 * intAB


CSET = [Fr(2), 2 * AL, 2 * LAM, 2 * H, 4 * AL, 2 * (1 - AL), 2 * (1 + AL)]


def breakpoints(lo, hi):
    pts = {Fr(lo), Fr(hi)}
    for c in CSET:
        # x = k / c in (lo, hi)
        k0 = fl(Fr(lo) * c) + 1
        k = k0
        while Fr(k) / c < hi:
            pts.add(Fr(k) / c)
            k += 1
    return sorted(pts)


def pieces(lo, hi, check=True):
    """yield (l, r, slope, intercept) of the affine R on each piece; verifies affinity at the midpoint."""
    bps = breakpoints(lo, hi)
    for l, r in zip(bps[:-1], bps[1:]):
        x1, x2 = (2 * l + r) / 3, (l + 2 * r) / 3
        R1, R2 = R(x1), R(x2)
        a = (R2 - R1) / (x2 - x1)
        b = R1 - a * x1
        if check:
            xm = (l + r) / 2
            assert R(xm) == a * xm + b, f"R not affine on [{l},{r}]"
            # extra: points near the ends
            for xx in (l + (r - l) / 97, r - (r - l) / 89):
                assert R(xx) == a * xx + b, f"R not affine near the ends of [{l},{r}]"
        yield l, r, a, b


def integral(lo, hi, exact=True):
    tot = Fr(0)
    n = 0
    for l, r, a, b in pieces(lo, hi):
        tot += a * (1 / l - 1 / r) + b / 2 * (1 / (l * l) - 1 / (r * r))
        n += 1
    return tot, n


if __name__ == "__main__":
    t0 = time.time()
    # ---- (1) [3,20] and Table 3
    I320, n = integral(3, 20)
    paper518 = Fr(322437603634266857629, 7535670527041937280000)
    print(f"(1) int_3^20 R/x^3 = {I320}  = {float(I320):.15f}   pieces = {n}")
    print(f"    equals Fauzan (5.18) exactly: {I320 == paper518}")
    table3 = {3: Fr(26807, 161280), 4: Fr(37383, 704000), 5: Fr(37523, 8236800), 6: Fr(-120923, 6552000),
              7: Fr(-10025233, 356428800), 8: Fr(-110029309, 3348864000), 9: Fr(-2278419487, 64465632000),
              10: Fr(-282415081, 7724640000), 11: Fr(-3236921227, 87524236800), 12: Fr(-3350220001, 90899827200),
              13: Fr(7424224373, 2217983040000), 14: Fr(16775764609, 955086612480), 15: Fr(369043847, 30810528000),
              16: Fr(651380108633, 86461373856000), 17: Fr(472851276229, 119820121344000),
              18: Fr(4930060867, 4724197793280), 19: Fr(-15199801, 11563552000)}
    bad = []
    for j, v in table3.items():
        mine, _ = integral(j, j + 1)
        if mine != v:
            bad.append((j, mine, v))
    print(f"    Table 3 unit integrals all match: {not bad}  {bad[:3]}")
    print(f"    sum of Table 3 == (5.18): {sum(table3.values()) == paper518}")
    print(f"    [t = {time.time()-t0:.1f}s]")

    # ---- (3) periodicity of Q = R - xF on [20, oo)
    t1 = time.time()
    # (3a) second difference in k vanishes (k^2-coefficient cancellation) -- numeric confirmation of the algebra
    import random
    random.seed(15)
    ok2 = all(R(y + 80) - 2 * R(y + 40) + R(y) == 0 and R(y + 4000) - R(y) == 100 * (R(y + 40) - R(y))
              for y in (Fr(random.randint(20 * 10 ** 6, 60 * 10 ** 6), 10 ** 6) for _ in range(2000)))
    print(f"(3a) R(y+80)-2R(y+40)+R(y)=0 and R(y+4000)-R(y)=100(R(y+40)-R(y)) at 2000 random y: {ok2}")
    # (3b) R(y+40) - R(y) - 40F(y) == 0 on every piece of [20,60): affine on pieces (breakpoint set is
    #      20-periodic), so 3 points per piece prove it on the whole piece.
    bps = breakpoints(20, 60)
    okD = True
    for l, r in zip(bps[:-1], bps[1:]):
        for w in (Fr(1, 3), Fr(1, 2), Fr(2, 3)):
            y = l + (r - l) * w
            if R(y + 40) - R(y) - 40 * F(y) != 0:
                okD = False
    print(f"(3b) R(y+40)-R(y) = 40F(y) on all {len(bps)-1} pieces of [20,60): {okD}")
    # (3c) exact range of Q on a period: Q = R - xF is quadratic on each piece
    Qmax, Qmin = None, None
    argmax = argmin = None
    nQ = 0
    for l, r, a, b in pieces(20, 60):
        # F affine on piece: F = p x + c0 ; get from two interior points
        x1, x2 = (2 * l + r) / 3, (l + 2 * r) / 3
        p = (F(x2) - F(x1)) / (x2 - x1); c0 = F(x1) - p * x1
        # Q(x) = a x + b - x (p x + c0) = -p x^2 + (a - c0) x + b   on the closed piece
        cand = [l, r]
        if p != 0:
            xv = (a - c0) / (2 * p)
            if l < xv < r:
                cand.append(xv)
        for x in cand:
            v = -p * x * x + (a - c0) * x + b
            if Qmax is None or v > Qmax:
                Qmax, argmax = v, x
            if Qmin is None or v < Qmin:
                Qmin, argmin = v, x
        # spot-check against direct evaluation and Fauzan's closed form at interior points
        for w in (Fr(1, 4), Fr(3, 4)):
            x = l + (r - l) * w
            qd = R(x) - x * F(x)
            assert qd == -p * x * x + (a - c0) * x + b
            if qd != Q_fauzan(x):
                nQ += 1
    print(f"(3c) exact range of Q on [20,60] (closure of pieces): [{Qmin} = {float(Qmin):.6f} at x={float(argmin):.4f},"
          f" {Qmax} = {float(Qmax):.6f} at x={float(argmax):.4f}]")
    print(f"     Fauzan's (5.14) [-1/2, 13/8] holds on the period: {Qmin >= Fr(-1,2) and Qmax <= Fr(13,8)}")
    print(f"(4)  Fauzan closed form (5.13) == R - xF at all checked interior points: {nQ == 0} (mismatches {nQ})")
    print(f"    [t = {time.time()-t1:.1f}s]")

    # ---- (2) int_20^4000 exactly (pieces), accumulated with fmpq for speed
    from flint import fmpq
    t2 = time.time()
    tot = fmpq(0)
    n = 0
    for l, r, a, b in pieces(20, 4000, check=True):
        L, Rr, Aa, Bb = fmpq(l.numerator, l.denominator), fmpq(r.numerator, r.denominator), \
            fmpq(a.numerator, a.denominator), fmpq(b.numerator, b.denominator)
        tot += Aa * (1 / L - 1 / Rr) + Bb / 2 * (1 / (L * L) - 1 / (Rr * Rr))
        n += 1
    from flint import arb, ctx
    ctx.prec = 200
    v = arb(tot)
    print(f"(2) int_20^4000 R/x^3 = {v}  pieces = {n}   [t = {time.time()-t2:.1f}s]")
    print(f"    denominator digits of exact value: {len(str(tot.q))}")
    open('./data/15-redteam-mu5/int_20_4000.txt', 'w').write(
        f"{tot.p}\n{tot.q}\n")
