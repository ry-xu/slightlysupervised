"""15-redteam-mu5: I_out from the definitions (5.8)-(5.10), rigorous tail bounds, A_M, and a high-precision
evaluation of the exact tail via periodicity + Hurwitz zeta (cross-check only).

Uses data/15-redteam-mu5/int_20_4000.txt (exact rational produced by exact_R.py).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import arb, fmpq, ctx
from exact_R import R, F, LAM, AL, breakpoints, pieces, fl, frac

ctx.prec = 200
A = lambda q: arb(fmpq(Fr(q).numerator, Fr(q).denominator))
pos = lambda u: u if u > 0 else Fr(0)

# ---------------- I_out directly from (5.8)-(5.10) ----------------
def R0(y):
    al = AL
    if Fr(1, 3) < y < Fr(1, 2):
        return 8 - 9 * y - 8 * al - 5 * min(al, 1 - 2 * y) - 5 * pos(1 + al - 3 * y)
    if Fr(1, 2) < y < 1:
        return 7 * (1 - y) - 6 * min(al, 1 - y) - 6 * pos(1 + al - 2 * y) + pos(1 + 4 * al - 2 * y)
    if y > 1:
        return Fr(0)
    raise ValueError(y)

def dcorr(y):
    al = AL
    if Fr(1, 3) < y < Fr(1, 2):
        return pos(1 + 4 * al - 3 * y - pos(1 + al - 3 * y))
    return Fr(0)

def Tout(y):
    return R0(y) - dcorr(y) - 2 * LAM * fl(1 / y) + sum(pos(2 * LAM - j * y) for j in range(1, 6))

def I_out():
    al, lam = AL, LAM
    cands = {Fr(1, 3), Fr(1, 2), Fr(1), 2 * lam, (1 - al) / 2, (1 + al) / 3, 1 - al, (1 + al) / 2,
             (1 + 4 * al) / 2, (1 + 4 * al) / 3} | {2 * lam / j for j in range(1, 6)}
    pts = sorted(p for p in cands if Fr(1, 3) <= p <= 2 * lam)
    tot = Fr(0)
    rows = []
    for l, r in zip(pts[:-1], pts[1:]):
        y1, y2, ym = (2 * l + r) / 3, (l + 2 * r) / 3, (l + r) / 2
        c = (Tout(y2) - Tout(y1)) / (y2 - y1); b = Tout(y1) - c * y1
        assert Tout(ym) == b + c * ym and Tout(l + (r - l) / 101) == b + c * (l + (r - l) / 101)
        tot += b * (r - l) + c * (r * r - l * l) / 2
        rows.append((l, r, b, c))
    # integral of d over (1/3,1/2) should be 9/640 per the paper
    dd = Fr(0)
    for l, r in zip(pts[:-1], pts[1:]):
        if r <= Fr(1, 2):
            y1, y2 = (2 * l + r) / 3, (l + 2 * r) / 3
            c = (dcorr(y2) - dcorr(y1)) / (y2 - y1); b = dcorr(y1) - c * y1
            dd += b * (r - l) + c * (r * r - l * l) / 2
    return tot, rows, dd


if __name__ == "__main__":
    Iout, rows, dd = I_out()
    print(f"I_out = {Iout} = {float(Iout):.12f};  equals 127751/96000: {Iout == Fr(127751, 96000)};  int d = {dd} (paper 9/640)")
    print("  affine pieces (l, r, b, c):")
    for l, r, b, c in rows:
        print(f"    [{l}, {r}]  b = {b}, c = {c}")

    I320 = Fr(322437603634266857629, 7535670527041937280000)
    p, q = open('./data/15-redteam-mu5/int_20_4000.txt').read().split()
    I20_4000 = fmpq(int(p), int(q))
    lam, al = A(LAM), A(AL)
    X = 4000
    # 09's tail bound: int_X^oo R/x^3 <= -lam/X + (3lam/(2al) + 13/16)/X^2  (needs Q <= 13/8 on [X,oo), G(X) = 0)
    assert frac(Fr(X)) == 0 and frac(AL * X) == 0
    tail09 = -lam / X + (3 * lam / (2 * al) + A(Fr(13, 16))) / X ** 2
    # my exact Q range (exact_R.py (3c)): Q <= 2790/1849 on [20,oo)
    tail_mine = -lam / X + (3 * lam / (2 * al) + A(Fr(2790, 1849)) / 2) / X ** 2
    Astar09 = A(Iout) + A(I320) + arb(I20_4000) + tail09
    print(f"\nint_3^20 = {float(I320):.13f};  int_20^4000 = {arb(I20_4000).str(15)}")
    print(f"tail bound at 4000 (09 formula) = {tail09.str(12)};  with exact sup Q = {tail_mine.str(12)}")
    print(f"A* (09 route) = {Astar09.str(15)}")
    for M in (200, 100000):
        AM = Astar09 + 7 * lam / M - (A(Fr(2923, 240)) - A(Fr(1, 4))) / M ** 2 + A(32) / M ** 3
        print(f"A_{M} <= {AM.str(15)}   (upper {AM.upper().str(12)})")
    # Fauzan's own A* for comparison
    AstarF = Fr(9928298118277006344769, 7535670527041937280000)
    print(f"Fauzan A* = {float(AstarF):.12f}; check A* = I_out + (5.18) - 2689/48000: {AstarF == Fr(127751,96000) + I320 - Fr(2689,48000)}")

    # ------- high-precision exact tail via periodicity (cross-check, not part of the bound) -------
    # For x = y + 40k, y in [Y0, Y0+40): R(x) = Q(y) + x F(y)  =>  sum_k R(x)/x^3 = Q(y) zeta(3,y/40)/40^3 + F(y) zeta(2,y/40)/40^2
    import mpmath as mp
    mp.mp.dps = 30
    def exact_tail(Y0):
        tot = mp.mpf(0)
        for l, r, a, b in pieces(Y0, Y0 + 40, check=False):
            x1, x2 = (2 * l + r) / 3, (l + 2 * r) / 3
            pF = (F(x2) - F(x1)) / (x2 - x1); cF = F(x1) - pF * x1
            # on this piece Q(y) = R(y) - y F(y) = (a - cF) y + b - pF y^2  (affine R, affine F)
            def g(y):
                Qy = (mp.mpf(a.numerator) / a.denominator - (mp.mpf(cF.numerator) / cF.denominator)) * y \
                    + mp.mpf(b.numerator) / b.denominator - mp.mpf(pF.numerator) / pF.denominator * y * y
                Fy = mp.mpf(pF.numerator) / pF.denominator * y + mp.mpf(cF.numerator) / cF.denominator
                return Qy * mp.zeta(3, y / 40) / 40 ** 3 + Fy * mp.zeta(2, y / 40) / 40 ** 2
            tot += mp.quad(g, [mp.mpf(l.numerator) / l.denominator, mp.mpf(r.numerator) / r.denominator])
        return tot
    T20 = exact_tail(20)
    T4000 = exact_tail(4000)
    T1e5 = exact_tail(100000)
    print(f"\n[cross-check, mpmath] int_20^oo R/x^3 = {mp.nstr(T20, 15)};  int_4000^oo = {mp.nstr(T4000, 15)};  int_1e5^oo = {mp.nstr(T1e5, 15)}")
    print(f"   consistency: int_20^4000 (exact) + int_4000^oo - int_20^oo = {mp.nstr(mp.mpf(float(arb(I20_4000).mid())) + T4000 - T20, 5)}")
    Ainf = float(Iout) + float(I320) + T20
    A1e5_true = Ainf + 6 * 37 / 40 / 100000 - T1e5
    print(f"   true A_oo = {mp.nstr(Ainf, 12)};  true A_100000 = I_out + 6lam/M + int_3^M = {mp.nstr(A1e5_true, 12)}")
