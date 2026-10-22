"""18-redteam-mu96: third, independent exact evaluation of int_20^4000 R(x)/x^3 dx via Fauzan's closed form
(5.12)-(5.13) + (B.1) (09 and 15 used the defining formulas (5.4)-(5.6)).
  R = xF + Q,  F = 4lam + 2lam f - 12lam g,
  Q = [tau(tau-sigma) - (tau-sigma)_+ + eta(1-eta)]/2 + 9 d0(g)(1-2d0(g)) - 3 e(f)e(g)(min(d0 f, d0 g) - 2 d0(f) d0(g)),
  f={x}, g={al x}, tau={2Hx}, sigma={2x}, eta={2 lam x}, d0(u)=min(u,1-u), e(u)=1 (u<=1/2), -1 (u>1/2).
R is affine on each piece between consecutive points of the union of (1/2)Z,(20/3)Z,(20/37)Z,(10/23)Z,(10/3)Z,(20/43)Z
(x^2 coefficients: 2lam(1-6al) + 2H(H-1) - 2lam^2 - 18al^2 + 6al = 0). int (b+cx)/x^3 = [-b/(2x^2) - c/x].
"""
import sys, time
from flint import fmpq, fmpz

lam, al, H = fmpq(37, 40), fmpq(3, 40), fmpq(23, 20)


def frac(x):
    return x - x.floor()


def d0(u):
    return min(u, 1 - u)


def e(u):
    return 1 if u <= fmpq(1, 2) else -1


def R(x):
    f, g = frac(x), frac(al * x)
    tau, sig, eta = frac(2 * H * x), frac(2 * x), frac(2 * lam * x)
    F = 4 * lam + 2 * lam * f - 12 * lam * g
    ts = tau - sig
    Q = (tau * (tau - sig) - (ts if ts > 0 else fmpq(0)) + eta * (1 - eta)) / 2 + 9 * d0(g) * (1 - 2 * d0(g)) \
        - 3 * e(f) * e(g) * (min(d0(f), d0(g)) - 2 * d0(f) * d0(g))
    return x * F + Q


def breakpoints(a, b):
    pts = set([a, b])
    for step in [fmpq(1, 2), fmpq(20, 3), fmpq(20, 37), fmpq(10, 23), fmpq(10, 3), fmpq(20, 43)]:
        k = (a / step).ceil()
        while True:
            x = step * k
            if x > b:
                break
            if x >= a:
                pts.add(x)
            k += 1
    return sorted(pts)


def integral(a, b, check_affine=True):
    pts = breakpoints(a, b)
    tot = fmpq(0)
    nbad = 0
    for l, r in zip(pts[:-1], pts[1:]):
        x1 = l + (r - l) / 3
        x2 = l + 2 * (r - l) / 3
        R1, R2 = R(x1), R(x2)
        c = (R2 - R1) / (x2 - x1)
        bb = R1 - c * x1
        if check_affine:
            xm = (l + r) / 2
            if R(xm) != bb + c * xm:
                nbad += 1
        tot += (-bb / (2 * r * r) - c / r) - (-bb / (2 * l * l) - c / l)
    return tot, len(pts) - 1, nbad


if __name__ == "__main__":
    t0 = time.time()
    v, n, nbad = integral(fmpq(20), fmpq(4000))
    print("int_20^4000 R/x^3 (closed form (5.12)-(5.13)+(B.1)) = %.16f ; pieces %d ; non-affine pieces %d [%.1fs]" % (
        float(v), n, nbad, time.time() - t0))
    num, den = open("./data/15-redteam-mu5/int_20_4000.txt").read().split()
    v15 = fmpq(int(num), int(den))
    print("15's exact value                                  = %.16f ; exactly equal: %s" % (float(v15), v == v15))
    # also (5.18): int_3^20
    v318, n2, nb2 = integral(fmpq(3), fmpq(20))
    print("int_3^20 R/x^3 via closed form = %s ; (5.18) = 322437603634266857629/7535670527041937280000 ; equal: %s ; non-affine %d" % (
        v318, v318 == fmpq(322437603634266857629, 7535670527041937280000), nb2))
