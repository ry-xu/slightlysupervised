"""Fauzan's arithmetic bookkeeping (Sections 4-5, App. B), generalized to (k, r, alpha, lam).

Construction: W = D_N^r / D_K (zeros of net order r-1 at -l^2, l<=N=alpha K; simple poles at -j^2, N<j<=K),
h = lam K rows (Fauzan: lam = 1-alpha), zeta(k), k odd.  Normalizer S_K = (K!)^{2h} 4^{h-1}/((N!)^{2rh} prod (2i)!^2).

log m_{K,M} / K^2 -> A_M,   A_inf = I_out + int_3^oo R(x)/x^3 dx   (M -> oo)
Inner range (x = K/p >= 3), Fauzan (4.2)-(4.8),(5.4)-(5.6) with 4 -> k-1, 6 -> r, 3 -> r/2:
   per-class valuation  L_a (L_a + 2 b_a - ell_K(a) - k),  b_a = (r/2) ell_N(a),  L_a = T - b_a + eps_a
   Gamma_k(x) = int_0^{1/2} (T-b)(T+b-ell-k) dz + s(2T-q-k) + (s-n_+)_+ ,  T = floor(2Hx), H = lam + r alpha/2
   Nf(x) = 2 lam x floor(x) - 2 r lam x floor(alpha x) - 2 J(lam x)
   R(x) = -Gamma_k(x) - Nf(x)
Outer range (y = p/K in (1/3,1)), Fauzan (4.12)-(4.14),(5.8)-(5.10):
   unremoved class with ell members: cost (k+2)(ell-2), 2 zero weights
   removed class (a<=N): weights min(0, i + r/2 - (ell+k-1)/2) for i<ell-2, plus one zero weight
   zero class: cost 1 (one pole p) or 4 (p,2p) -> O(1)
   rank correction min(r_p, Z), r_p/K = (2 lam + r alpha - 1 - (k-1) y/2)_+
   T(y) = costs/K + min(r_p,Z)/K - 2 lam floor(1/y) + sum_j (2 lam - j y)_+
Valid (as derived) for alpha < 1/6, ell <= 6 < k+5.
"""
import numpy as np
from math import floor
from fractions import Fraction as Fr


def ell(x, z):
    return floor(x - z) + floor(x + z) + 1


def d0(u):
    return min(u, 1 - u)


def frac(x):
    return x - floor(x)


def Gamma(x, k, r, al, lam):
    H = lam + r * al / 2
    T = floor(2 * H * x)
    s = H * x - T / 2
    q = floor(2 * x)
    npl = (2 * x - q) / 2
    bps = sorted(set([0.0, 0.5, d0(frac(x)), d0(frac(al * x))]))
    I = 0.0
    for a, b in zip(bps, bps[1:]):
        if b <= a:
            continue
        z = (a + b) / 2
        l = ell(x, z)
        bb = (r / 2) * ell(al * x, z)
        I += (T - bb) * (T + bb - l - k) * (b - a)
    return I + s * (2 * T - q - k) + max(s - npl, 0.0)


def J(u):
    m = floor(2 * u)
    return m * u - m * (m + 1) / 4


def Nf(x, r, al, lam):
    return 2 * lam * x * floor(x) - 2 * r * lam * x * floor(al * x) - 2 * J(lam * x)


def Rin(x, k, r, al, lam):
    return -Gamma(x, k, r, al, lam) - Nf(x, r, al, lam)


def inner_breakpoints(x0, x1, r, al, lam):
    H = lam + r * al / 2
    cs = [1.0, 2.0, al, 2 * al, 2 * lam, 2 * H, 2 * abs(H - 1), 2 * (1 - al), 2 * (1 + al), 4 * al]
    pts = {x0, x1}
    for c in cs:
        if c <= 0:
            continue
        k0 = int(np.ceil(x0 * c))
        k1 = int(np.floor(x1 * c))
        for kk in range(k0, k1 + 1):
            v = kk / c
            if x0 < v < x1:
                pts.add(v)
    return sorted(pts)


_gl3 = np.polynomial.legendre.leggauss(4)


def int_inner(x0, x1, k, r, al, lam, f=None):
    """int_{x0}^{x1} R(x)/x^3 dx exactly-ish: R is piecewise quadratic between breakpoints."""
    if f is None:
        f = lambda x: Rin(x, k, r, al, lam) / x ** 3
    pts = inner_breakpoints(x0, x1, r, al, lam)
    nodes, wts = _gl3
    tot = 0.0
    for a, b in zip(pts, pts[1:]):
        if b - a < 1e-15:
            continue
        # R piecewise polynomial (deg<=2) on (a,b); 1/x^3 smooth: subdivide for accuracy
        nsub = max(1, int((b - a) / max(a, 1.0) * 40) + 1)
        for i in range(nsub):
            aa = a + (b - a) * i / nsub
            bb = a + (b - a) * (i + 1) / nsub
            xm = (aa + bb) / 2
            hw = (bb - aa) / 2
            for t, w in zip(nodes, wts):
                tot += w * hw * f(xm + hw * t)
    return tot


def inner_tail(X, k, r, al, lam):
    """int_X^oo R/x^3 using R = x F + Q with mean F = lam (k - r) (oscillating parts O(1/X^2))."""
    Fbar = lam * (k - r)
    return Fbar / X


# ---------------- outer range ----------------
def removed_cost(l, k, r):
    c = 0.0
    zeros = 1
    for i in range(0, l - 2):
        w = i + r / 2 - (l + k - 1) / 2
        if w < 0:
            c += -2 * w
        else:
            zeros += 1
    return c, zeros


def unremoved_cost(l, k):
    return (k + 2) * (l - 2), (2 if l >= 2 else l)


def outer_T(y, k, r, al, lam, return_parts=False):
    """-L_p/K in the limit, y = p/K in (1/3, 2 lam)."""
    if y >= 1:
        val = sum(max(2 * lam - j * y, 0.0) for j in range(1, 8))
        return (val, 0, 0, 0) if return_parts else val
    if not (1 / 3 < y < 1):
        raise ValueError
    # class structure in z = a/p in (0,1/2)
    if y > 0.5:
        base = 2
        c1 = (1 - y) / y            # a <= K - p
        c2 = 2 - 1 / y              # a >= 2p - K
    else:
        base = 4
        c1 = 1 / y - 2              # a <= K - 2p
        c2 = 3 - 1 / y              # a >= 3p - K
    zr = al / y                     # removed: a <= N
    bps = sorted(set([0.0, 0.5] + [v for v in (c1, c2, zr) if 0 < v < 0.5]))
    cost = 0.0
    Z = 0.0
    for a, b in zip(bps, bps[1:]):
        z = (a + b) / 2
        l = base + (1 if z <= c1 else 0) + (1 if z >= c2 else 0)
        if z <= zr:
            c, zz = removed_cost(l, k, r)
        else:
            c, zz = unremoved_cost(l, k)
        cost += c * (b - a)
        Z += zz * (b - a)
    cost *= y
    Z *= y
    rp = max(0.0, 2 * lam + r * al - 1 - (k - 1) * y / 2)
    rank = min(rp, Z)
    scal = -2 * lam * floor(1 / y) + sum(max(2 * lam - j * y, 0.0) for j in range(1, 8))
    if return_parts:
        return cost + rank + scal, cost, rank, scal
    return cost + rank + scal


def int_outer(k, r, al, lam, n=40000):
    """int_{1/3}^{2 lam} T(y) dy by composite Gauss-Legendre (T piecewise affine)."""
    ys = np.linspace(1 / 3, 2 * lam, n + 1)
    nodes, wts = _gl3
    tot = 0.0
    parts = np.zeros(3)
    for a, b in zip(ys[:-1], ys[1:]):
        xm, hw = (a + b) / 2, (b - a) / 2
        for t, w in zip(nodes, wts):
            v, c, rk, sc = outer_T(xm + hw * t, k, r, al, lam, return_parts=True)
            tot += w * hw * v
            parts += w * hw * np.array([c, rk, sc])
    return tot, parts


def A_inf(k, r, al, lam, X=600.0, nout=20000, verbose=False):
    Iout, parts = int_outer(k, r, al, lam, n=nout)
    Iin = int_inner(3.0, X, k, r, al, lam)
    tail = inner_tail(X, k, r, al, lam)
    A = Iout + Iin + tail
    if verbose:
        print(f"k={k} r={r} al={al:.5f} lam={lam:.5f}: I_out={Iout:.6f} [classes {parts[0]:.5f}, rank {parts[1]:.5f}, "
              f"S_K {parts[2]:.5f}]  int_3^X R/x^3={Iin:.6f}  tail={tail:.6f}  A_inf={A:.6f}")
    return A, dict(Iout=Iout, Iin=Iin, tail=tail, parts=parts)


if __name__ == "__main__":
    k, r, al, lam = 5, 6, 3 / 40, 37 / 40
    # Fauzan checks
    I320 = int_inner(3.0, 20.0, k, r, al, lam)
    print("int_3^20 R/x^3 =", I320, "  paper (5.18) =", 322437603634266857629 / 7535670527041937280000)
    Iout, parts = int_outer(k, r, al, lam)
    print("I_out =", Iout, "  paper (5.10) =", 127751 / 96000, parts)
    I20 = int_inner(20.0, 2000.0, k, r, al, lam) + inner_tail(2000.0, k, r, al, lam)
    print("int_20^oo R/x^3 =", I20, "  paper bound (5.16) <=", -2689 / 48000)
    print("A* paper =", 9928298118277006344769 / 7535670527041937280000)
    A, info = A_inf(k, r, al, lam, verbose=True)
