"""Exact Hankel-determinant engine for the Fauzan-type construction at zeta(k), k odd >= 3.

Functional (Fauzan sec. 2, Calegari blog for general k):
   mu_X(t^e)          = (-1)^e B_{2e+2} (2e+k)! / ((k-1)! (2e+2)!)
   mu_X(1/(t+j^2))    = j^{k-1} (X - H_j^{(k)}) - 1/(k-1) + 1/(2j)
At X = zeta(k):  mu(R) = int_0^oo R(y^2) w_k(y) dy,  w_k(y) = 2(2pi)^{k-1} y^k/(k-1)! sum_l l^{k-1} e^{-2 pi l y} > 0.

Weight  W(t) = t^{e0} prod_{l>=1} (t+l^2)^{e_l},  e_l in {-1,0,1,2,...}  (poles simple, at -j^2 with e_j=-1).
Hankel  G(X) = [ mu_X( W(t) t^{i+j} ) ]_{0<=i,j<h},   Delta(X) = det G(X) in Q[X].
Fauzan:  e_l = r-1 (l<=N), -1 (N<l<=K), r=6, K=40n, N=3n, h=K-N.
"""
import time
from math import factorial
import flint
from flint import fmpz, fmpq, fmpz_poly, fmpq_poly, fmpq_mat, fmpz_mat, arb, ctx

_bern_cache = {}


def bern(n):
    if n not in _bern_cache:
        _bern_cache[n] = fmpq.bernoulli(n)
    return _bern_cache[n]


def mu_mono(k, L):
    """mu(t^e) for e = 0..L"""
    out = []
    fk1 = factorial(k - 1)
    for e in range(L + 1):
        v = bern(2 * e + 2) * fmpq(factorial(2 * e + k), fk1 * factorial(2 * e + 2))
        out.append(v if e % 2 == 0 else -v)
    return out


def harmonic(k, J):
    H = [fmpq(0)] * (J + 1)
    s = fmpq(0)
    for j in range(1, J + 1):
        s += fmpq(1, j ** k)
        H[j] = s
    return H


def fauzan_exps(K, N, r):
    """e_l for W = D_N^r / D_K"""
    e = {}
    for l in range(1, K + 1):
        e[l] = (r - 1) if l <= N else -1
    return e


def entries(k, e, h, e0=0):
    """Return lists a[m], b[m] (fmpq), m = 0..2h-2, with mu_X(W t^m) = a[m] + b[m] X."""
    poles = sorted(l for l, v in e.items() if v == -1)
    Lmax = max(e) if e else 0
    num = fmpz_poly([0] * e0 + [1])
    for l, v in e.items():
        if v > 0:
            num *= fmpz_poly([l * l, 1]) ** v
    den = fmpz_poly([1])
    for j in poles:
        den *= fmpz_poly([j * j, 1])
    Q0 = num // den  # polynomial part of W (exact division of the polynomial part)
    q0 = [fmpq(int(c)) for c in Q0.coeffs()] if Q0.degree() >= 0 else []
    # residues W_j
    Wj = {}
    for j in poles:
        x = -j * j
        nv = fmpz(int(num(x)))
        dv = fmpz(1)
        for i in poles:
            if i != j:
                dv *= (i * i - j * j)
        Wj[j] = fmpq(nv, 1) / dv
    H = harmonic(k, max(poles) if poles else 1)
    c1 = {j: Wj[j] * (j ** (k - 1)) for j in poles}
    c0 = {j: Wj[j] * (fmpq(1, 2 * j) - fmpq(1, k - 1) - (j ** (k - 1)) * H[j]) for j in poles}
    M = 2 * h - 1
    mom = mu_mono(k, M + len(q0) + 1)
    a, b, s = [], [], []
    cur1 = dict(c1)
    cur0 = dict(c0)
    curW = dict(Wj)
    for m in range(M):
        # residue sums
        bm = sum(cur1.values(), fmpq(0))
        am_res = sum(cur0.values(), fmpq(0))
        sm = sum(curW.values(), fmpq(0))
        # mu(Q_m) = sum_{e<m} s_{m-1-e} mu(t^e) + sum_e q0[e] mu(t^{e+m})
        muQ = fmpq(0)
        for ee in range(m):
            muQ += s[m - 1 - ee] * mom[ee]
        for ee, c in enumerate(q0):
            if c != 0:
                muQ += c * mom[ee + m]
        a.append(muQ + am_res)
        b.append(bm)
        s.append(sm)
        for j in poles:
            f = -j * j
            cur1[j] *= f
            cur0[j] *= f
            curW[j] *= f
    return a, b


def hankel_det(a, b, h):
    """Delta(X) = det(A + X B) as fmpq_poly (A, B Hankel from a, b)."""
    A = fmpq_mat(h, h, [a[i + j] for i in range(h) for j in range(h)])
    B = fmpq_mat(h, h, [b[i + j] for i in range(h) for j in range(h)])
    detB = B.det()
    if detB != 0:
        Mx = B.solve(A)
        cp = (-Mx).charpoly()  # det(x I + B^-1 A)
        return fmpq_poly([c for c in cp.coeffs()]) * detB
    # fallback: interpolation
    xs = list(range(h + 1))
    vals = [(A + B * x).det() for x in xs]
    return lagrange(xs, vals)


def lagrange(xs, ys):
    X = fmpq_poly([0, 1])
    P = fmpq_poly([0])
    for i, xi in enumerate(xs):
        L = fmpq_poly([1])
        d = fmpq(1)
        for j, xj in enumerate(xs):
            if j != i:
                L *= (X - xj)
                d *= (xi - xj)
        P += L * (ys[i] / d)
    return P


def content(P):
    cs = [c for c in P.coeffs() if c != 0]
    g = fmpz(0)
    L = fmpz(1)
    for c in cs:
        g = abs(c.p) if g == 0 else g.gcd(c.p)
        L = L * c.q // L.gcd(c.q)
    return fmpq(g, L)


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def vp_int(x, p):
    x = fmpz(x)
    if x == 0:
        return 10 ** 9
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def factor_content(c, B):
    """v_p(c) for primes p <= B; returns dict and leftover |num|,|den| after removal."""
    num, den = abs(c.p), c.q
    out = {}
    for p in primes_upto(B):
        v = 0
        while num % p == 0:
            num //= p
            v += 1
        while den % p == 0:
            den //= p
            v -= 1
        if v:
            out[p] = v
    return out, num, den


def log_fmpz(x):
    x = abs(fmpz(x))
    if x == 0:
        return float('-inf')
    nb = x.bit_length()
    if nb < 1000:
        import math
        return math.log(int(x))
    sh = nb - 64
    import math
    return math.log(int(x >> sh)) + sh * math.log(2)


def log_fmpq(c):
    return log_fmpz(c.p) - log_fmpz(c.q)


def zeta_arb(k, prec):
    old = ctx.prec
    ctx.prec = prec
    z = arb(k).zeta()
    ctx.prec = old
    return z


def eval_at_zeta(P, k):
    """log|P(zeta(k))| and sign for integer polynomial P (fmpz_poly), rigorous via arb."""
    import math
    cs = [fmpz(int(c)) if not isinstance(c, fmpz) else c for c in P.coeffs()]
    mx = max(log_fmpz(c) for c in cs if c != 0)
    prec = int(mx / math.log(2)) + 400
    while True:
        old = ctx.prec
        ctx.prec = prec
        z = arb(k).zeta()
        v = arb(0)
        for c in reversed(cs):
            v = v * z + arb(c)
        ok = not v.contains(0)
        if ok:
            lv = v.abs_lower().log()
            res = (float(v.abs_upper().log().mid()), 1 if v > 0 else -1, float(lv.mid()))
            ctx.prec = old
            return res
        ctx.prec = old
        prec *= 2
        if prec > 10 ** 7:
            raise RuntimeError("precision blowup")


def run(k, e, h, e0=0, B=None, verbose=False):
    t0 = time.time()
    a, b = entries(k, e, h, e0)
    t1 = time.time()
    D = hankel_det(a, b, h)
    t2 = time.time()
    c = content(D)
    P = D / c
    Pz = fmpz_poly([int(x) for x in P.coeffs()])
    lv, sg, _ = eval_at_zeta(Pz, k)
    t3 = time.time()
    res = dict(k=k, h=h, deg=D.degree(), log_cont=log_fmpq(c), logP=lv, sign=sg,
               logDelta=lv + log_fmpq(c), t_ent=t1 - t0, t_det=t2 - t1, t_eval=t3 - t2)
    if B is not None:
        fac, rn, rd = factor_content(c, B)
        res['vp'] = fac
        res['left'] = (int(rn) if rn.bit_length() < 64 else 'big', int(rd) if rd.bit_length() < 64 else 'big')
    if verbose:
        print({x: res[x] for x in res if x != 'vp'}, flush=True)
    return res, D, P


if __name__ == "__main__":
    import sys
    k = int(sys.argv[1]); K = int(sys.argv[2]); N = int(sys.argv[3]); r = int(sys.argv[4])
    h = int(sys.argv[5]) if len(sys.argv) > 5 else K - N
    res, D, P = run(k, fauzan_exps(K, N, r), h, B=4 * K, verbose=True)
