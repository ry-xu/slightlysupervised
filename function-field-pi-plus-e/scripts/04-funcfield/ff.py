"""
Core library for 04-funcfield: Laurent series in x = 1/theta over F_p.

Conventions
-----------
A = F_q[theta], k = F_q(theta), k_inf = F_q((x)), x = 1/theta, q = p^r.
|f| = q^{deg f}; for a series f = sum_{n>=v} c_n x^n with c_v != 0, deg f = -v.

All the constants we study have coefficients in F_p even when q = p^r > p:
  D_i = prod_{j<i} (theta^{q^i} - theta^{q^j})   (monic, in F_p[theta])
  e   = e_C(1) = sum_i 1/D_i,      1/D_i = x^{i q^i} prod_{j<i} (1 - x^{q^i - q^j})^{-1}
  u   = prod_{m>=1} (1 - theta^{1-q^m})^{-1} = prod_{m>=1} (1 - x^{q^m - 1})^{-1}
  pi~ = theta * (-theta)^{1/(q-1)} * u     (q = 2: pi~ = theta^2 * u)
so arithmetic mod p suffices (linear independence over F_p(theta) implies it over
F_q(theta) because F_q((x)) = F_q (x)_{F_p} F_p((x))).

A truncated series is a numpy int64 array c[0..N-1] (coefficients of x^0..x^{N-1}).
"""
import numpy as np
import flint


# ---------------------------------------------------------------- basic series ops
def div_one_minus_xa(c, a, p):
    """In place: c <- c / (1 - x^a)  (mod p, truncated)."""
    N = len(c)
    s = a
    while s < N:
        e = min(s + a, N)
        c[s:e] = (c[s:e] + c[s - a:e - a]) % p
        s = e
    return c


def mul_one_minus_xa(c, a, p):
    """In place: c <- c * (1 - x^a)."""
    if a < len(c):
        c[a:] = (c[a:] - c[:-a].copy()) % p
    return c


def to_poly(c, p):
    return flint.nmod_poly([int(v) for v in c], p)


def from_poly(P, N):
    out = np.zeros(N, dtype=np.int64)
    cs = P.coeffs()
    m = min(len(cs), N)
    if m:
        out[:m] = np.array([int(v) for v in cs[:m]], dtype=np.int64)
    return out


def mul(c1, c2, p):
    N = min(len(c1), len(c2))
    P = to_poly(c1[:N], p).mul_low(to_poly(c2[:N], p), N)
    return from_poly(P, N)


def inv(c, p):
    N = len(c)
    return from_poly(to_poly(c, p).inverse_series_trunc(N), N)


def val(c):
    """x-adic valuation (index of first nonzero coeff), or None if zero to precision."""
    nz = np.flatnonzero(c)
    return int(nz[0]) if len(nz) else None


# ---------------------------------------------------------------- the constants
def u_series(q, N, p):
    c = np.zeros(N, dtype=np.int64); c[0] = 1
    m = 1
    while q ** m - 1 < N:
        div_one_minus_xa(c, q ** m - 1, p)
        m += 1
    return c


def inv_D_series(q, i, N, p):
    """1/D_i as a series (zero array if its valuation i*q^i >= N)."""
    c = np.zeros(N, dtype=np.int64)
    v = i * q ** i
    if v >= N:
        return c
    c[v] = 1
    for j in range(i):
        div_one_minus_xa(c, q ** i - q ** j, p)
    return c


def e_series(q, N, p, upto=None):
    """e_C(1) = sum_{i} 1/D_i ; if upto is given, the partial sum S_upto."""
    c = np.zeros(N, dtype=np.int64)
    i = 0
    while i * q ** i < N and (upto is None or i <= upto):
        c = (c + inv_D_series(q, i, N, p)) % p
        i += 1
    return c


def U_M_series(q, M, N, p):
    """Partial product U_M = prod_{m<=M} (1 - x^{q^m-1})^{-1}."""
    c = np.zeros(N, dtype=np.int64); c[0] = 1
    for m in range(1, M + 1):
        if q ** m - 1 < N:
            div_one_minus_xa(c, q ** m - 1, p)
    return c


# ---------------------------------------------------------------- polynomials in theta
def theta_poly_from_xpoly(c, D, p):
    """Given a polynomial c(x) of degree <= D (array), return the theta-polynomial theta^D c(1/theta)."""
    coeffs = [0] * (D + 1)
    for n, v in enumerate(c[: D + 1]):
        if v % p:
            coeffs[D - n] = int(v % p)
    return flint.nmod_poly(coeffs, p)


def D_poly(q, i, p):
    """D_i as an nmod_poly in theta."""
    T = flint.nmod_poly([0, 1], p)
    P = flint.nmod_poly([1], p)
    for j in range(i):
        P *= T ** (q ** i) - T ** (q ** j)
    return P


def E_poly(q, M, p):
    """E_M = prod_{m<=M} (theta^{q^m - 1} - 1), the theta-denominator of U_M."""
    T = flint.nmod_poly([0, 1], p)
    P = flint.nmod_poly([1], p)
    for m in range(1, M + 1):
        P *= T ** (q ** m - 1) - 1
    return P


# ---------------------------------------------------------------- continued fractions
def cf_degrees(c, p, maxdeg=None):
    """Continued fraction of alpha = sum c_n x^n (alpha in F_p[[x]]) via Euclid on
    (A, theta^{N-1}), A = theta^{N-1} * alpha_trunc.  Returns list of (deg a_n, deg Q_n)
    where Q_n is the denominator of the n-th convergent; only reliable while
    2*deg Q_n + deg a_{n+1} < N-1 (we stop there)."""
    N = len(c)
    Npow = N - 1
    A = theta_poly_from_xpoly(c, Npow, p)
    B = flint.nmod_poly([0] * Npow + [1], p)
    out = []
    degQ_prev, degQ = None, 0     # deg Q_{-1} = -inf, deg Q_0 = 0
    a0, r = divmod(A, B)
    A, B = B, r
    out.append((a0.degree(), 0))
    while not B.is_zero():
        a, r = divmod(A, B)
        da = a.degree()
        newdegQ = degQ + da
        # reliability: the truncation error is |x^N| = q^{-N}; convergent n is reliable
        # while |alpha - P/Q| = q^{-2 deg Q_n - deg a_{n+1}} >> q^{-(N-1)}
        if 2 * degQ + da >= Npow - 2:
            break
        out.append((da, newdegQ))
        degQ = newdegQ
        if maxdeg is not None and degQ > maxdeg:
            break
        A, B = B, r
    return out
