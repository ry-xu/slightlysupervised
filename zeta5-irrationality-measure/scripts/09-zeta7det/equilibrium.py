"""Analytic (log-energy) exponent for the Fauzan-type Hankel determinant.

  Delta_K(zeta) = (1/h!) int prod_{i<j}(y_i^2-y_j^2)^2 prod W(y_i^2) w(y_i) dy_i,   y = K sqrt(t)
  log F_K(zeta)/K^2 -> U = E* + C*,
  E* = max_{sigma>=0, |sigma|=lam} [ I(sigma) - int V dsigma ],   I = intint log|s-t|,
  V(t) = 2 pi sqrt(t) - int_0^oo e(u) log(t+u^2) du    (e(u) = net exponent density of W per unit u = l/K)
  Fauzan profile: e(u) = r-1 on [0,alpha], -1 on (alpha,1]  =>  V = 2pi sqrt t + int_0^1 log - r int_0^alpha log.
  C* = -2lam + 2 r alpha lam (1-log alpha) + 3lam^2 - 2 lam^2 log(2 lam)    (from S_K, Fauzan (6.3) for r=6)

Equilibrium measure on one interval [a,b] with soft edges (Chebyshev method):
  2 pv int dsigma(s)/(t-s) = V'(t),  t = c + d x,  g(x) = V'(c+dx) = sum g_k T_k(x)
  conditions: g_0 = 0, mass d g_1 / 4 = lam;  density psi(c+dy) = (1/(2pi)) sqrt(1-y^2) sum_{k>=1} g_k U_{k-1}(y) (per unit y)
"""
import numpy as np
from numpy.polynomial import chebyshev as C
import mpmath as mp
from scipy.optimize import brentq, fsolve
from scipy.fft import dct


def make_profile(pieces):
    """pieces: list of (u0, u1, e) meaning exponent density e on [u0,u1] (in units of K)."""
    return pieces


def fauzan_profile(r, alpha, top=1.0):
    return [(0.0, alpha, r - 1.0), (alpha, top, -1.0)]


def _I(a, t):
    """int_0^a log(t+u^2) du"""
    st = np.sqrt(t)
    return a * np.log(t + a * a) - 2 * a + 2 * st * np.arctan(a / st)


def V(t, prof):
    t = np.asarray(t, dtype=float)
    out = 2 * np.pi * np.sqrt(t)
    for (u0, u1, e) in prof:
        out = out - e * (_I(u1, t) - _I(u0, t))
    return out


def dV(t, prof):
    t = np.asarray(t, dtype=float)
    st = np.sqrt(t)
    out = np.pi / st
    for (u0, u1, e) in prof:
        out = out - e * (np.arctan(u1 / st) - np.arctan(u0 / st)) / st
    return out


def cheb_coeffs(a, b, prof, n):
    c, d = (a + b) / 2, (b - a) / 2
    # Chebyshev points of the first kind
    j = np.arange(n)
    x = np.cos(np.pi * (j + 0.5) / n)
    g = dV(c + d * x, prof)
    co = dct(g, type=2) / n
    co[0] /= 2
    return co


def conditions(a, b, prof, lam, n):
    co = cheb_coeffs(a, b, prof, n)
    d = (b - a) / 2
    return co[0], d * co[1] / 4 - lam


def solve_support(prof, lam, n=4096, guess=(1e-3, 0.8)):
    def F(v):
        a, b = np.exp(v[0]), np.exp(v[1])
        if a >= b:
            return [1e3, 1e3]
        g0, mres = conditions(a, b, prof, lam, n)
        return [g0, mres]
    sol, info, ier, msg = fsolve(F, [np.log(guess[0]), np.log(guess[1])], full_output=True, xtol=1e-14)
    a, b = np.exp(sol)
    return a, b, ier, F(sol)


def density(a, b, prof, n, y):
    co = (b - a) / 2 * cheb_coeffs(a, b, prof, n)
    # psi(c+dy) per unit y = (1/2pi) sqrt(1-y^2) sum_{k>=1} g_k U_{k-1}(y); per unit t divide by d
    # U_{k-1}(y) = sin(k th)/sin th, y = cos th
    th = np.arccos(np.clip(y, -1, 1))
    k = np.arange(1, len(co))
    s = np.sin(np.outer(th, k)) @ co[1:]
    return s / (2 * np.pi)  # sqrt(1-y^2)*U_{k-1} = sin(k th); per unit y


def energy(a, b, prof, lam, n=4096):
    """E* = I(sigma) - int V dsigma for the soft-edge equilibrium measure on [a,b].
    sigma(dy) = psi(y) dy, psi(y) = (1/2pi) sum_k g_k sqrt(1-y^2) U_{k-1}(y), g = Chebyshev coeffs of d V'(c+dx).
    Potential: int_{-1}^1 log|x-y| sqrt(1-y^2) U_{k-1}(y) dy = (pi/2)(x^2-1/2-log2) (k=1),
                                                             (pi/2)[T_{k+1}/(k+1) - T_{k-1}/(k-1)] (k>=2).
    Quadrature: y = cos(th), th_j = pi j/n (trapezoid in theta, spectrally accurate)."""
    from scipy.fft import dst
    co = cheb_coeffs(a, b, prof, n)
    c, d = (a + b) / 2, (b - a) / 2
    g = d * co
    # Chebyshev-T coefficients of Q(x) = sum_{k>=2} g_k (T_{k+1}/(k+1) - T_{k-1}/(k-1))
    q = np.zeros(n + 2)
    for k in range(2, n):
        q[k + 1] += g[k] / (k + 1)
        q[k - 1] -= g[k] / (k - 1)
    def pot(x):
        return (g[1] * (np.pi / 2) * (x * x - 0.5 - np.log(2)) + (np.pi / 2) * C.chebval(x, q)) / (2 * np.pi)
    j = np.arange(1, n)
    th = np.pi * j / n
    y = np.cos(th)
    S = dst(g[1:], type=1) / 2          # S_j = sum_{k=1}^{n-1} g_k sin(pi j k / n)
    w = (np.pi / n) * S * np.sin(th) / (2 * np.pi)
    mass = w.sum()
    P = pot(y) + mass * np.log(d)       # U^sigma(c+dy)
    Vy = V(c + d * y, prof)
    I = P @ w
    ell = 2 * P - Vy
    sel = w > 1e-14 * w.max()
    return I - Vy @ w, mass, ell[sel].min(), ell[sel].max(), w


def Cstar(r, alpha, lam):
    return -2 * lam + 2 * r * alpha * lam * (1 - np.log(alpha)) + 3 * lam ** 2 - 2 * lam ** 2 * np.log(2 * lam)


def U_fauzan(r, alpha, lam, n=8192, guess=(1e-4, 0.75), verbose=False):
    prof = fauzan_profile(r, alpha)
    a, b, ier, resid = solve_support(prof, lam, n=n, guess=guess)
    E, mass, lo, hi, w = energy(a, b, prof, lam, n=n)
    if verbose:
        print(f"support [{a:.6g},{b:.6g}] ier={ier} resid={resid} mass={mass:.8f} ell in [{lo:.8f},{hi:.8f}] minw={w.min():.3e}")
    return E + Cstar(r, alpha, lam), dict(a=a, b=b, E=E, Cstar=Cstar(r, alpha, lam), mass=mass, ell=(lo, hi), minw=w.min())


if __name__ == "__main__":
    for n in (2048, 8192, 32768):
        U, info = U_fauzan(6, 3 / 40, 37 / 40, n=n, verbose=True)
        print(n, "U =", U, " (Fauzan bound -2733991/2000000 =", -2733991 / 2000000, ")")
