"""Fauzan-type positive functional for Dirichlet beta(k), k even.

Alternating Abel-Plana  sum_{n>=0} (-1)^n f(n) = f(0)/2 + i int_0^oo (f(iy)-f(-iy))/(2 sinh(pi y)) dy, applied to
f(z) = (z+a)^{-k} and (k-1) integrations by parts (boundary terms vanish) gives, with
   w(y) = (2 pi^{k-1}/(k-1)!) y^k sum_{m>=0} (2m+1)^{k-1} e^{-(2m+1) pi y}  > 0,
   int_0^oo w(y)/(y^2+a^2) dy = a^{k-1} sum_{n>=0} (-1)^n (n+a)^{-k} - 1/(2a).
At a = j+1/2:  sum_{n>=0}(-1)^n (n+j+1/2)^{-k} = (-1)^j 2^k (beta(k) - B_j),  B_j = sum_{m<j} (-1)^m/(2m+1)^k.
Moments:  int y^{2e} w = (-1)^e (4^{e+1}-1) B_{2e+2} (2e+k)! / ((k-1)! (2e+2)!)   [Genocchi: integral at odd p]
So   mu_X(t^e) as above,  mu_X(1/(t+(j+1/2)^2)) = (-1)^j 2^k (j+1/2)^{k-1} (X - B_j) - 1/(2j+1).
This script checks both formulas by quadrature."""
import mpmath as mp
from math import factorial

mp.mp.dps = 25


def w(k, y):
    q = mp.e ** (-mp.pi * y)
    # sum_{m>=0} (2m+1)^{k-1} q^{2m+1}
    s = mp.polylog(1 - k, q) - 2 ** (k - 1) * mp.polylog(1 - k, q * q)   # odd part of sum n^{k-1} q^n
    return 2 * mp.pi ** (k - 1) / mp.factorial(k - 1) * y ** k * s


def beta(k):
    return mp.nsum(lambda m: (-1) ** m / (2 * m + 1) ** k, [0, mp.inf])


for k in (2, 4, 6):
    bk = beta(k)
    for e in (0, 1, 2):
        I = mp.quad(lambda y: y ** (2 * e) * w(k, y), [0, 1, 4, mp.inf])
        pred = (-1) ** e * (4 ** (e + 1) - 1) * mp.bernoulli(2 * e + 2) * mp.factorial(2 * e + k) / (mp.factorial(k - 1) * mp.factorial(2 * e + 2))
        print(f"k={k} mu(t^{e}): quad {mp.nstr(I, 20)}  formula {mp.nstr(pred, 20)}")
    for j in (0, 1, 2, 5):
        a = j + mp.mpf(1) / 2
        I = mp.quad(lambda y: w(k, y) / (y * y + a * a), [0, 1, 4, mp.inf])
        Bj = sum(mp.mpf((-1) ** m) / (2 * m + 1) ** k for m in range(j))
        pred = (-1) ** j * 2 ** k * a ** (k - 1) * (bk - Bj) - 1 / (2 * a)
        print(f"k={k} j={j}: quad {mp.nstr(I, 20)}  formula {mp.nstr(pred, 20)}")
