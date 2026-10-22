"""18-redteam-mu96: explicit Theorem-2 bound vs true log||B_K||; and 16's Lemma 4 closed form for |r_j| (Theorem 1 only)."""
import math
from fractions import Fraction
from math import factorial
import sys
sys.path.insert(0, './scripts/18-redteam-mu96')
from true_norm_lanczos import r_coeffs
true = {40: 292.5613171408397, 80: 591.0372993642756, 160: 1188.2728511251, 240: 1785.6413312423,
        320: 2383.0664002111, 480: 3577.9992183777, 640: 4772.9897927749}
rho2 = 466851 / 62500
tS = 0.000071243158827
x1 = -1.069911
c0 = math.log((2 * math.pi) ** 4 / 24) - 2
lenG = 2 * (0.4 + 1 - x1); maxt2 = x1 * x1 + 0.04; dG = -1 - x1
print("Theorem 2 explicit: log||B_K|| <= K rho2 + 7logK + 2log h + E(K) + log(len max|tau|^2/2pi) + 7/d_Gamma")
for K in sorted(true):
    h = 37 * K // 40
    E = 3 * math.pi / (2 * K * math.sqrt(tS)) - 3 * math.log(K) - 2 * math.log(tS) - c0
    b = K * rho2 + 7 * math.log(K) + 2 * math.log(h) + E + math.log(lenG * maxt2 / (2 * math.pi)) + 7 / dG
    print("  K=%4d  bound %.3f  true %.3f  bound-true %.3f  (non-exponential part %.2f)" % (K, b, true[K], b - true[K], b - K * rho2))
# Lemma 4 closed form check
for K in [40, 80]:
    N = 3 * K // 40
    r = r_coeffs(K)
    ok = True
    for j in range(N + 1, K + 1):
        cf = Fraction(2 * factorial(j + N) ** 6, factorial(j - N - 1) ** 6 * factorial(K - j) * factorial(K + j)) * Fraction(K) ** (2 * K - 12 * N - 6)
        if abs(r[j]) != cf:
            ok = False
        if (r[j] > 0) != ((j - 1) % 2 == 0):
            ok = False
    print("Lemma 4 closed form |r_j| and sign (-1)^(j-1) exact for all N<j<=K at K=%d: %s" % (K, ok))
