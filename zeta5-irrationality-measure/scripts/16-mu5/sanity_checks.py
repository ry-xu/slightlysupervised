"""16-mu5: numerical sanity checks of the finite-K inequalities in the proof of Theorem 1 (notes/16-mu5.md).

 (S1) closed form |r_j| = 2 K^(2K-12N-6) (j+N)!^6 / ((j-N-1)!^6 (K-j)! (K+j)!)  vs the defining product (C.3);
      and the bound log|r_j| <= K phi(j/K) + log 2 + 1/2 + 3 log(4 pi K)            (Lemma 4)
 (S2) weight bound  log W_K(t) >= -K V(t) + 3 log K + 2 log t + log((2pi)^4/24) - 2   (Lemma 3)
 (S3) final bound  log||B_K|| <= K (M1 - m0) + 10 log K + C_rho + 3 pi/(2 K sqrt(t_min))   vs the exact
      ||B_K|| = 1/dist(zeta(5), roots of Delta_K) measured by 09 (root_distance.py) at K = 40, 80.
"""
import sys, math
import mpmath as mp
sys.path.insert(0, './scripts/16-mu5')
import evalmeasure as ev
import numpy as np

mp.mp.dps = 50
AL, LAM = mp.mpf(3) / 40, mp.mpf(37) / 40


def Gm(x):
    return x * mp.log(x) - x if x > 0 else mp.mpf(0)


def phi(s):
    return 6 * (Gm(s + AL) - Gm(s - AL)) - Gm(1 + s) - Gm(1 - s)


def S1(K):
    N = 3 * K // 40
    h = K - N
    worst = -mp.inf
    maxrel = 0
    for j in range(N + 1, K + 1):
        u = lambda i: mp.mpf(i) / K
        num = u(j) ** 4 * mp.fprod([(u(l) ** 2 - u(j) ** 2) for l in range(1, N + 1)]) ** 5
        den = mp.fprod([(u(k) ** 2 - u(j) ** 2) for k in range(N + 1, K + 1) if k != j])
        r_direct = abs(num / den)
        r_closed = 2 * mp.mpf(K) ** (2 * K - 12 * N - 6) * mp.factorial(j + N) ** 6 / (
            mp.factorial(j - N - 1) ** 6 * mp.factorial(K - j) * mp.factorial(K + j))
        maxrel = max(maxrel, abs(r_direct / r_closed - 1))
        bound = K * phi(u(j)) + mp.log(2) + mp.mpf(1) / 2 + 3 * mp.log(4 * mp.pi * K)
        worst = max(worst, mp.log(r_direct) - bound)
    return maxrel, worst


def logW(K, t):
    N = 3 * K // 40
    t = mp.mpf(t)
    s1 = 6 * mp.fsum(mp.log(t + (mp.mpf(l) / K) ** 2) for l in range(1, N + 1))
    s2 = mp.fsum(mp.log(t + (mp.mpf(j) / K) ** 2) for j in range(1, K + 1))
    y = K * mp.sqrt(t)
    w = (2 * mp.pi) ** 4 * y ** 5 / 12 * mp.nsum(lambda l: l ** 4 * mp.exp(-2 * mp.pi * l * y), [1, mp.inf])
    return s1 - s2 + mp.log(w) - mp.log(2 * mp.sqrt(t))


def Vm(t):
    t = mp.mpf(t)
    st = mp.sqrt(t)
    i1 = mp.log(1 + t) - 2 + 2 * st * mp.atan(1 / st)
    i2 = AL * mp.log(t + AL ** 2) - 2 * AL + 2 * st * mp.atan(AL / st)
    return 2 * mp.pi * st + i1 - 6 * i2


if __name__ == "__main__":
    for K in (40, 80, 400):
        rel, worst = S1(K)
        print(f"(S1) K={K}: max rel. error closed form vs product = {mp.nstr(rel, 3)};  max_j [log|r_j| - bound] = {mp.nstr(worst, 6)} (must be <= 0)")
    c0 = mp.log((2 * mp.pi) ** 4 / 24) - 2
    for K in (40, 80, 400):
        worst = -mp.inf
        for t in [7.2e-5, 1e-4, 1e-3, 5.9e-3, 0.02, 0.1, 0.3, 0.5, 0.75, 1.0, 1.5]:
            lhs = logW(K, t)
            rhs = -K * Vm(t) + 3 * mp.log(K) + 2 * mp.log(t) + c0
            worst = max(worst, rhs - lhs)
        print(f"(S2) K={K}: max_t [lower bound - log W_K(t)] = {mp.nstr(worst, 6)} (must be <= 0)")
    # (S3)
    rhoK = 75299683 / 10000000
    tmin = 7.12431588266478e-05
    Crho = 3 * math.log(4 * math.pi) - 2 * math.log(tmin) - math.log((2 * math.pi) ** 4 / 24) + 2 + math.log(2) + 0.5
    for K, logB in ((40, 292.5613171408397), (80, 591.0372993642756)):
        bnd = K * rhoK + 10 * math.log(K) + Crho + 3 * math.pi / (2 * K * math.sqrt(tmin))
        print(f"(S3) K={K}: exact log||B_K|| = {logB:.3f} (rate {logB/K:.4f});  Theorem-1 bound = {bnd:.3f} (rate {bnd/K:.4f});"
              f"  Fauzan's (C.1) without margin: kappa*K = {13.8658566788653*K:.1f}")
