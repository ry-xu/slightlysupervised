"""
Exact check of the q = 2 window data used in the proof of Theorem 2 (notes/04-funcfield.md).

For f monic in F_2[theta], F = deg f (f = 1 gives e = e_C(1)):
   S_N^f = sum_{i<=N} f^{-2^i}/D_i  has denominator f^{2^N} D_N of degree B = (N+F) 2^N,
   e_C(1/f) - S_N^f = x^{Y} (1 + O(x^{2^N})),   Y = (N+1+F) 2^{N+1}.
Choose N with N+F+1 = 2^k, M+1 = N+k, A = 2^{M+1}-1 = (N+F+1)2^N - 1. Then
   deg(Q_M D_N^f delta_u) = B - A = 1 - 2^N < 0,
   deg(Q_M D_N^f delta_e) = A + B - Y = -2^N - 1 < 0,
   2A + 2 = Y   (the 'tie' that forces the second-order analysis).
We verify the series facts exactly mod 2.
"""
import sys, os, time
import numpy as np
import flint
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ff import *

p = 2
T = flint.nmod_poly([0, 1], 2)

def inv_f_series(fpoly, Nprec):
    """1/f as a series in x (f monic deg F): 1/f = x^F * 1/frev(x)."""
    F = fpoly.degree()
    rev = [int(c) for c in reversed(fpoly.coeffs())]     # f(1/x) x^F
    r = np.zeros(Nprec, dtype=np.int64); r[:len(rev)] = rev
    ri = inv(r, p)
    out = np.zeros(Nprec, dtype=np.int64)
    out[F:] = ri[:Nprec - F]
    return out

def frob(c, j):
    """c(x)^(2^j) = c(x^(2^j)) for F_2 coefficients."""
    out = np.zeros(len(c), dtype=np.int64)
    step = 2 ** j
    m = (len(c) + step - 1) // step
    out[::step] = c[:m]
    return out

for fpoly, kmax in [(flint.nmod_poly([1], 2), 4), (T, 3), (T + 1, 3), (T**2 + T + 1, 3)]:
    F = fpoly.degree()
    print(f'=== f = {fpoly}  (F = {F})')
    for k in range(2, kmax + 1):
        N = 2 ** k - 1 - F
        if N < 1:
            continue
        M = N + k - 1
        A = 2 ** (M + 1) - 1
        B = (N + F) * 2 ** N
        Y = (N + 1 + F) * 2 ** (N + 1)
        Nprec = Y + 2 ** N + 50
        t0 = time.time()
        finv = inv_f_series(fpoly, Nprec)
        # e_C(1/f) - S_N^f = sum_{i>N} f^{-2^i}/D_i : compute leading term i=N+1 and i=N+2 contribution
        tail = np.zeros(Nprec, dtype=np.int64)
        i = N + 1
        while i * 2 ** i < Nprec:
            term = mul(frob(finv, i), inv_D_series(2, i, Nprec, p), p)
            tail = (tail + term) % p
            i += 1
        v = val(tail)
        # next nonzero after the leading one
        nz = np.flatnonzero(tail)
        gap = int(nz[1] - nz[0]) if len(nz) > 1 else None
        ok1 = (v == Y)
        ok2 = gap is None or gap >= 2 ** N
        # u - V_M
        if 2 * A + 4 < Nprec:
            u = u_series(2, 2 * A + 64, p)
            EM = np.zeros(2 * A + 64, dtype=np.int64); EM[0] = 1
            for m in range(1, M + 1):
                mul_one_minus_xa(EM, 2 ** m - 1, p)
            num = np.zeros(2 * A + 64, dtype=np.int64); num[0] = 1; num[A] = 1
            du = (u - mul(num, inv(EM, p), p)) % p
            vu = val(du)
        else:
            vu = 'skipped'
        print(f'  k={k}: N={N}, M={M}, A={A}, B={B}, Y={Y}; 2A+2==Y: {2*A+2==Y}; '
              f'val(e(1/f)-S_N)={v} (==Y: {ok1}); next term gap={gap} (>=2^N={2**N}: {ok2}); '
              f'val(u-V_M)={vu} (==2A: {vu==2*A}); margins: B-A={B-A}, A+B-Y={A+B-Y}; {time.time()-t0:.1f}s',
              flush=True)
