"""[KNOWN, Thakur 1992] q=2: e_C(1) = [1; [1],[2],[1],[3],[1],[2],[1],[4],...], a_n = [v2(n)+1],
[j] = theta^(2^j) + theta.  Check the exact partial quotients (not just degrees)."""
import sys, os
import flint
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ff import *
p, N = 2, 20000
e = e_series(2, N, p)
A = theta_poly_from_xpoly(e, N - 1, p)
B = flint.nmod_poly([0] * (N - 1) + [1], p)
T = flint.nmod_poly([0, 1], p)
a0, r = divmod(A, B); A, B = B, r
print('a0 =', a0)
ok, n, degQ = True, 0, 0
while True:
    a, r = divmod(A, B)
    if 2 * degQ + a.degree() >= N - 3:
        break
    n += 1
    v = (n & -n).bit_length() - 1
    expect = T ** (2 ** (v + 1)) + T
    if a != expect:
        ok = False; print('mismatch at n =', n); break
    degQ += a.degree(); A, B = B, r
print(f'exact match a_n = [v2(n)+1] for n = 1..{n}: {ok}  (deg Q_n up to {degQ})')
