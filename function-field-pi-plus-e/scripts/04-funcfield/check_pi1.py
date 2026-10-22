"""Check pi_1 := prod_{j>=1} (1 - [j]/[j+1]) (Thakur 2013 / Denis 1995 'pi') equals u * g with
g = prod_{j>=0} (1 - x^{(q-1) q^j})  (algebraic: g = (1 - x^{q-1}) g^q), so Denis's
'e and pi algebraically independent' transfers verbatim to e and u (hence to e and pi~)."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ff import *
for q, p in [(3, 3), (4, 2), (5, 5)]:
    N = 5000
    u = u_series(q, N, p)
    # 1 - [j]/[j+1] = [1]^{q^j}/[j+1] = (theta^q - theta)^{q^j} / (theta^{q^{j+1}} - theta)
    #   in x: x^{q^{j+1} - q^{j+1}} (1 - x^{q-1})^{q^j} / (1 - x^{q^{j+1}-1}) = (1 - x^{(q-1)q^j}) / (1 - x^{q^{j+1}-1})
    pi1 = np.zeros(N, dtype=np.int64); pi1[0] = 1
    j = 1
    while (q - 1) * q ** j < N or q ** (j + 1) - 1 < N:
        if (q - 1) * q ** j < N:
            mul_one_minus_xa(pi1, (q - 1) * q ** j, p)
        if q ** (j + 1) - 1 < N:
            div_one_minus_xa(pi1, q ** (j + 1) - 1, p)
        j += 1
    g = np.zeros(N, dtype=np.int64); g[0] = 1
    j = 0
    while (q - 1) * q ** j < N:
        mul_one_minus_xa(g, (q - 1) * q ** j, p); j += 1
    ug = mul(u, g, p)
    # note: u = prod_{m>=1} 1/(1-x^{q^m-1}) includes m=1 factor 1/(1-x^{q-1}); pi1's denominators start at q^2-1
    print(f'q={q}: pi_1 == u*g to x^{N}: {np.array_equal(pi1 % p, ug % p)}')
