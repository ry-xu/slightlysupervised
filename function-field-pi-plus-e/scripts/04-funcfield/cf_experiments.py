"""
Continued-fraction experiments for u, e = e_C(1) (and pi~ for q=2).

Outputs (data/04-funcfield/):
  cf_<name>_q<q>.txt : lines "n deg(a_n) deg(Q_n)"
Prints summary statistics:
  * Thakur 1992 pattern for e (q=2): a_n = [v2(n)+1] = theta^(2^(v2(n)+1)) - theta  (degree 2^(v2(n)+1))
  * running exponent estimate  E_n = 2 + deg a_{n+1} / deg Q_n  (max over the tail)
  * parity statistics of deg a_n for u (q=2)
Also checks the normalization e_C(pi~) = 0 for q = 2 (pi~ = theta^2 u).
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ff import *

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../data/04-funcfield')
os.makedirs(OUT, exist_ok=True)

def v2(n):
    c = 0
    while n % 2 == 0:
        n //= 2; c += 1
    return c

def save(name, q, cf):
    with open(os.path.join(OUT, f'cf_{name}_q{q}.txt'), 'w') as f:
        for n, (da, dQ) in enumerate(cf):
            f.write(f'{n} {da} {dQ}\n')

def exponent_profile(cf, label):
    # E = 2 + deg a_{n+1}/deg Q_n ; report the largest values with deg Q_n >= 20
    rows = []
    for n in range(1, len(cf) - 1):
        dQ = cf[n][1]
        da_next = cf[n + 1][0]
        if dQ >= 20:
            rows.append((2 + da_next / dQ, n, dQ, da_next))
    rows.sort(reverse=True)
    print(f'  [{label}] top approximation events (E_n = 2 + deg a_(n+1)/deg Q_n):')
    for r in rows[:8]:
        print(f'     E={r[0]:.4f}  n={r[1]}  degQ_n={r[2]}  deg a_(n+1)={r[3]}')
    # tail behaviour: max over deg Q_n in dyadic windows
    print(f'  [{label}] max E_n over windows of deg Q_n:')
    W = [(2**j, 2**(j+1)) for j in range(5, 20)]
    for lo, hi in W:
        vals = [r[0] for r in rows if lo <= r[2] < hi]
        if vals:
            print(f'     degQ in [{lo},{hi}): max E = {max(vals):.4f}  (#conv={len(vals)})')

t0 = time.time()
# ------------------------------------------------------------ normalization check, q = 2
q, p, N = 2, 2, 4096
u = u_series(q, N, p)
# pi~ = theta^2 u ; pi~^(2^i) = x^{-2^(i+1)} u(x^(2^i)) ; e_C(pi~) = sum_i pi~^(2^i)/D_i
# compute x^{S} * e_C(pi~) with S = max shift to stay in power series land
tot = np.zeros(N, dtype=np.int64)
i = 0
S = 64     # we multiply everything by x^S
while True:
    v = i * 2 ** i - 2 ** (i + 1) + S   # valuation of x^S * pi~^(2^i)/D_i
    if v >= N:
        break
    if v < 0:
        raise SystemExit('increase S')
    ui = np.zeros(N, dtype=np.int64)
    ui[::2 ** i] = u[: (N + 2 ** i - 1) // 2 ** i]          # u(x^(2^i)) = u^(2^i)
    # x^S * pi~^(2^i)/D_i = x^{i 2^i - 2^(i+1) + S} * u(x^{2^i}) * prod_j (1 - x^{2^i - 2^j})^{-1}
    c = np.zeros(N, dtype=np.int64)
    c[v:] = ui[: N - v]
    for j in range(i):
        div_one_minus_xa(c, 2 ** i - 2 ** j, p)
    tot = (tot + c) % p
    i += 1
vz = val(tot)
print(f'normalization: x^{S} e_C(pi~) has first nonzero coefficient at index {vz} (None = vanishes to precision {N})')

# ------------------------------------------------------------ CFs
NPREC = {2: 60000, 3: 40000, 4: 40000, 5: 30000}
for q, p in [(2, 2), (3, 3), (4, 2), (5, 5)]:
    N = NPREC[q]
    print(f'=== q={q}, precision x^{N}  (t={time.time()-t0:.1f}s)')
    u = u_series(q, N, p)
    e = e_series(q, N, p)
    for name, ser in [('e', e), ('u', u)]:
        cf = cf_degrees(ser, p)
        save(name, q, cf)
        degs = [d for d, _ in cf[1:]]
        print(f'  {name}: {len(cf)-1} partial quotients, last deg Q = {cf[-1][1]}')
        if name == 'e' and q == 2:
            ok = all(degs[n - 1] == 2 ** (v2(n) + 1) for n in range(1, len(degs) + 1))
            print(f'  Thakur 1992 pattern deg a_n = 2^(v2(n)+1) for all {len(degs)} computed n: {ok}')
        if name == 'u' and q == 2:
            from collections import Counter
            C = Counter(degs)
            odd_ge3 = sum(v for k, v in C.items() if k % 2 == 1 and k >= 3)
            print(f'  u (q=2): degree histogram (top): {sorted(C.items())[:12]} ... max deg = {max(degs)}; '
                  f'#odd degrees >=3: {odd_ge3}')
        exponent_profile(cf, f'{name}, q={q}')
print(f'done in {time.time()-t0:.1f}s')
