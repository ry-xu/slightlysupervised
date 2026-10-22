"""
Verify (exactly, by series arithmetic mod p) every degree claim used in the
Mobius-independence proofs (notes/04-funcfield.md, sections 3-4):

 (a) deg(u - U_M)  = 1 - q^{M+1}                     (all q)
 (b) deg(e - S_N)  = -(N+1) q^{N+1}, leading coeff 1   (all q)
 (c) D_i | D_N for i <= N, deg D_N = N q^N
 (d) E_M * U_M = theta^{A_M},  A_M = (q^{M+1}-q)/(q-1) - M
 (e) q = 2:  V_M := (1 + x^{2^{M+1}-1}) / E_M(x)   satisfies
        u - V_M = x^{2A} (1+x) u (1 + O(x^A)),   A = 2^{M+1}-1,
     i.e. Dirichlet-quality approximation: |Q u - P| = |Q|^{-1}, deg Q = A.
     Also: u - V_M - x^{2A}(1+x)u has valuation >= 3A  (checked)
 (f) for general q, V_M-type truncations: T_M = E_M u = sum_{j<q} x^{j p1} + O(x^{q p1}),
     p1 = q^{M+1} - 1.
"""
import sys
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ff import *

def check(cond, msg):
    print(("OK   " if cond else "FAIL ") + msg)
    if not cond:
        raise SystemExit(1)

for q, p, Nprec in [(2, 2, 40000), (3, 3, 40000), (4, 2, 40000), (5, 5, 40000), (7, 7, 40000), (8, 2, 40000), (9, 3, 40000)]:
    print(f"--- q={q} (p={p}), precision x^{Nprec}")
    u = u_series(q, Nprec, p)
    e = e_series(q, Nprec, p)
    # (a)
    M = 1
    while q ** (M + 1) - 1 < Nprec:
        UM = U_M_series(q, M, Nprec, p)
        v = val((u - UM) % p)
        check(v == q ** (M + 1) - 1, f"(a) M={M}: val(u-U_M) = {v} = q^(M+1)-1")
        M += 1
    # (b)
    N = 0
    while (N + 1) * q ** (N + 1) < Nprec:
        SN = e_series(q, Nprec, p, upto=N)
        d = (e - SN) % p
        v = val(d)
        check(v == (N + 1) * q ** (N + 1) and d[v] == 1, f"(b) N={N}: val(e-S_N) = {v} = (N+1)q^(N+1), lc=1")
        N += 1
    # (c) divisibility (small i,N only: polynomial degrees N q^N)
    for Nn in range(0, 4):
        if Nn * q ** Nn > 3000:
            break
        DN = D_poly(q, Nn, p)
        check(DN.degree() == Nn * q ** Nn, f"(c) deg D_{Nn} = {DN.degree()}")
        for i in range(Nn + 1):
            Di = D_poly(q, i, p)
            check((DN % Di).is_zero(), f"(c) D_{i} | D_{Nn}")
    # (d)
    for M in range(1, 4):
        if q ** (M + 1) > 3000:
            break
        EM = E_poly(q, M, p)
        AM = (q ** (M + 1) - q) // (q - 1) - M
        check(EM.degree() == AM, f"(d) deg E_{M} = {EM.degree()} = A_M = {AM}")
    # (f) structure of T_M = E_M u
    for M in range(1, 6):
        p1 = q ** (M + 1) - 1
        if q * p1 + 5 >= Nprec:
            break
        TM = u.copy()
        for m in range(1, M + 1):
            mul_one_minus_xa(TM, q ** m - 1, p)
        expect = np.zeros(q * p1, dtype=np.int64)
        for j in range(q):
            expect[j * p1] = 1
        check(np.array_equal(TM[: q * p1] % p, expect), f"(f) M={M}: T_M = sum_(j<q) x^(j p1) + O(x^(q p1))")

# (e) the q = 2 Dirichlet approximants
q, p, Nprec = 2, 2, 200000
u = u_series(q, Nprec, p)
print("--- (e) q=2 approximants V_M")
for M in range(1, 17):
    A = 2 ** (M + 1) - 1
    if 3 * A + 2 >= Nprec:
        break
    EM = np.zeros(Nprec, dtype=np.int64); EM[0] = 1
    for m in range(1, M + 1):
        mul_one_minus_xa(EM, 2 ** m - 1, p)       # E_M(x) = prod (1 + x^{2^m-1})
    num = np.zeros(Nprec, dtype=np.int64); num[0] = 1; num[A] = 1
    VM = mul(num, inv(EM, p), p)
    du = (u - VM) % p
    v = val(du)
    # predicted: x^{2A}(1+x)u
    pred = np.zeros(Nprec, dtype=np.int64)
    onepx_u = u.copy(); onepx_u[1:] = (onepx_u[1:] + u[:-1]) % p
    pred[2 * A:] = onepx_u[: Nprec - 2 * A]
    resid = (du - pred) % p
    vr = val(resid)
    degEM = int(np.flatnonzero(EM)[-1])
    check(v == 2 * A, f"(e) M={M}: A={A}, val(u-V_M) = {v} = 2A; deg_x E_M = {degEM} < A; "
                      f"val(u - V_M - x^(2A)(1+x)u) = {vr} (>= 3A = {3*A})")
    check(vr is None or vr >= 3 * A, "    second-order term is O(x^{3A})")
print("all checks passed")
