"""
Numerical verification for papers/ff-pi-plus-e-q2.{md,tex}  (writer 22-write-ff2).

Independent re-implementation: does NOT import scripts/04-funcfield/ff.py.
Everything is exact arithmetic in F_2[[x]] (x = 1/theta) truncated mod x^n, or in F_2[theta],
using python-flint nmod_poly (modulus 2).

Sections (numbers refer to the paper):
  [S0]  basic facts: deg D_i, D_i | D_N, x-expansion of 1/D_i, combinatorial form of 1/u
  [S1]  normalization: e_C(theta^2 u) = 0, e_C(theta u) = theta, log_C(1) = u/(1+x) = pi~/[1]
  [S2]  Lemma A (exact form), M = 1..14
  [S3]  Lemma B, N = 0..10
  [S4]  window lemma (exact integers) and the 'window is forced' remark
  [S5]  e_C(1/f): tails and windows for several f
  [S6]  full error expansion  (aVS+bV+cS+d) - (aue+bu+ce+d) = x^{2A}(gamma + O(x^{2^N+2-C}))
  [S7]  Thakur's continued fraction of e_C(1); q_{2^N-1} = D_N
  [S8]  iota-invariance (theta -> theta+1) of all partial quotients of pi~ and e_C(1);
        parity of partial quotients of h*pi~ and pi~/h (Proposition 8.5), incl. u = pi~/theta^2
  [S9]  Hermite-Pade: no relation a ue + b u + c e + d = 0 with deg a,b,c,d <= D

Usage:  python -I verify_q2_paper.py [n_cf] [D_hp]       (defaults 80000, 300; about 10 s)
"""
import sys, time, random
from collections import Counter
import flint

N_CF = int(sys.argv[1]) if len(sys.argv) > 1 else 80000
D_HP = int(sys.argv[2]) if len(sys.argv) > 2 else 300

def P(c):
    return flint.nmod_poly(c, 2)

ONE = P([1])
T = P([0, 1])             # theta (when used as a polynomial in theta) or x (as a series variable)
S_ = T ** 2 + T           # theta^2 + theta = [1]

def mono(k):
    return P([0] * k + [1])

def val(f):
    for i, v in enumerate(f.coeffs()):
        if int(v):
            return i
    return None

def tr(f, n):
    return f.truncate(n)

def inv(f, n):
    return f.inverse_series_trunc(n)

def mul(f, g, n):
    return f.mul_low(g, n)

def frob(g, j, n):
    """g(x)^(2^j) = g(x^(2^j)) mod x^n."""
    step = 2 ** j
    out = [0] * n
    for i, v in enumerate(g.coeffs()):
        if i * step >= n:
            break
        out[i * step] = int(v)
    return P(out)

def check(cond, msg):
    print(("OK   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        raise SystemExit(1)

def s2(m):
    return bin(m).count("1")

def v2(m):
    return (m & -m).bit_length() - 1

# ------------------------------------------------------------------ the constants as x-series
def E_ser(M, n):
    """E_M(x) = prod_{m=1}^{M} (1 + x^{2^m-1}) mod x^n."""
    f = ONE
    for m in range(1, M + 1):
        if 2 ** m - 1 < n:
            f = tr(f * (ONE + mono(2 ** m - 1)), n)
    return f

def u_ser(n):
    M = 1
    while 2 ** (M + 1) - 1 < n:
        M += 1
    return inv(E_ser(M, n), n)

def U_partial(M, n):
    return inv(E_ser(M, n), n)

def invD(i, n):
    """1/D_i = x^{i 2^i} prod_{j<i} (1 + x^{2^i - 2^j})^{-1} mod x^n."""
    v = i * 2 ** i
    if v >= n:
        return P([])
    g = ONE
    for j in range(i):
        g = tr(g * (ONE + mono(2 ** i - 2 ** j)), n - v)
    return tr(inv(g, n - v).left_shift(v), n)

def e_ser(n, upto=None):
    s = P([])
    i = 0
    while i * 2 ** i < n and (upto is None or i <= upto):
        s = s + invD(i, n)
        i += 1
    return s

def D_theta(i):
    d = ONE
    for j in range(i):
        d = d * (T ** (2 ** i) + T ** (2 ** j))
    return d

# ------------------------------------------------------------------ continued fractions
def cf(series, n, shift=0):
    """CF of alpha = x^{-shift} * series (series known mod x^n). Returns [(a_k, q_k)], keeping
    only certified entries: index k is accepted once 2 deg q_k + deg a_{k+1} < n - shift."""
    Nn = n - 1
    cs = [int(c) for c in series.coeffs()] + [0] * n
    A = P([cs[Nn - j] for j in range(Nn + 1)]).left_shift(shift)
    B = mono(Nn)
    prec = n - shift
    a0, r = divmod(A, B)
    out = [(a0, ONE)]
    q_prev, q = P([]), ONE
    X, Y = B, r
    pending = None
    while not Y.is_zero():
        a, r = divmod(X, Y)
        if pending is not None:
            if 2 * pending[1].degree() + a.degree() < prec:
                out.append(pending)
            else:
                break
        q_new = a * q + q_prev
        pending = (a, q_new)
        q_prev, q = q, q_new
        X, Y = Y, r
    return out

def iota(p):
    return p.compose(T + ONE)

t0 = time.time()
print(f"verify_q2_paper.py  (python-flint {flint.__version__}; n_cf = {N_CF}, D_hp = {D_HP})")

# ------------------------------------------------------------------ [S0]
print("=== [S0] basic facts")
n0 = 5000
comb = [0] * n0
j = 0
while 2 * j - s2(j) < n0:
    comb[2 * j - s2(j)] = 1
    j += 1
check(E_ser(20, n0) == P(comb), "1/u = sum_{j>=0} x^{2j - s_2(j)} mod x^5000 (exponents distinct)")
for i in range(0, 5):
    Di = D_theta(i)
    vv = i * 2 ** i
    rev = P([int(c) for c in reversed(Di.coeffs())])           # theta^{-deg} D_i as x-polynomial
    ok_ser = tr(mul(rev, invD(i, n0).right_shift(vv), n0 - vv), n0 - vv) == ONE
    ok_div = all((Di % D_theta(k2)).is_zero() for k2 in range(i + 1))
    ok_rec = (i == 0) or Di == (T ** (2 ** i) + T) * D_theta(i - 1) ** 2
    check(Di.degree() == vv and ok_ser and ok_div and ok_rec,
          f"D_{i}: degree {vv} = i 2^i; D_i = [i] D_(i-1)^2; D_j | D_{i} (j<=i); x-expansion of 1/D_{i}")

# ------------------------------------------------------------------ [S1]
print("=== [S1] normalization of pi~ and log_C(1)")
n1 = 20000
u = u_ser(n1)
acc, acc2 = P([]), P([])
i = 0
while (i - 2) * 2 ** i + 2 < n1:
    g = ONE
    for jj in range(i):
        g = tr(g * (ONE + mono(2 ** i - 2 ** jj)), n1)
    core = mul(frob(u, i, n1), inv(g, n1), n1)                  # u^{2^i} * x^{-i 2^i}/D_i
    acc = acc + tr(core.left_shift((i - 2) * 2 ** i + 2), n1)   # x^2 (theta^2 u)^{2^i}/D_i
    if (i - 1) * 2 ** i + 1 < n1:
        acc2 = acc2 + tr(core.left_shift((i - 1) * 2 ** i + 1), n1)   # x (theta u)^{2^i}/D_i
    i += 1
check(acc.is_zero(), f"e_C(theta^2 u) = 0 mod x^{n1 - 2}")
check(tr(acc2, n1) == ONE, f"e_C(theta u) = theta mod x^{n1 - 1}  (the nonzero theta-torsion point)")
logc = P([])
i = 0
while 2 ** (i + 1) - 2 < n1:
    logc = logc + tr(U_partial(i, n1).left_shift(2 ** (i + 1) - 2), n1)   # 1/L_i
    i += 1
check(tr(mul(logc, ONE + mono(1), n1), n1) == u, f"log_C(1) = sum 1/L_i = u/(1+x) = pi~/[1] mod x^{n1}")
Ct = lambda z: T * z + z * z                   # Carlitz C_theta(z) = theta z + z^2;  C_[1] = C_theta^2 + C_theta
check(Ct(Ct(ONE)) + Ct(ONE) == P([]), "C_{[1]}(1) = 0: 1 is a [1]-torsion point of the Carlitz module")

# ------------------------------------------------------------------ [S2]
print("=== [S2] Lemma A: V_M = (1+x^A)/E_M(x), A = 2^(M+1)-1")
for M in range(1, 15):
    A = 2 ** (M + 1) - 1
    n = 4 * A + 12
    uu = u_ser(n)
    EM = E_ser(M, n)
    QM = P([int(c) for c in reversed(EM.coeffs())]).left_shift(M + 1)     # theta^A E_M(1/theta)
    V = mul(ONE + mono(A), inv(EM, n), n)
    d = uu + V
    r = d + tr(mul(uu, ONE + mono(1), n).left_shift(2 * A), n)
    r2 = r + tr(uu.left_shift(4 * A + 1), n)
    rhs = mul(mul(ONE + mono(2 * A), ONE + mono(2 * A + 1), n), U_partial(M + 2, n), n)
    check(EM.degree() == A - M - 1 and QM.degree() == A and val(d) == 2 * A and val(r) == 4 * A + 1
          and val(r2) == 4 * A + 3 and rhs == tr(V, n),
          f"M={M:2d} A={A:5d}: deg_x E_M = A-M-1, deg Q_M = A, V_M = (1+x^2A)(1+x^(2A+1))U_(M+2), "
          f"u-V_M = x^2A(1+x)u + x^(4A+1)u + O(x^(4A+3)) [valuations {val(d)}, {val(r)}, {val(r2)}]")

# ------------------------------------------------------------------ [S3]
print("=== [S3] Lemma B: e - S_N = x^Y (1 + x^(2^N) + O(x^(2^N+1))), Y = (N+1) 2^(N+1)")
n3 = 24000
e = e_ser(n3)
for N in range(0, 11):
    Y = (N + 1) * 2 ** (N + 1)
    if Y + 2 ** N + 2 >= n3:
        break
    dd = e + e_ser(n3, upto=N)
    check(val(dd) == Y and val(dd + mono(Y)) == Y + 2 ** N,
          f"N={N:2d}: val(e - S_N) = {val(dd)} = Y, next term at Y + 2^N = {Y + 2 ** N}")

# ------------------------------------------------------------------ [S4]
print("=== [S4] window lemma")
okw = True
for k in range(1, 17):
    N = 2 ** k - 1; M = N + k - 1
    A = 2 ** (M + 1) - 1; B = N * 2 ** N; Y = (N + 1) * 2 ** (N + 1)
    okw &= (A == (N + 1) * 2 ** N - 1) and (2 * A + 2 == Y) and (B - A == 1 - 2 ** N) \
        and (A + B - Y == -1 - 2 ** N) and (2 * A + 1 >= 2 ** N + 2)
check(okw, "k=1..16: A = (N+1)2^N - 1, 2A+2 = Y, B-A = 1-2^N, A+B-Y = -1-2^N, 2A+1 >= 2^N+2")
forced = True
for N in range(2, 40):
    good = [M for M in range(0, N + 10)
            if (2 ** (M + 1) - 1 - N * 2 ** N) > 1 and ((N + 2) * 2 ** N - (2 ** (M + 1) - 1)) > 1]
    expect = [N + v2(N + 1) - 1] if (N + 1) & N == 0 else []
    forced &= (good == expect)
check(forced, "2<=N<40: B+1 < A < Y-B-1 has a solution A = 2^(M+1)-1 iff N+1 = 2^k (then M = N+k-1)")

# ------------------------------------------------------------------ [S5]
print("=== [S5] e_C(1/f): tail x^{Y_f}(1+O(x^{2^N})), Y_f = (N+1+F)2^(N+1), window 2A+2 = Y_f")
for fpoly in [T, T + 1, T ** 2, T ** 2 + T + 1, T ** 3 + T + 1, T ** 4 + T + 1]:
    F = fpoly.degree()
    for k in range(2, 5):
        N = 2 ** k - 1 - F
        if N < 1:
            continue
        M = N + k - 1
        A = 2 ** (M + 1) - 1
        B = (N + F) * 2 ** N
        Y = (N + 1 + F) * 2 ** (N + 1)
        n = Y + 2 ** N + 40
        if n > 300000:
            continue
        frev = P([int(c) for c in reversed(fpoly.coeffs())])
        finv = tr(inv(frev, n).left_shift(F), n)
        tail = P([])
        i = N + 1
        while (i + F) * 2 ** i < n:
            tail = tail + mul(frob(finv, i, n), invD(i, n), n)
            i += 1
        v2nd = val(tail + mono(Y))
        check(val(tail) == Y and (v2nd is None or v2nd >= Y + 2 ** N) and 2 * A + 2 == Y and B - A == 1 - 2 ** N,
              f"f = {str(fpoly).replace('x', 'theta'):18s} k={k}: N={N:2d} M={M:2d} A={A:6d} B={B:6d} Y={Y:6d}")

# ------------------------------------------------------------------ [S6]
print("=== [S6] error expansion with random a,b,c,d (deg <= C = 6): x^C-scaled valuations")
random.seed(1)
for k in (2, 3):
    N = 2 ** k - 1; M = N + k - 1; A = 2 ** (M + 1) - 1
    C = 6
    n = 2 * A + 2 ** N + 200
    uu = u_ser(n); ee = e_ser(n); Sn = e_ser(n, upto=N)
    V = mul(ONE + mono(A), inv(E_ser(M, n), n), n)
    for trial in range(3):
        def X(_):   # x^C * p(theta) for a random p of degree <= C, as an x-polynomial
            cs = [random.randint(0, 1) for _ in range(C + 1)]
            return P([cs[C - jj] for jj in range(C + 1)])
        aX, bX, cX, dX = X(0), X(0), X(0), X(0)
        lhs = mul(aX, mul(V, Sn, n), n) + mul(bX, V, n) + mul(cX, Sn, n) \
            + mul(aX, mul(uu, ee, n), n) + mul(bX, uu, n) + mul(cX, ee, n)
        alphaX = mul(aX, ee, n) + bX
        betaX = mul(aX, uu, n) + cX
        gammaX = mul(mul(ONE + mono(1), uu, n), alphaX, n) + tr(betaX.left_shift(2), n)
        vd = val(tr(lhs + tr(gammaX.left_shift(2 * A), n), n))
        check(vd is not None and vd == 2 * A + 2 ** N + 2 + val(betaX),
              f"k={k} trial {trial}: val[x^C((aVS+bV+cS+d)-(aue+bu+ce+d)) - x^2A x^C gamma] = {vd} "
              f"= 2A + 2^N + 2 + val(x^C beta)  (predicted next term x^(2A+2+2^N) beta)")

# ------------------------------------------------------------------ [S7]
print("=== [S7] Thakur's continued fraction of e_C(1)")
ne = 20000
Le = cf(e_ser(ne), ne)
okT = Le[0][0] == ONE and all(a == T ** (2 ** (v2(m) + 1)) + T for m, (a, q) in enumerate(Le) if m >= 1)
okQ = all(Le[2 ** N - 1][1] == D_theta(N) for N in range(0, 14) if 2 ** N - 1 < len(Le))
check(okT and okQ, f"e_C(1) = [1; [1],[2],[1],[3],...], a_m = [v_2(m)+1] for all {len(Le)-1} certified m "
                   f"(deg q up to {Le[-1][1].degree()}); q_(2^N-1) = D_N")

# ------------------------------------------------------------------ [S8]
print(f"=== [S8] iota-invariance and parity of partial quotients (precision x^{N_CF})")
ucf = u_ser(N_CF)
Lpi = cf(ucf, N_CF, shift=2)
okI = all(iota(a) == a for (a, q) in Lpi)
check(okI, f"pi~: all {len(Lpi)} certified partial quotients are fixed by theta -> theta+1 (lie in F_2[theta^2+theta])")
Le2 = cf(e_ser(8000), 8000)
check(all(iota(a) == a for (a, q) in Le2), f"e_C(1): all {len(Le2)} certified partial quotients are iota-fixed")
# Proposition 8.5: for h of degree j, partial quotients of pi~/h and h*pi~ have degree <= j or = j mod 2.
nn = 20000
uh = u_ser(nn)
for name, num_mult, shift, j in [("u = pi~/theta^2", ONE, 0, 2), ("theta u = pi~/theta", ONE, 1, 1),
                                  ("pi~/theta^3", ONE, -1, 3), ("pi~/(theta^2+theta+1)", None, 2, 2),
                                  ("theta pi~", ONE, 3, 1), ("theta^2 pi~", ONE, 4, 2)]:
    if num_mult is None:
        h = T ** 2 + T + ONE                                  # pi~/h = x^{-2} u * x^2/(1+x+x^2)
        ser = mul(uh, inv(P([1, 1, 1]), nn), nn)
        Lh = cf(ser, nn, shift=0)
    elif shift >= 0:
        Lh = cf(uh, nn, shift=shift)
    else:
        Lh = cf(uh.left_shift(-shift), nn, shift=0)           # x * u = pi~/theta^3
    degs = [a.degree() for (a, q) in Lh[1:]]
    bad = [d for d in degs if not (d <= j or (d - j) % 2 == 0)]
    hist = Counter(degs)
    check(not bad, f"{name:24s}: {len(degs):5d} partial quotients, degrees all <= {j} or = {j} mod 2; "
                   f"small-degree counts {sorted(hist.items())[:5]}")
Lu = cf(ucf, N_CF, shift=0)
du = Counter(a.degree() for (a, q) in Lu[1:])
check(all(d == 1 or d % 2 == 0 for d in du),
      f"u: {sum(du.values())} certified partial quotients (deg q up to {Lu[-1][1].degree()}): all of degree 1 or even; "
      f"max degree {max(du)}; deg-1 fraction {du[1] / sum(du.values()):.3f}")

# ------------------------------------------------------------------ [S9]
print(f"=== [S9] Hermite-Pade: a ue + b u + c e + d = 0 with deg a,b,c,d <= {D_HP}?")
D = D_HP
Lx = 3 * D + 100
n = Lx + D + 10
uu = u_ser(n); ee = e_ser(n); ue = mul(uu, ee, n)
cu = [int(c) for c in uu.coeffs()] + [0] * n
ce = [int(c) for c in ee.coeffs()] + [0] * n
cue = [int(c) for c in ue.coeffs()] + [0] * n
rows = []
for jx in range(-D, Lx + 1):                 # coefficient of x^jx in sum_i theta^i (a_i ue + b_i u + c_i e + d_i)
    rr = 0
    for i in range(D + 1):
        m = jx + i
        if m < 0:
            continue
        if cue[m]: rr |= 1 << i
        if cu[m]: rr |= 1 << (D + 1 + i)
        if ce[m]: rr |= 1 << (2 * (D + 1) + i)
        if m == 0: rr |= 1 << (3 * (D + 1) + i)
    rows.append(rr)
rank, piv = 0, {}
for rr in rows:
    while rr:
        hb = rr.bit_length() - 1
        if hb in piv:
            rr ^= piv[hb]
        else:
            piv[hb] = rr; rank += 1; break
check(rank == 4 * (D + 1), f"rank {rank} = number of unknowns {4 * (D + 1)} ({len(rows)} coefficient equations): no relation")

print(f"=== all checks passed ({time.time() - t0:.1f} s)")
