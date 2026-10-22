"""15-redteam-mu5: independent rigorous (arb ball arithmetic) certificate of the analytic constant.

Claim under test (09): for 09's 121-component nested arcsine measure rho,
     sup_{t>=0} (2U^rho(t) - V(t)) <= M0 = -6.6600532   and   U = lam*M0 - I(rho) + C* <= -1.3843989.

Method (my own implementation; same mathematical principle as Fauzan App. A.3, different code/bracketing):
  * U^rho is quasi-convex on [0,oo) for nested intervals (nonincreasing on [0,a_1], constant on [a_1,b_1],
    nondecreasing on [b_1,oo)), so  max_{[l,r]} U^rho = max(U^rho(l), U^rho(r)).
  * V(t) = Phi(sqrt t) is unimodal: Phi''(y) = -2/(1+y^2) + 12 al/(al^2+y^2) > 0 iff y^2 < 711/880,
    Phi'(0+) = -3 pi, Phi'(oo) = 2 pi > 0 and Phi' is decreasing on (sqrt(711/880), oo) towards 2 pi, so Phi' has
    exactly one zero y*.  I bracket y* in [y_lo, y_hi] (arb sign certificate) and bound
        min_{[l,r]} V >= V(r) if r <= y_lo^2,  V(l) if l >= y_hi^2,  else  lower(V(ball [y_lo^2, y_hi^2])).
    (V on a ball is evaluated by arb, i.e. a rigorous interval enclosure -- different from Fauzan's (A.6).)
  * Adaptive bisection of [0,2] starting from the breakpoints {0,2,a_j,b_j,y_lo^2,y_hi^2}.
  * t >= 2: since supp rho in (0,t): U^rho(t) <= lam log t;  V(t) >= 2 pi sqrt t + (1-6 al) log t - 2 al^3/t.
    So 2U-V <= (2lam-1+6al) log t + 2al^3/t - 2pi sqrt t =: g(t), decreasing on [2,oo); check g(2) < M0.
  * I(rho) by the nesting formula (A.2) (re-derived: U^{omega_j} = log(L_j/4) on [a_j,b_j] which contains
    all earlier intervals); cross-checked in explore_U.py by Gauss-Chebyshev quadrature (agreement 2e-8).
  * C* = -2lam + 12 al lam (1 - log al) + 3 lam^2 - 2 lam^2 log(2 lam)  (re-derived from Stirling for S_K).

usage: python u_certify.py [09|fauzan] [M0 as decimal string]
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import arb, fmpq, ctx
from measure_common import LAM, AL, load_09, load_fauzan, admissibility

ctx.prec = 128


def A(q):
    q = Fr(q)
    return arb(fmpq(q.numerator, q.denominator))


PI = arb.pi()
lamA, alA = A(LAM), A(AL)


def U_arc(t, a, b):
    """potential int log|t-u| d omega_[a,b](u) of the unit arcsine measure; t, a, b exact rationals."""
    if a <= t <= b:
        return (A(b - a) / 4).log()
    m = (a + b) / 2
    # |t-m| + sqrt((t-a)(t-b)), all exact rationals inside sqrt
    return ((A(abs(t - m)) + A((t - a) * (t - b)).sqrt()) / 2).log()


def make_U(a, b, c):
    cache = {}
    cA = [A(x) for x in c]

    def U(t):
        v = cache.get(t)
        if v is None:
            v = arb(0)
            for aj, bj, cj in zip(a, b, cA):
                v += cj * U_arc(t, aj, bj)
            cache[t] = v
        return v
    return U, cache


def V_ball(tb):
    """V on an arb ball tb (tb > 0 required); rigorous enclosure."""
    st = tb.sqrt()
    i1 = (1 + tb).log() - 2 + 2 * st * (1 / st).atan()
    i2 = alA * (tb + alA * alA).log() - 2 * alA + 2 * st * (alA / st).atan()
    return 2 * PI * st + i1 - 6 * i2


def V(t):
    if t == 0:
        return -2 - 12 * alA * alA.log() + 12 * alA
    return V_ball(A(t))


def dPhi(y):  # y exact rational > 0
    yA = A(y)
    return 2 * (PI + (1 / yA).atan() - 6 * (alA / yA).atan())


def bracket_min_V():
    # unimodality facts
    y0 = Fr(711, 880)
    assert dPhi(Fr(9, 10)) > 0  # 0.81 > 711/880 = 0.8079..., so this point is beyond the inflection; Phi' > 0 there
    # float root of Phi' then rational bracket
    import math
    al = 3 / 40
    f = lambda y: math.pi + math.atan(1 / y) - 6 * math.atan(al / y)
    lo, hi = 0.01, 0.5
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    ylo = Fr(int(lo * 10 ** 14) - 1, 10 ** 14)
    yhi = Fr(int(hi * 10 ** 14) + 2, 10 ** 14)
    assert dPhi(ylo) < 0 and dPhi(yhi) > 0, "bracket failed"
    qlo, qhi = ylo * ylo, yhi * yhi
    ball = A(qlo).union(A(qhi))
    Vmin_lb = V_ball(ball).lower()
    return qlo, qhi, Vmin_lb


def certify(a, b, c, M0, verbose=True, maxsplit=10 ** 6):
    t0 = time.time()
    adm = admissibility(a, b, c, verbose=verbose)
    assert adm['all_c_positive'] and adm['mass_exact'] and adm['nested'] and adm['len_gt_1_225'] and adm['support_in_(0,2)']
    U, cache = make_U(a, b, c)
    qlo, qhi, Vmin_lb = bracket_min_V()
    if verbose:
        print(f"   V-min bracket t in [{float(qlo):.15f}, {float(qhi):.15f}], min V >= {Vmin_lb}")
    M0A = A(M0)

    def Vmin(l, r):
        if r <= qlo:
            return V(r).lower()
        if l >= qhi:
            return V(l).lower()
        return Vmin_lb

    def bound(l, r):
        ul, ur = U(l).upper(), U(r).upper()
        um = ul if ul > ur else ur
        return 2 * um - Vmin(l, r)

    pts = sorted(set([Fr(0), Fr(2), qlo, qhi] + list(a) + list(b)))
    stack = list(zip(pts[:-1], pts[1:]))
    n_ok = n_split = 0
    worst = None
    minw = Fr(1)
    while stack:
        l, r = stack.pop()
        B = bound(l, r)
        if B < M0A:  # certain comparison in arb
            n_ok += 1
            if worst is None or B.upper() > worst:
                worst = B.upper()
            minw = min(minw, r - l)
            continue
        n_split += 1
        if n_split > maxsplit or r - l < Fr(1, 2 ** 60):
            raise RuntimeError(f"FAILED to certify on [{float(l)}, {float(r)}]: bound {B} vs M0 {float(M0)}")
        m = (l + r) / 2
        stack.append((l, m))
        stack.append((m, r))
    # t >= 2
    t2 = arb(2)
    g2 = (2 * lamA - 1 + 6 * alA) * t2.log() + 2 * alA ** 3 / t2 - 2 * PI * t2.sqrt()
    assert g2 < M0A, "t>=2 tail bound fails"
    assert (2 * lamA - 1 + 6 * alA) / arb(2).sqrt() < PI  # g'(t) < 0 on [2,oo): (2lam-1+6al)/t <= that/sqrt(2)/sqrt t < pi/sqrt t
    # energy (A.2) and C*
    S = Fr(0)
    I = arb(0)
    for aj, bj, cj in zip(a, b, c):
        S1 = S + cj
        I += A(S1 * S1 - S * S) * (A(bj - aj) / 4).log()
        S = S1
    Cs = -2 * lamA + 12 * alA * lamA * (1 - alA.log()) + 3 * lamA ** 2 - 2 * lamA ** 2 * (2 * lamA).log()
    Ub = lamA * M0A - I + Cs
    if verbose:
        print(f"   certified sup_(t>=0)(2U-V) < M0 = {float(M0):.10f}: {n_ok} final intervals, {n_split} splits,"
              f" {len(cache)} U-evaluations, min width {float(minw):.3g}, worst interval bound {worst}")
        print(f"   g(2) = {g2}")
        print(f"   I(rho) = {I}")
        print(f"   C*     = {Cs}")
        print(f"   U <= lam*M0 - I + C* = {Ub}   (upper end {Ub.upper()})")
        print(f"   time {time.time()-t0:.1f}s")
    return dict(M0=M0, I=I, Cs=Cs, U=Ub, n=n_ok, splits=n_split)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else '09'
    M0 = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr('-6.6600532')
    if which == 'fauzan':
        a, b, c = load_fauzan()
    else:
        a, b, c = load_09(which if which != '09' else '3p')
    print(f"== measure {which}, J = {len(a)}, target M0 = {M0} = {float(M0)}")
    certify(a, b, c, M0)
