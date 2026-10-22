"""Rigorous (arb ball arithmetic) certificate for Fauzan's analytic constant with a new nested-arcsine measure.
Lemma 6.1 form:  U = lam*M0 - I(rho) + C*,  valid if rho = sum_j c_j omega[a_j,b_j] (c_j > 0, nested intervals in (0,2),
sum c_j = lam = 37/40, lengths > 1/225) and 2U^rho(t) - V(t) <= M0 for all t >= 0.
Certification of the sup (Fauzan App. A.3): U^rho is nonincreasing on [0,a_1], constant on [a_1,b_1], nondecreasing on
[b_1,oo); V(t)=Phi(sqrt t) has a unique minimum q in (q-,q+) (Phi'' changes sign once at y^2 = 711/880, Phi'(0+)=-3pi,
Phi'(oo)=2pi, and we certify Phi'(sqrt q-) < 0 < Phi'(sqrt q+)).  On [l,r]:  2U - V <= 2 max(U(l),U(r)) - minV[l,r].
For t >= 2:  2U - V <= 1.3 log t + 6 al^3/t + 2/(3t) - 2 pi sqrt t  (decreasing; checked < M0 at t = 2).
All interval endpoints, a_j, b_j, c_j are exact rationals (denominator 10^12)."""
import sys, json
from fractions import Fraction as Fr
from flint import arb, ctx

ctx.prec = 160
LAM = Fr(37, 40); AL = Fr(3, 40)
lamA, alA = arb(37) / 40, arb(3) / 40
PI = arb.pi()


def A(fr):
    return arb(fr.numerator) / arb(fr.denominator)


def load(level):
    d = json.load(open('./data/09-zeta7det/arcsine_opt.json'))[str(level)]
    D = 10 ** 12
    a = [Fr(round(x * D), D) for x in d['a']]
    b = [Fr(round(x * D), D) for x in d['b']]
    c = [Fr(round(x * D), D) for x in d['c']]
    c[-1] += LAM - sum(c)  # exact mass
    assert all(x > 0 for x in c) and sum(c) == LAM
    for j in range(1, len(a)):
        assert a[j] < a[j - 1] and b[j] > b[j - 1], "not nested"
    assert all(bb - aa > Fr(1, 225) for aa, bb in zip(a, b)) and b[-1] < 2 and a[-1] > 0
    return a, b, c


def Uarc(t, a, b):
    """potential of the unit arcsine measure on [a,b] at t (arb); t is an exact rational."""
    if a <= t <= b:
        return (A(b - a) / 4).log()
    m = (a + b) / 2
    return (abs(A(t - m)) + (A((t - a) * (t - b))).sqrt()).log() - arb(2).log()


def Urho(t, a, b, c):
    s = arb(0)
    for aj, bj, cj in zip(a, b, c):
        s += A(cj) * Uarc(t, aj, bj)
    return s


def V(t):
    """(A.5): V(t) = log(1+t) - 6 al log(t+al^2) - 2 + 12 al + 2 sqrt t (pi + atan(1/sqrt t) - 6 atan(al/sqrt t))"""
    if t == 0:
        return -12 * alA * alA.log() - 2 + 12 * alA
    tA = A(t)
    st = tA.sqrt()
    return (1 + tA).log() - 6 * alA * (tA + alA * alA).log() - 2 + 12 * alA + 2 * st * (PI + (1 / st).atan() - 6 * (alA / st).atan())


def dPhi(y):
    return 2 * (PI + (1 / y).atan() - 6 * (alA / y).atan())


def certify(level, target_slack=Fr(1, 10 ** 5), maxdepth=30, verbose=True):
    a, b, c = load(level)
    # minimum of V: bracket q
    qm, qp = Fr(592050771, 10 ** 11), Fr(592050791, 10 ** 11)  # Fauzan's (A.6) bracket
    assert dPhi(A(qm).sqrt()) < 0 and dPhi(A(qp).sqrt()) > 0
    # lower bound V_* on [q-,q+] (Fauzan (A.6) form)
    qmA, qpA = A(qm), A(qp)
    Vstar = (1 + qmA).log() - 6 * alA * (qpA + alA * alA).log() - 2 + 12 * alA + \
        2 * qpA.sqrt() * (PI + (1 / qpA.sqrt()).atan() - 6 * (alA / qmA.sqrt()).atan())
    # the parenthesis is negative there, so using sqrt(q+) for the outer factor is the lower bound (as in Fauzan)
    assert (PI + (1 / qpA.sqrt()).atan() - 6 * (alA / qmA.sqrt()).atan()) < 0

    def minV(l, r):
        if r <= qm:
            return V(r)
        if l >= qp:
            return V(l)
        return Vstar

    a1, b1 = a[0], b[0]
    assert a1 < qm and qp < b1
    # numerical estimate of sup to set the target
    import numpy as np
    pts = sorted(set([Fr(0), Fr(2), qm, qp] + a + b))
    cache = {}

    def U(t):
        if t not in cache:
            cache[t] = Urho(t, a, b, c)
        return cache[t]

    def B(l, r):
        if a1 <= l and r <= b1:
            Uu = U(a1)
        else:
            Uu = U(l) if U(l) > U(r) else U(r)
            Uu = arb(max(U(l).upper(), U(r).upper()))
        return 2 * Uu - minV(l, r)

    # crude estimate of sup by dense sampling (float)
    samp = []
    for l, r in zip(pts[:-1], pts[1:]):
        for k in range(1, 8):
            t = l + (r - l) * k / 8
            samp.append(float((2 * U(t) - V(t)).mid()))
    est = max(samp)
    M0 = Fr(est).limit_denominator(10 ** 7) + target_slack
    if verbose:
        print(f"level {level}: J={len(a)}  sampled sup ~ {est:.7f}  target M0 = {float(M0):.7f}", flush=True)
    stack = list(zip(pts[:-1], pts[1:]))
    n_ok, n_split = 0, 0
    worst = arb(-100)
    while stack:
        l, r = stack.pop()
        val = B(l, r)
        if val.upper() <= A(M0).lower():
            n_ok += 1
            if val.upper() > worst.upper():
                worst = val
            continue
        if r - l < Fr(1, 2 ** maxdepth):
            raise RuntimeError(f"cannot certify on [{float(l)},{float(r)}]: B = {val}")
        m = (l + r) / 2
        stack += [(l, m), (m, r)]
        n_split += 1
    # tail t >= 2
    t2 = arb(2)
    f2 = arb(13) / 10 * t2.log() + 6 * alA ** 3 / t2 + arb(2) / (3 * t2) - 2 * PI * t2.sqrt()
    assert f2 < A(M0)
    # energy (A.2) and C* (6.3)
    S = [Fr(0)]
    for cj in c:
        S.append(S[-1] + cj)
    I = arb(0)
    for j in range(len(a)):
        I += A(S[j + 1] ** 2 - S[j] ** 2) * (A(b[j] - a[j]) / 4).log()
    Cs = -2 * lamA + 12 * alA * lamA * (1 - alA.log()) + 3 * lamA ** 2 - 2 * lamA ** 2 * (2 * lamA).log()
    Ucert = lamA * A(M0) - I + Cs
    if verbose:
        print(f"  certified: sup(2U-V) <= M0 = {float(M0):.7f} on {n_ok} intervals ({n_split} splits); worst B = {worst}")
        print(f"  I(rho) = {I}\n  C* = {Cs}\n  U <= {Ucert}")
    return dict(level=level, J=len(a), M0=str(M0), I=str(I), Cstar=str(Cs), U=Ucert, n_intervals=n_ok)


if __name__ == "__main__":
    lev = sys.argv[1] if len(sys.argv) > 1 else "2"
    slack = Fr(int(sys.argv[2]), 10 ** 6) if len(sys.argv) > 2 else Fr(1, 10 ** 5)
    res = certify(lev, target_slack=slack)
    print("U_cert upper bound:", res['U'].upper())
