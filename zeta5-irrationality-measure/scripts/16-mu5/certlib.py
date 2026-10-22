"""16-mu5: rigorous (arb ball arithmetic) certificates for nested arcsine measures.

Conventions (Fauzan, "zeta(5) is irrational", Zenodo 10.5281/zenodo.22826419, Sec. 6, App. A, App. C):
  alpha = 3/40, lambda = 37/40,
  V(t)    = 2 pi sqrt t + int_0^1 log(t+u^2) du - 6 int_0^alpha log(t+u^2) du          (6.1)
  U^rho(z) = int log|z-u| d rho(u)                                                    (Fauzan's sign)
  omega[a,b] = arcsine probability measure of [a,b];
  rho = sum_j c_j omega[a_j,b_j], intervals nested with index 0 innermost, c_j > 0.
  phi(s)  = 6[G(s+alpha) - G(s-alpha)] - G(1+s) - G(1-s),  G(x) = x log x - x  (= -Vtilde(-s^2)), notes/16-mu5.md.

Facts used (proved in notes/16-mu5.md, Sec. 3):
  (F1) U^{omega[a,b]}(t) = log((b-a)/4) on [a,b]; = log((|t-m| + sqrt((t-a)(t-b)))/2) off [a,b], m = (a+b)/2.
       It is concave on (-oo,a], on [a,b] and on [b,oo); non-increasing on (-oo,a], non-decreasing on [b,oo).
       Derivative off [a,b]: sign(t-m) / sqrt((t-a)(t-b)).
  (F2) Hence on any interval [l,r] containing no a_j, b_j in its interior, U^rho is concave (c_j > 0).
  (F3) V is decreasing on [0,t*] and increasing on [t*,oo) (unimodal; 15-redteam-mu5 sec. 2, re-proved in notes).
       V'(t)  = pi/sqrt t + arctan(1/sqrt t)/sqrt t - 6 arctan(al/sqrt t)/sqrt t,
       V''(t) = -pi/(2 t^1.5) - J(1,t) + 6 J(al,t),  J(a,t) = int_0^a du/(t+u^2)^2 = a/(2t(t+a^2)) + arctan(a/sqrt t)/(2 t^1.5).
  (F4) phi is concave on [alpha,1]: phi''(s) = -12 al/(s^2-al^2) - 2/(1-s^2) < 0.
"""
from fractions import Fraction as Fr
from bisect import bisect_left, bisect_right
import json, math, time
from flint import arb, fmpq, ctx

ctx.prec = 110

LAM = Fr(37, 40)
AL = Fr(3, 40)


def A(q):
    q = Fr(q)
    return arb(fmpq(q.numerator, q.denominator))


PI = arb.pi()
lamA, alA = A(LAM), A(AL)
CSTAR = -2 * lamA + 12 * alA * lamA * (1 - alA.log()) + 3 * lamA ** 2 - 2 * lamA ** 2 * (2 * lamA).log()   # (6.3)


# ------------------------------------------------------------------ the field V
def V_arb(t):
    """t: arb (point or ball) with t > 0, or exact 0"""
    if isinstance(t, (int, Fr)) and t == 0:
        return -2 - 12 * alA * alA.log() + 12 * alA
    if not isinstance(t, arb):
        t = A(t)
    st = t.sqrt()
    i1 = (1 + t).log() - 2 + 2 * st * (1 / st).atan()
    i2 = alA * (t + alA * alA).log() - 2 * alA + 2 * st * (alA / st).atan()
    return 2 * PI * st + i1 - 6 * i2


def dV_arb(t):
    if not isinstance(t, arb):
        t = A(t)
    st = t.sqrt()
    return (PI + (1 / st).atan() - 6 * (alA / st).atan()) / st


def _J(a, t):
    st = t.sqrt()
    return a / (2 * t * (t + a * a)) + (a / st).atan() / (2 * t * st)


def d2V_arb(t):
    """V'' on an arb ball t (t > 0)"""
    st = t.sqrt()
    return -PI / (2 * t * st) - _J(arb(1), t) + 6 * _J(alA, t)


def ball(l, r):
    """arb ball containing [l, r] (exact rationals)"""
    return A(l).union(A(r))


# ------------------------------------------------------------------ phi on the negative side
def G_arb(x):
    if isinstance(x, (Fr, int)) and x == 0:
        return arb(0)
    xa = A(x) if not isinstance(x, arb) else x
    return xa * xa.log() - xa


def phi_arb(s):
    """s exact rational in [AL, 1]"""
    s = Fr(s)
    return 6 * (G_arb(s + AL) - G_arb(s - AL)) - G_arb(1 + s) - G_arb(1 - s)


def dphi_arb(s):
    """phi'(s) = 6 log((s+al)/(s-al)) + log((1-s)/(1+s)), s in (AL, 1)"""
    s = Fr(s)
    return 6 * (A((s + AL) / (s - AL))).log() + (A((1 - s) / (1 + s))).log()


# ------------------------------------------------------------------ the measure
class Measure:
    def __init__(self, a, b, c, name='rho'):
        self.name = name
        self.a = [Fr(x) for x in a]
        self.b = [Fr(x) for x in b]
        self.c = [Fr(x) for x in c]
        self.J = len(self.a)
        self.aA = [A(x) for x in self.a]
        self.bA = [A(x) for x in self.b]
        self.cA = [A(x) for x in self.c]
        self.mA = [A((x + y) / 2) for x, y in zip(self.a, self.b)]
        self.kap = [(A(y - x) / 4).log() for x, y in zip(self.a, self.b)]       # log(L_j/4)
        # suffix sums of c_j * kappa_j  (components j >= j0 are "inside")
        self.suf = [arb(0)] * (self.J + 1)
        for j in range(self.J - 1, -1, -1):
            self.suf[j] = self.suf[j + 1] + self.cA[j] * self.kap[j]
        self.breaks = sorted(set(self.a + self.b))
        self.cache = {}

    def check(self, mass=LAM, verbose=True):
        ok = {}
        ok['positive'] = all(x > 0 for x in self.c)
        ok['mass_exact'] = (sum(self.c) == mass)
        ok['nested'] = all(self.a[j] < self.a[j - 1] for j in range(1, self.J)) and \
            all(self.b[j] > self.b[j - 1] for j in range(1, self.J)) and self.a[0] < self.b[0]
        ok['support_in_(0,2)'] = self.a[-1] > 0 and self.b[-1] < 2
        ok['min_len'] = float(min(y - x for x, y in zip(self.a, self.b)))
        ok['reg_const'] = sum(8 * float(cj) * math.sqrt(2) / (math.pi * math.sqrt(float(bj - aj)))
                              for aj, bj, cj in zip(self.a, self.b, self.c))
        if verbose:
            print(f"   [{self.name}] J={self.J} positive={ok['positive']} mass==lambda: {ok['mass_exact']} nested={ok['nested']}"
                  f" supp=[{float(self.a[-1]):.6g},{float(self.b[-1]):.8f}] in (0,2): {ok['support_in_(0,2)']}"
                  f" min_len={ok['min_len']:.3g} reg_const={ok['reg_const']:.3f}")
        return ok

    def j0(self, t):
        """first index j with a_j <= t <= b_j (components j >= j0 contain t); J if none"""
        # nested: containment is monotone in j
        lo, hi = 0, self.J
        while lo < hi:
            mid = (lo + hi) // 2
            if self.a[mid] <= t <= self.b[mid]:
                hi = mid
            else:
                lo = mid + 1
        return lo

    def U(self, t):
        """U^rho(t) at exact rational real t"""
        v = self.cache.get(t)
        if v is not None:
            return v
        t = Fr(t)
        j0 = self.j0(t)
        tA = A(t)
        v = self.suf[j0]
        for j in range(j0):
            d = tA - self.mA[j]
            q = A((t - self.a[j]) * (t - self.b[j]))
            v += self.cA[j] * ((abs(d) + q.sqrt()) / 2).log()
        self.cache[t] = v
        return v

    def dU(self, t):
        """U^rho'(t), t exact rational, not a breakpoint"""
        t = Fr(t)
        j0 = self.j0(t)
        v = arb(0)
        for j in range(j0):
            q = A((t - self.a[j]) * (t - self.b[j])).sqrt()
            if t > self.b[j]:
                v += self.cA[j] / q
            else:
                v -= self.cA[j] / q
        return v

    def dU_neg_upper_abs(self, x):
        """for x < a_last (x exact rational): |U'(x)| = sum c_j / sqrt((a_j-x)(b_j-x)) (exact enclosure)"""
        v = arb(0)
        for j in range(self.J):
            v += self.cA[j] / A((self.a[j] - x) * (self.b[j] - x)).sqrt()
        return v

    def energy(self):
        """(A.2): I(rho) = sum_j (S_j^2 - S_{j-1}^2) log(L_j/4), S_j cumulative from innermost"""
        S = Fr(0)
        I = arb(0)
        for j in range(self.J):
            S1 = S + self.c[j]
            I += A(S1 * S1 - S * S) * self.kap[j]
            S = S1
        return I

    @staticmethod
    def from_json(fn, key=None, denom=10 ** 15, mass=LAM, name=None):
        d = json.load(open(fn))
        if key is not None:
            d = d[key]
        a = [Fr(round(x * denom), denom) for x in d['a']]
        b = [Fr(round(x * denom), denom) for x in d['b']]
        c = [Fr(round(x * denom), denom) for x in d['c']]
        c[-1] += mass - sum(c)
        return Measure(a, b, c, name=name or fn.split('/')[-1])

    def to_json(self, fn):
        json.dump(dict(a=[str(x) for x in self.a], b=[str(x) for x in self.b], c=[str(x) for x in self.c]), open(fn, 'w'))

    @staticmethod
    def from_exact_json(fn, name=None):
        d = json.load(open(fn))
        return Measure([Fr(x) for x in d['a']], [Fr(x) for x in d['b']], [Fr(x) for x in d['c']], name=name or fn.split('/')[-1])


# ------------------------------------------------------------------ V-minimum bracket (unimodality)
def v_min_bracket():
    """t* = argmin V = y*^2 where Phi'(y*) = 0, Phi'(y) = 2(pi + arctan(1/y) - 6 arctan(al/y)).
    Returns (qlo, qhi, lower bound of min V)."""
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
    dPhi = lambda y: 2 * (PI + (1 / A(y)).atan() - 6 * (alA / A(y)).atan())
    assert dPhi(ylo) < 0 and dPhi(yhi) > 0
    qlo, qhi = ylo * ylo, yhi * yhi
    return qlo, qhi, V_arb(ball(qlo, qhi)).lower()
