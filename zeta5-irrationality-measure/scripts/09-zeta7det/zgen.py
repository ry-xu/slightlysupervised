"""General exact Hankel engine:  functional mu_X on rational functions of s with simple poles at s = -n^2
(n in a set of positive integers 'nodes'):
    mu_X(s^e) = mom[e] (rational),   mu_X(1/(s+n^2)) = c1(n) X + c0(n).
Weight W(s) = s^{e0} prod_n (s+n^2)^{e_n}, e_n in {-1,0,1,...};  Hankel G = [mu_X(W s^{i+j})], Delta = det G.

Families:
  zeta(k), k odd:  nodes n=j,  mom = (-1)^e B_{2e+2}(2e+k)!/((k-1)!(2e+2)!),
                   c1 = n^{k-1}, c0 = 1/(2n) - 1/(k-1) - n^{k-1} H_n^{(k)}.
  beta(k), k even (s = 4t, nodes n = 2j+1 = odd):  mom = 4^e (-1)^e (4^{e+1}-1) B_{2e+2}(2e+k)!/((k-1)!(2e+2)!),
                   c1 = (-1)^j n^{k-1}/2,  c0 = -(-1)^j n^{k-1} B_j/2 - 1/(4n),  B_j = sum_{m<j} (-1)^m/(2m+1)^k.
"""
import time, math
from math import factorial
from flint import fmpz, fmpq, fmpz_poly, fmpq_poly, fmpq_mat, arb, ctx
import sys
sys.path.insert(0, './scripts/09-zeta7det')
from zdet import bern, hankel_det, content, factor_content, log_fmpq, log_fmpz


class Family:
    def __init__(self, name, k):
        self.name, self.k = name, k

    def mom(self, L):
        k = self.k
        out = []
        for e in range(L + 1):
            v = bern(2 * e + 2) * fmpq(factorial(2 * e + k), factorial(k - 1) * factorial(2 * e + 2))
            if self.name == 'beta':
                v = v * (4 ** (e + 1) - 1) * (4 ** e)
            out.append(v if e % 2 == 0 else -v)
        return out

    def polevals(self, nodes):
        k = self.k
        res = {}
        if self.name == 'zeta':
            H = fmpq(0); last = 0
            for n in sorted(nodes):
                for v in range(last + 1, n + 1):
                    H += fmpq(1, v ** k)
                last = n
                res[n] = (fmpq(n ** (k - 1)), fmpq(1, 2 * n) - fmpq(1, k - 1) - (n ** (k - 1)) * H)
        elif self.name == 'beta':
            B = fmpq(0); m_done = 0
            for n in sorted(nodes):
                j = (n - 1) // 2
                while m_done < j:
                    B += fmpq((-1) ** m_done, (2 * m_done + 1) ** k); m_done += 1
                sg = 1 if j % 2 == 0 else -1
                res[n] = (fmpq(sg * n ** (k - 1), 2), fmpq(-sg * n ** (k - 1), 2) * B - fmpq(1, 4 * n))
        return res

    def value(self, prec):
        old = ctx.prec; ctx.prec = prec
        if self.name == 'zeta':
            v = arb(self.k).zeta()
        else:
            v = (arb(self.k).zeta(arb(1) / 4) - arb(self.k).zeta(arb(3) / 4)) / arb(4) ** self.k
        ctx.prec = old
        return v


class Comb35(Family):
    """L(R) = a mu^{(3)}(-t R) + b mu^{(5)}(R); under a zeta(3)+b zeta(5)=X:
       L(t^e) = -a mom3[e+1] + b mom5[e];  L(1/(t+j^2)) = j^4 X + c0(j)."""
    def __init__(self, a, b):
        self.name, self.k, self.a, self.b = 'comb35', 5, a, b
        self.f3, self.f5 = Family('zeta', 3), Family('zeta', 5)

    def mom(self, L):
        m3 = self.f3.mom(L + 1); m5 = self.f5.mom(L)
        return [-self.a * m3[e + 1] + self.b * m5[e] for e in range(L + 1)]

    def polevals(self, nodes):
        a, b = self.a, self.b
        res = {}
        H3 = fmpq(0); H5 = fmpq(0); last = 0
        m30 = self.f3.mom(0)[0]
        for n in sorted(nodes):
            for v in range(last + 1, n + 1):
                H3 += fmpq(1, v ** 3); H5 += fmpq(1, v ** 5)
            last = n
            c0 = a * (-m30 - (n ** 4) * H3 - fmpq(n * n, 2) + fmpq(n, 2)) + b * (-(n ** 4) * H5 - fmpq(1, 4) + fmpq(1, 2 * n))
            res[n] = (fmpq(n ** 4), c0)
        return res

    def value(self, prec):
        old = ctx.prec; ctx.prec = prec
        v = self.a * arb(3).zeta() + self.b * arb(5).zeta()
        ctx.prec = old
        return v


def entries_gen(fam, e, h, e0=0):
    poles = sorted(n for n, v in e.items() if v == -1)
    num = fmpz_poly([0] * e0 + [1])
    for n, v in e.items():
        if v > 0:
            num *= fmpz_poly([n * n, 1]) ** v
    den = fmpz_poly([1])
    for j in poles:
        den *= fmpz_poly([j * j, 1])
    Q0 = num // den
    q0 = [fmpq(int(c)) for c in Q0.coeffs()] if Q0.degree() >= 0 else []
    Wj = {}
    for j in poles:
        x = -j * j
        dv = fmpz(1)
        for i in poles:
            if i != j:
                dv *= (i * i - j * j)
        Wj[j] = fmpq(int(num(x)), 1) / dv
    pv = fam.polevals(poles)
    c1 = {j: Wj[j] * pv[j][0] for j in poles}
    c0 = {j: Wj[j] * pv[j][1] for j in poles}
    M = 2 * h - 1
    mom = fam.mom(M + len(q0) + 1)
    a, b, s = [], [], []
    for m in range(M):
        bm = sum(c1.values(), fmpq(0)); am = sum(c0.values(), fmpq(0)); sm = sum(Wj.values(), fmpq(0))
        muQ = fmpq(0)
        for ee in range(m):
            muQ += s[m - 1 - ee] * mom[ee]
        for ee, c in enumerate(q0):
            if c != 0:
                muQ += c * mom[ee + m]
        a.append(muQ + am); b.append(bm); s.append(sm)
        for j in poles:
            f = -j * j
            c1[j] *= f; c0[j] *= f; Wj[j] *= f
    return a, b


def eval_poly_at(P, fam):
    cs = [fmpz(int(c)) for c in P.coeffs()]
    mx = max(log_fmpz(c) for c in cs if c != 0)
    prec = int(mx / math.log(2)) + 400
    while True:
        old = ctx.prec; ctx.prec = prec
        z = fam.value(prec)
        v = arb(0)
        for c in reversed(cs):
            v = v * z + arb(c)
        if not v.contains(0):
            r = (float(v.abs_upper().log().mid()), 1 if v > 0 else -1)
            ctx.prec = old
            return r
        ctx.prec = old
        prec *= 2


def run_gen(fam, e, h, e0=0, B=None):
    t0 = time.time()
    a, b = entries_gen(fam, e, h, e0)
    D = hankel_det(a, b, h)
    c = content(D)
    P = D / c
    lv, sg = eval_poly_at(P, fam)
    res = dict(fam=fam.name, k=fam.k, h=h, deg=D.degree(), log_cont=log_fmpq(c), logP=lv, sign=sg,
               logDelta=lv + log_fmpq(c), t=time.time() - t0)
    if B:
        fac, rn, rd = factor_content(c, B)
        res['vp'] = {str(p): v for p, v in fac.items()}
    return res, D, P


def fauzan_like(fam, K, N, r):
    """zeros (order r-1) at the first N nodes, simple poles at the next K-N nodes."""
    e = {}
    for l in range(K):
        n = (l + 1) if fam.name in ('zeta', 'comb35') else (2 * l + 1)
        e[n] = (r - 1) if l < N else -1
    return e


if __name__ == "__main__":
    import json
    name = sys.argv[1]; k = int(sys.argv[2])
    Km, Nm, r, hm, nmax = float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), int(sys.argv[7])
    out = sys.argv[8] if len(sys.argv) > 8 else None
    fam = Family(name, k)
    for n in range(1, nmax + 1):
        K, N, h = int(round(Km * n)), int(round(Nm * n)), int(round(hm * n))
        res, D, P = run_gen(fam, fauzan_like(fam, K, N, r), h, B=4 * K + 10)
        res.update(n=n, K=K, N=N, r=r)
        print(json.dumps({x: res[x] for x in res if x != 'vp'}), flush=True)
        if out:
            with open(out, 'a') as fh:
                fh.write(json.dumps(res) + "\n")
