"""18-redteam-mu96: rigorous (arb) bounds for f(t) = 2U^rho(t) - V(t) on [0,2].

Own scheme (differs from 16's tangent-line/concavity scheme):
  On a sub-piece [p,q] of a breakpoint-free gap (l,r):
    f = g + (singular component term), g = sum_{nonsingular outside j} 2c_j U_j + 2*inside_const - V,
    g(t) in g(c) +/- sup_{[p,q]}|g'| (q-p)/2        (mean value theorem; g' enclosed by arb ball evaluation)
    singular component (endpoint exactly at p or q; only when p==l or q==r): monotone on the gap,
      so bounded above/below by its values at p and q (exact point evaluations).
  For [0,delta] with delta <= A_{J-1} < alpha^2: U decreasing (all components left of t) and V decreasing
  (V'(t) numerator pi + atan(1/sqrt t) - 6 atan(alpha/sqrt t) < 0 when sqrt t < alpha), so
  sup f <= 2U(0) - V(delta).
Everything is evaluated in arb with exact rational inputs; comparisons are arb-certain.
"""
import sys, time
from flint import arb, fmpq, ctx
sys.path.insert(0, './scripts/18-redteam-mu96')
from common import load_measure, ALPHA, LAM

ctx.prec = 128


class RealCert:
    def __init__(self, A, B, C):
        self.J = len(A)
        self.Aq, self.Bq, self.Cq = A, B, C
        self.A = [arb(x) for x in A]
        self.B = [arb(x) for x in B]
        self.C = [arb(x) for x in C]
        self.M = [(arb(A[j]) + arb(B[j])) / 2 for j in range(self.J)]
        self.two = arb(2)
        # inside constants: suffix sums of c_j log((b_j-a_j)/4)
        ins = [arb(0)] * (self.J + 1)
        for j in range(self.J - 1, -1, -1):
            ins[j] = ins[j + 1] + self.C[j] * (arb(B[j] - A[j]) / 4).log()
        self.ins = ins
        self.al = arb(ALPHA)
        self.pi = arb.pi()
        self.nevals = 0

    # ---- single components
    def Uj(self, j, t, side):
        """U_j(t) for t outside [a_j,b_j]; side = 'L' (t < a_j) or 'R' (t > b_j); t arb."""
        if side == 'L':
            return ((self.M[j] - t + ((self.A[j] - t) * (self.B[j] - t)).sqrt()) / 2).log()
        else:
            return ((t - self.M[j] + ((t - self.A[j]) * (t - self.B[j])).sqrt()) / 2).log()

    def dUj(self, j, t, side):
        if side == 'L':
            return -1 / ((self.A[j] - t) * (self.B[j] - t)).sqrt()
        else:
            return 1 / ((t - self.A[j]) * (t - self.B[j])).sqrt()

    # ---- V and V'
    def V(self, t):
        al = self.al
        st = t.sqrt()
        return (1 + t).log() - 6 * al * (t + al * al).log() - 2 + 12 * al + 2 * st * (
            self.pi + (1 / st).atan() - 6 * (al / st).atan())

    def V0(self):
        al = self.al
        return -12 * al * al.log() - 2 + 12 * al

    def dV(self, t):
        st = t.sqrt()
        return (self.pi + (1 / st).atan() - 6 * (self.al / st).atan()) / st

    # ---- gap classification
    def classify(self, l, r):
        """Return (side, k): components 0..k-1 are outside of type side, k..J-1 inside, for t in (l,r)."""
        Aq, Bq, J = self.Aq, self.Bq, self.J
        # left-type: A_j >= r ; right-type: B_j <= l
        kL = 0
        while kL < J and Aq[kL] >= r:
            kL += 1
        kR = 0
        while kR < J and Bq[kR] <= l:
            kR += 1
        assert not (kL > 0 and kR > 0)
        if kL > 0:
            # check inside for the rest
            assert all(Aq[j] <= l and Bq[j] >= r for j in range(kL, J))
            return 'L', kL
        if kR > 0:
            assert all(Aq[j] <= l and Bq[j] >= r for j in range(kR, J))
            return 'R', kR
        assert all(Aq[j] <= l and Bq[j] >= r for j in range(J))
        return None, 0

    def bounds_piece(self, p, q, l, r, side, k):
        """Rigorous (lower, upper) bounds of f on [p,q] subset [l,r] (gap with classification side,k)."""
        self.nevals += 1
        c = (p + q) / 2
        ca = arb(c)
        ball = arb(p).union(arb(q))
        hh = arb(q - p) / 2
        sing = None
        if side == 'R' and p == l and self.Bq[k - 1] == l:
            sing = k - 1
        if side == 'L' and q == r and self.Aq[k - 1] == r:
            sing = k - 1
        g = 2 * self.ins[k] - self.V(ca)
        dg = -self.dV(ball)
        for j in range(k):
            if j == sing:
                continue
            g += 2 * self.C[j] * self.Uj(j, ca, side)
            dg += 2 * self.C[j] * self.dUj(j, ball, side)
        err = dg.abs_upper() * hh
        lo = g - err
        hi = g + err
        if sing is not None:
            # value at its own endpoint is exactly log((b-a)/4) (continuity of the arcsine potential)
            endv = (arb(self.Bq[sing] - self.Aq[sing]) / 4).log()
            if side == 'R':
                up = endv
                uq = self.Uj(sing, arb(q), side)
            else:
                up = self.Uj(sing, arb(p), side)
                uq = endv
            if side == 'R':  # increasing on the gap
                lo += 2 * self.C[sing] * up
                hi += 2 * self.C[sing] * uq
            else:  # decreasing
                lo += 2 * self.C[sing] * uq
                hi += 2 * self.C[sing] * up
        return lo.lower(), hi.upper()

    def f_point(self, t):
        """f at an exact rational point t (classification by direct comparison)."""
        ta = arb(t)
        s = -self.V(ta) if t > 0 else -self.V0()
        for j in range(self.J):
            if t < self.Aq[j]:
                s += 2 * self.C[j] * self.Uj(j, ta, 'L')
            elif t > self.Bq[j]:
                s += 2 * self.C[j] * self.Uj(j, ta, 'R')
            else:
                s += 2 * self.C[j] * (arb(self.Bq[j] - self.Aq[j]) / 4).log()
        return s


def breakpoints(rc, lo, hi):
    pts = set([lo, hi])
    for x in rc.Aq + rc.Bq:
        if lo < x < hi:
            pts.add(x)
    return sorted(pts)


def certify(rc, lo, hi, target, mode, log=print, maxdepth=60):
    """mode 'sup': prove f < target on [lo,hi]; mode 'inf': prove f > target on [lo,hi].
    Returns (ok, npieces, worst)."""
    pts = breakpoints(rc, lo, hi)
    npieces = 0
    worst = None
    t0 = time.time()
    for gi in range(len(pts) - 1):
        l, r = pts[gi], pts[gi + 1]
        side, k = rc.classify(l, r)
        stack = [(l, r, 0)]
        while stack:
            p, q, dep = stack.pop()
            if mode == 'sup' and p == 0:
                # special handling near 0: U decreasing (all left), V decreasing on [0, q] if q < alpha^2
                assert side == 'L' and k == rc.J and q < ALPHA * ALPHA
                ub = (rc.f_point(fmpq(0)) + rc.V0() - rc.V(arb(q))).upper()  # 2U(0) - V(q)
                lb = None
            else:
                lb, ub = rc.bounds_piece(p, q, l, r, side, k)
            val = ub if mode == 'sup' else lb
            ok = (val < target) if mode == 'sup' else (val > target)
            if ok:
                npieces += 1
                if worst is None or (mode == 'sup' and val > worst) or (mode == 'inf' and val < worst):
                    worst = val
            else:
                if dep >= maxdepth:
                    log("FAIL at piece [%s,%s] bound %s" % (float(p), float(q), val))
                    return False, npieces, val
                m = (p + q) / 2
                stack.append((p, m, dep + 1))
                stack.append((m, q, dep + 1))
        if gi % 400 == 0:
            log("  gap %d/%d  pieces %d  evals %d  worst %s  [%.0fs]" % (gi, len(pts) - 1, npieces, rc.nevals,
                                                                       worst, time.time() - t0))
    return True, npieces, worst
