"""16-mu5: rigorous certificate (acb/arb ball arithmetic) for the CONTOUR form of the norm bound (Theorem 2 of
notes/16-mu5.md):

    limsup (1/K) log ||B_K||  <=  max_{tau in Gamma} H(tau)  -  inf_{t in S} (2U^rho(t) - V(t)),
    H(tau) = 2 U^rho(tau) - Vt(tau),   Vt(tau) = int_0^1 log|tau+u^2| du - 6 int_0^al log|tau+u^2| du,

for any closed contour Gamma, symmetric under conjugation, that encloses the pole segment [-1, -al^2] and stays at
positive distance from [-1, 0] and from supp rho.  Here Gamma^+ (upper half) is the polygon
    P0 = x1  ->  P1 = x1 + iY  ->  P2 = x2 + iY  ->  P3 = x2,       x1 < -1 (saddle of H on the real axis), x2 = 1,
and Gamma is Gamma^+ followed by its mirror image.

On a straight piece with midpoint c and half-length l (so the piece lies in the disc |tau - c| <= l), with
F = 2 Phi_rho - Vt_c analytic (H = Re F):
    H(tau) <= H(c) + |F'(c)| l + (1/2) l^2 max_disc |F''|,
    F'(tau)  = 2 sum_j c_j / (sqrt(tau-a_j) sqrt(tau-b_j)) - [A(1,tau) - 6 A(al,tau)],  A(a,tau) = atan(a/sqrt tau)/sqrt tau,
    |F''|   <= sum_j c_j [d_a^-3/2 d_b^-1/2 + d_a^-1/2 d_b^-3/2] + (1+6al)/d_0^2   on the disc,
    d_a = |c-a_j| - l, d_b = |c-b_j| - l, d_0 = dist(c,[-1,0]) - l.
usage: python contour_certify.py <measure.exact.json> <x1> <Y> <M1'>
"""
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import arb, acb, fmpq, ctx
from certlib import Measure, A, AL, LAM, alA, lamA, PI

ctx.prec = 110


def C(x, y):
    return acb(A(x), A(y))


class CMeasure:
    def __init__(self, meas):
        self.m = meas
        self.rA = [A((b - a) / 2) for a, b in zip(meas.a, meas.b)]
        self.mA = meas.mA
        self.cA = meas.cA
        self.kap = meas.kap          # log(L/4) = log(r/2)
        self.af = [float(x) for x in meas.a]
        self.bf = [float(x) for x in meas.b]
        self.cf = [float(x) for x in meas.c]

    def U(self, tau):
        """U^rho(tau) for complex tau (acb) off supp rho"""
        v = arb(0)
        for j in range(self.m.J):
            z = (tau - self.mA[j]) / self.rA[j]
            w = z + (z - 1).sqrt() * (z + 1).sqrt()
            v += self.cA[j] * (self.kap[j] + abs(w).log())
        return v

    def dPhi(self, tau):
        v = acb(0)
        for j in range(self.m.J):
            v += self.cA[j] / ((tau - self.m.aA[j]).sqrt() * (tau - self.m.bA[j]).sqrt())
        return v

    def d2Phi_bound(self, cx, cy, l):
        """upper bound for |Phi''| on the disc |tau - c| <= l (arb)"""
        v = arb(0)
        lA = A(l)
        for j in range(self.m.J):
            da = abs(C(cx - self.m.a[j], cy)) - lA
            db = abs(C(cx - self.m.b[j], cy)) - lA
            assert da > 0 and db > 0
            v += self.cA[j] / 2 * (1 / (da * da.sqrt() * db.sqrt()) + 1 / (da.sqrt() * db * db.sqrt()))
        return v


def I_log(a, c):
    """int_0^a log|u - c| du = Re[(a-c) Log(a-c) + c Log(-c)] - a   (c complex, Im c != 0)"""
    aa = acb(a)
    return ((aa - c) * (aa - c).log() + c * (-c).log()).real - a


def Vt(tau):
    st = tau.sqrt()
    ci = acb(0, 1) * st
    J = lambda a: I_log(a, ci) + I_log(a, -ci)
    return J(arb(1)) - 6 * J(alA)


def dVt(tau):
    st = tau.sqrt()
    Af = lambda a: (acb(a) / st).atan() / st
    return Af(arb(1)) - 6 * Af(alA)


def dist01(x, y):
    """dist(x+iy, [-1,0]) as arb"""
    if -1 <= x <= 0:
        return A(abs(y))
    if x < -1:
        return abs(C(x + 1, y))
    return abs(C(x, y))


def certify_path(cm, x1, Y, x2, M1p, verbose=True, maxpieces=200000):
    t0 = time.time()
    M1A = A(M1p)
    segs = [((x1, Fr(0)), (x1, Y)), ((x1, Y), (x2, Y)), ((x2, Y), (x2, Fr(0)))]
    stack = []
    for (xa, ya), (xb, yb) in segs:
        n0 = 64
        for k in range(n0):
            p = (xa + (xb - xa) * Fr(k, n0), ya + (yb - ya) * Fr(k, n0))
            q = (xa + (xb - xa) * Fr(k + 1, n0), ya + (yb - ya) * Fr(k + 1, n0))
            stack.append((p, q))
    nacc = nsplit = 0
    worst = None
    supp_lo, supp_hi = cm.m.a[-1], cm.m.b[-1]
    while stack:
        p, q = stack.pop()
        cx, cy = (p[0] + q[0]) / 2, (p[1] + q[1]) / 2
        l = max(abs(q[0] - p[0]), abs(q[1] - p[1])) / 2      # axis-parallel pieces
        # distance requirements (exact, conservative)
        d0 = dist01(cx, cy)
        dsupp = abs(C(cx - supp_lo, cy)) if cx < supp_lo else (abs(C(cx - supp_hi, cy)) if cx > supp_hi else A(abs(cy)))
        if not (d0 > 2 * A(l) and dsupp > 2 * A(l)):
            nsplit += 1
            m = ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
            stack.append((p, m)); stack.append((m, q))
            continue
        tau = C(cx, cy)
        H = 2 * cm.U(tau) - Vt(tau)
        dF = 2 * cm.dPhi(tau) - dVt(tau)
        d0l = d0 - A(l)
        D2 = 2 * cm.d2Phi_bound(cx, cy, l) + (1 + 6 * alA) / (d0l * d0l)
        UB = H + abs(dF) * A(l) + D2 * A(l) * A(l) / 2
        if UB < M1A:
            nacc += 1
            if worst is None or UB.upper() > worst:
                worst = UB.upper()
            continue
        nsplit += 1
        if nsplit > maxpieces or l < Fr(1, 10 ** 25):
            raise RuntimeError(f"contour certificate FAILED near {float(cx)}+{float(cy)}i (l={float(l):.3g}): UB {UB} vs {float(M1p)}")
        m = ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
        stack.append((p, m)); stack.append((m, q))
    if verbose:
        print(f"   max_Gamma H < M1' = {float(M1p):.10f} certified on Gamma+ = [{float(x1)}] -> [{float(x1)}+{float(Y)}i] -> "
              f"[{float(x2)}+{float(Y)}i] -> [{float(x2)}]: {nacc} pieces, {nsplit} splits, worst {worst}  [{time.time()-t0:.1f}s]")
    return dict(n=nacc, splits=nsplit, worst=worst)


if __name__ == "__main__":
    fn = sys.argv[1]
    x1 = Fr(sys.argv[2]); Y = Fr(sys.argv[3]); M1p = Fr(sys.argv[4])
    x2 = Fr(1)
    meas = Measure.from_exact_json(fn)
    meas.check()
    cm = CMeasure(meas)
    # value at the saddle point (real formula, cross-check of the complex one)
    s1 = None
    Hs_real = 2 * meas.U(x1)
    from certlib import G_arb
    # Vt(-s^2) for s > 1: G(s+1) - G(s-1) - 6[G(s+al) - G(s-al)]
    sA = A(-x1).sqrt()
    Gb = lambda z: z * z.log() - z
    Vt_real = Gb(sA + 1) - Gb(sA - 1) - 6 * (Gb(sA + alA) - Gb(sA - alA))
    Hs_real = Hs_real - Vt_real
    Hs_cplx = 2 * cm.U(C(x1, 0)) - Vt(C(x1, 0))
    print(f"   H(x1) real formula = {Hs_real.str(15)};  complex formula = {Hs_cplx.str(15)}")
    certify_path(cm, x1, Y, x2, M1p)
