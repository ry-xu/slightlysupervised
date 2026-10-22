"""18-redteam-mu96: rigorous (acb) bound for H(tau) = 2U^rho(tau) - Vt(tau) on 16's contour Gamma+.

Gamma+ : x1 -> x1 + i/5 -> 1 + i/5 -> 1   (x1 = -1069911/10^6); Gamma = Gamma+ u conj(Gamma+); H(conj z) = H(z).

Own scheme: on a sub-segment with centre c and half-length l,
     H(tau) <= H(c) + l * sup_{segment} |F'|,      H = Re F locally, F' = 2 sum_j c_j/P_j(w_j) - Vt'(tau),
with sup |F'| enclosed by acb evaluation on a box containing the segment (mean value theorem).

Branches (all analytic on the boxes used; asserted below):
 * arcsine potential: U_j(tau) = log|wt + P(wt)| - log 2, w = tau - m_j, P(w) = sqrt(w-r)sqrt(w+r) (principal),
   wt = -w on the left edge (Re w < 0 certainly; P odd so 1/P(w) = -1/P(-w)), wt = w on top/right edges
   (top: Im w = 1/5 != 0; right: Re(w - r_j) = 1 - b_j > 0).
 * J(a,tau) = int_0^a log|tau+u^2| du = Re[(s+a)Log(s+a) - (s-a)Log(s-a)] - 2a, s^2 = -tau,
   dJ/dtau = (Log(s-a) - Log(s+a))/(2s);  s = sqrt(-tau) on the left/top edges, s = i sqrt(tau) on the right edge.
   Validity: s +- u (0<=u<=a<=1) must avoid (-inf,0]: left edge Re s > 1; top edge Im s != 0; right edge Im s > 0.
"""
import sys, time
from flint import arb, acb, fmpq, ctx
sys.path.insert(0, './scripts/18-redteam-mu96')
from common import load_measure, ALPHA

ctx.prec = 128


class BranchError(Exception):
    pass


class ContourCert:
    def __init__(self, A, B, C):
        self.J = len(A)
        self.C = [arb(x) for x in C]
        self.M = [(arb(A[j]) + arb(B[j])) / 2 for j in range(self.J)]
        self.R = [(arb(B[j]) - arb(A[j])) / 2 for j in range(self.J)]
        self.al = arb(ALPHA)
        self.log2 = arb(2).log()
        self.nevals = 0

    def s_of(self, tau, edge):
        """Root s of s^2 = -tau on an analytic branch; raises BranchError if the validity condition
        (see module docstring) cannot be certified on this box (caller then subdivides)."""
        if edge == 'right':
            s = acb(0, 1) * tau.sqrt()
            ok = s.imag > 0
        else:
            s = (-tau).sqrt()
            if edge == 'left':
                ok = s.real > 1
            else:
                ok = (s.imag < 0) or (s.imag > 0)
        if not ok:
            raise BranchError(edge)
        return s

    def Jval(self, a, s):
        a = acb(a)
        return ((s + a) * (s + a).log() - (s - a) * (s - a).log()).real - 2 * a.real

    def dJ(self, a, s):
        a = acb(a)
        return ((s - a).log() - (s + a).log()) / (2 * s)

    def Hc(self, tau, edge):
        """Rigorous enclosure of H(tau) for a box tau (uses Re Log, which is an exact enclosure of log|.|)."""
        u = arb(0)
        for j in range(self.J):
            w = tau - self.M[j]
            wt = -w if edge == 'left' else w
            if edge == 'left':
                ok = wt.real > self.R[j]
            elif edge == 'right':
                ok = (wt - self.R[j]).real > 0
            else:
                ok = wt.imag > 0
            if not ok:
                raise BranchError(edge)
            P = (wt - self.R[j]).sqrt() * (wt + self.R[j]).sqrt()
            u += self.C[j] * ((wt + P).log().real - self.log2)
        s = self.s_of(tau, edge)
        Vt = self.Jval(1, s) - 6 * self.Jval(self.al, s)
        return 2 * u - Vt

    def dF(self, tau, edge):
        """Enclosure of F'(tau) = 2 sum c_j / P(w_j) - Vt'(tau) on a box."""
        self.nevals += 1
        d = acb(0)
        for j in range(self.J):
            w = tau - self.M[j]
            if edge == 'left':
                wt = -w
                if not (wt.real > self.R[j]):
                    raise BranchError(edge)
                P = (wt - self.R[j]).sqrt() * (wt + self.R[j]).sqrt()
                d -= self.C[j] / P
            else:
                if edge == 'top' and not (w.imag > 0):
                    raise BranchError(edge)
                if edge == 'right' and not ((w - self.R[j]).real > 0):
                    raise BranchError(edge)
                P = (w - self.R[j]).sqrt() * (w + self.R[j]).sqrt()
                d += self.C[j] / P
        s = self.s_of(tau, edge)
        dVt = self.dJ(1, s) - 6 * self.dJ(self.al, s)
        return 2 * d - dVt


X1 = fmpq(-1069911, 10 ** 6)
Y = fmpq(1, 5)


def edge_point(edge, u):
    """parametrisation by a rational u in [0,1]."""
    if edge == 'left':
        return X1, Y * u           # x1 + i*(u/5)
    if edge == 'top':
        return X1 + (1 - X1) * u, Y
    if edge == 'right':
        return fmpq(1), Y * u


def certify_contour(cc, target, log=print, maxdepth=50):
    worst = None
    npieces = 0
    lens = {'left': Y, 'top': 1 - X1, 'right': Y}
    for edge in ['left', 'top', 'right']:
        stack = [(fmpq(0), fmpq(1), 0)]
        while stack:
            u0, u1, dep = stack.pop()
            x0, y0 = edge_point(edge, u0)
            xx1, yy1 = edge_point(edge, u1)
            cx, cy = (x0 + xx1) / 2, (y0 + yy1) / 2
            c = acb(arb(cx), arb(cy))
            box = acb(arb(x0).union(arb(xx1)), arb(y0).union(arb(yy1)))
            ell = arb(lens[edge] * (u1 - u0)) / 2
            try:
                Hcv = cc.Hc(c, edge)
                dfa = cc.dF(box, edge).abs_upper()
                ub = (Hcv + ell * dfa).upper()
                certified = ub < target
            except BranchError:
                ub, certified = None, False
            if certified:
                npieces += 1
                if worst is None or ub > worst:
                    worst = ub
            else:
                if dep >= maxdepth:
                    log("FAIL on %s [%s,%s]: %s" % (edge, float(u0), float(u1), ub))
                    return False, npieces, ub
                um = (u0 + u1) / 2
                stack.append((u0, um, dep + 1))
                stack.append((um, u1, dep + 1))
        log("  edge %s done: pieces so far %d, worst %s" % (edge, npieces, worst))
    return True, npieces, worst
