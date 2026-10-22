"""16-mu5: rigorous certificates (arb) for a nested arcsine measure rho:

  (a) sup_{t >= 0} (2U^rho(t) - V(t)) < M0        -> Lemma 6.1 constant  U = lam*M0 - I(rho) + C*   (Fauzan (6.2)-(6.4))
  (b) inf_{t in S} (2U^rho(t) - V(t)) > m0         (S = supp rho = [a_last, b_last])
  (c) sup_{s in [al,1]} (2U^rho(-s^2) + phi(s)) < M1
      -> norm-rate constant  rho_K = M1 - m0  (notes/16-mu5.md, Theorem 1)

Interval bounds used on an interval [l,r] with no breakpoint a_j, b_j in its interior (h = r-l, c = (l+r)/2):
  sup:  min( 2 max(U(l),U(r)) - minV(l,r),                                  [quasi-convexity of U, unimodality of V]
             f(c) + |f'(c)| h/2 + max(0, -min V''[l,r]) h^2/8 )              [U concave (tangent), Taylor for V]
  inf:  max( 2 min(U(l),U(r)) - max(V(l),V(r)),                              [U concave => >= min of endpoints]
             U(l)+U(r)-V(c) - |2(U(r)-U(l))/h - V'(c)| h/2 - max(0, max V''[l,r]) h^2/8 )   [chord; Taylor for V]
  neg:  min( 2U(-s_r^2) + max_{[s_l,s_r]} phi (concavity tangents),
             F(c) + |F'(c)| h/2 + (4 lam / s_l^2) h^2/8 )                       [F = 2U(-s^2)+phi, F'' <= 4 lam/s^2]
usage: python certify.py <measure.json> <mode: U|rho|both> [M0] [m0] [M1]
"""
import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from flint import arb
from certlib import (Measure, A, V_arb, dV_arb, d2V_arb, ball, phi_arb, dphi_arb, v_min_bracket,
                     LAM, AL, lamA, alA, PI, CSTAR)


def certify_sup(meas, M0, maxint=4_000_000, verbose=True):
    """sup over [0,2] of 2U - V < M0 (rational M0); t >= 2 via the tail bound g(t)."""
    t0 = time.time()
    M0A = A(M0)
    qlo, qhi, Vmin_lb = v_min_bracket()

    def Vmin(l, r):
        if r <= qlo:
            return V_arb(r).lower()
        if l >= qhi:
            return V_arb(l).lower()
        return Vmin_lb

    pts = sorted(set([Fr(0), Fr(2), qlo, qhi] + [x for x in meas.breaks if 0 < x < 2]))
    stack = list(zip(pts[:-1], pts[1:]))
    nacc = nsplit = 0
    worst = None
    while stack:
        l, r = stack.pop()
        Ul, Ur = meas.U(l), meas.U(r)
        um = Ul if Ul.upper() > Ur.upper() else Ur
        B0 = 2 * um - Vmin(l, r)
        best = B0
        if not (B0 < M0A) and l > 0:
            c = (l + r) / 2
            h = A(r - l)
            fc = 2 * meas.U(c) - V_arb(c)
            dfc = 2 * meas.dU(c) - dV_arb(c)
            v2 = d2V_arb(ball(l, r)).lower()
            curv = -v2 if v2 < 0 else arb(0)
            B2 = fc + abs(dfc) * h / 2 + curv * h * h / 8
            if B2.upper() < best.upper():
                best = B2
        if best < M0A:
            nacc += 1
            if worst is None or best.upper() > worst:
                worst = best.upper()
            continue
        nsplit += 1
        if nsplit > maxint or r - l < Fr(1, 10 ** 30):
            raise RuntimeError(f"sup certificate FAILED near [{float(l)}, {float(r)}]: bound {best} vs M0 {float(M0)}")
        m = (l + r) / 2
        stack.append((l, m))
        stack.append((m, r))
    # t >= 2:  2U - V <= (2lam - 1 + 6al) log t + 2 al^3/t - 2 pi sqrt t  (decreasing on [2,oo))
    t2 = arb(2)
    g2 = (2 * lamA - 1 + 6 * alA) * t2.log() + 2 * alA ** 3 / t2 - 2 * PI * t2.sqrt()
    assert g2 < M0A and (2 * lamA - 1 + 6 * alA) / arb(2).sqrt() < PI
    # Fauzan's (6.8) form, as used verbatim in (6.7): 13/10 log t + 6 al^3/t - 2 pi sqrt t at t = 2 (plus sqrt2/K) < M0
    g2F = arb(13) / 10 * t2.log() + 6 * alA ** 3 / t2 - 2 * PI * t2.sqrt()
    assert g2F < M0A
    if verbose:
        print(f"   sup_(t>=0)(2U-V) < M0 = {float(M0):.10f} certified: {nacc} intervals, {nsplit} splits, "
              f"worst bound {worst}, tail g(2) = {g2.str(6)}, (6.8) at 2 = {g2F.str(6)}  [{time.time()-t0:.1f}s]")
    return dict(n=nacc, splits=nsplit, worst=worst)


def certify_inf(meas, m0, maxint=4_000_000, verbose=True):
    """inf over S = [a_last, b_last] of 2U - V > m0 (rational m0)."""
    t0 = time.time()
    m0A = A(m0)
    S0, S1 = meas.a[-1], meas.b[-1]
    pts = sorted(set([x for x in meas.breaks if S0 <= x <= S1]))
    stack = list(zip(pts[:-1], pts[1:]))
    nacc = nsplit = 0
    worst = None
    while stack:
        l, r = stack.pop()
        Ul, Ur = meas.U(l), meas.U(r)
        Vl, Vr = V_arb(l), V_arb(r)
        umin = Ul if Ul.lower() < Ur.lower() else Ur
        vmax = Vl if Vl.upper() > Vr.upper() else Vr
        B0 = 2 * umin - vmax
        best = B0
        if not (B0 > m0A):
            c = (l + r) / 2
            h = A(r - l)
            slope = (Ur - Ul) / h
            v2 = d2V_arb(ball(l, r)).upper()
            curv = v2 if v2 > 0 else arb(0)
            B2 = Ul + Ur - V_arb(c) - abs(2 * slope - dV_arb(c)) * h / 2 - curv * h * h / 8
            if B2.lower() > best.lower():
                best = B2
        if best > m0A:
            nacc += 1
            if worst is None or best.lower() < worst:
                worst = best.lower()
            continue
        nsplit += 1
        if nsplit > maxint or r - l < Fr(1, 10 ** 30):
            raise RuntimeError(f"inf certificate FAILED near [{float(l)}, {float(r)}]: bound {best} vs m0 {float(m0)}")
        m = (l + r) / 2
        stack.append((l, m))
        stack.append((m, r))
    if verbose:
        print(f"   inf_S(2U-V) > m0 = {float(m0):.10f} certified on S = [{float(S0):.6g}, {float(S1):.8f}]: {nacc} intervals,"
              f" {nsplit} splits, worst bound {worst}  [{time.time()-t0:.1f}s]")
    return dict(n=nacc, splits=nsplit, worst=worst)


def certify_neg(meas, M1, maxint=1_000_000, verbose=True):
    """sup over s in [al, 1] of F(s) = 2U(-s^2) + phi(s) < M1."""
    t0 = time.time()
    M1A = A(M1)
    lam4 = 4 * lamA
    nacc = nsplit = 0
    worst = None
    N0 = 64
    pts = [AL + (1 - AL) * Fr(k, N0) for k in range(N0 + 1)]
    stack = list(zip(pts[:-1], pts[1:]))

    def U_neg(s):
        return meas.U(-s * s)

    while stack:
        l, r = stack.pop()
        h = A(r - l)
        # zeroth order: U(-s^2) nondecreasing in s; phi concave
        cands = []
        if l > AL:
            dl = dphi_arb(l)
            cands.append(phi_arb(l) + (dl if dl > 0 else arb(0)) * h)
        if r < 1:
            dr = dphi_arb(r)
            cands.append(phi_arb(r) + (-dr if dr < 0 else arb(0)) * h)
        phimax = min(cands, key=lambda x: x.upper())
        B0 = 2 * U_neg(r) + phimax
        best = B0
        if not (B0 < M1A):
            c = (l + r) / 2
            Fc = 2 * U_neg(c) + phi_arb(c)
            dFc = 4 * A(c) * meas.dU_neg_upper_abs(-c * c) + dphi_arb(c)
            B2 = Fc + abs(dFc) * h / 2 + lam4 / A(l * l) * h * h / 8
            if B2.upper() < best.upper():
                best = B2
        if best < M1A:
            nacc += 1
            if worst is None or best.upper() > worst:
                worst = best.upper()
            continue
        nsplit += 1
        if nsplit > maxint or r - l < Fr(1, 10 ** 30):
            raise RuntimeError(f"neg certificate FAILED near [{float(l)}, {float(r)}]: bound {best} vs M1 {float(M1)}")
        m = (l + r) / 2
        stack.append((l, m))
        stack.append((m, r))
    if verbose:
        print(f"   sup_[al,1](2U(-s^2)+phi) < M1 = {float(M1):.10f} certified: {nacc} intervals, {nsplit} splits,"
              f" worst {worst}  [{time.time()-t0:.1f}s]")
    return dict(n=nacc, splits=nsplit, worst=worst)


def U_constant(meas, M0):
    I = meas.energy()
    Ub = lamA * A(M0) - I + CSTAR
    return I, Ub


if __name__ == "__main__":
    fn = sys.argv[1]
    mode = sys.argv[2]
    if fn.endswith('.exact.json'):
        meas = Measure.from_exact_json(fn)
    elif 'arcsine_opt' in fn:
        meas = Measure.from_json(fn, key='3p', denom=10 ** 12)
    else:
        meas = Measure.from_json(fn)
    ok = meas.check()
    assert ok['positive'] and ok['mass_exact'] and ok['nested'] and ok['support_in_(0,2)']
    out = {}
    if mode in ('U', 'both'):
        M0 = Fr(sys.argv[3])
        certify_sup(meas, M0)
        assert ok['reg_const'] <= 60
        I, Ub = U_constant(meas, M0)
        print(f"   I(rho) = {I.str(20)};  C* = {CSTAR.str(20)}")
        print(f"   ==> U <= lam*M0 - I + C* = {Ub.str(20)}  (upper {Ub.upper().str(15)})")
        out['U_upper'] = Ub.upper().str(20)
    if mode in ('rho', 'both'):
        k = 4 if mode == 'both' else 3
        m0, M1 = Fr(sys.argv[k]), Fr(sys.argv[k + 1])
        certify_inf(meas, m0)
        certify_neg(meas, M1)
        print(f"   ==> rho_K <= M1 - m0 = {M1 - m0} = {float(M1 - m0):.10f}")
        out['rho'] = str(M1 - m0)
    print(json.dumps(out))
