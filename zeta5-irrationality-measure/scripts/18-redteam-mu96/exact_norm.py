"""18-redteam-mu96: exact ||B_K|| for small K and finite-K checks of 16's Lemmas 3, 5, Cor 3' and Theorem 2.

||B_K|| = 1/min_j |X_j - zeta(5)| over the roots X_j of Delta_K(X) = det(G0 + X M) (all real, B_K symmetric).
Delta_K(X) = det(M) * det(X I + M^{-1} G0) = det(M) * charpoly(-M^{-1}G0)(X).
Extremal direction q*: exact solve (G0 + Xt M) v = e_1 at a rational Xt within 10^-(prec) of the nearest root.
"""
import sys, time, math
import mpmath as mp
from flint import fmpq, fmpz, fmpq_mat, fmpq_poly, arb, ctx
sys.path.insert(0, './scripts/18-redteam-mu96')
from exact_small_K import build, w_mp
from common import load_measure, FloatMeasure, V_float, Vt_float
import cypari2

pari = cypari2.Pari()
pari.allocatemem(2 * 10 ** 9)


def tompf(x):
    return mp.mpf(int(x.p)) / int(x.q)


def exact_norm(K, digits=None):
    t0 = time.time()
    D = build(K)
    N, h, G0, M = D['N'], D['h'], D['G0'], D['M']
    if digits is None:
        digits = int(8 * K) + 200
    Cm = M.solve(G0)  # M^{-1} G0
    negC = fmpq_mat(h, h, [-x for x in Cm.entries()])
    P = negC.charpoly()  # det(X I + C)
    print("K=%d: charpoly degree %d [%.1fs]" % (K, P.degree(), time.time() - t0), flush=True)
    pari.set_real_precision(digits)
    coeffs = [str(c) for c in P.coeffs()]
    pol = pari('Pol([' + ','.join(reversed(coeffs)) + '])')
    roots = pari.polrootsreal(pol, precision=int(digits * 3.33) + 64)
    nroots = len(roots)
    mp.mp.dps = digits
    z5 = mp.zeta(5)
    rts = [mp.mpf(str(r)) for r in roots]
    dist = sorted((abs(r - z5), r) for r in rts)
    dmin, Xstar = dist[0]
    print("K=%d: %d real roots (deg %d); nearest root distance %s ; second %s [%.1fs]" % (
        K, nroots, P.degree(), mp.nstr(dmin, 10), mp.nstr(dist[1][0], 10), time.time() - t0), flush=True)
    lognorm = -mp.log(dmin)
    return D, lognorm, Xstar, z5, dist


def extremal_vector(D, Xstar, z5, extra=60):
    """exact rational solve near the nearest root; returns rational q (list of fmpq) and Rayleigh quotient."""
    h, G0, M = D['h'], D['G0'], D['M']
    dps = mp.mp.dps
    # rational approximation of Xstar to dps-extra digits
    s = mp.nstr(Xstar, dps - extra, strip_zeros=False)
    Xt = fmpq(*[int(x) for x in mp.identify(0) and (0, 1)]) if False else None
    sign = -1 if s.startswith('-') else 1
    s = s.lstrip('-')
    ip, fp = s.split('.')
    Xt = fmpq(sign * int(ip + fp), 10 ** len(fp))
    A = fmpq_mat(h, h, [G0[i, k] + Xt * M[i, k] for i in range(h) for k in range(h)])
    e1 = fmpq_mat(h, 1, [fmpq(1)] + [fmpq(0)] * (h - 1))
    v = A.solve(e1)
    q = [v[i, 0] for i in range(h)]
    # normalise by max |coef|
    mx = max(abs(x) for x in q)
    q = [x / mx for x in q]
    # Rayleigh quotient |q^T M q| / q^T G q
    qM = sum((M[i, k] * q[i] * q[k] for i in range(h) for k in range(h)), fmpq(0))
    qG0 = sum((G0[i, k] * q[i] * q[k] for i in range(h) for k in range(h)), fmpq(0))
    qMf, qG0f = tompf(qM), tompf(qG0)
    qG = qG0f + z5 * qMf
    return q, qMf, qG


def finite_K_checks(K, D, q, qMf, qG, lognormB, consts):
    """Check 16's Lemma 3 / Cor 3' / Lemma 5 / Theorem 2 at finite K for the extremal p(tau) = q(K^2 tau)."""
    N, h = D['N'], D['h']
    Kf = mp.mpf(K)
    qpoly = [tompf(c) for c in q]
    def p(tau):
        r = mp.mpc(0) if isinstance(tau, mp.mpc) else mp.mpf(0)
        x = Kf ** 2 * tau
        for c in reversed(qpoly):
            r = r * x + c
        return r
    intW = qG / Kf ** (12 * N - 2 * K + 1)  # int p^2 W_K
    print("  check: log(K^5 |q^T M q|/K^(12N-2K+6) ... ) Rayleigh: log|qMq|-log qGq = %s vs log||B|| = %s" % (
        mp.nstr(mp.log(abs(qMf)) - mp.log(qG), 15), mp.nstr(lognormB, 15)))
    A, B, C = load_measure()
    fm = FloatMeasure(A, B, C)
    m0, M1c = consts['m0'], consts['M1c']
    tS = float(A[-1])
    c0 = math.log((2 * math.pi) ** 4 / 24) - 2
    EK = 3 * math.pi / (2 * K * math.sqrt(tS)) - 3 * math.log(K) - 2 * math.log(tS) - c0
    # Lemma 3 at points t in S: p(t)^2 <= h^2 K^2 exp(K V(t) + 3pi/(2K sqrt t) - 3 log K - 2 log t - c0) int p^2 W
    worst3 = -1e300
    for t in [float(A[-1]) * (1 + 1e-9)] + [0.0001 * 1.25 ** k for k in range(40) if 0.0001 * 1.25 ** k < float(B[-1])]:
        lhs = 2 * mp.log(abs(p(mp.mpf(t))))
        rhs = math.log(h * h * K * K) + K * V_float(t) + 3 * math.pi / (2 * K * math.sqrt(t)) - 3 * math.log(K) \
            - 2 * math.log(t) - c0 + mp.log(intW)
        worst3 = max(worst3, float(lhs - rhs))
    print("  Lemma 3 (Nikolskii) at 40 points of S: max(log lhs - log rhs) = %.3f (must be <= 0)" % worst3)
    # Cor 3' and Lemma 5 on Gamma
    x1 = -1.069911
    pts = [complex(x1, 0.2 * k / 20) for k in range(21)] + [complex(x1 + (1 - x1) * k / 40, 0.2) for k in range(41)] \
        + [complex(1.0, 0.2 * k / 20) for k in range(21)]
    worstC, worst5, worstT = -1e300, -1e300, -1e300
    u = [mp.mpf(j) / K for j in range(K + 1)]
    for z in pts:
        zm = mp.mpc(z.real, z.imag)
        lp2 = 2 * mp.log(abs(p(zm)))
        Uz = fm.Uc(z)
        rhs = math.log(h * h * K * K) + EK + K * (2 * Uz - m0) + mp.log(intW)
        worstC = max(worstC, float(lp2 - rhs))
        # Lemma 5: log|Omega| <= -K Vt + 7/d
        lo = mp.mpf(0)
        for l in range(1, N + 1):
            lo += 6 * mp.log(abs(zm + u[l] ** 2))
        for j in range(1, K + 1):
            lo -= mp.log(abs(zm + u[j] ** 2))
        d = min(abs(z - complex(max(-1, min(0, z.real)), 0)), 10)
        worst5 = max(worst5, float(lo - (-K * Vt_float(z) + 7 / d)))
        # integrand vs Theorem-2 pointwise bound: |tau^2 Omega p^2| <= |tau|^2 h^2K^2 e^{E+7/d} e^{K(H - m0)} int p^2W
        H = 2 * Uz - Vt_float(z)
        worstT = max(worstT, float(lp2 + lo - (math.log(h * h * K * K) + EK + 7 / d + K * (H - m0) + mp.log(intW))))
    print("  Cor 3' on Gamma (83 pts): max(log|p|^2 - log rhs) = %.3f (must be <= 0)" % worstC)
    print("  Lemma 5 on Gamma (83 pts): max(log|Omega| - (-K Vt + 7/d)) = %.3f (must be <= 0)" % worst5)
    print("  Theorem 2 integrand bound on Gamma: max(log|Omega p^2| - log bound) = %.3f (must be <= 0)" % worstT)
    # full Theorem-2 bound
    lenG, maxtau2, dG = 2 * (0.4 + 1 - x1), (x1 * x1 + 0.04), 1 + x1 if False else (-1 - x1)
    thm2 = K * (M1c - m0) + 7 * math.log(K) + 2 * math.log(h) + EK + math.log(lenG * maxtau2 / (2 * math.pi)) + 7 / dG
    print("  Theorem 2 explicit bound: log||B_K|| <= %.3f ; true log||B_K|| = %.6f ; rate true %.6f" % (
        thm2, float(lognormB), float(lognormB) / K))
    return dict(worst3=worst3, worstC=worstC, worst5=worst5, worstT=worstT, thm2=thm2)


if __name__ == "__main__":
    consts = dict(m0=-6.660814, M1c=0.808802)
    for K in [int(a) for a in sys.argv[1:]] or [40]:
        D, lognorm, Xstar, z5, dist = exact_norm(K)
        print("K=%d: log||B_K|| = %s ; (1/K) log||B_K|| = %s" % (K, mp.nstr(lognorm, 16), mp.nstr(lognorm / K, 10)),
              flush=True)
        q, qMf, qG = extremal_vector(D, Xstar, z5)
        finite_K_checks(K, D, q, qMf, qG, lognorm, consts)
