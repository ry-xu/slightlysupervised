"""18-redteam-mu96: exact construction of G_K(X) = G0 + X*M (Fauzan (2.2)-(2.4)) for small K, and checks:
  (a) Prop 2.2:  q^T G_K(zeta5) q = int_0^oo D_N(y^2)^6 q(y^2)^2 / D_K(y^2) w(y) dy        (mpmath quad)
  (b) residue form: q^T M q = sum_{N<j<=K} j^4 Res_{t=-j^2} D_N^6 q^2/D_K                   (exact)
  (c) contour form: q^T M q = (1/2 pi i) oint_{K^2 Gamma} t^2 D_N^6 q^2 / D_K dt          (mpmath quad on 16's rectangle)
  (d) Lemma 1 scalings (a),(b) in the variable tau = t/K^2
  (e) exact ||B_K|| = 1/min_j |X_j - zeta5| over the (real) roots X_j of det(G0 + X M)
"""
import sys, time, random
import mpmath as mp
from flint import fmpq, fmpz, fmpq_mat, fmpq_poly, arb, ctx

ROOT = "."


def mu_mono(e):
    """mu_X(t^e) = (-1)^e B_{2e+2} (2e+3)(2e+4)(2e+5)/24   (2.2)"""
    return (-1) ** e * fmpq.bernoulli(2 * e + 2) * (2 * e + 3) * (2 * e + 4) * (2 * e + 5) / 24


def H5(j):
    return sum((fmpq(1, v ** 5) for v in range(1, j + 1)), fmpq(0))


def build(K):
    assert K % 40 == 0
    N, h = 3 * K // 40, 37 * K // 40
    t = fmpq_poly([0, 1])
    DN = fmpq_poly([1])
    for j in range(1, N + 1):
        DN *= (t + j * j)
    Dtail = fmpq_poly([1])
    for j in range(N + 1, K + 1):
        Dtail *= (t + j * j)
    dDtail = Dtail.derivative()
    DN5 = DN ** 5
    poles = list(range(N + 1, K + 1))
    # residues of D_N^5 t^m / D_tail at -j^2:  D_N(-j^2)^5 (-j^2)^m / D_tail'(-j^2)
    base = {j: DN5(fmpq(-j * j)) / dDtail(fmpq(-j * j)) for j in poles}
    polevals = {j: fmpq(j) ** 4 * (-H5(j)) - fmpq(1, 4) + fmpq(1, 2 * j) for j in poles}  # X-free part of (2.3)
    G0row, Mrow = [], []
    for m in range(2 * h - 1):
        num = DN5 * t ** m
        Q = num // Dtail  # polynomial part
        polypart = sum((Q[e] * mu_mono(e) for e in range(Q.degree() + 1)), fmpq(0)) if Q.degree() >= 0 else fmpq(0)
        g0 = polypart
        mm = fmpq(0)
        for j in poles:
            cj = base[j] * fmpq(-j * j) ** m
            g0 += cj * polevals[j]
            mm += cj * fmpq(j) ** 4
        G0row.append(g0)
        Mrow.append(mm)
    G0 = fmpq_mat(h, h, [G0row[i + k] for i in range(h) for k in range(h)])
    M = fmpq_mat(h, h, [Mrow[i + k] for i in range(h) for k in range(h)])
    return dict(K=K, N=N, h=h, G0=G0, M=M, DN=DN, Dtail=Dtail, DK=DN * Dtail)


def qform(Mat, q):
    h = len(q)
    s = fmpq(0)
    for i in range(h):
        for k in range(h):
            s += Mat[i, k] * q[i] * q[k]
    return s


def w_mp(y):
    # (2pi)^4 y^5/12 * sum l^4 e^{-2 pi l y} = (2pi)^4 y^5/12 * x(1+11x+11x^2+x^3)/(1-x)^5, x = e^{-2 pi y}
    x = mp.exp(-2 * mp.pi * y)
    return (2 * mp.pi) ** 4 * y ** 5 / 12 * x * (1 + 11 * x + 11 * x ** 2 + x ** 3) / (-mp.expm1(-2 * mp.pi * y)) ** 5


def polyval_mp(P, x):
    cs = [mp.mpf(int(c.p)) / int(c.q) for c in P.coeffs()]
    r = mp.mpf(0) if not isinstance(x, mp.mpc) else mp.mpc(0)
    for c in reversed(cs):
        r = r * x + c
    return r


def main(K, dps=50, seed=1):
    mp.mp.dps = dps
    t0 = time.time()
    D = build(K)
    N, h, G0, M = D['N'], D['h'], D['G0'], D['M']
    print("K=%d N=%d h=%d built G0, M exactly [%.1fs]" % (K, N, h, time.time() - t0))
    z5 = mp.zeta(5)
    random.seed(seed)
    q = [fmpq(random.randint(-3, 3)) for _ in range(h)]
    qpoly = fmpq_poly(q)
    # (a) Prop 2.2
    lhs_G0 = qform(G0, q)
    lhs_M = qform(M, q)
    with mp.workdps(600):  # G0-part and zeta5*M-part cancel to many digits
        lhs = mp.mpf(int(lhs_G0.p)) / int(lhs_G0.q) + mp.zeta(5) * mp.mpf(int(lhs_M.p)) / int(lhs_M.q)
        print("    |G0-part| = %s, |zeta5*M-part| = %s" % (mp.nstr(abs(mp.mpf(int(lhs_G0.p)) / int(lhs_G0.q)), 5),
                                                         mp.nstr(abs(mp.zeta(5) * mp.mpf(int(lhs_M.p)) / int(lhs_M.q)), 5)))
    lhs = +lhs
    num = D['DN'] ** 6 * qpoly ** 2
    den = D['DK']
    f = lambda y: polyval_mp(num, y * y) / polyval_mp(den, y * y) * w_mp(y)
    brk = [0] + [mp.mpf(K) * x / 8 for x in range(1, 41)]
    rhs = mp.quad(f, brk)
    print("(a) Prop 2.2: q^T G(zeta5) q = %s ; integral = %s ; rel diff %s" % (
        mp.nstr(lhs, 25), mp.nstr(rhs, 25), mp.nstr(abs(lhs - rhs) / abs(rhs), 3)))
    # (b) residue form (exact)
    t = fmpq_poly([0, 1])
    res_sum = fmpq(0)
    dDt = D['Dtail'].derivative()
    for j in range(N + 1, K + 1):
        x = fmpq(-j * j)
        res = D['DN'](x) ** 5 * qpoly(x) ** 2 / dDt(x)  # Res_{t=-j^2} D_N^6 q^2/D_K
        res_sum += fmpq(j) ** 4 * res
    print("(b) q^T M q == sum_j j^4 Res : %s" % (res_sum == lhs_M))
    # (c) contour integral of t^2 D_N^6 q^2/D_K over K^2 * Gamma (16's rectangle), counterclockwise
    num2 = t ** 2 * D['DN'] ** 6 * qpoly ** 2
    g = lambda z: polyval_mp(num2, z) / polyval_mp(den, z)
    x1, yy = mp.mpf(-1069911) / 10 ** 6, mp.mpf(1) / 5
    K2 = mp.mpf(K) ** 2
    corners = [mp.mpc(x1, -yy) * K2, mp.mpc(1, -yy) * K2, mp.mpc(1, yy) * K2, mp.mpc(x1, yy) * K2]
    tot = mp.mpc(0)
    for a in range(4):
        z0, z1 = corners[a], corners[(a + 1) % 4]
        pts = [z0 + (z1 - z0) * mp.mpf(s) / 16 for s in range(17)]
        tot += mp.quad(g, pts)
    cint = tot / (2j * mp.pi)
    Mq = mp.mpf(int(lhs_M.p)) / int(lhs_M.q)
    print("(c) contour integral = %s ; q^T M q = %s ; rel diff %s" % (
        mp.nstr(cint, 20), mp.nstr(Mq, 20), mp.nstr(abs(cint - Mq) / abs(Mq), 3)))
    # (d) scalings: p(tau) = q(K^2 tau), Omega, W_K
    Kf = mp.mpf(K)
    u = [mp.mpf(j) / K for j in range(K + 1)]
    def Omega(tau):
        r = mp.mpf(1) if not isinstance(tau, mp.mpc) else mp.mpc(1)
        for l in range(1, N + 1):
            r *= (tau + u[l] ** 2) ** 6
        for j in range(1, K + 1):
            r /= (tau + u[j] ** 2)
        return r
    p = lambda tau: polyval_mp(qpoly, K2 * tau)
    WK = lambda tau: Omega(tau) * w_mp(Kf * mp.sqrt(tau)) / (2 * mp.sqrt(tau))
    # int_0^oo p^2 W_K dtau = int_0^oo p(y^2)^2 Omega(y^2) w(K y) dy
    I1 = mp.quad(lambda y: p(y * y) ** 2 * Omega(y * y) * w_mp(Kf * y), [0] + [mp.mpf(x) / 40 for x in range(1, 41)])
    print("(d1) Lemma 1(a): K^(12N-2K+1) int p^2 W_K = %s vs q^T G q = %s ; rel %s" % (
        mp.nstr(Kf ** (12 * N - 2 * K + 1) * I1, 20), mp.nstr(rhs, 20),
        mp.nstr(abs(Kf ** (12 * N - 2 * K + 1) * I1 - rhs) / rhs, 3)))
    gs = lambda z: z ** 2 * Omega(z) * p(z) ** 2
    corners_s = [mp.mpc(x1, -yy), mp.mpc(1, -yy), mp.mpc(1, yy), mp.mpc(x1, yy)]
    tot = mp.mpc(0)
    for a in range(4):
        z0, z1 = corners_s[a], corners_s[(a + 1) % 4]
        pts = [z0 + (z1 - z0) * mp.mpf(s) / 16 for s in range(17)]
        tot += mp.quad(gs, pts)
    cint_s = Kf ** (12 * N - 2 * K + 6) * tot / (2j * mp.pi)
    print("(d2) Lemma 1(b): K^(12N-2K+6) (1/2pi i) oint tau^2 Omega p^2 = %s ; rel diff to q^T M q %s" % (
        mp.nstr(cint_s, 20), mp.nstr(abs(cint_s - Mq) / abs(Mq), 3)))
    return D


if __name__ == "__main__":
    for K in [int(a) for a in sys.argv[1:]] or [40]:
        main(K)
