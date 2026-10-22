"""18-redteam-mu96: certify M0 (sup over t>=0 of 2U-V), m0 (inf over S of 2U-V), I(rho), C*, U for rho1600.
Usage: python run_real_certs.py [M0num M0den m0num m0den]
"""
import sys, time, json
from flint import arb, fmpq, ctx
sys.path.insert(0, './scripts/18-redteam-mu96')
from common import load_measure, admissibility, energy_A2, Cstar, ALPHA, LAM, ROOT
from certreal import RealCert, certify

ctx.prec = 128
OUT = ROOT + "/data/18-redteam-mu96/real_certs.log"
logf = open(OUT, "a")


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    logf.write(s + "\n")
    logf.flush()


M0 = fmpq(-66607665, 10 ** 7)
m0 = fmpq(-6660814, 10 ** 6)
if len(sys.argv) >= 5:
    M0 = fmpq(int(sys.argv[1]), int(sys.argv[2]))
    m0 = fmpq(int(sys.argv[3]), int(sys.argv[4]))

log("=" * 100)
log("18-redteam-mu96 real-axis certificates, %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
A, B, C = load_measure()
adm = admissibility(A, B, C)
log("admissibility:", json.dumps(adm))
assert adm["all_c_positive"] and adm["mass_exact"] and adm["strictly_nested"] and adm["support_in_(0,2)"] \
    and adm["reg_const_le_60"]
rc = RealCert(A, B, C)

# ---- tail t >= 2: 2U - V <= (2 lam + 6 al - 1) log t + 2 al^3/t - 2 pi sqrt t =: g(t), decreasing on [2,inf)
#      g'(t) = (2lam+6al-1)/t - 2al^3/t^2 - pi/sqrt t < 0 since (2lam+6al-1)/t <= 0.65 < pi/sqrt(t) needs t < ...
#      we check g'(t) < 0 for t >= 2 via: (2lam+6al-1)/t - pi/sqrt t = (1/sqrt t)((2lam+6al-1)/sqrt t - pi) < 0.
al, lam = arb(ALPHA), arb(LAM)
coef = 2 * lam + 6 * al - 1
assert coef > 0 and coef / arb(2).sqrt() < arb.pi()
g2 = coef * arb(2).log() + 2 * al ** 3 / 2 - 2 * arb.pi() * arb(2).sqrt()
fauzan68 = arb(13) / 10 * arb(2).log() + 6 * al ** 3 / 2 - 2 * arb.pi() * arb(2).sqrt()
log("tail: g(2) = %s ; Fauzan (6.8) rhs at 2 = %s ; both < M0 = %s : %s" % (g2, fauzan68, float(M0),
                                                                          bool(g2 < arb(M0) and fauzan68 < arb(M0))))
assert g2 < arb(M0) and fauzan68 < arb(M0)

t0 = time.time()
ok_sup, n_sup, w_sup = certify(rc, fmpq(0), fmpq(2), M0, 'sup', log=log)
log("SUP on [0,2]: f < M0 = %s : %s ; pieces %d ; worst upper bound %s ; [%.1fs]" % (
    M0, ok_sup, n_sup, w_sup, time.time() - t0))

t0 = time.time()
ok_inf, n_inf, w_inf = certify(rc, A[-1], B[-1], m0, 'inf', log=log)
log("INF on S=[%s,%s]: f > m0 = %s : %s ; pieces %d ; worst lower bound %s ; [%.1fs]" % (
    float(A[-1]), float(B[-1]), m0, ok_inf, n_inf, w_inf, time.time() - t0))

I = energy_A2(A, B, C)
Cs = Cstar()
U = lam * arb(M0) - I + Cs
log("I(rho) (A.2) = %s" % I)
log("C*           = %s" % Cs)
log("U = lam*M0 - I + C* = %s   (upper %s)" % (U, U.upper()))
log("claimed U <= -1.38549608255455 : %s" % bool(U < arb(fmpq(-138549608255455, 10 ** 14))))
log("RESULT ok_sup=%s ok_inf=%s" % (ok_sup, ok_inf))
