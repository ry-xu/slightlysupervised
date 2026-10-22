"""
Final exact ledger for a certified second trial.            [rigorous: exact rationals + log 2 bounds]
usage: final_ledger.py K2 MX MY R2      (certified caps, decimal strings -> exact rationals)
Optimises the target U (1e-8 grid) and the cost c2 (1e-5 grid) jointly, then kappa (1e-7 grid),
and checks every inequality of Litsey's generic ledger (ledger.check).  Prints the exponent mu = rho*kappa.
"""
import sys, json
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from ledger import check, best_c2, best_kappa, LO2, HI2, CHI, W

K2, MX, MY, R2 = (Fr(s) for s in sys.argv[1:5])
F2 = MX + MY
delta = Fr(1, 10**7)
Umax = Fr(7, 10) * LO2 - K2 - F2           # h0 <= 0 guaranteed for U <= Umax (log2 >= LO2)
best = None
for slack_e8 in range(0, 2000, 10):        # U = Umax - slack, slack in [0, 2e-5]
    U = Fr(int(Umax * 10**8) - slack_e8, 10**8)
    c2 = best_c2(U, K2, F2, R2)
    gap_lo = 400 * (U - Fr(8795, 4800)) - 242 * HI2
    eta = Fr(int(gap_lo * 10**7) - 1, 10**7)
    E = eta + CHI
    kap = best_kappa(E, delta)
    mu = (20 * c2 + W) * kap
    if best is None or mu < best[0]:
        best = (mu, U, c2, kap)
mu, U, c2, kap = best
print('chosen: U =', U, '=', float(U), ' c2 =', c2, ' kappa =', kap)
ok, out = check(U, K2, F2, R2, c2, kap, delta)
print(json.dumps(dict(ok=ok, U=str(U), c2=str(c2), kappa=str(kap), mu=float(mu), mu_exact=str(mu),
                      rho=out['rho'], eta=out['eta'], eta_plus_chi=out['eta_plus_chi'],
                      limit=out['mu_limit_18rho_over_E'])))
