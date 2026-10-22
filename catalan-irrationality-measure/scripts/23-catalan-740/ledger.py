"""
Exact replay of the rational inequalities in Litsey's generic 740-ledger, for new constants.
(Lean sources: Parametric/StrongBudgets.lean, GapBudgets.lean, GapFinalTransfer.lean, U6Main.lean,
 GapContentTransfer.lean `measure_of_content_on`.)               [rigorous: exact rationals + log 2 bounds]

Inputs (all exact rationals):
  U    : new dual target (strongTarget)
  K2   : normCap of the second trial (>= 2||p||^2 + ||v||^2/2, certified)
  F2   : fieldCap = MX + MY of the second trial (certified upper bounds)
  R2   : marked residual cap (gapResidualCap, certified) on |x| >= 401/2500
  c2   : new second marked cost (gapSecondCost)
  kappa, delta : PNT-window parameters;  mu = rho * kappa.
First trial unchanged: K1 = 67761/100000, F1 = -247913/100000, lam = 168493357/100000000, c1 = 421/100.
Weight cost w = 17913/1000; content chi = 507/2000 + 3331/10000 + 624/10000.
log 2 in [LO2, HI2] with LO2 = 0.6931471803, HI2 = 0.6931471808 (Real.log_two_gt_d9 / lt_d9).
"""
from fractions import Fraction as Fr

LO2 = Fr(6931471803, 10**10)
HI2 = Fr(6931471808, 10**10)
K1, F1, LAM, C1 = Fr(67761, 100000), Fr(-247913, 100000), Fr(168493357, 100000000), Fr(421, 100)
W = Fr(17913, 1000)
CHI = Fr(507, 2000) + Fr(3331, 10000) + Fr(624, 10000)


def upper_of(expr_lo, expr_hi):
    return max(expr_lo, expr_hi)


def lin_in_log2(a, b):
    """bounds of a + b*log2."""
    v1, v2 = a + b * LO2, a + b * HI2
    return min(v1, v2), max(v1, v2)


def check(U, K2, F2, R2, c2, kappa, delta, verbose=True):
    out = {}
    # (L2) second budget:  h0 + alpha a + beta t - 2 log2 a^2 - log2 t^2 <= 0  for a,t >= 0
    h0_lo, h0_hi = lin_in_log2(K2 + F2 + U, Fr(-7, 10))
    al_lo, al_hi = lin_in_log2(-F2 + R2 - c2, Fr(13, 5))
    be_lo, be_hi = lin_in_log2(-R2 - Fr(3, 4), Fr(1))
    al = max(al_hi, 0); be = max(be_hi, 0)
    worst2 = h0_hi + al * al / (8 * LO2) + be * be / (4 * LO2)
    out['L2_second_budget'] = (worst2 < 0, float(worst2))
    # (L3) first budget at U (all real a): need h0 <= 0, alpha >= 0, h0 + alpha^2/(8 log2) <= 0
    g_lo, g_hi = lin_in_log2(K1 + F1 + U, Fr(-7, 10))
    a1_lo, a1_hi = lin_in_log2(-F1 + LAM / 5 - C1, Fr(13, 5))
    worst1 = g_hi + max(a1_hi, 0) ** 2 / (8 * LO2)
    out['L3_first_budget'] = (g_hi <= 0 and a1_lo >= 0 and worst1 < 0, float(worst1))
    # (L4) c1 <= c2
    out['L4_c1_le_c2'] = (C1 <= c2, float(c2 - C1))
    rho = 20 * c2 + W
    out['rho'] = float(rho)
    # (L6) eta < 400 * (U - (121/200) log2 - 8795/4800); take eta as a rational just below
    gap_lo = 400 * (U - Fr(8795, 4800)) - 242 * HI2
    eta = Fr(int(gap_lo * 10**7) - 1, 10**7)     # rational strictly below the certified lower bound
    out['eta'] = float(eta)
    out['L6_eta'] = (eta < gap_lo, float(gap_lo))
    E = eta + CHI
    out['eta_plus_chi'] = float(E)
    # (L7) PNT window conditions
    out['L7'] = (2 < kappa and 0 < delta and 1 < 2 * (1 - delta * (2 * kappa - 1)) and E < 36
                 and 18 + 4 * 18 * delta * kappa < E * kappa,
                 float(E * kappa - 18 - 72 * delta * kappa))
    mu = rho * kappa
    out['mu'] = float(mu)
    out['mu_exact'] = str(mu)
    out['mu_limit_18rho_over_E'] = float(18 * rho / E)
    ok = all(v[0] for k, v in out.items() if k.startswith('L'))
    out['ALL'] = ok
    if verbose:
        for k, v in out.items():
            print(k, v)
    return ok, out


def best_c2(U, K2, F2, R2):
    """smallest c2 (in 1/10^5 steps) for which the second budget closes."""
    c = Fr(400000, 100000)
    step = Fr(1, 100000)
    lo, hi = Fr(4), Fr(5)
    for _ in range(60):
        mid = (lo + hi) / 2
        h0_lo, h0_hi = lin_in_log2(K2 + F2 + U, Fr(-7, 10))
        al_lo, al_hi = lin_in_log2(-F2 + R2 - mid, Fr(13, 5))
        be_lo, be_hi = lin_in_log2(-R2 - Fr(3, 4), Fr(1))
        worst = h0_hi + max(al_hi, 0) ** 2 / (8 * LO2) + max(be_hi, 0) ** 2 / (4 * LO2)
        if worst < 0:
            hi = mid
        else:
            lo = mid
    return Fr(int(hi * 10**5) + 1, 10**5)


def best_kappa(eta_plus_chi, delta):
    """smallest kappa (rational, 1e-7 grid) with 18 + 72 delta kappa < E kappa."""
    k = Fr(18) / (eta_plus_chi - 72 * delta)
    return Fr(int(k * 10**7) + 1, 10**7)


if __name__ == "__main__":
    # control: Litsey's own constants must pass and reproduce rho*kappa = 740
    U = Fr(22563367, 10000000)
    K2 = Fr(3096999, 4000000); F2 = Fr(-127269173, 50000000); R2 = Fr(-131, 2500)
    c2 = Fr(13421, 3125)
    print('--- control: Litsey constants, kappa = 3700000/519037, delta = 1e-7')
    check(U, K2, F2, R2, c2, Fr(3700000, 519037), Fr(1, 10**7))
