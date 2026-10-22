"""15-redteam-mu5: common definitions for the analytic-constant check (independent of 09's code).

Conventions taken from Fauzan, "zeta(5) is irrational" (Zenodo 10.5281/zenodo.22826419), Sec. 6 and App. A:
  V(t)   = 2 pi sqrt t + int_0^1 log(t+u^2) du - 6 int_0^alpha log(t+u^2) du          (6.1)
  U^rho  = int log|t-u| drho(u),  I(rho) = int int log|t-u| drho drho
  U_bound = lam*M0 - I(rho) + C*,  M0 >= sup_{t>=0} (2U^rho - V)                         (6.4)
  C* = -2 lam + 12 alpha lam (1 - log alpha) + 3 lam^2 - 2 lam^2 log(2 lam)              (6.3)
lam = 37/40, alpha = 3/40.

The measure under test is 09's 121-component nested arcsine measure, read from
data/09-zeta7det/arcsine_opt.json (key '3p'), with the conversion convention stated in 09's notes:
a_j, b_j, c_j rounded to 10^-12 Z, last weight adjusted so that sum c_j = 37/40 exactly.
Only the DATA is taken from 09; all formulas below are re-derived and re-implemented here.
"""
import json
from fractions import Fraction as Fr
import math

LAM = Fr(37, 40)
AL = Fr(3, 40)
DATA09 = './data/09-zeta7det/arcsine_opt.json'

# Fauzan Table 1 (10^12 a_j, 10^12 b_j, 10^12 c_j), transcribed from the PDF (page 22)
FAUZAN_T1 = [
    (3906748086, 8992695531, 10515596180),
    (2312248264, 15340997855, 29471737793),
    (1402286665, 25730180724, 42934365099),
    (881725356, 41909578246, 58204231966),
    (578197906, 65851089563, 69037621310),
    (396324613, 99481037884, 78873099189),
    (283911191, 144325727458, 84856120711),
    (212206188, 201105762729, 88396082127),
    (165097686, 269345996903, 88303382125),
    (133347132, 347089554156, 85472321255),
    (111522114, 430806704415, 78899184238),
    (96349355, 515561896511, 70353471918),
    (85815639, 595448778546, 58838976615),
    (78667711, 664241383483, 44421321106),
    (74129565, 716160577112, 30462865791),
    (71741310, 746637295669, 5959622577),
]


def load_fauzan():
    D = 10 ** 12
    a = [Fr(x, D) for x, _, _ in FAUZAN_T1]
    b = [Fr(y, D) for _, y, _ in FAUZAN_T1]
    c = [Fr(z, D) for _, _, z in FAUZAN_T1]
    return a, b, c


def load_09(key='3p'):
    d = json.load(open(DATA09))[key]
    D = 10 ** 12
    a = [Fr(round(x * D), D) for x in d['a']]
    b = [Fr(round(x * D), D) for x in d['b']]
    c = [Fr(round(x * D), D) for x in d['c']]
    c[-1] += LAM - sum(c)
    return a, b, c


def admissibility(a, b, c, verbose=True):
    """Checks every hypothesis of Lemma 6.1 / Sec. 6.1 / App. A.3 that the measure must satisfy.
    Returns dict of booleans and the actual numbers."""
    J = len(a)
    res = {}
    res['J'] = J
    res['all_c_positive'] = all(x > 0 for x in c)
    res['min_c'] = float(min(c))
    res['mass_exact'] = (sum(c) == LAM)
    # nesting with index 1 innermost: a_J < ... < a_1 < b_1 < ... < b_J
    res['nested'] = all(a[j] < a[j - 1] for j in range(1, J)) and all(b[j] > b[j - 1] for j in range(1, J)) and a[0] < b[0]
    res['min_len'] = float(min(bb - aa for aa, bb in zip(a, b)))
    res['len_gt_1_225'] = all(bb - aa > Fr(1, 225) for aa, bb in zip(a, b))
    res['support_in_(0,2)'] = (min(a) > 0) and (max(b) < 2)
    res['support'] = (float(min(a)), float(max(b)))
    # regularization constant of Sec. 6.1:  sum_j 8 c_j sqrt(2)/(pi sqrt(L_j))  must be <= 60
    reg = sum(8 * float(cj) * math.sqrt(2) / (math.pi * math.sqrt(float(bj - aj))) for aj, bj, cj in zip(a, b, c))
    res['reg_const'] = reg
    res['reg_le_60'] = reg <= 60
    # all denominators divide 10^12 ?
    res['denoms_10^12'] = all((10 ** 12) % x.denominator == 0 for x in a + b + c)
    if verbose:
        for k, v in res.items():
            print(f"   {k}: {v}")
    return res


# ---------------- float versions (exploration only, not used in the certificate) ----------------
def U_arc_f(t, a, b):
    if a <= t <= b:
        return math.log((b - a) / 4)
    m = (a + b) / 2
    return math.log((abs(t - m) + math.sqrt((t - a) * (t - b))) / 2)


def V_f(t):
    al = 3 / 40
    if t == 0:
        return -2 - 12 * al * math.log(al) + 12 * al
    st = math.sqrt(t)
    i1 = math.log(1 + t) - 2 + 2 * st * math.atan(1 / st)            # int_0^1 log(t+u^2) du
    i2 = al * math.log(t + al * al) - 2 * al + 2 * st * math.atan(al / st)  # int_0^al log(t+u^2) du
    return 2 * math.pi * st + i1 - 6 * i2
