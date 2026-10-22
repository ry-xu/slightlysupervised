"""Reproduce Fauzan's analytic constant U = lam*M0 - I(rho) + C* <= -2733991/2000000 (Lemma 6.1, App. A)
from his Table 1 (16 arcsine measures), and check 2U^rho - V <= M0 on a fine grid; compare with the
true equilibrium value from equilibrium.py."""
import sys
sys.path.insert(0, './scripts/09-zeta7det')
import numpy as np
import mpmath as mp
from equilibrium import V, fauzan_profile, U_fauzan, Cstar

mp.mp.dps = 30
T1 = [  # 1e12 a_j, 1e12 b_j, 1e12 c_j
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
a = [mp.mpf(x[0]) / 10 ** 12 for x in T1]
b = [mp.mpf(x[1]) / 10 ** 12 for x in T1]
c = [mp.mpf(x[2]) / 10 ** 12 for x in T1]
lam = mp.mpf(37) / 40
alpha = mp.mpf(3) / 40
print("sum c_j =", mp.nstr(sum(c), 15), " (37/40 =", mp.nstr(lam, 15), ")")
# (A.2): I(rho) = sum_j (S_j^2 - S_{j-1}^2) log((b_j-a_j)/4)
S = [mp.mpf(0)]
for cj in c:
    S.append(S[-1] + cj)
I = sum((S[j + 1] ** 2 - S[j] ** 2) * mp.log((b[j] - a[j]) / 4) for j in range(16))
Cs = -2 * lam + 12 * alpha * lam * (1 - mp.log(alpha)) + 3 * lam ** 2 - 2 * lam ** 2 * mp.log(2 * lam)
M0 = mp.mpf(-1329) / 200
U = lam * M0 - I + Cs
print("I(rho) =", mp.nstr(I, 15), " C* =", mp.nstr(Cs, 15))
print("U = lam M0 - I + C* =", mp.nstr(U, 15), "  paper bound -2733991/2000000 =", -2733991 / 2000000)


def Uarc(t, aa, bb):
    if aa <= t <= bb:
        return mp.log((bb - aa) / 4)
    m = (aa + bb) / 2
    return mp.log(abs(t - m) + mp.sqrt((t - aa) * (t - bb))) - mp.log(2)


def Urho(t):
    return sum(c[j] * Uarc(t, a[j], b[j]) for j in range(16))


prof = fauzan_profile(6, 3 / 40)
ts = np.concatenate([np.geomspace(1e-7, 1e-2, 3000), np.linspace(1e-2, 2.5, 5000)])
worst = -1e9
for t in ts:
    val = float(2 * Urho(mp.mpf(t))) - float(V(t, prof))
    worst = max(worst, val)
print("max over grid of 2U^rho - V =", worst, "  (paper M0 = -1329/200 = -6.645, App.A bound -6.6450021)")
Ut, info = U_fauzan(6, 3 / 40, 37 / 40, n=8192)
print("true equilibrium U* =", Ut, info)
print("Fauzan's measure is within", float(U) - Ut, "of the optimum")
