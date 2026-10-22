"""Optimize the weights of a nested-arcsine test measure (Fauzan App. A form) for
   U_bound(c) = lam * sup_t (2 U^rho(t) - V(t)) - I(rho) + C*,     rho = sum_j c_j omega[a_j,b_j]
(signed c allowed: Lemma 6.2 needs only mass(sigma - rho) = 0).  I(rho) = sum_{j,k} c_j c_k L_{max(j,k)},
L_j = log((b_j-a_j)/4) (nested intervals, innermost j=0).  Sequential LP on a fine grid (HiGHS):
   minimize lam*tau - [2 c0^T Q c - c0^T Q c0]   s.t.  2 (M c)_i - tau <= V_i,  sum c = lam,  |c - c0| <= trust.
Intervals: Fauzan's Table 1, optionally refined by inserting geometric/arith. midpoints."""
import sys, json
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, './scripts/09-zeta7det')
import equilibrium as eq

r, al, lam = 6, 3 / 40, 37 / 40
POSMIN = 1e-6   # keep weights positive (monotone-potential certificate)
prof = eq.fauzan_profile(r, al)
T1 = [(3906748086, 8992695531, 10515596180), (2312248264, 15340997855, 29471737793), (1402286665, 25730180724, 42934365099),
      (881725356, 41909578246, 58204231966), (578197906, 65851089563, 69037621310), (396324613, 99481037884, 78873099189),
      (283911191, 144325727458, 84856120711), (212206188, 201105762729, 88396082127), (165097686, 269345996903, 88303382125),
      (133347132, 347089554156, 85472321255), (111522114, 430806704415, 78899184238), (96349355, 515561896511, 70353471918),
      (85815639, 595448778546, 58838976615), (78667711, 664241383483, 44421321106), (74129565, 716160577112, 30462865791),
      (71741310, 746637295669, 5959622577)]


def Uarc(t, a, b):
    t = np.asarray(t, dtype=float)
    out = np.empty_like(t)
    inside = (t >= a) & (t <= b)
    out[inside] = np.log((b - a) / 4)
    o = ~inside
    tt = t[o]
    out[o] = np.log(np.abs(tt - (a + b) / 2) + np.sqrt((tt - a) * (tt - b))) - np.log(2)
    return out


def refine(a, b, levels):
    """insert intermediate nested intervals (geometric mean of a's, arithmetic of b's)"""
    for _ in range(levels):
        na, nb = [a[0]], [b[0]]
        for j in range(1, len(a)):
            na.append(np.sqrt(a[j - 1] * a[j])); nb.append(0.5 * (b[j - 1] + b[j]))
            na.append(a[j]); nb.append(b[j])
        a, b = np.array(na), np.array(nb)
    return a, b


def grid(a, b):
    pts = [np.geomspace(1e-9, a.min(), 400), np.linspace(b.max(), 3.0, 600)]
    edges = np.unique(np.concatenate([a, b]))
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi / lo > 1.5 and lo > 0:
            pts.append(np.geomspace(lo, hi, 60))
        else:
            pts.append(np.linspace(lo, hi, 60))
    return np.unique(np.concatenate(pts))


def evaluate(a, b, c, tg):
    U = sum(c[j] * Uarc(tg, a[j], b[j]) for j in range(len(a)))
    val = 2 * U - eq.V(tg, prof)
    L = np.log((b - a) / 4)
    J = len(a)
    Q = np.array([[L[max(j, k)] for k in range(J)] for j in range(J)])
    I = c @ Q @ c
    return val.max(), I, lam * val.max() - I + eq.Cstar(r, al, lam), tg[val.argmax()]


def optimize(a, b, c0, iters=30, trust=0.02, verbose=True):
    tg = grid(a, b)
    J = len(a)
    M = np.stack([2 * Uarc(tg, a[j], b[j]) for j in range(J)], axis=1)
    Vg = eq.V(tg, prof)
    L = np.log((b - a) / 4)
    Q = np.array([[L[max(j, k)] for k in range(J)] for j in range(J)])
    c = np.array(c0, dtype=float)
    best = None
    for it in range(iters):
        g = 2 * Q @ c                       # gradient of I
        # variables x = (c_1..c_J, tau); minimize lam*tau - g.c
        cost = np.concatenate([-g, [lam]])
        A_ub = np.hstack([M, -np.ones((len(tg), 1))])
        b_ub = Vg
        A_eq = np.concatenate([np.ones(J), [0.0]])[None, :]
        bounds = [(max(ci - trust, POSMIN), ci + trust) for ci in c] + [(None, None)]
        res = linprog(cost, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[lam], bounds=bounds, method='highs')
        if res.status != 0:
            trust /= 2; continue
        cn = res.x[:J]
        M0, I, Ub, targ = evaluate(a, b, cn, tg)
        if best is None or Ub < best[0] - 1e-12:
            best = (Ub, cn.copy(), M0, I)
            c = cn
        else:
            trust /= 2
        if verbose:
            print(f"  it {it}: U_bound={Ub:.7f} sup={M0:.7f} I={I:.7f} trust={trust:.2e}", flush=True)
        if trust < 1e-7:
            break
    return best


if __name__ == "__main__":
    a0 = np.array([x[0] for x in T1]) / 1e12
    b0 = np.array([x[1] for x in T1]) / 1e12
    c0 = np.array([x[2] for x in T1]) / 1e12
    tg = grid(a0, b0)
    print("Fauzan Table 1 on our grid: sup, I, U_bound, argmax =", evaluate(a0, b0, c0, tg))
    out = {}
    for lev in (0, 1, 2):
        a, b = refine(a0, b0, lev)
        cinit = np.interp(np.arange(len(a)), np.linspace(0, len(a) - 1, len(c0)), c0)
        cinit *= lam / cinit.sum()
        print(f"--- refinement level {lev}: J={len(a)}")
        best = optimize(a, b, cinit, iters=40, trust=0.01)
        print(f"level {lev}: best U_bound = {best[0]:.7f}  (sup {best[2]:.7f}, I {best[3]:.7f}); U* = -1.385510")
        out[lev] = dict(a=a.tolist(), b=b.tolist(), c=best[1].tolist(), U_bound=best[0], M0=best[2], I=best[3])
        json.dump(out, open('./data/09-zeta7det/arcsine_opt.json', 'w'))
