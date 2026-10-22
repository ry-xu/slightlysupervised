"""Fit a nested-arcsine test measure rho = sum_j c_j omega[a_j,b_j] (Fauzan App. A form, signed c_j allowed)
to the weighted equilibrium problem of Fauzan's zeta(5) construction, by collocation:
   2 U^rho(t_i) - V(t_i) = ell   at collocation points t_i in the support,  sum c_j = lam.
Then U_bound = lam * sup_t (2U^rho - V) - I(rho) + C*   (Lemma 6.1 / (6.4) with M0 = sup).
I(rho) for nested intervals: sum_j (S_j^2 - S_{j-1}^2) log((b_j-a_j)/4)   (Fauzan (A.2)).
Output: intervals and weights (float) + the numerical U_bound; certification is done in arcsine_certify.py."""
import sys, json
import numpy as np
sys.path.insert(0, './scripts/09-zeta7det')
import equilibrium as eq

r, al, lam = 6, 3 / 40, 37 / 40
prof = eq.fauzan_profile(r, al)


def Uarc(t, a, b):
    t = np.asarray(t, dtype=float)
    out = np.empty_like(t)
    inside = (t >= a) & (t <= b)
    out[inside] = np.log((b - a) / 4)
    m = (a + b) / 2
    o = ~inside
    tt = t[o]
    out[o] = np.log(np.abs(tt - m) + np.sqrt((tt - a) * (tt - b))) - np.log(2)
    return out


def I_nested(a, b, c):
    S = np.concatenate([[0.0], np.cumsum(c)])
    return float(np.sum((S[1:] ** 2 - S[:-1] ** 2) * np.log((b - a) / 4)))


def fit(J, n_eq=8192, verbose=True):
    U, info = eq.U_fauzan(r, al, lam, n=n_eq)
    A, B = info['a'], info['b']
    # nested intervals: a_j decreasing to ~A, b_j increasing to ~B, innermost around the field minimum
    # choose them as level sets of the equilibrium density's "arcsine decomposition": geometric toward A, linear toward B
    s = np.linspace(0, 1, J + 1)[1:]
    q = 0.0059  # approx location of min of V
    a_j = q * (A / q) ** s * 1.0
    b_j = q + (B - q) * s ** 1.0
    a_j[-1], b_j[-1] = A, B
    # collocation points: Chebyshev-distributed in log scale near A and linear near B
    xs = np.cos(np.pi * (np.arange(4 * J) + 0.5) / (4 * J))
    tcol = (A + B) / 2 + (B - A) / 2 * xs
    tcol = np.unique(np.concatenate([tcol, np.geomspace(A * 1.0001, q, 2 * J)]))
    Vc = eq.V(tcol, prof)
    M = np.stack([2 * Uarc(tcol, a_j[j], b_j[j]) for j in range(J)], axis=1)
    # unknowns c (J) and ell; minimize ||M c - V - ell||^2 subject to sum c = lam  (least squares with constraint)
    Aug = np.hstack([M, -np.ones((len(tcol), 1))])
    # constraint via elimination: add heavily weighted row
    w = 1e6
    Aaug = np.vstack([Aug, w * np.concatenate([np.ones(J), [0.0]])])
    baug = np.concatenate([Vc, [w * lam]])
    sol, *_ = np.linalg.lstsq(Aaug, baug, rcond=None)
    c, ell = sol[:J], sol[J]
    # evaluate sup of 2U - V on a fine grid over [0, 3]
    tg = np.unique(np.concatenate([np.geomspace(1e-9, 1e-2, 20000), np.linspace(1e-2, 3, 20000)]))
    Ug = sum(c[j] * Uarc(tg, a_j[j], b_j[j]) for j in range(J))
    val = 2 * Ug - eq.V(tg, prof)
    M0 = val.max()
    I = I_nested(a_j, b_j, c)
    Cs = eq.Cstar(r, al, lam)
    Ub = lam * M0 - I + Cs
    if verbose:
        print(f"J={J}: ell(fit)={ell:.6f}  sup(2U-V)={M0:.6f} at t={tg[val.argmax()]:.3g}  I={I:.6f}  U_bound={Ub:.6f}  (U*={U:.6f}, gap {Ub-U:.2e})  min c={c.min():.2e}")
    return dict(J=J, a=a_j.tolist(), b=b_j.tolist(), c=c.tolist(), M0=M0, I=I, U_bound=Ub, Ustar=U)


if __name__ == "__main__":
    res = []
    for J in (16, 24, 32, 48):
        res.append(fit(J))
    json.dump(res, open('./data/09-zeta7det/arcsine_fit.json', 'w'))
