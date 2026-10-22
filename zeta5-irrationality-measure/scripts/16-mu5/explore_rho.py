"""16-mu5: float exploration of the potential-theoretic norm bound for ||B_K|| (Fauzan App. C).

Bound (proved in notes/16-mu5.md, Thm 1):
   limsup (1/K) log ||B_K||  <=  sup_{s in [alpha,1]} [2U^rho(-s^2) + phi(s)]  -  inf_{t in S} [2U^rho(t) - V(t)]
for any (signed) measure rho with mass >= lambda, supp rho = S compact in (0,oo), continuous potential.
   V(t)   = 2 pi sqrt t + int_0^1 log(t+u^2) du - 6 int_0^alpha log(t+u^2) du        (Fauzan (6.1))
   phi(s) = 6[G(s+a) - G(s-a)] - G(1+s) - G(1-s),   G(x) = x log x - x            (= -Vtilde(-s^2))
   U^rho(z) = int log|z-u| drho(u)                                                  (Fauzan's sign)
The equilibrium measure gives the sharp constant rho* (it is also the true growth rate of ||B_K||).
"""
import sys, json, math
import numpy as np
sys.path.insert(0, './scripts/09-zeta7det')
import equilibrium as eq

AL, LAM = 3 / 40, 37 / 40


def G(x):
    x = np.asarray(x, float)
    out = np.where(x > 0, x * np.log(np.where(x > 0, x, 1.0)) - x, 0.0)
    return out


def phi(s):
    return 6 * (G(s + AL) - G(s - AL)) - G(1 + s) - G(1 - s)


def V(t):
    return eq.V(t, eq.fauzan_profile(6, AL))


def U_arc(t, a, b):
    """potential of unit arcsine measure on [a,b] at real t (vectorized)"""
    t = np.asarray(t, float)
    m = (a + b) / 2
    inside = (t >= a) & (t <= b)
    d = np.abs(t - m)
    out = np.log((d + np.sqrt(np.maximum((t - a) * (t - b), 0))) / 2)
    return np.where(inside, np.log((b - a) / 4), out)


def load09(key='3p'):
    d = json.load(open('./data/09-zeta7det/arcsine_opt.json'))[key]
    return np.array(d['a']), np.array(d['b']), np.array(d['c'])


def bound_for_arcsine(a, b, c, ngrid=200001):
    S0, S1 = a.min(), b.max()
    t = np.concatenate([np.linspace(S0, S1, ngrid), np.geomspace(S0, S1, ngrid)])
    t.sort()
    U = sum(cj * U_arc(t, aj, bj) for aj, bj, cj in zip(a, b, c))
    f = 2 * U - V(t)
    m0 = f.min()
    tmin = t[f.argmin()]
    s = np.linspace(AL, 1, 20001)
    x = -s * s
    Ux = sum(cj * U_arc(x, aj, bj) for aj, bj, cj in zip(a, b, c))
    g = 2 * Ux + phi(s)
    M1 = g.max()
    smax = s[g.argmax()]
    # also sup over [0,oo) of 2U-V for reference
    tt = np.concatenate([np.linspace(0, 2, 400001)])
    UU = sum(cj * U_arc(tt, aj, bj) for aj, bj, cj in zip(a, b, c))
    sup_all = (2 * UU - V(tt)).max()
    return dict(m0=m0, tmin=tmin, M1=M1, smax=smax, rho=M1 - m0, sup_f=f.max(), sup_all=sup_all, S=(S0, S1))


def equilibrium_bound(n=8192):
    prof = eq.fauzan_profile(6, AL)
    a, b, ier, resid = eq.solve_support(prof, LAM, n=n, guess=(1e-4, 0.75))
    E, mass, lo, hi, w = eq.energy(a, b, prof, LAM, n=n)
    # quadrature nodes for the measure: y_j = cos(pi j/n), weights w (from eq.energy)
    c, d = (a + b) / 2, (b - a) / 2
    j = np.arange(1, n)
    y = np.cos(np.pi * j / n)
    tn = c + d * y
    F = 0.5 * (lo + hi)
    s = np.linspace(AL, 1, 4001)
    x = -s * s
    Ux = np.array([np.sum(w * np.log(tn - xx)) for xx in x])
    g = 2 * Ux + phi(s)
    return dict(a=a, b=b, mass=mass, F_lo=lo, F_hi=hi, M1=g.max(), smax=s[g.argmax()], rho=g.max() - F,
                g_at_alpha=g[0], g_at_1=g[-1], U=E + eq.Cstar(6, AL, LAM))


if __name__ == "__main__":
    for key in ['0', '1', '2', '3p']:
        a, b, c = load09(key)
        r = bound_for_arcsine(a, b, c)
        print(f"09 measure {key:>2} (J={len(a)}): inf_S(2U-V) = {r['m0']:.6f} at t={r['tmin']:.3e}; sup_S = {r['sup_f']:.6f};"
              f" sup_[0,2] = {r['sup_all']:.6f}; M1 = {r['M1']:.6f} at s={r['smax']:.4f};  rho_bound = {r['rho']:.5f}")
    for n in (2048, 8192):
        r = equilibrium_bound(n)
        print(f"equilibrium n={n}: support [{r['a']:.6g},{r['b']:.6g}] mass {r['mass']:.8f} F in [{r['F_lo']:.8f},{r['F_hi']:.8f}]"
              f"  M1 = {r['M1']:.6f} at s = {r['smax']:.5f} (g(alpha)={r['g_at_alpha']:.5f}, g(1)={r['g_at_1']:.5f});"
              f"  rho* = {r['rho']:.6f};  U* = {r['U']:.7f}")
    # Fauzan's crude kappa for comparison
    a_, b_ = 1 / 18, 1 / 3
    Y = (2 + a_ + b_) / (b_ - a_)
    R = Y + math.sqrt(Y * Y - 1)
    kap = 2 * math.pi * math.sqrt(b_) + math.log(1 + b_) - 6 * AL * math.log(a_) + LAM * (1 - math.log(AL * LAM)) + 2 * LAM * math.log(R)
    print("Fauzan kappa(1/18,1/3) =", kap)
