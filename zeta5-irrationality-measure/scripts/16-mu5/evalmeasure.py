"""16-mu5: float evaluation of a nested arcsine measure (a_j, b_j, c_j), index 0 innermost.
Reports
  * the Lemma 6.1 objective  U(rho) = lam * sup_{t>=0}(2U^rho - V) - I(rho) + C*      (Fauzan (6.4), (A.2))
  * the norm-bound objective rho(rho) = sup_{s in [al,1]}(2U^rho(-s^2) + phi(s)) - inf_{t in supp}(2U^rho - V)
"""
import sys, json, math
import numpy as np

AL, LAM = 3 / 40, 37 / 40
CSTAR = -2 * LAM + 12 * AL * LAM * (1 - math.log(AL)) + 3 * LAM ** 2 - 2 * LAM ** 2 * math.log(2 * LAM)


def V(t):
    t = np.asarray(t, float)
    st = np.sqrt(np.maximum(t, 1e-300))
    i1 = np.log(1 + t) - 2 + 2 * st * np.arctan(1 / st)
    i2 = AL * np.log(t + AL * AL) - 2 * AL + 2 * st * np.arctan(AL / st)
    return 2 * np.pi * st + i1 - 6 * i2


def G(x):
    x = np.asarray(x, float)
    xs = np.where(x > 0, x, 1.0)
    return np.where(x > 0, xs * np.log(xs) - xs, 0.0)


def phi(s):
    return 6 * (G(s + AL) - G(s - AL)) - G(1 + s) - G(1 - s)


def U_arc(t, a, b):
    t = np.asarray(t, float)
    m = (a + b) / 2
    inside = (t >= a) & (t <= b)
    q = np.maximum((t - a) * (t - b), 0)
    with np.errstate(divide='ignore'):
        out = np.log((np.abs(t - m) + np.sqrt(q)) / 2)
    return np.where(inside, np.log((b - a) / 4), out)


def Upot(t, a, b, c):
    t = np.asarray(t, float)
    out = np.zeros_like(t)
    for aj, bj, cj in zip(a, b, c):
        out += cj * U_arc(t, aj, bj)
    return out


def energy(a, b, c):
    """(A.2): I(rho) = sum_j (S_j^2 - S_{j-1}^2) log(L_j/4), S_j cumulative from innermost"""
    S = 0.0
    I = 0.0
    for aj, bj, cj in zip(a, b, c):
        S1 = S + cj
        I += (S1 * S1 - S * S) * math.log((bj - aj) / 4)
        S = S1
    return I


def tgrid(a, b, n=20000, extra=()):
    lo, hi = min(a), max(b)
    pts = np.concatenate([np.geomspace(1e-9, 2, n), np.linspace(0, 2, n), np.array(list(a) + list(b) + list(extra))])
    return np.unique(pts)


def evaluate(a, b, c, n=20000, ns=4000):
    a, b, c = map(np.asarray, (a, b, c))
    t = tgrid(a, b, n)
    f = 2 * Upot(t, a, b, c) - V(t)
    sup_all = f.max()
    t_sup = t[f.argmax()]
    S0, S1 = a.min(), b.max()
    ts = np.unique(np.concatenate([np.linspace(S0, S1, n), np.geomspace(S0, S1, n), a, b]))
    fs = 2 * Upot(ts, a, b, c) - V(ts)
    m0 = fs.min()
    t_m0 = ts[fs.argmin()]
    s = np.linspace(AL, 1, ns)
    g = 2 * Upot(-s * s, a, b, c) + phi(s)
    M1 = g.max()
    I = energy(a, b, c)
    return dict(U=LAM * sup_all - I + CSTAR, M0=sup_all, t_sup=t_sup, I=I, m0=m0, t_m0=t_m0, M1=M1,
                s_M1=s[g.argmax()], rho=M1 - m0, mass=c.sum(), support=(S0, S1))


if __name__ == "__main__":
    for fn in sys.argv[1:]:
        d = json.load(open(fn))
        if 'a' not in d:
            d = d['3p']
        r = evaluate(d['a'], d['b'], d['c'])
        print(f"{fn.split('/')[-1]}: J={len(d['a'])} mass={r['mass']:.12f} supp=[{r['support'][0]:.4g},{r['support'][1]:.6f}]"
              f"\n   U = {r['U']:.7f} (sup 2U-V = {r['M0']:.7f} at t={r['t_sup']:.4g}; I={r['I']:.7f})"
              f"\n   rho = {r['rho']:.6f} (M1 = {r['M1']:.6f} at s={r['s_M1']:.4f}; inf_S(2U-V) = {r['m0']:.6f} at t={r['t_m0']:.4g})")
