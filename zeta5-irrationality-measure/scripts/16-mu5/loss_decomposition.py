"""16-mu5: where Fauzan's norm constant kappa = 13.866 loses (float; all quantities in units of K in the exponent).

kappa(a,b) = F(a,b) + lam(1 - log(al lam)) + 2 lam log R     (Fauzan (C.6)), a = 1/18, b = 1/3.
In the language of Theorem 1 (notes/16-mu5.md) Fauzan's proof uses the comparison measure rho_F = lam*omega[a,b] and
three decoupled crude bounds.  Remove them one at a time:
  A: (C.2) weight bound F(a,b) = 2pi sqrt b + log(1+b) - 6 al log a   ->  exact max_{[a,b]} V
  B: (C.4) |r_j| <= exp(lam K (1 - log al lam))                      ->  exact exponent max_s phi(s)
  C: decoupled  max_j |r_j| * max_{[-1,0]} |p|^2                      ->  coupled sup_s (2 U^{rho_F}(-s^2) + phi(s))
  D: best single interval [a,b]                                       (still rho = lam * arcsine)
  E: equilibrium measure (sharp constant rho*)
"""
import sys, math, json
import numpy as np
from scipy.optimize import minimize
sys.path.insert(0, './scripts/16-mu5')
import evalmeasure as ev
import explore_rho as ex

AL, LAM = 3 / 40, 37 / 40
s = np.linspace(AL, 1, 200001)
phis = ev.phi(s)
phimax = phis.max()


def parts(a, b):
    L = b - a
    m, r = (a + b) / 2, L / 2
    Y = (2 + a + b) / (b - a)
    R = Y + math.sqrt(Y * Y - 1)
    F = 2 * math.pi * math.sqrt(b) + math.log(1 + b) - 6 * AL * math.log(a)
    t = np.linspace(a, b, 200001)
    Vmax = ev.V(t).max()
    kF = F + LAM * (1 - math.log(AL * LAM)) + 2 * LAM * math.log(R)
    kA = Vmax + LAM * (1 - math.log(AL * LAM)) + 2 * LAM * math.log(R)
    kB = Vmax + phimax + 2 * LAM * math.log(R)
    Uneg = LAM * ev.U_arc(-s * s, a, b)
    kC = Vmax - 2 * LAM * math.log(L / 4) + (2 * Uneg + phis).max()
    return dict(F=F, Vmax=Vmax, R=R, kF=kF, kA=kA, kB=kB, kC=kC)


if __name__ == "__main__":
    a, b = 1 / 18, 1 / 3
    p = parts(a, b)
    print(f"Fauzan [a,b] = [1/18, 1/3]:  F(a,b) = {p['F']:.4f} vs max V = {p['Vmax']:.4f};  lam(1-log al lam) = "
          f"{LAM*(1-math.log(AL*LAM)):.4f} vs max phi = {phimax:.4f} (at s = {s[phis.argmax()]:.4f});  2 lam log R = {2*LAM*math.log(p['R']):.4f}")
    print(f"  kappa_F (Fauzan)                       = {p['kF']:.4f}")
    print(f"  A: exact max_[a,b] V                   = {p['kA']:.4f}   (gain {p['kF']-p['kA']:.4f})")
    print(f"  B: + exact |r_j| exponent max phi     = {p['kB']:.4f}   (gain {p['kA']-p['kB']:.4f})")
    print(f"  C: + coupling sup_s(2U(-s^2)+phi)     = {p['kC']:.4f}   (gain {p['kB']-p['kC']:.4f})")
    # D: optimize a, b for kC
    f = lambda v: parts(math.exp(v[0]), math.exp(v[0]) + math.exp(v[1]))['kC']
    best = minimize(f, [math.log(0.01), math.log(0.5)], method='Nelder-Mead', options=dict(xatol=1e-6, fatol=1e-8, maxiter=4000))
    aD, bD = math.exp(best.x[0]), math.exp(best.x[0]) + math.exp(best.x[1])
    print(f"  D: best single arcsine [a,b] = [{aD:.5f},{bD:.5f}] -> {best.fun:.4f}   (gain {p['kC']-best.fun:.4f})")
    r = ex.equilibrium_bound(8192)
    print(f"  E: equilibrium measure (sharp)        = {r['rho']:.6f}   (gain {best.fun-r['rho']:.4f})")
    print(f"  F: certified 1600-arcsine measure     = {75299683/10000000:.6f}   (loss vs E: {75299683/10000000 - r['rho']:.2e})")
    # also: Fauzan's own kappa optimized over (a,b) (09: 13.8635)
    g = lambda v: parts(math.exp(v[0]), math.exp(v[0]) + math.exp(v[1]))['kF']
    bF = minimize(g, [math.log(1 / 18), math.log(1 / 3 - 1 / 18)], method='Nelder-Mead', options=dict(xatol=1e-7, fatol=1e-9))
    print(f"  (Fauzan's kappa minimized over (a,b): {bF.fun:.4f} at a={math.exp(bF.x[0]):.5f}, b={math.exp(bF.x[0])+math.exp(bF.x[1]):.5f})")
    print(f"  equilibrium: support [{r['a']:.6g},{r['b']:.6f}], F = {r['F_lo']:.8f}, M1 = {r['M1']:.6f} at s* = {r['smax']:.5f}")
