"""16-mu5: [HEURISTIC] how the irrationality-measure exponent of Fauzan's method depends on the zero fraction alpha
(r = 6 zeros, lambda = 1 - alpha).  For each alpha:
   rho(alpha) = min over contours of max H - F   (Theorem 2 constant with the equilibrium measure; grid minimax)
   eps(alpha) = -(A_inf + U*)                    (09's generalized bookkeeping model, bookkeeping.A_inf, + U*)
   mu(alpha)  = rho * lambda / eps.
Only alpha = 3/40 is backed by Fauzan's verified arithmetic; other alpha would need Secs. 3-5 re-proved."""
import sys, json, math
import numpy as np
sys.path.insert(0, './scripts/16-mu5')
sys.path.insert(0, './scripts/09-zeta7det')
import equilibrium as eq
from bookkeeping import A_inf
import contour_explore as ce


def run(al, r=6, n=4096):
    lam = 1 - al
    prof = eq.fauzan_profile(r, al)
    a, b, ier, resid = eq.solve_support(prof, lam, n=n, guess=(1e-4, 0.8))
    E, mass, lo, hi, w = eq.energy(a, b, prof, lam, n=n)
    U = E + eq.Cstar(r, al, lam) if r == 6 else None
    c, d = (a + b) / 2, (b - a) / 2
    j = np.arange(1, n)
    tn = c + d * np.cos(np.pi * j / n)
    F = 0.5 * (lo + hi)
    Vt = lambda z: ce.int_log_abs(1.0, z) - r * ce.int_log_abs(al, z)
    xs = np.linspace(-2.0, 0.3, 231)
    ys = np.linspace(0.0, 1.0, 101)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    Z = X + 1j * Y
    H = 2 * ce.U_nodes(Z, tn, w) - Vt(Z)
    blocked = (np.abs(Y) < 1e-12) & (X >= -1.0) & (X <= -al * al)
    left_ok = (np.abs(Y) < 1e-12) & (X < -1.0)
    right_ok = (np.abs(Y) < 1e-12) & (X > -al * al)
    v, where, best = ce.minimax(xs, ys, H, left_ok, right_ok, blocked)
    # refine the real-axis minimum left of -1
    x = np.linspace(-1.6, -1.00001, 6001)
    hx = 2 * ce.U_nodes(x + 0j, tn, w) - Vt(x + 0j)
    s = np.linspace(al, 1, 2001)
    hs = 2 * ce.U_nodes(-s * s + 0j, tn, w) - Vt(-s * s + 0j)
    A5, _ = A_inf(5, r, al, lam, X=400.0, nout=6000)
    eps = -(A5 + U)
    rho_c = max(v, hx.min()) - F
    rho_seg = hs.max() - F
    return dict(al=al, lam=lam, U=U, A=A5, eps=eps, F=F, rho_contour=rho_c, rho_segment=rho_seg,
                mu_contour=rho_c * lam / eps if eps > 0 else None, mu_segment=rho_seg * lam / eps if eps > 0 else None)


if __name__ == "__main__":
    for al in [float(x) for x in sys.argv[1].split(',')]:
        print(json.dumps(run(al)), flush=True)
