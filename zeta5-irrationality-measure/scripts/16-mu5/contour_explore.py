"""16-mu5: exploration of the CONTOUR (sign-aware) bound for ||B_K||.  [float, exploration only]

q^T M_K q = sum_j Res_{t=-j^2} [t^2 Omega(t) q(t)^2] = (1/2 pi i) oint_Gamma t^2 Omega(t) q(t)^2 dt
for any closed Gamma enclosing all poles -j^2, N < j <= K (and nothing else is singular).  After scaling, the
integrand has exponential rate  H(tau) = 2U^rho(tau) - Vt(tau)  on Gamma (Bernstein-Walsh for |p|^2, continuum limit
of |Omega|), with Vt(tau) = int_0^1 log|tau+u^2| du - 6 int_0^al log|tau+u^2| du.  So
      rate <= min_Gamma max_{tau in Gamma} H(tau)  -  inf_S (2U^rho - V).
Gamma must separate the pole segment [-1,-al^2] from oo: it crosses R once left of -1 and once right of -al^2.
This script computes the minimax ("mountain pass") value on a grid with a bottleneck Dijkstra.
"""
import sys, json, heapq
import numpy as np
sys.path.insert(0, './scripts/16-mu5')
sys.path.insert(0, './scripts/09-zeta7det')
import equilibrium as eq

AL, LAM = 3 / 40, 37 / 40


def I_log(a, c):
    """int_0^a log|u - c| du for complex array c"""
    c = np.asarray(c, complex)
    with np.errstate(divide='ignore', invalid='ignore'):
        t1 = np.where(np.abs(a - c) > 0, (a - c) * np.log(a - c + 0j), 0)
        t2 = np.where(np.abs(c) > 0, c * np.log(-c + 0j), 0)
    return np.real(t1 + t2) - a


def int_log_abs(a, tau):
    """int_0^a log|tau + u^2| du (complex tau)"""
    r = np.sqrt(np.asarray(tau, complex))
    return I_log(a, 1j * r) + I_log(a, -1j * r)


def Vt(tau):
    return int_log_abs(1.0, tau) - 6 * int_log_abs(AL, tau)


def equilibrium_nodes(n=8192):
    prof = eq.fauzan_profile(6, AL)
    a, b, ier, resid = eq.solve_support(prof, LAM, n=n, guess=(1e-4, 0.75))
    E, mass, lo, hi, w = eq.energy(a, b, prof, LAM, n=n)
    c, d = (a + b) / 2, (b - a) / 2
    j = np.arange(1, n)
    tn = c + d * np.cos(np.pi * j / n)
    return tn, w, 0.5 * (lo + hi)


def U_nodes(z, tn, w):
    z = np.asarray(z, complex)
    out = np.zeros(z.shape)
    for t, ww in zip(tn, w):
        out += ww * np.log(np.abs(z - t))
    return out


def minimax(xs, ys, Hgrid, left_ok, right_ok, blocked):
    """bottleneck path on grid (4-neighbour) from any left_ok node to any right_ok node avoiding blocked."""
    nx, ny = len(xs), len(ys)
    best = np.full((nx, ny), np.inf)
    pq = []
    for i in range(nx):
        for k in range(ny):
            if left_ok[i, k] and not blocked[i, k]:
                best[i, k] = Hgrid[i, k]
                heapq.heappush(pq, (Hgrid[i, k], i, k))
    while pq:
        v, i, k = heapq.heappop(pq)
        if v > best[i, k]:
            continue
        if right_ok[i, k]:
            return v, (xs[i], ys[k]), best
        for di, dk in ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)):
            ii, kk = i + di, k + dk
            if 0 <= ii < nx and 0 <= kk < ny and not blocked[ii, kk]:
                nv = max(v, Hgrid[ii, kk])
                if nv < best[ii, kk]:
                    best[ii, kk] = nv
                    heapq.heappush(pq, (nv, ii, kk))
    return np.inf, None, best


if __name__ == "__main__":
    tn, w, F = equilibrium_nodes()
    xs = np.linspace(-2.0, 0.3, 461)
    ys = np.linspace(0.0, 1.2, 241)
    X, Y = np.meshgrid(xs, ys, indexing='ij')
    Z = X + 1j * Y
    H = 2 * U_nodes(Z, tn, w) - Vt(Z)
    # reference: values on the segment
    s = np.linspace(AL, 1, 2001)
    Hseg = 2 * U_nodes(-s * s + 0j, tn, w) - Vt(-s * s + 0j)
    print(f"F = {F:.6f};  max_segment H = {Hseg.max():.6f} at s = {s[Hseg.argmax()]:.4f}  -> segment rate {Hseg.max()-F:.6f}")
    blocked = (np.abs(Y) < 1e-12) & (X >= -1.0) & (X <= -AL * AL)
    left_ok = (np.abs(Y) < 1e-12) & (X < -1.0)
    right_ok = (np.abs(Y) < 1e-12) & (X > -AL * AL)
    v, where, best = minimax(xs, ys, H, left_ok, right_ok, blocked)
    print(f"minimax over contours (grid): {v:.6f}  ->  contour rate {v - F:.6f}   (reached right axis at {where})")
    # where is the bottleneck? find the max along the optimal path region: report H at a few points
    i_min = np.unravel_index(np.argmin(np.where(left_ok, H, np.inf)), H.shape)
    print("min H on left axis (x<-1):", H[i_min], "at x =", xs[i_min[0]])
    for yy in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5):
        k = np.argmin(np.abs(ys - yy))
        row = H[:, k]
        msk = (xs >= -1.0) & (xs <= -AL * AL)
        print(f"  y = {ys[k]:.3f}: max_x in [-1,-al^2] H = {row[msk].max():.5f} at x = {xs[msk][row[msk].argmax()]:.3f};  min over x<0 = {row[xs<0].min():.5f}")
    json.dump(dict(F=F, segment=float(Hseg.max()), contour=float(v)), open('./data/16-mu5/contour_explore.json', 'w'))
