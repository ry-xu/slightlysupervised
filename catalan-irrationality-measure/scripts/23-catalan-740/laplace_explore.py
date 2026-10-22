"""
Explore whether the Laplace (reflected-row) sheet bound improves OpenAI's kappa=2 energy.

Functional family (theta in [0,1]):  theta * Phi_2(sigma+tau) + (1-theta) * Psi(sigma,tau)
  Phi_2 : OpenAI kappa=2 majorant functional (005 (78)/(80), primal form)
  Psi   : Laplace bound  |R_n| <= max_S |det_diff(S)| * |det_rem(S^c)|,
          sigma = S-nodes (positive, mass g/n), tau = the rest (mass 1-g/n).
usage: laplace_explore.py theta Mx Ms [design n0 b0 q0 g0 h0]      [NUMERICAL]
"""
import sys, json, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from energy import *


def build(theta, Mx, Ms, design=(20, 3, 1, 2, 1)):
    n0, b0, q0, g0, h0 = design
    al, be, gg, gq, et = (b0 + q0) / n0, b0 / n0, g0 / n0, q0 / n0, h0 / n0
    ex, exp_, es = grids(Mx, Ms)
    W2 = lambda x: (al + 2 * gg) * np.log(np.abs(x)) + 2 * et * np.log(1 - x) - (1 + al + gg + et) * np.log(1 + x * x)
    Ws = lambda s: be * np.log(s) + gq * np.log(1 - s)
    # populations: 0 = sigma on (0,1), 1 = tau on (-1,1), 2 = nu on (0,1)
    sizes = [len(exp_) - 1, len(ex) - 1, len(es) - 1]
    masses = [gg, 1 - gg, 1.0]
    P = QP(sizes, masses)
    P.const = (1 + al + gg) * np.log(2)
    fsig = lambda x: W2(x) + (1 - theta) * (1 - 2 * gg) * np.log(x)
    P.add_lin(0, cell_avg_1d(exp_, fsig))
    P.add_lin(1, cell_avg_1d(ex, W2))
    P.add_lin(2, cell_avg_1d(es, Ws))
    Kss = loglog_cells(exp_, exp_); Ktt = loglog_cells(ex, ex); Kst = loglog_cells(exp_, ex)
    Jss = smooth_cells(exp_, exp_, J_ker); Jtt = smooth_cells(ex, ex, J_ker); Jst = smooth_cells(exp_, ex, J_ker)
    Xs = smooth_cells(exp_, es, XS_ker); Xt = smooth_cells(ex, es, XS_ker)
    Knn = loglog_cells(es, es)
    P.add_quad(0, 0, Kss, 1.0)
    P.add_quad(1, 1, Ktt, 1.0)
    P.add_quad(0, 1, Kst, 1.0 + theta)
    P.add_quad(0, 0, Jss, theta / 2 + (1 - theta))
    P.add_quad(1, 1, Jtt, 0.5)
    P.add_quad(0, 1, Jst, 1.0)
    P.add_quad(2, 2, Knn, 1.0)
    P.add_quad(0, 2, Xs, -1.0)
    P.add_quad(1, 2, Xt, -1.0)
    return P, (ex, exp_, es)


if __name__ == "__main__":
    theta = float(sys.argv[1])
    Mx = int(sys.argv[2]); Ms = int(sys.argv[3])
    design = tuple(int(v) for v in sys.argv[4:9]) if len(sys.argv) > 8 else (20, 3, 1, 2, 1)
    t0 = time.time()
    P, (ex, exp_, es) = build(theta, Mx, Ms, design)
    val, z = P.solve(iters=3000, step=0.5, polish=True)
    sig = z[P.blk(0)]; tau = z[P.blk(1)]
    xm = (ex[:-1] + ex[1:]) / 2; xpm = (exp_[:-1] + exp_[1:]) / 2
    out = dict(theta=theta, Mx=Mx, Ms=Ms, design=design, Phi=float(val), physical=float(val + 0.5 * np.log(2)),
               sigma_supp=[float(xpm[sig > 1e-7].min()), float(xpm[sig > 1e-7].max())],
               tau_neg_mass=float(tau[xm < 0].sum()), time=time.time() - t0)
    print(json.dumps(out), flush=True)
