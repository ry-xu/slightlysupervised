"""16-mu5: Buyarov-Rakhmanov nested-arcsine test measures.

For the field V of Fauzan (6.1), let S_m = [a(m), b(m)] be the support of the equilibrium measure of mass m
(maximize I(s) - int V ds over s >= 0 of mass m).  Buyarov-Rakhmanov: d sigma_m / dm = omega_{S_m} (arcsine of S_m),
so sigma_lambda = int_0^lambda omega_{S_m} dm.  A midpoint discretization gives a nested arcsine measure
    rho_J = sum_j (lambda/J) omega_{S_{m_j}},   m_j = (j - 1/2) lambda / J,
which is admissible for Fauzan's Lemma 6.1 scheme (nested intervals) and for the norm bound of notes/16-mu5.md.
This script computes S_m (float, Chebyshev soft-edge solver of 09) and writes the measures to data/16-mu5/.
"""
import sys, json, math
import numpy as np
sys.path.insert(0, './scripts/09-zeta7det')
import equilibrium as eq

AL, LAM = 3 / 40, 37 / 40
PROF = eq.fauzan_profile(6, AL)


def support(m, guess, n=4096):
    a, b, ier, resid = eq.solve_support(PROF, m, n=n, guess=guess)
    return a, b, ier, resid


T_STAR = 0.005920507796065   # argmin V (bracket of 15-redteam)
def _Vpp(t, h=1e-7):
    return (eq.V(t + h, PROF) - 2 * eq.V(t, PROF) + eq.V(t - h, PROF)) / h ** 2
VPP = float(_Vpp(T_STAR))


def _ok(a, b, resid):
    return (0 < a < b) and max(abs(r) for r in resid) < 1e-9


def trace_supports(ms, n=4096):
    """continuation from large m to small m; fallback semicircle guess r = 2 sqrt(m/V'') (quadratic field)"""
    out = {}
    guess = (7.12e-5, 0.754)
    for m in sorted(ms, reverse=True):
        r = 2 * (m / VPP) ** 0.5
        sc = (max(T_STAR - r, 1e-7), T_STAR + r)
        tries = [guess, sc, (guess[0] * 0.9, guess[1] / 0.9), (guess[0] / 0.9, guess[1] * 0.9)]
        for gs in tries:
            a, b, ier, resid = support(m, gs, n=n)
            if _ok(a, b, resid):
                break
        else:
            raise RuntimeError(f"support solve failed for m = {m}: {a}, {b}, {resid}")
        out[m] = (a, b)
        guess = (a, b)
    return out


def br_measure(J, n=4096, grading=None):
    """grading: optional function mapping uniform u in (0,1) to mass fraction (to refine near m=lambda)"""
    if grading is None:
        edges = np.linspace(0, LAM, J + 1)
    else:
        edges = LAM * grading(np.linspace(0, 1, J + 1))
    mids = 0.5 * (edges[1:] + edges[:-1])
    w = np.diff(edges)
    S = trace_supports(list(mids), n=n)
    a = [S[m][0] for m in mids]
    b = [S[m][1] for m in mids]
    return dict(a=a, b=b, c=list(w), m=list(mids))


if __name__ == "__main__":
    J = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    g = sys.argv[2] if len(sys.argv) > 2 else 'uniform'
    if g == 'uniform':
        grading = None
    elif g == 'sqrt':
        grading = lambda u: 1 - (1 - u) ** 2     # finer near m = lambda (outer intervals)
    else:
        raise ValueError(g)
    d = br_measure(J, grading=grading)
    fn = f'./data/16-mu5/br_{g}_{J}.json'
    json.dump(d, open(fn, 'w'))
    print("wrote", fn, "innermost", d['a'][0], d['b'][0], "outermost", d['a'][-1], d['b'][-1])
