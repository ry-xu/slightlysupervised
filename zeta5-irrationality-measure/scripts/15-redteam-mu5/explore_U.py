"""15-redteam-mu5: float exploration (NOT a certificate) of 09's measure and Fauzan's measure."""
import sys, math
import numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from measure_common import *

def Urho_vec(t, a, b, c):
    t = np.asarray(t, dtype=float)
    out = np.zeros_like(t)
    for aj, bj, cj in zip(a, b, c):
        aj, bj, cj = float(aj), float(bj), float(cj)
        m = (aj + bj) / 2
        inside = (t >= aj) & (t <= bj)
        val = np.empty_like(t)
        val[inside] = math.log((bj - aj) / 4)
        o = ~inside
        val[o] = np.log((np.abs(t[o] - m) + np.sqrt((t[o] - aj) * (t[o] - bj))) / 2)
        out += cj * val
    return out

def V_vec(t):
    al = 3 / 40
    t = np.asarray(t, dtype=float)
    st = np.sqrt(t)
    with np.errstate(divide='ignore', invalid='ignore'):
        i1 = np.log1p(t) - 2 + 2 * st * np.arctan(1 / st)
        i2 = al * np.log(t + al * al) - 2 * al + 2 * st * np.arctan(al / st)
    v = 2 * math.pi * st + i1 - 6 * i2
    v = np.where(t == 0, -2 - 12 * al * math.log(al) + 12 * al, v)
    return v

def energy_A2(a, b, c):
    S = 0.0; I = 0.0
    for aj, bj, cj in zip(a, b, c):
        S1 = S + float(cj)
        I += (S1 * S1 - S * S) * math.log(float(bj - aj) / 4)
        S = S1
    return I

def energy_quad(a, b, c, n=4000):
    """I(rho) = sum_i c_i int U^rho d omega_i, by Gauss-Chebyshev quadrature (independent of (A.2))."""
    k = np.arange(1, n + 1)
    cosk = np.cos((2 * k - 1) * math.pi / (2 * n))
    I = 0.0
    for ai, bi, ci in zip(a, b, c):
        m = float(ai + bi) / 2; r = float(bi - ai) / 2
        pts = m + r * cosk
        I += float(ci) * Urho_vec(pts, a, b, c).mean()
    return I

def Cstar():
    lam, al = 37 / 40, 3 / 40
    return -2 * lam + 12 * al * lam * (1 - math.log(al)) + 3 * lam ** 2 - 2 * lam ** 2 * math.log(2 * lam)

for name, (a, b, c) in [('fauzan16', load_fauzan()), ('09-3p', load_09('3p'))]:
    print(f"== {name}")
    admissibility(a, b, c)
    # dense grid on [0,2] plus extra density near the support, plus all endpoints
    grid = np.unique(np.concatenate([np.linspace(0, 2, 200001), np.geomspace(1e-9, 0.8, 400001),
                                     np.array([float(x) for x in a + b])]))
    f = 2 * Urho_vec(grid, a, b, c) - V_vec(grid)
    i = np.argmax(f)
    print(f"   float sup(2U-V) on grid = {f[i]:.10f} at t = {grid[i]:.6g}")
    Ia = energy_A2(a, b, c); Iq = energy_quad(a, b, c)
    print(f"   I(rho): (A.2) = {Ia:.12f}   quadrature = {Iq:.12f}   diff = {Ia-Iq:.2e}")
    lam = 37 / 40
    print(f"   C* = {Cstar():.12f};  U(float) = lam*sup - I + C* = {lam*f[i] - Ia + Cstar():.10f}")
    # where is 2U-V within 1e-5 of its sup?
    near = grid[f > f[i] - 1e-5]
    print(f"   points within 1e-5 of sup: {len(near)}, range [{near.min():.4g}, {near.max():.4g}]")
    # 'flatness': E(rho)=I - int V drho  vs J(rho) = lam*sup - I
