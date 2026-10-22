"""18-redteam-mu96: common definitions (own code; 16's code only read for the data format).

Conventions (Fauzan, Zenodo 10.5281/zenodo.22826419):
  alpha = 3/40, lambda = 37/40, K in 40Z, N = 3K/40, h = 37K/40.
  V(t) = 2 pi sqrt t + int_0^1 log(t+u^2) du - 6 int_0^alpha log(t+u^2) du        (6.1)
  U^rho(z) = int log|z-u| d rho(u)                                                   (Fauzan's sign)
  Vt(tau) = int_0^1 log|tau+u^2| du - 6 int_0^alpha log|tau+u^2| du  (complex tau; V = 2 pi sqrt t + Vt on t>0)

Measure file data/16-mu5/rho1600.exact.json: lists a, b, c of exact rationals, index 0 = innermost interval.
"""
import json
import math
import numpy as np
from flint import fmpq, fmpz, arb, acb, ctx

ROOT = "."
ALPHA = fmpq(3, 40)
LAM = fmpq(37, 40)


def load_measure(path=ROOT + "/data/16-mu5/rho1600.exact.json"):
    d = json.load(open(path))
    A = [fmpq(*map(int, s.split("/"))) if "/" in s else fmpq(int(s)) for s in d["a"]]
    B = [fmpq(*map(int, s.split("/"))) if "/" in s else fmpq(int(s)) for s in d["b"]]
    C = [fmpq(*map(int, s.split("/"))) if "/" in s else fmpq(int(s)) for s in d["c"]]
    return A, B, C


def admissibility(A, B, C):
    """Exact checks of the hypotheses of Fauzan's Lemma 6.1 proof (items 1,2,4,5)."""
    J = len(A)
    out = {}
    out["J"] = J
    out["all_c_positive"] = all(c > 0 for c in C)
    out["mass_exact"] = (sum(C, fmpq(0)) == LAM)
    out["mass"] = str(sum(C, fmpq(0)))
    # strict nesting with index 0 innermost: A[k+1] < A[k] < B[k] < B[k+1]
    nest = all(A[k + 1] < A[k] for k in range(J - 1)) and all(B[k] < B[k + 1] for k in range(J - 1)) \
        and all(A[k] < B[k] for k in range(J))
    out["strictly_nested"] = nest
    out["support"] = (float(A[-1]), float(B[-1]))
    out["support_in_(0,2)"] = (A[-1] > 0) and (B[-1] < 2)
    L = [B[k] - A[k] for k in range(J)]
    out["min_len"] = float(min(L))
    out["max_len"] = float(max(L))
    # regularization constant sum_j 8 c_j sqrt2 / (pi sqrt(L_j)), rigorous upper bound in arb
    old = ctx.prec
    ctx.prec = 128
    s = arb(0)
    for k in range(J):
        s += 8 * arb(C[k]) * arb(2).sqrt() / (arb.pi() * arb(L[k]).sqrt())
    out["reg_const_upper"] = float(s.upper())
    out["reg_const_le_60"] = bool(s < 60)
    # max number of components with length < 1/225 (Fauzan's original sufficient condition)
    out["n_len_below_1/225"] = sum(1 for x in L if x <= fmpq(1, 225))
    ctx.prec = old
    return out


def energy_A2(A, B, C):
    """I(rho) = sum_j (S_j^2 - S_{j-1}^2) log((b_j - a_j)/4), index 0 innermost (Fauzan (A.2))."""
    S = fmpq(0)
    I = arb(0)
    for k in range(len(A)):
        Sn = S + C[k]
        I += arb(Sn * Sn - S * S) * (arb(B[k] - A[k]) / 4).log()
        S = Sn
    return I


def Cstar():
    a = arb(ALPHA)
    l = arb(LAM)
    return -2 * l + 12 * a * l * (1 - a.log()) + 3 * l * l - 2 * l * l * (2 * l).log()


# ---------------------------------------------------------------- float (numpy) exploration only
class FloatMeasure:
    def __init__(self, A, B, C):
        self.a = np.array([float(x) for x in A])
        self.b = np.array([float(x) for x in B])
        self.c = np.array([float(x) for x in C])
        self.m = 0.5 * (self.a + self.b)
        self.r = 0.5 * (self.b - self.a)

    def U(self, t):
        """U^rho(t), t real (vectorised over components; scalar t)."""
        a, b, c, m, r = self.a, self.b, self.c, self.m, self.r
        out = np.where((t >= a) & (t <= b), np.log(r / 2.0),
                       np.log((np.abs(t - m) + np.sqrt(np.abs((t - a) * (t - b)))) / 2.0))
        return float(np.dot(c, out))

    def Uc(self, z):
        """U^rho(z) for complex z (any z)."""
        w = z - self.m
        ws = np.where(w.real >= 0, w, -w)
        P = np.sqrt(ws - self.r) * np.sqrt(ws + self.r)
        return float(np.dot(self.c, np.log(np.abs(ws + P) / 2.0)))


def V_float(t):
    al = 3 / 40
    if t == 0:
        return -12 * al * math.log(al) - 2 + 12 * al
    st = math.sqrt(t)
    return math.log(1 + t) - 6 * al * math.log(t + al * al) - 2 + 12 * al + 2 * st * (
        math.pi + math.atan(1 / st) - 6 * math.atan(al / st))


def J_float(a, tau):
    """int_0^a log|tau+u^2| du for complex tau not on [-a^2,0], via s = sqrt(-tau) or c = i sqrt(tau)."""
    tau = complex(tau)
    if tau.real < 0:
        s = np.sqrt(-tau)  # principal, Re s > 0
    else:
        s = 1j * np.sqrt(tau)  # Im s = Re sqrt(tau) > 0 unless tau on (-inf,0]
    val = (s + a) * np.log(s + a) - (s - a) * np.log(s - a)
    return float(val.real) - 2 * a


def Vt_float(tau):
    return J_float(1.0, tau) - 6 * J_float(3 / 40, tau)
