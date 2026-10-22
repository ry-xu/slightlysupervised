"""
Export an optimised theta (Param layout) to exact data in Litsey's format and choose caps.  [exact + NUMERICAL]
  u_k = cs_{k-1}/10^8 + Re sum r z^k,  cs in Z, z,r in Q+Qi.
usage: export_trial.py in.npz out.json
Writes: Pc, P, Vc, V (exact), float diagnostics, proposed caps MX, MY, RM (float max + margins).
"""
import sys, json
import numpy as np
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from trial import Trial, dual_value, xgrid, sgrid, Mfield, LOG2
from optdual import Param


def frac_of(xf, den=10**10):
    return Fr(int(round(xf * den)), den)


def zfrac(z):
    """exact rational base: Litsey's bases are short decimals; recover them."""
    return Fr(str(round(z.real, 12))).limit_denominator(10**7), Fr(str(round(z.imag, 12))).limit_denominator(10**7)


def export(npz):
    zf = np.load(npz, allow_pickle=True)
    th = zf['theta']; Kp = int(zf['Kp']); Kv = int(zf['Kv'])
    tp = [complex(z) for z in zf['tp']]; tv = [complex(z) for z in zf['tv']]
    par = Param(Kp, tp, Kv, tv)
    p, v = par.to_trials(th)
    def pack(tr):
        cs = [int(round(c * 1e8)) for c in tr.c]
        tails = []
        for z, r in tr.tails:
            zr, zi = zfrac(z)
            rr, ri = frac_of(r.real), frac_of(r.imag)
            tails.append([str(zr), str(zi), str(rr), str(ri)])
        return cs, tails
    Pc, P = pack(p); Vc, V = pack(v)
    return Pc, P, Vc, V


if __name__ == "__main__":
    Pc, P, Vc, V = export(sys.argv[1])
    p = Trial.from_litsey(Pc, P); v = Trial.from_litsey(Vc, V)
    xg = xgrid(800001); sg = sgrid(800001)
    out = dual_value(p, v, xg, sg)
    sel = np.abs(xg) >= 401 / 2500
    out['R'] = float(Mfield(p, xg[sel]).max())
    out['U'] = -out['dual']
    print(json.dumps({k: float(val) for k, val in out.items()}))
    d = dict(Pc=Pc, P=P, Vc=Vc, V=V, float_diag={k: float(val) for k, val in out.items()})
    json.dump(d, open(sys.argv[2], 'w'))
