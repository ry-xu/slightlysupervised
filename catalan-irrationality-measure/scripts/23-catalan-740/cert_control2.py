"""Control: certify Litsey's StrongTrial caps with the angle-variable arb certifier.
mode pass: Litsey's caps (must succeed); mode fail: caps lowered by 1e-8 (must fail).  [rigorous]"""
import sys, json, time
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from certify2 import *
SP = '/private/tmp/claude-501/-Users-ryanxu-Documents-dev/313d84c6-3fcb-46b6-9f6f-56b316f7c448/scratchpad/'
d = json.load(open(SP + 'litsey/strongtrial.json'))
p = Tr(d['Pc'], d['P']); v = Tr(d['Vc'], d['V'])
mode = sys.argv[1]
t0 = time.time()
MX = float(Fr(d['MX'])); MY = float(Fr(d['MY'])); RM = -131 / 2500
if mode == 'fail':
    MX -= 1e-8; MY -= 1e-8; RM -= 1e-5
ok, res = cert_all(p, v, MX, MY, RM)
print('ALL OK' if ok else 'FAILED', round(time.time() - t0, 1), 's')
