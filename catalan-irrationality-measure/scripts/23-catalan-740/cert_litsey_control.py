"""Control: certify Litsey's own StrongTrial caps with my independent arb code (must pass),
and check that slightly smaller caps fail (negative control).  [rigorous ball arithmetic]"""
import sys, json, time
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from certify import *
SP = '/private/tmp/claude-501/-Users-ryanxu-Documents-dev/313d84c6-3fcb-46b6-9f6f-56b316f7c448/scratchpad/'
d = json.load(open(SP + 'litsey/strongtrial.json'))
p = ExactTrial(d['Pc'], d['P']); v = ExactTrial(d['Vc'], d['V'])
t0 = time.time()
NP = p.norm2(); NV = v.norm2()
print('NP', NP, ' Litsey', d['NP']); print('NV', NV, ' Litsey', d['NV'])
K = 2 * NP + NV / 2
print('2NP+NV/2 upper', upper(K), ' <= 3096999/4000000 ?', upper(K) <= 3096999 / 4000000)
mode = sys.argv[1] if len(sys.argv) > 1 else 'pass'
MX = float(Fr(d['MX'])); MY = float(Fr(d['MY']))
if mode == 'fail':
    MX -= 1e-8; MY -= 1e-8
ok, r = cert_trial(p, v, MX, MY, -131 / 2500 if mode == 'pass' else -131 / 2500 - 1e-5)
print('ALL OK' if ok else 'FAILED', time.time() - t0)
