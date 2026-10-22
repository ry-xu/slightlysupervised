"""18-redteam-mu96: negative controls - each certificate must FAIL when its target is moved past the float extremum."""
import sys, bisect
from flint import fmpq, arb, acb, ctx
sys.path.insert(0, './scripts/18-redteam-mu96')
from common import load_measure
from certreal import RealCert, certify, breakpoints
from certcontour import ContourCert, certify_contour
ctx.prec = 128
A, B, C = load_measure()
rc = RealCert(A, B, C)
quiet = lambda *a: None
pts = breakpoints(rc, fmpq(0), fmpq(2))
# sup: gap containing float argmax t = 0.3359919
i = bisect.bisect(pts, fmpq(335991898628, 10 ** 12))
l, r = pts[i - 1], pts[i]
for tgt in [fmpq(-66607669, 10 ** 7), fmpq(-666076682, 10 ** 8), fmpq(-6660766817, 10 ** 9)]:
    ok, n, w = certify(rc, l, r, tgt, 'sup', log=quiet, maxdepth=40)
    print("SUP negative control, gap [%.9f,%.9f], target %.10f (float sup -6.66076681785): certified=%s (expect False)" % (
        float(l), float(r), float(tgt), ok))
ok, n, w = certify(rc, l, r, fmpq(-66607665, 10 ** 7), 'sup', log=quiet)
print("SUP positive control same gap, target -6.6607665: certified=%s pieces=%d" % (ok, n))
# inf: gaps adjacent to breakpoint 0.334225221 (float argmin over S)
j = bisect.bisect(pts, fmpq(334225221, 10 ** 9))
l, r = pts[j - 2], pts[j]
for tgt in [fmpq(-66608135, 10 ** 7), fmpq(-666081359, 10 ** 8)]:
    ok, n, w = certify(rc, l, r, tgt, 'inf', log=quiet, maxdepth=40)
    print("INF negative control, [%.9f,%.9f], target %.10f (float inf -6.66081360140): certified=%s (expect False)" % (
        float(l), float(r), float(tgt), ok))
ok, n, w = certify(rc, l, r, fmpq(-6660814, 10 ** 6), 'inf', log=quiet)
print("INF positive control same range, target -6.660814: certified=%s pieces=%d" % (ok, n))
# contour: H(x1) = 0.80880185128
cc = ContourCert(A, B, C)
for tgt in [fmpq(8088018, 10 ** 7), fmpq(80880185, 10 ** 8)]:
    ok, n, w = certify_contour(cc, tgt, log=quiet, maxdepth=30)
    print("CONTOUR negative control target %.10f (H(x1)=0.80880185128): certified=%s (expect False)" % (float(tgt), ok))
