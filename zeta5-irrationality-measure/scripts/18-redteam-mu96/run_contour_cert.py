"""18-redteam-mu96: certify max_{Gamma} H < M1c (16's contour), plus point value at x1 and the true saddle."""
import sys, time
from flint import arb, acb, fmpq, ctx
sys.path.insert(0, './scripts/18-redteam-mu96')
from common import load_measure, ROOT, FloatMeasure, Vt_float
from certcontour import ContourCert, certify_contour, X1, Y
ctx.prec = 128
logf = open(ROOT + "/data/18-redteam-mu96/contour_cert.log", "a")
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); logf.write(s + "\n"); logf.flush()
target = fmpq(808802, 10 ** 6)
if len(sys.argv) > 2:
    target = fmpq(int(sys.argv[1]), int(sys.argv[2]))
log("=" * 100)
log("contour certificate %s ; target M1c = %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), target))
A, B, C = load_measure()
cc = ContourCert(A, B, C)
hx1 = cc.Hc(acb(arb(X1), 0), 'left')
log("H(x1) = %s, x1 = %s" % (hx1, X1))
# distance of Gamma from [-1,0]: left edge Re = x1 < -1, top edge Im = 1/5, right edge Re = 1
log("d_Gamma = min(|x1+1|, 1/5, 1) = %s" % float(min(-1 - X1, Y, fmpq(1))))
t0 = time.time()
ok, n, worst = certify_contour(cc, target, log=log)
log("RESULT max_Gamma+ H < %s : %s ; pieces %d ; worst upper bound %s ; F' evals %d [%.1fs]" % (
    target, ok, n, worst, cc.nevals, time.time() - t0))
