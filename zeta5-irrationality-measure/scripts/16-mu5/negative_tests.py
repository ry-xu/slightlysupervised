"""16-mu5: negative controls -- each certificate must FAIL when the target is beyond the float extremum."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from certlib import Measure
from certify import certify_sup, certify_inf, certify_neg

meas = Measure.from_json('./data/16-mu5/br_mix_800.json')
# float extrema (evalmeasure, 200k grid): sup f = -6.66074088737815, inf_S f = -6.660873212729079, M1 = 0.8691542147431187
for name, fn, target in [('sup', certify_sup, Fr('-6.6607410')), ('inf', certify_inf, Fr('-6.6608730')),
                         ('neg', certify_neg, Fr('0.8691541'))]:
    try:
        fn(meas, target, maxint=20000)
        print(f"{name}: UNEXPECTED SUCCESS at target {float(target)}  <-- BUG")
    except RuntimeError as e:
        print(f"{name}: failed as expected at target {float(target)}: {str(e)[:150]}")

# contour certificate: must fail below H(x1) = 0.80880185128...
from contour_certify import CMeasure, certify_path
meas16 = Measure.from_exact_json('./data/16-mu5/rho1600.exact.json')
cm = CMeasure(meas16)
try:
    certify_path(cm, Fr('-1.069911'), Fr(1, 5), Fr(1), Fr('0.8088018'), maxpieces=3000)
    print("contour: UNEXPECTED SUCCESS  <-- BUG")
except RuntimeError as e:
    print("contour: failed as expected:", str(e)[:160])
