"""16-mu5: stability check of true_norm.py (finer quadrature, higher precision) at one K."""
import sys, json, time
import mpmath as mp
sys.path.insert(0, './scripts/16-mu5')
import true_norm as tn
K = int(sys.argv[1])
mp.mp.dps = 60
panels = [mp.mpf(i) / 800 for i in range(0, 81)] + [mp.mpf(1) / 10 + mp.mpf(i) / 100 for i in range(1, 141)]
t0 = time.time()
r = tn.run(K, dps=60, panels=panels, nper=64, eig=False)
r['time'] = round(time.time() - t0, 1); r['note'] = 'dps60, 220 panels x 64 nodes'
print(json.dumps(r))
