import sys
sys.path.insert(0, './scripts/18-redteam-mu96')
from true_norm_mp import run
for K in [int(a) for a in sys.argv[1:]]:
    res = run(K, n=96, fine=1.5)
    print("K=%d nodes=%d (mp Stieltjes, fine) : log||B_K|| = %.10f ; rate = %.8f ; orth %s [%.0fs]" % (
        K, res['nodes'], res['lognorm'], res['lognorm'] / K, {str(k): '%.1e' % v for k, v in res['diag'].items()}, res['t']), flush=True)
