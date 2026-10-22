import sys
sys.path.insert(0, './scripts/18-redteam-mu96')
from true_norm_mp import run
for spec in sys.argv[1:]:
    K, n, fine = spec.split(':'); K = int(K); n = int(n); fine = float(fine)
    res = run(K, n=n, fine=fine, dps=60)
    print("K=%d n=%d fine=%.1f nodes=%d (mp Stieltjes dps60) : log||B_K|| = %.10f ; rate = %.8f ; orth %s [%.0fs]" % (
        K, n, fine, res['nodes'], res['lognorm'], res['lognorm'] / K, {str(k): '%.1e' % v for k, v in res['diag'].items()}, res['t']), flush=True)
