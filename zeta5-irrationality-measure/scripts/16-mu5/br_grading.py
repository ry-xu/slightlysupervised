"""16-mu5: try mass gradings for the Buyarov-Rakhmanov discretization (float exploration)."""
import sys, json
import numpy as np
sys.path.insert(0, './scripts/16-mu5')
import br_measure as br
import evalmeasure as ev

gradings = {
    'smooth3': lambda u: 3 * u ** 2 - 2 * u ** 3,
    'sq0': lambda u: u ** 2,
    'cub0': lambda u: u ** 3,
    'mix': lambda u: 0.5 * (u ** 3) + 0.5 * (3 * u ** 2 - 2 * u ** 3),
}
J = int(sys.argv[1])
names = sys.argv[2].split(',')
for name in names:
    d = br.br_measure(J, grading=gradings[name])
    r = ev.evaluate(d['a'], d['b'], d['c'])
    json.dump(d, open(f'./data/16-mu5/br_{name}_{J}.json', 'w'))
    print(f"{name} J={J}: inner [{d['a'][0]:.6g},{d['b'][0]:.6g}] c1={d['c'][0]:.2e} outer [{d['a'][-1]:.4g},{d['b'][-1]:.6f}]"
          f"  U={r['U']:.7f} (sup at t={r['t_sup']:.4g}, M0={r['M0']:.7f})  rho={r['rho']:.6f} (M1={r['M1']:.6f}, m0={r['m0']:.6f} at {r['t_m0']:.4g})", flush=True)
