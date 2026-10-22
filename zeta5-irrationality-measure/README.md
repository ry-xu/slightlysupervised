# μ(ζ(5)) ≤ 95.91

**Claim.** The irrationality exponent of ζ(5) is at most 95.9083454190. For every ε > 0, all sufficiently large b and all integers a,
|ζ(5) − a/b| > b^(−95.91).

This improves the bound μ(ζ(5)) ≤ 260 that A. Fauzan obtained alongside his proof that ζ(5) is irrational
([Zenodo, 17 Sep 2026](https://doi.org/10.5281/zenodo.22826419)).

## Status: please read first

- **AI-generated.** The mathematics, the certified computations and the manuscript were produced by Claude (Anthropic) agents,
  run and curated by Ryan Xu.
- **Checked by AI, not yet by a human expert.** Two separate agents re-derived the argument from Fauzan's paper and wrote their
  own certification code. Both reproduce every certified constant. No human mathematician has checked the argument yet.
- **Builds on Fauzan's paper.** Several statements from Fauzan's paper are used unchanged. His irrationality proof has been
  formalized in Lean ([Romik](https://github.com/danromik/zeta5-irrationality), [mo271](https://github.com/mo271/zeta5)).
  The new arguments here are **not** formalized.
- Corrections and criticism are very welcome: ryan@ryanxu.net.

## What is new

Fauzan's bound comes from a perturbation argument whose cost is controlled by the norm of a certain matrix. The paper proves an
essentially sharp bound for that norm, log‖B_K‖ ≤ 7.469616·K + O(log K), in place of Fauzan's 13.866·K. It uses a contour
integral around the poles of the moment functional and a weighted Bernstein–Walsh inequality. It also sharpens Fauzan's decay
rate to ε > 0.0720416. Details are in `paper/`.

## Contents

| Path | What |
|---|---|
| `paper/mu-zeta5.tex`, `paper/mu-zeta5.md` | The manuscript (LaTeX and Markdown versions of the same text) |
| `scripts/16-mu5/` | Main certification. `mu_final16.py` certifies every constant and assembles the bound |
| `scripts/15-redteam-mu5/`, `scripts/18-redteam-mu96/` | Two independent re-verifications, each written separately |
| `scripts/09-zeta7det/` | Earlier work (the intermediate bound 181) and a third computation of the exact integral (`exact_A.py`) |
| `data/16-mu5/rho1600.exact.json` | The certified test measure. SHA-256 `0529b8fd…8474068f` |
| `data/15-redteam-mu5/int_20_4000.txt` | The exact integral I₂₀,₄₀₀₀. SHA-256 `84c23d27…5cce5f956` |
| `data/*/*.log` | Logs of the certification runs |

Some exploration scripts in these folders depend on intermediate data that is not included. They are kept for transparency. The
certification path is the one listed below.

## Reproduce

Python 3.12 with `python-flint==0.9.0`, `mpmath==1.3.0` and `cypari2` (see `requirements.txt`). Run from this folder:

```bash
python scripts/16-mu5/mu_final16.py
```

It takes about 2.5 minutes on one core, and its last line is `==> mu(zeta(5)) <= 95.91`. Its output matches
`data/16-mu5/mu_final16.log` apart from timings. This was re-checked from this folder on 23 Oct 2026.

Independent re-checks (Table 6 of the paper):

```bash
python scripts/18-redteam-mu96/run_real_certs.py
python scripts/18-redteam-mu96/run_contour_cert.py
python scripts/18-redteam-mu96/int_R_closed_form.py
python scripts/15-redteam-mu5/exact_R.py
```

## Credits

- A. Fauzan, for the irrationality proof and the framework this builds on.
- W. Zudilin's determinantal criterion (arXiv:1507.05697).
- D. Romik and mo271, for the Lean formalizations used to locate formal statements of Fauzan's results.

## License

Text: CC BY 4.0. Code: MIT.
