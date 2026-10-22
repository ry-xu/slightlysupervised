# μ(G) ≤ 727.03 for Catalan's constant

**Claim.** For every integer a and all sufficiently large b, |G − a/b| > b^(−727.03), where G is Catalan's constant
(exponent 18175647944007/25000000000 ≈ 727.026). The previous best bound is 740, from S. Litsey's Lean development
([slitsey/catalan-irrationality-measure](https://github.com/slitsey/catalan-irrationality-measure)). That repo still states
740 as of 23 Oct 2026.

## Status: please read first

- **Built on Litsey's development, with one input changed.** The only new ingredient is a better second (κ = 2) trial for
  the energy dual of his record design (20,3,1,2,1). Everything else is his argument and his ledger.
- **Not checked in Lean.** The new trial was certified in ball arithmetic (python-flint/arb) and the ledger was replayed
  exactly in rationals. It has **not** been run through Litsey's Lean build, so it's a candidate until his pipeline checks it.
- **AI-generated** by a Claude (Anthropic) agent, run and curated by Ryan Xu, who is not a specialist. Corrections welcome:
  ryan@ryanxu.net.

## The new trial

- 400 + 1000 finite modes, Litsey's 45 geometric tail bases, plus 11 new bases at the column-support edge frequencies. The data
  are exact rationals in Litsey's StrongTrial format: `data/23-catalan-740/trial_cand2.json`.
- The certified dual value is **2.2564477**, against 2.2563367 for Litsey's trial. The method's ceiling for this design is about
  2.2564510, so this route bottoms out near μ ≈ 726.66.
- Certified in arb/acb ball arithmetic: the row, column and marked-residual field suprema, both norms and p₁. The ledger
  inequalities are replayed exactly and re-checked with 50-digit mpmath.
- Side result: re-balancing Litsey's own certified constants, with no new trial, already gives 739.72.

Full notes, including routes that failed and why: `notes/23-catalan-740.md`.

## Reproduce

Run from this folder (Python 3.12, `python-flint`, `mpmath`, `numpy`, `scipy`).

Ledger (seconds):
```bash
python scripts/23-catalan-740/final_ledger.py 0.774561586 -0.97760061 -1.56820567 -0.0521786
```
It prints `ALL True` and `mu_exact 18175647944007/25000000000`. This was re-checked from this folder on 23 Oct 2026.

Field certificates (slow, several hours on one core; logs in `data/23-catalan-740/cert_cand2_*.log`):
```bash
python scripts/23-catalan-740/cert_new2.py data/23-catalan-740/trial_cand2.json X -0.97760061
python scripts/23-catalan-740/cert_new2.py data/23-catalan-740/trial_cand2.json Y -1.56820567
```

Controls: `control_litsey_pass.log` reproduces Litsey's 740 from his own trial, and the `control_*_fail.log` files show the
certifier rejecting deliberately broken inputs.

## Credits and third-party data

- S. Litsey: the Lean development, the ledger, and the trial format this builds on.
  `data/23-catalan-740/litsey_strongtrial.json` is his original trial, exported from his repository (Apache-2.0) and
  included only for the control runs.
- OpenAI, for the irrationality proof of G (family 005 in [openai/math](https://github.com/openai/math)).

## License

Text: CC BY 4.0. Code: MIT. Litsey's trial data: Apache-2.0, as in his repository.
