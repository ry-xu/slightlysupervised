# The Carlitz analogue of "π + e is irrational", in characteristic 2

**Claim.** Let π̃ be the Carlitz period and e_C(1) the value at 1 of the Carlitz exponential for F₂[θ]. Then
1, π̃, e_C(1) and π̃·e_C(1) are linearly independent over F₂(θ). In particular π̃ + e_C(1), π̃·e_C(1) and π̃/e_C(1)
are not in F₂(θ).

This is the function-field ("toy world") version of the open problem of whether π + e, πe and π/e are irrational. For
q ≥ 3, Denis (1995) proved much more: π̃ and e_C(1) are algebraically independent. His argument needs a non-trivial
element of F_q^×, which does not exist when q = 2. As far as we could find, the q = 2 case was open.

## Status: please read first

- **AI-generated.** The proof and manuscript were produced by Claude (Anthropic) agents, run and curated by Ryan Xu, who is
  not a specialist.
- **The proof is elementary** (a Liouville-type argument), and it was checked line by line by a second agent. Every series
  identity and degree count is also checked by exact arithmetic in two independent sets of scripts.
- **Not yet checked by a human expert, and novelty not confirmed with specialists.** The main risk is that a general
  theorem in the literature already covers q = 2. Corrections welcome: ryan@ryanxu.net.

## Also in the paper

- Why Denis's method fails at q = 2, and which single ingredient is missing in characteristic 0.
- π̃ is fixed by θ ↦ θ + 1, so every partial quotient of π̃'s continued fraction has even degree.

## Contents

| Path | What |
|---|---|
| `paper/ff-pi-plus-e-q2.tex`, `.md` | The manuscript |
| `scripts/22-write-ff2/verify_q2_paper.py` | An independent check of every identity and degree count in the paper (about 8 s) |
| `scripts/04-funcfield/` | The original study's scripts: windows, degrees, Thakur's continued fraction, rank checks |
| `data/` | Logs and continued-fraction data |

## Reproduce

From this folder (Python 3.12, `python-flint`, `numpy`):

```bash
python scripts/22-write-ff2/verify_q2_paper.py
```

It ends with `=== all checks passed`. This was re-checked from this folder on 23 Oct 2026.

## License

Text: CC BY 4.0. Code: MIT.
