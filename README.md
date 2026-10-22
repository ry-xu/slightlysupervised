# slightlysupervised

Experiments and results from projects run largely by AI agents, with a human (Ryan Xu) choosing problems, checking output
and deciding what to publish. Each folder is self-contained, and its README opens with the result's status.

| Folder | What | Status |
|---|---|---|
| [`zeta5-irrationality-measure/`](zeta5-irrationality-measure/) | μ(ζ(5)) ≤ 95.91, improving the bound of 260 | AI-generated; independently re-checked by AI; not yet verified by a human expert |
| [`catalan-irrationality-measure/`](catalan-irrationality-measure/) | μ(G) ≤ 727.03 for Catalan's constant (previous: 740) | Builds on S. Litsey's Lean development; new trial certified in ball arithmetic, not yet run through Lean |
| [`function-field-pi-plus-e/`](function-field-pi-plus-e/) | The q = 2 Carlitz analogue of "π + e, πe, π/e are irrational" | AI-generated; elementary proof, checked by a second agent and by scripts; novelty not yet confirmed with specialists |
| [`wilhelm-1583-cipher/`](wilhelm-1583-cipher/) | A 1583 cipher letter to Duke Wilhelm V of Bavaria, read for the first time | Key recovered, 86% of cipher words read directly; not yet checked against the archive images |
| [`mnist_how_many_pixels.ipynb`](mnist_how_many_pixels.ipynb) | How many pixels do you need to solve MNIST? | Notebook |

Also see [Scattered Threads](https://ryanxu.net/scattered-threads/), a search for objects divided between museums.

Corrections welcome: ryan@ryanxu.net. Unless a folder says otherwise, text is CC BY 4.0 and code is MIT (see `LICENSE`).
