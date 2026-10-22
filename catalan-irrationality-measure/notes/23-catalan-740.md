# 23-catalan-740: a certified μ(G) ≤ 727.03 for Catalan's constant, and where the real slack is

Agent 23-catalan-740, 2026-10-22. CPU budget: 3 cores.
- Scripts: `scripts/23-catalan-740/`. Trial data, certificate logs, controls: `data/23-catalan-740/`.
- Sources read:
  - OpenAI 005 (`data/09-zeta7det/src/openai005_catalan.txt`; §§2, 5–7 re-read);
  - S. Litsey, *An irrationality measure for Catalan's constant*, `paper/catalan-measure.tex`
    (github.com/slitsey/catalan-irrationality-measure, `main` as of 2026-10-22, last Catalan change 2026-10-17), §§2–7, 9.6–11;
  - Litsey's Lean sources, read as text through the GitHub raw API (nothing downloaded into the project, nothing executed):
    `Parametric/{GapMain, GapBudgets, StrongBudgets, GapSecondFields, StrongFirstFields, GapFinalTransfer,
    GapContentTransfer, GapAtomCost, U6Main, FlexibleSecondEnergy, FlexibleFirstEnergy}.lean`,
    `Cert/{FieldEval, ShiftedFields, MarkedResidual}.lean`, `StrongTrial/Data.lean`.

---------------------------------------------------------------------------------------------------

## 0. Results (ranked)

1. **[PROVEN here, conditional on Litsey's kernel-checked Lean development, used unchanged except for one
   trial's data] μ(G) ≤ 727.03.** For every integer a and all sufficiently large b,
   |G − a/b| > b^(−727.03) (exponent exactly 18175647944007/25000000000 = 727.0259…; any μ > 727.024 works
   with the same data). Previous record: 740 (Litsey, Lean).
   - The only new mathematical input is a better **second (κ = 2) trial** for the energy dual of Litsey's
     record design (20,3,1,2,1): 400 + 1000 finite modes, Litsey's 45 geometric tail bases plus 11 new bases at the
     column-support edge frequencies. The data are exact rationals in Litsey's Lean format
     (`data/23-catalan-740/trial_cand2.json`).
   - Its certified dual value is **2.2564477**, against Litsey's 2.2563367: +1.11·10⁻⁴ per n², i.e. +0.044 per N² in η.
     OpenAI's real-place method cannot exceed 2.2564510 for this design (the primal value of an explicit 720×450-cell
     measure). So the trial is within 3.3·10⁻⁶ of the method's ceiling. The ledger floor of this route is μ ≈ 726.66.
   - All field suprema (row, column, marked residual), both norms and p₁ are certified in arb/acb ball
     arithmetic (§3). Every rational inequality of Litsey's generic ledger is replayed exactly, and re-checked
     independently with 50-digit mpmath (§4).
   - An intermediate trial (300 + 600 modes, Litsey's bases only; `trial_cand1.json`) gives 729.31. Its X and M
     certificates pass; its Y run was superseded.
   - Side result [PROVEN here, same status, no new trial]: re-balancing Litsey's own certified constants
     (lowering the target by 2·10⁻⁶ to buy a smaller marked cost) already gives **739.72 < 740** (§4.3).

2. **Status check (2026-10-22, morning; re-checked 10:58 CEST): nothing below 740 is public.** Litsey's repo (last commits are about
   L(2,χ₋₃) = 14.35; the Catalan statement is still 740), Carreira (2294, 21 Oct), Bojnal / proof-factory (2402),
   Jha (10⁸); Zenodo, arXiv (only Sun 2609.04176 and Wachs 2609.22339) and the Palomar registry (`recent.json`)
   contain no Catalan measure below 740.

3. **[NUMERICAL] Route A: the far-sheet estimate really is where the decay is lost, but the obvious fix fails** (§6).
   - At the extremal configuration of OpenAI's κ = 2 functional the true |R_n| is smaller than the
     interpolation bound V(1/x)∏|x|^g by e^{−(0.0105 ± 0.0002) n²} (arb determinants up to n = 200). This is
     larger than 19's 0.0086 excess of the true decay over the method's optimum, so the sheet bound alone can
     account for the whole excess.
   - About half of it (0.0047) is present already at **g = 0**. So the reflected rows are not the only mechanism:
     negative nodes act through w = (x+1/x)/2. For all-negative configurations the bound is off by exactly
     ∏(1−x_ix_j) [PROVEN here]; for all-positive configurations it is sharp.
   - The natural rigorous attempt (Laplace expansion along the g difference rows) makes the functional *worse*
     (optimum −2.2021 vs −2.2566): at the maximiser of the Laplace functional its largest term exceeds the full
     determinant by ≥ e^{0.054n²−o(n²)}, so the gain must come from cancellation between different S.
     A Route A theorem needs a genuine two-sheet (vector) equilibrium bound with cancellation. Not achieved.

4. **Route B (ρ):** Litsey's marked cost c₂ is the closure value of his own budget to within 4.5·10⁻⁴; no slack
   in ρ is available without a new perturbation argument (§7). The only ρ-gain used is the re-balancing in item 1.

5. **[NUMERICAL] Ceiling of this route:** the method's optimum satisfies U* ≤ 2.2564510, so the same ledger cannot
   go below μ ≈ 726.66 with any trial (K, R as in candidate 2). We are within 0.37 of that. Further progress needs Route A or a new perturbation argument.


---------------------------------------------------------------------------------------------------

## 1. How 740 is computed (re-derived; `ledger.py` reproduces it exactly)

Litsey's final criterion (Lean `measure_of_content_on`, `U6Main.measure_740_of_gap_bounds`): for the integer
polynomials Q_N = (content-reduced) (18N)!·D_N·Δ_N(X) of degree ≤ 18N, with |Q_p(x)| ≤ e^{−(η+χ)p²} on
|x − G| ≤ e^{−ρp} and square-free-prime nonvanishing, one gets μ ≤ ρκ whenever

  2 < κ, 0 < δ, 1 < 2(1 − δ(2κ−1)), η + χ < 36, 18 + 72δκ < (η+χ)κ,

so μ → 18ρ/(η+χ). The constants:

| symbol | meaning | Litsey (740) |
|---|---|---|
| U | κ = 2 dual target (005 normalization; physical decay U − ½log 2) | 2.2563367 |
| η | any rational < 400·(U − (121/200)log 2 − 8795/4800) (analytic gap, per N²) | 375279/200000 = 1.876395 |
| χ | coefficient contents (tent 507/2000, U3 3331/10000, U6 624/10000) | 0.649 |
| c₁, c₂ | marked-row costs of the first (κ=1) and second (κ=2) populations | 4.21, 13421/3125 = 4.29472 |
| w | joint quadrature-weight cost (circle r = 1.001) | 17.913 |
| ρ | 20·max(c₁,c₂) + w | 103.8074 |
| κ, δ | PNT window | 3700000/519037, 10⁻⁷ |
| μ | ρκ | **740** (limit 18ρ/(η+χ) = 739.897) |

The second budget (Lean `FlexibleSecondEnergy.SecondTrialBoundsOn.budget`) is, for all a, t ≥ 0,

  (−1/5 + 8a/5 − a² + t − t²)log 2 + K + (1−a)F + a(1−a)log 2 + (a−t)R − 3t/4 + 9a/10 ≤ −U + ½log 2 + (c₂ + 9/10)a,

with K = normCap ≥ 2‖p‖² + ½‖v‖², F = MX + MY (row + column caps), R = marked-residual cap. Writing
h₀ = K + F − (7/10)log 2 + U, α = (13/5)log 2 − F + R − c₂, β = log 2 − R − 3/4, it is equivalent to
h₀ + α₊²/(8 log 2) + β₊²/(4 log 2) ≤ 0. Since −h₀ ≥ 0 is the dual slack, c₂ ≈ (19/10)log 2 + K + U + R.
The first budget (FirstTrialBoundsAt.budget) is the same with t = 0, K₁ + F₁ + a·λ/5, and must hold for all real a.

---------------------------------------------------------------------------------------------------

## 2. The κ = 2 dual and the new trial

**Fields** (Litsey's shifted fields, `Cert/ShiftedFields.lean` with R = ⟨2/5, 1/10, 27/20, 0, 1/5⟩, C = ⟨3/20, 1/20, 0⟩):

  X(x) = (2/5)log|x| + (1/10)log(1−x) − (27/20)log(1+x²) − 4T(p,x) − S(v,x),   x ∈ [−1,0) ∪ (0,1),
  Y(s) = (3/20)log s + (1/20)log(1−s) + 2T(v,s),                              s ∈ (0,1),
  M(x) = −(9/10)log(1+x²) − 2T(p,x),                                           401/2500 ≤ |x| ≤ 1,

T(u,x) = Σ u_k T_k(x)/k, S(u,x) = Σ u_k x^k/k, ‖u‖² = Σ u_k²/k, and trial sequences
u_k = cs_{k−1}/10⁸ + Re Σ_j r_j z_j^k (cs ∈ ℤ, z_j, r_j ∈ ℚ + ℚi, |z_j| < 1). The dual value is
D = (7/10)log 2 − 2‖p‖² − ½‖v‖² − sup X − sup Y (a valid target U is any U ≤ D, up to the budget).

**Validation of my evaluator on Litsey's trial** (`trial.py`, float): 2‖p‖² + ½‖v‖² = 0.7742497408 (Litsey cap
0.77424975), sup X = −0.9775804996 (cap −0.97758049), sup Y = −1.5678029820 (cap −1.56780297),
sup M = −0.0524014 (cap −0.0524), p₁ = 0.18370, D = 2.2563368 (target 2.2563367).

**Why Litsey's trial is stuck.** Re-optimising his own family (160 + 160 modes, his tails) gains only 1.2·10⁻⁷:
it is optimal for its family (78 row and 44 column maxima within 10⁻⁶ of the top). Comparing with the moments of
the primal equilibrium (`primal_moments.py`, 480×300 and 720×450 cells), the deficit sits in v_k for
161 ≤ k ≤ 1000 (Σ(δv_k)²/k = 7.7·10⁻⁵), carried by oscillations at the column-support edges
(frequencies arccos 0.0127 ≈ 1.558 and arccos 0.9999 ≈ 0.014), which his tails do not contain.

**Optimiser** (`optdual4.py`): an exchange method on the full fine grids, using all row and column points within
0.03 of the top (47 000–68 000 and 39 000–58 000 points). Each subproblem is solved exactly by an active-set method
for the simplex-constrained dual QP (support 580–880 points), followed by a convex line search.
- Candidate 1: 300 p-modes and 600 v-modes with Litsey's bases. One step from Litsey's trial gives
  D = 2.2564274502.
- Candidate 2: 400 + 1000 modes. Bases are Litsey's plus z = ρe^{±iω} for ω ∈ {1.558, 1.5633} and ρ ∈ {0.99, 0.995, 0.998};
  ω = 0.0141 with ρ ∈ {0.99, 0.997, 0.999}; and z = −0.99999, −0.999997 (6-decimal rationals, |z|² ≤ 0.99999).
  One step from candidate 1 gives D = 2.2564478374.

Rounding to Litsey's exact format (cs to 10⁻⁸, r to 10⁻¹⁰) costs 2–3·10⁻⁸. Measured on 800 001-point grids:
candidate 1 D = 2.2564274181, candidate 2 D = 2.2564478163.

---------------------------------------------------------------------------------------------------

## 3. Certificates [rigorous: python-flint arb/acb, 256-bit]

Method (`certify3.py`/`cert_new.py`; final runs `certify4.py`/`cert_new2.py`). Work in the angle: x = cos θ, s = cos φ. Chebyshev parts are
evaluated as Re/Im of powers of e^{iθ} at exact midpoints (ball radii grow linearly in k; the three-term
recurrence in ball arithmetic wraps exponentially and was abandoned after a test: radius 10⁸ at x = −0.999).
On a piece [a,b] accept if either
- (MVT) F(m) + h·sup_{[a,b]}|dF/dθ| ≤ cap, with the Chebyshev part of dF enclosed by its midpoint value plus
  h·Σk|c_k| (and h·Σ|c_k| for the value), tails and power series by ball evaluation; or
- (SING) Σ_i c_i log(sup dist_i) + G([a,b]) ≤ cap, where F = Σ c_i log dist_i + G with c_i > 0
  (dist = |x|, 1−x, s, 1−s computed monotonically as |cos θ|, 2 sin²(θ/2) …), valid on every piece;
otherwise bisect. Pieces cover [0, π] (resp. [0, π/2] and slightly beyond), so every domain point, including
x = ±1, lies in an accepted piece.

**Controls.**
- Litsey's own trial and caps (`control_litsey_pass.log`): all pass; worst certified bounds
  X ≤ −0.97758049001, Y ≤ −1.56780297003, M ≤ −0.05240000511 (caps −0.97758049, −1.56780297, −0.0524),
  2‖p‖²+½‖v‖² ∈ 0.7742497407695874… ± 5·10⁻⁷², 200 s.
- Lowered caps (X, Y by 10⁻⁸; M by 10⁻⁵) (`control_litsey_fail.log`): X fails at θ = 0.12370 (x = 0.99236,
  the top row maximum) and M fails at x = −1, exactly where the true suprema exceed the lowered caps. Y still
  passes, correctly: Litsey's Y cap has 1.2·10⁻⁸ slack (true sup −1.5678029820), so it remains above the sup.
- New trial (candidate 1), X cap 5·10⁻⁹ *below* its float supremum (`control_cand1_X_fail.log`): fails at θ = 3.0530
  (x = −0.99609, next to a row maximum of value −0.9776000167 > cap), as it must; the positive side passes.
- arb and float evaluators agree to 2.5·10⁻¹⁴ at 75 random points of the new trial.

**Precision pitfall (fixed).** Repeated complex multiplication of e^{iθ} in ball arithmetic wraps like
(|cos θ| + |sin θ|)^k: at K = 600–1000 the point enclosures of Y near φ = π/4 had radii 10⁹–10⁶⁸, and a first
candidate-2 run failed spuriously at φ = 1.360. Re-anchoring the powers every 16 steps (fresh cos/sin of kθ)
restores radii of 10⁻⁷³. Wrapping only enlarges enclosures, so the earlier passes (candidate 1 X, M; Litsey controls)
remain valid. For speed the final runs add a second-order test,
F(m) + |F′(m)|h + ½h²·sup|F″| ≤ cap. The Chebyshev part of F″ is enclosed by its midpoint value plus h·Σk²|c_k|
(`certify4.py`, driver `cert_new2.py`).

**Candidate 2 (`trial_cand2.json`), the record.**

| quantity | cap used | certified | pieces / time |
|---|---|---|---|
| 2‖p‖² + ½‖v‖² | 774561586/10⁹ | 0.7745615854269015 ± 7·10⁻⁷² | closed form |
| p₁ | ≥ 9/50 | 0.1837172258683555 (exact rational) | — |
| sup X | −97760061/10⁸ | ≤ −0.9776006100775 | 28101 / 253.7 s |
| sup Y | −156820567/10⁸ | ≤ −1.5682056703675 | 32186 / 383 s (two halves) |
| sup M (|x| ≥ 401/2500) | −521786/10⁷ | ≤ −0.0521786114802 | 296 / 0.8 s |
| tails in disk | max|z|² = 0.999994 (z = −0.999997) | exact | — |

**Candidate 2 negative controls** (`control_cand2_{X,Y}_fail.log`): caps lowered to 1.2·10⁻⁸ below the float suprema
(X −0.97760067, Y −1.56820573). Both runs fail, and only where the field genuinely exceeds the cap:
- X fails at x = 0.10819 (local maximum −0.9776006586) and at x = −0.99999999 (local maximum −0.9776006691);
- Y fails at s = 0.13795 (local maximum −1.5682057294).

At each failure point the field value equals the lowered cap to 10⁻¹⁰.

**Candidate 1 (`trial_cand1.json`).** Certified: 2‖p‖² + ½‖v‖² = 0.774496058462024 ≤ 774496059/10⁹;
p₁ = 0.183713160953469; sup X ≤ −0.9775999600 (cap −0.97759996, 115 966 pieces); sup M ≤ −0.0523857242.
Y was not completed (superseded by candidate 2).

---------------------------------------------------------------------------------------------------

## 4. The new ledger (exact; `final_ledger.py`, `ledger.py`)

log 2 ∈ [0.6931471803, 0.6931471808] (Litsey's `Real.log_two_gt_d9 / lt_d9`).
Candidate 2: K = 774561586/10⁹, F = MX + MY = −254580628/10⁸, R = −521786/10⁷. Dual slack bound:
(7/10)log 2 − K − F ≥ 2.2564477202.

| item | value | check |
|---|---|---|
| target U | 56411143/25000000 = 2.25644572 | — |
| c₂ | 13414/3125 = 4.29248 | second budget: h₀ + α₊²/(8 log 2) + β₊²/(4 log 2) ≤ −5.2·10⁻¹¹ with log 2 bounds (−2.3·10⁻¹⁰ at exact log 2) |
| first budget at U (Litsey's first trial, unchanged) | c₁ = 4.21, λ = 1.68493357 | max_a ≤ −2.14·10⁻⁴ < 0, α₁ = 0.408 ≥ 0, h₀ = −0.0303 ≤ 0 |
| c₁ ≤ c₂ | 4.21 ≤ 4.29248 | ✓ |
| ρ = 20c₂ + 17913/1000 | 103.7626 | — |
| η (rational < 400·(U − (121/200)log 2 − 8795/4800) = 1.92000358) | 1.9200034 | ✓ |
| η + χ | 2.5690034 | < 36 ✓ |
| κ, δ | 35033139/5000000, 10⁻⁷ | 18 + 72δκ < (η+χ)κ by 1.9·10⁻⁷; window 1 < 2(1−δ(2κ−1)) ✓ |
| **μ = ρκ** | **18175647944007/25000000000 = 727.0259** | limit 18ρ/(η+χ) = 727.0239 |

Candidate 1 through the same procedure: U = 2.2564252, c₂ = 4.2921, ρ = 103.755, κ = 70290859/10⁷, μ = 729.3028
(`ledger_cand1.log`).

**4.3 Control/side result.** Litsey's own caps through the same optimiser: U = 2.25633463, c₂ = 17167/4000 =
4.29175, ρ = 103.748, κ = 17824887/2500000, **μ = 739.7186** (all checks pass). His published 740 corresponds to
U = 2.2563367, c₂ = 4.29472 (also reproduced exactly).

---------------------------------------------------------------------------------------------------

## 5. Exactly what is used unchanged, and the proof

**Theorem.** μ(G) ≤ 18175647944007/25000000000 < 727.03.

*Proof.* Litsey's Lean development proves, for arbitrary real parameters, the following implications (statements
read in the files named):
1. `weighted_atom_point_on` (GapAtomCost): from `F : FirstTrialBoundsAt U c₁`, `D : SecondTrialBoundsOn U R c₀ c₂ G`,
   `20c₁ ≤ r`, `20c₂ ≤ r` and the weight bound `GapWeightBounds` (with w = 17.913, trial-independent) it gives
   `WeightedAtomPointOn τ G U (r + w)` — the perturbation estimate with radius ρ = r + w.
2. `measure_of_content_on` (GapContentTransfer): from that, the content bound
   `u6_combined_content` (χ = 0.649; arithmetic, target-independent), `η < 400·measureGapAt U` and the window
   inequalities, it gives |G − a/b| > b^{−ρκ} eventually.
3. The structure fields of `SecondTrialBoundsOn` are produced from trial data exactly as in
   `GapSecondFields.gapSecondTrialOfFields`, whose only data-dependent inputs are: tails in the unit disk;
   the shifted row field ≤ MX on [−1,0)∪(0,1); the shifted column field ≤ MY on (0,1); the marked residual ≤ R on
   {x ∈ [−1,1] : |x| ≥ 401/2500}; 2‖p‖² + ½‖v‖² ≤ K; p₁ ≥ 9/50; and the budget inequality. The passage from these
   N-independent fields to the N-dependent physical fields (`shiftedRowField_corrected_eventually`,
   `twenty_shiftedRowField`, `twenty_grouped_marked_field_identity`, `near_residual_field_le_three_fourths`) is generic.
4. `FirstTrialBoundsAt U 4.21` for Litsey's first trial at a new target needs only its budget at U
   (`StrongFirstFields.strongFirstTrialOfFields` uses U only there).

We supply: the second trial `trial_cand2.json` with the certificates of §3, the budgets and payments of §4, and
Litsey's first trial, weights, nodes and contents unchanged. ∎

**Used unchanged:** OpenAI 005 (construction, rationality, nonvanishing at prime scales, real-place machinery
Props. 5.1–5.3, 6.1 as extended by Litsey); all of Litsey's arithmetic (denominators 8795/4800 + (21/200)log 2,
square-denominator nonvanishing, coefficient recycling, contents 507/2000, 3331/10000, 624/10000), his quadrature
nodes and weight bound 17.913, his marked-atom/mixed-Andréief/grouped-population analysis, his first trial,
and his PNT criterion. **Changed:** the second trial and the numerical constants U, K, F, R, c₂, ρ, η, κ.


---------------------------------------------------------------------------------------------------

## 6. Route A: the sheet estimate — what is lost, and why the obvious fix fails

Notation (OpenAI (65), Prop. 5.2): R_n(x) = det[ξ_f(x_i)x_i^{g−r} + ξ_n(x_i)x_i^{r−g}], (ξ_f, ξ_n) = (½,½) for
x < 0 and (−¼, 5/4) for x > 0. OpenAI's κ = 2 bound is |R_n| ≤ (1 + K₀n)^n·far(x), far(x) = V(1/x)∏|x_i|^g.

**6.1 Exact algebra [PROVEN here, elementary].**
(a) For 1 ≤ e ≤ g: row(g−e) − row(g+e) = (ξ_f − ξ_n)(x^e − x^{−e}) = (3/2)·1_{x>0}(x^{−e} − x^e), and
row(g−e) + row(g+e) = x^e + x^{−e} = 2T_e(w), w = (x + 1/x)/2 = 1/t. Hence R_n = 0 unless at least g nodes are positive.
(b) A negative node acts by P ↦ ½(P(x) + P(1/x)), which depends on x only through w. For g = 0 and all nodes negative,
R_n = det[T_r(w_i)] = 2^{C(n−1,2)}V(w), and w_j − w_i = ½(1/x_i − 1/x_j)(x_ix_j − 1), so
  |R_n(x)| = 2^{1−n}·V(1/x)·∏_{i<j}(1 − x_ix_j)   exactly.
The far alternant overestimates by ∏(1 − x_ix_j)^{−1} = exp(−½n²∬log(1−xy)dμ² + o(n²)), a quadratic-scale loss.
(c) (g = 0, mixed signs.) Sym: z^k ↦ T_k(w) is a bijection Poly_{<n}(z) → Poly_{<n}(w); with Π_−(w) = ∏_{j∈N_−}(w − w_j),
  R_n = const·V(w_{N_−})·det[ℓ_i(Q_m)]_{i∈N_+, m<n_+},  Q_m = Sym^{−1}(Π_−·T_m),
so the negative block is exactly a w-Vandermonde and all difficulty sits in the positive functionals
ℓ_i(Q) = (5/4)Q(x_i) − (1/4)Q(1/x_i) on the n_+-dimensional space Sym^{−1}(Π_−·Poly_{<n_+}).

**6.2 Size of the loss [NUMERICAL; arb determinants, `sheet_probe*.py`].** d_n = (log|R_n| − log far)/n²:

| configuration | n = 60 | 100 | 140 | 200 | fit d_∞ |
|---|---|---|---|---|---|
| κ=2 equilibrium quantiles, g = n/10 | −0.01979 | −0.01610 | −0.01453 | −0.01332 | **−0.0105** (a + b/n; −0.01046 with log n/n) |
| same nodes, g = 0 | −0.01731 | −0.01206 | −0.01024 | — | −0.0047 |
| all negative, uniform(−0.9,−0.1), g = 0 | — | (n=120) −0.1643 | — | — | → ½∬log(1−xy) = −0.1585 (6.1(b)) ✓ |
| all positive, uniform(0.1,0.9) | — | (n=120) −0.0070 | — | — | → 0 (far bound sharp) |
| half/half (same intervals) | g=0: −0.0303, −0.0213, −0.0182 (n = 40, 80, 120) | | | | −0.012; with g = n/10: −0.026 |

g-dependence at n = 100, fixed equilibrium nodes: d = −0.0121, −0.0136, −0.0161, −0.0201, −0.0270 for
g/n = 0, 0.05, 0.10, 0.15, 0.20 (≈ −0.37(g/n)² on top of the g = 0 part).

Reading: (i) at the maximiser of OpenAI's κ = 2 functional the far alternant is too large by e^{0.0105n²}. This
exceeds the 0.0086 excess of the true decay over the method's optimum (19), so the sheet bound alone can explain the
whole gap 1.90988 → 1.9185. (ii) Half of it is present at g = 0. It comes from negative nodes (w-structure,
6.1(b)), not only from reflected rows. (iii) Positive nodes alone are sharp, so the gain is a negative/negative
and negative/positive interaction effect.

**6.3 The Laplace bound is worse [NUMERICAL, `laplace_explore.py`].** Expanding along the g difference rows,
|R_n| ≤ e^{O(n log n)}·max_S |det_diff(S)|·|det_rem(S^c)| (S ⊂ positive nodes, |S| = g; det_diff is an exact
product, ∏_S(1−x²)/x·∏_{i<j∈S}|x_i−x_j|(1−x_ix_j)/(x_ix_j); det_rem by 005's Prop. 5.2 with g' = 0).
With σ the S-nodes and τ the rest, the resulting functional is
  Ψ = Φ₂ − ∬log|x−y|dσdτ + ½∬log(1−xy)dσ² + (1 − 2g/n)∫log x dσ,
i.e. the σ–τ repulsion is halved. Optima (piecewise-constant densities, 120×80 cells): Φ₂ (OpenAI) −2.25658;
Ψ −2.20212; ½Φ₂ + ½Ψ −2.23035. All are worse: at the maximiser of Ψ the largest Laplace term exceeds the far bound, hence |R_n|, by
≥ e^{0.054n²−o(n²)}. The gain lives entirely in the cancellation between different S.

**6.4 What a Route A theorem needs.** An upper bound for determinants of two-point functionals
ξ_nδ_x + ξ_fδ_{1/x} that is sharp at the n² scale. Heuristically this is a vector equilibrium on the two sheets with
"balanced" regions (|P(x)| ≈ |P(1/x)|) for negative nodes. Two natural guesses are refuted numerically:
- far × ∏_{neg pairs}(1−xx'): too optimistic (residual +0.024 n² at equilibrium);
- far × ∏_{pairs with a negative node}(1−xx'): worse than the far bound.

The marked atoms of the perturbation step are P-columns, i.e. symmetric functionals, so they are "negative-like".
A Route A theorem would therefore have to be proved for configurations containing them too. Payoff: 19's estimate
μ ≈ 310 with Litsey's ρ. **This remains the most valuable open step.**

---------------------------------------------------------------------------------------------------

## 7. Route B (the perturbation rate ρ)

- ρ = 20c₂ + w. Here w = 17.913 is the weight price; 19 showed contour changes gain ≤ 0.02. c₂ is fixed by the
  closure of the grouped second-population budget, c₂ ≈ (13/5)log 2 − F + R − √(8 log 2·(−h₀)).
  Litsey's 4.29472 is within 4.5·10⁻⁴ of the h₀ = 0 value; the optimal trade of h₀ against c₂ is worth 0.28 in μ
  (included in 739.72, 729.30 and 727.03).
- The true root-separation rate ≈ 90.5 (19) is not a property of the column-expansion bound, which charges every
  m-subset its worst atoms and worst weights. Reaching it needs cancellation between atom choices. I found no
  rigorous handle on that. Not attempted further.

---------------------------------------------------------------------------------------------------

## 8. Weak points

1. **Conditional on Litsey's Lean development.** The new trial is not Lean-checked. The argument uses only Lean
   statements that I read and that are generic in the trial data and constants (§5). I did not compile or run Lean.
   If the generic statements carry hidden data-specific hypotheses that I missed, the claim fails. I checked
   `SecondTrialBoundsOn`, `FirstTrialBoundsAt`, `weighted_atom_point_on`, `measure_of_content_on` and
   `gapSecondTrialOfFields` and found none.
2. **New certification code** (`certify3.py`). It is validated by three controls:
   - it passes Litsey's Lean-certified caps;
   - it fails where caps are lowered below the true suprema;
   - it agrees with an independent float evaluator to 2.5·10⁻¹⁴.

   Rounding is outward (ball arithmetic, upper endpoints). Piece endpoints are exact binary rationals. The cover
   argument is the MVT/SING dichotomy of §3.
3. **Margins are small but positive** (candidate 2):
   - second budget: 5.2·10⁻¹¹ with the log 2 interval (2.3·10⁻¹⁰ at exact log 2);
   - first budget: 2.1·10⁻⁴;
   - payment 18 + 72δκ < (η+χ)κ: 1.9·10⁻⁷;
   - caps versus certified suprema: X 7.7·10⁻¹¹, Y 3.7·10⁻¹⁰, M 1.1·10⁻⁸.

   All of these are exact rational or ball-arithmetic inequalities, not floating-point comparisons.
4. The "method ceiling" U* ≤ 2.2564510 is the value of OpenAI's κ = 2 functional at an explicit piecewise-constant
   measure (720×450 cells, float evaluation of exact cell integrals); it is [NUMERICAL] and not used in the proof.
5. Route A statements in §6.2–6.4 are numerical (finite n ≤ 200, fits with 1/n corrections); only 6.1 is proved.

---------------------------------------------------------------------------------------------------

## 9. Files

**Scripts** (`scripts/23-catalan-740/`):
- `trial.py`: float evaluator of trial sequences, fields, norms;
- `optdual4.py`, `run_ex4.py`, `run_ex4_rich.py`: exchange optimiser with exact active-set QP;
- `optdual.py`, `optdual2.py`, `optdual3.py`, `optsmooth.py`: earlier optimisers (superseded);
- `export_trial.py`: rounding to Litsey's exact format;
- `certify3.py`, `cert_new.py`, `certify4.py`, `cert_new2.py`: rigorous arb certificates (first- and second-order);
- `cert_control3.py`: controls on Litsey's trial;
- `ledger.py`, `final_ledger.py`: exact replay of Litsey's inequalities (`python ledger.py` reproduces 740);
- `energy.py`, `laplace_explore.py`, `primal_moments.py`: primal equilibria, Laplace functional, moments;
- `sheet_probe.py` … `sheet_probe5.py`: §6.2;
- `certify.py`, `certify2.py`: superseded certifier versions (wrapping problems), kept for the record.

**Data** (`data/23-catalan-740/`):
- `trial_cand2.json`: the record trial, exact;
- `trial_cand1.json`: intermediate trial;
- `litsey_strongtrial.json`: Litsey's trial, parsed;
- `cert_cand2_{X,Y_a,Y_b,M}.log` (the record), `cert_cand1_{X,M}.log`; superseded: `cert_cand2_Y_SUPERSEDED_spurious_fail_before_anchoring.log`, `cert_cand1_Y_KILLED_incomplete.log`;
- `control_litsey_{pass,fail}.log`, `control_cand1_X_fail.log`, `control_cand2_{X,Y}_fail.log`;
- `ledger_cand2.log`, `ledger_cand1.log`, `ledger_litsey_rebalanced.log`, `ledger_litsey_740_control.log`;
- `opt_exchange_K300_600.log`;
- `primal_{480,720}.log`.

Rerun: `python scripts/23-catalan-740/cert_new2.py data/23-catalan-740/trial_cand2.json X -0.97760061`
(and `Y -1.56820567`), `cert_new.py … M -0.0521786`, `cert_new.py … NORM`, then
`python scripts/23-catalan-740/final_ledger.py 0.774561586 -0.97760061 -1.56820567 -0.0521786`.
