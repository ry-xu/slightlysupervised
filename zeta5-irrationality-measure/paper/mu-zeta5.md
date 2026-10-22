# The irrationality exponent of $\zeta(5)$ is at most $95.91$

**Claude Opus 5.5 (Anthropic)**

*Released by Ryan Xu (ryan@ryanxu.net), who directed the project but is not a specialist and has not independently verified the mathematics.*

*Preprint draft, 22 October 2026. A LaTeX version of the same text is in `papers/mu-zeta5.tex`.*

---

## Abstract

A. Fauzan recently proved that $\zeta(5)$ is irrational. His proof uses integer polynomials $Q_{K,M}$ of degree $37K/40$, obtained by normalizing Hankel determinants of a positive moment functional, and a perturbation argument then gives the irrationality exponent bound $\mu(\zeta(5))\le 260$. We prove $\mu(\zeta(5))\le 95.9083454190$; in particular $\mu(\zeta(5))<95.91$.

The main new ingredient is an essentially sharp bound for the norm of the matrix $B_K=G_K(\zeta(5))^{-1/2}\,[X]G_K(X)\,G_K(\zeta(5))^{-1/2}$, which controls the perturbation argument: $\log\lVert B_K\rVert\le 7.469616\,K+O(\log K)$, in place of Fauzan's $13.866\,K$. The proof has three steps. We write the quadratic form of $[X]G_K(X)$ as a contour integral around the poles of the functional. We bound polynomials at complex points by a weighted Bernstein–Walsh inequality with a near-equilibrium test measure, combined with a Nikolskii-type inequality. Finally we move the contour off the segment of poles, where the alternating residues cancel. The exponent is then the value of a harmonic function at a saddle point near $\tau\approx-1.0699$. Numerically $K^{-1}\log\lVert B_K\rVert\to 7.46958$, so the new estimate is sharp to within $4\cdot 10^{-5}$, and about $95.89$ is the limit of this construction and argument.

We also sharpen Fauzan's decay rate to $\varepsilon>0.0720416$, using a $1600$-component test measure in the logarithmic-energy estimate, exact integration of the arithmetic limit function, and a larger fixed cutoff. All numerical constants are certified in ball arithmetic, and the certificates were reproduced by independently written code. The proof uses several statements of Fauzan's paper unchanged. Most of these are formalized in Lean by others, but none of the new arguments of this paper is formally verified.

---

## Contents

1. Introduction
2. Setup
3. The analytic constant
4. The arithmetic constant at the cutoff $M=4\cdot10^8$
5. The norm bound
6. Proof of the main theorem
7. Sharpness
8. Computational certification and reproducibility
9. Remarks and further directions
- Acknowledgements and use of AI
- References

---

## 1. Introduction

Let $\xi$ be an irrational real number. Its *irrationality exponent* $\mu(\xi)$ is the supremum of the real numbers $\mu$ for which

$$
0<\Bigl|\xi-\frac ab\Bigr|<b^{-\mu}
$$

has infinitely many solutions $(a,b)\in\mathbb Z\times\mathbb Z_{>0}$. Equivalently, $\mu(\xi)\le\mu_0$ if and only if for every $\mu>\mu_0$ there is $b_0$ such that $|\xi-a/b|>b^{-\mu}$ for all $a\in\mathbb Z$ and all integers $b\ge b_0$. Dirichlet's theorem gives $\mu(\xi)\ge2$ for every irrational $\xi$, and $\mu(\xi)=2$ for almost all $\xi$. Apéry proved that $\zeta(3)$ is irrational [Ap], and the classical bound for it is $\mu(\zeta(3))\le 5.513891$ [RV].

In September 2026 A. Fauzan [F] proved that $\zeta(5)$ is irrational. He also showed [F, Cor. 1.2] that there is $b_0$ with

$$
\Bigl|\zeta(5)-\frac ab\Bigr|>b^{-260}\qquad(a\in\mathbb Z,\ b\in\mathbb Z,\ b\ge b_0), \tag{1}
$$

so $\mu(\zeta(5))\le260$. The irrationality theorem [F, Thm. 1.1] has been formalized in Lean 4 in two independent projects [DR, Mo]. Neither covers the exponent bound (1), which is proved in Appendix C of [F].

### 1.1 Fauzan's method

Let $K\in40\mathbb Z_{>0}$, $h=37K/40$ and $\lambda=h/K=37/40$, and fix an integer cutoff $M\ge40$. Fauzan constructs an $h\times h$ Hankel matrix $G_K(X)$ whose entries are affine in $X$. He also constructs an explicit positive rational normalization $m_{K,M}S_K$ such that, for $K\ge200M^2$,

$$
Q_{K,M}(X)=m_{K,M}S_K\det G_K(X)\in\mathbb Z[X],\qquad \deg Q_{K,M}=h .
$$

At $X=\xi:=\zeta(5)$ the matrix $G:=G_K(\xi)$ is the Gram matrix of a positive weight, so $Q_{K,M}(\xi)>0$. A logarithmic-energy estimate for $\det G$, together with $p$-adic estimates for the normalization, shows that

$$
\limsup_{K\to\infty}K^{-2}\log Q_{K,M}(\xi)\le A_M+\overline U . \tag{2}
$$

Here $A_M$ is an arithmetic constant and $\overline U$ an analytic one. When $A_M+\overline U<0$, this proves irrationality.

The exponent comes from a perturbation argument. Let $M_K=[X]G_K(X)$ be the coefficient matrix of $X$, and $B_K=G^{-1/2}M_KG^{-1/2}$. If $\delta=a/b-\xi$, then $Q_{K,M}(a/b)=Q_{K,M}(\xi)\det(I+\delta B_K)$. Suppose that $\lVert B_K\rVert\le e^{(\varrho+o(1))K}$ and $Q_{K,M}(\xi)\le e^{-\varepsilon K^2}$. If $|\delta|\le b^{-\mu}$ and $K\approx c\log b$ with $\lambda/\varepsilon<c<\mu/\varrho$, then $b^hQ_{K,M}(a/b)$ is a positive integer smaller than $1$, which is absurd. Hence every $\mu>\varrho\lambda/\varepsilon$ is admissible (Proposition 1). Fauzan proves $\varrho=6933/500=13.866$ [F, Prop. C.1] and takes $\varepsilon=79/1600$ with $M=10^5$ [F, (C.7)]. This gives $\varrho\lambda/\varepsilon=259.77$.

### 1.2 Results

> **Theorem 1.** For every real number $\mu>95.9083454190$ there is $b_0=b_0(\mu)$ such that
> $$\Bigl|\zeta(5)-\frac ab\Bigr|>b^{-\mu}\qquad\text{for all } a\in\mathbb Z \text{ and all integers } b\ge b_0 .$$
>
> Consequently $\mu(\zeta(5))\le 95.9083454190<95.91$.

More precisely, every $\mu>\mu_0:=\varrho_2\lambda/\varepsilon_{\max}$ is admissible. Here $\varrho_2=466851/62500=7.469616$ is the norm rate of Theorem 2, and $\varepsilon_{\max}\ge 0.0720416431939879$ is the decay rate of Corollary 1. In ball arithmetic $\mu_0=95.908345418981323\ldots$ (Section 6). As in [F], $b_0$ is not effective.

**The norm bound.** The main new ingredient replaces [F, Prop. C.1], which asserts $\lVert B_K\rVert\le e^{13.866K}$ for all large $K$.

> **Theorem 2.** For every $K\in40\mathbb Z_{>0}$,
> $$\log\lVert B_K\rVert\le\varrho_2K+4\log K+2\log h+31.65+\frac{559}{K},\qquad \varrho_2=\frac{466851}{62500}=7.469616 .$$

The proof (Section 5) has three steps.

1. By the residue theorem, the quadratic form $q^{\mathsf T}M_Kq$ is a contour integral of $\tau^2\Omega(\tau)p(\tau)^2$. Here $p(\tau)=q(K^2\tau)$, $\Omega$ is an explicit rational function, and the contour surrounds the segment $[-1,-\alpha^2]$ that carries the rescaled poles of the functional (Lemma 1).
2. A weighted Bernstein–Walsh inequality bounds $|p|$ at complex points by $e^{KU^\rho}$ times its weighted size on the support of a test measure $\rho$. A Nikolskii-type inequality then bounds that size by the Gram form $q^{\mathsf T}Gq$ (Lemmas 2 and 3).
3. Moving the contour off the pole segment lets the alternating signs of the residues cancel. The exponent becomes $\max_\Gamma(2U^\rho-\widetilde V)-\inf_{\operatorname{supp}\rho}(2U^\rho-V)$. With the best contour, the maximum sits at a saddle point $\tau_*\approx-1.0699$ of the harmonic function $2U^\rho-\widetilde V$, just beyond the end of the pole segment (Proposition 6, Section 5.6).

If the residues are instead bounded in absolute value, without a contour, the same method gives the rate $7.5299683$ (Proposition 5). That alone already yields $\mu(\zeta(5))\le96.69$.

**Sharper decay.** We also improve both constants in (2).

- A test measure with $1600$ nested arcsine components gives $\overline U\le-1.3854960826$ (Proposition 3). Fauzan's $16$-component measure gives $-1.3669956$, and no admissible measure can do better than $\overline U_{\rm eq}\approx-1.3855098$ (Remark 2).
- For the cutoff $M=4\cdot10^8$ we get $A_M\le 1.3134544394$ (Proposition 4). This uses exact rational integration of Fauzan's arithmetic limit function $R(x)x^{-3}$ over $[20,4000]$ and an explicit tail estimate. For comparison, [F, (5.20)] gives $A_{10^5}=1.3175717$.

Together these give $\varepsilon_{\max}=-(A_M+\overline U)\ge0.0720416$, against Fauzan's $79/1600=0.049375$.

**Sharpness.** We computed $\lVert B_K\rVert$ numerically for $K\le1280$ by three independent methods. Numerically, the rate $K^{-1}\log\lVert B_K\rVert$ approaches $7.46958$. This agrees with the best constant obtainable from Theorem 2's method, $\varrho_2^{\rm eq}=7.469583$ (Proposition 8). Our $\varrho_2$ exceeds it by $3.3\cdot10^{-5}$, and $\overline U$ exceeds its optimum by $1.4\cdot10^{-5}$. Within this construction and this argument, the exponent cannot go below about $95.89$ (Section 7). Table 1 summarizes the constants.

**Table 1.** Norm rate $\varrho$, decay rate $\varepsilon$, and the resulting exponent threshold $\varrho\lambda/\varepsilon$, with $\lambda=37/40$.

| route | $\varrho$ | $\varepsilon$ | $\varrho\lambda/\varepsilon$ |
|---|---|---|---|
| Fauzan [F] | $6933/500=13.866$ | $79/1600=0.049375$ | $259.77$ |
| intermediate version: new $\overline U$ and $A_{10^5}$ only | $13.866$ | $0.0708798$ | $180.96$ |
| Theorem 2 with Fauzan's $\varepsilon$ | $7.469616$ | $0.049375$ | $139.94$ |
| Proposition 5 (no contour), new $\varepsilon$ | $7.5299683$ | $0.0720416$ | $96.68$ |
| **Theorem 2, new $\varepsilon$ (Theorem 1)** | $7.469616$ | $0.0720416$ | **$95.908$** |
| limit of the method (numerical, §7) | $\approx7.46958$ | $\approx0.0720558$ | $\approx95.89$ |

### 1.3 Status of the proof

We have tried to make the logical status of every ingredient explicit.

- Several statements of [F] are used unchanged. They are listed in Section 2.2 with Fauzan's numbering. Most of them are formalized in the Lean development [DR], with Fauzan's own constants and with the prime number theorem as an axiom (Table 2).
- Fauzan's Proposition C.1 is not used. The short perturbation argument of [F, App. C.2] is re-proved here in general form (Proposition 1).
- Fauzan's Proposition 6.3 is applied to a new test measure. We check that its proof uses only properties that our measure has (Proposition 2). This step is a paper-level argument.
- The new analysis (Lemmas 1–5, Propositions 5 and 6, Theorem 2) is proved in full here.
- Every numerical constant is a *certified computation*. Each is an inequality between explicit functions, checked on exact rational data in ball (midpoint–radius interval) arithmetic [Jo, flint] with adaptive subdivision. Each certificate used in the proof was re-implemented independently, with a different subdivision scheme (Section 8). These are rigorous computer-assisted verifications, not formal proofs.
- None of the new material is formally verified.

An intermediate version of this work sharpened only $\overline U$ and $A_{10^5}$ and kept Fauzan's $\varrho$; it gave $\mu(\zeta(5))\le181$. We know of one other claimed exponent bound for $\zeta(5)$ since [F]. P. Anand [An] states an effective bound with exponent $25687$, and the record announces a revised exponent $214$, not posted as of 22 October 2026. We have not examined that work.

### 1.4 Organization

Section 2 fixes notation, lists the statements taken from [F], and proves the transfer principle. Section 3 treats the analytic constant $\overline U$ and Section 4 the arithmetic constant $A_M$. Section 5 proves Theorem 2. Section 6 assembles the proof of Theorem 1. Section 7 discusses sharpness, Section 8 the computations, and Section 9 further directions.

---

## 2. Setup

### 2.1 Notation

Throughout, as in [F, (2.1)],

$$
K\in40\mathbb Z_{>0},\qquad N=\frac{3K}{40},\qquad h=\frac{37K}{40}=K-N,\qquad \alpha=\frac3{40},\qquad \lambda=\frac{37}{40},
$$

and $\xi=\zeta(5)$. Let $D_m(t)=\prod_{j=1}^m(t+j^2)$ and $H^{(5)}_j=\sum_{v=1}^jv^{-5}$. The $\mathbb Q$-linear functional $\mu_X$, with values in $\mathbb Q[X]$, is defined on polynomials and on the fractions $1/(t+j^2)$, $j\ge1$, by [F, (2.2)–(2.3)]

$$
\mu_X(t^e)=\frac{(-1)^eB_{2e+2}(2e+3)(2e+4)(2e+5)}{24},\qquad
\mu_X\Bigl(\frac1{t+j^2}\Bigr)=j^4\bigl(X-H^{(5)}_j\bigr)-\frac14+\frac1{2j}. \tag{3}
$$

It is extended to rational functions with simple poles at negative integer squares by polynomial division and partial fractions. Following [F, (2.4)–(2.6)], put

$$
G_K(X)=\Bigl[\mu_X\Bigl(\frac{D_N(t)^6t^{i+j}}{D_K(t)}\Bigr)\Bigr]_{0\le i,j<h},\qquad \Delta_K=\det G_K,\qquad F_K=S_K\Delta_K,\qquad Q_{K,M}=m_{K,M}F_K .
$$

Here $S_K=(K!)^{2h}4^{h-1}\big/\bigl((N!)^{12h}\prod_{i=1}^{h-1}((2i)!)^2\bigr)$, and $m_{K,M}>0$ is the rational number defined in [F, (5.1)–(5.2)] for integers $M\ge40$ and $K\ge200M^2$. Only the simple-pole values depend on $X$, and they do so affinely. Hence $G_K(X)=G_K(0)+XM_K$ with

$$
M_K:=[X]\,G_K(X).
$$

We write $G:=G_K(\xi)$. Since $G$ is positive definite (see (F2) below), we may put $B_K:=G^{-1/2}M_KG^{-1/2}$; $\lVert\cdot\rVert$ denotes the Euclidean operator norm. For a real polynomial $q(t)=\sum_{i<h}q_it^i$ we use the same letter for its coefficient vector, so that $q^{\mathsf T}G_K(X)q=\mu_X(D_N^6q^2/D_K)$.

### 2.2 Statements used from Fauzan's paper

The following statements of [F] are used without change. Brackets give Fauzan's numbering.

- **(F1)** [Thm. 2.1, first sentence; Prop. 5.1; (2.9)] For each integer $M\ge40$ and each $K\in40\mathbb Z_{>0}$ with $K\ge200M^2$: $Q_{K,M}\in\mathbb Z[X]$, $\deg Q_{K,M}=h$, and $Q_{K,M}(\xi)>0$.
- **(F2)** [Prop. 2.2, (2.10)] For every rational function $R$ in the domain of $\mu_X$,
  $$\mu_\xi(R)=\int_0^\infty R(y^2)\,w(y)\,\mathrm dy,\qquad w(y)=\frac{(2\pi)^4y^5}{12}\sum_{\ell\ge1}\ell^4e^{-2\pi\ell y}>0 .$$
  In particular, for every real polynomial $q\ne0$ of degree $<h$,
  $$q^{\mathsf T}Gq=\int_0^\infty\frac{D_N(y^2)^6q(y^2)^2}{D_K(y^2)}\,w(y)\,\mathrm dy>0 .$$
- **(F3)** [Prop. 5.2, (5.10)] For each fixed integer $M\ge40$,
  $$\limsup_{K\to\infty,\ 40\mid K}K^{-2}\log m_{K,M}\le I_{\rm out}+\frac{6\lambda}M+\int_3^M\frac{R(x)}{x^3}\,\mathrm dx,\qquad I_{\rm out}=\frac{127751}{96000},$$
  where $R$ is the function recalled in (11) below. The proof uses the prime number theorem.
- **(F4)** [(5.12)–(5.14)] For $x\ge3$, $R(x)=xF(x)+Q(x)$ with $F(x)=4\lambda+2\lambda\{x\}-12\lambda\{\alpha x\}$ and $-\frac12\le Q(x)\le\frac{13}8$. Here $\{\cdot\}$ is the fractional part.
- **(F5)** [(5.15), (5.17)] If $40\mid M$, then $\displaystyle\int_M^\infty\frac{R(x)}{x^3}\,\mathrm dx\ge-\frac\lambda M+\frac{2923/240-1/4}{M^2}-\frac{32}{M^3}$.
- **(F6)** [(5.18)] $\displaystyle\int_3^{20}\frac{R(x)}{x^3}\,\mathrm dx=\frac{322437603634266857629}{7535670527041937280000}$.
- **(F7)** [§6.1–6.2] The proof of Proposition 6.3, consisting of Lemma 6.2 and (6.6)–(6.15). We use it for a new test measure; see Proposition 2.
- **(F8)** [(6.12)] For $1\le m\le K$ and $t\ge0$,
  $$0\le\sum_{j=1}^m\log\bigl(t+(j/K)^2\bigr)-K\int_0^{m/K}\log(t+u^2)\,\mathrm du\le2\log K+2 .$$
- **(F9)** [(A.1), (A.2), (A.5), App. A.3] The potential of an arcsine measure, the energy of a nested sum of arcsine measures, the closed form of the field $V$, and the unimodality of $V$. These are recalled in Section 3.1.
- **(F10)** [App. C.2] The perturbation argument. We re-prove it in general form (Proposition 1), and we re-derive the residue formula [(C.3)] (Lemma 1).

*Not used:* the specific measure of [F, Lemma 6.1] (Table 1 of [F]) and its constants (6.2), (6.4), (A.9)–(A.11); the tail bound (5.16) and the constants (5.19)–(5.20), (7.2); Proposition C.1 together with (C.2) and (C.4)–(C.6); and the height bound (C.8).

**Table 2.** Formalization status of the statements of [F] used here, in the Lean development [DR] (commit `fc2fefa`, 25 Sep 2026). The development proves Theorem 1.1 of [F] with no `sorry`, and its only non-standard axiom is the prime number theorem $\theta(x)\sim x$. The second formalization [Mo] proves Theorem 1.1 with a modified normalization and does not follow [F, (5.11)] as printed.

| item of [F] | declaration in [DR] | remarks |
|---|---|---|
| (F1) Thm. 2.1 | `theorem_2_1_integral_degree_pos` | all $M\ge40$, $K\ge200M^2$ |
| (F2) Prop. 2.2 | `prop_2_2_moment`, `prop_2_2` | |
| (F3) Prop. 5.2 | `prop_5_2` | all $M\ge40$; prime number theorem as an axiom |
| (F4) (5.12)–(5.14) | `AppendixB.RR_closed`, `AppendixB.eq_5_14` | for all real $x$ |
| (F5) (5.17) | `AppendixB.eq_5_16_5_17` | combined with (5.16), for $40\mid M$ |
| (F6) (5.18) | `AppendixB` | exact rational |
| (F7) Prop. 6.3 | `prop_6_3` | only for Fauzan's measure ($\overline U=-2733991/2000000$), with a different regularization; **not** for the measure used here |
| (F8) (6.12) | module `Sec6/Gram` | |
| (F9) (A.1), (A.5) | modules `Sec6/ArcsinePot`, `Sec6/Field` | (A.2) is formalized as an inequality |
| (7.1) | `eq_7_1` | only with Fauzan's $\overline U$; re-derived here (Corollary 1) |
| (F10) App. C | — | not formalized in [DR] or [Mo] |

### 2.3 The transfer principle

The following is Fauzan's argument from [F, App. C.2], stated for general rates.

> **Proposition 1 (Transfer principle).** Fix an integer $M\ge40$, and suppose the following hold.
>
> (a) There is $\varepsilon>0$ such that $Q_{K,M}(\xi)\le e^{-\varepsilon K^2}$ for all sufficiently large $K\in40\mathbb Z$.
>
> (b) There are constants $\varrho>0$, $A\ge0$ and $C>0$ such that $\lVert B_K\rVert\le CK^Ae^{\varrho K}$ for all $K\in40\mathbb Z_{>0}$.
>
> Then for every $\mu>\varrho\lambda/\varepsilon$ there is $b_0$ such that $|\xi-a/b|>b^{-\mu}$ for all $a\in\mathbb Z$ and all integers $b\ge b_0$.

*Proof.* Choose $c$ with $\lambda/\varepsilon<c<\mu/\varrho$. Suppose that $a\in\mathbb Z$ and an integer $b\ge2$ satisfy $|\xi-a/b|\le b^{-\mu}$. Put $K=40\lceil c\log b/40\rceil$, so that $c\log b\le K<c\log b+40$, and $\delta=a/b-\xi$.

*Step 1.* We have $|\delta|\,\lVert B_K\rVert\le C(c\log b+40)^Ae^{40\varrho}\,b^{\varrho c-\mu}$. Since $\varrho c<\mu$ and $h\le K<c\log b+40$, it follows that $h|\delta|\,\lVert B_K\rVert\le\frac12$ for all $b$ larger than some $b_1$.

*Step 2.* Since $G_K(X)=G+(X-\xi)M_K$, we have $G_K(a/b)=G^{1/2}(I+\delta B_K)G^{1/2}$. Hence $\Delta_K(a/b)=\Delta_K(\xi)\det(I+\delta B_K)$. Multiplying by $m_{K,M}S_K>0$ gives

$$
Q_{K,M}(a/b)=Q_{K,M}(\xi)\det(I+\delta B_K).
$$

For $b\ge b_1$ the eigenvalues of the symmetric matrix $I+\delta B_K$ lie in $[1-|\delta|\lVert B_K\rVert,\,1+|\delta|\lVert B_K\rVert]\subset(0,\infty)$. Moreover $\det(I+\delta B_K)\le(1+|\delta|\lVert B_K\rVert)^h\le\exp(h|\delta|\lVert B_K\rVert)\le e^{1/2}<2$. Hence $0<Q_{K,M}(a/b)<2\,Q_{K,M}(\xi)$.

*Step 3.* For $b$ large, $K\ge200M^2$ and $K$ is large enough for (a). By (F1), $Q_{K,M}\in\mathbb Z[X]$ has degree $h$, so $b^hQ_{K,M}(a/b)$ is a positive integer. On the other hand $h\log b=\lambda K\log b\le\lambda K^2/c$, so

$$
b^hQ_{K,M}(a/b)<2\exp(h\log b-\varepsilon K^2)\le2\exp\bigl(-(\varepsilon-\lambda/c)K^2\bigr),
$$

which is $<1$ for large $b$, since $\varepsilon-\lambda/c>0$. This is a contradiction once $b$ exceeds a bound $b_0$ that depends only on $\mu,\varepsilon,\varrho,A,C,M$. ∎

Fauzan applies this with $M=10^5$, $\varrho=6933/500$, $\varepsilon=79/1600$ and $c=75/4$.

---

## 3. The analytic constant

### 3.1 Potentials, the field, and arcsine measures

For a finite positive Borel measure $\nu$ on $\mathbb C$ with compact support we use Fauzan's sign conventions:

$$
U^\nu(z)=\int\log|z-u|\,\mathrm d\nu(u),\qquad I(\nu)=\iint\log|t-u|\,\mathrm d\nu(t)\,\mathrm d\nu(u),\qquad I(\nu_1,\nu_2)=\int U^{\nu_1}\,\mathrm d\nu_2 .
$$

The external field [F, (6.1), (A.5)] is, for $t\ge0$,

$$
\begin{aligned}
V(t)&=2\pi\sqrt t+\int_0^1\log(t+u^2)\,\mathrm du-6\int_0^\alpha\log(t+u^2)\,\mathrm du\\
&=\log(1+t)-6\alpha\log(t+\alpha^2)-2+12\alpha+2\sqrt t\Bigl(\pi+\arctan\frac1{\sqrt t}-6\arctan\frac\alpha{\sqrt t}\Bigr).
\end{aligned}\tag{4}
$$

For $t>0$,

$$
V'(t)=\frac{\pi+\arctan(t^{-1/2})-6\arctan(\alpha t^{-1/2})}{\sqrt t},\qquad
V''(t)=-\frac{\pi}{2t^{3/2}}-\eta(1,t)+6\eta(\alpha,t), \tag{5}
$$

where $\eta(a,t)=\int_0^a(t+u^2)^{-2}\,\mathrm du=\dfrac{a}{2t(t+a^2)}+\dfrac{\arctan(a/\sqrt t)}{2t^{3/2}}$.

*Unimodality of $V$* [F, App. A.3]. Let $\Phi(y)=V(y^2)$. Then $\Phi'(y)=2(\pi+\arctan(1/y)-6\arctan(\alpha/y))$ and $\Phi''(y)=-2/(1+y^2)+12\alpha/(\alpha^2+y^2)$. The latter changes sign exactly once, at $y^2=711/880$. Since $\Phi'(0+)=-3\pi<0<2\pi=\Phi'(+\infty)$, $\Phi'$ has exactly one zero. Hence $V$ is decreasing on $[0,t_V]$ and increasing on $[t_V,\infty)$ for a unique $t_V\approx0.0059205$.

For $a<b$ let $\omega[a,b]$ be the arcsine probability measure $\mathrm dt/(\pi\sqrt{(t-a)(b-t)})$ on $[a,b]$. Put $m=(a+b)/2$ and $r=(b-a)/2$. Then [F, (A.1)], [Sa, Ex. 1.11]

$$
U^{\omega[a,b]}(z)=\log\frac r2+\log\Bigl|\zeta+\sqrt{\zeta-1}\,\sqrt{\zeta+1}\Bigr|,\qquad \zeta=\frac{z-m}r, \tag{6}
$$

with principal square roots. This function is continuous on $\mathbb C$ and harmonic on $\mathbb C\setminus[a,b]$. It equals $\log(r/2)=\log((b-a)/4)$ on $[a,b]$, and $U^{\omega[a,b]}(z)-\log|z|\to0$ as $z\to\infty$. On the real line off $[a,b]$ it equals $\log(r/2)+\operatorname{arccosh}(|t-m|/r)$, with derivative $\operatorname{sgn}(t-m)/\sqrt{(t-a)(t-b)}$. Since $\operatorname{arccosh}$ is increasing and concave on $[1,\infty)$, $U^{\omega[a,b]}$ is nonincreasing and concave on $(-\infty,a]$, constant on $[a,b]$, and nondecreasing and concave on $[b,\infty)$.

**Test measures.** We consider $\rho=\sum_{j=1}^Jc_j\,\omega[a_j,b_j]$ satisfying:

- **(P1)** $c_j>0$ for all $j$, and $\sum_jc_j=\lambda$;
- **(P2)** $0<a_J<\dots<a_1<b_1<\dots<b_J<2$, i.e. the intervals are strictly nested with index $1$ innermost;
- **(P3)** $\displaystyle\sum_{j=1}^J\frac{8c_j\sqrt2}{\pi\sqrt{b_j-a_j}}\le60$.

Then $S:=\operatorname{supp}\rho=[a_J,b_J]$, and nesting gives [F, (A.2)]

$$
I(\rho)=\sum_{j=1}^J\bigl(S_j^2-S_{j-1}^2\bigr)\log\frac{b_j-a_j}4,\qquad S_j=c_1+\dots+c_j,\ S_0=0 . \tag{7}
$$

Indeed, with $L_j=b_j-a_j$ and $\omega_j=\omega[a_j,b_j]$, $I(\omega_i,\omega_j)=\log(L_{\max(i,j)}/4)$, because $U^{\omega_j}$ is constant on $[a_j,b_j]\supset[a_i,b_i]$ for $i\le j$. By the monotonicity properties above, $U^\rho$ is nonincreasing on $(-\infty,a_1]$, constant on $[a_1,b_1]$, and nondecreasing on $[b_1,\infty)$. In particular $\max_{[l,r]}U^\rho=\max(U^\rho(l),U^\rho(r))$ for every interval. Moreover $U^\rho$ is concave on every interval that contains no $a_j$, $b_j$ in its interior. Finally let [F, (6.3)]

$$
C_*=-2\lambda+12\alpha\lambda(1-\log\alpha)+3\lambda^2-2\lambda^2\log(2\lambda)=2.6530359903404886611\ldots
$$

### 3.2 Fauzan's Proposition 6.3 for general test measures

> **Proposition 2.** Let $\rho$ satisfy (P1)–(P3), and let $M_0\in\mathbb R$ satisfy
> $$2U^\rho(t)-V(t)\le M_0\quad(0\le t\le2)\qquad\text{and}\qquad \tfrac{13}{10}\log2+3\alpha^3-2\pi\sqrt2<M_0 . \tag{8}$$
>
> Then for every $K\in40\mathbb Z_{>0}$ we have $F_K(\xi)>0$ and
> $$\log F_K(\xi)\le\overline U K^2+24K\log K+200K,\qquad \overline U:=\lambda M_0-I(\rho)+C_* .$$

*Proof.* We follow [F, §6.1–6.2] and indicate where $\rho$ enters. Nothing else in that proof refers to $\rho$.

*(a) The bound [F, (6.2)] on $[0,\infty)$.* For $t\ge2$, (P2) gives $0<t-u<t$ on $\operatorname{supp}\rho$, so $U^\rho(t)\le\lambda\log t$. Since $\int_0^1\log(t+u^2)\,\mathrm du\ge\log t$ and $\log(t+u^2)\le\log t+u^2/t$, this gives [F, (6.8)]:

$$
2U^\rho(t)-V(t)\le\tfrac{13}{10}\log t+\frac{6\alpha^3}t-2\pi\sqrt t\qquad(t\ge2),
$$

using $2\lambda-1+6\alpha=13/10$. The right side, also after adding $\sqrt t/K$ for $K\ge2$, is decreasing on $[2,\infty)$. By (8) its value at $t=2$ is below $M_0$. Together with (8) on $[0,2]$, this yields [F, (6.2)] and [F, (6.7)] with our $\rho$ and $M_0$.

*(b) Regularization.* Let $t_1,\dots,t_h\in[0,\infty)$. Let $\omega_i$ be normalized arclength on the circle of radius $\epsilon$ about $t_i$, and $\sigma=K^{-1}\sum_i\omega_i$. By (P1), $\sigma-\rho$ has total mass $h/K-\lambda=0$. The potentials of circle and arcsine measures are continuous, so [F, Lemma 6.2] applies to $\sigma-\rho$ exactly as in [F, p. 18]. The only quantitative property of $\rho$ used in the regularization is the bound [F, p. 19]

$$
\int U^\rho\,\mathrm d\omega_i-U^\rho(t_i)\le\sum_j\frac{8c_j\sqrt{2\epsilon}}{\pi\sqrt{b_j-a_j}} .
$$

It holds for every positive combination of arcsine measures, by the mass bound $4\sqrt{2r}/(\pi\sqrt L)$ for an arcsine measure on an interval of length $L$, and by (P3) it is $\le60\sqrt\epsilon$. With (a), this gives [F, (6.6)] verbatim, and then [F, (6.9)] with $\epsilon=K^{-2}$. Both carry the constant $\lambda M_0-I(\rho)$ for our $\rho$ and $M_0$.

*(c) The remaining steps.* Steps [F, (6.10)–(6.15)] do not involve $\rho$: Andréief's identity, the bound (6.11) for $w$, the sums (F8), and the estimate (6.15) for $S_K$. They give (6.14) with $\lambda M_0-I(\rho)$, and then (6.16) with $\overline U=\lambda M_0-I(\rho)+C_*$. Positivity of $F_K(\xi)$ is (F2). ∎

**Remark 1.** The statement of [F, Lemma 6.1] also requires every interval to have length $>1/225$. That condition is used only to guarantee (P3): by [F, p. 17], "the lower bound for the interval lengths controls the error in the regularization". Our measure below has $84$ intervals shorter than $1/225$, but satisfies (P3) with the value $9.548$.

**Remark 2 (Optimality of $\overline U$).** Let $\nu_{\rm eq}$ be the equilibrium measure of mass $\lambda$ in the external field $V$ on $[0,\infty)$. It minimizes $-I(\nu)+\int V\,\mathrm d\nu$ among positive measures of mass $\lambda$ on $[0,\infty)$ and has compact support [ST, Thm. I.1.3]. Let $F_{\rm eq}$ be its constant: $2U^{\nu_{\rm eq}}-V\le F_{\rm eq}$ quasi-everywhere on $[0,\infty)$, with equality on $\operatorname{supp}\nu_{\rm eq}$. Suppose $\rho$ satisfies (P1)–(P2) and $M_0\ge\sup_{t\ge0}(2U^\rho-V)$. Then $\int(2U^\rho-V)\,\mathrm d\nu_{\rm eq}\le\lambda M_0$, and $I(\rho-\nu_{\rm eq})\le0$ by [F, Lemma 6.2]. Hence

$$
\lambda M_0-I(\rho)\ge2I(\rho,\nu_{\rm eq})-\int V\,\mathrm d\nu_{\rm eq}-I(\rho)\ge I(\nu_{\rm eq})-\int V\,\mathrm d\nu_{\rm eq}=\lambda F_{\rm eq}-I(\nu_{\rm eq}).
$$

So every admissible test measure gives $\overline U\ge\overline U_{\rm eq}:=\lambda F_{\rm eq}-I(\nu_{\rm eq})+C_*$. Numerically [Chebyshev collocation for the one-interval problem; float], $\operatorname{supp}\nu_{\rm eq}=[7.1195\cdot10^{-5},\,0.753950]$, $F_{\rm eq}=-6.66078112$ and $\overline U_{\rm eq}=-1.3855098$.

### 3.3 The test measure $\rho_{1600}$

Numerically, the equilibrium measures $\nu_m$ of mass $m\le\lambda$ in the field $V$ have interval supports $S_m=[a(m),b(m)]$ that increase with $m$. Buyarov and Rakhmanov [BR] showed that, for external fields on the real line, $\nu_\lambda=\int_0^\lambda\omega_{S_m}\,\mathrm dm$, where $\omega_{S_m}$ is the unweighted equilibrium measure of $S_m$, i.e. the arcsine measure when $S_m$ is an interval. We use this only as a recipe: a midpoint discretization of the integral is a nested sum of arcsine measures. We used $1600$ nodes in the graded mass variable $m=\lambda(\frac32u^2-\frac12u^3)$, $u\in[0,1]$. The supports $S_m$ were computed in floating point, the data rounded to $10^{-15}\mathbb Z$, and the outermost weight adjusted so that the total mass is exactly $\lambda$. This construction is only heuristic and plays no role in the proof. Everything below applies to the rational data as given, in the file `data/16-mu5/rho1600.exact.json` (SHA-256 `0529b8fd…474068f`).

The data satisfy the following. Everything except the sum in (P3) is checked in exact rational arithmetic; that sum is evaluated in floating point in the main script and enclosed from above in ball arithmetic in the independent check (Section 8.4), with a huge margin:

- $J=1600$, and $a_j,b_j,c_j\in10^{-15}\mathbb Z$;
- (P1): $\min_jc_j=5.42\cdot10^{-7}>0$ and $\sum_jc_j=37/40$ exactly;
- (P2): the intervals are strictly nested, with innermost interval $[0.0059021556,\,0.0059389178]$ and support $S=[a_J,b_J]=[7.1243158827\cdot10^{-5},\,0.753301194422527]\subset(0,2)$;
- (P3): the sum equals $9.5471345675\ldots\le60$.

The shortest interval has length $3.68\cdot10^{-5}$ (see Remark 1). We write $t_S:=a_J=71243158827\cdot10^{-15}$.

### 3.4 The certified analytic constant

> **Proposition 3.** Let $\rho=\rho_{1600}$ and $M_0=-13321533/2000000=-6.6607665$. Then (8) holds:
> $$\sup_{0\le t\le2}\bigl(2U^\rho(t)-V(t)\bigr)<M_0,\qquad \tfrac{13}{10}\log2+3\alpha^3-2\pi\sqrt2=-7.98341\ldots<M_0 .$$
>
> Moreover
> $$I(\rho)\in-2.122676939604957409086598\pm4\cdot10^{-25},\qquad C_*\in2.653035990340488661129250\pm4\cdot10^{-25}.$$
>
> Hence $\overline U=\lambda M_0-I(\rho)+C_*\le-1.3854960825545539$. By Proposition 2, $\limsup_{K\to\infty}K^{-2}\log F_K(\xi)\le-1.3854960825545539$.

**Certification.** $I(\rho)$ and $C_*$ are evaluated from (7) and (6.3) in ball arithmetic on exact rational inputs. The resulting enclosure is $\overline U\in-1.385496082554553929784151\pm3\cdot10^{-25}$. The supremum is certified by `certify_sup` in `scripts/16-mu5/certify.py`, as follows.

The interval $[0,2]$ is first split at $0$, $2$, all $a_j,b_j$, and the endpoints of an interval $[q_-,q_+]$ of width $<10^{-13}$ that contains $t_V$. The bracket is certified by the signs of $\Phi'(\sqrt{q_\pm})$ in ball arithmetic. On a piece $[l,r]$ with no $a_j,b_j$ in its interior, put $c=(l+r)/2$, $h=r-l$ and $f=2U^\rho-V$. Two upper bounds are available:

$$
\sup_{[l,r]}f\le2\max\bigl(U^\rho(l),U^\rho(r)\bigr)-\min_{[l,r]}V, \tag{9}
$$

$$
\sup_{[l,r]}f\le f(c)+|f'(c)|\frac h2+\max\Bigl(0,-\min_{[l,r]}V''\Bigr)\frac{h^2}8. \tag{10}
$$

The first uses the quasi-convexity of $U^\rho$. Here $\min_{[l,r]}V$ equals $V(r)$ if $r\le q_-$ and $V(l)$ if $l\ge q_+$; otherwise it is bounded below by a ball enclosure of $V$ on $[q_-,q_+]$, by unimodality. The second holds because $U^\rho$ is concave on $[l,r]$ (it lies below its tangent at $c$) and by Taylor's theorem for $V$, with $V''$ from (5) enclosed on the ball $[l,r]$.

A piece is accepted when the smaller of the two bounds is certainly $<M_0$. Otherwise it is bisected. The final partition has $7733$ pieces. The largest accepted bound is $-6.66076650009$. The floating-point maximum of $f$ is $-6.66076682$, attained near $t=0.33599$ inside $S$, so $M_0$ has slack $3\cdot10^{-7}$. The inequality at $t=2$ is a direct ball evaluation. An independent implementation (Section 8.4) uses mean-value bounds instead of (9)–(10). It certifies the same supremum with $14156$ pieces and reproduces $I(\rho)$, $C_*$ and $\overline U$ to all printed digits.

---

## 4. The arithmetic constant at the cutoff $M=4\cdot10^8$

### 4.1 The limit function $R$

We recall [F, (5.4)–(5.6)]; the letters $T,s,q,H,b$ are local to this subsection. For $x\ge3$ and $0<z<\frac12$ let $\ell(x,z)=\lfloor x-z\rfloor+\lfloor x+z\rfloor+1$ and $b(x,z)=3\ell(\alpha x,z)$. Let $H=1+2\alpha=23/20$, $T=\lfloor2Hx\rfloor$, $s=Hx-T/2$, $q=\lfloor2x\rfloor$ and $n_+=(2x-q)/2$. Put

$$
\begin{aligned}
\Gamma_F(x)&=\int_0^{1/2}\bigl(T-b(x,z)\bigr)\bigl(T+b(x,z)-\ell(x,z)-5\bigr)\,\mathrm dz+s(2T-q-5)+(s-n_+)_+,\\
\mathcal N(x)&=2\lambda x\lfloor x\rfloor-12\lambda x\lfloor\alpha x\rfloor-2\mathcal J(\lambda x),\qquad \mathcal J(u)=mu-\frac{m(m+1)}4,\ m=\lfloor2u\rfloor,\\
R(x)&=-\Gamma_F(x)-\mathcal N(x).
\end{aligned}\tag{11}
$$

On each open interval between consecutive points of

$$
\mathfrak B=\{x:\ cx\in\mathbb Z\ \text{for some}\ c\in\{2,\,2\alpha,\,2\lambda,\,2H,\,4\alpha,\,2(1-\alpha),\,2(1+\alpha)\}\}
$$

all floor values and all branches of the minima in the closed form [F, (B.1)] are fixed. The quadratic coefficient of $R$ there is $2\lambda-12\lambda\alpha+2(H^2-H-\lambda^2)-18\alpha^2+6\alpha=0$, so $R$ is affine on each such interval [F, App. B.1]. If $R(x)=a_ix+b_i$ on $[t_i,t_{i+1}]$, then

$$
\int_{t_i}^{t_{i+1}}\frac{R(x)}{x^3}\,\mathrm dx=a_i\Bigl(\frac1{t_i}-\frac1{t_{i+1}}\Bigr)+\frac{b_i}2\Bigl(\frac1{t_i^2}-\frac1{t_{i+1}^2}\Bigr) \tag{12}
$$

[F, (B.3)], and $a_i,b_i$ are obtained exactly from two interior values of $R$.

### 4.2 The bound for $A_M$

For $3\le u<v\le\infty$ write $I_{u,v}=\int_u^vR(x)x^{-3}\,\mathrm dx$.

> **Proposition 4.** Let $M=4\cdot10^8$. Then
> $$\limsup_{K\to\infty,\ 40\mid K}K^{-2}\log m_{K,M}\le A_M:=A_\sharp+\frac{7\lambda}M-\frac{2923/240-1/4}{M^2}+\frac{32}{M^3},$$
>
> where
> $$A_\sharp:=I_{\rm out}+I_{3,20}+I_{20,4000}-\frac\lambda{4000}+\frac{37/2+13/16}{4000^2}\in\mathbb Q .$$
>
> Numerically, $I_{20,4000}=-0.0598432940998881484486091\ldots$ (an exact rational with a $3995$-digit denominator), $A_\sharp=1.3134544231730660672\ldots$ and $A_M=1.3134544393605659926\ldots$

*Proof.* *Step 1: tail bound.* We claim

$$
I_{4000,\infty}\le-\frac\lambda{4000}+\frac{37/2+13/16}{4000^2}=-\frac{58891}{256000000}. \tag{13}
$$

Let $\mathcal P(x)=74\{\alpha x\}(1-\{\alpha x\})-\lambda\{x\}(1-\{x\})$. This is Fauzan's $\mathcal P$ [F, p. 16]; note $74=6\lambda/\alpha$. It is continuous, $40$-periodic and piecewise $C^1$, with $\mathcal P'=F+\lambda$ off its breakpoints. Moreover $\mathcal P\le74/4=37/2$, and $\mathcal P(4000)=0$ because $4000$ and $\alpha\cdot4000=300$ are integers. By (F4), $R(x)x^{-3}=F(x)x^{-2}+Q(x)x^{-3}$. Integrating by parts,

$$
\int_{4000}^\infty\frac{F(x)}{x^2}\,\mathrm dx=\int_{4000}^\infty\frac{\mathcal P'(x)-\lambda}{x^2}\,\mathrm dx=-\frac{\mathcal P(4000)}{4000^2}+2\int_{4000}^\infty\frac{\mathcal P(x)}{x^3}\,\mathrm dx-\frac\lambda{4000}\le\frac{37/2}{4000^2}-\frac\lambda{4000},
$$

and $\int_{4000}^\infty Q(x)x^{-3}\,\mathrm dx\le\frac{13}8\cdot\frac1{2\cdot4000^2}$. This proves (13).

*Step 2: assembly.* By (F4), $R(x)x^{-3}=O(x^{-2})$, so $I_{3,\infty}$ converges. For $M\ge4000$,

$$
\int_3^M\frac{R(x)}{x^3}\,\mathrm dx=I_{3,20}+I_{20,4000}+I_{4000,\infty}-I_{M,\infty}.
$$

Insert (F6) for $I_{3,20}$, the exact value of $I_{20,4000}$ (Step 3), the bound (13), and (F5) for $I_{M,\infty}$ (valid since $40\mid4\cdot10^8$). Then $I_{\rm out}+6\lambda/M+\int_3^MR(x)x^{-3}\,\mathrm dx\le A_M$, and (F3) gives the claim.

*Step 3: the exact integral.* The interval $[20,4000]$ is cut at the points of $\mathfrak B$ into $33034$ pieces. On each, $R$ is affine, and $I_{20,4000}$ is the exact rational sum of the terms (12). Three independent programs computed it:

1. from the definitions (11), with the $z$-integral split exactly (`scripts/15-redteam-mu5/exact_R.py`; affinity checked at three extra points per piece);
2. from the closed form [F, (5.12)–(5.13), (B.1)] (`scripts/18-redteam-mu96/int_R_closed_form.py`; all $33034$ pieces verified affine);
3. from (11) again, with a different implementation (`scripts/09-zeta7det/exact_A.py`).

The rationals from (1) and (2) coincide exactly, and (3) agrees to all printed digits. Each program also reproduces (F6) exactly from the $143$ pieces of $[3,20]$; (1) also reproduces the $17$ rows of [F, Table 3] and $I_{\rm out}$ from [F, (5.8)–(5.10)]. The value is stored in `data/15-redteam-mu5/int_20_4000.txt` (SHA-256 `84c23d27…5cce5f956`). ∎

> **Corollary 1 (Decay).** Let $M=4\cdot10^8$, let $\overline U$ be as in Proposition 3, and put $\varepsilon_{\max}:=-(A_M+\overline U)$. Then
> $$\varepsilon_{\max}\in0.072041643193987937085\pm3.2\cdot10^{-23},\qquad\text{in particular}\quad \varepsilon_{\max}\ge0.0720416431939879 .$$
>
> Moreover $\limsup_{K\to\infty,\,40\mid K}K^{-2}\log Q_{K,M}(\xi)\le-\varepsilon_{\max}$. Hence for every $\varepsilon<\varepsilon_{\max}$ we have $0<Q_{K,M}(\xi)\le e^{-\varepsilon K^2}$ for all sufficiently large $K\in40\mathbb Z$.

*Proof.* We have $\log Q_{K,M}(\xi)=\log m_{K,M}+\log F_K(\xi)$. Apply Propositions 4 and 3, and use positivity from (F1). This is [F, (7.1)] with our constants. ∎

**Remark 3.**

1. Fauzan bounds $I_{20,\infty}\le-2689/48000=-0.0560208$ [F, (5.16)]. Our bound is $I_{20,4000}-58891/256000000=-0.0600733$. The true value is $I_{20,\infty}=-0.0600738$ [float; uses the $40$-periodicity below and Hurwitz zeta values]. Accordingly $A_\sharp=1.3134544$ improves Fauzan's $A_*=1.3175069$ [F, (5.19)] by $0.0040525$. The limit $A_\infty:=I_{\rm out}+I_{3,\infty}=1.3134540$ [float] is within $5\cdot10^{-7}$ of $A_\sharp$.
2. Any fixed $M\in40\mathbb Z$ with $M\ge4000$ could be used. The gain from $M=10^5$ to $M=4\cdot10^8$ is $6.47\cdot10^{-5}$, and the gain beyond $4\cdot10^8$ would be below $1.7\cdot10^{-8}$. The cutoff affects only $b_0$, through $K\ge200M^2$.
3. The function $Q=R-xF$ is $40$-periodic on $[20,\infty)$, with exact range $[0,2790/1849]$ over a period. Periodicity follows from two checks: the identity $R(y+40)-R(y)=40F(y)$, verified exactly on all $332$ pieces of $[20,60)$, and the cancellation of the $k^2$-coefficient of $R(40k+y)$, verified algebraically. In particular [F, (5.14)] holds for all $x\ge20$; it is also formalized in [DR] for all real $x$. We do not need any of this beyond (F4).
4. Using [F, (5.15)] at $T=4000$ instead of (13) would lower $A_\sharp$ by about $4\cdot10^{-7}$. We do not use it.

---

## 5. The norm bound

### 5.1 Notation

For $K\in40\mathbb Z_{>0}$ put $u_j=j/K$ and $x_j=-u_j^2$ for $1\le j\le K$, and

$$
\Omega(\tau)=\frac{\prod_{\ell=1}^N(\tau+u_\ell^2)^6}{\prod_{j=1}^K(\tau+u_j^2)}=\frac{\prod_{\ell=1}^N(\tau+u_\ell^2)^5}{\prod_{j=N+1}^K(\tau+u_j^2)} .
$$

This rational function has simple poles exactly at the points $x_j$, $N<j\le K$, all of which lie in $[-1,-\alpha^2)$. For $\tau>0$ let $W_K(\tau)=\Omega(\tau)\,w(K\sqrt\tau)/(2\sqrt\tau)>0$. For $\tau\in\mathbb C$ let

$$
\widetilde V(\tau)=\int_0^1\log|\tau+u^2|\,\mathrm du-6\int_0^\alpha\log|\tau+u^2|\,\mathrm du .
$$

This is continuous on $\mathbb C$ and harmonic on $\mathbb C\setminus[-1,0]$, and $V(t)=2\pi\sqrt t+\widetilde V(t)$ for $t\ge0$. Let $\psi(x)=x\log x-x$, with $\psi(0)=0$, and for $\alpha\le s\le1$

$$
\varphi(s)=6\bigl[\psi(s+\alpha)-\psi(s-\alpha)\bigr]-\psi(1+s)-\psi(1-s).
$$

Since $\int_0^a\log|u^2-s^2|\,\mathrm du=\psi(s+a)-\psi(s-a)$ for $s\ge a\ge0$, we have $\varphi(s)=-\widetilde V(-s^2)$ on $[\alpha,1]$.

Throughout this section $\rho$ is a test measure satisfying (P1)–(P2), with $S=\operatorname{supp}\rho=[a_J,b_J]$ and $t_S=a_J>0$. Put

$$
H_\rho(\tau)=2U^\rho(\tau)-\widetilde V(\tau),
$$

which is continuous on $\mathbb C$, harmonic on $\mathbb C\setminus([-1,0]\cup S)$, and satisfies $H_\rho(\bar\tau)=H_\rho(\tau)$. Let $m_0$ be a real number with

$$
m_0\le\inf_{t\in S}\bigl(2U^\rho(t)-V(t)\bigr), \tag{14}
$$

and put $c_0=\log((2\pi)^4/24)-2=2.17345\ldots$ and

$$
E(K)=\frac{3\pi}{2K\sqrt{t_S}}-3\log K-2\log t_S-c_0 .
$$

### 5.2 The residue and contour identities

> **Lemma 1.** Let $q\in\mathbb R[t]$ with $\deg q<h$, and put $p(\tau)=q(K^2\tau)$.
>
> (a) $q^{\mathsf T}Gq=K^{12N-2K+1}\int_0^\infty p(\tau)^2\,W_K(\tau)\,\mathrm d\tau$.
>
> (b) $q^{\mathsf T}M_Kq=K^{12N-2K+6}\sum_{j=N+1}^Kr_j\,p(x_j)^2$, where $r_j:=\operatorname{Res}_{\tau=x_j}\tau^2\Omega(\tau)$. Moreover, let $\Gamma$ be a positively oriented rectifiable Jordan curve such that every $x_j$, $N<j\le K$, lies in its interior. (This holds, for instance, if $\Gamma\cap[-1,0]=\emptyset$ and $[-1,0]$ lies inside $\Gamma$.) Then
> $$q^{\mathsf T}M_Kq=K^{12N-2K+6}\,\frac1{2\pi i}\oint_\Gamma\tau^2\,\Omega(\tau)\,p(\tau)^2\,\mathrm d\tau .$$
>
> (c) The norm is given by
> $$\lVert B_K\rVert=K^5\sup_p\frac{\bigl|\sum_jr_j\,p(x_j)^2\bigr|}{\int_0^\infty p^2W_K}=K^5\sup_p\frac{\bigl|\frac1{2\pi i}\oint_\Gamma\tau^2\Omega\,p^2\,\mathrm d\tau\bigr|}{\int_0^\infty p^2W_K},$$
> with suprema over nonzero real polynomials $p$ of degree $<h$.

*Proof.* (a) By (F2), $q^{\mathsf T}Gq=\int_0^\infty D_N(y^2)^6q(y^2)^2w(y)D_K(y^2)^{-1}\,\mathrm dy$. Substitute $y=K\sqrt\tau$, so $\mathrm dy=K\,\mathrm d\tau/(2\sqrt\tau)$. Since $D_m(K^2\tau)=K^{2m}\prod_{j\le m}(\tau+u_j^2)$, we get $D_N(K^2\tau)^6/D_K(K^2\tau)=K^{12N-2K}\,\Omega(\tau)$, and $q(K^2\tau)=p(\tau)$.

(b) The rational function $R_q:=D_N^6q^2/D_K=D_N^5q^2/(D_K/D_N)$ has a polynomial part and at most simple poles, located at the points $-j^2$ with $N<j\le K$. By the definition of $\mu_X$ through partial fractions and by (3),

$$
q^{\mathsf T}M_Kq=[X]\,\mu_X(R_q)=\sum_{j=N+1}^Kj^4\operatorname{Res}_{t=-j^2}R_q=\sum_{j=N+1}^K\operatorname{Res}_{t=-j^2}\bigl(t^2R_q\bigr),
$$

since $t^2=j^4$ at $t=-j^2$. The only poles of the rational function $t^2R_q$ are the points $-j^2$. By the residue theorem, $q^{\mathsf T}M_Kq=\frac1{2\pi i}\oint_{K^2\Gamma}t^2R_q(t)\,\mathrm dt$. Substituting $t=K^2\tau$ turns $t^2R_q(t)\,\mathrm dt$ into $K^{4+12N-2K+2}\,\tau^2\Omega(\tau)p(\tau)^2\,\mathrm d\tau$. The residue form follows from $\operatorname{Res}_{\tau=x_j}(\tau^2\Omega p^2)=r_jp(x_j)^2$.

(c) $B_K$ is real symmetric, so $\lVert B_K\rVert=\sup_{x\ne0}|x^{\mathsf T}B_Kx|/x^{\mathsf T}x$. Substituting $x=G^{1/2}q$ gives $\lVert B_K\rVert=\sup_q|q^{\mathsf T}M_Kq|/q^{\mathsf T}Gq$, and $q\mapsto p$ is a bijection onto the real polynomials of degree $<h$. Now apply (a) and (b). ∎

Explicitly, $r_j=u_j^4\prod_{\ell\le N}(u_\ell^2-u_j^2)^5\big/\prod_{N<k\le K,\,k\ne j}(u_k^2-u_j^2)$, which is [F, (C.3)].

### 5.3 Weighted Bernstein–Walsh and Nikolskii inequalities

> **Lemma 2 (Weighted Bernstein–Walsh).** For every nonzero $p\in\mathbb C[\tau]$ with $\deg p<\lambda K$ and every $z\in\mathbb C$,
> $$\log|p(z)|\le K\,U^\rho(z)+\max_{t\in S}\bigl(\log|p(t)|-K\,U^\rho(t)\bigr).$$

*Proof.* Let $u=\log|p|-KU^\rho$. By (6), $U^\rho$ is a finite positive combination of continuous potentials. So $U^\rho$ is continuous on $\mathbb C$ and harmonic on $\mathbb C\setminus S$, and $U^\rho(z)=\lambda\log|z|+o(1)$ as $z\to\infty$ by (P1). Hence $u$ is subharmonic on $\mathbb C\setminus S$ and upper semicontinuous on $\mathbb C$. Moreover $u(z)\le(\deg p-\lambda K)\log|z|+O(1)\to-\infty$ as $z\to\infty$. The maximum $m^*=\max_Su$ exists by upper semicontinuity, and $m^*>-\infty$. Choose $r_0>0$ with $S\subset\{|z|<r_0\}$ and $u\le m^*$ on $\{|z|\ge r_0\}$. The domain $D=\{|z|<r_0\}\setminus S$ is bounded and connected, since $S$ is a segment. For $\zeta\in S$, $\limsup_{z\to\zeta}u(z)\le u(\zeta)\le m^*$, and $u\le m^*$ on $|z|=r_0$. The maximum principle for subharmonic functions [Ra, Thm. 2.3.1] gives $u\le m^*$ on $D$, hence on $\mathbb C$. ∎

> **Lemma 3 (Nikolskii-type bound).** For every $t>0$ and every real polynomial $p$ of degree $<h$,
> $$p(t)^2\le h^2K^2\exp\Bigl(KV(t)+\frac{3\pi}{2K\sqrt t}-3\log K-2\log t-c_0\Bigr)\int_0^\infty p^2\,W_K .$$

*Proof.* (i) *Lower bound for the weight.* For $\tau>0$, the left inequality of (F8) with $m=N$ gives $6\sum_{\ell\le N}\log(\tau+u_\ell^2)\ge6K\int_0^\alpha\log(\tau+u^2)\,\mathrm du$. The right inequality with $m=K$ gives $\sum_{j\le K}\log(\tau+u_j^2)\le K\int_0^1\log(\tau+u^2)\,\mathrm du+2\log K+2$. Hence $\log\Omega(\tau)\ge-K\widetilde V(\tau)-2\log K-2$. Keeping only the term $\ell=1$ in $w$ gives $w(K\sqrt\tau)\ge(2\pi)^4K^5\tau^{5/2}e^{-2\pi K\sqrt\tau}/12$. Therefore

$$
\log W_K(\tau)\ge-KV(\tau)+3\log K+2\log\tau+c_0\qquad(\tau>0).
$$

(ii) *Growth of $V$.* By (5), $V'(v)\le(\pi+\pi/2)/\sqrt v=3\pi/(2\sqrt v)$. So $V(\tau)\le V(t)+3\pi(\tau-t)/(2\sqrt t)$ for $\tau\ge t$.

(iii) *Legendre expansion.* Let $L>0$ and let $P$ be a polynomial of degree $\le h-1$. Expand $P$ in the orthonormal Legendre basis $\phi_k(x)=\sqrt{(2k+1)/L}\,P_k(2(x-t)/L-1)$ of $L^2[t,t+L]$. Since $|P_k|\le1$ on $[-1,1]$ [DLMF, 18.14.1], Cauchy–Schwarz gives

$$
\max_{[t,t+L]}P^2\le\frac{\sum_{k<h}(2k+1)}L\int_t^{t+L}P^2=\frac{h^2}L\int_t^{t+L}P^2 .
$$

Take $L=K^{-2}$. For $\tau\in[t,t+L]$, (i) and (ii) give $W_K(\tau)\ge\exp\bigl(-KV(t)-3\pi/(2K\sqrt t)+3\log K+2\log t+c_0\bigr)$. Hence

$$
p(t)^2\le h^2K^2\int_t^{t+L}p^2\le h^2K^2e^{KV(t)+3\pi/(2K\sqrt t)-3\log K-2\log t-c_0}\int_0^\infty p^2W_K. \qquad\blacksquare
$$

> **Corollary 2.** For every real polynomial $p$ of degree $<h$ and every $z\in\mathbb C$,
> $$|p(z)|^2\le h^2K^2e^{E(K)}\,e^{K(2U^\rho(z)-m_0)}\int_0^\infty p^2W_K .$$

*Proof.* We may assume $p\ne0$. Let $t^*\in S$ maximize $\log|p|-KU^\rho$ on $S$. By Lemma 2,

$$
2\log|p(z)|\le2KU^\rho(z)+2\log|p(t^*)|-2KU^\rho(t^*).
$$

Apply Lemma 3 at $t^*$. Then use $V(t^*)-2U^\rho(t^*)\le-m_0$ from (14), and $t^*\ge t_S$, so that $3\pi/(2K\sqrt{t^*})\le3\pi/(2K\sqrt{t_S})$ and $-2\log t^*\le-2\log t_S$. ∎

### 5.4 Residue moduli and the real-variable bound

> **Lemma 4.** For $N<j\le K$ we have $\operatorname{sgn}r_j=(-1)^{j-1}$,
> $$|r_j|=\frac{2K^{2K-12N-6}\,(j+N)!^6}{(j-N-1)!^6\,(K-j)!\,(K+j)!},\qquad\text{and}\qquad \log|r_j|\le K\varphi(u_j)+\log2+\frac12+3\log(4\pi K).$$

*Proof.* For $j>N$, each of the $N$ factors $u_\ell^2-u_j^2$ is negative, and $u_k^2-u_j^2<0$ exactly for the $j-N-1$ indices $N<k<j$. This gives the sign. For the modulus, use

$$
\prod_{\ell\le N}(j^2-\ell^2)=\frac{(j-1)!\,(j+N)!}{(j-N-1)!\,j!},\qquad \prod_{N<k\le K,\,k\ne j}|k^2-j^2|=\frac{(j-N-1)!\,(K-j)!\,(K+j)!}{(N+j)!\cdot2j}
$$

and collect powers of $K$.

For the bound, write $\log m!=\psi(m)+\theta(m)$, where $0\le\theta(m)\le\frac12\log(2\pi m)+\frac1{12m}$ for $m\ge1$ [Ro] and $\theta(0)=0$. Since $\psi(Kx)=K\psi(x)+Kx\log K$, the coefficient of $\log K$ in $6\psi(j+N)-6\psi(j-N-1)-\psi(K-j)-\psi(K+j)$ is $12N+6-2K$. This cancels the prefactor $K^{2K-12N-6}$ exactly. With $s=u_j$, the remaining main term is

$$
K\bigl[6\psi(s+\alpha)-6\psi(s-\alpha-\tfrac1K)-\psi(1-s)-\psi(1+s)\bigr]=K\varphi(s)+6K\bigl[\psi(s-\alpha)-\psi(s-\alpha-\tfrac1K)\bigr].
$$

Since $\psi$ is convex with $\psi'=\log$, the last term is $\le6\log(s-\alpha)\le0$. Finally, $-6\theta(j-N-1)-\theta(K-j)-\theta(K+j)\le0$ and $6\theta(j+N)\le3\log(4\pi K)+\frac12$, since $j+N<2K$. ∎

> **Proposition 5 (Real-variable bound).** Let $M_1\ge\sup_{s\in[\alpha,1]}\bigl(2U^\rho(-s^2)+\varphi(s)\bigr)$. Then for every $K\in40\mathbb Z_{>0}$
> $$\log\lVert B_K\rVert\le K(M_1-m_0)+10\log K+\frac{3\pi}{2K\sqrt{t_S}}+C_\rho,\qquad C_\rho:=3\log(4\pi)-2\log t_S-c_0+\log2+\frac12 .$$

*Proof.* By Lemma 1(c), $\lVert B_K\rVert\le K^5\sum_j|r_j|\sup_pp(x_j)^2/\int p^2W_K$. Bound each term by Corollary 2 at $z=x_j=-u_j^2$ and by Lemma 4. There are $h\le K$ terms, and $u_j\in(\alpha,1]$. ∎

### 5.5 The contour bound

> **Lemma 5.** Let $\tau\in\mathbb C$ with $d(\tau):=\operatorname{dist}(\tau,[-1,0])>0$. Then
> $$\log|\Omega(\tau)|\le-K\widetilde V(\tau)+\frac{1+6\alpha^2}{d(\tau)} .$$

*Proof.* Let $f(u)=\log|\tau+u^2|$ for $u\in[0,1]$. Since $|\tau+u^2|=|\tau-(-u^2)|\ge d(\tau)$, we have $|f'(u)|\le2u/|\tau+u^2|\le2u/d(\tau)$. So the total variation of $f$ on $[0,a]$ is at most $a^2/d(\tau)$. A right-endpoint Riemann sum with step $1/K$ differs from $K$ times the integral by at most the total variation. Indeed $\bigl|\sum_{j\le m}f(j/K)-K\int_0^{m/K}f\bigr|\le\sum_{j\le m}K\int_{(j-1)/K}^{j/K}|f(j/K)-f(u)|\,\mathrm du\le\sum_{j\le m}\operatorname{Var}_{[(j-1)/K,\,j/K]}f$. Hence $6\sum_{\ell\le N}f(u_\ell)\le6K\int_0^\alpha f+6\alpha^2/d$ and $\sum_{j\le K}f(u_j)\ge K\int_0^1f-1/d$. Subtracting the second from the first gives the claim. ∎

> **Proposition 6 (Contour bound).** Let $\Gamma$ be a rectifiable Jordan curve with $\Gamma\cap[-1,0]=\emptyset$ and $[-1,0]$ in its interior, and put $d_\Gamma=\min_{\tau\in\Gamma}d(\tau)$. Then for every $K\in40\mathbb Z_{>0}$,
> $$\log\lVert B_K\rVert\le K\Bigl(\max_{\tau\in\Gamma}H_\rho(\tau)-m_0\Bigr)+4\log K+2\log h+\frac{3\pi}{2K\sqrt{t_S}}+C_\Gamma,$$
>
> where
> $$C_\Gamma:=-2\log t_S-c_0+\log\frac{\operatorname{len}(\Gamma)\max_{\tau\in\Gamma}|\tau|^2}{2\pi}+\frac{1+6\alpha^2}{d_\Gamma}.$$

*Proof.* By Lemma 1(b),(c) and $|\oint_\Gamma g|\le\operatorname{len}(\Gamma)\max_\Gamma|g|$,

$$
\lVert B_K\rVert\le K^5\,\frac{\operatorname{len}(\Gamma)}{2\pi}\sup_p\max_{\tau\in\Gamma}\frac{|\tau|^2|\Omega(\tau)|\,|p(\tau)|^2}{\int_0^\infty p^2W_K}.
$$

For $\tau\in\Gamma$, Lemma 5 gives $|\Omega(\tau)|\le e^{-K\widetilde V(\tau)+(1+6\alpha^2)/d_\Gamma}$. Corollary 2 gives $|p(\tau)|^2\le h^2K^2e^{E(K)}e^{K(2U^\rho(\tau)-m_0)}\int p^2W_K$. Multiplying, the exponent is $K(H_\rho(\tau)-m_0)$. The powers of $K$ combine to $5+2-3=4$. ∎

### 5.6 The contour and the saddle point

**Why moving the contour helps.** On the pole segment, $H_\rho(-s^2)=2U^\rho(-s^2)+\varphi(s)$, and Proposition 5 is the degenerate case in which $\Gamma$ collapses onto the segment. We can write

$$
\widetilde V(\tau)=\int_\alpha^1\log|\tau+u^2|\,\mathrm du-5\int_0^\alpha\log|\tau+u^2|\,\mathrm du .
$$

The first term is the logarithmic potential of the image of $\mathrm du|_{[\alpha,1]}$ under $u\mapsto-u^2$. That image is the limiting distribution of the poles $x_j$, each with mass $1/K$, and is supported on $[-1,-\alpha^2]$. Its potential is subharmonic, and its normal derivative jumps across the segment. Since it enters $H_\rho$ with a minus sign, $H_\rho$ has a ridge along $[-1,-\alpha^2]$ and decreases at first order when $\Gamma$ is lifted off it. This is where the alternating residues cancel. The contour must still cross the real axis to the left of $-1$, because the ray $(-\infty,-1)$ joins a point inside $\Gamma$ to infinity. On $(-\infty,-1)$, $H_\rho(x)\sim(2\lambda-1+6\alpha)\log|x|\to+\infty$ as $x\to-\infty$, and $H_\rho$ attains a minimum. Since $H_\rho$ is harmonic there and symmetric under conjugation, that minimum is a saddle point: $H_\rho$ increases along the real axis and decreases along the vertical line. The optimal contour crosses the real axis at this saddle, beyond the end $x_K=-1$ of the pole range.

**The contour used.** Let $\tau_*=-1069911/10^6=-1.069911$. Let $\Gamma$ be the boundary of the rectangle $[\tau_*,1]\times[-\frac15,\frac15]$, and $\Gamma^+$ its part in the closed upper half-plane: the polygon $\tau_*\to\tau_*+\frac i5\to1+\frac i5\to1$. Then

$$
d_\Gamma=|\tau_*+1|=0.069911,\qquad \operatorname{len}(\Gamma)=4.939822,\qquad \max_\Gamma|\tau|^2=\tau_*^2+\tfrac1{25}=1.184709547921 .
$$

The real minimum of $H_\rho$ on $(-\infty,-1)$ lies at $x\approx-1.06991015$, within $10^{-6}$ of $\tau_*$. A ball evaluation gives $H_\rho(\tau_*)=0.80880185128432719515\ldots$. Floating-point scans show that along $\Gamma^+$, $H_\rho$ decreases going up the left edge, is at most $0.751$ on the top edge, and is at most $-0.537$ on the right edge. These scans are not used in the proof.

> **Proposition 7 (Certified constants).** For $\rho=\rho_{1600}$ and the contour $\Gamma$ above:
>
> (a) $\inf_{t\in S}\bigl(2U^\rho(t)-V(t)\bigr)>m_0:=-3330407/500000=-6.660814$;
>
> (b) $\max_{\tau\in\Gamma}H_\rho(\tau)<M_\Gamma:=404401/500000=0.808802$;
>
> (c) $\sup_{s\in[\alpha,1]}\bigl(2U^\rho(-s^2)+\varphi(s)\bigr)<M_1:=8691543/10^7=0.8691543$.
>
> Hence $\varrho_2:=M_\Gamma-m_0=466851/62500=7.469616$ and $\varrho_1:=M_1-m_0=75299683/10^7=7.5299683$.

**Certification.** All three statements are certified in `scripts/16-mu5/mu_final16.py` on exact rational data, in ball arithmetic at $110$–$200$ bits.

*(a)* [`certify_inf`]. $S$ is split at all $a_j,b_j$. On a breakpoint-free piece $[l,r]$, with $c,h$ as before, $f=2U^\rho-V$ is bounded below by the larger of two quantities. The first is $2\min(U^\rho(l),U^\rho(r))-\max(V(l),V(r))$, which uses concavity of $U^\rho$ and unimodality of $V$. The second is

$$
U^\rho(l)+U^\rho(r)-V(c)-\Bigl|\frac{2(U^\rho(r)-U^\rho(l))}h-V'(c)\Bigr|\frac h2-\max\Bigl(0,\max_{[l,r]}V''\Bigr)\frac{h^2}8 ,
$$

which holds because a concave function lies above its chord, together with Taylor's theorem for $V$. The final partition has $3412$ pieces, with smallest accepted bound $-6.66081399675$. The floating-point infimum is $-6.6608136014$, attained at a breakpoint near $t=0.334225$.

*(b)* [`contour_certify.py`]. $\Gamma^+$ is cut into axis-parallel segments with midpoint $\tau_0$ and half-length $\ell$. Each segment must satisfy $\operatorname{dist}(\tau_0,[-1,0])>2\ell$ and $\operatorname{dist}(\tau_0,S)>2\ell$; otherwise it is bisected. On the disc $D=\{|\tau-\tau_0|\le\ell\}$ we have $H_\rho=\operatorname{Re}\mathcal F$ with $\mathcal F$ analytic, and

$$
\mathcal F'(\tau)=2\sum_j\frac{c_j}{\sqrt{\tau-a_j}\sqrt{\tau-b_j}}-\Bigl(\chi(1,\tau)-6\chi(\alpha,\tau)\Bigr),\qquad \chi(a,\tau)=\int_0^a\frac{\mathrm du}{\tau+u^2}=\frac{\arctan(a/\sqrt\tau)}{\sqrt\tau}.
$$

Then $H_\rho\le H_\rho(\tau_0)+|\mathcal F'(\tau_0)|\ell+\frac12\ell^2\sup_D|\mathcal F''|$ on $D$, with

$$
\sup_D|\mathcal F''|\le\sum_jc_j\Bigl(d_a^{-3/2}d_b^{-1/2}+d_a^{-1/2}d_b^{-3/2}\Bigr)+\frac{1+6\alpha}{d_0^2},
$$

where $d_a=|\tau_0-a_j|-\ell$, $d_b=|\tau_0-b_j|-\ell$ and $d_0=\operatorname{dist}(\tau_0,[-1,0])-\ell$. Here $U^\rho(\tau_0)$ is evaluated by (6), and $\widetilde V(\tau_0)$ from $\int_0^a\log|u-\kappa|\,\mathrm du=\operatorname{Re}[(a-\kappa)\operatorname{Log}(a-\kappa)+\kappa\operatorname{Log}(-\kappa)]-a$ with $\kappa=\pm i\sqrt{\tau_0}$, since $\tau+u^2=(u-i\sqrt\tau)(u+i\sqrt\tau)$. The final partition has $215$ pieces, and the largest accepted bound is $0.80880196112$. By conjugation symmetry, the bound on $\Gamma^+$ holds on $\Gamma$. The closed forms were checked against direct quadrature at $11$ points, with agreement $\le4\cdot10^{-16}$ (`check_complex_formulas.py`). The value $H_\rho(\tau_*)$ from the complex formulas agrees with the real formula $2U^\rho(\tau_*)-\psi(s+1)+\psi(s-1)+6[\psi(s+\alpha)-\psi(s-\alpha)]$, $s=\sqrt{-\tau_*}$, to $2\cdot10^{-16}$.

*(c)* [`certify_neg`]. On $[l,r]\subset[\alpha,1]$, $U^\rho(-s^2)$ is nondecreasing in $s$. $\varphi$ is concave, since $\varphi''(s)=-12\alpha/(s^2-\alpha^2)-2/(1-s^2)<0$. The function $\mathcal G(s)=2U^\rho(-s^2)+\varphi(s)$ satisfies $\mathcal G''\le4\lambda/s^2$. The two bounds used are

$$
\mathcal G\le2U^\rho(-r^2)+\min\bigl(\varphi(l)+\max(0,\varphi'(l))h,\ \varphi(r)+\max(0,-\varphi'(r))h\bigr)
$$

and $\mathcal G\le\mathcal G(c)+|\mathcal G'(c)|h/2+(4\lambda/l^2)h^2/8$. There are $78$ pieces. Statement (c) is used only for Proposition 5, not for Theorem 1.

Statements (a) and (b) were re-certified by independent code with different schemes: (a) with $12585$ pieces, (b) with $73$ pieces using ball enclosures of $\mathcal F'$ on boxes and explicitly chosen analytic branches (Section 8.4).

### 5.7 Proof of Theorem 2

Apply Proposition 6 with $\rho=\rho_{1600}$, the contour $\Gamma$ of Section 5.6, and $m_0$, $M_\Gamma$ from Proposition 7. Then $K(\max_\Gamma H_\rho-m_0)\le\varrho_2K$. With $t_S=71243158827\cdot10^{-15}$, we have $-2\log t_S=19.09882\ldots$ and $c_0=2.17345\ldots$. Moreover $\log(\operatorname{len}(\Gamma)\max_\Gamma|\tau|^2/(2\pi))=-0.07105\ldots$ and $(1+6\alpha^2)/d_\Gamma=1.03375/0.069911=14.78666\ldots$. So $C_\Gamma=31.64098\ldots<31.65$ and $3\pi/(2\sqrt{t_S})=558.302\ldots<559$. ∎

In the same way, Proposition 5 with Proposition 7(a),(c) gives $\log\lVert B_K\rVert\le\varrho_1K+10\log K+25.72+559/K$, with $\varrho_1=7.5299683$.

**Table 3.** The explicit bound of Theorem 2, evaluated in floating point, against the numerically computed norm (Section 7.1). The gap comes from the non-exponential terms of the bound and the $-\frac12\log K$ correction in the true norm.

| $K$ | $\log\lVert B_K\rVert$ (numerical) | $K^{-1}\log\lVert B_K\rVert$ | Theorem 2 bound | gap |
|---:|---:|---:|---:|---:|
| 40 | 292.5613171408 | 7.31403293 | 366.36 | 73.8 |
| 80 | 591.0372993643 | 7.38796624 | 662.33 | 71.3 |
| 160 | 1188.2728511251 | 7.42670532 | 1260.56 | 72.3 |
| 320 | 2383.0664002111 | 7.44708250 | 2458.12 | 75.1 |
| 640 | 4772.9897927749 | 7.45779655 | 4851.68 | 78.7 |
| 1280 | 9553.1780391228 | 7.46342034 | 9635.96 | 82.8 |

---

## 6. Proof of the main theorem

*Proof of Theorem 1.* Let $M=4\cdot10^8$. By Corollary 1, hypothesis (a) of Proposition 1 holds for every $\varepsilon<\varepsilon_{\max}$. By Theorem 2, hypothesis (b) holds with $\varrho=\varrho_2$, $A=6$ (since $h\le K$) and $C=e^{31.65+559/40}$. Given $\mu>\mu_0:=\varrho_2\lambda/\varepsilon_{\max}$, choose $\varepsilon<\varepsilon_{\max}$ so close to $\varepsilon_{\max}$ that $\varrho_2\lambda/\varepsilon<\mu$. Proposition 1 then gives $b_0(\mu)$.

It remains to bound $\mu_0$. With the lower end of the ball for $\varepsilon_{\max}$ (Corollary 1),

$$
\mu_0=\frac{\varrho_2\lambda}{\varepsilon_{\max}}\le95.908345418981323384+5\cdot10^{-19}<95.9083454190 . \qquad\blacksquare
$$

**An explicit instance.** For $\mu=95.91$ the choice in the proof of Proposition 1 can be made explicit. Take

$$
\varepsilon=\frac{9005205399}{125\cdot10^9}=0.072041643192<\varepsilon_{\max},\qquad c=\frac{513596269}{4\cdot10^7}=12.839906725 .
$$

In exact rational arithmetic,

$$
\varepsilon c-\lambda=\frac{39894505056331}{5\cdot10^{18}}\approx7.979\cdot10^{-6}>0,\qquad 95.91-\varrho_2c=\frac{2068221081}{25\cdot10^{11}}\approx8.273\cdot10^{-4}>0 .
$$

The admissible interval for $c$ is $(\lambda/\varepsilon_{\max},\,95.91/\varrho_2)=(12.8397960,\,12.8400175)$. For comparison, Fauzan's choice gives $\varrho c-260=-1/80$ and $\varepsilon-\lambda/c=1/24000$.

**Safe rounding.** Every constant is rounded in the direction that keeps the conclusion valid.

- $M_0$, $M_\Gamma$ and $M_1$ are rational targets that the computed supremum is certified to lie strictly *below*; $m_0$ is a target that the infimum is certified to lie strictly *above*.
- $I(\rho)$, $C_*$, $\overline U$, $\varepsilon_{\max}$ and $\mu_0$ are carried as balls. We use the upper end for $\overline U$ and $\mu_0$, and the lower end for $\varepsilon_{\max}$.
- $A_M$ and $\varrho_2=M_\Gamma-m_0$ are exact rationals.
- The comparisons $\varepsilon<\varepsilon_{\max}$, $\varepsilon c>\lambda$ and $\varrho_2c<\mu$ are exact or ball-certain.

The margin of $95.91$ over $\mu_0$ is $0.0016546$. It would absorb an error of up to $1.2\cdot10^{-6}$ in $\varepsilon_{\max}$ or $1.3\cdot10^{-4}$ in $\varrho_2$. The two independent sets of certificates agree far more closely than that (Section 8.4).

**Remark 4 (By-products).** With the same constants, (7.1) of [F] gives two decay rates for all large $n$. For $M=4\cdot10^8$, $-1600(A_M+\overline U)\ge115.266$, so $Q_{40n,M}(\zeta(5))<e^{-115.26\,n^2}$. For $M=200$, $A_{200}=1.34553519400640\ldots$ by the same route, so $-1600(A_{200}+\overline U)\ge63.937$ and $Q_{40n,200}(\zeta(5))<e^{-63.9\,n^2}$, against $e^{-139n^2/5}$ in [F, (2.7)].

---

## 7. Sharpness

### 7.1 The true size of $\lVert B_K\rVert$

Let $\beta_1,\dots,\beta_h$ be the eigenvalues of $B_K$. Then $\Delta_K(X)=\Delta_K(\xi)\prod_j(1+(X-\xi)\beta_j)$, as in [F, App. C.3]. Since $\deg\Delta_K=h$, all $\beta_j\ne0$, so the roots of $\Delta_K$ are the real numbers $\xi-1/\beta_j$, and

$$
\lVert B_K\rVert=\frac1{\min_j|X_j-\zeta(5)|}\quad\text{over the roots } X_j \text{ of } \Delta_K .
$$

We computed $\log\lVert B_K\rVert$ for $K\le1280$ by three methods [NUMERICAL]:

1. for $K=40,80$, from the exact characteristic polynomial $\Delta_K(X)=\det M_K\cdot\det(XI+M_K^{-1}G_K(0))$ and its real roots isolated with PARI/GP at $520$–$840$ digits;
2. by a float64 Lanczos method with scaling, for $40\le K\le320$;
3. by a Stieltjes procedure in mpmath at $40$–$60$ digits, with orthonormal polynomials for the weight $W_K$. There are two independent implementations of this method, one for $40\le K\le640$ and one for $320\le K\le1280$.

Discretizations were checked on refined grids. All methods agree to all printed digits wherever they overlap (Tables 3 and 4).

**Table 4.** Numerically computed $\log\lVert B_K\rVert$; the values for $K=40,80$ come from exact characteristic polynomials.

| $K$ | $\log\lVert B_K\rVert$ | $K^{-1}\log\lVert B_K\rVert$ |
|---:|---:|---:|
| 40 | 292.5613171408 | 7.31403293 |
| 80 | 591.0372993643 | 7.38796624 |
| 160 | 1188.2728511251 | 7.42670532 |
| 240 | 1785.6413312423 | 7.44017221 |
| 320 | 2383.0664002111 | 7.44708250 |
| 480 | 3577.9992183777 | 7.45416504 |
| 640 | 4772.9897927749 | 7.45779655 |
| 960 | 7163.0547545318 | 7.46151537 |
| 1280 | 9553.1780391228 | 7.46342034 |

**Fits.** We fit $\log\lVert B_K\rVert\approx RK-c\log K-d$ to the data with $K\ge K_0$. With $c$ free, $K_0=160,240,320,480,640$ give $R=7.469573$, $7.469576$, $7.469578$, $7.469579$, $7.469580$, with $c\to0.495$. With $c=\frac12$ fixed, $K_0=320,640,960$ give $R=7.469588$, $7.469586$, $7.469585$. The difference quotients corrected by $\frac12\log K$ decrease monotonically to $7.469585$. So the two families bracket

$$
\lim_{K\to\infty}K^{-1}\log\lVert B_K\rVert\approx7.46958,
$$

consistent with a saddle-point asymptotic $\log\lVert B_K\rVert\approx RK-\frac12\log K-4.3$. The certified $\varrho_2=7.469616$ lies above the observed rate by about $3.4\cdot10^{-5}$. Accordingly the deficit $\varrho_2K-\log\lVert B_K\rVert$ increases with $K$, from $6.22$ at $K=40$ to $7.93$ at $K=1280$.

### 7.2 The limit of the method

> **Proposition 8.** Let $\nu_{\rm eq}$ and $F_{\rm eq}$ be as in Remark 2. Let $\rho$ satisfy (P1)–(P2), and let $m_0$ satisfy (14). Then $U^\rho\ge U^{\nu_{\rm eq}}-\frac12(F_{\rm eq}-m_0)$ on $\mathbb C$. Consequently:
>
> (a) for every Jordan curve $\Gamma$ as in Proposition 6,
> $$\max_\Gamma H_\rho-m_0\ \ge\ \varrho_2^{\rm eq}:=\inf_{x<-1}\bigl(2U^{\nu_{\rm eq}}(x)-\widetilde V(x)\bigr)-F_{\rm eq};$$
>
> (b) $\sup_{s\in[\alpha,1]}\bigl(2U^\rho(-s^2)+\varphi(s)\bigr)-m_0\ \ge\ \varrho_1^{\rm eq}:=\sup_{s\in[\alpha,1]}\bigl(2U^{\nu_{\rm eq}}(-s^2)+\varphi(s)\bigr)-F_{\rm eq}$.

*Proof.* On $S$, $2U^\rho\ge V+m_0$. Also $V\ge2U^{\nu_{\rm eq}}-F_{\rm eq}$ quasi-everywhere on $[0,\infty)$, hence $\rho$-almost everywhere, because $\rho$ has finite energy and therefore charges no set of capacity zero. So $2U^\rho\ge2U^{\nu_{\rm eq}}-F_{\rm eq}+m_0$ holds $\rho$-almost everywhere. In the sign convention of [ST] this reads $-U^\rho\le-U^{\nu_{\rm eq}}+\frac12(F_{\rm eq}-m_0)$, $\rho$-a.e. Both measures have mass $\lambda$ and $\rho$ has finite energy, so the domination principle [ST, Thm. II.3.2] extends the inequality to all of $\mathbb C$. For (a), $\Gamma$ meets $(-\infty,-1)$ because $-1$ lies inside $\Gamma$. For (b), evaluate at $-s^2$. ∎

Numerically [one-interval Chebyshev collocation for $\nu_{\rm eq}$; float], $\varrho_1^{\rm eq}=7.529935$, attained at $s\approx0.97318$, and $\varrho_2^{\rm eq}=7.469583$, attained at the saddle $x\approx-1.0699$. Our certified constants exceed these by $3.3\cdot10^{-5}$: $\varrho_1=7.5299683$ and $\varrho_2=7.469616$. The numerical true rate $\approx7.46958$ (Section 7.1) agrees with $\varrho_2^{\rm eq}$ and not with $\varrho_1^{\rm eq}$. So, numerically, Theorem 2 is sharp, and the sign cancellation captured by the contour is real: it is worth $0.060$ in the rate. The lower bound of Proposition 8 is proved, but the values of $\varrho_1^{\rm eq}$ and $\varrho_2^{\rm eq}$ are floating-point computations.

### 7.3 Consequence for the exponent

Within this construction and this argument, the norm rate cannot be pushed below about $7.46958$: Proposition 8 limits the method, and the numerics limit the norm itself. By Remark 2, $\overline U\ge\overline U_{\rm eq}\approx-1.3855098$. The right side of (F3) equals $A_\infty+6\lambda/M-I_{M,\infty}$, where $A_\infty=I_{\rm out}+I_{3,\infty}\approx1.3134540$; it exceeds $A_\infty$ (by about $7\lambda/M$, since $I_{M,\infty}\approx-\lambda/M$). Hence $\varepsilon\le-(A_\infty+\overline U_{\rm eq})\approx0.0720558$, and

$$
\frac{\varrho\lambda}\varepsilon\ \gtrsim\ \frac{7.46958\cdot0.925}{0.0720558}\approx95.89 .
$$

Our $95.908$ is within $0.02$ of this limit. Going further requires a different construction, or a sharper arithmetic estimate than (F3), not better analysis of $\lVert B_K\rVert$ or $\overline U$. Nothing here bears on the true value of $\mu(\zeta(5))$, which is presumably $2$.

### 7.4 Where Fauzan's constant is lost

Fauzan's (C.2)–(C.6) can be read as Proposition 5 with the comparison measure $\lambda\,\omega[\frac1{18},\frac13]$, and with the residue weights and the polynomial growth bounded separately by crude estimates. Table 5 removes the losses one at a time [floating point; `scripts/16-mu5/loss_decomposition.py`; not part of the proof and not independently re-verified].

**Table 5.** Decomposition of the improvement $13.866\to7.470$ of the norm rate.

| step | constant | gain |
|---|---:|---:|
| Fauzan's $\kappa(\frac1{18},\frac13)$ | $13.866$ | |
| (C.2): exact $\max_{[a,b]}V$ in place of $\mathcal F(a,b)$ | $12.266$ | $1.600$ |
| (C.4): exact $\max\varphi$ in place of $\lambda(1-\log\alpha\lambda)$ | $10.027$ | $2.239$ |
| couple the residue weight with the polynomial growth | $9.520$ | $0.507$ |
| best single arcsine interval | $9.074$ | $0.446$ |
| equilibrium measure ($\varrho_1^{\rm eq}$) | $7.530$ | $1.544$ |
| contour, sign cancellation ($\varrho_2^{\rm eq}$) | $7.470$ | $0.060$ |

---

## 8. Computational certification and reproducibility

### 8.1 Software and data

All computations use Python 3.12.13, python-flint 0.9.0 [flint] (Arb ball arithmetic [Jo]), mpmath 1.3.0 [mpmath], and, for the exact roots in Section 7.1, PARI/GP [PARI] through cypari2. All inputs to certificates are exact rationals. Every comparison is "ball-certain": it is accepted only if the entire ball lies on the required side. Paths are relative to the root of the accompanying repository. The two data inputs to the proof are

- `data/16-mu5/rho1600.exact.json`: the test measure, SHA-256 `0529b8fda6f1ccb5b31e32f5638651ae7228544d7fe7da8bc7638f8f8474068f`;
- `data/15-redteam-mu5/int_20_4000.txt`: the exact $I_{20,4000}$, SHA-256 `84c23d2736d19deef5ef8e69b51ae586532fbe66018d0c2caeb33755cce5f956`.

### 8.2 The main script

Table 6 lists the certified quantities and the scripts. The main script `scripts/16-mu5/mu_final16.py` uses the helpers `certlib.py`, `certify.py` and `contour_certify.py` in the same directory. In one run it:

- checks (P1)–(P3);
- certifies $\sup_{[0,2]}(2U^\rho-V)<M_0$ and the inequality at $t=2$;
- computes $I(\rho)$, $C_*$ and $\overline U$;
- certifies (a)–(c) of Proposition 7;
- computes $A_\sharp$ and $A_M$ in exact rationals, with the cross-check $A_{10^5}=1.3135191719802$;
- computes $\varepsilon_{\max}$ and $\mu_0$, and checks the final inequalities in exact rationals.

For this manuscript it was re-run on one core (22 October 2026, `.venv/bin/python scripts/16-mu5/mu_final16.py`). It finished in $139$ s, of which $132$ s is the supremum certificate, and its output was identical to the archived log `data/16-mu5/mu_final16.log`. Its final lines are:

```
eps_max = -(A_M + U) = [0.072041643193987937085 +/- 3.16e-23]
[Theorem 2] threshold rho_K*lam/eps_max = [95.9083454190 +/- 1.87e-11];
  eps = 9005205399/125000000000, c = 513596269/40000000:
  eps*c - lam = 7.979e-06 > 0, mu - rho_K*c = 8.273e-04 > 0
  ==> mu(zeta(5)) <= 95.91
```

(In this log, "Theorem 2" is the contour bound, Proposition 6. The displayed radius $1.87\cdot10^{-11}$ is a printing artefact of the 12-digit output; the ball computed at higher precision is $95.90834541898132338\ldots\pm2\cdot10^{-54}$.)

**Table 6.** Certified quantities and the scripts that certify them. The labels 09, 15, 16, 18 refer to the directories `scripts/09-zeta7det`, `scripts/15-redteam-mu5`, `scripts/16-mu5` and `scripts/18-redteam-mu96`.

| statement | certified | primary (16) | independent re-check (18) |
|---|---|---|---|
| $\sup_{[0,2]}(2U^\rho-V)<M_0$ | $M_0=-6.6607665$ | `certify_sup`, 7733 pieces | `run_real_certs.py`, 14156 pieces |
| $I(\rho)$, $C_*$, $\overline U$ | balls, Prop. 3 | `certlib.py` | identical digits |
| $\inf_S(2U^\rho-V)>m_0$ | $m_0=-6.660814$ | `certify_inf`, 3412 pieces | `run_real_certs.py`, 12585 pieces |
| $\max_\Gamma H_\rho<M_\Gamma$ | $M_\Gamma=0.808802$ | `contour_certify.py`, 215 pieces | `run_contour_cert.py`, 73 pieces |
| $\sup(2U^\rho(-s^2)+\varphi)<M_1$ | $M_1=0.8691543$ | `certify_neg`, 78 pieces | not re-checked (not needed for Thm. 1) |
| $I_{20,4000}$ | exact rational | `exact_R.py` (15), `exact_A.py` (09) | `int_R_closed_form.py`, exactly equal |
| $A_M$, $\varepsilon_{\max}$, $\mu_0$, final inequalities | Cor. 1, §6 | `mu_final16.py` | `final_arith.py` |

### 8.3 Negative controls

A certificate that always succeeds proves nothing. We therefore checked that each certificate *fails* when its target is moved past the floating-point extremum.

- `scripts/16-mu5/negative_tests.py`: with an $800$-component measure of the same type, the supremum, infimum and segment certificates fail at targets $10^{-7}$ beyond the float extremum. With $\rho_{1600}$, the contour certificate fails at target $0.8088018<H_\rho(\tau_*)$.
- `scripts/18-redteam-mu96/negative_controls.py`, with $\rho_{1600}$. The supremum certificate fails at targets $-6.6607669$ and $-6.66076682$ (float maximum $-6.66076681785$) and succeeds at $M_0$ on the same subinterval. The infimum certificate fails at $-6.6608135$ and $-6.66081359$ (float minimum $-6.66081360140$). The contour certificate fails at $0.8088018$ and $0.80880185$.

### 8.4 Independent verifications

The work was checked twice. In each case a separate team member (an AI agent; see the Acknowledgements) re-derived the arguments from [F] and wrote their own code, reading the original code only for data formats.

**For the intermediate bound $181$** (`scripts/15-redteam-mu5`):

- exact reproduction of (F6), of the $17$ rows of [F, Table 3] and of $I_{\rm out}$;
- an independent exact $I_{20,4000}$;
- the $40$-periodicity of $Q$ (Remark 3);
- a re-derivation of the tail bound (13);
- an independent certificate for the $121$-component measure of the intermediate version;
- a calibration: on Fauzan's own $16$-component measure, the code reproduces [F, (A.10)] and [F, (6.4)] exactly ($I=-2.12659344514705$, $\overline U\le-1.3669955645$).

The `scripts/16-mu5` code in turn reproduces the $121$-component value $\overline U\le-1.3843989249340$ to all $13$ printed digits.

**For the bound of this paper** (`scripts/18-redteam-mu96`):

- admissibility of $\rho_{1600}$ in exact arithmetic;
- certificates for $M_0$ and $m_0$ by mean-value bounds with ball enclosures of $f'$ ($14156$ and $12585$ pieces), and for $M_\Gamma$ by ball enclosures of $\mathcal F'$ on boxes with explicitly chosen branches ($73$ pieces);
- $H_\rho(\tau_*)=0.80880185128432719515721067\pm3\cdot10^{-36}$;
- $I(\rho)$, $C_*$ and $\overline U$ to all printed digits;
- $I_{20,4000}$ from the closed form [F, (5.12)–(5.13), (B.1)], exactly equal to the other computation;
- $A_M$, $\varepsilon_{\max}$ and $\mu_0$;
- at $K=40$, an exact check of Lemma 1: the integral representation to all $25$ printed digits, the residue identity as an exact rational equality, and the contour identity on $\Gamma$ to relative error $8\cdot10^{-51}$;
- at $K=40$ and $80$, a test of every inequality in the proof of Theorem 2 with the exact extremal polynomial (obtained from an exact rational solve at the root of $\Delta_K$ nearest $\zeta(5)$): Lemma 3 at $40$ points of $S$, Corollary 2 and Lemma 5 at $83$ points of $\Gamma$, and the integrand bound of Proposition 6, all holding with margins of at least $7$ in the logarithm;
- an exact check of the closed form and signs of $r_j$ (Lemma 4) for all $j$ at $K=40,80$;
- the values of $\lVert B_K\rVert$ in Table 4, by three methods.

The finite-$K$ tests of Lemma 5 used the cruder constant $7/d$. The constant $1+6\alpha^2$ stated here comes from the proof above. Finally, `scripts/16-mu5/sanity_checks.py` checks Lemma 4 against the defining products, and the weight bound of Lemma 3(i), at $K=40,80,400$.

### 8.5 What is and is not machine-checked

- *Formally verified (Lean, by others):* the statements of [F] marked in Table 2, with Fauzan's own constants, modulo the prime number theorem as an axiom. This does not include Appendix C of [F].
- *Certified numerical computation (this work):* the inequalities of Propositions 3 and 7, the exact rational $I_{20,4000}$, and the final arithmetic. These are rigorous, provided the Arb library and our scripts are correct. Two independent implementations reduce, but do not eliminate, this dependence. They are not formal proofs.
- *Paper-level proofs (this work):* Propositions 1, 2 and 4, Lemmas 1–5, Propositions 5, 6 and 8, and the assembly. Each was re-derived by an independent reader, but none is formalized.
- *Floating-point only (not used in the proof):* the construction of $\rho_{1600}$ and of the contour, the equilibrium quantities, $\lVert B_K\rVert$ for $K>80$, the fits, Table 5, and the heuristics of Section 9.

---

## 9. Remarks and further directions

**Retuning the construction.** Section 7 shows that, with $\alpha=3/40$, the analysis is essentially exhausted. A heuristic model varies $\alpha=N/K$ with $\lambda=1-\alpha$ and zero multiplicity $6$. In it, $\varepsilon(\alpha)$ comes from a direct generalization of Fauzan's local estimates combined with the equilibrium value of $\overline U$, and $\varrho_2^{\rm eq}(\alpha)$ from the equilibrium measure. The model gives approximately the following [HEURISTIC]:

| $\alpha$ | $0.07$ | $0.075$ | $0.08$ | $0.0875$ | $0.1$ |
|---|---|---|---|---|---|
| $\mu(\alpha)\approx$ | $100.2$ | $95.6$ | $93.2$ | $92.2$ | $99.4$ |

The minimum is near $\alpha=7/80$ ($K=80n$, $N=7n$). At $\alpha=3/40$ the model's $\varepsilon$ is slightly optimistic ($0.0723$ against our certified $0.0720$). Realizing this gain would require re-proving the local estimates of [F, §§3–5] for the new parameters, a substantial task that we have not attempted.

**Formalization.** The natural next step is to formalize the new parts in Lean, building on [DR]:

1. Proposition 2 for general nested arcsine measures. The development of [DR] is written for Fauzan's measure, with a different regularization.
2. Lemmas 1–5 and Proposition 6. These need the residue theorem for a rectangle, the maximum principle for subharmonic functions, and elementary estimates.
3. The certificates. Following [DR], these could be done by interval arithmetic checked by `decide +kernel`. Our partitions have a few thousand cells, but the contour certificate needs complex logarithms and square roots.
4. Appendix C.2 of [F], i.e. Proposition 1.

This would make Theorem 1 fully machine-checked, modulo the prime number theorem.

**Effectivity.** Corollary 1 rests on the prime number theorem in limit form, so $b_0$ is ineffective. Even with an effective prime number theorem, the condition $K\ge200M^2$ with $M=4\cdot10^8$ forces $\log b_0\gtrsim2.5\cdot10^{18}$. A small cutoff costs only $O(1/M)$ in $\varepsilon$. For example, $M=10^4$ would lower $\varepsilon$ by about $6.5\cdot10^{-4}$, which costs about $0.9$ in the exponent.

**Other constructions.** The proof of Theorem 2 uses only three features: a positive weight, simple poles on a segment, and a test measure for the logarithmic energy. It should apply to other Hankel-determinant constructions of this type, such as the variant with $(K,N,h)=(12n,n,11n)$ discussed by Calegari [Ca]. We expect it to give essentially sharp norm rates there as well [HEURISTIC]; we have not checked this.

---

## Acknowledgements and use of AI

The research reported here was carried out by Claude (Anthropic) agents with separate roles. This includes the mathematical arguments of Sections 3–6, the test measure, all certified computations, the independent re-verifications, and the drafting of this manuscript.

- One agent obtained the intermediate bound $181$.
- One agent developed the norm bound and the certificates for the bound of this paper.
- Two further agents each re-derived one of these two results independently, working from [F] with their own code.
- One agent wrote this text from the research notes.

**Authorship and responsibility.** The mathematics, computations and text are the work of Claude Opus 5.5 agents. Ryan Xu set up the project, chose what to release, and is the contact for questions and corrections (ryan@ryanxu.net). No human expert has verified the arguments.

We relied on the Lean developments [DR, Mo] and on the audit [Au] to locate the formalized statements of [F].

---

## References

- **[An]** P. Anand, *Spectral theory and Diophantine bounds for Hankel determinants of odd zeta values*, Zenodo, doi:10.5281/zenodo.23528910, 21 October 2026.
- **[Ap]** R. Apéry, *Irrationalité de $\zeta2$ et $\zeta3$*, Astérisque **61** (1979), 11–13.
- **[Au]** Claude Opus 5 (Anthropic), prepared for D. Romik, *Referee audit of "ζ(5) is irrational" (A. Fauzan, 17 September 2026)*, 22 September 2026; file `docs/zeta5-audit.tex` in [DR].
- **[BR]** V. S. Buyarov and E. A. Rakhmanov, *Families of equilibrium measures in an external field on the real axis*, Mat. Sb. **190** (1999), no. 6, 11–22; English transl. Sb. Math. **190** (1999), 791–802.
- **[Ca]** F. Calegari, *zeta(5) is irrational*, blog post, Persiflage (galoisrepresentations.org), 24 September 2026.
- **[DLMF]** NIST Digital Library of Mathematical Functions, https://dlmf.nist.gov/, §18.14.
- **[DR]** D. Romik (director), *zeta5-irrationality: a Lean 4 formalization of A. Fauzan's proof*, https://github.com/danromik/zeta5-irrationality, commit `fc2fefa` (25 September 2026). The Lean code was written by Claude (Anthropic).
- **[F]** A. Fauzan, *ζ(5) is irrational*, preprint, Zenodo, doi:10.5281/zenodo.22826419, version 1, 17 September 2026.
- **[flint]** python-flint, version 0.9.0, Python bindings for FLINT and Arb, https://github.com/flintlib/python-flint.
- **[Jo]** F. Johansson, *Arb: efficient arbitrary-precision midpoint-radius interval arithmetic*, IEEE Trans. Comput. **66** (2017), 1281–1292.
- **[Mo]** mo271, *Zeta5: a Lean 4 formalization of the irrationality of ζ(5)*, https://github.com/mo271/Zeta5, commit `7fe7367` (27 September 2026).
- **[mpmath]** The mpmath development team, *mpmath: a Python library for arbitrary-precision floating-point arithmetic*, version 1.3.0, https://mpmath.org/.
- **[PARI]** The PARI Group, *PARI/GP*, Univ. Bordeaux, https://pari.math.u-bordeaux.fr/.
- **[Ra]** T. Ransford, *Potential Theory in the Complex Plane*, London Math. Soc. Student Texts 28, Cambridge Univ. Press, 1995.
- **[Ro]** H. Robbins, *A remark on Stirling's formula*, Amer. Math. Monthly **62** (1955), 26–29.
- **[RV]** G. Rhin and C. Viola, *The group structure for ζ(3)*, Acta Arith. **97** (2001), 269–293.
- **[Sa]** E. B. Saff, *Logarithmic potential theory with applications to approximation theory*, Surv. Approx. Theory **5** (2010), 165–200, arXiv:1010.3760.
- **[ST]** E. B. Saff and V. Totik, *Logarithmic Potentials with External Fields*, Grundlehren der mathematischen Wissenschaften 316, Springer, 1997.
