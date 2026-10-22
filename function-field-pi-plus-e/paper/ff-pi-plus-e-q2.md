# On the Carlitz analogue of $\pi+e$ and $\pi e$ in characteristic two

**Claude Opus 5.5 (Anthropic)**

*Released by Ryan Xu (ryan@ryanxu.net), who directed the project but is not a specialist and has not independently verified the mathematics.*

*Preprint, October 2026*

*MSC 2020:* 11J72 (primary); 11J93, 11J70, 11J61.
*Keywords:* Carlitz module, Carlitz period, Carlitz exponential, linear independence, continued fractions, function fields.

---

> **Abstract.** Let $\tilde\pi$ be the Carlitz period and $e_C(1)=\sum_{i\ge 0}1/D_i$ the value at $1$ of the Carlitz exponential for $\mathbb F_2[\theta]$; both lie in $\mathbb F_2((1/\theta))$. We prove that $1,\ \tilde\pi,\ e_C(1),\ \tilde\pi\,e_C(1)$ are linearly independent over $\mathbb F_2(\theta)$. In particular $\tilde\pi+e_C(1)$, $\tilde\pi\,e_C(1)$ and $\tilde\pi/e_C(1)$ do not lie in $\mathbb F_2(\theta)$, the exact analogue of the open problems for $\pi+e$, $\pi e$ and $\pi/e$. The same holds with $e_C(1)$ replaced by $e_C(1/f)$ for any nonzero $f\in\mathbb F_2[\theta]$. It also holds over every finite extension of $\mathbb F_2(\theta)$ that is linearly disjoint from $\mathbb F_2((1/\theta))$.
>
> For $q\ge 3$, Denis proved in 1995 that $\tilde\pi$ and $e_C(1)$ are algebraically independent over $\mathbb F_q(\theta)$. His argument uses the substitution $\theta\mapsto\xi\theta$ with $\xi\in\mathbb F_q^\times\smallsetminus\{1\}$, which does not exist when $q=2$.
>
> Our proof is an elementary Liouville-type argument. It combines a new explicit family of approximations to $\tilde\pi$ of Dirichlet quality with the partial sums of $e_C(1)$, which are the convergents of Thakur's continued fraction that precede its large partial quotients. Balancing the two families forces the two error terms to have the same size. The exact second-order term of the approximations to $\tilde\pi$ breaks this tie.
>
> We isolate the single ingredient that is missing in characteristic zero. We also show that $\tilde\pi$ is fixed by $\theta\mapsto\theta+1$. Consequently all partial quotients of $\tilde\pi$ have even degree, and those of $\tilde\pi/\theta^2$ have degree $1$ or even degree.

---

## 1. Introduction

### 1.1. The classical problem

It is a famous open problem whether $\pi+e$, $\pi e$ or $\pi/e$ is irrational. All three would follow from the linear independence of $1,\pi,e,\pi e$ over $\mathbb Q$. That statement is equivalent to saying that $\pi$ is not of the form $(ae+b)/(ce+d)$ with $a,b,c,d\in\mathbb Q$.

Because $\pi$ and $e$ are transcendental, at least one of $\pi+e$ and $\pi e$ is transcendental: both are roots of $X^2-(\pi+e)X+\pi e$. Beyond this, essentially nothing is known. The algebraic independence of $\pi$ and $e$, which would settle all these questions, is a special case of Schanuel's conjecture.

### 1.2. The Carlitz setting

Let $q$ be a power of a prime. Write $k=\mathbb F_q(\theta)$, $k_\infty=\mathbb F_q((1/\theta))$, and let $\mathbb C_\infty$ be the completion of an algebraic closure of $k_\infty$. The well-known analogy is
$$\mathbb Z\leftrightarrow\mathbb F_q[\theta],\quad \mathbb Q\leftrightarrow k,\quad \mathbb R\leftrightarrow k_\infty,\quad \mathbb C\leftrightarrow\mathbb C_\infty .$$
It is completed by Carlitz's exponential [Car35]
$$e_C(z)=\sum_{i\ge 0}\frac{z^{q^i}}{D_i},\qquad D_i=\prod_{j=0}^{i-1}\bigl(\theta^{q^i}-\theta^{q^j}\bigr).$$
This is an entire $\mathbb F_q$-linear function satisfying $e_C(\theta z)=\theta\,e_C(z)+e_C(z)^q$.

- Its kernel is $\mathbb F_q[\theta]\,\tilde\pi$, where the *Carlitz period* $\tilde\pi\in\mathbb C_\infty^\times$ is unique up to $\mathbb F_q^\times$.
- $D_i$ is the product of the monic polynomials of degree $i$, and plays the role of $n!$.
- $\tilde\pi$ is the analogue of $2\pi i$, and $e_C(1)=\sum_i 1/D_i$ is the analogue of $e$ [Tha92, §0], [Tha04].

Both $\tilde\pi$ and $e_C(1)$ are transcendental over $k$. For $\tilde\pi$ this is due to Wade [Wad41]. The Carlitz analogue of the Hermite–Lindemann theorem, which gives the transcendence of $e_C(1)$, also goes back to Wade (see [Pel17, §2]).

The function-field analogue of the problem in §1.1 asks whether $1,\tilde\pi,e_C(1),\tilde\pi\,e_C(1)$ are linearly independent over $k$. For $q\ge3$ the right coefficient field is $k\bigl((-\theta)^{1/(q-1)}\bigr)$; see §1.5.

### 1.3. What is known

Positive characteristic offers several deep algebraic independence theorems.

- **Thiery [Thi92]** proved the analogue of the Lindemann–Weierstrass theorem. If $\alpha_1,\dots,\alpha_n\in\bar k$ are linearly independent over $k$, then $e_C(\alpha_1),\dots,e_C(\alpha_n)$ are algebraically independent over $k$ (see [Pel17, §2]).
- **Papanikolas [Pap08]**, using the linear independence criterion of Anderson, Brownawell and Papanikolas [ABP04], proved the following. Carlitz logarithms of algebraic points that are linearly independent over $k$ are algebraically independent over $\bar k$; up to a factor in $k$, $\tilde\pi$ is one of them.

The value $e_C(1)$, however, is an exponential value at an algebraic point, not a logarithm of one. So the "mixed" question about $\tilde\pi$ and $e_C(1)$, the true analogue of $\pi$ versus $e$, is not covered by these results. For $q\ge 3$ it was settled by Denis:

> **Theorem (Denis [Den95, Cor. 2(a)]).** *If $q\ge 3$, then $\tilde\pi$ and $e_C(1)$ are algebraically independent over $k$.*

The proof rests on a symmetry. For $\xi\in\mathbb F_q^\times$, the substitution $\theta\mapsto\xi\theta$ multiplies $D_i$ by $\xi^i$. Hence it maps $e_C(1)=\sum_iD_i^{-1}$ to
$$\sum_i\xi^{-i}D_i^{-1}=e_C(\mu)/\mu,\qquad \mu^{q-1}=\xi^{-1},$$
while it fixes a normalization of the period that differs from $\tilde\pi$ by an algebraic factor.

When $\xi\neq1$ we have $\mu\notin\mathbb F_q$, so $1$ and $\mu$ are $k$-linearly independent. An algebraic relation between $\tilde\pi$ and $e_C(1)$ would then produce one between $e_C(1)$ and $e_C(\mu)$, contradicting Thiery's theorem.

For $q=2$ there is no $\xi\ne 1$. In fact, the only nontrivial automorphism of $\mathbb F_2[\theta]$, namely $\theta\mapsto\theta+1$, fixes both $\tilde\pi$ and $e_C(1)$ (Proposition 8.1), so no argument of this kind is available. Accordingly, Pellarin's discussion [Pel17, §2] states Denis's result under the hypothesis $q\ge3$.

Denis also formulated a Carlitz analogue of Schanuel's conjecture ([Pel17, Conj. 2]). It would imply the algebraic independence of $\tilde\pi$ and $e_C(1)$ for every $q$, including $q=2$. We are not aware of any unconditional result for $q=2$ relating $\tilde\pi$ and $e_C(1)$ beyond the transcendence of each.

### 1.4. Results

From now on $q=2$. Then
$$\tilde\pi=\theta^2\prod_{m\ge1}\bigl(1+\theta^{1-2^m}\bigr)^{-1}$$
lies in $k_\infty=\mathbb F_2((1/\theta))$ (see §2), as does $e_C(1)$.

> **Theorem 1.1.** *Let $q=2$. The four elements $1,\ \tilde\pi,\ e_C(1),\ \tilde\pi\,e_C(1)$ of $\mathbb F_2((1/\theta))$ are linearly independent over $\mathbb F_2(\theta)$.*

> **Corollary 1.2.** *Let $q=2$. None of $\tilde\pi+e_C(1)$, $\tilde\pi\,e_C(1)$, $\tilde\pi/e_C(1)$ lies in $\mathbb F_2(\theta)$. More generally, $\tilde\pi\neq\dfrac{a\,e_C(1)+b}{c\,e_C(1)+d}$ for all $a,b,c,d\in\mathbb F_2(\theta)$ with $c\,e_C(1)+d\neq0$.*

> **Theorem 1.3.** *Let $q=2$.*
> 1. *For every nonzero $f\in\mathbb F_2[\theta]$, the elements $1,\ \tilde\pi,\ e_C(1/f),\ \tilde\pi\,e_C(1/f)$ are linearly independent over $\mathbb F_2(\theta)$.*
> 2. *Theorem 1.1 and part 1 remain true over any finite extension $L\subset\mathbb C_\infty$ of $\mathbb F_2(\theta)$ that is linearly disjoint from $\mathbb F_2((1/\theta))$ over $\mathbb F_2(\theta)$. This holds, for instance, if $L/\mathbb F_2(\theta)$ is separable with a single place above $\infty$. In particular, both statements hold over $\overline{\mathbb F}_2(\theta)$.*

The elements $e_C(1/f)$ are analogues of $e^{1/n}$. Thakur studies their Möbius images as analogues of Hurwitz numbers [Tha11, §3].

### 1.5. Why $q=2$ is the most classical case

Four features make $q=2$ the closest analogue of the classical problem.

1. **$\tilde\pi$ is "real".** For $q\ge3$ we have $\tilde\pi\in\lambda k_\infty$ with $\lambda=(-\theta)^{1/(q-1)}\notin k_\infty$, just as $2\pi i\in i\mathbb R$. The analogue of "$1,\pi,e,\pi e$ over $\mathbb Q$" is then "$1,\tilde\pi,e_C(1),\tilde\pi e_C(1)$ over $k(\lambda)$". For $q=2$, $\lambda=\theta$ and $\tilde\pi\in k_\infty$, just as $\pi\in\mathbb R$.
2. **No extra Diophantine room.** For $q\ge3$, Frobenius makes both constants unusually well approximable.
   - $e_C(1)$ has approximation exponent $q$ [Tha11, Thm. 1].
   - $u=\tilde\pi/(\theta\lambda)\in k_\infty$ has exponent $\ge q-1$, with equality for $q\ge5$ [Tha11, Thm. 7].

   This Liouville-like behaviour has no classical counterpart, and it makes an elementary proof easy for $q\ge 3$ (§8.2). For $q=2$, the exponent of $e_C(1)$ is $2$, as for $e$, and numerically the exponent of $\tilde\pi$ is also $2$. Our argument works exactly at the threshold (§8.1).
3. **Hurwitz-type continued fraction.** For $q=2$, Thakur [Tha92] found
   $$e_C(1)=[1;[1],[2],[1],[3],[1],[2],[1],[4],\dots],$$
   whose $n$-th partial quotient is $[v_2(n)+1]$, where $[j]=\theta^{2^j}+\theta$. Compare Euler's $e=[2;1,2,1,1,4,1,1,6,\dots]$. In both, the partial quotients grow along an explicit subsequence, and the preceding convergents are explicit; for $e_C(1)$ they are the partial sums of $\sum1/D_i$ (Remark 4.2). These are exactly the approximations our proof uses.
4. **No symmetry.** Classically, the only continuous automorphisms of $\mathbb C$ fix both $\pi$ and $e$. For $q=2$, the only nontrivial automorphism of $\mathbb F_2[\theta]$ fixes both $\tilde\pi$ and $e_C(1)$ (Proposition 8.1). Denis's mechanism is absent in both settings.

Theorem 1.1 thus shows that, in the most faithful function-field model, the $\pi+e$ problem can be solved by an explicit Diophantine argument. In §8.3 we isolate the one input that has no known classical counterpart.

### 1.6. Idea of the proof

Write $\tilde\pi=\theta^2u$, and suppose $aue+bu+ce+d=0$ with $a,b,c,d\in\mathbb F_2[\theta]$, where $e=e_C(1)$. We form
$$Z_k=Q_MD_N\,(aV_MS_N+bV_M+cS_N+d)\in\mathbb F_2[\theta].$$
Here $S_N=\sum_{i\le N}1/D_i$, and $V_M=P_M/Q_M$ is a new explicit approximation to $u$ of Dirichlet quality (Lemma 3.1). We show that $\deg Z_k<0$ while $Z_k\ne0$, which is absurd. The argument has four steps.

1. **The window.** The parameters must lie in a narrow window, $N=2^k-1$ and $M=N+k-1$ (Lemma 5.1). This window is forced, and inside it the error terms from $u-V_M$ and from $e-S_N$ have the same size.
2. **The leading term.** The tie is broken by the exact second-order term
   $$u-V_M=x^{2A}(1+x)u+O(x^{4A+1}),\qquad x=1/\theta.$$
   The leading coefficient of $Z_k$ is then an element $\gamma\in k_\infty$ independent of $k$, so $\gamma=0$ for every relation.
3. **A derived relation.** Substituting the relation into $\gamma=0$ yields a new relation of a special shape.
4. **A second round.** Applying the same argument to the derived relation finishes the proof.

**Organization.** §2 fixes normalizations. §§3–5 contain the three lemmas: approximations to $u$, the tail of $e_C(1)$, and the window. §6 proves Theorem 1.1, and §7 proves the corollaries. §8 contains remarks on four topics:
- the limits of the method;
- the case $q\ge3$ and the situation in characteristic zero;
- a symmetry, with its consequences for continued fractions;
- open problems.

§9 describes the numerical checks.

---

## 2. Setup and normalizations

Throughout, $q=2$ and $x=1/\theta$, so $k=\mathbb F_2(\theta)$ and $k_\infty=\mathbb F_2((x))$.

**Degree.** For $0\ne\xi=\sum_{n\ge n_0}c_nx^n$ with $c_{n_0}\ne0$, set $\deg\xi=-n_0$ and $|\xi|=2^{\deg\xi}$; also $\deg0=-\infty$. This extends the degree of polynomials. We have $\deg(\xi\eta)=\deg\xi+\deg\eta$, and $\deg(\xi+\eta)\le\max(\deg\xi,\deg\eta)$, with equality if the degrees differ. The *Liouville principle* is
$$\text{a nonzero element of }\mathbb F_2[\theta]\text{ has degree}\ \ge 0. \tag{2.1}$$
We write $O(x^m)$, for $m\in\mathbb Z$, for an unspecified element of $x^m\,\mathbb F_2[[x]]$. In characteristic $2$ signs play no role.

**Carlitz factorials.** Let $[i]=\theta^{2^i}+\theta$, $D_0=1$, and $D_i=\prod_{j=0}^{i-1}(\theta^{2^i}+\theta^{2^j})$ for $i\ge1$.

> **Lemma 2.1.**
> 1. *$D_i=[i]\,D_{i-1}^2$ for $i\ge1$. Hence $D_i\mid D_N$ for $i\le N$, and $\deg D_i=i\,2^i$.*
> 2. *$1/D_i=x^{i2^i}\prod_{j<i}\bigl(1+x^{2^i-2^j}\bigr)^{-1}$. In particular $1/D_i\in x^{i2^i}\bigl(1+x^{2^{i-1}}\mathbb F_2[[x]]\bigr)$ for $i\ge1$.*

*Proof.* (1) We have $D_{i-1}^2=\prod_{j=0}^{i-2}(\theta^{2^i}+\theta^{2^{j+1}})=\prod_{j=1}^{i-1}(\theta^{2^i}+\theta^{2^j})$. The remaining factor, $j=0$, is $[i]$. So $D_{i-1}\mid D_i$, and by induction $\deg D_i=2^i+2(i-1)2^{i-1}=i2^i$. (2) Factor $\theta^{2^i}$ out of each of the $i$ factors. The smallest exponent among $2^i-2^j$ with $j<i$ is $2^{i-1}$. $\square$

**The constants.** Put
$$e:=e_C(1)=\sum_{i\ge0}\frac1{D_i}\ \in\ 1+x^2\,\mathbb F_2[[x]],\qquad u:=\prod_{m\ge1}\bigl(1+x^{2^m-1}\bigr)^{-1}\ \in\ 1+x\,\mathbb F_2[[x]].$$
From now on $e$ denotes $e_C(1)$. The exceptions are §8.3 and §8.5, where $e$ is Euler's number and the Carlitz value is written $e_C(1)$ in full.

The classical product formula for the Carlitz period [Car35], [Gos96], [Pel17, (2)] is
$$\tilde\pi=-(-\theta)^{q/(q-1)}\prod_{m\ge1}(1-\theta^{1-q^m})^{-1}.$$
For $q=2$ it reads
$$\tilde\pi=\theta^2u\in k_\infty,\qquad \deg\tilde\pi=2 . \tag{2.2}$$
For $q=2$ there is no choice of root or sign, and $\tilde\pi$ is uniquely determined because $\mathbb F_2^\times=\{1\}$. We also checked numerically, to precision $x^{20000}$ (§9), that $e_C(\theta^2u)=0$ and that $e_C(\theta u)=\theta$, the nonzero $\theta$-torsion point.

Since $\theta^2\in k^\times$, Theorem 1.1 is equivalent to:
$$\text{if } a,b,c,d\in k \text{ satisfy } aue+bu+ce+d=0,\ \text{ then } a=b=c=d=0. \tag{2.3}$$

We also use the partial products and the polynomials
$$U_n:=\prod_{m=1}^n\bigl(1+x^{2^m-1}\bigr)^{-1},\qquad E_M(x):=\prod_{m=1}^M\bigl(1+x^{2^m-1}\bigr)\in\mathbb F_2[x],$$
so that $U_M=1/E_M$.

> **Remark 2.2 ($\tilde\pi$ as a logarithm).** For $q=2$ we have $C_{[1]}(1)=0$, that is, $1$ is a $[1]$-torsion point of the Carlitz module. Indeed $C_\theta(1)=\theta+1$ and $C_\theta(\theta+1)=\theta+1$. Moreover
> $$\tilde\pi=[1]\,\log_C(1),\qquad \log_C(z)=\sum_{i\ge0}\frac{z^{2^i}}{L_i},\quad L_i=[i][i-1]\cdots[1],$$
> by [Tha92, Remark (3)]; we verified $\log_C(1)=u/(1+x)$ numerically. This is the analogue of $2\pi i=2\log(-1)$. Theorem 1.1 may thus be read as the $k$-linear independence of $1$, $\log_C(1)$, $\exp_C(1)$ and $\log_C(1)\exp_C(1)$.

---

## 3. Lemma A: Dirichlet-quality approximations to $\tilde\pi$

> **Lemma 3.1.** *Let $M\ge1$ and $A:=2^{M+1}-1$, and put*
> $$V_M:=\frac{1+x^{A}}{E_M(x)},\qquad P_M:=\theta^A+1,\qquad Q_M:=\theta^A\,E_M(1/\theta).$$
> 1. *$\deg_xE_M=A-M-1$. Consequently $Q_M=\theta^{M+1}\prod_{m=1}^{M}(\theta^{2^m-1}+1)\in\mathbb F_2[\theta]$ has degree $A$, and $Q_MV_M=P_M$.*
> 2. *$V_M=(1+x^{2A})(1+x^{2A+1})\,U_{M+2}$.*
> 3. *$u-V_M=x^{2A}(1+x)\,u+x^{4A+1}u+O(x^{4A+3})$. In particular $u-V_M=x^{2A}(1+x)\,u+O(x^{4A+1})$.*
> 4. *$|Q_Mu-P_M|=2^{-A}=|Q_M|^{-1}$.*

*Proof.*

1. We have $\sum_{m=1}^M(2^m-1)=2^{M+1}-2-M=A-M-1$. Multiplying $E_M(1/\theta)$ by $\theta^A=\theta^{M+1}\theta^{A-M-1}$ gives the stated product, which is monic of degree $A$. Also $Q_MV_M=\theta^A(1+\theta^{-A})=P_M$.
2. In characteristic $2$, $(1+x^A)^2=1+x^{2A}$, and $2^{M+2}-1=2A+1$. Hence
   $$V_M=\frac{(1+x^A)^2(1+x^{2A+1})}{E_M(x)\,(1+x^{A})(1+x^{2A+1})}=(1+x^{2A})(1+x^{2A+1})\,U_{M+2}.$$
3. Write $u=U_{M+2}R$ with $R=\prod_{m\ge M+3}(1+x^{2^m-1})^{-1}$. Since $2^{M+3}-1=4A+3$, we have $R\in1+x^{4A+3}\mathbb F_2[[x]]$, so $U_{M+2}=u+O(x^{4A+3})$. By part 2,
   $$V_M=\bigl(1+x^{2A}+x^{2A+1}+x^{4A+1}\bigr)\,u+O(x^{4A+3}),$$
   which is part 3, since signs are irrelevant.
4. $Q_Mu-P_M=Q_M(u-V_M)$ has degree $A-2A=-A$ by part 3, because $\deg((1+x)u)=0$. $\square$

> **Remark 3.2.** The partial products $U_M=\theta^{A-M-1}/\prod_{m\le M}(\theta^{2^m-1}+1)$ are useless for $q=2$. Their denominators have degree $A-M-1$ while $u-U_M$ has degree $-A$, so $|\text{den}\cdot u-\text{num}|=2^{-M-1}$, far from Dirichlet quality.
>
> The approximation $V_M=(1+x^A)\,U_M$ replaces the next factor $(1+x^A)^{-1}$ of the product by its first-order truncation $1+x^A$. The denominator degree grows by only $M+1$, while the error drops from $x^{A}$ to $x^{2A}$; equivalently, $V_M=(1+x^{2A})U_{M+1}$.
>
> What the proof needs is not only the size of $u-V_M$ but the exact shape of its leading term, part 3 of Lemma 3.1.

---

## 4. Lemma B: the tail of $e_C(1)$

> **Lemma 4.1.** *For $N\ge0$ let $S_N:=\sum_{i=0}^N1/D_i$, $B:=\deg D_N=N2^N$ and $Y:=(N+1)2^{N+1}$. Then $D_NS_N\in\mathbb F_2[\theta]$ and*
> $$e-S_N=x^{Y}\bigl(1+x^{2^N}+O(x^{2^N+1})\bigr);\qquad\text{in particular}\quad e-S_N=x^{Y}\bigl(1+O(x^{2^N})\bigr).$$

*Proof.* First, $D_NS_N=\sum_{i\le N}D_N/D_i\in\mathbb F_2[\theta]$ by Lemma 2.1(1).

By Lemma 2.1(2), $1/D_{N+1}=x^{Y}\prod_{j\le N}(1+x^{2^{N+1}-2^j})^{-1}$. The exponent $2^{N+1}-2^j$ equals $2^N$ for $j=N$ and is at least $3\cdot2^{N-1}\ge2^N+1$ for $j<N$. Hence $1/D_{N+1}=x^Y(1+x^{2^N}+O(x^{2^N+1}))$.

For $i\ge N+2$, $1/D_i=O(x^{(N+2)2^{N+2}})$, and $(N+2)2^{N+2}-Y=(N+3)2^{N+1}>2^N+1$. $\square$

> **Remark 4.2 (Thakur's continued fraction).** By [Tha92, Thm. 1 and its proof], for $q=2$ the convergent of index $2^N-1$ of $e=[1;[1],[2],[1],[3],\dots]$ is $S_N$, with denominator $D_N$. The next partial quotient is $[N+1]$, of degree $2^{N+1}$. This is consistent with Lemma 4.1, since $Y=2B+2^{N+1}$.
>
> Thus $|D_Ne-D_NS_N|=|D_N|^{-1-2/N}$: the approximations $S_N$ beat Dirichlet's exponent by the factor $|[N+1]|=|D_N|^{2/N}$. These "spikes" play the role of the partial quotients $2m$ in Euler's continued fraction of $e$. We checked $q_{2^N-1}=D_N$, and all partial quotients, numerically up to denominators of degree $10^4$ (§9).

---

## 5. The window

The approximations of §§3–4 must be balanced: in the proof, the error $u-V_M$ is multiplied by $D_N$ and the error $e-S_N$ by $Q_M$.

> **Lemma 5.1 (window).** *Let $k\ge1$, $N:=2^k-1$ and $M:=N+k-1$. Let $A=2^{M+1}-1$, $B=N2^N$ and $Y=(N+1)2^{N+1}$ be as in Lemmas 3.1 and 4.1. Then*
> $$A=(N+1)2^N-1,\qquad 2A+2=Y,\qquad B-A=1-2^N,\qquad A+B-Y=-1-2^N .$$

*Proof.* $2^{M+1}=2^{N+k}=2^k\cdot 2^N=(N+1)2^N$. The rest is arithmetic. $\square$

> **Remark 5.2 (the window is forced, and the two terms tie).**
>
> *The window is forced.* The argument of §6 needs both
> $$\deg(Q_MD_N\,x^{2A})=B-A\qquad\text{and}\qquad \deg(Q_MD_N\,x^{Y})=A+B-Y$$
> to tend to $-\infty$. That is, $N2^N<A<(N+2)2^N$ with both margins unbounded. Here $A+1$ is a power of $2$ lying in $(N2^N,(N+2)2^N]$, so it is a multiple of $2^N$ and equals $(N+1)2^N$ or $(N+2)2^N$. The second candidate leaves the margin $Y-A-B=1$, so $A+1=(N+1)2^N$. Hence $N+1=2^k$ and $M=N+k-1$, so Lemma 5.1 describes the only possible choice.
>
> *The terms tie.* For this choice $Y=2A+2$. The two error terms $x^{2A}(1+x)u\,\alpha$ and $x^{Y}\beta$ of §6 therefore have the same order of magnitude, up to a factor independent of $k$, and neither can be shown to dominate. This tie is why the exact form of Lemma 3.1(3) is needed.

---

## 6. Proof of Theorem 1.1

Let $a,b,c,d\in k$ satisfy
$$a\,u\,e+b\,u+c\,e+d=0 . \tag{6.1}$$
Clearing denominators, we may assume $a,b,c,d\in\mathbb F_2[\theta]$. Put
$$\alpha:=ae+b,\qquad \beta:=au+c,\qquad \gamma:=(1+x)\,u\,\alpha+x^2\beta\ \in k_\infty .$$

### 6.1. The characteristic-2 identity

For $k\ge1$, let $N,M,A,B,Y$ be as in Lemma 5.1. Write $V=V_M$, $S=S_N$, $\delta_u:=u-V$, $\delta_e:=e-S$, and
$$Z_k:=Q_MD_N\,\bigl(aVS+bV+cS+d\bigr).$$

> **Lemma 6.1.** *$Z_k\in\mathbb F_2[\theta]$, and $Z_k=Q_MD_N\,\bigl(\delta_u\alpha+\delta_e\beta+a\,\delta_u\delta_e\bigr)$.*

*Proof.* Since $Q_MV=P_M$ and $D_NS\in\mathbb F_2[\theta]$, each of the four products below is a polynomial:
$$Q_MD_N\cdot aVS=a\,P_M\,(D_NS),\quad Q_MD_N\cdot bV,\quad Q_MD_N\cdot cS,\quad Q_MD_N\,d.$$
Substituting $V=u-\delta_u$ and $S=e-\delta_e$ gives
$$aVS+bV+cS+d=(aue+bu+ce+d)-\delta_u(ae+b)-\delta_e(au+c)+a\,\delta_u\delta_e .$$
The first bracket vanishes by (6.1), and in characteristic $2$ the signs disappear. $\square$

### 6.2. The leading term $\gamma$

> **Proposition 6.2.** *If $a,b,c,d\in\mathbb F_2[\theta]$ satisfy (6.1), then $\gamma=0$.*

*Proof.* Let $C:=\max(\deg a,\deg b,\deg c,\deg d)$. Since $\deg u=\deg e=0$, we have $\deg\alpha\le C$ and $\deg\beta\le C$. Lemma 3.1(3), Lemma 4.1 and $Y=2A+2$ give
$$\delta_u\alpha=x^{2A}(1+x)\,u\,\alpha+O(x^{4A+1-C}),\qquad \delta_e\beta=x^{2A+2}\beta+O(x^{2A+2+2^N-C}),\qquad a\,\delta_u\delta_e=O(x^{4A+2-C}).$$
Since $2A+1=(N+1)2^{N+1}-1\ge 2^N+2$, it follows that
$$\delta_u\alpha+\delta_e\beta+a\,\delta_u\delta_e=x^{2A}\,(\gamma+\varepsilon_k),\qquad \varepsilon_k=O\bigl(x^{2^N+2-C}\bigr).$$

Suppose $\gamma\ne0$. As $k\to\infty$ we have $N=2^k-1\to\infty$, so for large $k$ we get $\deg\varepsilon_k<\deg\gamma$. Hence $\gamma+\varepsilon_k\neq0$ and $\deg(\gamma+\varepsilon_k)=\deg\gamma$. By Lemmas 6.1 and 5.1, using $\deg Q_M=A$ and $\deg D_N=B$,
$$Z_k\neq 0,\qquad \deg Z_k=A+B-2A+\deg\gamma=1-2^N+\deg\gamma ,$$
which is negative for large $k$. This contradicts (2.1), so $\gamma=0$. $\square$

The element $\gamma$ does not depend on $k$. This is exactly where the precise leading terms of Lemmas 3.1 and 4.1 enter.

### 6.3. The derived relation

By (6.1), $u\alpha=aue+bu=ce+d$. Hence $\gamma=0$ reads
$$(1+x)(ce+d)+x^2(au+c)=0 .$$
Multiply by $\theta^2$, noting that $\theta^2(1+x)=\theta^2+\theta=[1]$:
$$a\,u+[1]\,c\,e+\bigl(c+[1]\,d\bigr)=0 . \tag{6.2}$$
This is again a relation of the form (6.1), now with coefficients
$$(a',b',c',d')=\bigl(0,\ a,\ [1]c,\ c+[1]d\bigr)\in\mathbb F_2[\theta]^4 .$$

### 6.4. The second round

Apply Proposition 6.2 to (6.2). There $\alpha'=a'e+b'=a$ and $\beta'=a'u+c'=[1]c$. Since $x^2[1]=1+x$,
$$0=\gamma'=(1+x)\,u\,a+x^2[1]\,c=(1+x)\,(au+c).$$
Hence $au+c=0$. Since $u\notin k$ (Lemma 6.3), $a=c=0$. Then (6.2) reduces to $[1]\,d=0$, so $d=0$, and (6.1) reduces to $bu=0$, so $b=0$. This proves (2.3), and hence Theorem 1.1. $\square$

> **Lemma 6.3.** *$u\notin k$.*

*Proof.* Of course $u$ is even transcendental [Wad41]; here is a self-contained argument. Suppose $u=r/s$ with $r,s\in\mathbb F_2[\theta]$ and $s\neq0$. Then $s\,(Q_Mu-P_M)=rQ_M-sP_M\in\mathbb F_2[\theta]$. It is nonzero by Lemma 3.1(4), and for $M$ large it has degree $\deg s-A<0$. This contradicts (2.1). $\square$

The proof never uses $e\notin k$; this follows from Theorem 1.1 itself (take $a=b=0$).

---

## 7. Corollaries

### 7.1. Möbius images (Corollary 1.2)

Suppose $\tilde\pi\,(ce+d)=ae+b$ with $a,b,c,d\in k$. Then $c\,\tilde\pi e+d\,\tilde\pi+ae+b=0$, and Theorem 1.1 gives $a=b=c=d=0$, which is impossible if $ce+d\ne0$.

The relations $\tilde\pi+e=r$, $\tilde\pi e=r$ and $\tilde\pi/e=r$ with $r\in k$ correspond to the matrices
$$\begin{pmatrix}1&r\\0&1\end{pmatrix},\qquad \begin{pmatrix}0&r\\1&0\end{pmatrix},\qquad \begin{pmatrix}r&0\\0&1\end{pmatrix}.$$
Symmetrically, $e$ is not a $k$-Möbius image of $\tilde\pi$. $\square$

### 7.2. Proof of Theorem 1.3(1): $e_C(1/f)$

Let $0\ne f\in\mathbb F_2[\theta]$, $F:=\deg f$, and
$$e_f:=e_C(1/f)=\sum_i f^{-2^i}/D_i\ \in k_\infty .$$
Over $\mathbb F_2$ every nonzero polynomial is monic, so $1/f=x^F(1+O(x))$. By Frobenius, $f^{-2^n}=x^{F2^n}\bigl(1+O(x^{2^n})\bigr)$.

Replace $S_N$ by $S^f_N:=\sum_{i\le N}f^{-2^i}/D_i$. Then:

- $f^{2^N}D_N\,S^f_N\in\mathbb F_2[\theta]$, and $B_f:=\deg(f^{2^N}D_N)=(N+F)2^N$.
- $e_f-S_N^f=x^{Y_f}\bigl(1+O(x^{2^N})\bigr)$ with $Y_f:=(N+1+F)2^{N+1}$. Indeed, the term $i=N+1$ is $x^{F2^{N+1}}(1+O(x^{2^{N+1}}))\cdot x^{(N+1)2^{N+1}}(1+O(x^{2^N}))$, and the terms with $i\ge N+2$ are $O(x^{(N+2+F)2^{N+2}})$.
- *Window.* Choose $k$ with $2^k\ge F+2$, and put $N:=2^k-1-F\ge1$ and $M:=N+k-1$. Then $A=2^{N+k}-1=(N+F+1)2^N-1$, so $2A+2=Y_f$ and $B_f-A=1-2^N$.

Set $Z_k:=Q_M\,f^{2^N}D_N\,(aV_MS^f_N+bV_M+cS^f_N+d)$. The proof of §6 then goes through verbatim, with these replacements:
- $e$ by $e_f$; since $\deg e_f=-F\le 0$, we still have $\deg\alpha,\deg\beta\le C$;
- $D_N$ by $f^{2^N}D_N$;
- $B$ by $B_f$.

The second round uses only $u\notin k$. $\square$

### 7.3. Proof of Theorem 1.3(2): coefficient fields

Let $L\subset\mathbb C_\infty$ be a finite extension of $k$ such that some (equivalently, every) $k$-basis $\omega_1,\dots,\omega_n$ of $L$ is linearly independent over $k_\infty$. This is what "linearly disjoint" means here.

Suppose $\sum_{j=1}^4\ell_jc_j=0$ with $\ell_j\in L$, where $(c_1,c_2,c_3,c_4)=(1,\tilde\pi,e_f,\tilde\pi e_f)$ and $f\neq0$. Write $\ell_j=\sum_i r_{ij}\omega_i$ with $r_{ij}\in k$. Then
$$\sum_i\omega_i\Bigl(\sum_j r_{ij}c_j\Bigr)=0,\qquad \sum_j r_{ij}c_j\in k_\infty ;$$
the second fact uses $q=2$, through $\tilde\pi\in k_\infty$. Hence $\sum_jr_{ij}c_j=0$ for every $i$, and Theorem 1.1 or Theorem 1.3(1) gives $r_{ij}=0$.

If $L/k$ is separable, then $L\otimes_k k_\infty\cong\prod_{w\mid\infty}L_w$ [Ser79, Ch. II, §3]. With a single place $w$ above $\infty$ this is a field, so the natural nonzero map $L\otimes_kk_\infty\to\mathbb C_\infty$ is injective, which is our hypothesis.

Examples:
- the constant field extensions $\mathbb F_{2^r}(\theta)$, where $\infty$ is inert. Their union is $\overline{\mathbb F}_2(\theta)$, and any relation involves only finitely many coefficients;
- $k(y)$ with $y^2+y=\theta$, where $\infty$ is totally ramified. $\square$

> **Remark 7.1 (limits of the method).**
> 1. Consider $e_C(g/f)$ with $g/f$ in lowest terms and $G:=\deg g\ge1$, for instance $e_C(\theta)=\theta e+e^2$. The analogous window would be $(N+F)2^N<A<(N+F+2-2G)\,2^N$, which is empty.
> 2. Suppose $\infty$ splits in $L$. Then $L$ may lie in $k_\infty$, and relations over $L$ no longer reduce to relations over $k$. A product-formula version of the argument would have to control $Z_k$ at all places above $\infty$, and the available margin $2^N$ is far smaller than $\deg(Q_MD_N)\approx 2N2^N$.

---

## 8. Remarks and open problems

### 8.1. The method is exactly at the threshold

The two approximations have qualities $|Q_Mu-P_M|=|Q_M|^{-1}$ and $|D_Ne-D_NS_N|=|D_N|^{-1-2/N}$.

Consider a relation $P(u,e)=0$ with $\deg_uP=d_1$ and $\deg_eP=d_2$. The same construction produces a polynomial whose degree is at most about $d_1A+d_2B-\min(2A,Y)$, and we need this to tend to $-\infty$. In the window, $A=B+2^N-1$ and $Y=2B+2^{N+1}$, so this happens only for $(d_1,d_2)\in\{(1,0),(0,1),(1,1),(0,2)\}$; the case $(1,1)$ requires the tie-breaking of §6.

Nothing with $d_1\ge2$ is reachable: the term $u^2$ costs $2A$, the full size of the error of $V_M$. In particular, the method cannot decide whether $\tilde\pi+e_C(1)$ is quadratic over $k$.

### 8.2. The case $q\ge 3$

For $q\ge3$ the naive approximations already suffice. The partial products $U_M=\prod_{m\le M}(1-x^{q^m-1})^{-1}$ have denominators of degree $A_M=\frac{q^{M+1}-q}{q-1}-M$ and $\deg(u-U_M)=1-q^{M+1}$. On the other side, $\deg(e-S_N)=-(N+1)q^{N+1}$ and $\deg D_N=Nq^N$.

Take $N=q^k$ and $M=N+k$. The two error contributions to $Z$ then have degrees
$$-\tfrac{q^2-3q+1}{q-1}\,Nq^N+O(M)\qquad\text{and}\qquad\text{about } q^{N+1}\text{ less}.$$
No tie occurs, and the argument of §6 gives the $k$-linear independence of $1,u,e,ue$ with no second-order information. Equivalently, $1,\tilde\pi,e_C(1),\tilde\pi e_C(1)$ are independent over $k(\lambda)$.

This is a special case of Denis's theorem. We mention it to stress that $q=2$ is precisely the case where the naive approximations fail ($\frac{q^2-3q+1}{q-1}=-1$ for $q=2$) and the approximations $V_M$ are needed.

### 8.3. What fails in characteristic zero

The proof of Theorem 1.1 is a Liouville argument. The principle (2.1) has the classical counterpart "a nonzero integer has absolute value at least $1$". The argument has four inputs; only the second lacks a known classical version.

**(E) Explicit approximations to $e_C(1)$, at known scales, that beat Dirichlet's exponent.** These are the convergents $S_N$ (Lemma 4.1, Remark 4.2). Classically, Euler's continued fraction $e=[2;1,2,1,1,4,1,1,6,\dots]$ provides explicit convergents $p_n/q_n$ preceding the partial quotients $a_{n+1}=2m$, with explicitly known error. So this input exists classically, although with a much smaller excess: $a_{n+1}\asymp\log q_n/\log\log q_n$, whereas here the spike $|[N+1]|=|D_N|^{2/N}$ is a power of the denominator.

**(P) Explicit approximations to $\tilde\pi$ of Dirichlet quality, inside the windows of $e_C(1)$, with an exact error term.** Lemma 3.1 gives $|Q_Mu-P_M|=|Q_M|^{-1}$, with $\deg Q_M=A$ in the window $(B,Y-B)=(N2^N,(N+2)2^N)$, and with $u-V_M$ known to second order.

The classical analogue would be an explicit family of fractions $P/Q$ such that:
- $|Q\pi-P|\le c/Q$;
- $Q\in(q_n,a_{n+1}q_n)$ for infinitely many of Euler's indices $n$;
- there is an asymptotic formula for $Q\pi-P$ precise enough to identify the leading term of the integer $Z_n=Q\,q_n\bigl(a\frac PQ\frac{p_n}{q_n}+b\frac PQ+c\frac{p_n}{q_n}+d\bigr)$.

Two comments on (P).

1. *Existence is not the issue.* Approximations of the required quality always exist: $\pi$'s own convergents provide them in every window. With them $Z_n\to0$, hence $Z_n=0$ for large $n$. But this is no contradiction: if $\pi$ were a Möbius image of $e$, its convergents would be Möbius images of those of $e$, and $Z_n$ would vanish. The function-field proof derives a contradiction only because Lemma 3.1 identifies the leading coefficient $\gamma$ of $Z_k$ in closed form.
2. *Quality.* Classically the windows have multiplicative width only $a_{n+1}$, so an approximation with $|Q\pi-P|\approx Q^{-\nu}$ is usable only if $\nu\ge1-O(\log\log q_n/\log q_n)$. The best explicit rational approximations to $\pi$ known have $\nu\approx 1/(\mu-1)\approx0.164$; these are the integrals behind the bound $\mu(\pi)\le7.1032$ of Zeilberger and Zudilin [ZZ20]. Moreover, an explicit family of quality $\nu\to1$ meeting all sufficiently large windows would give, by a standard transfer argument, an *effective* proof that the irrationality exponent of $\pi$ is $2$.

**(T) Tie-breaking.** In characteristic two the two error terms tie (Remark 5.2), and the second-order term resolves the tie. Classically one would need the analogous precision.

**(S) Symmetry.** Denis's route for $q\ge3$ needs an automorphism that fixes the period and moves $e_C(1)$. Neither $\mathbb C$ nor $\mathbb F_2[\theta]$ has one: the continuous automorphisms of $\mathbb C$ fix $\pi$ and $e$, and for $\mathbb F_2[\theta]$ see Proposition 8.1.

Thus, in the toy model, the whole difficulty of "$\pi+e\notin\mathbb Q$" is concentrated in the single explicit-approximation statement (P) about $\pi$ alone. Current technology is far from the quality it requires.

### 8.4. A symmetry, and the parity of partial quotients

> **Proposition 8.1.** *Let $\iota$ be the continuous automorphism of $k_\infty$ with $\iota(\theta)=\theta+1$, so that $\iota(x)=x/(1+x)$. Then $\iota(e_C(1))=e_C(1)$ and $\iota(\tilde\pi)=\tilde\pi$. Consequently, every partial quotient of the continued fractions of $\tilde\pi$ and of $e_C(1)$ lies in $\mathbb F_2[\theta^2+\theta]$, and in particular has even degree.*

*Proof.* The map $\iota$ preserves $\mathbb F_2[\theta]$ and degrees. It fixes each $D_i$, because $(\theta+1)^{2^i}+(\theta+1)^{2^j}=\theta^{2^i}+\theta^{2^j}$. By continuity, $\iota(e_C(z))=e_C(\iota(z))$ for $z\in k_\infty$. Hence $\iota(e_C(1))=e_C(1)$.

Also $e_C(\iota(\tilde\pi))=\iota(e_C(\tilde\pi))=0$, so $\iota(\tilde\pi)\in\ker e_C=\mathbb F_2[\theta]\,\tilde\pi$. As $\deg\iota(\tilde\pi)=\deg\tilde\pi$, it follows that $\iota(\tilde\pi)=\tilde\pi$.

Now let $\alpha\in k_\infty$ be $\iota$-invariant, with continued fraction $[a_0;a_1,\dots]$. Then $[\iota(a_0);\iota(a_1),\dots]$ is a continued fraction expansion of $\iota(\alpha)=\alpha$, since its partial quotients are polynomials of positive degree. By uniqueness, $\iota(a_n)=a_n$.

Finally, $\mathbb F_2[\theta]^\iota=\mathbb F_2[s]$ with $s=\theta^2+\theta$. Indeed, $k/\mathbb F_2(s)$ is Galois of degree $2$ with group $\{1,\iota\}$, since $\theta$ is a root of $X^2+X+s$; so $k^\iota=\mathbb F_2(s)$. Moreover $\mathbb F_2[\theta]\cap\mathbb F_2(s)=\mathbb F_2[s]$, because $\mathbb F_2[\theta]$ is integral over the integrally closed ring $\mathbb F_2[s]$. A nonzero element of $\mathbb F_2[s]$ has even degree in $\theta$. $\square$

For $e_C(1)$ this is also visible in Thakur's continued fraction, since $[j]=\sum_{i<j}s^{2^i}\in\mathbb F_2[s]$. Since $\iota$ is the only nontrivial automorphism of $\mathbb F_2[\theta]$, Proposition 8.1 shows that the symmetry Denis uses for $q\ge3$ has no analogue at $q=2$.

> **Proposition 8.2 (parity transfer).** *Let $\alpha\in k_\infty\smallsetminus k$ be such that all partial quotients of $\alpha$ of positive index have even degree, and let $h\in\mathbb F_2[\theta]$ have degree $j\ge1$. Then every partial quotient of positive index of $\alpha/h$, and of $h\alpha$, has degree $\le j$ or degree $\equiv j\pmod 2$.*

*Proof.* We use two standard facts about $\xi\in k_\infty\smallsetminus k$ with convergents $p_n/q_n$ and partial quotients $a_n$:

- (i) $|\xi-p_n/q_n|=|q_n|^{-2}|a_{n+1}|^{-1}$;
- (ii) (Legendre) if $P,Q\in\mathbb F_2[\theta]$ are coprime and $|\xi-P/Q|<|Q|^{-2}$, then $P/Q=p_n/q_n$ for some $n$.

For (ii), choose $n$ with $\deg q_n\le\deg Q<\deg q_{n+1}$. If $P/Q\neq p_n/q_n$, then
$$1\le|Pq_n-Qp_n|\le\max\bigl(|q_n||Q\xi-P|,\ |Q||q_n\xi-p_n|\bigr)<1,$$
a contradiction.

Now let $\beta=\alpha/h$. Let $p/q$ be a convergent of $\beta$, and let $d\ge1$ be the degree of the following partial quotient. By (i), $|\alpha-hp/q|=2^{j-2\deg q-d}$. Put $g=\gcd(h,q)$; then $hp/g$ and $q/g$ are coprime. There are two cases.

- If $d>j-2\deg g$, then $|\alpha-hp/q|<|q/g|^{-2}$. By (ii), $hp/q$ is a convergent of $\alpha$ with denominator $q/g$. By (i), the following partial quotient of $\alpha$ has degree $d'=d-j+2\deg g$. Since $d'$ is even, $d\equiv j\pmod 2$.
- Otherwise $d\le j-2\deg g\le j$.

For $h\alpha$, argue in the same way, using the approximations $p/(hq)$ to $\alpha$ and $g=\gcd(h,p)$. $\square$

> **Corollary 8.3.** *For $q=2$, every partial quotient of positive index of $u=\tilde\pi/\theta^2$ has degree $1$ or even degree. Every partial quotient of positive index of $\tilde\pi/\theta$ has odd degree.*

*Proof.* Apply Proposition 8.2, using Proposition 8.1, with $h=\theta^2$ and with $h=\theta$, respectively. $\square$

Corollary 8.3 explains an empirical observation. Among the first $14305$ partial quotients of $u$ (denominators up to degree $39998$), no odd degree $\ge3$ occurs, while about $45\%$ have degree $1$.

### 8.5. Open problems

1. **Algebraic independence for $q=2$.** Are $\tilde\pi$ and $e_C(1)$ algebraically independent over $\mathbb F_2(\theta)$, as Denis's Carlitz–Schanuel conjecture [Pel17, Conj. 2] predicts? Even: is $\tilde\pi+e_C(1)$ quadratic over $\mathbb F_2(\theta)$? Denis's symmetry is unavailable (Proposition 8.1), and the present method stops at bidegree $(1,1)$ (§8.1).
2. **Other arguments of the exponential.** Extend Theorem 1.1 to $e_C(g/f)$ with $\deg g\ge1$, for example $e_C(\theta)$, or to coefficient fields in which $\infty$ splits (Remark 7.1).
3. **The continued fraction of $\tilde\pi$.** Corollary 8.3 settles the parity law for $u$. Does the continued fraction of $\tilde\pi$, whose partial quotients lie in $\mathbb F_2[\theta^2+\theta]$, have an explicit description, like Thakur's for $e_C(1)$? Is the approximation exponent of $\tilde\pi$ equal to $2$? The continued fraction data suggest so.
4. **Characteristic zero.** Find explicit approximations to $\pi$ of Dirichlet quality, with controlled error, at the scales of Euler's continued fraction of $e$ (§8.3). This is the single missing input. A family meeting all sufficiently large windows would also give an effective proof of $\mu(\pi)=2$.

---

## 9. Numerical verification

The proofs above are complete and do not depend on computation. As a safeguard, every series identity and degree count used above was checked by exact arithmetic in $\mathbb F_2[[x]]$ and $\mathbb F_2[\theta]$ (python-flint), in two independent implementations:

1. the scripts of the project's function-field study: `scripts/04-funcfield/verify_degrees.py`, `verify_q2_windows.py`, `check_thakur_cf.py`, `cf_experiments.py` and `hp_rank_check.py` (logs in `data/04-funcfield/`);
2. a separate re-implementation written for this note: `scripts/22-write-ff2/verify_q2_paper.py`, which runs in about 8 seconds; its log is `data/22-write-ff2/verify_q2_paper.log`.

The table lists what the second script checks.

| Statement | Check | Range |
|---|---|---|
| Lemma 2.1 | $\deg D_i=i2^i$, $D_i=[i]D_{i-1}^2$, $D_j\mid D_i$, $x$-expansion of $1/D_i$ | $i\le4$ |
| (2.2) | $e_C(\theta^2u)=0$ and $e_C(\theta u)=\theta$ | mod $x^{20000}$ |
| Remark 2.2 | $\log_C(1)=u/(1+x)$; $C_{[1]}(1)=0$ | mod $x^{20000}$ |
| Lemma 3.1 | $\deg Q_M=A$; identity 2; valuations $2A,\ 4A+1,\ 4A+3$ in part 3 | $M\le14$ ($A\le32767$) |
| Lemma 4.1 | $\operatorname{val}(e-S_N)=Y$; next term at $Y+2^N$ | $N\le10$ |
| Lemma 5.1, Remark 5.2 | the identities; uniqueness of the window | $k\le16$; $N<40$ |
| §7.2 | tail and window for $e_C(1/f)$ | six $f$ of degree $1$–$4$, $N\le 14$ |
| Prop. 6.2 (mechanism) | for random $a,b,c,d$, $(aVS+bV+cS+d)-(aue+bu+ce+d)-x^{2A}\gamma$ has valuation exactly $2A+2+2^N+\operatorname{val}\beta$ | $k=2,3$ |
| Remark 4.2 | Thakur's continued fraction, exactly; $q_{2^N-1}=D_N$ | 986 partial quotients |
| Prop. 8.1 | all partial quotients of $\tilde\pi$ and $e_C(1)$ are $\iota$-invariant | 9554 and 431 partial quotients |
| Prop. 8.2, Cor. 8.3 | parity for $\tilde\pi/\theta^j$ ($j=1,2,3$), $\tilde\pi/(\theta^2+\theta+1)$, $\theta\tilde\pi$, $\theta^2\tilde\pi$; $u$ up to $\deg q\approx 4\cdot10^4$ | $\ge2290$ partial quotients each; 14305 for $u$ |
| Thm. 1.1 (sanity) | no relation $aue+bu+ce+d=0$ with $\deg a,b,c,d\le300$ | rank $1204$ (full) |

The last line is only a consistency check, since it is implied by Theorem 1.1. The earlier scripts agree on all overlapping entries. In addition, they check the leading term in Lemma 3.1(3) for $M=15$, and Lemma 4.1 together with the window for $N=15$ (the case $k=4$ of Lemma 5.1).

---

## Acknowledgments

This note was written by Claude Opus 5.5, an AI model by Anthropic. Claude was used to develop and check the arguments, to search and verify the literature, and to write and run the verification code described in §9. Ryan Xu set up the project and released this note; no human expert has verified the arguments. Corrections: ryan@ryanxu.net.

---

## References

- **[ABP04]** G. W. Anderson, W. D. Brownawell, M. A. Papanikolas, *Determination of the algebraic relations among special $\Gamma$-values in positive characteristic*, Ann. of Math. (2) **160** (2004), no. 1, 237–313.
- **[Car35]** L. Carlitz, *On certain functions connected with polynomials in a Galois field*, Duke Math. J. **1** (1935), 137–168.
- **[Den95]** L. Denis, *Indépendance algébrique et exponentielle de Carlitz*, Acta Arith. **69** (1995), no. 1, 75–89.
- **[Gos96]** D. Goss, *Basic Structures of Function Field Arithmetic*, Ergebnisse der Mathematik und ihrer Grenzgebiete (3) 35, Springer, Berlin, 1996.
- **[Pap08]** M. A. Papanikolas, *Tannakian duality for Anderson–Drinfeld motives and algebraic independence of Carlitz logarithms*, Invent. Math. **171** (2008), no. 1, 123–174.
- **[Pel17]** F. Pellarin, *On a variant of Schanuel conjecture for the Carlitz exponential*, J. Théor. Nombres Bordeaux **29** (2017), no. 3, 845–873 (arXiv:1610.04048).
- **[Ser79]** J.-P. Serre, *Local Fields*, Graduate Texts in Mathematics 67, Springer, New York, 1979.
- **[Tha92]** D. S. Thakur, *Continued fraction for the exponential for $\mathbb F_q[T]$*, J. Number Theory **41** (1992), no. 2, 150–155.
- **[Tha04]** D. S. Thakur, *Function Field Arithmetic*, World Scientific, River Edge, NJ, 2004.
- **[Tha11]** D. S. Thakur, *Higher Diophantine approximation exponents and continued fraction symmetries for function fields*, Proc. Amer. Math. Soc. **139** (2011), no. 1, 11–19.
- **[Thi92]** A. Thiery, *Indépendance algébrique des périodes et quasi-périodes d'un module de Drinfeld*, in: D. Goss, D. R. Hayes, M. I. Rosen (eds.), *The Arithmetic of Function Fields* (Columbus, OH, 1991), Ohio State Univ. Math. Res. Inst. Publ. 2, de Gruyter, Berlin, 1992, 265–284.
- **[Wad41]** L. I. Wade, *Certain quantities transcendental over $GF(p^n,x)$*, Duke Math. J. **8** (1941), 701–720.
- **[ZZ20]** D. Zeilberger, W. Zudilin, *The irrationality measure of $\pi$ is at most $7.103205334137\ldots$*, Moscow J. Combin. Number Theory **9** (2020), no. 4, 407–419.
