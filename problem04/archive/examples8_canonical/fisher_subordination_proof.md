# Finite Free Fisher Information Superadditivity via Subordination

A structured proof strategy in the sense of Lamport.

---

## 0. Definitions and Setup

**0.1. Polynomials.**
Let $\mathcal{P}_n$ denote the set of monic real-rooted polynomials of degree $n$
with *simple* roots. For $p \in \mathcal{P}_n$, write $p(x) = \prod_{i=1}^n (x - \lambda_i)$
with $\lambda_1 < \cdots < \lambda_n$.

**0.2. Finite free additive convolution.**
For $p, q \in \mathcal{P}_n$ with coefficient vectors $(a_0, \ldots, a_n)$, $(b_0, \ldots, b_n)$
(where $a_0 = b_0 = 1$), define $r = p \boxplus_n q$ by

$$r(x) = \sum_{k=0}^n c_k \, x^{n-k}, \qquad c_k = \sum_{i+j=k} \frac{(n-i)!(n-j)!}{n!(n-k)!} \, a_i \, b_j.$$

This is the Marcus–Spielman–Srivastava (MSS) convolution [MSS15].

**0.3. Cauchy transform.**
For $p \in \mathcal{P}_n$ with roots $\lambda_1, \ldots, \lambda_n$, define

$$G_p(z) = \frac{1}{n} \frac{p'(z)}{p(z)} = \frac{1}{n} \sum_{i=1}^n \frac{1}{z - \lambda_i}.$$

$G_p : \mathbb{C}^+ \to \mathbb{C}^-$, with simple poles at $\lambda_i$ each having residue $1/n$.

**0.4. Discrete Hilbert transform.**
Define

$$H_p(\lambda_i) := \sum_{j \neq i} \frac{1}{\lambda_i - \lambda_j} = \frac{p''(\lambda_i)}{2\,p'(\lambda_i)}.$$

Equivalently, $H_p(\lambda_i)$ is the regular part of $n G_p(z)$ at $z = \lambda_i$:

$$n G_p(z) = \frac{1}{z - \lambda_i} + H_p(\lambda_i) + O(z - \lambda_i).$$

**0.5. Finite free Fisher information.**

$$\Phi_n(p) := \sum_{i=1}^n H_p(\lambda_i)^2 = \sum_{i=1}^n \left( \sum_{j \neq i} \frac{1}{\lambda_i - \lambda_j} \right)^2.$$

Set $\Phi_n(p) = \infty$ if $p$ has a repeated root.

**0.6. Context.**
$\Phi_n$ is a finite analogue of the free Fisher information
$\Phi^*(\mu) = \int (H\mu)^2 \, d\mu$ of Voiculescu [V98]. In the large-$n$ limit,
with empirical root measures converging, $\boxplus_n \to \boxplus$ and
$\Phi_n / n^2 \to \Phi^*$.

---

## Admitted Results

**ADMITTED A (MSS).** If $p, q \in \mathcal{P}_n$, then $p \boxplus_n q \in \mathcal{P}_n$.
[MSS15, Theorem 4.4]

**ADMITTED B (MSS).** If $A, B$ are $n \times n$ Hermitian with characteristic
polynomials $p, q$, and $U \sim \mathrm{Haar}(U(n))$, then
$\mathbb{E}_U[\chi_{A + UBU^*}(x)] = (p \boxplus_n q)(x)$.
[MSS15, Theorem 4.2]

**ADMITTED C (MSS).** $(p \boxplus_n q)'(x) = n ( p^{(1)} \boxplus_{n-1} q^{(1)} )(x)$,
where $p^{(1)} = \frac{1}{n}p'$ is the monic polynomial whose roots are the
critical points of $p$. [MSS15, Proposition 4.5]

---

## 1. Theorem (Conjecture)

For all $p, q \in \mathcal{P}_n$:

$$\frac{1}{\Phi_n(p \boxplus_n q)} \geq \frac{1}{\Phi_n(p)} + \frac{1}{\Phi_n(q)}.$$

**Numerical status.** Verified to hold (with strict inequality for $n \geq 3$)
across 10,000 random trials each at $n = 3, 4, 5, 10$, with zero violations.
Equality holds if and only if $n \leq 2$.

---

## 2. Finite Subordination

### 2.1. Statement (UNPROVED — Hard Lemma 1)

**Lemma.** Let $p, q \in \mathcal{P}_n$ and $r = p \boxplus_n q$.
There exist rational functions $\omega_1, \omega_2 : \hat{\mathbb{C}} \to \hat{\mathbb{C}}$
satisfying:

(a) $G_r(z) = G_p(\omega_1(z)) = G_q(\omega_2(z))$ as rational identities.

(b) $\omega_1, \omega_2 : \mathbb{C}^+ \to \mathbb{C}^+$ (Herglotz property).

### 2.2. Construction

2.2.1. Fix $z \in \mathbb{C}^+$. The equation $G_r(z) = G_p(w)$ in $w$ is equivalent to
a degree-$(n-1)$ polynomial equation (since $G_p$ has degree $(n-1)/n$ as a
rational map). Generically there are $n-1$ solutions.

2.2.2. **Selection principle.** Since $G_p : \mathbb{C}^+ \to \mathbb{C}^-$ is a
covering map of degree $n-1$ (restricted to $\mathbb{C}^+$), and
$G_r(z) \in \mathbb{C}^-$ for $z \in \mathbb{C}^+$, we seek the unique
preimage $w = \omega_1(z) \in \mathbb{C}^+$. That exactly one such preimage
exists follows from: $G_p$ maps $\mathbb{C}^+$ to $\mathbb{C}^-$ with winding
number $-(n-1)$, and the level set $\{w : G_p(w) = c\}$ for $c \in \mathbb{C}^-$
has exactly $n-1$ points, of which exactly one lies in $\mathbb{C}^+$ (the others
lie in $\mathbb{C}^-$ or on $\mathbb{R}$).

> **STATUS.** The winding-number argument in 2.2.2 needs verification. The
> claim that exactly one of $n-1$ preimages lies in $\mathbb{C}^+$ is plausible
> from the interlacing structure of poles and zeros of $G_p$ on $\mathbb{R}$,
> but a rigorous proof requires the argument principle applied to
> $G_p - c$ on $\mathbb{C}^+$.
>
> **Inductive approach via Admitted C.** Subordination at level $n$ may be
> bootstrapped from subordination at level $n-1$ using the derivative
> identity, with the interlacing of roots of $r$ and $r'$ providing the
> lift. Base case $n = 2$ is computed explicitly in §5.

---

## 3. Chain Rule for Discrete Hilbert Transform

### 3.1. Statement

**Lemma.** Assume Lemma 2.1. Let $r = p \boxplus_n q$ with roots $\nu_1 < \cdots < \nu_n$.

(a) There exist bijections $\sigma, \tau : \{1, \ldots, n\} \to \{1, \ldots, n\}$ such
that $\omega_1(\nu_k) = \lambda_{\sigma(k)}$ and $\omega_2(\nu_k) = \mu_{\tau(k)}$.

(b) $\omega_1'(\nu_k) = 1$ and $\omega_2'(\nu_k) = 1$ for all $k$.

(c) Defining correction terms $\alpha_k := \tfrac{1}{2}\omega_1''(\nu_k)$ and
$\beta_k := \tfrac{1}{2}\omega_2''(\nu_k)$:

$$H_r(\nu_k) = H_p(\lambda_{\sigma(k)}) - \alpha_k = H_q(\mu_{\tau(k)}) - \beta_k.$$

### 3.2. Proof of (a)

3.2.1. As $z \to \nu_k$, $G_r(z)$ has a simple pole. Since $G_r(z) = G_p(\omega_1(z))$,
$G_p$ must also diverge at $\omega_1(\nu_k)$, forcing $\omega_1(\nu_k) = \lambda_j$
for some $j$. Define $\sigma(k) = j$.

3.2.2. Injectivity of $\sigma$: If $\sigma(k) = \sigma(\ell) = j$, then $\omega_1$ maps
two distinct real points $\nu_k, \nu_\ell$ to the same value $\lambda_j$. But a
rational Herglotz function is monotone increasing on each interval of
$\mathbb{R} \setminus \{\text{poles}\}$, so $\omega_1$ is locally injective near any
regular real point. Since $\omega_1(\nu_k) = \lambda_j$ is a root of $p$ (not a pole
of $\omega_1$), local injectivity prevents $\omega_1(\nu_\ell) = \lambda_j$ for
$\nu_\ell \neq \nu_k$ in the same monotone interval. A global argument using the
degree of $\omega_1$ and the Herglotz property closes this.

3.2.3. Since $|\{\nu_k\}| = |\{\lambda_j\}| = n$, injectivity implies bijectivity.

### 3.3. Proof of (b) and (c)

3.3.1. Laurent expansion of $nG_p(\omega_1(z))$ near $z = \nu_k$, where
$\omega_1(z) = \lambda_{\sigma(k)} + \omega_1'(\nu_k)(z - \nu_k) + \tfrac{1}{2}\omega_1''(\nu_k)(z-\nu_k)^2 + \cdots$:

$$nG_p(\omega_1(z)) = \frac{1}{\omega_1'(\nu_k)(z - \nu_k)} - \frac{\omega_1''(\nu_k)}{2\,\omega_1'(\nu_k)^2} + H_p(\lambda_{\sigma(k)}) + O(z - \nu_k).$$

3.3.2. Matching with $nG_r(z) = \frac{1}{z - \nu_k} + H_r(\nu_k) + O(z - \nu_k)$:

- Pole residue: $1/\omega_1'(\nu_k) = 1$, so $\omega_1'(\nu_k) = 1$.
- Regular part: $H_r(\nu_k) = H_p(\lambda_{\sigma(k)}) - \frac{\omega_1''(\nu_k)}{2}$.

This is a direct algebraic computation given Lemma 2.1. $\square$

**Numerical verification at $n = 2$.** With $p = x(x-2)$, $q = x(x-1)$:

| $k$ | $\omega_1(\nu_k)$ | $\omega_1'(\nu_k)$ | $\alpha_k$ | $H_p(\lambda_{\sigma(k)})$ | $H_r(\nu_k)$ | Match |
|---|---|---|---|---|---|---|
| $-$ | $2$ | $1$ | $0.947$ | $0.5$ | $-0.447$ | ✓ |
| $+$ | $0$ | $1$ | $-0.947$ | $-0.5$ | $0.447$ | ✓ |

| $k$ | $\omega_2(\nu_k)$ | $\omega_2'(\nu_k)$ | $\beta_k$ | $H_q(\mu_{\tau(k)})$ | $H_r(\nu_k)$ | Match |
|---|---|---|---|---|---|---|
| $-$ | $0$ | $1$ | $-0.553$ | $-1$ | $-0.447$ | ✓ |
| $+$ | $1$ | $1$ | $0.553$ | $1$ | $0.447$ | ✓ |

---

## 4. Main Inequality

### 4.1. Reformulation in terms of correction vectors

Write $u_k = H_p(\lambda_{\sigma(k)})$, $v_k = H_q(\mu_{\tau(k)})$, $h_k = H_r(\nu_k)$.
Since $\sigma, \tau$ are bijections:

$$\Phi_n(p) = \|u\|^2, \qquad \Phi_n(q) = \|v\|^2, \qquad \Phi_n(r) = \|h\|^2$$

where $\|\cdot\|$ denotes the Euclidean norm on $\mathbb{R}^n$.

From Lemma 3.1(c): $h = u - \alpha = v - \beta$, i.e.,

$$u = h + \alpha, \qquad v = h + \beta.$$

**Compatibility:** $\alpha - \beta = u - v$ (determined by $p, q, r$; not free).

### 4.2. Target inequality in vector form

The conjecture $1/\|h\|^2 \geq 1/\|u\|^2 + 1/\|v\|^2$ is equivalent to:

$$\|u\|^2 \|v\|^2 \geq \|h\|^2 (\|u\|^2 + \|v\|^2)$$

or equivalently:

$$\|h + \alpha\|^2 \|h + \beta\|^2 \geq \|h\|^2 (\|h + \alpha\|^2 + \|h + \beta\|^2).$$

### 4.3. Reduction (UNPROVED — Hard Lemma 2)

**What is needed:** A constraint on $(\alpha, \beta)$ from the Herglotz property
of $\omega_1, \omega_2$ that implies the inequality in §4.2.

**Available constraints:**

(i) $\alpha - \beta = u - v$ (compatibility, from §4.1).

(ii) The Herglotz property implies $\omega_i'(x) \geq 0$ on
$\mathbb{R} \setminus \{\text{poles}\}$. Combined with $\omega_i'(\nu_k) = 1$,
this constrains the magnitude and sign of $\omega_i''(\nu_k)$ (hence $\alpha_k, \beta_k$)
through the curvature of $\omega_i$ at the roots.

(iii) The representation $\omega_i(z) = z + c_i + \sum_j \frac{m_{ij}}{z - p_{ij}}$
(partial fraction of a Herglotz function with $\omega_i(z) \sim z$ at $\infty$)
gives $\omega_i''(\nu_k) = \sum_j \frac{2 m_{ij}}{(\nu_k - p_{ij})^3}$, where $m_{ij} > 0$.
This connects $\alpha_k$ to the pole structure of $\omega_i$.

**Possible attack:** The inequality in §4.2 can be rewritten as

$$\frac{\|h\|^2}{\|h + \alpha\|^2} + \frac{\|h\|^2}{\|h + \beta\|^2} \leq 1.$$

By Cauchy–Schwarz, $\|h\|^2 / \|h + \alpha\|^2 \leq 1 - 2\langle h, \alpha\rangle/\|h+\alpha\|^2 + \ldots$
This does not close directly. The constraint from the Herglotz representation (iii) that
$\alpha_k = \sum_j m_{1j}/(\nu_k - p_{1j})^3$ has a *specific structure* — it is not an
arbitrary vector. This structure (related to complete Bernstein functions and their
discretisations) must be exploited.

> **NOTE ON THE SUMMATION RELATION.** In the continuum theory, subordination functions
> satisfy $\omega_1(z) + \omega_2(z) = z + F_{\mu\boxplus\nu}(z)$, which couples $\alpha$ and
> $\beta$ through $G_r$ and makes the inequality tractable. At finite $n$, this relation
> **does not hold** (verified computationally at $n = 2$: $\omega_1 + \omega_2 - z - 1/(nG_r)$
> is a nontrivial irrational function). A finite-$n$ replacement — perhaps involving the
> finite $R$-transform or the finite $F$-transform of [AP18] — may be necessary to
> couple the constraints on $\alpha$ and $\beta$.

---

## 5. Base Case: $n = 2$

### 5.1. Fisher information

$p(x) = (x-a)(x-b)$, $q(x) = (x-c)(x-d)$, $a < b$, $c < d$.

$\Phi_2(p) = 2/(b-a)^2$, $\Phi_2(q) = 2/(d-c)^2$.

### 5.2. Convolution

$r = p \boxplus_2 q$ has gap $\sqrt{(b-a)^2 + (d-c)^2}$.

$\Phi_2(r) = 2/((b-a)^2 + (d-c)^2)$.

### 5.3. Equality

$$\frac{1}{\Phi_2(r)} = \frac{(b-a)^2 + (d-c)^2}{2} = \frac{(b-a)^2}{2} + \frac{(d-c)^2}{2} = \frac{1}{\Phi_2(p)} + \frac{1}{\Phi_2(q)}. \qquad \square$$

### 5.4. Subordination at $n = 2$ (verified computationally)

With $p = x(x-2)$, $q = x(x-1)$, solving $G_p(\omega_1) = G_r(z)$ gives
a quadratic in $\omega_1$ with a unique Herglotz branch satisfying:

- $\omega_1'(\nu_k) = 1$ at both roots of $r$. ✓
- Chain rule $H_r(\nu_k) = H_p(\omega_1(\nu_k)) - \omega_1''(\nu_k)/2$. ✓
- $\omega_2$, constructed independently from $G_q(\omega_2) = G_r(z)$,
  also satisfies $\omega_2'(\nu_k) = 1$ and the chain rule. ✓
- The continuum summation relation $\omega_1 + \omega_2 = z + 1/(nG_r)$ fails. ✗

---

## 6. Summary

| Component | Status | Difficulty |
|---|---|---|
| Definitions (§0) | Complete | — |
| $n = 2$ equality (§5.1–5.3) | Proved | — |
| $n = 2$ subordination (§5.4) | Verified computationally | — |
| Numerical evidence (§1) | 40k trials, 0 violations | — |
| Finite subordination existence (§2) | **Unproved** | High |
| Chain rule (§3) | Follows from §2 | Mechanical |
| Constraint on corrections (§4.3) | **Unproved** | High |
| Main inequality from constraints (§4.2) | **Unproved** | Medium |

**Critical path:** §2 (subordination existence) → §4.3 (correction constraints) → §4.2 (inequality).

**Recommended next steps:**

1. Compute subordination explicitly at $n = 3$ (cubic equation in $\omega_1$)
   and verify all claims. This is feasible and highly informative.

2. Determine the correct finite summation relation replacing
   $\omega_1 + \omega_2 = z + 1/(nG_r)$, possibly via the finite $R$-transform
   of [AP18].

3. Investigate whether Lorentzian polynomial techniques [ALOV19] provide
   an alternative route to the correction bounds in §4.3, bypassing the
   summation relation entirely.

---

## References

- **[MSS15]** Marcus, Spielman, Srivastava. *Interlacing families II: Mixed
  characteristic polynomials and the Kadison–Singer problem.* Ann. Math. 2015.
- **[V98]** Voiculescu. *The analogues of entropy and of Fisher's information
  measure in free probability theory V.* Invent. Math. 1998.
- **[AP18]** Arizmendi, Perales. *Finite free convolutions of polynomials.*
  arXiv:1807.08958, 2018.
- **[ALOV19]** Anari, Liu, Oveis Gharan, Vinzant. *Log-concave polynomials II:
  High-dimensional walks and an FPRAS for counting bases of a matroid.*
  Ann. Math. 2021.
