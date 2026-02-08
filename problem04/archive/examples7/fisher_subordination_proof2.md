# Finite Free Fisher Information Superadditivity via Subordination

A structured proof strategy in the sense of Lamport.

---

## 0. Definitions and Setup

**0.1. Polynomials.**
Let $\mathcal{P}_n$ denote the set of monic real-rooted polynomials of degree $n$
with *simple* roots. For $p \in \mathcal{P}_n$, write $p(x) = \prod_{i=1}^n (x - \lambda_i)$
with $\lambda_1 < \cdots < \lambda_n$. Note that the result should later be extended to real-rooted polynomials with repeated roots.

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
