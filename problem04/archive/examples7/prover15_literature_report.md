# PROVER-15: Literature Search and Framework Discovery Report

**Date:** 2026-02-08
**Role:** Literature search + novel framework identification
**Target:** Fisher superadditivity conjecture for finite free convolution

---

## Executive Summary

The conjecture `1/Phi_n(p boxplus_n q) >= 1/Phi_n(p) + 1/Phi_n(q)` is the
**exact finite analog of Voiculescu's free Stam inequality**, which is a known
theorem in free probability theory. The infinite-dimensional version was proved
by Voiculescu (1998) and the proof mechanism has been extensively studied. The
key insight is that `Phi_n(p) = sum_i H_p(lambda_i)^2` is the finite analog of
Voiculescu's free Fisher information `Phi*(mu) = integral (H mu)^2 d mu`, and
the MSS finite free convolution `boxplus_n` converges to free additive convolution
`boxplus`. The conjecture therefore asks: does the free Stam inequality hold at
finite n, not just in the limit?

**No existing paper was found that proves this exact finite result.** However,
the structural framework for proving it exists and has been partially developed.
The most promising approach uses the **subordination method** combined with
**polynomial differentiation as free heat flow**.

---

## 1. The Free Stam Inequality (Known Theorem)

### Statement
For freely independent self-adjoint random variables X, Y with finite free
Fisher information:

```
1/Phi*(X + Y) >= 1/Phi*(X) + 1/Phi*(Y)
```

where `Phi*(X) = ||J(X)||^2` is the squared L^2 norm of the **conjugate variable**
(free score) `J(X)`.

### Connection to Our Conjecture
- Voiculescu's free Fisher information: `Phi*(mu) = integral (H_mu(x))^2 d mu(x)`
  where `H_mu` is the Hilbert transform of the measure mu.
- Our finite free Fisher information: `Phi_n(p) = sum_i H_p(lambda_i)^2`
  where `H_p(lambda_i) = sum_{j != i} 1/(lambda_i - lambda_j)`.
- The quantity `H_p(lambda_i)` is the discrete Hilbert transform of the
  empirical root measure `mu_p = (1/n) sum delta_{lambda_i}` evaluated at
  the atoms.
- As n -> infinity with empirical measures converging, `Phi_n/n^2 -> Phi*`.
- The MSS convolution `boxplus_n` converges to free convolution `boxplus`.
- Therefore our conjecture is the exact finite-n version of the free Stam
  inequality.

### Proof Technique (Infinite Case)
Voiculescu (1998, Inventiones Math.) proved the free Stam inequality using:
1. **Conjugate variables / free score**: J(X) is defined via the relation
   `tau(J(X) * p(X, X*)) = (d/dt)|_{t=0} tau(p(X + tS, (X+tS)*))` where S
   is a free semicircular.
2. For absolutely continuous measures, `J(X) = 2*pi*H_f(X)` where H_f is the
   Hilbert transform of the density f.
3. **Subordination**: For X + Y with X, Y freely independent, the Cauchy
   transform satisfies `G_{X+Y}(z) = G_X(omega_Y(z)) = G_Y(omega_X(z))` where
   omega_X, omega_Y are subordination functions with `omega_X(z) + omega_Y(z) = z`.
4. The free score of X+Y decomposes via conditional expectations:
   `J(X+Y) = E[J(X) | X+Y] + E[J(Y) | X+Y]`.
5. By Pythagoras (L^2 projection contracts): `||J(X+Y)||^2 <= ||J(X)||^2` and
   similarly for Y. The Stam inequality follows from the harmonic mean property
   of projections.

**Key references:**
- [V98] Voiculescu, "The analogues of entropy and of Fisher's information
  measure in free probability theory V," Invent. Math. 132 (1998), 189-227.
- [BS01] Biane, Speicher, "Free diffusions, free entropy and free Fisher
  information," Ann. I.H.P. Prob. Stat. 37 (2001), 581-606.

---

## 2. The Shlyakhtenko-Tao Framework (Key Paper)

### Paper
Shlyakhtenko, Tao (with appendix by Jekel), "Fractional free convolution
powers," arXiv:2009.01882 (2020).

### Main Results
1. **Monotonicity of free entropy** along normalized free convolution powers
   `mu^{boxplus k}` for k >= 1: `chi(k^{-1/2} * mu^{boxplus k})` is
   non-decreasing in k.
2. **Monotonicity of free Fisher information**: the Fisher information is
   non-increasing.
3. **Variational description** of the free convolution power process.

### Two Proofs
1. **Free probability proof**: Uses the free score function J(X). Key identity:
   `J([pXp] : [pBp]) = k * E([p J(X:B) p] | [pXp], [pBp])`.
   Pythagoras's theorem (conditional expectations contract in L^2 norm) gives
   monotonicity.
2. **Analytic proof**: Self-contained, using integration by parts and contour
   integration. Computes partial derivatives in parameter k using a Burgers-type
   equation and discovers a positive semi-definite kernel K(z,w) whose properties
   establish monotonicity.

### Relevance to Our Conjecture
The analytic proof technique is potentially adaptable to the finite setting
because:
- It works with Cauchy transforms and contour integrals
- The polynomial Cauchy transform `G_p(z) = (1/n) p'(z)/p(z)` is well-defined
- The positive semi-definite kernel approach may have a finite-n analog
- The connection to the minor process (which IS finite-dimensional) provides
  a bridge

---

## 3. Polynomial Differentiation as Free Heat Flow

### Key Connection (MSS Admitted C)
The MSS derivative formula states:
```
(p boxplus_n q)'(x) = n * (p^{(1)} boxplus_{n-1} q^{(1)})(x)
```
where `p^{(1)} = (1/n) p'` is the monic polynomial whose roots are the critical
points of p.

### Differentiation = Free Convolution with Semicircle
In the large-n limit, repeated differentiation of a polynomial is equivalent to
free convolution with a semicircular distribution (Steinerberger 2020,
Shlyakhtenko-Tao 2020). The empirical root measure after taking the (epsilon*n)-th
derivative converges to the original empirical measure raised to a free
convolution power.

### The Heat Flow Conjecture
Hall, Ho (arXiv:2202.09660): When polynomial coefficients evolve according to
the heat equation, the roots' large-N limiting behavior follows a first-order
nonlinear PDE (the free Burgers equation). This is the polynomial analog of the
classical de Bruijn identity.

### Finite de Bruijn Identity Candidate
The derivative formula `(p boxplus_n q)' = n (p^{(1)} boxplus_{n-1} q^{(1)})`
is the finite analog of "free convolution commutes with the free heat semigroup."
If we can express `Phi_n` in terms of how `p` changes under differentiation
(i.e., under convolution with the "polynomial semicircle"), this would give
a finite de Bruijn identity:
```
d/dt [entropy(p boxplus_n gamma_t)] = Phi_n(p boxplus_n gamma_t)
```
where gamma_t is a Gaussian polynomial. This would reduce the Stam inequality
to convexity of entropy.

---

## 4. Gribinski's Polynomial Entropy (Directly Relevant)

### Paper
Gribinski, "A notion of entropy on the roots of polynomials,"
arXiv:1907.12826 (2019, revised 2023).

### Main Result
Defines `Dis(p)` (related to the discriminant) as an entropy for polynomials
and proves:
```
Dis(p - t*p') > Dis(p)   for t != 0
```
This is **entropy increases under differentiation**, which is the polynomial
analog of "entropy increases along free heat flow."

### Connection to Our Problem
- Gribinski proves entropy increases under the finite free heat semigroup
- The free Stam inequality in the classical setting follows from:
  (a) entropy is concave along the heat flow (de Bruijn identity), AND
  (b) Fisher information is the derivative of entropy along the heat flow
- If Gribinski's result can be extended to show that `Phi_n` is the
  "derivative" of `Dis` along the polynomial heat flow, then concavity
  of `Dis` under free convolution would imply the Stam inequality

### Critical Gap
Gribinski does NOT define a finite free Fisher information or prove the finite
Stam inequality. His entropy monotonicity result is one of the two ingredients
needed; the other (the de Bruijn identity relating Phi_n to d/dt of entropy)
appears to be open.

---

## 5. Free Probability via Entropic Optimal Transport

### Paper
Boucheron, Gaetan, Guionnet (arXiv:2309.12196, 2023).

### Main Results
- Express logarithmic potentials of free convolution via entropic optimal
  transport on the symmetric group
- Derive inequalities relating free and classical operations
- Apply large deviation principles to MSS quadrature formulas

### Relevance
This paper works directly at the intersection of MSS convolution and
information theory. The entropic optimal transport formulation might provide
the right variational framework for proving finite information inequalities.
The connection to the symmetric group provides a finite-dimensional setting.

---

## 6. The 2025 Paper: Free Probabilistic Denoising Diffusion

### Paper
arXiv:2510.22778 (October 2025), "A Free Probabilistic Framework for Denoising
Diffusion Models: Entropy, Transport, and Reverse Processes."

### Main Results
- Develops free analogs of diffusion models using Voiculescu's theory
- Proves free logarithmic Sobolev, Talagrand, and HWI inequalities
- Establishes gradient-flow structure in noncommutative Wasserstein space
- Uses free Malliavin calculus and Clark-Ocone representation
- The conjugate variable (free score) drives the reverse-time SDE

### Relevance
This is the most recent paper connecting free Fisher information to dynamical
processes. The functional inequalities proved here (free log-Sobolev, HWI)
are stronger than the Stam inequality and hold in the infinite-dimensional
setting. Their finite-dimensional analogs would imply our conjecture.

---

## 7. The Szarek-Voiculescu Free Entropy Power Inequality

### Paper
Szarek, Voiculescu, "Volumes of restricted Minkowski sums and the free
analogue of the entropy power inequality," Comm. Math. Phys. 178 (1996), 563-570.

### Result
For freely independent X, Y:
```
exp(2 * chi(X + Y)) >= exp(2 * chi(X)) + exp(2 * chi(Y))
```
where chi is the free entropy. This is proved using the Brascamp-Lieb-Luttinger
rearrangement inequality applied to microstates (matricial approximants).

### Connection
The entropy power inequality is STRONGER than the Stam inequality (classically,
EPI implies Stam via differentiation). If a finite version of the free EPI could
be proved for polynomials, our conjecture would follow. The Szarek-Voiculescu
proof uses matricial microstates, which ARE finite-dimensional objects.

---

## 8. Recommended Proof Framework

Based on this literature survey, I recommend the following framework, in order
of promise:

### Framework A: Finite Subordination + L^2 Projection (MOST PROMISING)

**Idea:** Adapt Voiculescu's original proof of the free Stam inequality to the
finite polynomial setting.

**Steps:**
1. Define finite subordination functions omega_p, omega_q for
   `p boxplus_n q`: find analytic functions such that
   `G_{p boxplus_n q}(z) = G_p(omega_q(z))`.
   (These exist finitely because G_p is a rational function of degree n.)

2. Express `H_{p boxplus_n q}(lambda_i)` in terms of `H_p` and `H_q` via
   these subordination functions.

3. Prove the finite L^2 contraction:
   ```
   sum_i H_{p boxplus_n q}(rho_i)^2 <= sum_j H_p(lambda_j)^2
   ```
   where rho_i are roots of p boxplus_n q and lambda_j are roots of p.
   This is the finite Pythagoras/projection inequality.

4. Combine for both p and q to get the Stam inequality.

**Why this might work:** The Cauchy transform of a polynomial is just a
rational function. Subordination for rational functions is concrete algebra.
The key inequality is an L^2 contraction that may follow from the interlacing
properties that MSS already proved.

**Key obstacle:** Defining the "right" finite subordination. In the infinite
case, subordination functions are unique; for polynomials, we need to identify
the correct finite analog.

### Framework B: Polynomial de Bruijn + Gribinski Entropy

**Idea:** Complete Gribinski's program by:
1. Defining `Phi_n` as the derivative of `Dis(p)` (Gribinski's polynomial
   entropy) along the polynomial heat flow.
2. Proving concavity of `Dis(p boxplus_n gamma_t)` as a function of t.
3. Deriving the Stam inequality from concavity + the de Bruijn identity.

**Why this might work:** Gribinski already proved entropy monotonicity.
The missing piece is the quantitative relationship between Phi_n and d/dt Dis.

**Key obstacle:** Gribinski's entropy `Dis(p)` may not be exactly
`log(discriminant(p))` or may require modification to get the right de Bruijn
identity.

### Framework C: Analytic Kernel Method (Shlyakhtenko-Tao)

**Idea:** Adapt the self-contained analytic proof from the Shlyakhtenko-Tao
paper on fractional free convolution powers.

**Steps:**
1. Parameterize `r_t = p boxplus_n (t * q)` (if such fractional finite free
   convolution makes sense).
2. Compute `d/dt Phi_n(r_t)` using the Cauchy transform equation.
3. Show this derivative has a sign that implies the Stam inequality.

**Why this might work:** The Shlyakhtenko-Tao analytic proof already works
with Cauchy transforms and contour integrals, which are well-defined for
polynomials.

**Key obstacle:** Fractional finite free convolution powers may not be
well-defined for all parameters.

### Framework D: Random Matrix / Matricial Microstates

**Idea:** Use the random matrix representation (MSS Admitted B):
`E_U[chi_{A + UBU*}] = p boxplus_n q`

to reduce the polynomial inequality to a random matrix inequality.

**Steps:**
1. Express `Phi_n(E[chi_{A+UBU*}])` in terms of A, B, U.
2. Use Jensen's inequality or convexity arguments on the random matrix side.
3. The deterministic matrices A, B carry the "separated" information that
   should add up.

**Why this might work:** The random matrix representation gives a concrete
probabilistic interpretation. Jensen's inequality applied to convex functions
of random matrices is well-studied.

**Key obstacle:** `Phi_n` is defined on the expected polynomial, not on
individual realizations. The relationship between `Phi_n(E[chi])` and
`E[Phi_n(chi)]` is unclear.

---

## 9. Key Open Questions

1. **Does finite subordination exist for MSS convolution?** I.e., for
   `r = p boxplus_n q`, do there exist functions omega_p, omega_q such that
   `G_r(z) = G_p(omega_q(z))`? If so, what are their properties?

2. **Is Phi_n = d/dt Dis(p boxplus_n gamma_t)?** Where gamma_t is a
   Gaussian polynomial with variance t and Dis is Gribinski's entropy.

3. **Does the positive semi-definite kernel from Shlyakhtenko-Tao have a
   finite polynomial analog?**

4. **Can the Szarek-Voiculescu proof be "de-limited"?** Their proof uses
   volumes of sets in matrix space. Can this be made finite without taking
   n -> infinity?

---

## 10. Papers Directly Relevant (Annotated Bibliography)

| # | Paper | Year | Relevance |
|---|-------|------|-----------|
| 1 | Voiculescu, "Analogues of entropy and Fisher info V" | 1998 | **Proves free Stam inequality** (our conjecture's limit) |
| 2 | Szarek-Voiculescu, "Volumes of restricted Minkowski sums" | 1996 | **Free entropy power inequality** (stronger than Stam) |
| 3 | Shlyakhtenko-Tao, "Fractional free convolution powers" | 2020 | **Two proofs of Fisher monotonicity**, analytic method adaptable |
| 4 | Biane-Speicher, "Free diffusions, free entropy, free Fisher" | 2001 | **Free de Bruijn identity** (entropy derivative = Fisher info) |
| 5 | Gribinski, "A notion of entropy on roots of polynomials" | 2019 | **Polynomial entropy increases under differentiation** |
| 6 | Marcus, "Polynomial convolutions and finite free probability" | 2021 | **Finite free transforms** (K-transform, R-transform analogs) |
| 7 | Marcus-Spielman-Srivastava, "Finite free convolutions" | 2015/2022 | **Foundation**: derivative formula, interlacing |
| 8 | Arizmendi-Perales, "Cumulants for finite free convolution" | 2018 | **Finite free cumulants** (additive under boxplus_n) |
| 9 | Steinerberger, "Free convolution powers via roots" | 2020 | **Differentiation = free convolution** at polynomial level |
| 10 | Hall-Ho, "Heat flow conjecture for polynomials" | 2022 | **Polynomial heat flow** and limiting PDE |
| 11 | arXiv:2510.22778, "Free probabilistic denoising diffusion" | 2025 | **Free log-Sobolev, HWI inequalities** (stronger than Stam) |
| 12 | arXiv:2309.12196, "Free probability via entropic OT" | 2023 | **Entropic optimal transport** for MSS formulas |
| 13 | Arizmendi-Fujie-Perales-Ueda, "S-transform in finite free" | 2024 | **Finite free S-transform**, cumulant-differentiation |
| 14 | arXiv:2508.21483, "Finite N precursors of free cumulants" | 2025 | **Finite-N corrections** to free cumulant theory |
| 15 | Shlyakhtenko-Tao, "Flow of polynomial roots under diff." | 2020 | **PDE for root dynamics** = free Burgers equation |

---

## 11. Verdict

**The conjecture is NOT already known to be true in the literature.** It is the
finite-n analog of a known infinite-dimensional result (Voiculescu's free Stam
inequality), but the passage from n=infinity to finite n is non-trivial and has
not been carried out.

**The conjecture is very likely true** because:
1. It holds in the limit n -> infinity (free Stam inequality)
2. It is numerically verified for all tested n
3. The structural framework (subordination, score functions) exists
4. Related finite results (Gribinski's entropy monotonicity) have been proved

**The most promising proof strategy is Framework A (finite subordination):**
1. The free Stam inequality proof via subordination + L^2 contraction
   is the simplest proof in the infinite case
2. All ingredients (Cauchy transforms, Hilbert transforms at roots,
   subordination functions) make sense for polynomials
3. The interlacing theory of MSS may provide the missing L^2 contraction
   estimate

**The key technical step** is to show that for `r = p boxplus_n q` with
roots rho_1, ..., rho_n:
```
sum_i [H_r(rho_i)]^2  <=  (n_p / n) * sum_j [H_p(lambda_j)]^2
                         + (n_q / n) * sum_k [H_q(mu_k)]^2
```
or more precisely, that `H_r` is an "average" of `H_p` and `H_q` evaluated
at subordinated points, with the averaging being an L^2 contraction.

---

## 12. Immediate Next Steps for the Proof Team

1. **Numerical experiment**: For small n (n=3,4), compute the subordination
   functions omega_p, omega_q explicitly for `p boxplus_n q` and verify
   that `G_r(z) = G_p(omega_q(z))`.

2. **Read Voiculescu 1998 carefully**: The proof of the free Stam inequality
   (Proposition 5.8 and surrounding material) should be studied line by line.

3. **Read Shlyakhtenko-Tao Section 4**: The self-contained analytic proof
   of Fisher monotonicity. Identify which steps are infinite-dimensional
   and which can be made finite.

4. **Investigate Gribinski's entropy**: Is `Dis(p)` related to `log disc(p)`?
   If so, is `d/dt log disc(p boxplus_n gamma_t)` related to `Phi_n`?

5. **Check**: Does the finite free K-transform (Marcus) provide subordination?
   I.e., does `K_{p boxplus_n q}` factor through `K_p` and `K_q`?

---

## Sources

- [Voiculescu 1998](https://link.springer.com/article/10.1007/s002220050222)
- [Szarek-Voiculescu 1996](https://link.springer.com/article/10.1007/BF02108815)
- [Shlyakhtenko-Tao 2020](https://arxiv.org/abs/2009.01882)
- [Tao's blog post on fractional free convolution](https://terrytao.wordpress.com/2020/09/07/free-fractional-convolution-powers/)
- [Biane-Speicher 2001](https://www.numdam.org/item/AIHPB_2001__37_5_581_0/)
- [Gribinski 2019](https://arxiv.org/abs/1907.12826)
- [Marcus 2021](https://arxiv.org/abs/2108.07054)
- [MSS 2015](https://arxiv.org/abs/1504.00350)
- [Arizmendi-Perales 2018](https://arxiv.org/abs/1611.06598)
- [Steinerberger 2020](https://arxiv.org/abs/2009.03869)
- [Hall-Ho 2022](https://arxiv.org/abs/2202.09660)
- [Free probabilistic diffusion 2025](https://arxiv.org/abs/2510.22778)
- [Free probability via entropic OT 2023](https://arxiv.org/abs/2309.12196)
- [S-transform in finite free 2024](https://arxiv.org/abs/2408.09337)
- [Finite N precursors 2025](https://arxiv.org/abs/2508.21483)
- [Flow of polynomial roots 2020](https://arxiv.org/abs/2012.09080)
