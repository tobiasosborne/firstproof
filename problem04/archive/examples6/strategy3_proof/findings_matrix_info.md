# Findings: Matrix / Information-Theoretic Approach to Fisher Superadditivity

**Node:** 1.10.3
**Owner:** matrix-prover (prover)
**Date:** 2026-02-08
**Status:** CONJECTURE REFUTED (for n >= 4)

---

## Executive Summary

The conjecture `1/Phi_n(p boxplus_n q) >= 1/Phi_n(p) + 1/Phi_n(q)` is
**FALSE** for n >= 4, even for centered polynomials. It holds exactly for n=2
(trivially) and n=3 (centered), but fails with explicit counterexamples
at n=4. The free probability motivation (Voiculescu's exact additivity
`1/Phi*(mu boxplus nu) = 1/Phi*(mu) + 1/Phi*(nu)`) does NOT produce a
one-sided finite inequality.

---

## 1. Definitions Recap

- `Phi_n(p) = sum_i H_p(lambda_i)^2` where `H_p(lambda_i) = sum_{j != i} 1/(lambda_i - lambda_j)`
- `r = p boxplus_n q` via Marcus-Spielman-Srivastava finite free convolution
- `boxplus_n` defined by: `hat_e_k(r) = sum_{j=0}^k hat_e_j(p) * hat_e_{k-j}(q)` where `hat_e_k = e_k / C(n,k)`

## 2. Random Matrix Connection

If A, B are n x n Hermitian with characteristic polynomials p, q, and U is
Haar-distributed on U(n), then:

    E_U[det(xI - (A + UBU*))] = p boxplus_n q

This is the Marcus-Spielman-Srivastava (2015) identity. The finite free
convolution IS the expected characteristic polynomial.

## 3. Free Probability Connection

Voiculescu's free Fisher information for a compactly-supported measure mu:

    Phi*(mu) = integral |H_mu(x)|^2 dmu(x)

where H_mu is the Hilbert transform. For a discrete (empirical) measure:

    Phi*(mu_p) = Phi_n(p) / n^3

**Voiculescu's theorem** (for free convolution in the n -> infinity limit):

    1/Phi*(mu boxplus nu) = 1/Phi*(mu) + 1/Phi*(nu)

This is an EXACT EQUALITY for free convolution of freely independent
variables. The finite free convolution boxplus_n converges to the free
convolution as n -> infinity.

**The conjecture** was that the finite version satisfies the inequality:

    1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q)

i.e., `superadditivity` of `1/Phi_n` under finite free convolution.

## 4. Matrix Formulation

### Cauchy Matrix Expression
Let C be the off-diagonal Cauchy matrix: `C_{ij} = 1/(lambda_i - lambda_j)` for `i != j`, `C_{ii} = 0`.

Then:
- `H_i = (C @ 1)_i` where `1 = (1,...,1)^T`
- `Phi_n = ||C @ 1||^2 = 1^T C^T C 1`
- C is antisymmetric: `C^T = -C`, so `C^T C = -C^2`
- Therefore: `Phi_n = -sum_{i,j} (C^2)_{ij}`

This is verified numerically for all tested n.

### Translation Invariance
Phi_n is translation invariant: shifting all roots by the same constant c
does not change Phi_n. This is immediate from the definition since
`lambda_i - lambda_j` is invariant.

## 5. Exact Results

### n=2

For p with roots a, b and q with roots c, d:

    1/Phi_n(r) - 1/Phi_n(p) - 1/Phi_n(q) = (a+b)(c+d)/2

This is:
- **Positive** when both means have the same sign
- **Zero** when at least one polynomial is centered (sum of roots = 0)
- **Negative** when means have opposite signs

**Verified exactly.** The n=2 case gives equality for centered polynomials.

### n=3

**Centered symmetric** `p = {-a, 0, a}`, `q = {-b, 0, b}`:

    Phi_n(p) = 9/(2a^2)
    r has roots {-sqrt(a^2+b^2), 0, sqrt(a^2+b^2)}
    Phi_n(r) = 9/(2(a^2+b^2))

    1/Phi_n(r) = 2(a^2+b^2)/9 = 2a^2/9 + 2b^2/9 = 1/Phi_n(p) + 1/Phi_n(q)

**EXACT EQUALITY** for all a, b > 0. Verified numerically.

**Centered asymmetric** (roots sum to 0, but not symmetric):

316 valid random tests. **ALL have gap > 0** (strictly positive).
- min gap: 7.28e-03
- max gap: 1.78e+01
- mean gap: 2.29e+00
- failures: 0

**The superadditivity inequality HOLDS for n=3 centered polynomials.**

### n=4 (COUNTEREXAMPLES FOUND)

Even for centered polynomials with well-separated roots:

**Counterexample 1** (centered symmetric):
- p = {-5, -1, 1, 5}, q = {-5, -1, 1, 5}
- gap = -4.49e-01 (NEGATIVE)

**Counterexample 2** (centered symmetric):
- p = {-3, -1, 1, 3}, q = {-3, -1, 1, 3}
- gap = -9.37e-02 (NEGATIVE)

**Counterexample 3** (centered, well-separated integers):
- p = {-35, 0, 12, 23}, q = {-19, -6, -4, 29}
- gap = -2.16e+01 (NEGATIVE)

**Large-scale test** (n=4 centered, well-separated):
- pass=262, fail=56, skip=31
- min gap = -4.36e+01

**The superadditivity inequality FAILS for n=4.**

### n >= 5

Even more failures proportionally. The conjecture is false.

## 6. Convexity Analysis

### Phi_n as function of root vector

Phi_n appears to be **convex** along linear interpolations of root vectors
(20/20 random tests confirmed). This makes sense: it is a sum of squares
of rational functions, composed with specific linear operations.

### 1/Phi_n

1/Phi_n is **neither convex nor concave** in general along linear root
interpolations (only 3/20 convex, 11/20 concave in random tests).

### Jensen's Inequality Approach

Since r = E_U[chi_{A+UBU*}] and Phi_n appears convex in roots:

    Phi_n(r) <= E_U[Phi_n(A + UBU*)]   (by Jensen, if Phi_n is convex)

Monte Carlo confirms this (`Phi_n(r) < E_U[Phi_n(A+UBU*)]` for n=3,4).

But this gives the WRONG direction for the superadditivity inequality
(it says Phi_n(r) is SMALL, hence 1/Phi_n(r) is LARGE, which would
support the conjecture). The problem is that E_U[Phi_n(A+UBU*)] does
NOT equal Phi_n(A) + Phi_n(B), so the Jensen approach does not close.

### 1/Phi_n in Polynomial Coefficients (n=2)

For n=2: `1/Phi_n = (s^2 - 4p)/2` where s = sum of roots, p = product.
This is linear in p (product) and convex in s (sum). So 1/Phi_n is
convex in the coefficient space for n=2. However, this does not
generalize usefully.

## 7. Literature Search Results

### Voiculescu's Free Fisher Information
- Voiculescu (1993), "The analogues of entropy and Fisher's information
  measure in free probability theory. I", Comm. Math. Phys. 155, 71-92.
- Defines Phi*(X) for noncommutative random variables.
- Additivity: for X, Y freely independent, `Phi*(X,Y) = Phi*(X) + Phi*(Y)`.
- Free Cramer-Rao inequality established.

### Superadditivity of Classical Fisher Information
- Carlen (1991), "Superadditivity of Fisher's information and logarithmic
  Sobolev inequalities", J. Funct. Anal.
- The CLASSICAL Fisher information satisfies superadditivity under
  classical convolution (not free convolution).

### Marcus-Spielman-Srivastava Finite Free Probability
- Marcus, Spielman, Srivastava (2015), "Finite free convolutions of
  polynomials", Probab. Theory Related Fields.
- Marcus (2021), "Polynomial convolutions and (finite) free probability",
  arXiv:2108.07054.
- Defines finite free convolution via hat-elementary symmetric functions.
- Proves real-rootedness preservation.
- Develops finite analogues of R-transform and Cauchy/Stieltjes transforms.

### Missing: Finite Free Fisher Information
No literature was found specifically addressing:
- A "finite free Fisher information" functional
- Superadditivity inequalities for Phi_n under finite free convolution
- Finite analogues of Voiculescu's exact additivity

This appears to be a novel direction that has not been explored in the
literature.

## 8. Why the Conjecture Fails

### Root cause: finite free convolution is NOT free convolution

Voiculescu's exact additivity holds for the FREE convolution of FREELY
INDEPENDENT variables in the large-n limit. The finite free convolution
boxplus_n is a POLYNOMIAL operation that only approximates free convolution
as n -> infinity.

At finite n, boxplus_n introduces correction terms of order 1/n that can
go in either direction. The exact n=2 formula

    gap = (sum_p)(sum_q) / 2

shows that even at the simplest level, the "correction" depends on
the means and can be negative.

### Structural issue: Phi_n depends on higher-order root statistics

Phi_n involves sums of squares of Cauchy-type sums over roots. The finite
free convolution preserves the first few hat-elementary symmetric functions
exactly, but the effect on Phi_n (which involves ALL pairwise root
differences) is not controlled by a simple inequality.

## 9. What IS True

1. **n=2 centered:** Exact equality `1/Phi_2(r) = 1/Phi_2(p) + 1/Phi_2(q)`.

2. **n=3 centered:** Superadditivity holds with strict inequality for
   asymmetric roots. This is the one case where the conjecture is
   non-trivially true.

3. **n=3 centered symmetric:** Exact equality (proven analytically).

4. **Large n (asymptotic):** The gap should converge to zero by
   Voiculescu's theorem (the finite free convolution converges to
   free convolution).

5. **Same-sign means:** The inequality tends to hold when the means of
   the root sets have the same sign (both positive or both negative),
   even without centering.

## 10. Recommendation

**REFUTE node 1.10.3.** The matrix trace / info-theoretic approach
does not yield a proof because the underlying conjecture is false for
n >= 4. The node should be refuted with the n=4 counterexamples.

**Possible salvage:** If the conjecture is MODIFIED to apply only to n=3,
or if a different functional (not Phi_n) is used, there might be a viable
approach. But the current formulation of the superadditivity inequality
is false.

---

## Files

- `investigate_matrix_approach.py` — Full investigation (Parts 1-12)
- `investigate_matrix_approach_v2.py` — Focused on well-separated roots (Parts A-I)
- `findings_matrix_info.md` — This document
