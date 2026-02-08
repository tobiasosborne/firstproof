# VERIFIER-13 Report: Verification of PROVER-14 Claims

## Summary

Five claims from PROVER-14 were subjected to adversarial verification using symbolic algebra (sympy), exact algebraic reasoning, and intensive numerical testing. Three claims are VALID, one is PARTIALLY VALID, and one is MOSTLY VALID with corrections needed.

---

## Claim 1: Phi_n = 2 * Sm2

**VERDICT: VALID**

### Statement
Phi_n(p) = sum_i H_p(lambda_i)^2 = 2 * Sm2, where Sm2 = sum_{i<j} 1/(lambda_i - lambda_j)^2.

### Verification

**Triple identity (algebraic).** The core identity is:
```
1/((a-b)(a-c)) + 1/((b-a)(b-c)) + 1/((c-a)(c-b)) = 0
```
This was verified symbolically (sympy returns 0) and by manual common-denominator calculation: with D = (a-b)(b-c)(c-a), the numerators sum to -(b-c)-(c-a)-(a-b) = 0.

**Cross-term cancellation.** Expanding Phi_n = sum_i H_i^2 gives:
- Diagonal terms: sum_{i!=j} 1/(lambda_i - lambda_j)^2 = 2*Sm2
- Cross terms: B = 2 * sum_i sum_{j<k, j!=i, k!=i} 1/((lambda_i - lambda_j)(lambda_i - lambda_k))

Each cross term belongs to exactly one unordered triple {i,j,k}. The combinatorial count confirms this: n*(n-1)*(n-2)/2 total cross terms, C(n,3) = n*(n-1)*(n-2)/6 triples, exactly 3 terms per triple. By the triple identity, each triple sums to 0, so B = 0.

**Factor of 2 confirmed.** Numerical tests for n=3,...,7 (500 trials each) give max relative error below 8e-16. Symbolic verification for n=3 and n=4 confirms Phi - 2*Sm2 = 0 exactly.

### Issues in PROVER-14's presentation
PROVER-14's initial proof attempt (lines 9-22 of key_identity.py) is circular. The partial fraction attempt (lines 24-47) has sign confusion. The correct proof via the triple identity (lines 186-217) is sound, though the common-denominator calculation at lines 196-217 contains some false starts before arriving at the right answer.

---

## Claim 2: S2 Additivity under MSS Convolution

**VERDICT: VALID**

### Statement
S2(p boxplus_n q) = S2(p) + S2(q), where S2 = sum_{i<j} (lambda_i - lambda_j)^2.

### Verification

**Exact algebraic proof.** Using S2 = (n-1)*a_1^2 - 2n*a_2 (where a_k are the polynomial coefficients) and the MSS formula:
- c_1 = a_1 + b_1
- c_2 = a_2 + b_2 + ((n-1)/n)*a_1*b_1

We get:
```
S2(r) = (n-1)*(a_1+b_1)^2 - 2n*(a_2 + b_2 + ((n-1)/n)*a_1*b_1)
      = (n-1)*a_1^2 + 2(n-1)*a_1*b_1 + (n-1)*b_1^2 - 2n*a_2 - 2n*b_2 - 2(n-1)*a_1*b_1
      = [(n-1)*a_1^2 - 2n*a_2] + [(n-1)*b_1^2 - 2n*b_2]
      = S2(p) + S2(q)
```

The cross-term 2(n-1)*a_1*b_1 from expanding c_1^2 cancels exactly with -2(n-1)*a_1*b_1 from the c_2 term. This is an exact algebraic identity, not approximate.

**Numerical confirmation.** For n=3,...,6, max relative errors are below 3e-14 (both centered and non-centered polynomials).

### Notes on cumulants
PROVER-14 mentions "S2 = n * k_2" where k_2 is the second finite free cumulant. This is imprecise: for centered polynomials, S2 = n^2 * kappa_2 (not n * kappa_2). For non-centered polynomials, the relationship involves correction terms. However, the S2 additivity result does NOT depend on the cumulant interpretation -- it follows directly from the MSS coefficient formula.

---

## Claim 3: Universal Identity sum_i H_i * lambda_i = n(n-1)/2

**VERDICT: VALID**

### Statement
For any monic degree-n polynomial with simple roots lambda_1, ..., lambda_n:
sum_i H_p(lambda_i) * lambda_i = n(n-1)/2.

### Verification

**Algebraic proof.**
```
sum_{i} H_i * lambda_i = sum_{i!=j} lambda_i / (lambda_i - lambda_j)
                       = sum_{i<j} [lambda_i/(lambda_i-lambda_j) + lambda_j/(lambda_j-lambda_i)]
                       = sum_{i<j} [lambda_i/(lambda_i-lambda_j) - lambda_j/(lambda_i-lambda_j)]
                       = sum_{i<j} (lambda_i - lambda_j)/(lambda_i - lambda_j)
                       = sum_{i<j} 1
                       = C(n,2) = n(n-1)/2
```

This is a completely elementary pairing argument. Confirmed symbolically for n=3 (gives 3) and numerically for n=2,...,10 with errors below 2e-15, including edge cases with nearly-equal roots.

---

## Claim 4: Q = Sm2 * S2 Scale Invariance and Q_r <= max(Q_p, Q_q)

**VERDICT: PARTIALLY VALID**

### Sub-claim 4a: Q is scale-invariant
**VALID.** Under lambda_i -> c*lambda_i: Sm2 -> Sm2/c^2, S2 -> c^2*S2, so Q -> Q. Verified numerically for c in {0.1, 0.5, 2.0, 10.0, -3.0}.

### Sub-claim 4b: Cauchy-Schwarz bound Q >= [n(n-1)/2]^2
**VALID.** By the Cauchy-Schwarz inequality applied to vectors (1/|d_{ij}|) and (|d_{ij}|) over all C(n,2) pairs. Equality requires all pairwise distances equal, which is impossible for n >= 3 distinct real numbers.

### Sub-claim 4c: Q_r <= max(Q_p, Q_q)
**NUMERICALLY SUPPORTED BUT UNPROVED.** Tested with:
- 5000 random trials per n for n=3,4,5: zero violations
- 3000 adversarial trials (extreme gap ratios) per n for n=3,4: zero violations
- Targeted self-convolution tests: Q_r/Q_p < 1 always

However, worst ratios Q_r/max(Q_p,Q_q) approach 1 for n=3 self-convolution (0.9999996), suggesting the bound could be tight.

**IMPORTANT: Q_r <= min(Q_p, Q_q) is FALSE.** Approximately 40% of random trials for n=3 violate this stronger bound.

### Sub-claim 4d: Q_r <= max(Q_p, Q_q) implies the main conjecture
**INVALID.** PROVER-14 initially claims this implication (Part 6, lines 157-177) but then correctly identifies the error: knowing Q_r <= Q_p (WLOG Q_p >= Q_q) gives LHS >= (S2_p+S2_q)/Q_p and RHS >= (S2_p+S2_q)/Q_p, so both sides are bounded below by the same quantity. This is inconclusive.

The correct equivalent reformulation is:
```
Q_r <= Q_p * Q_q * (S2_p + S2_q) / (S2_p * Q_q + S2_q * Q_p)
```
This is the weighted harmonic mean of Q_p, Q_q with weights S2_p, S2_q. It is EQUIVALENT to (not a consequence of) the main conjecture 1/Phi_r >= 1/Phi_p + 1/Phi_q.

---

## Claim 5: n=3 Trigonometric Parametrization

**VERDICT: MOSTLY VALID**

### Sub-claim 5a: Parametrization
**VALID.** For centered cubic x^3 + sigma*x + tau = 0 with sigma < 0 and discriminant Delta = -4*sigma^3 - 27*tau^2 > 0, the roots are:
```
x_k = 2*sqrt(-sigma/3) * cos(phi + 2k*pi/3),  k = 0, 1, 2
```
where cos(3*phi) = -3*sqrt(3)*tau / (2*(-sigma)^{3/2}) and phi in (0, pi/3).

This is the standard Cardano-Vieta trigonometric solution. Verified numerically.

### Sub-claim 5b: G(c) = -9/(16c^6 - 24c^4 + 9c^2 - 1)
**VALID.** Where c = cos(phi) and G(c) = F(phi) = csc^2(phi) + csc^2(phi+pi/3) + csc^2(phi+2pi/3).

Key discovery (verified but not stated clearly by PROVER-14):
```
16c^6 - 24c^4 + 9c^2 - 1 = (4c^3 - 3c)^2 - 1 = cos^2(3phi) - 1 = -sin^2(3phi)
```
Therefore G(c) = 9/sin^2(3*phi).

This follows from the beautiful identity:
```
csc^2(x) + csc^2(x + pi/3) + csc^2(x + 2pi/3) = 9 * csc^2(3x)
```
which itself follows by differentiating the product identity sin(x)*sin(x+pi/3)*sin(x+2pi/3) = sin(3x)/4 twice (via the logarithmic derivative).

### Sub-claim 5c: Q = 3*F(phi)/2
**VALID.** Since Q = Sm2*S2 and Sm2 = F(phi)/(12R^2), S2 = 18R^2, we get Q = F(phi)*18R^2/(12R^2) = 3F(phi)/2 = 27/(2*sin^2(3phi)).

### Sub-claim 5d: Convexity
**VALID but imprecisely stated.** F(phi) = 9*csc^2(3phi) is convex on (0, pi/3), which is immediate since csc^2 is convex on any interval where it is defined.

More importantly, Q as a function of u = cos(3phi) is:
```
Q(u) = 27/(2(1 - u^2))
```
This is convex in u on (-1, 1) since Q''(u) = 27(2+6u^2)/(1-u^2)^3 > 0.

However, **convexity of Q in u does NOT directly yield the conjecture**, because the MSS formula gives cos(3phi_r) as a CONTRACTED (not convex-combination) average of cos(3phi_p) and cos(3phi_q). The contraction factor (a_p^{3/2} + a_q^{3/2})/(a_p+a_q)^{3/2} < 1 pushes cos(3phi_r) toward 0, which pushes Q_r toward 13.5 (the minimum).

The unresolved n=3 conjecture reduces to:
```
F(phi_r) <= (a_p + a_q) / (a_p/F(phi_p) + a_q/F(phi_q))
```
where phi_r is determined by cos(3phi_r) = [a_p^{3/2}*cos(3phi_p) + a_q^{3/2}*cos(3phi_q)] / (a_p+a_q)^{3/2}.

This remains unproved despite the elegant reformulation.

---

## Errors and Gaps Found

1. **Circular reasoning (Claim 1, lines 9-22).** PROVER-14's first proof attempt defines cross_i in terms of H_i^2 and then tries to derive H_i^2 from cross_i, which is circular. The correct proof (via triple identity) appears later.

2. **False implication (Claim 4, Part 6).** The claim "Q_r <= max(Q_p, Q_q) implies the main conjecture" is wrong. PROVER-14 eventually corrects this.

3. **Missing proof for n=3.** Despite extensive analysis, the n=3 case remains UNPROVED. The reduction to F(phi_r) <= weighted_harmonic_mean is correct but the final inequality is not established.

4. **Imprecise cumulant relation (Claim 2).** "S2 = n * k_2" should be "S2 = n^2 * kappa_2" for centered polynomials, where kappa_2 is the second finite free cumulant.

5. **Variable confusion (Claim 5).** PROVER-14 sometimes writes "c = cos(phi)" and sometimes implies "c" relates to cos(3phi). The natural MSS variable is u = cos(3phi), not c = cos(phi). The formula G(c) = -9/(16c^6-24c^4+9c^2-1) is in terms of c = cos(phi), but the MSS addition law operates on tau and sigma, which translates to a transformation of cos(3phi).

---

## Counterexamples

No counterexamples to any of the main claims were found despite intensive adversarial testing. Specifically:
- Phi_n = 2*Sm2: no violations (exact identity)
- S2 additivity: no violations (exact identity)
- sum H_i*lambda_i = n(n-1)/2: no violations (exact identity)
- Q_r <= max(Q_p, Q_q): no violations in ~50,000 trials across n=3,4,5 with both random and adversarial inputs
- Main conjecture 1/Phi_r >= 1/Phi_p + 1/Phi_q: no violations found

The strongest near-violation found for Q_r <= max(Q_p, Q_q) was Q_r/max(Q_p, Q_q) = 0.9999996 for n=3 self-convolution with phi near the boundary (phi close to 0 or pi/3).

---

## Correct Statements of Results

For reference, here are the precise versions of all validated claims:

**Theorem 1.** For distinct real numbers lambda_1, ..., lambda_n, define H_i = sum_{j!=i} 1/(lambda_i - lambda_j). Then sum_i H_i^2 = 2 * sum_{i<j} 1/(lambda_i - lambda_j)^2.

**Theorem 2.** For monic degree-n polynomials p, q and r = p boxplus_n q (MSS convolution), S2(r) = S2(p) + S2(q) where S2(f) = sum_{i<j} (mu_i - mu_j)^2 for roots mu_i of f.

**Theorem 3.** For distinct real numbers lambda_1, ..., lambda_n, sum_i lambda_i * H_i = n(n-1)/2 where H_i = sum_{j!=i} 1/(lambda_i - lambda_j).

**Theorem 4 (partial).** Q = Sm2 * S2 is invariant under affine scaling lambda_i -> c*lambda_i + d. By Cauchy-Schwarz, Q >= [n(n-1)/2]^2 with strict inequality for n >= 3.

**Theorem 5.** For n=3, Q(phi) = 27/(2*sin^2(3phi)) where phi parametrizes the root shape via the Cardano-Vieta trigonometric solution. Equivalently, G(c) = -9/(16c^6-24c^4+9c^2-1) = 9/sin^2(3*arccos(c)) for c = cos(phi) in (1/2, 1).
