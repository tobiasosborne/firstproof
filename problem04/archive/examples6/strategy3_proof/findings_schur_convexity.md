# Findings: Schur Convexity / Majorization Approach to Fisher Superadditivity

**Agent:** prover (Schur convexity investigation)
**Date:** 2026-02-08
**Scripts:** `investigate_schur_convexity.py`, `investigate_schur_followup.py`, `investigate_schur_ostrowski_detail.py`

---

## Executive Summary

The Schur convexity / majorization approach was investigated as a potential proof strategy for Fisher superadditivity:

```
1/Phi_n(p boxplus_n q) >= 1/Phi_n(p) + 1/Phi_n(q)
```

**Bottom line:** The approach yields valuable structural insights but does NOT provide a direct proof. The core obstacle is that superadditivity is a statement about THREE objects (p, q, r) while majorization-based Schur arguments relate only TWO objects at a time.

**Most significant discovery:** The centered root vector of r = p boxplus q ALWAYS majorizes the centered root vectors of both p and q individually. This appears to be a theorem (100% of 3500+ trials for n = 2 through 8). Combined with other properties, this gives 1/Phi(r) >= max(1/Phi(p), 1/Phi(q)), which is HALF the desired inequality.

---

## Detailed Findings

### Finding 1: Centered Root Majorization (the strongest result)

**Statement:** For the MSS convolution r = p boxplus_n q:
```
centered_r majorizes centered_p    (100% of all trials)
centered_r majorizes centered_q    (100% of all trials)
```
where centered_f = roots(f) - mean(roots(f)).

**Tested:** n = 2, 3, 4, 5, 6, 7, 8 with 500 random trials each. Zero exceptions.

**Significance:** This means the roots of r are "more spread out" than those of either p or q individually, in the precise sense of majorization. This is likely a consequence of the random matrix interpretation (A + UBU* averages over unitary conjugation, spreading eigenvalues).

**Limitation:** This only gives Phi(r) <= min(Phi(p), Phi(q)) (verified at 100%), which yields 1/Phi(r) >= max(1/Phi(p), 1/Phi(q)). But the target inequality 1/Phi(r) >= 1/Phi(p) + 1/Phi(q) requires controlling BOTH terms simultaneously.

---

### Finding 2: Phi_n is Convex in the Root Vector

**Statement:** Phi_n(t*x + (1-t)*y) <= t*Phi_n(x) + (1-t)*Phi_n(y) for all root vectors x, y.

**Tested:** n = 3, 4, 5 with 2000 trials each. Zero violations.

**Note:** The MSS convolution r is NOT a convex combination of p and q in root space, so Jensen's inequality does not directly apply. However, this convexity is an interesting structural fact.

---

### Finding 3: Phi_n is NOT Schur-convex in Roots

Phi_n is NEITHER Schur-convex NOR Schur-concave as a function of the centered root vector:
- n=3: ~98.5% of pairs show Phi(x) < Phi(y) when x majorizes y (looks Schur-concave), but ~1.5% go the other way.
- n=4, 5: similar pattern with more mixed behavior.

This means we cannot combine Finding 1 with a Schur argument on roots.

---

### Finding 4: Phi_n Schur-convexity in Gaps is PARTIAL

**n = 3:** Phi_n is genuinely Schur-convex in the gap vector (0 violations in all tests, Schur-Ostrowski condition verified).

**n >= 4:** Phi_n is NOT strictly Schur-convex in the gap vector. There exist counterexamples where g1 majorizes g2 but Phi(g1) < Phi(g2). However:
- The violations are small (ratio ~ 0.9998-0.9999).
- They occur primarily when two ADJACENT gaps are nearly equal.
- When using T-transforms on the MAX vs MIN gap, the violations disappear.

**Important nuance:** 1/Phi_n is Schur-concave in gaps for n = 3 (verified), but NOT for n >= 4 (small counterexamples exist). So a direct Schur-concavity argument for 1/Phi in gaps is only viable for n = 3.

---

### Finding 5: Gap Vector Majorization Under MSS Convolution

**Raw gap vectors:** The gap vectors of r, p, q have different total sums, so ordinary majorization does not apply.

**Normalized gap vectors (sum to 1):**
- ng_p majorizes ng_r in ~84% of cases (n = 3-6).
- ng_r majorizes ng_p in ~8% of cases.
- Neither in ~0% of cases (they are almost always comparable).

This means the normalized gap structure of r is MORE UNIFORM than that of p. The MSS convolution "equalizes" the gaps. This is consistent with the fact that equally spaced roots give the tightest case for the Fisher inequality.

---

### Finding 6: Variance (p_2) Super-additivity

**Statement:** For centered roots, sum(roots_r^2) >= sum(roots_p^2) + sum(roots_q^2).

**Tested:** n = 3, 4, 5, 6 with 500 trials each. Ratio always >= 1.000000 (exact to machine precision).

This is a known consequence of the MSS convolution: e_2(r) = e_2(p) + e_2(q) - (n-1)/n * e_1(p)*e_1(q) plus correction terms, which simplifies under centering.

**Higher power sums also super-additive:**
- p_4(r) >= p_4(p) + p_4(q): min ratio ~ 1.016 (n=3), 1.031 (n=4), 1.041 (n=5)
- p_6(r) >= p_6(p) + p_6(q): min ratio ~ 1.024 (n=3), 1.038 (n=4), 1.066 (n=5)

These all follow from the centered root majorization (Finding 1), since power sums are Schur-convex.

---

### Finding 7: Shape Factor Decomposition

By scale invariance, Phi(c * roots) = Phi(roots) / c^2. So:
```
1/Phi(roots) = var(roots) / (Phi(roots) * var(roots)) = var / sf
```
where sf = Phi * var is the dimensionless "shape factor" depending only on gap ratios.

The superadditivity becomes:
```
var(r) / sf(r) >= var(p) / sf(p) + var(q) / sf(q)
```

Since var(r) >= var(p) + var(q) (Finding 6), this would follow from:
```
(var_p + var_q) / sf(r) >= var_p / sf_p + var_q / sf_q
```
i.e., sf(r) <= harmonic mean of sf(p) and sf(q) weighted by variances.

**Observation:** sf ranges from about 3.0 (uniform gaps) to > 100 (very unequal gaps) for n = 3-5. The shape factor of r tends to be SMALLER than those of p and q (because r has more uniform gaps), which is the right direction. But quantifying this precisely is the challenge.

---

### Finding 8: Derivative Identity Verification

The MSS derivative identity (p boxplus q)' = n * (p^{(1)} boxplus_{n-1} q^{(1)}) was verified numerically. The ratio Phi(crit_r) / Phi(roots_r) varies between 0.18 and 0.54. There is no simple universal relationship that would support an induction argument directly.

---

### Finding 9: Pythagorean Property

For equally spaced roots, gaps_r = sqrt(gaps_p^2 + gaps_q^2) (element-wise). For non-equally-spaced roots, this is only approximately true (within 10% in ~22% of cases). The Pythagorean property is exact only for the variance (Finding 6).

---

## Assessment: Why Schur/Majorization Does Not Suffice

The fundamental obstacle is structural:

1. **Majorization gives one-sided bounds:** centered_r majorizes centered_p gives Phi(r) <= Phi(p), i.e., 1/Phi(r) >= 1/Phi(p). But we need 1/Phi(r) >= 1/Phi(p) + 1/Phi(q), which is the SUM of two terms.

2. **No "additive" majorization:** There is no relation of the form "roots_r majorizes roots_p + roots_q" that would convert the superadditivity to a single Schur bound.

3. **Phi is not Schur-convex/concave in roots:** So even the one-sided bound from Finding 1 requires a different argument (not pure Schur theory).

4. **Gap majorization has wrong direction:** The normalized gaps of r are MORE uniform (LESS majorized) than those of p, q. Combined with partial Schur-convexity of Phi in gaps, this gives Phi(r) <= Phi(p) by a different route, but still only one-sided.

5. **Two-polynomial interaction:** The key difficulty is that superadditivity couples p and q. Majorization compares r to p OR r to q, but never handles both simultaneously in a way that produces the sum 1/Phi(p) + 1/Phi(q).

---

## What MIGHT Work (Honest Assessment)

### Possible path via shape factor (Finding 7):
If one could prove sf(r) <= min(sf(p), sf(q)), then:
```
1/Phi(r) = var(r)/sf(r) >= (var(p)+var(q))/min(sf(p),sf(q))
         >= var(p)/sf(p) + var(q)/sf(q) = 1/Phi(p) + 1/Phi(q)
```
This would work if sf(r) <= min(sf(p), sf(q)). Numerically, sf(r) is often much smaller than both sf(p) and sf(q), but I have NOT verified this universally.

### Possible path via normalized gap Schur bound (n=3 only):
For n = 3, Phi is genuinely Schur-convex in gaps, and the normalized gaps of r are majorized by those of p. This, combined with the variance super-additivity, could give a complete proof for n = 3.

### Possible path via log-convexity or information geometry:
The natural metric on the space of polynomials with distinct roots may have curvature properties that interact well with the MSS convolution. This is unexplored.

---

## Literature Search Summary

1. **Carlen (1991):** Proved superadditivity of classical Fisher information using scaling properties and conditional entropy. The proof technique is specific to the Gaussian setting and does not obviously transfer to the polynomial case.

2. **Voiculescu (1993-1997):** Introduced free entropy and free Fisher information. The free entropy power inequality was proved using restricted Minkowski sums and the Brascamp-Lieb-Luttinger rearrangement inequality. This is a fundamentally infinite-dimensional argument.

3. **Marcus-Spielman-Srivastava (2015, 2022):** The finite free convolution preserves real-rootedness and provides strong root bounds. The barrier method used in the Kadison-Singer resolution is related but focused on interlacing, not majorization of the kind needed here.

4. **Marcus (2021):** Survey on finite free probability. Does not address entropy or Fisher information directly. The "finite free cumulants are additive" property is the key structural fact.

5. **Gribinski (2019):** Defines a notion of entropy for polynomial roots (arXiv:1907.12826) and proves it increases under finite free addition. This is conceptually related but the entropy notion is different from Fisher information.

6. **Schur-Horn inequalities for hyperbolic polynomials (2025):** Recent work uses Garding's concavity for hyperbolicity cones. The MSS convolution preserves hyperbolicity, so this could be relevant, but the connection to Fisher information is unclear.

**No existing result in the literature directly addresses finite free Fisher information superadditivity.**

---

## Recommendations

1. **Do NOT pursue pure Schur/majorization approach** for the general proof. It fundamentally cannot handle the additive structure of the inequality.

2. **DO use the centered root majorization** as a building block. It provides the crucial bound Phi(r) <= min(Phi(p), Phi(q)) which is necessary (though not sufficient).

3. **Investigate the shape factor bound** sf(r) <= min(sf(p), sf(q)) as a separate conjecture. If true, it closes the inequality.

4. **For n = 3 specifically,** the Schur approach has a better chance due to genuine Schur-convexity of Phi in gaps. A specialized n = 3 proof via this route may be achievable.

5. **Consider information-theoretic approaches** modeled on Carlen's classical proof, adapted to the finite polynomial setting.

---

## Appendix: Definitions

- **Majorization:** x majorizes y (x >- y) if sum of k largest entries of x >= sum of k largest entries of y for all k, with equal total sums.
- **Weak majorization:** Same but without the equal total sum constraint.
- **Schur-convex:** f is Schur-convex if x >- y implies f(x) >= f(y).
- **Schur-concave:** f is Schur-concave if x >- y implies f(x) <= f(y).
- **Schur-Ostrowski criterion:** f is Schur-convex iff (x_i - x_j)(df/dx_i - df/dx_j) >= 0 for all i, j (for differentiable f).
- **Phi_n(p):** Sum of H_p(lambda_i)^2 where H_p(lambda_i) = sum_{j != i} 1/(lambda_i - lambda_j).
- **MSS convolution:** (p boxplus_n q)(x) = E_U[det(xI - A - UBU*)] where U is Haar unitary.
