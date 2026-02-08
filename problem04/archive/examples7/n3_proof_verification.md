# Verification Report: n=3 Superadditivity Proof

## Status: VERIFIED WITH CAVEATS

**Verifier:** VERIFIER-8 (adversarial)
**Date:** 2026-02-08
**Scripts:** `verifier8_n3_check.py`, `verifier8_deep_checks.py`, `verifier8_equality_check.py`, `verifier8_matrix_connection.py`, `verifier8_convention_check.py`, `verifier8_convention_detail.py`

---

## Executive Summary

The claimed n=3 proof of finite free Fisher information superadditivity is **CORRECT**. Every step has been verified both algebraically (via sympy symbolic computation) and numerically (over thousands of random trials). No errors were found.

Three minor caveats are noted, none of which affect the mathematical correctness of the proof.

---

## Step 1: Formula verification

**Claim:** For n=3, `1/Phi_3 = (2/9)*k2 - (2/27)*k3^2/k2^2`

**Derivation verified:**
1. N_3 := Phi_3 * disc = 18 * e2^2 -- **VERIFIED SYMBOLICALLY** (sympy: difference = 0, using centered roots a, b, -(a+b))
2. disc = -4*e2^3 - 27*e3^2 -- standard discriminant formula, verified numerically
3. 1/Phi_3 = disc / N_3 = (-4*e2^3 - 27*e3^2) / (18*e2^2)
4. Cumulant substitution: e2 = -k2, e3 = (2/9)*k3 -- **VERIFIED** (kappa_2 = -e2 and kappa_3 = (9/2)*e3 for centered n=3, confirmed over 30 trials)
5. Resulting formula: (2/9)*k2 - (2/27)*k3^2/k2^2 -- **VERIFIED SYMBOLICALLY** (sympy: difference = 0)

**Numerical verification:** 50 random centered cubic polynomials, max error 7.11e-15. PASS.

**Additional finding:** The formula holds for ALL polynomials in P_3, not just centered ones, because k2 and k3 are both translation-invariant (proved symbolically). The centering step in the proof is therefore unnecessary, though harmless.

**Verdict: CORRECT.**

---

## Step 2: Centering argument

**Claim:** WLOG we may assume both p and q are centered (kappa_1 = 0).

**Verification:**

1. **Phi_n is translation-invariant:** H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j). Under roots -> roots + c, the differences lambda_i - lambda_j are unchanged. Therefore Phi_n(p(.−c)) = Phi_n(p). **VERIFIED ALGEBRAICALLY** (exact cancellation). Numerically confirmed; small discrepancies at large shifts (c = 10000) are purely floating-point (relative error ~2e-10).

2. **MSS convolution commutes with translation:** (p(.−a)) ⊞_n (q(.−b)) = (p ⊞_n q)(.−(a+b)). **VERIFIED** numerically over 30 trials (max root difference < 1e-6).

3. **k2, k3 are translation-invariant:** Under roots -> roots + t, the e_k change but the cumulants kappa_2 and kappa_3 do not. **VERIFIED SYMBOLICALLY** via sympy (k2_shifted - k2_original = 0, k3_shifted - k3_original = 0).

4. **Centering argument:** Since Phi_n, k2, k3 are all translation-invariant, and MSS commutes with translation, the inequality for general (p, q) is equivalent to the inequality for centered (p_c, q_c).

**Note:** Since the formula 1/Phi_3 = (2/9)*k2 - (2/27)*k3^2/k2^2 holds for ALL polynomials (not just centered), the centering step is logically unnecessary. This is a minor expositional point, not an error.

**Verdict: CORRECT (centering is valid but unnecessary).**

---

## Step 3: Domain verification

**Claim:** k2 > 0 for all p in P_3 (simple roots).

**Algebraic proof:** For centered n=3 with roots a, b, -(a+b):
- k2 = a^2 + ab + b^2 = (a + b/2)^2 + (3/4)*b^2 >= 0
- Equality iff a = b = 0, which means all three roots equal zero (not simple)
- Therefore k2 > 0 for all p in P_3. **QED.**

Since k2 is translation-invariant, this holds for all (not just centered) polynomials.

**Numerical verification:** 1000 random polynomials, minimum k2 observed = 0.000032, all positive. PASS.

**Domain of (k2, k3):** For simple roots, disc > 0 implies k2^3 - k3^2/3 > 0, i.e., |k3| < sqrt(3) * k2^(3/2). Verified over 500 trials. PASS.

**Verdict: CORRECT.**

---

## Step 4: Superadditivity reduction

**Claim:** The superadditivity inequality reduces to:
`k3_p^2/k2_p^2 + k3_q^2/k2_q^2 >= (k3_p+k3_q)^2/(k2_p+k2_q)^2`

**Verification of the reduction:**
- 1/Phi_3(r) = (2/9)*k2_r - (2/27)*k3_r^2/k2_r^2
- By cumulant additivity: k2_r = k2_p + k2_q, k3_r = k3_p + k3_q
- The (2/9)*k2 terms cancel exactly (linearity)
- Remaining: -(2/27) * [(k3_p+k3_q)^2/(k2_p+k2_q)^2 - k3_p^2/k2_p^2 - k3_q^2/k2_q^2]
- This is >= 0 iff (k3_p+k3_q)^2/(k2_p+k2_q)^2 <= k3_p^2/k2_p^2 + k3_q^2/k2_q^2

**Cumulant additivity:** Verified **SYMBOLICALLY** via sympy for n=3 using the exact MSS coefficient formula. All three cumulants (k1, k2, k3) satisfy k_i(p ⊞ q) = k_i(p) + k_i(q) exactly. This is an algebraic identity, not a numerical approximation.

**Key inequality verification:** The substitution x = k3_p/k2_p, y = k3_q/k2_q, s = k2_p/(k2_p+k2_q), t = 1-s transforms the inequality to (s*x + t*y)^2 <= x^2 + y^2, which follows from Jensen's inequality (convexity of z^2) and the trivial bound s*x^2 + t*y^2 <= x^2 + y^2.

**Numerical verification:** 10,000 random (k2_p, k2_q, k3_p, k3_q) with k2_p, k2_q > 0. Zero violations. PASS.

**Verdict: CORRECT.**

---

## Step 5: Matrix positive-definiteness argument

**Claim:** The matrix M = [[t^3(2s+t), -s^2*t^2], [-s^2*t^2, s^3(2t+s)]] with s=k2_p, t=k2_q is positive definite, which implies the key inequality.

**Verification:**

1. **Connection to the inequality:** The quadratic form v^T * Q * v where Q = M/(s^2*t^2*(s+t)^2) equals exactly k3_p^2/s^2 + k3_q^2/t^2 - (k3_p+k3_q)^2/(s+t)^2. **VERIFIED** numerically (ratio = 1.0000000000 over 20 trials).

2. **det(M) = 2*s^3*t^3*(s+t)^2:** Derived algebraically from det(Q) = 2/(s*t*(s+t)^2) and M = s^2*t^2*(s+t)^2 * Q. **VERIFIED** numerically over 100 trials.

3. **trace(M) = t^3(2s+t) + s^3(2t+s) > 0:** Both diagonal entries are positive for s, t > 0.

4. **M positive definite:** trace > 0 and det > 0 implies both eigenvalues positive. **VERIFIED** numerically (all eigenvalues > 0 over 20 trials).

**Assessment:** The matrix argument is CORRECT but unnecessarily complex. The Jensen/Cauchy-Schwarz proof (Step 4) is simpler and more transparent:
- (sx+ty)^2 <= s*x^2 + t*y^2 (Jensen, since z^2 is convex and s+t=1)
- s*x^2 + t*y^2 <= x^2 + y^2 (trivial: x^2+y^2 - sx^2 - ty^2 = tx^2 + sy^2 >= 0)

Both proofs are valid. The Jensen argument has the advantage of being a two-line proof.

**Verdict: CORRECT (but unnecessarily complex).**

---

## Step 6: Edge cases

### k3 = 0 (equally-spaced roots)
- Roots {-1, 0, 1}: k2 = 1, k3 = 0, 1/Phi_3 = 2/9 = 0.222... Formula gives 0.222... MATCH.
- Convolution of two equally-spaced: {-1,0,1} ⊞ {-2,0,2} gives EXACT equality (margin < 1e-15).

### k2_p >> k2_q
- Tested ratios 10, 100, 1000, 10000. All give RHS > LHS (inequality satisfied). No issues.

### Near boundary (nearly-repeated roots)
- Tested eps from 0.1 down to 1e-5. Formula and direct computation agree to ~1e-17. As roots approach coalescence, 1/Phi_3 -> 0 (disc -> 0 while N_3 stays bounded). No singularity issues.

### Adversarial counterexample search
- 5000 trials with various distributions (random, near-coalescence, extreme spread, skew). Zero violations. Minimum margin: 7.2e-7.
- 2000 additional direct trials (no cumulants). Zero violations. Minimum margin: 8.5e-5.

**Verdict: All edge cases pass.**

---

## Equality characterization

**Equality holds iff k3_p = k3_q = 0**, i.e., both p and q have equally-spaced roots.

**Proof:** The Jensen step gives equality iff x = y (i.e., k3_p/k2_p = k3_q/k2_q). The second step gives equality iff x = y = 0. Combined: equality iff k3_p = k3_q = 0.

For centered n=3: k3 = 0 iff e3 = abc = 0 iff one root is zero iff roots are {-a, 0, a} (equally spaced).

**Verified numerically:**
- Equally-spaced pairs: margin = 0 (exact equality)
- Non-equally-spaced: margin > 0 (strict inequality, tested 5 pairs)
- Same k3/k2 ratio but nonzero: strict inequality (tested 3 configurations)

This matches the HANDOFF claim.

---

## Caveats

**C1 (Minor, expositional):** The matrix M argument is correct but the Jensen/Cauchy-Schwarz proof is simpler. The final writeup should prefer the two-line Jensen argument.

**C2 (Minor, completeness):** The formula N_3 = Phi_3 * disc = 18*e2^2 is verified symbolically but could benefit from a standalone derivation in the proof document. (The symbolic verification using sympy constitutes a proof, but a human-readable derivation would be clearer.)

**C3 (Minor, dependency):** Cumulant additivity under MSS is verified symbolically for n=3 from the explicit MSS coefficient formula. This is a self-contained algebraic identity for degree 3. The proof does not depend on any general theory of finite free cumulants -- it can be checked by direct computation for n=3 specifically.

---

## Conclusion

The n=3 superadditivity proof is **mathematically correct** with a complete, gap-free proof chain:

1. Express 1/Phi_3 as (2/9)*k2 - (2/27)*k3^2/k2^2 (verified symbolically)
2. Cumulant additivity under MSS: k_i(p ⊞ q) = k_i(p) + k_i(q) (verified symbolically)
3. Linear term cancels; correction term superadditivity reduces to (sx+ty)^2 <= x^2+y^2
4. This follows from Jensen's inequality (two lines)

The only external dependency is MSS Theorem 4.4 (real-rootedness preservation), which is well-established in the literature.

**No challenges filed.** The proof is correct.

---

## Appendix: Convention check

The proof uses the Arizmendi-Perales (AP) convention for finite free cumulants:
- kappa_1 = tilde_a_1
- kappa_2 = -n * (tilde_a_2 - tilde_a_1^2)
- kappa_3 = (n^2/2) * (tilde_a_3 - 3*tilde_a_2*tilde_a_1 + 2*tilde_a_1^3)

where tilde_a_k = (-1)^k * a_k / C(n,k) and a_k are the polynomial coefficients.

An alternative convention without the n-prefactors (k2' = -(tilde_a_2 - tilde_a_1^2), k3' = (1/2)*(tilde_a_3 - ...)) also gives additive cumulants, since k2_AP = n * k2' and k3_AP = n^2 * k3'. Any scalar multiple of an additive quantity is additive.

The formula `1/Phi_3 = (2/9)*k2 - (2/27)*k3^2/k2^2` uses the AP convention. With the alternative convention, the formula would be `1/Phi_3 = (2/3)*k2' - (2/3)*k3'^2/k2'^2`. Both are correct; what matters is internal consistency, which is maintained throughout the proof.

**Verified:** The AP convention gives Phi_2 * k2 = 1 exactly for n=2 (10 trials, ratio = 1.0000000000). This serves as a sanity check on the convention.
