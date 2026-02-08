# Proof Attempt: <h, alpha> >= 0 via Herglotz Convexity

## Node: 1.7.2
## Author: herglotz-prover
## Date: 2026-02-08
## Status: PROOF NOT FOUND -- Herglotz approach definitively blocked

---

## 1. Setup

Let p, q be monic real-rooted polynomials of degree n with simple roots. Let r = p boxplus_n q (MSS finite free convolution). Define:

- Roots: nu_1 < ... < nu_n (of r), lambda_1 < ... < lambda_n (of p), mu_1 < ... < mu_n (of q)
- Subordination: omega_1 determined by G_r(z) = G_p(omega_1(z)), with omega_1(nu_k) = lambda_k
- From the implicit function theorem: omega_1'(nu_k) = 1 for all k
- h_k = H_r(nu_k) = sum_{j != k} 1/(nu_k - nu_j)
- u_k = H_p(lambda_k) = sum_{j != k} 1/(lambda_k - lambda_j)
- alpha_k = omega_1''(nu_k)/2 = u_k - h_k

**Target:** sum_k h_k * alpha_k >= 0.

**Equivalently:** sum_k H_r(nu_k) * H_p(lambda_k) >= sum_k H_r(nu_k)^2 = Phi_n(r).

---

## 2. The Herglotz Approach

### 2.1 The Idea

Define phi(z) = omega_1(z) - z. Then:
- phi(nu_k) = lambda_k - nu_k
- phi'(nu_k) = 0 (since omega_1'(nu_k) = 1)
- phi''(nu_k) = 2 * alpha_k

If omega_1 were a "strict Herglotz function" meaning Im(omega_1(z)) >= Im(z) for z in C^+, then phi would map C^+ to closure(C^+), and the Nevanlinna representation theorem would give:

    phi(z) = c_0 + sum_{j=1}^{n-1} c_j / (d_j - z),   with c_j > 0

This would force phi'(z) = sum c_j/(d_j - z)^2 > 0 for all real z away from poles, which would give alpha_k = phi''(nu_k)/2 having a specific sign structure that could yield <h, alpha> >= 0.

### 2.2 Why This Approach Fails

**The fundamental obstruction:** omega_1 is NOT a Herglotz function in the standard sense.

Through extensive numerical computation (5 test cases, multiple complex points):

1. **omega_1 maps C^+ to C^+:** Im(omega_1(z)) > 0 when Im(z) > 0 (verified numerically).

2. **But Im(omega_1(z)) < Im(z) sometimes:** The minimum ratio Im(omega_1(z))/Im(z) ranges from 0.10 to 0.50, confirming Im(omega_1(z)) can be substantially less than Im(z). Specifically, we found points where Im(phi(z)) = Im(omega_1(z)) - Im(z) < 0.

3. **omega_1 is a BRANCH of an algebraic function**, not a single-valued rational function. It is defined by the degree-n algebraic equation F(z,w) = r(z)*p'(w) - r'(z)*p(w) = 0. For n >= 3, this is irreducible and omega_1 is not globally rational.

4. **The Nevanlinna contradiction:** If omega_1 had the standard Nevanlinna form (z + constant + sum of positive simple poles), then phi' > 0 everywhere on R, contradicting phi'(nu_k) = 0. Since we proved algebraically that omega_1'(nu_k) = 1 (from the implicit function theorem applied to F(z,w) = 0), the Nevanlinna form cannot hold.

**Conclusion:** The Herglotz approach is **definitively blocked**. The subordination function omega_1 in the finite MSS setting is an algebraic function, not a rational Nevanlinna function, and phi = omega_1 - id does not have the monotonicity properties needed.

---

## 3. Alternative Reformulations Discovered

### 3.1 Cauchy Matrix Decomposition

**Identity (verified algebraically and numerically):**

    <h, alpha> = sum_{k < l} (h_k - h_l)(phi_k - phi_l) / [(nu_k - nu_l)(lambda_k - lambda_l)]

where phi_k = nu_k - lambda_k.

**Properties:**
- The weights 1/[(nu_k - nu_l)(lambda_k - lambda_l)] are always POSITIVE (since both root sequences are sorted in the same order).
- Individual terms (h_k - h_l)(phi_k - phi_l) can be negative (7.7% of tested pairs).
- But the weighted sum is ALWAYS positive (verified in 500+ MSS trials).
- For arbitrary (non-MSS) sorted point pairs, the sum is negative 71% of the time, confirming the MSS structure is essential.

### 3.2 Interpolation Monotonicity

Define x_k(t) = (1-t)*nu_k + t*lambda_k for t in [0,1].

**Observation 1:** Phi(t) = sum_k H(x_k(t))^2 is ALWAYS monotone increasing along this interpolation path (0/500 violations). This gives the weaker result Phi_r <= Phi_p.

**Observation 2:** F(t) = <h_r, h_t> where h_t = H-values at x(t) satisfies F(t) >= F(0) = ||h_r||^2 for all t in [0,1] (0/500 violations). In particular F(1) = <h_r, u> >= ||h_r||^2, which IS the statement <h, alpha> >= 0. However, F is not always monotone -- F'(t) can be negative at some intermediate t values (4/1000 cases), but F never dips below F(0).

### 3.3 Schur-concavity and Majorization

**Observation 3:** The Fisher information Phi(x) = sum H(x_k)^2 is SCHUR-CONCAVE (tested: Schur-convexity fails in 500/500 cases; Schur-concavity direction confirmed).

**Observation 4:** The centered r-roots ALWAYS majorize the centered p-roots (500/500 tests), reflecting the "root-spreading" property of free convolution. Combined with Schur-concavity, this gives Phi_r <= Phi_p, but this is weaker than <h, alpha> >= 0.

### 3.4 Universal Identities

- sum_k alpha_k = 0 (antisymmetry)
- sum_k x_k * H(x_k) = n(n-1)/2 for ANY n distinct points (universal; independent of the specific point set)
- sum_k x_k^2 * H(x_k) = (n-1) * sum x_k (also universal)

These identities are too weak to prove <h, alpha> >= 0.

---

## 4. What Would Be Needed for a Proof

### 4.1 Via the Cauchy Matrix Decomposition (Section 3.1)

One would need to show that the specific structure of MSS convolution forces
sum_{k<l} (h_k-h_l)(phi_k-phi_l) / [(nu_k-nu_l)(lambda_k-lambda_l)] >= 0.

This is a weighted correlation between the "H-value increments" (h_k - h_l) and the "shift increments" (phi_k - phi_l). The MSS convolution creates a specific correlation structure between these two sequences that keeps the sum positive. Proving this likely requires understanding the fine structure of how the MSS coefficient formula transforms root configurations.

### 4.2 Via the Interpolation (Section 3.2)

The observation F(t) >= F(0) for all t in [0,1] is equivalent to <h_r, h_t> >= ||h_r||^2 for all t. This is stronger than <h, alpha> >= 0 (which only needs t = 1). A proof could proceed by showing:

1. At t = 0: F(0) = ||h_r||^2.
2. F is continuous on [0,1].
3. If F(t_0) = ||h_r||^2 for some t_0 > 0, then F'(t_0) > 0 (the derivative is strictly positive whenever F touches the baseline).

This "barrier argument" would prevent F from crossing below ||h_r||^2. It requires understanding the relationship between F(t) and F'(t) at points where F is at its minimum, which involves the second-order structure of the subordination.

### 4.3 Via Schur-concavity + Majorization (Section 3.3)

This gives only Phi_r <= Phi_p (weaker). To upgrade this to <h, alpha> >= 0, one would need a refined version of the Schur-concavity argument that takes into account the CORRELATION between h and u (not just their norms).

---

## 5. What the n = 2 Case Reveals

For n = 2, the proof is explicit:
- <h, alpha> = 2(R - s)/(R^2 * s) where s = gap of p, R = gap of r = sqrt(s^2 + t^2), t = gap of q.
- Since R > s (because t > 0), we have <h, alpha> > 0 with equality iff q is trivial.

This uses the explicit formula for n = 2 where omega_1 is a Mobius transformation. For n >= 3, omega_1 is genuinely algebraic and no such closed form exists.

---

## 6. Conclusion

**The Herglotz convexity approach to proving <h, alpha> >= 0 FAILS.** The subordination function omega_1 in the finite MSS setting is an algebraic function (branch of an n-th degree algebraic curve), not a Nevanlinna-class rational function. While omega_1 does map C^+ to C^+, it does NOT satisfy Im(omega_1(z)) >= Im(z), which means phi = omega_1 - id does NOT have the monotonicity structure needed for the Herglotz argument.

**New structural insights discovered:**
- A Cauchy matrix decomposition of <h, alpha> as a weighted sum of concordance terms
- Interpolation F(t) = <h_r, h_t> satisfies F(t) >= F(0) (numerically, stronger than needed)
- Phi(t) = ||h_t||^2 is monotone increasing along the interpolation
- Phi is Schur-concave and the MSS convolution enforces majorization

**Gaps remaining:** The inequality <h, alpha> >= 0 is NOT a consequence of any single general principle (not Herglotz, not Schur-convexity, not majorization alone). It requires the specific MSS structure. The most promising path is the Cauchy matrix decomposition combined with properties of the MSS coefficient formula, or the interpolation barrier argument. Both require new ideas beyond the scope of this investigation.
