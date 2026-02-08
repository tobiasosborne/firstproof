# R_4 Superadditivity: Proof Results

**Author:** PROVER-11
**Date:** 2026-02-08
**Status:** PARTIAL PROOF (k3=0 case complete; full case numerically verified)

---

## 1. Main Results

### Result 1: k3=0 Case (PROVED ANALYTICALLY)

**Theorem.** For s, t > 0 and u in (-4s^2, 4s^2), v in (-4t^2, 4t^2) with u + v in (-4(s+t)^2, 4(s+t)^2):

R_4(s+t, 0, u+v) >= R_4(s, 0, u) + R_4(t, 0, v)

where R_4(K2, 0, K4) = -K4^2 / (24 * K2 * (4*K2^2 - K4)).

**Proof.** Let A = 4s^2 - u > 0, B = 4t^2 - v > 0, C = 4(s+t)^2 - (u+v) > 0.

The inequality is equivalent to (subadditivity of -R_4):

  (u+v)^2 / ((s+t)*C) <= u^2/(s*A) + v^2/(t*B)      ... (*)

**Step 1 (Key Lemma): s*A + t*B <= (s+t)*C.**

Compute:
```
(s+t)*C - s*A - t*B
= (s+t)*[4(s+t)^2 - u - v] - s*(4s^2 - u) - t*(4t^2 - v)
= 4(s+t)^3 - (s+t)(u+v) - 4s^3 + su - 4t^3 + tv
= 4[(s+t)^3 - s^3 - t^3] - (s+t)(u+v) + su + tv
= 4*[3s^2*t + 3st^2] - su - tu - sv - tv + su + tv
= 12st(s+t) - tu - sv
```

Since u < 4s^2 and v < 4t^2:
```
12st(s+t) - tu - sv > 12st(s+t) - 4s^2*t - 4st^2 = 8st(s+t) > 0.  QED
```

**Step 2 (Cauchy-Schwarz):** With a_1 = sqrt(s*A), b_1 = u/sqrt(s*A), a_2 = sqrt(t*B), b_2 = v/sqrt(t*B), the Cauchy-Schwarz inequality gives:

```
(a_1*b_1 + a_2*b_2)^2 <= (a_1^2 + a_2^2)*(b_1^2 + b_2^2)
(u + v)^2 <= (s*A + t*B) * (u^2/(s*A) + v^2/(t*B))
```

**Combining Steps 1 and 2:**

```
(u+v)^2/((s+t)*C) <= (u+v)^2/(s*A+t*B)    [by Step 1: s*A+t*B <= (s+t)*C]
                   <= u^2/(s*A) + v^2/(t*B) [by Step 2: Cauchy-Schwarz]
```

This proves (*) and hence the theorem. QED

**Numerical verification:** 0 violations in 500,000 valid random trials.

---

### Result 2: Full Case (NUMERICALLY VERIFIED)

**Conjecture (verified numerically).** For all cumulant vectors (k2_p, k3_p, k4_p) and (k2_q, k3_q, k4_q) in the domain (both denominator factors positive), the superadditivity gap

```
R_4(k2_p + k2_q, k3_p + k3_q, k4_p + k4_q) - R_4(k2_p, k3_p, k4_p) - R_4(k2_q, k3_q, k4_q) >= 0
```

**Numerical verification:** 0 violations in 447,000 valid random trials.

---

### Result 3: Partial Fraction Decomposition (VERIFIED)

The function -R_4 admits the partial fraction decomposition in K4:

```
-R_4 = (4*K2^3 - K3^2) / (6*(4*K2^2 - K4))
     + K3^2*(4*K2^3 - K3^2) / (6*K2^2*(4*K2^3 + K2*K4 - 2*K3^2))
     - K4 / (24*K2)
     - (2*K2^3 + K3^2) / (12*K2^2)
```

**Verified symbolically** (sympy: difference = 0).

Key structural fact: `4*K2^3 - K3^2 > 0` on the entire domain. Proof: from D2 > 0 we get `2*K3^2 < 4*K2^3 + K2*K4`, and from D1 > 0 we get `K4 < 4*K2^2`, so `K3^2 < 2*K2^3 + K2*K4/2 < 2*K2^3 + 2*K2^3 = 4*K2^3`.

---

## 2. Failed Approaches (Full Case)

### 2.1 Joint Concavity of f(u,v)
R_4 = K2 * f(K3/K2^{3/2}, K4/K2^2). If f were concave, the perspective argument would give superadditivity. But the Hessian of f is **indefinite** at points like (u, v) = (0, 1) and (1, 0). FAILED.

### 2.2 Term-by-Term Subadditivity
The partial fraction decomposition has 4 terms. The term -K4/(24*K2) is NOT individually subadditive (it is a ratio, and (K4_p+K4_q)/(K2_p+K2_q) is a weighted average, not a sum). So we cannot prove subadditivity term by term. FAILED.

### 2.3 Discriminant Structure
The numerator of -R_4, viewed as a quadratic in K3^2, has discriminant `16*(K2^2 - 2K4)*(4K2^2 - K4)^2`. This factors nicely but does not directly yield a proof.

---

## 3. Promising Directions (Not Yet Attempted or Incomplete)

1. **SDP-based SOS Certificate:** The gap numerator is a specific polynomial. Use DSOS/SDSOS tools to find an explicit sum-of-squares certificate for non-negativity.

2. **Perturbation from k3=0:** The k3=0 case is proved. For small k3, can we show the gap remains positive by a continuity/monotonicity argument? The gap at k3=0 is strictly positive (except at degenerate points), so small perturbations should preserve it.

3. **Schur-Convexity:** Fix scaled cumulant ratios and study the gap as a function of (k2_p, k2_q). The symmetry under (p <-> q) is suggestive of Schur-convexity methods.

4. **Generalized Cauchy-Schwarz with Vector Weights:** Instead of scalar weights (s*A, t*B), use a matrix-valued weight that captures all three cumulants simultaneously. For n=3, a 2x2 positive-definite matrix worked. For n=4, a 3x3 or 4x4 matrix may be needed.

---

## 4. Scripts

- `R4_prover11.py` - Initial setup and formula verification
- `R4_prover11_part3.py` - k3=0 Cauchy-Schwarz proof development
- `R4_prover11_part5.py` - Clean k3=0 numerical verification
- `R4_prover11_part7.py` - Discriminant factorization
- `R4_prover11_part8.py` - Partial fraction decomposition + full numerical
- `R4_prover11_final.py` - Consolidated verification of all results

---

## 5. Summary Table

| Component | Status |
|-----------|--------|
| R_4 formula | VERIFIED (prior work) |
| k3=0 analytical proof | **COMPLETE** (Cauchy-Schwarz + Key Lemma) |
| Full case numerical | **VERIFIED** (0 violations in 447K trials) |
| Partial fraction decomposition | **VERIFIED** (sympy) |
| Joint concavity of f(u,v) | **FAILED** (Hessian indefinite) |
| Term-by-term subadditivity | **FAILED** (K4/K2 term not subadditive) |
| Full analytical proof | **OPEN** |
