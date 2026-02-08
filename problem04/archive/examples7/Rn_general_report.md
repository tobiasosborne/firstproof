# C_n and R_n General Report

## Author: PROVER-10b
## Date: 2026-02-08
## Script: `examples7/Rn_prover10b.py`

---

## 1. Main Result: C_n Formula

**Theorem (verified numerically for n=2,...,10):**

```
C_n = 4 / (n^2 * (n-1))
```

in the convention where `their_kappa_2 = -n * ta_2` (the convention used in the
original conjecture statement where C_2 = 1, C_3 = 2/9, C_4 = 1/12, C_5 = 1/25).

Equivalently, in the Arizmendi-Perales normalized coefficient convention (ta_k):

```
C_n(ta) = -4 / (n * (n-1))
```

### Verification Table

| n  | C_n(their) predicted | C_n(their) computed | Match |
|----|---------------------|--------------------:|:-----:|
| 2  | 1                   | 1.0000000000        |  Y    |
| 3  | 2/9                 | 0.2222222222        |  Y    |
| 4  | 1/12                | 0.0833333333        |  Y    |
| 5  | 1/25                | 0.0400000000        |  Y    |
| 6  | 1/45                | 0.0222222222        |  Y    |
| 7  | 2/147               | 0.0136054422        |  Y    |
| 8  | 1/112               | 0.0089285714        |  Y    |
| 9  | 1/162               | 0.0061728395        |  Y    |
| 10 | 1/225               | 0.0044444444        |  Y    |

All ratios are 1.0000000000 to 10 decimal places.

### Correction to Original Conjecture

The original conjecture stated C_n = 2/(n(n-1)). The actual formula is:

```
C_n = 4/(n^2(n-1)) = (2/n) * 2/(n(n-1))
```

The original conjecture was off by a factor of 2/n. This discrepancy arose from
an incorrect normalization of the second cumulant kappa_2.

---

## 2. Proof of C_n (Hermite Polynomial Argument)

### Setup

At the Gaussian locus (kappa_3 = kappa_4 = ... = kappa_n = 0), the normalized
coefficients satisfy ta_{2m} = (2m-1)!! * K2^m and ta_{2m+1} = 0, where K2 = kappa_2 = ta_2.

### Proof Steps

1. **Hermite connection:** At the Gaussian locus, the polynomial is
   p_n(x) = s^{n/2} He_n(x/sqrt(s)) where s = -K2 > 0 and He_n is the
   probabilist's Hermite polynomial.

2. **Roots:** lambda_i = sqrt(s) * xi_i where xi_i are the roots of He_n.

3. **Key identity:** H_i = xi_i / 2 for all Hermite roots.
   - Proof: He_n'(x) = n * He_{n-1}(x) (derivative identity)
   - He_n''(x) = n(n-1) * He_{n-2}(x) (second derivative)
   - Three-term recurrence at root xi_i: 0 = xi_i * He_{n-1}(xi_i) - (n-1) * He_{n-2}(xi_i)
   - Therefore He_{n-2}(xi_i) = xi_i * He_{n-1}(xi_i) / (n-1)
   - H_i = He_n''(xi_i) / (2 * He_n'(xi_i)) = n(n-1) * He_{n-2}(xi_i) / (2n * He_{n-1}(xi_i)) = xi_i/2

4. **Sum of squares:** sum_i xi_i^2 = n(n-1) by Vieta's formula for He_n.
   - He_n(x) = x^n - C(n,2)*x^{n-2} + ..., so sum xi_i = 0, sum_{i<j} xi_i*xi_j = -C(n,2)
   - sum xi_i^2 = (sum xi_i)^2 - 2*sum_{i<j} xi_i*xi_j = 0 + 2*C(n,2) = n(n-1)

5. **Phi_n computation:**
   Phi_n = sum_i H_i^2 = (1/s) * sum_i (xi_i/2)^2 = n(n-1)/(4s)

6. **C_n extraction:**
   1/Phi_n = 4s/(n(n-1)) = -4*K2/(n(n-1))
   So C_n(ta) = -4/(n(n-1)).
   In "their" convention: C_n = 4/(n^2(n-1)).

---

## 3. Finite Free Cumulants

### Definition

The finite free cumulants are defined via the EGF logarithm:
- T(x) = 1 + sum_{k>=2} ta_k * x^k/k!  (EGF of normalized coefficients)
- log T(x) = sum_{k>=2} kappa_k * x^k/k!

Under MSS convolution: T_r = T_p * T_q, so log is additive:
kappa_k(p conv q) = kappa_k(p) + kappa_k(q).

### Explicit Formulas (centered, kappa_1 = 0)

| k | kappa_k in terms of ta_j |
|---|--------------------------|
| 2 | ta_2 |
| 3 | ta_3 |
| 4 | ta_4 - 3*ta_2^2 |
| 5 | ta_5 - 10*ta_2*ta_3 |
| 6 | ta_6 - 15*ta_2*ta_4 - 10*ta_3^2 + 30*ta_2^3 |

### Additivity Verification

All cumulants verified additive under MSS convolution to machine precision:

| n | kappa_2 error | kappa_3 error | kappa_4 error | kappa_5 error | kappa_6 error |
|---|:---:|:---:|:---:|:---:|:---:|
| 3 | 1e-14 | 4e-14 | -- | -- | -- |
| 4 | 2e-14 | 3e-14 | 3e-13 | -- | -- |
| 5 | 1e-14 | 2e-14 | 2e-13 | 1e-12 | -- |
| 6 | 7e-15 | 2e-14 | 7e-14 | 4e-13 | 2e-12 |

---

## 4. R_n Structure

### n=2: R_2 = 0

Trivial: 1/Phi_2 = kappa_2 exactly.

### n=3: R_3 = -(1/6) * K3^2 / K2^2

- Always <= 0 (since K2^2 > 0 and K3^2 >= 0)
- Superadditive by Cauchy-Schwarz / Jensen's inequality
- **Proof:** Need K3^2/K2^2 + K3'^2/K2'^2 >= (K3+K3')^2/(K2+K2')^2.
  Set x = K3/K2, y = K3'/K2', t = K2/(K2+K2'). Then RHS = (xt + y(1-t))^2.
  By convexity of z^2: (xt+y(1-t))^2 <= t*x^2+(1-t)*y^2 <= x^2+y^2 = LHS. QED.

### n=4: R_4 (exact rational function)

```
R_4 = (-54*K2^3*K3^2 + 6*K2^2*K4^2 - 45*K2*K3^2*K4 + 27*K3^4 - K4^3)
      / (9*(6*K2^2 + K4)*(6*K2^3 - K2*K4 + 3*K3^2))
```

Special cases:
- At K4=0: R_4 = -K3^2*(2*K2^3 - K3^2) / (6*K2^2*(2*K2^3 + K3^2))
- At K3=0: R_4 = K4^2 / (9*K2*(6*K2^2 + K4))
- At K3=K4=0: R_4 = 0 (consistent with C_n extraction)

### n=5: R_5 (numerical approximation)

Leading terms (polynomial regression, relative residual 0.38):
```
R_5 ~ -0.220 * K3^2/K2^2 + 0.013 * K4/K2 + 0.017 * K3*K5/K2^3
     + 0.034 * K4^2/K2^3 - 0.226 * K3^2*K4/K2^4 + 0.162 * K3^4/K2^5
```

Note: R_5 is a rational function (not polynomial) in cumulants, so the polynomial
fit has limited accuracy. The exact formula requires symbolic computation of
Phi_5 * disc_5, which has weight 18.

### n=6: R_6 (numerical approximation)

Leading terms (2-feature fit, relative residual 0.77):
```
R_6 ~ -0.100 * K3^2/K2^2 + 0.006 * K4/K2 + (higher order)
```

---

## 5. Superadditivity Results

### Full 1/Phi_n superadditivity

| n | Trials | Violations | Min margin |
|---|-------:|-----------:|:-----------|
| 3 | 3371   | 0          | 1.89e-04   |
| 4 | 2830   | 0          | 4.11e-03   |
| 5 | 2243   | 0          | 3.93e-03   |
| 6 | 1698   | 0          | 4.49e-03   |

### R_n superadditivity

| n | Trials | Violations | Min margin |
|---|-------:|-----------:|:-----------|
| 3 | 3371   | 0          | 1.89e-04   |
| 4 | 2830   | 0          | 4.11e-03   |
| 5 | 2243   | 0          | 3.93e-03   |
| 6 | 1698   | 0          | 4.49e-03   |

Note: The margins for R_n and 1/Phi_n superadditivity are identical, which
is expected since kappa_2 is additive and thus the C_n*kappa_2 term cancels.

### R_n Sign

R_n <= 0 in ALL observed cases for ALL n:

| n | Positive | Negative | Total |
|---|:--------:|:--------:|:-----:|
| 3 | 0        | 442      | 442   |
| 4 | 0        | 753      | 753   |
| 5 | 0        | 1304     | 1304  |
| 6 | 0        | 1467     | 1467  |

**Conjecture:** R_n <= 0 for all n and all valid (simple-root) polynomials.
This means the Gaussian locus MAXIMIZES 1/Phi_n among polynomials with the
same kappa_2 value.

---

## 6. Structural Patterns

1. **Weight-2 homogeneity:** R_n(t^2*K2, t^3*K3, ...) = t^2 * R_n(K2, K3, ...)

2. **Gaussian locus is extremal:** R_n = 0 at kappa_3 = ... = kappa_n = 0,
   and R_n < 0 away from this locus.

3. **Dominant term is K3^2/K2^2:** For all n >= 3, the leading contribution
   to R_n involves K3^2/K2^2 with a negative coefficient.

4. **Superadditivity reduces to Cauchy-Schwarz-type inequality:** For n=3,
   superadditivity of R_3 is exactly the Cauchy-Schwarz inequality. For
   higher n, the rational function structure makes the proof more complex
   but the underlying mechanism appears similar.

5. **Rational function structure:** 1/Phi_n = disc_n / (Phi_n * disc_n)
   where Phi_n * disc_n is a polynomial in the elementary symmetric functions.
   The rational function structure persists in cumulant coordinates.

---

## 7. Open Problems

1. **Prove R_n <= 0 for general n.** This would establish that the Gaussian
   locus is the maximum of 1/Phi_n at fixed kappa_2.

2. **Prove R_n is superadditive for general n.** This is the KEY step that
   would complete the proof of the main superadditivity conjecture.

3. **Find exact R_n formulas for n >= 5.** This requires symbolic computation
   of Phi_n * disc_n (weight n(n-1)-2 polynomial in e_2, ..., e_n).

4. **Understand the connection to free probability.** The cumulants are the
   LOG of the EGF of normalized coefficients, matching the classical
   moment-cumulant relation in free probability. The W_k(i,j) = C(k,i)
   structure (binomial convolution) connects to exponential generating functions.

---

## 8. Convention Summary

| Quantity | "ta" convention | "their" convention |
|----------|:---------------:|:------------------:|
| kappa_2  | ta_2 (< 0)      | -n*ta_2 (> 0)      |
| C_n      | -4/(n(n-1))     | 4/(n^2(n-1))       |
| 1/Phi_n  | C_n(ta)*ta_2 + R_n | C_n(their)*their_k2 + R_n |
| R_n      | same in both conventions | same |
