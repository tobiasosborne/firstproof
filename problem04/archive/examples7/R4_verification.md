# R_4 Verification Report (VERIFIER-9)

**Date:** 2026-02-08
**Agent:** VERIFIER-9 (adversarial)
**Scripts:** `R4_verify.py`, `R4_verify_part2.py`

---

## Executive Summary

Three claims were verified: the R_4 formula, the constant C_4, and the superadditivity of R_4. Two are fully confirmed. A major error was found in the general C_n conjecture.

| Claim | Verdict | Details |
|-------|---------|---------|
| R_4 formula | **CONFIRMED** | Algebraically verified via sympy |
| C_4 = 1/12 | **CONFIRMED** | Correct: 4/(4^2 * 3) = 1/12 |
| C_n = 2/(n(n-1)) general | **REFUTED** | Correct formula: C_n = 4/(n^2(n-1)) |
| R_4 superadditivity | **SUPPORTED** | 0 violations in 17,000+ adversarial trials |

**Challenge filed:** `ch-3e848e240f7d92e7` (major) against node 1.5.2.3 for wrong C_n general formula.

---

## Part A: Derivation of 1/Phi_4 from First Principles

### Method

For centered degree-4 polynomial `p(x) = x^4 + a_2 x^2 + a_3 x + a_4`:

1. Computed N_4 = Phi_4 * disc as a polynomial in (a_2, a_3, a_4) by numerical fitting with rational reconstruction, confirmed by symbolic computation in sympy.

2. Used the weight-homogeneity constraint: under root scaling by t, we have a_k -> t^k a_k, Phi_4 -> t^{-2} Phi_4, disc -> t^{12} disc, so N_4 is homogeneous of weight 10 with weights (2, 3, 4).

3. This restricts N_4 to exactly 5 monomials with 2a + 3b + 4c = 10.

### Results

**N_4 = Phi_4 * disc:**
```
N_4 = -8*e2^5 - 64*e2^3*e4 - 36*e2^2*e3^2 + 384*e2*e4^2 - 432*e3^2*e4
```

**Discriminant (centered quartic):**
```
disc = 256*e4^3 - 128*e2^2*e4^2 + 144*e2*e3^2*e4 - 27*e3^4 + 16*e2^4*e4 - 4*e2^3*e3^2
```

**Cumulant-coefficient relations (n=4 centered):**
```
a_2 = -3*k2/2
a_3 = -k3/2
a_4 = 3*k2^2/16 - 3*k4/32
```

Substituting into disc and N_4:

```
disc = (27/128) * (32*k2^6 - 32*k2^3*k3^2 - 6*k2^2*k4^2 + 24*k2*k3^2*k4 - 8*k3^4 - k4^3)

N_4 = (81/16) * (4*k2^2 - k4) * (4*k2^3 + k2*k4 - 2*k3^2)
```

**1/Phi_4 = disc / N_4:**
```
1/Phi_4 = (32*k2^6 - 32*k2^3*k3^2 - 6*k2^2*k4^2 + 24*k2*k3^2*k4 - 8*k3^4 - k4^3)
         / (384*k2^5 - 192*k2^2*k3^2 - 24*k2*k4^2 + 48*k3^2*k4)
```

At k3 = k4 = 0: 1/Phi_4 = 32*k2^6 / (384*k2^5) = k2/12, confirming C_4 = 1/12.

### R_4 extraction

```
R_4 = 1/Phi_4 - k2/12
    = (-16*k2^3*k3^2 - 4*k2^2*k4^2 + 20*k2*k3^2*k4 - 8*k3^4 - k4^3)
      / (384*k2^5 - 192*k2^2*k3^2 - 24*k2*k4^2 + 48*k3^2*k4)
```

This is **identical** to the claimed formula (the claimed denominator is `24*(16*k2^5 - ...)` which expands to `384*k2^5 - ...`). **Algebraically verified by sympy: difference = 0.**

### Key factorization

The denominator of 1/Phi_4 factors as:
```
24 * (4*k2^2 - k4) * (4*k2^3 + k2*k4 - 2*k3^2)
```

This factorization is important for understanding the domain structure.

---

## Part B: The C_n Discrepancy

### Error in HANDOFF.md

HANDOFF.md line 75 states: `C_n = 2/(n(n-1))` with values `C_2=1, C_3=2/9, C_4=1/12`.

**Check:**
- 2/(2*1) = 1 = C_2 -- matches
- 2/(3*2) = 1/3 -- but C_3 = 2/9, so **DOES NOT MATCH**
- 2/(4*3) = 1/6 -- but C_4 = 1/12, so **DOES NOT MATCH**

The HANDOFF table lists the correct individual VALUES (C_2=1, C_3=2/9, C_4=1/12) but gives the wrong GENERAL FORMULA.

### Correct formula

Computed C_n empirically for n = 2, 3, 4, 5:

| n | C_n (computed) | C_n = 4/(n^2(n-1)) | C_n = 2/(n(n-1)) | Match correct? | Match wrong? |
|---|---|---|---|---|---|
| 2 | 1 | 1 | 1 | Yes | Yes |
| 3 | 2/9 | 2/9 | 1/3 | Yes | **No** |
| 4 | 1/12 | 1/12 | 1/6 | Yes | **No** |
| 5 | 1/25 | 1/25 | 1/10 | Yes | **No** |

The correct formula is:
```
C_n = 4 / (n^2 * (n-1))
```

Equivalently: `C_n * C(n,2) = 2/n`, or `C_n = (2/n) / C(n,2)`.

**Impact:** The specific claimed values C_4 = 1/12 and C_3 = 2/9 are correct. The error is only in the general formula, which would affect attempts to prove the conjecture for arbitrary n.

---

## Part C: Superadditivity Testing

### 1/Phi_4 superadditivity (via MSS convolution)

| Test category | Trials | Violations |
|---|---|---|
| Random Gaussian | 10,000 | 0 |
| Extreme ratio (100:1) | 500 | 0 |
| Extreme ratio (1000:1) | 500 | 0 |
| Extreme ratio (1:100) | 500 | 0 |
| Extreme ratio (1:1000) | 500 | 0 |
| Large k3/k4 vs k2 | 500 | 0 |
| Sign (1,1) | 500 | 0 |
| Sign (1,-1) | 500 | 0 |
| Sign (-1,1) | 500 | 0 |
| Sign (-1,-1) | 500 | 0 |
| **Total** | **~14,000** | **0** |

Minimum gap: 2.60e-03 (strictly positive in all valid trials).

### Direct R_4 superadditivity (using cumulant vectors)

7,036 valid tests, 0 violations, minimum gap 2.52e-03.

Additional check: all 1,000 tested summed cumulant vectors correspond to valid polynomials with real simple roots (guaranteed by MSS real-rootedness theorem).

---

## Part D: Domain Analysis

### Valid domain
The domain of valid (k2, k3, k4) is defined by:
- k2 > 0 (positive variance)
- disc > 0 (four distinct real roots)

In kappa coordinates: disc = (27/128) * (32*k2^6 - 32*k2^3*k3^2 - 6*k2^2*k4^2 + 24*k2*k3^2*k4 - 8*k3^4 - k4^3) > 0.

### Denominator positivity
The denominator `24*(4*k2^2 - k4)*(4*k2^3 + k2*k4 - 2*k3^2)` was tested on 9,170 valid samples. It is **always positive** on the valid domain (0 negative cases found).

This means 1/Phi_4 is well-defined and finite on the entire valid domain.

### Sign of R_4
R_4 is **always negative** on the valid domain (0 positive cases in 10,000 samples). This makes sense: R_4 is the negative correction from higher-order cumulants, reducing 1/Phi_4 below its "pure Gaussian" value of k2/12.

### Boundary behavior
- As k3, k4 -> 0: R_4 -> 0 (quadratically, as expected)
- The discriminant boundary (disc = 0) is not approached in the R_4 formula since both numerator and denominator of 1/Phi_4 vanish there proportionally.

---

## Part E: Structural Observations

### N_4 factorization
```
N_4 = (81/16) * (4*k2^2 - k4) * (4*k2^3 + k2*k4 - 2*k3^2)
```

The factor `(4*k2^2 - k4)` is always positive on the valid domain. The factor `(4*k2^3 + k2*k4 - 2*k3^2)` is also always positive.

### R_4 numerator
The numerator `-16*k2^3*k3^2 - 4*k2^2*k4^2 + 20*k2*k3^2*k4 - 8*k3^4 - k4^3` does NOT factor over the rationals (confirmed by sympy). This makes analytical proofs of superadditivity more challenging.

### Superadditivity mechanism
Since R_4 is always negative, superadditivity `R_4(k_p + k_q) >= R_4(k_p) + R_4(k_q)` means the correction is **less negative** for the sum than for the parts. Equivalently, |R_4| is subadditive. This is the analog of the n=3 case where R_3 = -2*k3^2/(27*k2^2) and the Cauchy-Schwarz inequality ensures subadditivity of |R_3|.

---

## Challenges Filed

| ID | Node | Severity | Description |
|---|---|---|---|
| ch-3e848e240f7d92e7 | 1.5.2.3 | major | General formula C_n = 2/(n(n-1)) is wrong. Correct: C_n = 4/(n^2(n-1)). |

---

## Files

- `/home/tobiasosborne/Projects/af-tests/examples7/R4_verify.py` -- Main verification script
- `/home/tobiasosborne/Projects/af-tests/examples7/R4_verify_part2.py` -- Algebraic verification and C_n pattern
- `/home/tobiasosborne/Projects/af-tests/examples7/R4_verification.md` -- This report
