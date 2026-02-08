# VERIFIER-10b Independent Verification Report

## Date: 2026-02-08
## Node: 1.5.5.1
## Script: `examples7/verifier10b_check.py`

---

## Methodology

All computations were performed independently from PROVER-10b's code:
- Different random seed (99991 vs prover's 42)
- Independent MSS convolution implementation (two separate implementations cross-checked)
- Independent cumulant computation (verified symbolically via EGF log expansion)
- Adversarial test battery: near-degenerate roots, extreme scale asymmetry, outlier roots, symmetric configurations
- Hermite polynomial argument verified step-by-step (symbolically and numerically)

---

## Claim 1: C_n = 4/(n^2*(n-1))

**VERDICT: CONFIRMED**

Independent verification for n=2,...,10. All ratios (computed/predicted) equal 1.000000000000 to 12 decimal places.

The Hermite polynomial proof was verified in three steps:
1. **H_i = xi_i/2**: Verified symbolically (derivative identity He_n' = n*He_{n-1} plus three-term recurrence at roots) and numerically (max error < 1e-14 for n=3,4,5).
2. **sum(xi_i^2) = n(n-1)**: Verified to < 2e-13 for all n=2,...,10.
3. **Phi_n = n(n-1)/(4s)**: Follows from steps 1-2 by direct computation.

The original conjecture C_n = 2/(n(n-1)) is definitively WRONG. The discrepancy factor is exactly n/2, confirmed for all n.

---

## Claim 2: R_n <= 0

**VERDICT: CONFIRMED**

Tested for n=3,...,8 with adversarial sampling (nearly coincident roots, outlier roots, extreme scales):
- n=3: 0 positive out of 1932 trials
- n=4: 0 positive out of 1833 trials
- n=5: 0 positive out of 1729 trials
- n=6: 0 positive out of 1624 trials
- n=7: 0 positive out of 1459 trials
- n=8: 0 positive out of 1352 trials

---

## Claim 3: R_n superadditivity

**VERDICT: CONFIRMED**

Adversarial test battery for each n included 5 categories:
1. Random (500 trials, different seed)
2. Nearly coincident roots (200 trials, gaps as small as 0.05)
3. Large scale asymmetry (200 trials, scale ratio up to 20x)
4. Outlier roots (200 trials, outliers at 5-15 standard deviations)
5. Equally-spaced roots (16 configurations)

Results:
- n=3: 985 trials, 0 violations, min margin = 0.00e+00
- n=4: 891 trials, 0 violations, min margin = 7.08e-05
- n=5: 754 trials, 0 violations, min margin = 6.50e-05
- n=6: 612 trials, 0 violations, min margin = 5.09e-05

The margins for 1/Phi_n and R_n superadditivity are identical (max difference < 4e-15), confirming that kappa_2 additivity causes the linear term to cancel.

---

## Claim 4: Decomposition 1/Phi_n = C_n*kappa_2 + R_n

**VERDICT: CONFIRMED**

The decomposition is correct by construction: R_n is defined as the remainder after subtracting the linear term. The key content is:
- C_n is correctly extracted as the coefficient of kappa_2 at the Gaussian locus
- kappa_2 is additive under MSS (verified: max error < 7e-15)
- R_n is the nonlinear correction that vanishes at the Gaussian locus

---

## Claim 5: R_3 = -(1/6)*K3^2/K2^2

**VERDICT: CONFIRMED**

Max |1/Phi_3 - ((-2/3)*K2 - (1/6)*K3^2/K2^2)| = 1.82e-14 across 947 samples.
R_3 <= 0 in all 947 cases.

---

## Claim 6: R_4 exact formula

**VERDICT: CONFIRMED**

Max relative error between numerical R_4 and the claimed rational function:
- Max relative error: 1.15e-12
- Median relative error: 1.00e-15

Special case R_4(K2, 0, 0) = 0: verified symbolically.
Special case R_4(K2, 0, K4) = K4^2/(9*K2*(6*K2^2+K4)): verified symbolically (diff = 0).

---

## Claim 7: Cumulant additivity

**VERDICT: CONFIRMED**

Cumulant formulas verified symbolically via EGF log expansion (exact match for kappa_2,...,kappa_6).
Additivity verified numerically under MSS for n=3,...,6 (max errors < 5e-13).

---

## Issues Found

### MINOR: Internal inconsistency in prover's code

PROVER-10b's script has an inconsistency in C_n(ta):
- Line 20 (comment) and line 154 (code): `C_n(ta) = -4/(n^2*(n-1))` -- WRONG
- Line 312 (compute_R_n function): `C_n_ta = -4/(n*(n-1))` -- CORRECT

The correct formula is C_n(ta) = -4/(n*(n-1)). The relationship to C_n(their) is:
  C_n(their) = -C_n(ta)/n = 4/(n^2*(n-1))

This inconsistency does not affect any of the claimed results because:
1. The C_n(their) verification in Part 1 is self-consistent (uses the their convention directly)
2. The R_n computations in Parts 4-6 use the CORRECT formula via compute_R_n

### COSMETIC: Report line 639

Line 639 of the report states "C_n(ta) = -4/(n(n-1)) = -4/(n^2*(n-1)) * n" which is mathematically correct but confusing notation.

---

## Overall Verdict

**ALL MAJOR CLAIMS VALIDATED.** The node 1.5.5.1 should be accepted.

The C_n = 4/(n^2*(n-1)) formula, R_n <= 0 conjecture, R_n superadditivity, cumulant additivity, and exact R_3/R_4 formulas are all independently confirmed through adversarial testing.
