# VERIFIER-11: Adversarial Verification Report

**Target:** PROVER-11's k3=0 case of R_4 superadditivity
**Node:** 1.5.2.6.1
**Date:** 2026-02-08
**Verdict:** VALID

---

## 1. Summary

PROVER-11 claims an analytical proof that for s, t > 0 and u in (-4s^2, 4s^2), v in (-4t^2, 4t^2):

R_4(s+t, 0, u+v) >= R_4(s, 0, u) + R_4(t, 0, v)

where R_4(K2, 0, K4) = -K4^2 / (24 * K2 * (4*K2^2 - K4)).

The proof proceeds in two steps:
1. Key Lemma: sA + tB <= (s+t)C (where A = 4s^2-u, B = 4t^2-v, C = 4(s+t)^2-(u+v))
2. Cauchy-Schwarz inequality application

**I have verified this proof is CORRECT. No gaps, errors, or missing steps were found.**

---

## 2. Detailed Verification

### Check 0: Formula correctness
The R_4 formula at k3=0 was verified symbolically using SymPy. Starting from the full R_4 formula:
```
R_4 = (-16*K2^3*K3^2 - 4*K2^2*K4^2 + 20*K2*K3^2*K4 - 8*K3^4 - K4^3) /
      (24*(4*K2^2-K4)*(4*K2^3+K2*K4-2*K3^2))
```
Setting K3=0 and simplifying gives:
```
-K4^2*(4K2^2+K4) / (24*K2*(4K2^2-K4)*(4K2^2+K4)) = -K4^2 / (24*K2*(4K2^2-K4))
```
**PASS.** Formula is correct.

### Check 1: Key Lemma algebra
(s+t)C - sA - tB was expanded symbolically. Result: 12*s^2*t + 12*s*t^2 - s*v - t*u = 12st(s+t) - tu - sv. **PASS.** Algebra is correct.

### Check 2: Lower bound
From u < 4s^2 (and t > 0): tu < 4s^2*t, so -tu > -4s^2*t.
From v < 4t^2 (and s > 0): sv < 4st^2, so -sv > -4st^2.
Therefore: 12st(s+t) - tu - sv > 12st(s+t) - 4s^2*t - 4st^2 = 8st(s+t) > 0.
**PASS.** The bound direction is correct; the proof only needs the *upper* bounds on u, v.

### Check 3: Cauchy-Schwarz application
With a_1 = sqrt(sA), b_1 = u/sqrt(sA), a_2 = sqrt(tB), b_2 = v/sqrt(tB):
- a_1*b_1 + a_2*b_2 = u + v
- a_1^2 + a_2^2 = sA + tB
- b_1^2 + b_2^2 = u^2/(sA) + v^2/(tB)

CS gives: (u+v)^2 <= (sA+tB)*(u^2/(sA) + v^2/(tB)). The substitution is valid since sA > 0 and tB > 0 on the domain. **PASS.**

### Check 4: Logical chain
- Step 1: sA+tB <= (s+t)C, all positive. So 1/((s+t)C) <= 1/(sA+tB). Multiplying by (u+v)^2 >= 0: (u+v)^2/((s+t)C) <= (u+v)^2/(sA+tB).
- Step 2: (u+v)^2/(sA+tB) <= u^2/(sA) + v^2/(tB).
- Combined: (u+v)^2/((s+t)C) <= u^2/(sA) + v^2/(tB), which is (*).
- (*) is equivalent to R_4(s+t,0,u+v) >= R_4(s,0,u) + R_4(t,0,v) because R_4 = -(1/24)*K4^2/(K2*(4K2^2-K4)) and the -(1/24) factor flips the inequality.

**PASS.** Direction is correct throughout.

### Check 5: Domain constraints
The domain for R_4(K2, 0, K4) requires:
- D1: 4K2^2 - K4 > 0  (i.e., K4 < 4K2^2)
- D2: 4K2^3 + K2*K4 > 0, i.e., K4 > -4K2^2 (since K2 > 0)

Combined: K4 in (-4K2^2, 4K2^2).

For the sum (s+t, u+v):
- D1: 4(s+t)^2 - (u+v) > 4(s+t)^2 - 4s^2 - 4t^2 = 8st > 0. **AUTOMATIC.**
- D2: 4(s+t)^2 + (u+v) > 4(s+t)^2 - 4s^2 - 4t^2 = 8st > 0. **AUTOMATIC.**

The sum domain constraints are redundant, automatically implied by the individual constraints. **PASS.**

### Check 6: Edge/boundary cases
Tested numerically:
- u -> 4s^2 (A -> 0): gap -> +infinity. OK.
- u -> -4s^2 (D2 -> 0): gap remains positive. OK.
- u near +boundary, v near -boundary: gap large and positive. OK.
- Very asymmetric s >> t: gap positive. OK.
- Both u, v near lower boundary: gap positive (converges to finite limit). OK.

**PASS.** No boundary issues.

### Check 7-8: Numerical verification
- 0 violations in 999,324 valid trials with aggressive sampling (including boundary-focused, asymmetric, and tiny-value strategies).
- Minimum gap: ~6.7e-11 (near machine epsilon, at u, v near zero).

**PASS.**

### Check 9-10: Bound analysis
The lower bound 8st(s+t) is tight at the upper boundary (u=4s^2, v=4t^2) and loose at the lower boundary (margin 16st(s+t)). The proof does not need the lower bounds u > -4s^2, v > -4t^2 for Steps 1-2, only for the domain of R_4 itself. **PASS.**

### Check 11: Sign handling
The equivalence chain:
```
R_4(s+t,0,u+v) >= R_4(s,0,u) + R_4(t,0,v)
iff -(1/24)[(u+v)^2/((s+t)C)] >= -(1/24)[u^2/(sA)] + -(1/24)[v^2/(tB)]
iff -(u+v)^2/((s+t)C) >= -u^2/(sA) - v^2/(tB)
iff (u+v)^2/((s+t)C) <= u^2/(sA) + v^2/(tB)
```
**PASS.** The negative sign in R_4 correctly flips the inequality.

### Check 12-14: Independent verifications
- Step-by-step numerical verification: 0 failures in 492,686 valid trials. **PASS.**
- Trivial case u=v=0: gap=0. **PASS.**
- Symmetric case s=t, u=v: gap = u^2*(12s^2-u)/(24s*(8s^2-u)*(4s^2-u)) >= 0 on domain. **PASS.**

---

## 3. Issues Found

**None.** The proof is complete and correct.

### Minor notes (not affecting validity):
- The theorem statement includes the hypothesis u+v in (-4(s+t)^2, 4(s+t)^2), but this is automatically implied by the individual domain constraints. Not an error, just a redundancy.
- The proof is actually slightly stronger than stated: Step 1 shows (s+t)C - sA - tB >= 8st(s+t) > 0, not just > 0.

---

## 4. Verification Script

See `examples7/verifier11_check.py` for the full verification code (15 independent checks, ~1.5M total numerical trials across all tests).

---

## 5. Verdict

**VALID.** The k3=0 case of R_4 superadditivity is correctly proved by PROVER-11's two-step argument (Key Lemma + Cauchy-Schwarz). The algebra is correct, the inequality directions are correct, the domain constraints are complete, and no counterexamples exist.
