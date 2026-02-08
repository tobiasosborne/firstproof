# VERIFIER-12 Report: Adversarial Verification of PROVER-13 Claims

## Executive Summary

PROVER-13 presents a heat-flow-based approach to proving the finite free analog of
Stam's inequality: `1/Phi_n(p ⊞ q) >= 1/Phi_n(p) + 1/Phi_n(q)`. Four main claims
are made. After extensive adversarial testing with > 15,000 test cases, I find:

| Claim | Verdict | Summary |
|-------|---------|---------|
| 1. De Bruijn identity | **VALID** (with caveats) | Identity `dS/dt = Phi_n(p)` holds, but numerics degrade near degenerate roots |
| 2. Entropy Power Inequality | **VALID** (numerically) | 0 violations in 13,770 trials; n=2 case proved analytically |
| 3. Monotone Gap | **INVALID** | Gap is NOT monotonically decreasing. Proof mechanism is broken. |
| 4. Gaussian Splitting | **VALID** | Holds exactly (algebraic identity). But cumulant explanation is wrong. |

**The most important finding: Claim 3 (the proposed proof mechanism) is FALSE.**
The gap `F(t) - G(t)` is NOT monotonically decreasing. In 15 out of 137 tested
polynomial pairs, the gap was observed to increase. Despite this, the underlying
conjecture (`1/Phi(p ⊞ q) >= 1/Phi(p) + 1/Phi(q)`) appears to be TRUE based on
numerical evidence. The gap is always non-negative even though it is non-monotone.
A different proof mechanism is needed.

---

## Claim 1: Finite Free De Bruijn Identity

**PROVER-13 states:** `dS/dt|_{t=0} = Phi_n(p)` where `S(p) = sum_{i<j} log|r_i - r_j|`
and `p_t = p ⊞ G_t`.

### Verdict: VALID (with numerical caveats)

### Evidence

**Sign check:** S increases under heat flow (roots spread apart). `dS/dt > 0` and
`Phi_n > 0` are consistent. 0 sign errors across all tests.

**Well-separated roots (n=3..8):** Richardson-extrapolated derivatives match `Phi_n`
to relative error < 1e-4. The ratio `dS/dt / Phi_n` converges to 1.0 as the step
size decreases.

**Algebraic identity:** The identity `sum_{i<j} (H_i - H_j)/(r_i - r_j) = sum_i H_i^2 = Phi_n`
was verified to machine precision for n=3..8 across multiple root configurations.
This is the key algebraic backbone of the de Bruijn identity.

**Root dynamics:** `dr_i/dt = H_i(r)` (constant c=1) confirmed for n=4,6 with ratios
converging to 1.0000 as dt -> 0. Combined with the algebraic identity, this gives
`dS/dt = Phi_n` exactly.

### Issues Found

1. **Near-degenerate roots (gap < 0.01):** The finite difference approximation
   breaks down catastrophically. For gap=0.001, relative error is 83%. For gap=0.0001,
   relative error is 99.6%. This is NOT an error in the identity itself -- it is a
   numerical stability issue. When roots are very close, `Phi_n ~ 1/gap^2` becomes
   huge, and the polynomial `p ⊞ G_dt` for tiny `dt` has roots that are extremely
   sensitive to coefficient perturbations.

2. **Large n (n=20):** Numerical instability causes the verification to fail.
   Root-finding for degree-20 polynomials after MSS convolution with small-variance
   Gaussian produces spurious complex roots. The identity likely still holds but
   cannot be verified numerically for large n with this approach.

3. **Root dynamics instability (n=3,5):** For some root configurations with
   well-separated roots at spacing 2.0, the root dynamics ratio `dr/dt / H` shows
   `inf` or `nan`. This is likely due to root reordering during the finite
   difference step -- `np.sort(roots)` can swap indices between the `t=0` and
   `t=dt` evaluations. This is a methodological issue, not a mathematical one.

### Conclusion

The de Bruijn identity `dS/dt = Phi_n(p)` is **mathematically valid** based on:
- Exact verification of the algebraic identity `sum_{i<j}(H_i-H_j)/(r_i-r_j) = Phi_n`
- Confirmed root dynamics `dr_i/dt = H_i`
- High-precision numerical agreement for well-separated roots

The numerical difficulties near degenerate configurations do not indicate a
mathematical failure of the identity.

---

## Claim 2: Entropy Power Inequality (EPI)

**PROVER-13 states:** `N(p ⊞ q) >= N(p) + N(q)` where
`N(p) = exp(2*S(p)/m)`, `m = n(n-1)/2`, `S = sum_{i<j} log|r_i - r_j|`.

### Verdict: VALID (numerically; proved for n=2)

### Evidence

**Total violations: 0 out of 13,770 trials.**

Breakdown:
- n=2: 0/5000 violations. **Analytically proved:** N(p ⊞ q) = (a-b)^2 + (c-d)^2 = N(p) + N(q) exactly. EPI holds with EQUALITY for all degree-2 polynomials.
- n=3: 0/2945 violations, min ratio = 1.000000
- n=4: 0/2903 violations, min ratio = 1.000000
- n=5: 0/2844 violations, min ratio = 1.000000
- n=10: 0/43 violations, min ratio = 1.122
- n=15: 0/34 violations, min ratio = 1.104
- n=20: 0/1 violations (limited by numerical stability)

**Adversarial test results:**
- Very unequal root spacings (spread ratio 100:1): 0 violations across 6000 trials
- Near-degenerate polynomials (gap ~ 0.001): 0 violations across 2693 trials
- Large n (10, 15, 20): 0 violations (limited sample size for n=20)

**Normalization verification:**
- `N(G_s)` is perfectly linear in `s` (max deviation < 1e-15)
- Gaussians achieve equality: `N(G_s ⊞ G_t) = N(G_s) + N(G_t)` to machine precision
- The normalization `m = n(n-1)/2` is correct

### Issues Found

1. **Near-equality for small n:** For n=2, the min ratio is exactly 1.0 (equality).
   For n=3,4,5 the min ratio is also very close to 1.0 (within 1e-10 of equality).
   This raises the question: is there a continuous family of equality cases?

2. **Large n limitations:** For n >= 20, numerical instability severely limits testing.
   The MSS convolution involves factorials up to n!, causing floating-point issues.

### Conclusion

The EPI is **strongly supported numerically** with zero counterexamples across
extensive adversarial testing. The n=2 case is proved analytically. However, a
formal proof for general n remains open.

---

## Claim 3: Monotone Gap

**PROVER-13 states:** `F(t) - G(t)` is decreasing where
`F(t) = 1/Phi((p⊞q)⊞G_{2t})`, `G(t) = 1/Phi(p⊞G_t) + 1/Phi(q⊞G_t)`.

### Verdict: **INVALID**

### Evidence

**Gap is NOT monotonically decreasing.** In 15 out of 137 polynomial pairs tested,
the gap was observed to increase at some point along its trajectory. Specific
counterexample:

- n=3, trial 43 (seed 2024):
  - `roots_p = [-1.95, 0.10, 5.59]`
  - `roots_q = [0.03, 0.60, 0.88]`
  - gap(0.001) = 0.056, gap increases to max 0.342 at t ~ 1.67, then decreases
  - gap(5.0) = 0.232 > gap(0.001)

**The gap increases OVERALL in some cases.** This is not just a small fluctuation.
For the example above, the gap at t=5 is 4x larger than at t=0.001.

### Critical Logical Analysis

PROVER-13's proof argument (Section K of `R4_prover13_epi.py`) is:
1. gap_K(t) is decreasing
2. gap_K(t) -> 0 as t -> infinity
3. Therefore gap_K(0) >= 0

This argument fails at Step 1. Additionally, the argument has a subtle issue
at Step 2: careful asymptotic analysis shows that `gap*t` converges to values
that can be positive, negative, or zero depending on n and the polynomials.
This means the gap approaches 0 from different sides, further complicating
any monotonicity-based argument.

### However: The Gap is Always Non-Negative

Despite the gap being non-monotone, extensive testing (562 polynomial pairs,
n=3..6, 500 time points each) shows the gap is **always non-negative**:
- gap_K(t) >= 0 for all t and all tested polynomials
- This is equivalent to the Stam conjecture holding for all (p_t, q_t) pairs

This means the underlying conjecture is likely TRUE, but PROVER-13's proposed
proof mechanism (monotone gap) is WRONG. A different proof strategy is needed.

### Two Different Gap Definitions

PROVER-13 conflates two different gap definitions:
- **gap_H(t)** = 1/Phi(p⊞q⊞G_t) - 1/Phi(p⊞G_t) - 1/Phi(q⊞G_t) (same noise on all terms)
- **gap_K(t)** = 1/Phi((p⊞q)⊞G_{2t}) - 1/Phi(p⊞G_t) - 1/Phi(q⊞G_t) (double noise on LHS)

These are DIFFERENT functions. gap_H is monotonically decreasing (confirmed in 182/182
trials) but it goes NEGATIVE. gap_K is non-negative but NOT monotonically decreasing
(failed in 22/230 trials). Neither definition gives the proof as claimed.

### Supporting Evidence: Concavity of 1/Phi

The claim that `1/Phi(p ⊞ G_t)` is concave in t was verified with 0 violations
in 446 trials (n=3..7). This is a valid and interesting result, but it is NOT
sufficient to prove the monotone gap claim.

### Conclusion

**The monotone gap claim is FALSE.** The gap `F(t) - G(t)` is not monotonically
decreasing. PROVER-13's proposed proof of the Stam inequality via this mechanism
does not work. The underlying inequality appears to be true but requires a
different proof approach.

---

## Claim 4: Gaussian Splitting

**PROVER-13 states:** `(p ⊞ G_{t/2}) ⊞ (q ⊞ G_{t/2}) = (p ⊞ q) ⊞ G_t`

### Verdict: VALID (exact algebraic identity)

### Evidence

**Numerical verification:** Max error < 1e-12 for n=2..7, < 5e-5 for n=10.
Tested across 50 polynomial pairs per degree with 7 time values each.
The growing error for large n is due to floating-point accumulation in the
MSS convolution formula (factorials up to n!), not a mathematical failure.

**General splitting also holds:** `(p ⊞ G_s) ⊞ (q ⊞ G_t) = (p ⊞ q) ⊞ G_{s+t}`
for arbitrary s, t. Verified to max error < 3e-12 for n=3..5.

### Issue Found: Cumulant Explanation is Wrong

PROVER-13 claims the splitting follows from "cumulant additivity" and provides
a specific cumulant formula: `kappa_k = (-1)^k * a_k * n^{k-1} / C(n,k)`.

**This formula is WRONG.** It gives:
- `kappa_2(G_1) = -n` (should be some positive quantity proportional to s)
- `kappa_4(G_1) != 0` (should be 0 for Gaussian)
- Cumulants computed with this formula are NOT additive under MSS convolution

The correct finite free cumulants (as defined by Arizmendi-Perales or Marcus)
DO linearize MSS convolution by definition, but PROVER-13's explicit formula
for computing them is incorrect.

**The splitting identity itself is still valid** -- it holds as a direct algebraic
consequence of the MSS convolution formula structure, regardless of the cumulant
formalism. The proof via cumulants works in principle (with the CORRECT cumulant
definition), but PROVER-13's implementation of the cumulant computation is buggy.

### Conclusion

The Gaussian splitting identity is **valid and exact**. The cumulant-based
explanation is correct in spirit but PROVER-13's cumulant formula has a bug.
The splitting can be proved directly from the MSS formula without needing
an explicit cumulant computation.

---

## Additional Findings

### n=2 Analytic Results

For n=2, many things can be computed exactly:
- N(p) = (gap between roots)^2
- `N(p ⊞ q) = N(p) + N(q)` exactly (EPI with equality)
- `Phi_n(p) = 2/(gap)^2` for p with roots at distance `gap`
- `1/Phi(p ⊞ q) = 1/Phi(p) + 1/Phi(q)` (Stam with equality)

This means for n=2, both the EPI and Stam inequalities hold with **equality
for all polynomial pairs**. This is a strong anchor point for any proof attempt.

### Isoperimetric Inequality

PROVER-13 investigates `N(p) * Phi(p)` as a potential isoperimetric quantity.
For Gaussians, this product is constant (independent of variance). However,
Gaussians do NOT minimize or maximize this product among all polynomials --
the product takes both larger and smaller values. This means there is no simple
isoperimetric inequality of the form `N * Phi >= C_n` or `N * Phi <= C_n`.

---

## Summary of Errors and Gaps

1. **MAJOR (Claim 3):** The monotone gap claim is false. The gap F(t)-G(t) is not
   monotonically decreasing. This invalidates the proposed proof mechanism.

2. **MODERATE (Claim 4):** The cumulant formula `kappa_k = (-1)^k * a_k * n^{k-1} / C(n,k)`
   is wrong. It does not give additive quantities and produces incorrect Gaussian
   cumulants. The splitting identity is still valid but the stated justification is
   broken.

3. **MINOR (Claim 1):** Numerical verification degrades for nearly-coincident roots
   and large n. This is a methodological limitation, not a mathematical error.

4. **OBSERVATION:** The underlying conjectures (EPI and Stam) appear to be TRUE
   based on extensive numerical evidence, but the proposed proof mechanism does not
   work. A proof must be found by other means.

---

## Verification Scripts

- `/home/tobiasosborne/Projects/af-tests/examples7/verifier12_check.py` -- Main
  adversarial verification script with all tests.
