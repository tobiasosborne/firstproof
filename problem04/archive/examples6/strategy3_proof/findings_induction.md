# Findings: Inductive Approach to Fisher Superadditivity

## Investigation: Node 1.10.1

**Investigator:** induction-prover
**Date:** 2026-02-08

---

## 1. Critical Bug Fix: The Correct Boxplus Formula

The MSS finite free additive convolution formula used in earlier scripts was
**WRONG for non-centered polynomials**. The correct formula is:

```
g_k = sum_{i+j=k} C(n-j, i) / C(n, i) * e_i(p) * e_j(q)
```

where `e_k` are the elementary symmetric polynomials of the roots. This was
verified against Monte Carlo simulation of `E[det(xI - A - UBU*)]` over
Haar-random unitaries.

The previously used formula with `w = (n-i)!(n-j)!/(n!(n-k)!)` applied to
*normalized* coefficients `a_k = e_k/C(n,k)` gives the SAME result for
centered polynomials (mean of roots = 0), but DIVERGES for non-centered ones.

**Impact:** All previous scripts in this directory that used `boxplus_mss()`
with the old formula gave incorrect results for non-centered polynomials.
The conjecture appeared to be violated when it was actually valid.

## 2. Conjecture Status: NUMERICALLY CONFIRMED

With the corrected boxplus formula:

```
1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q)
```

where `r = p boxplus_n q`, is **confirmed for n = 2, 3, 4, 5, 6** over
500+ random trials each. Zero violations found.

| n | Trials | Violations | Min excess |
|---|--------|-----------|------------|
| 2 | 500    | 0         | ~0 (equality) |
| 3 | 500    | 0         | ~0 (near equality for specific cases) |
| 4 | 500    | 0         | 0.0004     |
| 5 | 500    | 0         | 0.0007     |
| 6 | 500    | 0         | 0.0005     |

## 3. Base Case n=2: Exact Equality

For n=2 with roots a < b and c < d:

- `Phi_2(p) = 2/(b-a)^2`
- `1/Phi_2(p) = (b-a)^2/2`
- The correct boxplus gives `gap_r^2 = gap_p^2 + gap_q^2` **exactly** (Pythagorean).
- So `1/Phi_2(r) = 1/Phi_2(p) + 1/Phi_2(q)` with **equality**.

Proof sketch:
- `g_2 = e_2(p) + e_1(p)*e_1(q)/2 + e_2(q)` = ab + (a+b)(c+d)/2 + cd
- disc = `(a+b+c+d)^2 - 4*g_2 = (a-b)^2 + (c-d)^2`

## 4. MSS Derivative Identity: Verified

The identity `r'(x) = n * (p^{(1)} boxplus_{n-1} q^{(1)})(x)` where
`p^{(1)} = (1/n)*p'` is verified to **machine precision** (~1e-15) for
n = 3, 4, 5 with the correct boxplus formula.

(The old formula gave errors of order 1e-2 for this identity, which was
another sign that it was wrong.)

## 5. Induction Decomposition

The exact identity:
```
F_n(p,q) = F_{n-1}(p^{(1)}, q^{(1)}) + [delta_r - delta_p - delta_q]
```
where:
- `F_n(p,q) = 1/Phi_n(r) - 1/Phi_n(p) - 1/Phi_n(q)` (excess)
- `delta_p = 1/Phi_n(p) - 1/Phi_{n-1}(p^{(1)})` (level drop)
- `r^{(1)}` roots = `p^{(1)} boxplus_{n-1} q^{(1)}` roots (by MSS)

### Key Findings:

**(a) delta_p is ALWAYS NEGATIVE:**
- For all n tested (3-6), `1/Phi_n(p) < 1/Phi_{n-1}(p^{(1)})` always.
- Intuition: higher degree means more H_p terms, so Phi_n > Phi_{n-1}.

**(b) The correction `delta_r - delta_p - delta_q` is NOT always >= 0:**
- For n=4: 193/500 samples have negative correction.
- For n=5: 232/500 samples have negative correction.
- So direct induction via "delta is superadditive" FAILS.

**(c) BUT: F_n/F_{n-1} is ALWAYS POSITIVE** (the key finding!):
- For n=4: `F_4/F_3` ranges from 0.25 to 218,000 -- always positive.
- For n=5: `F_5/F_4` ranges from 0.38 to 15.4 -- always positive.
- This means: when the correction is negative, `F_{n-1}` more than compensates.

**(d) The conjecture holds at EVERY derivative level:**
- For p, q of degree 5, the excess `F_k` is >= 0 for k = 5, 4, 3, 2.
- This is consistent with the inequality being TRUE for all degrees.

## 6. Phi Ratio for Uniform Spacing

For uniformly spaced roots with gap d:

| n | Phi_n/Phi_{n-1} |
|---|-----------------|
| 3 | 3.0000          |
| 4 | 2.0062          |
| 5 | 1.6746          |
| 6 | 1.5086          |
| 7 | 1.4089          |
| 8 | 1.3422          |

The ratio is NOT constant for general root configurations, but the
uniform-spacing ratio converges to approximately 1 as n grows.

## 7. Assessment of the Inductive Approach

### What works:
- Base case n=2: proved (equality).
- MSS derivative identity: exact, connects degree n to degree n-1.
- The decomposition `F_n = F_{n-1} + correction` is exact.
- Numerically, `F_n/F_{n-1} > 0` always (strong evidence).

### What does NOT work:
- The correction term `delta_r - delta_p - delta_q` is not always >= 0.
- There is no simple formula relating `Phi_n(p)` to `Phi_{n-1}(p^{(1)})`.
- The Phi ratio depends on the full root configuration, not just gaps.

### The gap:
The inductive approach requires proving ONE of:
1. `F_n/F_{n-1} > 0` directly (i.e., that F_n and F_{n-1} have the same sign).
2. A lower bound on `F_{n-1}` that dominates the negative correction.
3. A different decomposition where the correction IS non-negative.

None of these are achieved by the current analysis. The ratio `F_n/F_{n-1}`
being positive is numerically robust but conceptually unclear -- it says
"the excess at degree n is positive whenever the excess at degree n-1 is
positive," which is essentially restating the conjecture at two levels.

## 8. Conclusion

**The inductive approach via MSS derivatives does NOT close.** The key
structural fact (derivative convolution = convolution of derivatives) is
verified and exact, but the relationship between `Phi_n` and `Phi_{n-1}`
is too complex for simple bounding.

The most promising lead is the ratio `F_n/F_{n-1} > 0` (always positive
numerically). If one could show this ratio is bounded below by some positive
function of the root data, the induction would close. But this seems as
hard as the original conjecture.

**Alternative approaches suggested by the data:**
1. The n=3 case has `F_3 >= 0` with equality exactly when uniform spacing
   is preserved. A direct proof of n=3 (without induction) might be tractable.
2. The "convolution telescoping" (Section 8 of the script) shows that `F_k >= 0`
   at ALL derivative levels simultaneously. This suggests a proof that works
   at all levels at once, possibly via a monotonicity argument on the derivative
   chain, rather than level-by-level induction.
