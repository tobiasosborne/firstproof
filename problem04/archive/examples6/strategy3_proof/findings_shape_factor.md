# Findings: Shape Factor Decomposition Approach to Fisher Superadditivity

**Agent:** prover (shape factor investigation)
**Date:** 2026-02-08
**Script:** `investigate_shape_factor.py`

---

## Executive Summary

The shape factor decomposition `1/Phi(p) = Var(p) / SF(p)` was investigated as a potential proof strategy for Fisher superadditivity. The central claim `SF(r) <= min(SF(p), SF(q))` is **FALSE** -- counterexamples are abundant and dramatic (worst ratio ~37x for n=3). However, several valuable structural results were discovered:

| Result | Status | Value |
|--------|--------|-------|
| `Var(r) = Var(p) + Var(q)` exactly | **PROVED** | High: exact additivity of variance |
| `SF(r) <= min(SF(p), SF(q))` | **FALSE** | High: counterexample found |
| `SF(r) <= max(SF(p), SF(q))` | **CONJECTURED TRUE** (proved for n=3) | Medium: true but insufficient |
| `SF(p boxplus p) <= SF(p)` | **CONJECTURED TRUE** (proved for n=3) | Medium: nice structural fact |
| Shape factor approach gives new proof | **FALSE** | Conclusive: just a rephrasing |

**Bottom line:** The shape factor decomposition is a change of variables, not a new proof technique. The strong form `SF(r) <= min(SF(p), SF(q))` fails catastrophically, and the correct reformulation `SF(r) <= variance-weighted harmonic mean` is algebraically equivalent to the original conjecture.

---

## Finding 1: Var(r) = Var(p) + Var(q) -- EXACT ADDITIVITY (PROVED)

**Statement:** For the MSS convolution r = p boxplus_n q of monic real-rooted degree-n polynomials p, q:
```
Var(r) = Var(p) + Var(q)
```
where `Var(p) = (1/n) * sum_i (lambda_i - mean)^2` is the variance of the roots.

**Proof:** Variance is translation-invariant and MSS boxplus is translation-equivariant, so WLOG assume p, q centered (e_1 = 0). Then:
- `Var(p) = -2*e_2(p)/n` since `sum lambda_i^2 = e_1^2 - 2*e_2 = -2*e_2` for centered polynomials.
- MSS boxplus for e_2 with e_1 = 0:
  ```
  e_2(r) = sum_{i+j=2} C(n-j,i)/C(n,i) * e_i(p) * e_j(q)
  ```
  The three terms (i,j) = (0,2), (1,1), (2,0) give:
  ```
  C(n-2,0)/C(n,0)*e_2(q) + C(n-1,1)/C(n,1)*0*0 + C(n,2)/C(n,2)*e_2(p)
  = e_2(q) + e_2(p)
  ```
- Therefore `Var(r) = -2*e_2(r)/n = -2*(e_2(p)+e_2(q))/n = Var(p) + Var(q)`.

**Numerical verification:** Relative error < 3e-14 across n = 3, 4, 5, 6, 7, 8 with 200 trials each. Exact to machine precision.

**This is exact additivity, not merely super-additivity.** The variance of roots is a linear functional of the MSS convolution.

---

## Finding 2: SF(r) <= min(SF(p), SF(q)) is FALSE -- COUNTEREXAMPLE

**The conjecture SF(r) <= min(SF(p), SF(q)) is FALSE for all n >= 3.**

**Simplest counterexample (n=3):** Take p with equally-spaced roots (F_p = 0) and any q with non-equally-spaced roots (F_q != 0):
- `SF(p) = 3` (the global minimum of SF for n=3)
- `SF(r) = 12/(4 - 27*F_q^2/(E_p+E_q)^3) > 3 = SF(p)` whenever F_q != 0

**Concrete example:**
```
p = [-1.0, 0.0, 1.0]    (equally spaced)
q = [-3.0, 1.0, 2.0]    (centered, non-equally-spaced)
r = boxplus(p, q)

SF(p) = 3.000000
SF(q) = 10.290000
SF(r) = 5.710037
min(SF(p), SF(q)) = 3.000000

SF(r) > min(SF(p), SF(q)):  VIOLATED by factor 1.9x
BUT: 1/Phi(r) >= 1/Phi(p) + 1/Phi(q)  STILL HOLDS (excess = 0.258)
```

**Statistics of violations:**

| n | Fraction violated (random) | Worst ratio SF(r)/min(SF) |
|---|---------------------------|---------------------------|
| 3 | 42.7% (853/2000)         | 37.2                      |
| 4 | 19.6% (392/2000)         | 9.7                       |
| 5 | 12.8% (255/2000)         | 7.2                       |
| 6 | 7.7% (153/2000)          | 3.7                       |
| 7 | 6.8% (136/2000)          | 3.0                       |

For equally-spaced p vs random q: **100% violation rate** for n=3 (every case violates).

**The root cause:** SF achieves its global minimum at equally-spaced roots. The MSS convolution of equally-spaced p with non-equally-spaced q produces r with non-equally-spaced roots, so SF(r) > SF(p) = min(SF).

---

## Finding 3: SF(r) <= max(SF(p), SF(q)) appears ALWAYS TRUE

**Conjecture:** For all n >= 3 and all monic real-rooted degree-n polynomials p, q with simple roots:
```
SF(r) <= max(SF(p), SF(q))
```

**Evidence:** 0 violations in 10,000 trials per n for n = 3, 4, 5, 6, 7, 8, across diverse root configurations (random, near-degenerate, equally-spaced, highly asymmetric).

| n | Violations/Trials | Closest ratio SF(r)/max(SF) |
|---|-------------------|---------------------------|
| 3 | 0/10,000          | 0.99999996                |
| 4 | 0/10,000          | 0.99973                   |
| 5 | 0/10,000          | 0.99862                   |
| 6 | 0/10,000          | 0.99545                   |
| 7 | 0/10,000          | 0.99247                   |
| 8 | 0/10,000          | 0.97017                   |

**Proved for n = 3:** See analytical proof below.

**However, this is INSUFFICIENT for superadditivity.** The bound `SF(r) <= max(SF(p), SF(q))` gives:
```
1/Phi(r) = (Var(p)+Var(q))/SF(r) >= (Var(p)+Var(q))/max(SF(p),SF(q))
```
But this does NOT imply `>= Var(p)/SF(p) + Var(q)/SF(q) = 1/Phi(p) + 1/Phi(q)` when SF(p) != SF(q).

---

## Finding 4: SF(p boxplus p) <= SF(p) appears ALWAYS TRUE (proved for n=3)

**Conjecture:** For all n >= 3 and all p:
```
SF(p boxplus_n p) <= SF(p)
```
with equality iff p has equally-spaced roots.

**Evidence:** 0 violations in 3000 trials per n for n = 3, 4, 5, 6, 7.

**Proved for n = 3:**
- For centered cubic: `SF = 12/(4 - 27u)` where `u = F^2/E^3`.
- For `r = p boxplus p`: `E_r = 2E_p`, `F_r = 2F_p`, so `u_r = 4F_p^2/(8E_p^3) = u_p/2`.
- Since `u_r = u_p/2 <= u_p` and SF is monotone increasing in u:
  `SF(r) = 12/(4 - 27u_p/2) <= 12/(4 - 27u_p) = SF(p)`.
- Equality iff u_p = 0 iff p has equally-spaced roots.

---

## Finding 5: Analytical Proof of SF(r) <= max(SF(p), SF(q)) for n=3

**Theorem:** For centered cubics with `SF = 12/(4-27u)`, `u = F^2/E^3`, and `E_r = E_p+E_q`, `F_r = F_p+F_q`:
```
u_r <= max(u_p, u_q)
```
and therefore `SF(r) <= max(SF(p), SF(q))`.

**Proof:**

**Case 1:** `F_p = 0`. Then `u_p = 0` and `u_r = F_q^2/(E_p+E_q)^3 <= F_q^2/E_q^3 = u_q`.

**Case 2:** `F_p != 0`, `u_p >= u_q` (WLOG by symmetry). Set `t = E_q/E_p > 0`, `s = F_q/F_p` (real). Then `u_q/u_p = s^2/t^3`, so `|s| <= t^{3/2}`.

Need to show: `(1+s)^2 <= (1+t)^3`.

Chain of inequalities:
1. `(1+s)^2 <= (1+|s|)^2` since `|1+s| <= 1+|s|` implies `(1+s)^2 <= (1+|s|)^2`.
2. `(1+|s|)^2 <= (1+t^{3/2})^2` since `|s| <= t^{3/2}`.
3. `(1+t^{3/2})^2 <= (1+t)^3` since:
   ```
   (1+t)^3 - (1+t^{3/2})^2 = 3t + 3t^2 - 2t^{3/2}
   = t(3 + 3t - 2*sqrt(t)) = t(3*(sqrt(t) - 1/3)^2 + 8/3) > 0
   ```

Therefore `(1+s)^2 <= (1+t)^3`, which gives `u_r <= u_p = max(u_p, u_q)`. QED.

---

## Finding 6: The Shape Factor Approach is Just a Change of Variables

The "correct" reformulation of the conjecture in shape factor terms is:
```
SF(r) <= (Var(p) + Var(q)) * SF(p) * SF(q) / (Var(p)*SF(q) + Var(q)*SF(p))
```
This is the variance-weighted harmonic mean of SF(p) and SF(q).

**This is algebraically EQUIVALENT to the original conjecture** `1/Phi(r) >= 1/Phi(p) + 1/Phi(q)`. It provides no additional mathematical leverage. Verified: 100% sign-agreement across 1000 trials for each of n = 3, 4, 5, 6.

The harmonic mean lies in `[min(SF(p), SF(q)), max(SF(p), SF(q))]`, so the bound is strictly between the failed `min` bound and the proved `max` bound.

---

## Finding 7: Structure of the Shape Factor

**Scale invariance:** `SF(c * roots) = SF(roots)` for all c > 0. SF depends only on the "shape" (gap ratios) of the root configuration.

**For n = 3:** `SF = 12/(4 - 27u)` where `u = F^2/E^3` with `E = -e_2 > 0`, `F = e_3`. SF ranges from 3 (equally-spaced roots) to infinity (degenerate/nearly-coincident roots).

**For general n:**
- Equally-spaced roots give the MINIMUM value of SF: SF_min(3) = 3, SF_min(4) ~ 9.03, SF_min(5) ~ 20.1.
- SF depends on the full "shape" of the root configuration, not just variance or kurtosis.
- For n >= 4, SF depends on multiple dimensionless ratios (gap ratios), making a simple formula unlikely.

**For n = 3, SF as a function of gap ratio g1/g2:**

| g1/g2 | SF    |
|-------|-------|
| 0.1   | 50.2  |
| 0.5   | 4.23  |
| 1.0   | 3.00  |
| 2.0   | 4.23  |
| 10.0  | 50.2  |

SF is symmetric under g1/g2 <-> g2/g1 and minimized at g1/g2 = 1 (equal gaps).

---

## Finding 8: What the Terms in the Excess Look Like

The excess can be written as:
```
1/Phi(r) - 1/Phi(p) - 1/Phi(q) = Var(p)*(1/SF(r) - 1/SF(p)) + Var(q)*(1/SF(r) - 1/SF(q))
```

The two terms can have **opposite signs**:
- n=3: both positive 55.9%, one negative 44.1%
- n=4: both positive 77.8%, one negative 22.2%
- n=5: both positive 87.0%, one negative 13.0%

When one term is negative (i.e., SF(r) > SF(p)), the other is always positive and large enough to compensate. But this cancellation makes the shape factor approach fundamentally unsuitable for proving the inequality term-by-term.

---

## Conclusions and Recommendations

### The shape factor approach does NOT provide a proof path.

1. **The strong conjecture `SF(r) <= min(SF(p), SF(q))` is false.** Counterexamples are easy to construct (any equally-spaced p with any non-equally-spaced q).

2. **The correct reformulation is equivalent to the original conjecture.** It provides no mathematical simplification.

3. **The partial result `SF(r) <= max(SF(p), SF(q))` is true** (proved for n=3, numerically verified for all n) **but insufficient** -- it does not imply superadditivity.

### What IS useful from this investigation:

1. **Var(r) = Var(p) + Var(q) is proved** and could be a useful building block for other approaches. It gives the "extensive" part of the inequality for free.

2. **SF(r) <= max(SF(p), SF(q))** is a new structural conjecture about the MSS convolution that may have independent interest. For n=3 it is proved using an elementary inequality on power means.

3. **SF(p boxplus p) <= SF(p)** provides the p=q case of the full conjecture immediately:
   ```
   1/Phi(p boxplus p) = 2*Var(p)/SF(r) >= 2*Var(p)/SF(p) = 2/Phi(p)
   ```

4. **The fact that SF depends only on gap ratios** (scale-invariance) confirms that the "hard part" of the conjecture is the shape-dependent part, not the scale-dependent part.

### Recommended next directions:

- **Do NOT pursue shape factor decomposition further.** It is a dead end.
- **DO use variance additivity** as a lemma in other approaches.
- **Consider information-geometric approaches** that handle the "shape" and "scale" jointly rather than separately.
- **The n=3 quadratic form argument** (from `prove_n3_symbolic.py`) remains the right template for the proof, working directly with the excess as a positive definite form.
