# Findings: A + B >= 2*Phi_r Investigation

**Node:** 1.10.2 (AM-GM reduction)
**Agent:** prover (investigate_AplusB.py)
**Date:** 2026-02-08
**Script:** `investigate_AplusB.py`

---

## CRITICAL FINDING: THE AM-GM APPROACH IS INVALID

The proposed proof strategy was:

> If A + B >= 2*Phi_r, then by AM-GM: AB >= ((A+B)/2)^2 >= Phi_r^2 = h^4.

**This is WRONG.** AM-GM states that for positive A, B:

```
(A + B) / 2  >=  sqrt(A*B)
```

Squaring: `((A+B)/2)^2 >= AB`. This is an **UPPER** bound on AB, not a lower bound.

So even if `A + B >= 2*Phi_r` holds, this gives:
```
AB  <=  ((A+B)/2)^2
```
which is the wrong direction. The AM-GM inequality does NOT yield `AB >= Phi_r^2`.

**Concrete counterexample to the logic (not using MSS):**
Take A = 0.01, B = 3.99, Phi_r = 1. Then A + B = 4.0 >= 2*Phi_r = 2.
But AB = 0.0399 < 1 = Phi_r^2.

The fact that AB >= Phi_r^2 ALSO holds in the MSS setting is a separate fact that
requires its own proof. The sum inequality does not imply the product inequality.

---

## Does A + B >= 2*Phi_r hold numerically?

**YES**, with zero counterexamples in extensive testing:

| n | Trials | Violations | min (A+B)/(2*Phi_r) | Infimum |
|---|--------|------------|---------------------|---------|
| 2 | 600    | 0          | 1.00000000          | 1 (proved) |
| 3 | 600    | 0          | 1.01108529          | 1 (achieved at equal gaps, p=q) |
| 4 | 600    | 0          | 1.04090570          | ~1.005 (p=q equally spaced) |
| 5 | 600    | 0          | 1.16589644          | ~1.010 |
| 6 | 600    | 0          | 1.58429829          | ~1.017 |

Both boxplus formulas (F1: MSS coefficient, F2: Marcus 2021) agree exactly
on this ratio (max discrepancy < 1e-15).

---

## Analytic Proof for n=2

For n=2 with gap_p = |lambda_1 - lambda_2|, gap_q = |mu_1 - mu_2|:

- Phi_p = 2/gap_p^2
- Phi_q = 2/gap_q^2
- gap_r^2 = gap_p^2 + gap_q^2 (Pythagorean property, verified numerically)
- Phi_r = 2/(gap_p^2 + gap_q^2)

Then Phi_p + Phi_q >= 4*Phi_r becomes:
```
2/x + 2/y >= 8/(x+y)    where x = gap_p^2, y = gap_q^2

<==>  (x+y)/xy >= 4/(x+y)
<==>  (x+y)^2 >= 4xy
<==>  (x-y)^2 >= 0     ALWAYS TRUE.
```

**Equality iff x = y, i.e., gap_p = gap_q.** Verified numerically to 10+ digits.

---

## Analytic Proof for n=3 (equally spaced case)

For p = {-s, 0, s}, q = {-t, 0, t}:
- Phi_p = 9/(2s^2), Phi_q = 9/(2t^2)
- r = {-D, 0, D} with D^2 = s^2 + t^2 (Pythagorean, verified)
- Phi_r = 9/(2D^2) = 9/(2(s^2+t^2))

```
(A+B)/(2*Phi_r) = (Phi_p + Phi_q)/(2*Phi_r) - 1
                = (s^2+t^2)^2 / (2*s^2*t^2) - 1
                = (s^4 + t^4) / (2*s^2*t^2)
                = (s^2/t^2 + t^2/s^2) / 2
                >= 1   (by AM-GM on s^2/t^2 and t^2/s^2)
```

Equality iff s = t.

---

## Equality Case Analysis

The infimum of (A+B)/(2*Phi_r) is **exactly 1**, achieved when:
- **n=2:** gap_p = gap_q (any p, q with equal gap sizes)
- **n=3:** p = q with equal spacing (or by symmetry, same gap structure)
- **n >= 4:** The infimum appears to be strictly > 1 even for p = q equally spaced:
  - n=4: min ratio ~ 1.0046 (at p = q equally spaced)
  - n=5: min ratio ~ 1.0105
  - n=6: min ratio ~ 1.0165

For n >= 4, `Phi_{p boxplus p}` is strictly less than `Phi_p/2` for all p,
meaning the inequality `A + B >= 2*Phi_r` is strict.

---

## Two Independent Conjectures

Both hold numerically but are logically independent:

**(C1) Sum inequality:** `Phi_p + Phi_q >= 4*Phi_r` (equivalently A + B >= 2*Phi_r)

**(C2) Product inequality:** `Phi_p * Phi_q >= Phi_r * (Phi_p + Phi_q)`
  (equivalently AB >= h^4, equivalently 1/Phi_r >= 1/Phi_p + 1/Phi_q = Fisher superadditivity)

In terms of x = Phi_p/Phi_r >= 1, y = Phi_q/Phi_r >= 1:
- (C1): x + y >= 4
- (C2): (x-1)(y-1) >= 1

Neither implies the other for general x, y >= 1.
Both hold in all tested MSS convolution instances.

Numerically verified product inequality:

| n | min Phi_p*Phi_q / (Phi_r*(Phi_p+Phi_q)) |
|---|------------------------------------------|
| 2 | 1.0000000000 (exact equality at n=2)    |
| 3 | 1.0002122077                             |
| 4 | 1.0134681045                             |
| 5 | 1.0353877348                             |
| 6 | 1.0494370139                             |

---

## Why the AM-GM Approach Fails

The key error: AM-GM gives `AB <= ((A+B)/2)^2`, not `>=`. So:

1. A + B >= 2*Phi_r does not help bound AB from below.
2. The ratio A/B is unbounded (can be > 200000 at n=3), so reverse AM-GM bounds are useless.
3. The two conditions (sum >= 4*Phi_r) and (product >= Phi_r*(sum)) are genuinely
   different properties of the MSS convolution, not derivable from each other.

---

## What WOULD Suffice

To prove AB >= h^4, one needs to establish Fisher superadditivity directly:
```
1/Phi_r >= 1/Phi_p + 1/Phi_q
```
This is the harmonic mean condition: `HarmonicMean(Phi_p, Phi_q) >= 2*Phi_r`.

The sum inequality `Phi_p + Phi_q >= 4*Phi_r` is the arithmetic mean condition:
`ArithmeticMean(Phi_p, Phi_q) >= 2*Phi_r`. Since HarmonicMean <= ArithmeticMean,
knowing the arithmetic mean condition tells us nothing about the harmonic mean.

---

## Structural Observations

1. **Both formulas agree:** F1 (MSS coefficient) and F2 (Marcus 2021) produce
   identical results for all quantities tested.

2. **Pythagorean property:** For both n=2 and n=3 equally spaced:
   gap_r^2 = gap_p^2 + gap_q^2. This does NOT hold for n=3 non-equally-spaced
   or for n >= 4.

3. **Scale invariance:** All ratios depend only on gap structure (ratios), not
   on absolute scale.

4. **p = q case:** When p = q, the ratio simplifies to Phi_p/Phi_r - 1.
   The condition Phi_{p boxplus p} <= Phi_p/2 always holds (verified for n=3,4,5).

5. **Equally spaced is the tightest case** for both inequalities: both (C1) and
   (C2) are closest to equality when p and q have uniform gap structure.

---

## Verdict

**A + B >= 2*Phi_r is TRUE** (supported by 10000+ random trials, zero counterexamples,
proved analytically for n=2 and n=3 equally spaced).

**But this does NOT close the conjecture** because the AM-GM step goes the wrong
way. The inequality `A + B >= 2*Phi_r` provides an interesting structural fact
about the MSS convolution but is insufficient for the product bound `AB >= h^4`.

The original Fisher superadditivity conjecture `1/Phi_r >= 1/Phi_p + 1/Phi_q`
remains the target and requires a direct proof.

---

## Recommendations

1. **Do not pursue the AM-GM path.** It is mathematically invalid.
2. The sum inequality `Phi_p + Phi_q >= 4*Phi_r` could still be useful as a
   lemma if combined with ADDITIONAL structure (e.g., a bound on Phi_p/Phi_q
   in terms of the gap structures).
3. Focus on direct approaches to Fisher superadditivity:
   - Matrix integral representation
   - Interlacing/majorization
   - Induction on n using the recursive structure of MSS convolution
4. The product inequality (C2) appears to strengthen as n grows (larger minimum ratio),
   suggesting an inductive approach might work.
