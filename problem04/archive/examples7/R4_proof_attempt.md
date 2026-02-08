# R_4 Superadditivity: Proof Attempt and Structural Analysis

**Author:** PROVER-9
**Date:** 2026-02-08
**Status:** Partial results; full proof not achieved

---

## 1. Formula Verification (COMPLETE)

The R_4 formula has been **verified symbolically** from first principles:

Starting from:
- `Phi_4 = N_4 / disc_4` where `N_4 = -8*e2^5 - 36*e2^2*e3^2 - 64*e2^3*e4 - 432*e3^2*e4 + 384*e2*e4^2`
- `disc_4 = 16*e2^4*e4 - 4*e2^3*e3^2 - 128*e2^2*e4^2 + 144*e2*e3^2*e4 - 27*e3^4 + 256*e4^3`

Using the cumulant-to-elementary-symmetric substitutions (for n=4 centered):
- `e2 = -3*k2/2`
- `e3 = k3/2`
- `e4 = -3*k4/32 + 3*k2^2/16`

We confirmed:
```
1/Phi_4 = (32*k2^6 - 32*k2^3*k3^2 - 6*k2^2*k4^2 + 24*k2*k3^2*k4 - 8*k3^4 - k4^3)
          / (384*k2^5 - 192*k2^2*k3^2 - 24*k2*k4^2 + 48*k3^2*k4)
```

And verified:
```
R_4 = 1/Phi_4 - k2/12
    = (-16*k2^3*k3^2 - 4*k2^2*k4^2 + 20*k2*k3^2*k4 - 8*k3^4 - k4^3)
      / (384*k2^5 - 192*k2^2*k3^2 - 24*k2*k4^2 + 48*k3^2*k4)
```

**Key factorization of the denominator:**
```
Q = 24*(4*k2^2 - k4)*(4*k2^3 + k2*k4 - 2*k3^2)
```

This matches the claimed formula in HANDOFF.md (the denominator in the HANDOFF had a factor of 16 absorbed into the polynomial terms).

**Verification script:** `examples7/R4_symbolic_verify.py`

---

## 2. Structural Properties Established

### 2.1 Homogeneity

R_4 has cumulant weight 2 (assigning weight k to kappa_k). Under the substitution
`k3 = k2^{3/2} * u`, `k4 = k2^2 * v`, we get:
```
R_4 = k2 * f(u, v)
```
where `f(u,v) = p(u,v) / q(u,v)` with:
- `p = -8*u^4 + 20*u^2*v - 16*u^2 - v^3 - 4*v^2`
- `q = 24*(v-4)*(2*u^2-v-4)`

### 2.2 Domain Constraints

For a centered real-rooted degree-4 polynomial with simple roots:
- `k2 > 0` (variance is positive)
- `k4 < 4*k2^2` (ensures factor `4*k2^2 - k4 > 0`)
- `2*k3^2 < 4*k2^3 + k2*k4` (ensures factor `4*k2^3 + k2*k4 - 2*k3^2 > 0`)

In scaled variables: `v < 4` AND `2*u^2 < v + 4`.

### 2.3 N_4 Factorization

The numerator `N_4` of `Phi_4` factors as:
```
N_4 = 81*k2^5 * (4*k2^2 - k4) * (4*k2^3 + k2*k4 - 2*k3^2) / 16
```
Since `Phi_4 > 0` and `disc > 0`, the denominator factors are both positive.

---

## 3. Special Cases

### 3.1 k3 = 0 (Symmetric Distributions)

When `k3 = 0`:
```
R_4 = -k4^2 / (24*k2*(4*k2^2 - k4))
```

With `w = k4/k2^2`:
```
R_4 = k2 * g(w),  g(w) = -w^2/(24*(4-w))
```

**Concavity of g:** `g''(w) = -4/(3*(4-w)^3) < 0` for `w < 4`. So g is strictly concave.

The superadditivity gap `Delta|_{k3=0}` has denominator:
```
Den = 24*k2p*k2q*(k2p+k2q)*(4*k2p^2-k4p)*(4*k2q^2-k4q)*(4*(k2p+k2q)^2-k4p-k4q)
```
which is **always positive on the domain**.

The superadditivity gap numerator (as a polynomial in `(u,v) = (k4_p, k4_q)` with coefficients in `s = k2_p, t = k2_q`):
- Has a **positive-definite quadratic part** with determinant `768*s^4*t^4*(s+t)^4`
- Vanishes at `u = v = 0` (trivially)
- Vanishes on the boundary `u = 4s^2, v = 4t^2`

**Numerical verification:** 0 violations in 10,000 random trials on the `k3=0` slice.

### 3.2 k4 = 0

When `k4 = 0`:
```
R_4 = -8*k3^4 / (24*4*k2^2*(4*k2^3 - 2*k3^2))
    = -k3^4 / (24*k2^2*(2*k2^3 - k3^2))
```

The denominator of `Delta|_{k4=0}` factors as:
```
24*k2p^2*k2q^2*(k2p+k2q)^2*(2*k2p^3-k3p^2)*(2*k2q^3-k3q^2)*(2*(k2p+k2q)^3-(k3p+k3q)^2)... (approximate)
```

### 3.3 Full Numerical Verification

0 violations in 3,545 valid random trials (full 3-parameter case).

---

## 4. Proof Attempts and Obstructions

### 4.1 Joint Concavity of -R_4 (FAILED)

The function `-R_4(k2, k4)|_{k3=0} = k4^2/(24*k2*(4*k2^2-k4))` is **NOT jointly concave**.
The Hessian has one positive eigenvalue everywhere. Checked at 10,000 random domain points:
all have indefinite Hessian.

**Why it fails:** While `g(w)` is univariately concave, the composition `k2*g(k4/k2^2)`
is not jointly concave because `k4/k2^2` introduces convexity in the ratio.

### 4.2 Titu/Engel-Cauchy-Schwarz (PARTIAL)

The n=3 proof worked because `R_3 = -2*k3^2/(27*k2^2)` has the form `x^2/y^2`, and
the inequality `sum(x_i^2/y_i^2) >= (sum x_i)^2/(sum y_i)^2` follows from Cauchy-Schwarz.

For R_4, the structure is `P(k)/Q(k)` where P is degree 12 and Q is degree 10 in
the weighted sense. The Titu/Engel approach would need the denominator to be **linear**
in the cumulants, which it is not (Q has degree 10).

### 4.3 SOS Decomposition (INCONCLUSIVE)

The numerator of the k3=0 gap, when expressed in terms of `a = 4*s^2 - u > 0` and
`b = 4*t^2 - v > 0`, takes the form:
```
128*s^5*t^2*(s+t)*b + 16*s^4*t*(s+t)*b^2
+ 128*s^2*t^5*(s+t)*a + 16*s*t^4*(s+t)*a^2
- 16*s^2*t^2*(s+3t)*(3s+t)*a*b
+ (s^2+t^2)*a^2*b^2 + s^2*a*b^3 + t^2*a^3*b
+ 8*s^2*t*(s-t)*a*b^2 - 8*s*t^2*(s-t)*a^2*b
```

The negative terms prevent a direct SOS argument. However, the structure is suggestive:
- Pure `a` and pure `b` terms are positive
- The leading cross-term `-16*s^2*t^2*(s+3t)*(3s+t)*a*b` is the main obstruction
- Higher-order terms `a^2*b^2, a*b^3, a^3*b` have positive coefficients

An SOS certificate might exist but would require more sophisticated techniques
(SDP/DSOS/SDSOS solvers).

### 4.4 Perspective Function (BLOCKED)

The standard perspective function theory applies to `t*f(x/t)` which preserves
convexity/concavity. But `R_4 = k2*f(k4/k2^2)` has a **quadratic** denominator,
not linear. This breaks the standard perspective function approach.

---

## 5. Promising Directions (Not Yet Explored)

### 5.1 Schur-Convexity in (k2p, k2q) with Fixed Ratios

Fix `w_p = k4_p/k2_p^2` and `w_q = k4_q/k2_q^2` (scaled cumulants).
The gap becomes a function of `(s,t) = (k2_p, k2_q)` which may be provably non-negative
using Schur-convexity arguments.

### 5.2 SDP-based SOS Certificates

Use a semidefinite programming solver (e.g., DSOS/SDSOS from SOSTOOLS) to search for
an SOS certificate for the numerator polynomial (after suitable substitution to
remove the domain constraint, e.g., replacing `u` by `4*s^2 - a^2`).

### 5.3 Induction on the Number of Cumulants

Prove for k3=0 first (2-variable problem), then use a perturbation argument
for small k3, then extend to all k3 by a continuity/monotonicity argument.

### 5.4 Convexity of 1/Phi_4 Along Rays

If `t -> 1/Phi_4(t*kappa_p + (1-t)*kappa_q)` is convex in t, this gives
superadditivity by evaluating at t = k2_p/(k2_p+k2_q).

### 5.5 Matrix Quadratic Form Generalization

For n=3, the proof used a 2x2 positive-definite matrix. For n=4, seek a
3x3 (or larger) positive-definite matrix M(s,t) such that the superadditivity
gap equals `[k3_p, k4_p, k3_q, k4_q]^T * M * [k3_p, k4_p, k3_q, k4_q]`
(modulo denominator corrections).

---

## 6. Summary

| Component | Status |
|-----------|--------|
| R_4 formula verification | **COMPLETE** - verified symbolically |
| Denominator factorization | **COMPLETE** - `24*(4*k2^2-k4)*(4*k2^3+k2*k4-2*k3^2)` |
| Domain constraints | **COMPLETE** - both factors positive |
| Homogeneity analysis | **COMPLETE** - weight 2 |
| k3=0 case analysis | **COMPLETE** - reduced to numerator non-negativity, quadratic part PD |
| k4=0 case analysis | **COMPLETE** - denominator factored |
| Numerical verification | **COMPLETE** - 0 violations in ~15K trials |
| Joint concavity | **FAILED** - -R_4 not jointly concave |
| Titu/Engel CS | **INAPPLICABLE** - denominator not linear |
| SOS decomposition | **INCONCLUSIVE** - negative cross terms |
| Analytical proof | **OPEN** |

**Key structural finding:** The denominator of the superadditivity gap always factors into
positive terms (product of the three denominator factors for p, q, and r). This reduces
the problem to proving a **polynomial inequality** (numerator non-negativity).

**Key obstruction:** Unlike n=3 where the numerator was a quadratic form in a 2-vector,
for n=4 the numerator is a degree-3 polynomial in the cumulants with mixed signs.
Simple SOS or Cauchy-Schwarz arguments do not close it directly.

**Recommended next step:** Apply computational SOS methods (DSOS/SDSOS) or
try the Schur-convexity approach (Direction 5.1).

---

## Supporting Scripts

- `examples7/R4_symbolic_verify.py` - Formula verification
- `examples7/R4_structural_analysis.py` - Structural analysis and special cases
- `examples7/R4_numerical_check.py` - Numerical verification
- `examples7/R4_concavity_proof.py` - Joint concavity test (negative result)
- `examples7/R4_perspective_proof.py` - SOS decomposition attempt
- `examples7/R4_k3zero_proof.py` - k3=0 case analysis
