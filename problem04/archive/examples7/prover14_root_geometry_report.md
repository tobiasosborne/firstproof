# PROVER-14: Root Geometry Report on Fisher Superadditivity

**Author:** PROVER-14
**Date:** 2026-02-08
**Status:** Major structural results proved; full proof reduced to a clean inequality

---

## Executive Summary

This report establishes several key structural results about the Fisher superadditivity conjecture using root geometry. The main achievements are:

1. **Proved** the identity `Phi_n = 2 * sum_{i<j} 1/(lambda_i - lambda_j)^2` for all n.
2. **Reformulated** the conjecture as superadditivity of the inverse-square interaction energy.
3. **Derived** a clean trigonometric parametrization for n=3 that reduces the conjecture to a single inequality about a function G(psi) = F(psi/3).
4. **Numerically verified** G is **convex in cos(psi)**, which if proved would give a route to the full n=3 proof.
5. **Found** the symbolic formula G(c) = -9/(16c^6 - 24c^4 + 9c^2 - 1) where c = cos(phi).

---

## 1. Key Identity: Phi_n = 2 * Sm2

**THEOREM (Proved analytically and symbolically for all n):**

```
Phi_n(p) = sum_i H_p(lambda_i)^2 = 2 * sum_{i<j} 1/(lambda_i - lambda_j)^2
```

where `lambda_1 < ... < lambda_n` are the roots of p.

**Proof:** Expand `Phi_n = sum_i (sum_{j!=i} 1/d_{ij})^2` into diagonal terms `sum_{i!=j} 1/d_{ij}^2 = 2*Sm2` and cross terms. The cross terms, summed over all unordered triples {i,j,k}, give zero because for each triple:

```
1/((a-b)(a-c)) + 1/((b-a)(b-c)) + 1/((c-a)(c-b)) = 0
```

This follows from putting over a common denominator: the numerator is `(b-c) - (a-c) + (a-b) = 0`. QED.

**Verified:** Symbolically for n=3,4 and numerically for n=2,...,7.

---

## 2. S2 Additivity Under MSS

**FACT:** `S2(p boxplus_n q) = S2(p) + S2(q)` where `S2 = sum_{i<j} (lambda_i - lambda_j)^2`.

This holds because `S2 = n * (n-1) * k_2` (proportional to the second finite free cumulant), and cumulants are additive under MSS convolution.

**Verified:** Numerically to machine precision for n=2,3,4,5.

---

## 3. Reformulated Conjecture

Using `Phi_n = 2*Sm2`, the conjecture `1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q)` becomes:

```
1/Sm2(r) >= 1/Sm2(p) + 1/Sm2(q)
```

where `Sm2 = sum_{i<j} 1/(lambda_i - lambda_j)^2` is the **total inverse-square interaction energy**.

**Physical interpretation:** Sm2 is the total potential energy of unit charges on a line interacting via the inverse-square potential V(r) = 1/r^2. The conjecture says this energy satisfies the **harmonic mean bound** under MSS convolution:

```
Sm2(r) <= Sm2(p) * Sm2(q) / (Sm2(p) + Sm2(q))
```

---

## 4. Universal Identities

**Proved:**
- `sum_i H_p(lambda_i) = 0` (well-known)
- `sum_i H_p(lambda_i) * lambda_i = n(n-1)/2` (universal, independent of roots)
- `sum_i H_p(lambda_i) * lambda_i^2 = (n-1) * sum_j lambda_j` (depends only on mean)

---

## 5. Scale-Invariant Shape Factor

Define `Q = Sm2 * S2` (dimensionless, scale-invariant).

**Properties:**
- For n=2: Q = 1 always (trivial).
- For n=3: Q depends only on the gap ratio r = t/s, with minimum Q = 13.5 at r = 1 (equal gaps).
- Q is Schur-convex in the gap vector (more equal gaps = smaller Q).
- **Numerically verified:** Q_r <= max(Q_p, Q_q) under MSS convolution (0 violations in 5000 trials).

---

## 6. n=2: Exact Equality

For n=2: `Sm2 = 1/d^2`, `S2 = d^2`, so `1/Sm2 = d^2 = S2`.

By S2 additivity: `1/Sm2(r) = S2(r) = S2(p) + S2(q) = 1/Sm2(p) + 1/Sm2(q)`. **Equality!**

This is equivalent to the **Pythagorean theorem** on root gaps: `gap_r^2 = gap_p^2 + gap_q^2`.

---

## 7. Trigonometric Parametrization for n=3

For centered n=3 polynomials `p(x) = x^3 + sigma*x + tau`:

- MSS convolution ADDS (sigma, tau) coordinate-wise.
- Set `a = -sigma > 0`, and `phi in (0, pi/3)` with `cos(3*phi) = -3*sqrt(3)*tau / (2*a^{3/2})`.
- Roots: `2*sqrt(a/3) * cos(phi + 2k*pi/3)` for k=0,1,2.

Then:
```
Sm2 = F(phi) / (12*a/3) = F(phi) / (4a)
1/Sm2 = 4a / F(phi)
```
where:
```
F(phi) = csc^2(phi) + csc^2(phi + pi/3) + csc^2(phi + 2*pi/3)
```

**Symbolic formula:** `F(phi) = G(cos(phi))` where `G(c) = -9 / (16c^6 - 24c^4 + 9c^2 - 1)`.

Equivalently, `Q = 3*F(phi)/2`, so F determines the shape factor.

### MSS Addition in (a, phi) Coordinates

```
a_r = a_p + a_q
cos(3*phi_r) = [a_p^{3/2} * cos(3*phi_p) + a_q^{3/2} * cos(3*phi_q)] / (a_p + a_q)^{3/2}
```

**Key contraction property:** Since `a_p^{3/2} + a_q^{3/2} < (a_p + a_q)^{3/2}` (strict concavity of x^{3/2}), the formula for cos(3*phi_r) is a **contracted average** -- it pushes phi_r toward pi/6 (equal gaps).

### Reduced Conjecture for n=3

The conjecture reduces to:

```
F(phi_r) <= (a_p + a_q) / (a_p/F(phi_p) + a_q/F(phi_q))
```

i.e., F(phi_r) is bounded by the **weighted harmonic mean** of F(phi_p) and F(phi_q) with weights a_p, a_q.

**Numerically verified:** 0 violations in 100,000 trials.

---

## 8. Convexity of G in cos(psi)

Let `psi = 3*phi`, so `psi in (0, pi)`. Define `G(psi) = F(psi/3)`.

**Numerically verified:** G is **convex in cos(psi)**.

If proved, this convexity combined with the contraction property of cos(3*phi_r) would give a route to the n=3 proof, because:
- cos(psi_r) is a contracted weighted average of cos(psi_p) and cos(psi_q)
- G convex in cos means G(phi_r) <= appropriate combination of G(phi_p) and G(phi_q)

However, the exact argument requires careful handling of the different weights (a^{3/2} in the cos formula vs a in the harmonic mean).

---

## 9. Near-Equality Analysis

The gap `1/Sm2(r) - 1/Sm2(p) - 1/Sm2(q)` approaches zero when:
- phi_p and phi_q are both close to pi/6 (both p and q have nearly equal gaps), OR
- One of a_p, a_q is very small (one polynomial dominates)

---

## 10. Key Findings Not Yet Exploited

1. **Q_r <= max(Q_p, Q_q):** The shape factor Q monotonically decreases under MSS convolution. This is a STRONG structural property but does not directly imply the conjecture.

2. **Hessian structure:** The Hessian of the log-potential W (where H = -grad W) is positive semidefinite on the hyperplane sum=0, making W convex in the root configuration.

3. **Critical point representation:** `H_p(lambda_i) = (1/2) sum_j 1/(lambda_i - mu_j')` where mu_j' are critical points of p. Combined with ADMITTED C (MSS respects the critical-point tower), this gives a recursive structure.

4. **Cauchy matrix:** The matrix `C_{ij} = 1/(lambda_i - mu_j')` is a Cauchy matrix with known determinant and eigenvalue properties. Phi_n = (1/4)||C * 1||^2.

---

## 11. Recommended Next Steps

1. **Prove G convex in cos(psi)** for n=3. This is a single-variable calculus problem: show `d^2G/du^2 >= 0` where `u = cos(psi)` and `G(c) = -9/(16c^6 - 24c^4 + 9c^2 - 1)` with `c = cos(psi/3)`.

2. **Complete the n=3 proof** using convexity of G plus the contraction property. The gap between "convexity in cos" and "harmonic mean bound with different weights" needs to be bridged.

3. **Generalize to n=4:** Find the analogous trigonometric parametrization and shape factor. For n=4, there are 2 angular parameters instead of 1.

4. **Prove Q_r <= max(Q_p, Q_q)** as a standalone result about MSS convolution. This would be interesting in its own right.

---

## 12. Scripts

| Script | Purpose |
|--------|---------|
| `R4_prover14_roots.py` | Initial root geometry exploration, MSS implementation, basic verification |
| `R4_prover14_electrostatic.py` | Electrostatic interpretation, Hessian analysis, subordination |
| `R4_prover14_cauchy_schwarz.py` | Cauchy-Schwarz approaches, n=2 equality proof, Cauchy matrix |
| `R4_prover14_S2_additive.py` | S2 additivity, Phi_n=2*Sm2 discovery, Q=Sm2*S2 analysis |
| `R4_prover14_key_identity.py` | Formal proof of Phi_n=2*Sm2, reformulated conjecture |
| `R4_prover14_Q_bound.py` | Q-bound analysis, weighted harmonic mean, n=3 symbolic MSS |
| `R4_prover14_n3_proof.py` | Trigonometric parametrization, G(c) formula, reduced conjecture |

---

## Summary Table

| Result | Status |
|--------|--------|
| Phi_n = 2*Sm2 identity | **PROVED** (all n) |
| S2 additivity under MSS | **VERIFIED** (known from cumulant theory) |
| n=2 equality | **PROVED** |
| Q_r <= max(Q_p, Q_q) | **NUMERICALLY VERIFIED** (5000 trials) |
| G(c) = -9/(16c^6-24c^4+9c^2-1) | **DERIVED** |
| G convex in cos(psi) | **NUMERICALLY VERIFIED** |
| n=3 reduced conjecture | **FORMULATED** (0 violations in 100K trials) |
| Full n=3 proof | **OPEN** (reduced to convexity argument) |
| Full n >= 4 proof | **OPEN** |
