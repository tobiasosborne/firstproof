# Findings: Algebraic Identity for AB - ||h||^4

**Node:** 1.10.2 (algebraic identity / manifestly non-negative expression)
**Agent:** algebraic-prover
**Date:** 2026-02-08
**Script:** `investigate_algebraic_identity.py`

---

## Summary of Findings

### The Core Identity

AB - ||h||^4 admits the following clean factorization:

```
AB - h^4 = Phi_p * Phi_q * Phi_r * (1/Phi_r - 1/Phi_p - 1/Phi_q)
```

where Phi_f = ||H_f||^2 is the finite Fisher information of polynomial f.

**Equivalently:** AB - h^4 = ||u||^2 * ||v||^2 - ||h||^2 * (||u||^2 + ||v||^2)

So **AB >= h^4** is exactly equivalent to **1/Phi_r >= 1/Phi_p + 1/Phi_q** (Fisher superadditivity).

This is a tautology -- it reformulates the target inequality in clean notation but does not prove it.

---

### Key Structural Finding: AM-GM Reduction

**MAIN RESULT:** A stronger inequality appears to hold:

```
A + B >= 2 * Phi_r     (i.e., Phi_p + Phi_q >= 4 * Phi_r)
```

This was verified numerically for 1000 random trials each at n = 2, 3, 4, 5, 6 with **zero counterexamples**. If this can be proved, then by AM-GM:

```
AB >= ((A + B) / 2)^2 >= Phi_r^2 = h^4
```

This reduces the product inequality AB >= h^4 to the **SUM inequality** A + B >= 2*Phi_r, which may be more tractable.

**Status:** NOT PROVED. This is a conjecture supported by extensive numerics.

---

### Equality Cases

1. **n=2:** AB = h^4 exactly, always. Proved algebraically:
   - For p = (a,b), q = (c,d): D^2 = (a-b)^2 + (c-d)^2
   - AB = 4/D^4 = h^4 (perfect cancellation, Pythagorean structure)

2. **n=3, equally spaced:** AB = h^4 exactly.
   - For p = (-s, 0, s), q = (-t, 0, t): r = (-D, 0, D) where D = sqrt(s^2 + t^2)
   - Phi_p = 9/(2s^2), Phi_q = 9/(2t^2), Phi_r = 9/(2D^2)
   - AB = 81/(4D^4) = h^4 (same Pythagorean structure as n=2)

3. **n=3, unequal gaps:** AB > h^4 strictly.
   - Minimum AB/h^4 over all n=3 configurations: exactly 1, achieved at equal spacing.
   - The excess AB - h^4 comes from perpendicular components of alpha, beta relative to h.

---

### Gram Matrix Analysis

The 3x3 Gram matrix G = [<e_i, e_j>] for {h, alpha, beta}:

```
G = [[h^2,      <h,a>,     <h,b>  ],
     [<h,a>,    ||a||^2,   <a,b>  ],
     [<h,b>,    <a,b>,     ||b||^2]]
```

**Critical finding for n=3:** det(G) = 0 always (to machine precision, ~1e-15).
This means h, alpha, beta are **always linearly dependent** for n=3.
Specifically, h lies in span(alpha, beta).

For n >= 4, h, alpha, beta are generically independent (det(G) != 0).

---

### The <h, alpha> Sign Issue

**CONFIRMED:** <h, alpha> CAN be negative for n >= 3.
- n=3: negative in 15/2000 random trials, worst value = -1.164
- n=4: negative in 15/2000 random trials, worst value = -2.997
- n=5: negative in 5/2000 random trials, worst value = -1.620

Despite <h, alpha> being negative, AB - h^4 remains non-negative in all tested cases.
**The proof CANNOT assume <h, alpha> >= 0.**

---

### Decomposition: Parallel vs Perpendicular

Using the notation:
- p = <h, alpha> / ||h||^2, q = <h, beta> / ||h||^2
- s = ||alpha_perp||^2 / ||h||^2, t = ||beta_perp||^2 / ||h||^2

Then AB/h^4 = (p(p+2) + s)(q(q+2) + t) = parallel + cross + perp

Where:
- parallel = p(p+2) * q(q+2)  -- can be < 1 (in 13/2000 cases)
- cross = p(p+2)*t + q(q+2)*s -- always >= 0 when parallel < 1 (compensates)
- perp = s*t >= 0

The parallel term alone does NOT suffice; the perpendicular contributions are essential.

---

### Approaches That Do NOT Work

1. **<h, alpha> >= 0 assumption:** FALSE (counterexamples at n >= 3)
2. **min(A, B) >= h^2:** FALSE (holds only in 956/2000 n=3 cases)
3. **Cauchy-Schwarz on <u,v>:** Gives AB >= <u,v>^2 - h^2*(u^2+v^2), which is not useful
4. **Lagrange identity decomposition:** AB - h^4 = <u,v>^2 + ||u x v||^2 - <h,u>^2 - ||h x u||^2 - <h,v>^2 - ||h x v||^2, but individual terms can have either sign
5. **Gram determinant approach:** det(Gram(h,u,v)) >= 0 is TRUE but gives a different quantity from AB - h^4

---

### Additional Structural Observations

1. **Scale invariance:** AB/h^4 depends only on the SHAPE of root configuration (gap ratios), not on overall scale.

2. **Translation invariance:** AB/h^4 is unchanged under shifting all roots by a constant.

3. **Antisymmetric matrix representation:** h = M_r * e where M_r[k,l] = 1/(nu_k - nu_l) is antisymmetric and e = (1,...,1). Phi_r = -e^T M_r^2 e.

4. **sum(nu_k) = sum(lambda_k) + sum(mu_k):** The MSS convolution preserves trace (e_1 is additive).

5. **lambda_k + mu_k - nu_k is NOT constant** across k (the subordination identity omega_1 + omega_2 = id + c does not hold in finite n).

---

### Recommended Proof Strategies (Ordered by Promise)

**Strategy 1 (AM-GM reduction):** Prove A + B >= 2*Phi_r, then AB >= h^4 by AM-GM.
Reformulated: prove Phi_p + Phi_q >= 4*Phi_r.
This is equivalent to: (1/F_p + 1/F_q) >= 4/F_r where F = 1/Phi.
By the target inequality F_r >= F_p + F_q, we get 1/F_r <= 1/(F_p + F_q).
So 4/F_r <= 4/(F_p + F_q). And 1/F_p + 1/F_q >= (F_p + F_q)/(F_p*F_q) = ... this is CIRCULAR.
**Verdict:** Cannot derive A + B >= 2h2 from the target. This is a SEPARATE inequality that may need its own proof.

**Strategy 2 (Convexity/majorization):** Use the specific structure of MSS convolution to bound Phi_r in terms of Phi_p, Phi_q. The roots of r are "more spread" than those of p,q individually, which forces Phi_r to be smaller.

**Strategy 3 (n=3 direct):** For n=3, use the linear dependence h = a*alpha + b*beta to reduce AB - h^4 to a 2D problem. The constraints on a, b from the MSS structure may directly force the inequality.

---

### Gaps and Honest Assessment

1. **No manifestly non-negative expression found** for AB - h^4 in terms of h, alpha, beta, or u, v, h.

2. The AM-GM path (A+B >= 2h^2) is a CONJECTURE, not a theorem. It would suffice but has not been proved.

3. The linear dependence of {h, alpha, beta} for n=3 is a strong structural fact, but it is not clear how to exploit it for the proof.

4. The key difficulty: AB - h^4 involves a PRODUCT of two differences (Phi_p - Phi_r)(Phi_q - Phi_r), and comparing this product to Phi_r^2 requires understanding the CORRELATION between how much p and q individually spread relative to r. This correlation comes from the MSS convolution structure, which couples p and q.

5. For n=2, the proof is trivial (Pythagorean). For n >= 3, the proof genuinely requires properties of the MSS convolution that go beyond simple algebraic manipulation of h, alpha, beta.
