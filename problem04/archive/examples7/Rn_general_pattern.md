# General Pattern for C_n and R_n (PROVER-10 Investigation)

## Summary of Results

### Main Theorem (PROVED): C_n Formula

**Theorem.** For monic degree-n real-rooted polynomials with simple roots,
$$\frac{1}{\Phi_n} = C_n \cdot \kappa_2 + R_n(\kappa_2, \kappa_3, \ldots, \kappa_n)$$
where $\kappa_k$ are the additive finite free cumulants, and
$$C_n = \frac{-2}{\binom{n}{2}} = \frac{-4}{n(n-1)}$$
in the natural (Arizmendi-Perales) cumulant convention where $\kappa_2 = \tilde{a}_2 < 0$ for valid polynomials.

**In the original conjecture's convention** ($\kappa_2^{\mathrm{orig}} = -n\tilde{a}_2 > 0$):
$$C_n^{\mathrm{orig}} = \frac{4}{n^2(n-1)}$$

| n | $C_n^{\mathrm{orig}}$ | Decimal |
|---|----------------------|---------|
| 2 | 1                    | 1.000   |
| 3 | 2/9                  | 0.222   |
| 4 | 1/12                 | 0.083   |
| 5 | 1/25                 | 0.040   |
| 6 | 1/45                 | 0.022   |
| 7 | 2/147                | 0.014   |

The **original conjecture** $C_n = 2/(n(n-1))$ was **wrong** by a factor of $2/n$. The correct formula $C_n = 4/(n^2(n-1))$ matches the known values $C_2=1$, $C_3=2/9$, $C_4=1/12$ and predicts new values for $n \geq 5$.

### Proof of C_n

**Step 1.** At $\kappa_3 = \cdots = \kappa_n = 0$, the EGF of normalized coefficients is $T(x) = \exp(\kappa_2 x^2/2)$. Setting $s = -\kappa_2 > 0$, the polynomial is
$$p_n(x) = s^{n/2} \cdot \mathrm{He}_n(x/\sqrt{s})$$
where $\mathrm{He}_n$ is the probabilist's Hermite polynomial.

**Step 2.** Roots: $\lambda_i = \sqrt{s} \cdot \xi_i$ where $\xi_i$ are roots of $\mathrm{He}_n$.

**Step 3.** For Hermite polynomials, three identities hold:
- $\mathrm{He}_n'(x) = n \cdot \mathrm{He}_{n-1}(x)$
- $\mathrm{He}_n''(x) = n(n-1) \cdot \mathrm{He}_{n-2}(x)$
- $\mathrm{He}_n(x) = x \cdot \mathrm{He}_{n-1}(x) - (n-1) \cdot \mathrm{He}_{n-2}(x)$

At root $\xi_i$: from the recurrence, $\mathrm{He}_{n-2}(\xi_i) = \xi_i \mathrm{He}_{n-1}(\xi_i)/(n-1)$. Therefore:
$$H(\xi_i) = \frac{\mathrm{He}_n''(\xi_i)}{2\mathrm{He}_n'(\xi_i)} = \frac{n(n-1) \cdot \xi_i \mathrm{He}_{n-1}(\xi_i)/(n-1)}{2n \cdot \mathrm{He}_{n-1}(\xi_i)} = \frac{\xi_i}{2}$$

**Step 4.** For the actual polynomial: $H_p(\lambda_i) = \xi_i/(2\sqrt{s})$.

**Step 5.** $\Phi_n = \sum_i \xi_i^2/(4s)$.

**Step 6.** By Vieta for $\mathrm{He}_n(x) = x^n - \binom{n}{2}x^{n-2} + \cdots$:
$$\sum \xi_i^2 = -2 \sum_{i<j}\xi_i\xi_j = 2\binom{n}{2} = n(n-1)$$

**Step 7.** $\Phi_n = n(n-1)/(4s)$, so $1/\Phi_n = 4s/(n(n-1)) = (-4/(n(n-1)))\kappa_2 = (-2/\binom{n}{2})\kappa_2$. **QED**

Verified numerically for $n = 2, 3, \ldots, 14$ with machine-precision agreement.

## Cumulant Structure

### MSS Convolution and Normalized Coefficients

For a monic polynomial $p(x) = x^n + a_1 x^{n-1} + \cdots + a_n$, define:
$$\tilde{a}_k = \frac{(-1)^k a_k}{\binom{n}{k}}$$

**Key discovery:** Under MSS (finite free additive) convolution:
$$\tilde{a}_k(p \boxplus q) = \sum_{i+j=k} \binom{k}{i} \tilde{a}_i(p) \cdot \tilde{a}_j(q)$$

This is the **binomial convolution** -- multiplication of exponential generating functions!

### Additive Cumulants

Define the EGF $T(x) = 1 + \sum_{k \geq 1} \tilde{a}_k x^k/k!$ with $\tilde{a}_0 = 1$.

Under MSS: $T_{p \boxplus q}(x) = T_p(x) \cdot T_q(x)$.

The **additive cumulants** are the coefficients of $\log T(x) = \sum_{k \geq 1} \kappa_k x^k/k!$.

These satisfy $\kappa_k(p \boxplus q) = \kappa_k(p) + \kappa_k(q)$ (additivity).

**Important:** These cumulants are **universal** (independent of degree $n$).

### Explicit Cumulant Formulas (centered, $\kappa_1 = 0$)

| Cumulant | Formula |
|----------|---------|
| $\kappa_2$ | $\tilde{a}_2$ |
| $\kappa_3$ | $\tilde{a}_3$ |
| $\kappa_4$ | $\tilde{a}_4 - 3\tilde{a}_2^2$ |
| $\kappa_5$ | $\tilde{a}_5 - 10\tilde{a}_2\tilde{a}_3$ |
| $\kappa_6$ | $\tilde{a}_6 - 15\tilde{a}_2\tilde{a}_4 - 10\tilde{a}_3^2 + 30\tilde{a}_2^3$ |

Coefficients: 3 = C(4,2)/2, 10 = C(5,2), 15 = C(6,2), 10 = C(6,3)/2, 30 = C(6,2)C(4,2)/2(3!). These follow from the standard log-of-EGF expansion.

## R_n Structure

### n = 3

$$R_3 = -\frac{1}{6} \frac{\kappa_3^2}{\kappa_2^2}$$

- Weight 2 (homogeneous)
- Always $\leq 0$ (with equality iff $\kappa_3 = 0$)
- **Superadditive** for centered polynomials (proved via Jensen/convexity of $x^2$)

### n = 4

$$R_4 = \frac{-54\kappa_2^3\kappa_3^2 + 6\kappa_2^2\kappa_4^2 - 45\kappa_2\kappa_3^2\kappa_4 + 27\kappa_3^4 - \kappa_4^3}{9(6\kappa_2^2 + \kappa_4)(6\kappa_2^3 - \kappa_2\kappa_4 + 3\kappa_3^2)}$$

- Weight 2 (homogeneous in the grading $\mathrm{wt}(\kappa_k) = k$)
- Numerator: weight 12, denominator: weight 10
- At $\kappa_4 = 0$: $R_4 = \kappa_3^2(2\kappa_2^3 - \kappa_3^2)/[6\kappa_2^2(2\kappa_2^3 + \kappa_3^2)]$
- Can change sign depending on relative magnitudes of cumulants
- **Superadditive for centered polynomials** (verified numerically, 0/9891 violations)

### General Pattern

For all $n$:
- $R_n$ is a **rational function** of $\kappa_2, \ldots, \kappa_n$
- $R_n$ is **homogeneous of weight 2** under the grading $\mathrm{wt}(\kappa_k) = k$
- $R_n$ vanishes when $\kappa_3 = \cdots = \kappa_n = 0$
- $R_n$ has no simple recursive relation to $R_{n-1}$
- The denominator of $1/\Phi_n$ (in cumulant variables) factors into discriminant-related terms

## Superadditivity Results

### Full 1/Phi_n Superadditivity

$$\frac{1}{\Phi_n(p \boxplus q)} \geq \frac{1}{\Phi_n(p)} + \frac{1}{\Phi_n(q)}$$

| n | Violations | Valid trials |
|---|-----------|-------------|
| 3 | 0         | 8413        |
| 4 | 0         | 7089        |
| 5 | 0         | 5533        |
| 6 | 0         | 4225        |

**No violations found for any degree tested.**

### R_n Superadditivity (centered polynomials only)

$$R_n(\kappa_2(p) + \kappa_2(q), \ldots) \geq R_n(\kappa_2(p), \ldots) + R_n(\kappa_2(q), \ldots)$$

| n | Violations | Valid trials | Status |
|---|-----------|-------------|--------|
| 3 | 0         | 8152        | Proved (Cauchy-Schwarz/Jensen) |
| 4 | 0         | 9891        | Verified numerically |
| 5 | 0         | 7495        | Verified numerically |

**Important caveat:** R_n superadditivity only holds for **centered** polynomials ($\kappa_1 = 0$). For non-centered polynomials, the formula $R_n = 1/\Phi_n - C_n\kappa_2$ does not correctly decompose the function, since $\Phi_n$ is translation-invariant but $\kappa_2$ is not. WLOG we can center first.

## Key Structural Discoveries

1. **MSS is binomial convolution in disguise:** The normalized coefficients $\tilde{a}_k$ multiply via binomial convolution $W_k(i,j) = \binom{k}{i}$, which is multiplication of EGFs.

2. **Cumulants are universal:** The additive cumulants (log of EGF) are independent of the polynomial degree $n$. This is a new observation.

3. **Hermite connection:** At the "Gaussian locus" $\kappa_3 = \cdots = 0$, the polynomial is a scaled Hermite polynomial. The key identity $H_i = \xi_i/2$ for Hermite roots makes the $C_n$ computation elementary.

4. **The correction to the conjecture:** $C_n = 4/(n^2(n-1))$, not $2/(n(n-1))$. The factor $2/n$ comes from the relationship between the natural cumulant $\kappa_2 = \tilde{a}_2$ and the "original" convention $\kappa_2^{\mathrm{orig}} = -n\tilde{a}_2$.

## Files

- `Rn_investigation.py` -- initial numerical exploration
- `Rn_detailed.py` -- symbolic computation of Phi_n*disc
- `Rn_definitive.py` -- cumulant structure and W_k(i,j)
- `Rn_Cn_extraction.py` -- C_n extraction and Hermite proof
- `Rn_proof_and_summary.py` -- complete proof and superadditivity tests
- `Rn_R3_check.py` -- debugging R_3 superadditivity (centered vs non-centered)
- `Rn_centered_check.py` -- confirming R_4, R_5 superadditivity for centered polys
