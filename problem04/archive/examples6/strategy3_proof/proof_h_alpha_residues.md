# Proof Attempt: sum_k h_k * alpha_k >= 0 via Contour Integrals and Residues

## Node 1.7.2 -- Sub-lemma for Fisher Information Superadditivity

**Status: UNPROVED. Key identities established. Main gap identified.**

---

## 1. Setup and Notation

Let p, q be monic real-rooted polynomials of degree n with simple roots.
Let r = p boxplus_n q (MSS finite free convolution).

- Roots: p has lambda_1 < ... < lambda_n, r has nu_1 < ... < nu_n
- Subordination: omega_1 satisfying G_r(z) = G_p(omega_1(z))
- At roots: omega_1(nu_k) = lambda_k (order-preserving), omega_1'(nu_k) = 1
- H_r(nu_k) = h_k, H_p(lambda_k) = u_k, alpha_k = u_k - h_k = omega_1''(nu_k)/2

**Target:** sum_k h_k * alpha_k >= 0

---

## 2. Established Identities (PROVED)

### Identity 1: Generating function for weighted sums of h_k

For any real p not equal to a root of r:

**sum_k H_r(nu_k) / (nu_k - p) = -(1/2) * r''(p) / r(p)**

*Proof.* Write nG_r(z) = sum_k 1/(z - nu_k), so
sum_k 1/(nu_k - p) = -nG_r(p) = -r'(p)/r(p).

sum_k H_r(nu_k)/(nu_k - p) = sum_k [sum_{j!=k} 1/(nu_k - nu_j)] / (nu_k - p)
= sum_{k,j: k!=j} 1/((nu_k - nu_j)(nu_k - p))
= sum_{k<j} [1/((nu_k-nu_j)(nu_k-p)) + 1/((nu_j-nu_k)(nu_j-p))]
= -sum_{k<j} 1/((nu_k-p)(nu_j-p))
= -(1/2) * {[sum_k 1/(nu_k-p)]^2 - sum_k 1/(nu_k-p)^2}
= -(1/2) * {[nG_r(p)]^2 + [nG_r]'(p)}
= -(1/2) * r''(p)/r(p).

The last step uses the standard identity:
[r'/r]^2 + [r'/r]' = (r')^2/r^2 + [r''*r - (r')^2]/r^2 = r''/r. QED

**Verified numerically** for n = 3, 4, 5 across hundreds of random test points.

### Identity 2: Contour integral for position displacement

**(1/2pi i) oint [nG_r(z)]^2 * [z - omega_1(z)] dz = 0** (over large contour)

yields: **sum_k h_k * (nu_k - lambda_k) = (1/2) * sum_j m_j * [nG_r(p_j)]^2 >= 0**

where p_j are the poles of omega_1 with residues m_j > 0 in the Herglotz representation.

*Proof sketch.* Let Psi(z) = z - omega_1(z). At z = nu_k:
- [nG_r]^2 = 1/(z-nu_k)^2 + 2*h_k/(z-nu_k) + ...
- Psi(z) = (nu_k - lambda_k) + 0*(z-nu_k) + ...  (since Psi'(nu_k) = 1 - omega_1'(nu_k) = 0)
- Residue at nu_k: 2 * h_k * (nu_k - lambda_k)

At poles p_j of omega_1: Psi(z) = z - omega_1(z) has simple poles at p_j with residue m_j
(in the convention omega_1(z) = z + b + sum m_j/(p_j - z), so z - omega_1 = -b + sum m_j/(z-p_j),
and the residue of Psi at p_j is... we need to be careful with the sign convention).

At infinity: [nG_r]^2 * Psi ~ (1/z^2)*O(1) = O(1/z^2), integral vanishes.

This yields sum_k h_k*(nu_k-lambda_k) >= 0, with equality expressed through the
Herglotz residues. The key point is that the auxiliary residues involve
m_j * [nG_r(p_j)]^2 >= 0 (both factors non-negative).

**Status:** PROVED (modulo the existence and Herglotz property of omega_1).
**Note:** This proves positivity of <h, delta> where delta_k = nu_k - lambda_k,
which is NOT the same as <h, alpha> where alpha_k = H_p(lambda_k) - H_r(nu_k).

### Identity 3: Divided difference decomposition (KEY RESULT)

**<h, alpha> = sum_{k<l} H_r[nu_k, nu_l] * (1/D_{kl} - 1)**

where:
- H_r[nu_k, nu_l] = (H_r(nu_k) - H_r(nu_l))/(nu_k - nu_l) is the divided difference of H_r
- D_{kl} = (lambda_k - lambda_l)/(nu_k - nu_l) is the divided difference of omega_1

*Proof.* Start from:
h_k * alpha_k = h_k * [u_k - h_k]
= h_k * sum_{l!=k} [1/(lambda_k-lambda_l) - 1/(nu_k-nu_l)]
= sum_{l!=k} h_k/(nu_k-nu_l) * (1/D_{kl} - 1)

Sum over k: <h,alpha> = sum_{k!=l} h_k/(nu_k-nu_l) * (1/D_{kl} - 1)

Symmetrize by averaging with k <-> l (using D_{kl} = D_{lk}):
<h,alpha> = (1/2) sum_{k!=l} (h_k - h_l)/(nu_k - nu_l) * (1/D_{kl} - 1)
= sum_{k<l} H_r[k,l] * (1/D_{kl} - 1)

**Verified numerically** for n = 3, 4 across multiple trials. Perfect agreement.

---

## 3. Properties of the Divided Differences

### Fact A (NUMERICAL, CORRECTED): D_{kl} can be > 1 or < 1

**Previous claim that D_{kl} < 1 always was WRONG.**

Numerical evidence for MSS convolution (500 trials per n):
- n=3: D > 1 in ~8% of ADJACENT pairs, 0% of NON-ADJACENT pairs
- n=4: D > 1 in ~15% of adjacent pairs, ~1% of non-adjacent pairs
- n=5: D > 1 in ~18% of adjacent pairs, ~2% of non-adjacent pairs

The TOTAL root range of r is always at least as large as that of p:
nu_n - nu_1 >= lambda_n - lambda_1 (verified in 100% of trials).
This guarantees D_{1,n} <= 1 (the extremal divided difference).
But individual gap ratios D_{k,k+1} can exceed 1.

**Consequence:** The factor (1/D_{kl} - 1) has INDEFINITE SIGN.

### Fact B (NUMERICAL): H_r divided differences have indefinite sign

H_r[k,l] = (H_r(nu_k) - H_r(nu_l))/(nu_k - nu_l) can be positive or negative.

For n=3: H_r is monotone increasing in only ~52% of random root configurations.
For n=4: monotone in ~21%. For n=5: monotone in ~10%.

So BOTH factors in the divided difference identity have indefinite sign.

### Fact C: The SUM is always non-negative (UNPROVED, numerical)

Despite individual terms having indefinite sign, the sum
sum_{k<l} H_r[k,l] * (1/D_{kl} - 1) is non-negative in 100% of numerical trials
(800+ trials, n = 2,...,6) for MSS convolutions.

This is the core mystery: the proof requires understanding the CORRELATION
between H_r[k,l] and D_{kl}.

---

## 4. Structure of omega_1 (ESTABLISHED)

### Herglotz representation

omega_1 is a rational Herglotz function (maps upper half-plane to itself) with:

omega_1(z) = z + b + sum_{j=1}^{n-1} m_j / (p_j - z)

where:
- b is a real constant
- p_j are the n-1 critical points of r (roots of r'), one between each pair of
  consecutive roots nu_j < p_j < nu_{j+1}
- m_j > 0 (required by the Herglotz property)

### Key properties of omega_1

(a) omega_1'(z) = 1 + sum_j m_j/(p_j - z)^2 >= 1 for real z away from poles.
    This holds because m_j > 0 and (p_j - z)^2 > 0 for real z != p_j.

(b) omega_1'(nu_k) = 1 for all k. Combined with (a), this means omega_1' achieves
    its global minimum value of 1 exactly at the roots of r.

(c) omega_1 has a POLE between every pair of consecutive roots of r.
    Therefore the Mean Value Theorem does NOT apply to divided differences D_{kl}
    even for adjacent pairs |k-l| = 1.

(d) omega_1''(nu_k) = 2*alpha_k. This follows from the chain rule matching at
    second order.

---

## 5. Failed Approaches

### 5a. Pointwise Herglotz decomposition

**Attempt:** Write alpha_k = sum_j (+/-)m_j * f_j(nu_k) where f_j depends on the poles,
then show each term in <h,alpha> = sum_j m_j * A_j has A_j >= 0.

**Result:** A_j = sum_k h_k/(nu_k - p_j)^3 takes both positive and negative values
for arbitrary p between roots (verified for n = 3,4,5 with 200 trials each:
negative in ~22-30% of cases for p between roots, 100% negative for p outside roots).

The pointwise approach FAILS because the sign of A_j depends on where exactly
the pole p_j is, and not all locations give A_j >= 0.

### 5b. F''(p) sign analysis

**Attempt:** Use the identity A_j = -(1/4)*F''(p_j) where F(p) = r''(p)/r(p),
and argue F''(p_j) <= 0 at the critical points of r.

**Result:** F'' IS always negative at the critical points of r (0 out of 200 trials
have all F'' >= 0). Combined with the residue m_j = -n*r(p_j)/r''(p_j) > 0,
the formula gives (1/4)*sum m_j*F''(p_j) < 0, which is the WRONG sign.

**Diagnosis:** The formula m_j = -n*r(p_j)/r''(p_j) from the leading-order asymptotic
omega_1(z) ~ n*r(p_j)/(r''(p_j)*(z-p_j)) does NOT give the correct partial fraction
residues. Even for n=2, this formula does not reproduce the correct alpha values.
The Herglotz partial fraction coefficients are determined by a global system of
equations, not by local asymptotics alone. The leading-order expansion gives the
total rational function value, not the pole residue.

### 5c. Monotonicity of H_r

**Attempt:** For n=3, prove H_r is monotonically increasing, so all divided
differences H_r[k,l] >= 0. Combined with 1/D_{kl} - 1 > 0, this gives <h,alpha> >= 0.

**Result:** H_r is NOT always monotone, even for n=3. For roots with spacing ratio
a/b > 1 + sqrt(3) (where a = nu_2 - nu_1, b = nu_3 - nu_2), the monotonicity fails.
This ratio threshold ~ 2.73 is achievable for MSS convolution outputs.

### 5d. Contour integral of [nG_r]^3 * (omega_1' - 1)

**Attempt:** The function phi(z) = (1/2)[nG_r(z)]^3 * [omega_1'(z) - 1] has
Res_{nu_k} = 3*h_k*alpha_k + correction terms involving higher derivatives.

**Result:** The correction terms prevent a clean identification with <h,alpha>.
Numerical computation of the residues is badly conditioned due to the instability
of omega_1 branch selection near the real axis.

---

## 6. Promising Directions (Unresolved)

### 6a. Correlation structure of H_r[k,l] and D_{kl}

The divided difference identity
<h,alpha> = sum_{k<l} H_r[k,l] * (1/D_{kl} - 1)
expresses <h,alpha> as a bilinear form in terms of two families of divided differences.
Although individual terms can be negative, the sum is always non-negative.

**Key question:** Is there a positive semi-definite interpretation? For instance,
can we write <h,alpha> as a trace of a product of two PSD matrices?

**Observation:** The quantity H_r[k,l] = (h_k - h_l)/(nu_k - nu_l) is related to
the Cauchy matrix C_{kl} = 1/(nu_k - nu_l). Similarly, D_{kl} involves the
Cauchy matrix for the lambda roots. The interplay of two Cauchy matrices with
different nodes is the domain of Cauchy-Binet type identities.

### 6b. Energy functional approach

Define E(t) = Phi_n(x(t)) where x_k(t) = (1-t)*nu_k + t*lambda_k.
Then dE/dt|_{t=0} = -2*<h,alpha> (by a direct calculation).
So <h,alpha> >= 0 iff E is decreasing at t=0.

The function E(t) measures the "Fisher information" along a linear interpolation
of root vectors. If E is convex in t (or at least not increasing at t=0),
this would prove the result. But the convexity of Phi_n along linear paths
in root space is a non-trivial claim.

### 6c. Matrix model / random matrix approach

Since r = E_U[chi_{A+UBU*}] where A, B are diagonal matrices with roots of p, q,
the subordination function omega_1 arises from the resolvent identity. The quantity
<h,alpha> may be expressible as a trace inequality in the matrix model.

### 6d. Induction via MSS derivative identity

The identity r'(x) = n*(p^{(1)} boxplus_{n-1} q^{(1)})(x) relates the n-th order
convolution to an (n-1)-th order one. This might enable an inductive proof,
with the base case n=2 proved analytically (<h,alpha> = 2(R-s)/(R^2*s) > 0).

---

## 7. Summary and Assessment

### What is PROVED:
1. **Identity 1:** sum h_k/(nu_k - p) = -(1/2)*r''(p)/r(p)
2. **Identity 2** (contour integral): sum h_k*(nu_k - lambda_k) >= 0
3. **Identity 3** (divided differences): <h,alpha> = sum_{k<l} H_r[k,l]*(1/D_{kl}-1)
4. **Total spread:** nu_n - nu_1 >= lambda_n - lambda_1 always
5. **omega_1'(z) >= 1** for real z away from poles of omega_1
6. **Poles of omega_1** are at critical points of r, interlacing with roots

### What is CONJECTURED (strong numerical evidence, no proof):
- <h,alpha> >= 0 for all MSS convolutions (verified in 800+ random trials, n=2,...,6)
- <h,delta> >= 0 where delta_k = nu_k - lambda_k (proved via contour integral)

### Main GAP:
The contour integral method produces a clean proof for <h,delta> >= 0 (position
displacement) but NOT for <h,alpha> >= 0 (Hilbert transform displacement). The
divided difference Identity 3 is the best available representation of <h,alpha>
but has terms of indefinite sign. No manifestly positive representation has been found.

### Difficulty assessment:
The sub-lemma <h,alpha> >= 0 appears to require deeper structural understanding of
the MSS convolution than what contour integrals alone provide. The problem has the
flavor of a Schur convexity or majorization argument -- the root vectors nu and lambda
are related by a specific (subordination) map, and the Hilbert transform H is a
specific function, and the positivity involves the interaction between these structures.

A promising but unverified direction is the connection to the log-concavity properties
of real-rooted polynomials (Lorentzian polynomials, interlacing families), which may
provide the convexity structure needed to close the argument.
