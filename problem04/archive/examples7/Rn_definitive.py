"""
PROVER-10: Definitive computation of C_n and R_n.

KEY FINDINGS FROM PRIOR ANALYSIS:
1. ta_k = (-1)^k * a_k / C(n,k) are the normalized coefficients.
2. ta_2 and ta_3 are additive under MSS convolution. ta_4 is NOT.
3. The FINITE FREE CUMULANTS kappa_k (which ARE all additive) differ
   from the simple ta_k at order 4+.
4. For n=3: 1/Phi_3 = (2/9)*k2 - (2/27)*k3^2/k2^2 where k2=-e2, k3=(9/2)*e3
5. The conjecture C_n = 2/(n(n-1)) is WRONG for n=3 (gives 1/3, actual 2/9).

THIS SCRIPT:
- Correctly computes the finite free cumulants (verifying additivity)
- Expresses 1/Phi_n in cumulants for n=2,3,4,5
- Extracts C_n for each n
- Tests superadditivity numerically
"""

import numpy as np
from itertools import combinations
from math import comb, factorial
from fractions import Fraction
import sympy as sp
from sympy import (symbols, Rational, expand, simplify, factor, cancel,
                   together, Poly, Integer, Symbol)
np.random.seed(42)

# ============================================================
# CORRECT FINITE FREE CUMULANTS
# ============================================================
# Reference: Arizmendi, Perales (2018), "Cumulants for Finite Free Convolution"
#
# For a monic polynomial p(x) = x^n + c_1*x^{n-1} + ... + c_n,
# define the NORMALIZED COEFFICIENTS:
#   a_k = c_k / (-1)^k / C(n,k)   [Arizmendi-Perales notation]
# or equivalently a_k = (-1)^k * c_k / C(n,k)
#
# The finite free cumulants kappa_1, ..., kappa_n are defined by the
# MOMENT-CUMULANT FORMULA on the non-crossing partition lattice:
#
#   a_k = sum_{pi in NC(k)} prod_{V in pi} kappa_{|V|} / prod_{V in pi} n^{|V|-1} ...
#
# Actually, the CORRECT relation from Arizmendi-Perales Theorem 3.6 is:
# Define m_k = a_k (the normalized "moments"). Then the FREE cumulants
# kappa_k are defined via the R-transform or equivalently:
#
# m_k = sum_{pi in NC(k)} (n-k+|pi|)! * n^{k-2*|pi|} / (n-k)! * prod_{V} kappa_{|V|}
#
# This is complex. Let me instead COMPUTE kappa_k by REQUIRING ADDITIVITY.

# ============================================================
# STRATEGY: Find the cumulants by requiring additivity under MSS
# ============================================================

def phi_n_num(roots):
    n = len(roots)
    total = 0.0
    for i in range(n):
        H = sum(1.0/(roots[i]-roots[j]) for j in range(n) if j != i)
        total += H**2
    return total

def mss_convolve(roots_p, roots_q):
    n = len(roots_p)
    assert len(roots_q) == n
    p_coeffs = np.poly(roots_p)
    q_coeffs = np.poly(roots_q)
    r_coeffs = np.zeros(n + 1)
    r_coeffs[0] = 1.0
    for k in range(1, n + 1):
        ck = 0.0
        for i in range(k + 1):
            j = k - i
            if i <= n and j <= n:
                w = factorial(n - i) * factorial(n - j) / (factorial(n) * factorial(n - k))
                ck += w * p_coeffs[i] * q_coeffs[j]
        r_coeffs[k] = ck
    r_roots = np.roots(r_coeffs)
    return np.sort(np.real(r_roots))

def get_ta(roots):
    """Get normalized coefficients ta_k = (-1)^k * a_k / C(n,k)."""
    n = len(roots)
    coeffs = np.poly(roots)
    ta = {}
    for k in range(1, n+1):
        ta[k] = ((-1)**k * coeffs[k]) / comb(n, k)
    return ta

# ============================================================
# PART 1: Verify which quantities are additive
# ============================================================
print("=" * 70)
print("PART 1: Testing additivity of various quantities under MSS")
print("=" * 70)

# For n=4, test what IS additive.
# We know: ta_2, ta_3 are additive. ta_4 is NOT.
# The finite free cumulants require a correction at order 4:
# kappa_4 = some_function(ta_2, ta_3, ta_4) that IS additive.

# The correction must involve ta_2^2 (since ta_2 is additive,
# ta_2(p⊞q) = ta_2(p)+ta_2(q), so ta_2(p⊞q)^2 = (ta_2(p)+ta_2(q))^2
# ≠ ta_2(p)^2 + ta_2(q)^2 in general).

# So kappa_4 = alpha*ta_4 + beta*ta_2^2 with alpha, beta chosen so that
# alpha*(ta_4(p)+ta_4(q)+CROSS_TERM) + beta*(ta_2(p)+ta_2(q))^2 is additive.
# The CROSS_TERM from ta_4 needs to cancel beta*2*ta_2(p)*ta_2(q).

# Actually: MSS convolution has ta_k(p⊞q) = ta_k(p) + ta_k(q) ONLY for k<=n-1??
# Let me check more carefully.

n = 4
print(f"\nn = {n}: Testing ta_k additivity")
errors_ta = {k: [] for k in range(1, n+1)}

for trial in range(200):
    rp = np.random.randn(n) * 2
    rp = rp - np.mean(rp)
    rp = np.sort(rp)
    rq = np.random.randn(n) * 2
    rq = rq - np.mean(rq)
    rq = np.sort(rq)

    if min(np.diff(rp)) < 0.2 or min(np.diff(rq)) < 0.2:
        continue

    rr = mss_convolve(rp, rq)

    ta_p = get_ta(rp)
    ta_q = get_ta(rq)
    ta_r = get_ta(rr)

    for k in range(1, n+1):
        errors_ta[k].append(ta_r[k] - (ta_p[k] + ta_q[k]))

for k in range(1, n+1):
    errs = np.array(errors_ta[k])
    print(f"  ta_{k}: max |error| = {np.max(np.abs(errs)):.2e}, mean = {np.mean(errs):.2e}")

# Now test n=5
n = 5
print(f"\nn = {n}: Testing ta_k additivity")
errors_ta_5 = {k: [] for k in range(1, n+1)}

for trial in range(200):
    rp = np.random.randn(n) * 2
    rp = rp - np.mean(rp)
    rp = np.sort(rp)
    rq = np.random.randn(n) * 2
    rq = rq - np.mean(rq)
    rq = np.sort(rq)

    if min(np.diff(rp)) < 0.2 or min(np.diff(rq)) < 0.2:
        continue

    rr = mss_convolve(rp, rq)

    ta_p = get_ta(rp)
    ta_q = get_ta(rq)
    ta_r = get_ta(rr)

    for k in range(1, n+1):
        errors_ta_5[k].append(ta_r[k] - (ta_p[k] + ta_q[k]))

for k in range(1, n+1):
    errs = np.array(errors_ta_5[k])
    print(f"  ta_{k}: max |error| = {np.max(np.abs(errs)):.2e}, mean = {np.mean(errs):.2e}")

# ============================================================
# PART 2: Find the correct additive cumulants
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Finding correct finite free cumulants")
print("=" * 70)

# From the MSS convolution formula, the coefficient a_k of x^{n-k} in p⊞q is:
# (p⊞q)_k = sum_{i+j=k} C(n-i,j)*C(n-j,i)/C(n,k)^2 * C(n,i)*C(n,j) * p_i * q_j
# Wait, let me re-derive from the MSS formula.
#
# The MSS convolution: if p(x) = x^n + p_1*x^{n-1}+...+p_n and q similarly,
# then (p⊞q)(x) = x^n + r_1*x^{n-1}+...+r_n where:
# r_k = sum_{i+j=k} w(n,i,j) * p_i * q_j
# with w(n,i,j) = (n-i)!(n-j)! / (n! * (n-k)!)
#
# In terms of ta_k = (-1)^k p_k / C(n,k):
# p_k = (-1)^k * C(n,k) * ta_k(p)
# r_k = sum_{i+j=k} w(n,i,j) * (-1)^i C(n,i) ta_i(p) * (-1)^j C(n,j) ta_j(q)
# (-1)^k C(n,k) ta_k(r) = sum_{i+j=k} w(n,i,j) (-1)^{i+j} C(n,i) C(n,j) ta_i(p) ta_j(q)
# ta_k(r) = sum_{i+j=k} w(n,i,j) C(n,i) C(n,j) / C(n,k) * ta_i(p) * ta_j(q)

# Define W_k(i,j) = w(n,i,j) * C(n,i) * C(n,j) / C(n,k)
#                  = (n-i)!(n-j)! / (n!(n-k)!) * C(n,i)*C(n,j)/C(n,k)
# For i=0: w(n,0,k) = n!(n-k)!/(n!(n-k)!) = 1, and C(n,0)=1
# So W_k(0,k) = 1*1*C(n,k)/C(n,k) = 1
# Similarly W_k(k,0) = 1.
# For i=0 or j=0 terms: ta_k(r) = ta_k(p) + ta_k(q) + CROSS_TERMS
# So the CROSS_TERMS for k=4 are:
# W_4(1,3)*ta_1(p)*ta_3(q) + W_4(3,1)*ta_3(p)*ta_1(q) + W_4(2,2)*ta_2(p)*ta_2(q)

# For CENTERED polys (ta_1 = 0):
# ta_4(r) = ta_4(p) + ta_4(q) + W_4(2,2)*ta_2(p)*ta_2(q)
# This confirms ta_4 is NOT additive, and the correction is:
# kappa_4 := ta_4 - (W_4(2,2)/2) * ta_2^2
# so that kappa_4(r) = kappa_4(p) + kappa_4(q).

# Let's compute W_4(2,2) for various n.
print("Computing W_k(i,j) = w(n,i,j)*C(n,i)*C(n,j)/C(n,k)")
print("For centered polys, cross term for ta_4 is W_4(2,2)*ta_2(p)*ta_2(q)")

for n in [4, 5, 6]:
    # W_4(2,2) = w(n,2,2) * C(n,2)^2 / C(n,4)
    # w(n,2,2) = (n-2)!(n-2)! / (n!(n-4)!)
    w = factorial(n-2)**2 / (factorial(n) * factorial(n-4))
    W = w * comb(n,2)**2 / comb(n,4)
    print(f"  n={n}: W_4(2,2) = {w:.6f} * {comb(n,2)**2} / {comb(n,4)} = {W:.6f}")
    frac_W = Fraction(W).limit_denominator(1000)
    print(f"         = {frac_W}")

# For n=4: W_4(2,2) = (2!*2!)/(4!*0!) * C(4,2)^2/C(4,4)
# = 4/(24*1) * 36/1 = (4*36)/24 = 144/24 = 6
# So kappa_4 = ta_4 - 3*ta_2^2 (divide by 2 since there are 2 cross terms ta_2(p)*ta_2(q))
# Wait: ta_4(r) = ta_4(p) + ta_4(q) + 6*ta_2(p)*ta_2(q)
# kappa_4 = ta_4 - alpha*ta_2^2
# kappa_4(r) = ta_4(p) + ta_4(q) + 6*ta_2(p)*ta_2(q) - alpha*(ta_2(p)+ta_2(q))^2
#            = (ta_4(p)-alpha*ta_2(p)^2) + (ta_4(q)-alpha*ta_2(q)^2) + (6-2*alpha)*ta_2(p)*ta_2(q)
# For additivity: 6 - 2*alpha = 0, so alpha = 3.
# kappa_4 = ta_4 - 3*ta_2^2

print("\nFor n=4 centered:")
print("  kappa_4 = ta_4 - 3*ta_2^2")

# Verify numerically
print("\n  Verifying kappa_4 = ta_4 - 3*ta_2^2 additivity for n=4:")
n = 4
errors_k4 = []
for trial in range(100):
    rp = np.random.randn(n) * 2; rp -= np.mean(rp); rp = np.sort(rp)
    rq = np.random.randn(n) * 2; rq -= np.mean(rq); rq = np.sort(rq)
    if min(np.diff(rp)) < 0.2 or min(np.diff(rq)) < 0.2:
        continue
    rr = mss_convolve(rp, rq)
    ta_p = get_ta(rp); ta_q = get_ta(rq); ta_r = get_ta(rr)
    k4_p = ta_p[4] - 3*ta_p[2]**2
    k4_q = ta_q[4] - 3*ta_q[2]**2
    k4_r = ta_r[4] - 3*ta_r[2]**2
    errors_k4.append(k4_r - (k4_p + k4_q))

print(f"  Max |error|: {max(abs(e) for e in errors_k4):.2e}")

# For n=5: compute W_4(2,2), and also need W_5 cross terms
print("\nFor n=5 centered:")
# ta_4(r) = ta_4(p) + ta_4(q) + W_4(2,2)*ta_2(p)*ta_2(q)
# For n=5: W_4(2,2) = w(5,2,2)*C(5,2)^2/C(5,4) = (3!*3!)/(5!*1!)*100/5
# = 36/120 * 20 = 36*20/120 = 720/120 = 6
# Hmm, same as n=4? Let me compute more carefully.
n = 5
w_422 = factorial(n-2)**2 / (factorial(n) * factorial(n-4))
W_422 = w_422 * comb(n,2)**2 / comb(n,4)
print(f"  W_4(2,2) for n=5: {W_422}")
# So kappa_4 = ta_4 - 3*ta_2^2 for n=5 too? Let me verify.

print("  Verifying kappa_4 = ta_4 - 3*ta_2^2 additivity for n=5:")
errors_k4_5 = []
for trial in range(100):
    rp = np.random.randn(n) * 2; rp -= np.mean(rp); rp = np.sort(rp)
    rq = np.random.randn(n) * 2; rq -= np.mean(rq); rq = np.sort(rq)
    if min(np.diff(rp)) < 0.2 or min(np.diff(rq)) < 0.2:
        continue
    rr = mss_convolve(rp, rq)
    ta_p = get_ta(rp); ta_q = get_ta(rq); ta_r = get_ta(rr)
    k4_p = ta_p[4] - 3*ta_p[2]**2
    k4_q = ta_q[4] - 3*ta_q[2]**2
    k4_r = ta_r[4] - 3*ta_r[2]**2
    errors_k4_5.append(k4_r - (k4_p + k4_q))

print(f"  Max |error|: {max(abs(e) for e in errors_k4_5):.2e}")

# Now for kappa_5:
# ta_5(r) = ta_5(p) + ta_5(q) + W_5(2,3)*ta_2(p)*ta_3(q) + W_5(3,2)*ta_3(p)*ta_2(q)
# (no ta_1*ta_4 terms since ta_1=0 for centered)
# W_5(2,3) = w(5,2,3)*C(5,2)*C(5,3)/C(5,5)
# w(5,2,3) = (5-2)!(5-3)!/(5!(5-5)!) = 3!*2!/(120*1) = 12/120 = 1/10
# W_5(2,3) = (1/10)*10*10/1 = 10
# Similarly W_5(3,2) = w(5,3,2)*C(5,3)*C(5,2)/C(5,5) = same = 10

n = 5
w_523 = factorial(n-2)*factorial(n-3)/(factorial(n)*factorial(n-5))
W_523 = w_523 * comb(n,2)*comb(n,3)/comb(n,5)
print(f"\n  W_5(2,3) for n=5: {W_523}")
# So ta_5(r) = ta_5(p) + ta_5(q) + 10*ta_2(p)*ta_3(q) + 10*ta_3(p)*ta_2(q)
# kappa_5 = ta_5 - alpha*ta_2*ta_3
# kappa_5(r) = ta_5(p) + ta_5(q) + 10*(ta_2(p)*ta_3(q) + ta_3(p)*ta_2(q))
#              - alpha*(ta_2(p)+ta_2(q))*(ta_3(p)+ta_3(q))
#            = (ta_5(p)-alpha*ta_2(p)*ta_3(p)) + (ta_5(q)-alpha*ta_2(q)*ta_3(q))
#              + (10-alpha)*(ta_2(p)*ta_3(q) + ta_3(p)*ta_2(q))
# For additivity: alpha = 10
# kappa_5 = ta_5 - 10*ta_2*ta_3

print("  kappa_5 = ta_5 - 10*ta_2*ta_3")
print("  Verifying additivity:")
errors_k5_5 = []
for trial in range(100):
    rp = np.random.randn(n) * 2; rp -= np.mean(rp); rp = np.sort(rp)
    rq = np.random.randn(n) * 2; rq -= np.mean(rq); rq = np.sort(rq)
    if min(np.diff(rp)) < 0.2 or min(np.diff(rq)) < 0.2:
        continue
    rr = mss_convolve(rp, rq)
    ta_p = get_ta(rp); ta_q = get_ta(rq); ta_r = get_ta(rr)
    k5_p = ta_p[5] - 10*ta_p[2]*ta_p[3]
    k5_q = ta_q[5] - 10*ta_q[2]*ta_q[3]
    k5_r = ta_r[5] - 10*ta_r[2]*ta_r[3]
    errors_k5_5.append(k5_r - (k5_p + k5_q))

print(f"  Max |error|: {max(abs(e) for e in errors_k5_5):.2e}")

# For n=6: need kappa_6 as well
# ta_6(r) = ta_6(p) + ta_6(q) + cross terms involving:
# W_6(2,4), W_6(4,2), W_6(3,3) (since ta_1=0 for centered)
# Also W_6(1,5) etc. but ta_1=0.

n = 6
print(f"\nFor n={n} centered:")
# Compute all W_6(i,j) with i+j=6, i>=2, j>=2
for i_val in range(2, 5):
    j_val = 6 - i_val
    if j_val < 2 or j_val > n:
        continue
    w_ij = factorial(n-i_val)*factorial(n-j_val)/(factorial(n)*factorial(n-6))
    W_ij = w_ij * comb(n,i_val)*comb(n,j_val)/comb(n,6)
    print(f"  W_6({i_val},{j_val}) = {W_ij}")

# For kappa_6 = ta_6 - W_6(2,4)*ta_2*ta_4 - W_6(3,3)/2*ta_3^2 - ... hmm, this gets complex.
# Also need to subtract kappa_4*ta_2 type terms.
# Actually the formula is more nuanced for k=6 since we need to consider
# all NC(6) partitions.

# Let me just compute the general cumulant-from-ta formula.
# The MSS convolution formula gives:
# ta_k(p⊞q) = sum_{i+j=k} W_k(i,j) * ta_i(p) * ta_j(q)
# where W_k(0,k) = W_k(k,0) = 1.
# This is a POLYNOMIAL algebra: the ta_k multiply via a convolution.
# The cumulants are the PRIMITIVE elements of this algebra
# (additive under the comultiplication ta_k -> sum W_k(i,j) ta_i ⊗ ta_j).

# For centered polys, the cumulants are:
# kappa_2 = ta_2
# kappa_3 = ta_3
# kappa_4 = ta_4 - (W_4(2,2)/2)*ta_2^2 = ta_4 - 3*ta_2^2
# kappa_5 = ta_5 - W_5(2,3)*ta_2*ta_3
# kappa_6 = ta_6 - W_6(2,4)*ta_2*ta_4 - (W_6(3,3)/2)*ta_3^2
#           + correction for "double" non-primitivity

# Actually, this is exactly the LOG of the "exponential" generating function!
# In the algebra of sequences (ta_k) with MSS multiplication, the cumulants
# are the logarithm of the "character".

# Let me compute W_k(i,j) systematically.
def compute_W(n_val, k, i_val, j_val):
    """W_k(i,j) for MSS convolution of degree-n polynomials."""
    if i_val + j_val != k:
        return 0.0
    if i_val < 0 or j_val < 0 or i_val > n_val or j_val > n_val or k > n_val:
        return 0.0
    w = factorial(n_val-i_val)*factorial(n_val-j_val)/(factorial(n_val)*factorial(n_val-k))
    W = w * comb(n_val, i_val) * comb(n_val, j_val) / comb(n_val, k)
    return W

# For centered case (ta_1=0), the cross terms in ta_k(p⊞q) are:
# sum_{i=2}^{k-2} W_k(i,k-i) * ta_i(p)*ta_{k-i}(q)
# (since ta_0=1 gives W_k(0,k)=1 and W_k(k,0)=1, which are the diagonal terms)
# (and ta_1=0 eliminates the i=1 and j=1 terms)

# The cumulants (primitive elements) are obtained by the LOG formula:
# kappa_k = ta_k - sum over non-trivial "forest" decompositions...
# This is the classical logarithm in the connected graded Hopf algebra.

# For our purposes, the cumulant formula for centered case is:
# kappa_2 = ta_2
# kappa_3 = ta_3
# kappa_4 = ta_4 - sum_{i+j=4, i,j>=2} (W_4(i,j)/2) * kappa_i * kappa_j
#         = ta_4 - (W_4(2,2)/2)*kappa_2^2
# kappa_5 = ta_5 - sum_{i+j=5, i,j>=2} W_5(i,j) * kappa_i * kappa_j / 2
#         Wait, W_5(2,3) ≠ W_5(3,2) in general? Let me check.

n = 5
print(f"\nW_5(2,3) = {compute_W(5, 5, 2, 3)}")
print(f"W_5(3,2) = {compute_W(5, 5, 3, 2)}")
# These should be equal for the centered case.
# If they are, then kappa_5 = ta_5 - W_5(2,3)*kappa_2*kappa_3

# For kappa_6, we need:
# ta_6(p⊞q) = ta_6(p)+ta_6(q) + W_6(2,4)*[ta_2(p)*ta_4(q)+ta_4(p)*ta_2(q)]
#              + W_6(3,3)*ta_3(p)*ta_3(q)
# But ta_4 = kappa_4 + 3*kappa_2^2 (since kappa_4 = ta_4 - 3*ta_2^2 = ta_4 - 3*kappa_2^2)
# So the cross terms in terms of kappa are more complex.

# SIMPLER APPROACH: Use the recursive logarithm formula directly.
# Define c_k = sum_{i+j=k, i,j>=1} W_k(i,j)*ta_i*ta_j (cross terms only, excluding i=0,j=0)
# For centered: c_k = sum_{i+j=k, i,j>=2} W_k(i,j)*ta_i*ta_j (since ta_1=0)
# Then the cumulants are defined by:
# ta_k = kappa_k + sum over ordered forests ...
# This is getting very complex. Let me just numerically determine the cumulants.

print("\n" + "=" * 70)
print("NUMERICAL DETERMINATION OF CUMULANTS")
print("=" * 70)

# For n=6 centered, I need kappa_4, kappa_5, kappa_6.
# Strategy: use the CROSS TERM structure.
# ta_k(p⊞q) = ta_k(p) + ta_k(q) + C_k(ta(p), ta(q))
# where C_k is the cross term (sum over i+j=k, i,j>=2).
# The cumulant kappa_k must satisfy:
# kappa_k(p⊞q) = kappa_k(p) + kappa_k(q)
# kappa_k = ta_k + polynomial_in_lower_ta

# Compute all W_k(i,j) for n=6
n = 6
print(f"\nMSS cross-term coefficients W_k(i,j) for n={n}:")
for k in range(2, n+1):
    for i_val in range(2, k-1):
        j_val = k - i_val
        if j_val >= 2:
            W = compute_W(n, k, i_val, j_val)
            if abs(W) > 1e-10:
                frac = Fraction(W).limit_denominator(10000)
                print(f"  W_{k}({i_val},{j_val}) = {frac}")

# Now determine kappa_4 for general n (centered):
# kappa_4 = ta_4 - (W_4(2,2)/2) * ta_2^2
# For the cross term: W_4(2,2)/2 (divided by 2 because W_4(2,2) appears for
# both i=2,j=2 orderings, but the bilinear form has symmetry)
# Actually: in ta_4(p⊞q), the (2,2) cross term gives W_4(2,2)*ta_2(p)*ta_2(q).
# For kappa_4 = ta_4 - alpha*ta_2^2:
# kappa_4(p⊞q) = ta_4(p⊞q) - alpha*(ta_2(p)+ta_2(q))^2
#               = ta_4(p)+ta_4(q)+W_4(2,2)*ta_2(p)*ta_2(q) - alpha*ta_2(p)^2 - alpha*ta_2(q)^2 - 2*alpha*ta_2(p)*ta_2(q)
#               = kappa_4(p)+kappa_4(q) + (W_4(2,2)-2*alpha)*ta_2(p)*ta_2(q)
# So alpha = W_4(2,2)/2.

print(f"\n" + "=" * 70)
print("GENERAL CUMULANT FORMULAS (centered)")
print("=" * 70)

for n in [3, 4, 5, 6]:
    W422 = compute_W(n, 4, 2, 2)
    alpha_4 = W422 / 2
    print(f"\nn={n}: kappa_4 = ta_4 - {Fraction(alpha_4).limit_denominator(1000)} * ta_2^2")

    if n >= 5:
        W523 = compute_W(n, 5, 2, 3)
        W532 = compute_W(n, 5, 3, 2)
        # ta_5(p⊞q) has cross terms W_5(2,3)*ta_2(p)*ta_3(q) + W_5(3,2)*ta_3(p)*ta_2(q)
        # kappa_5 = ta_5 - beta*ta_2*ta_3
        # Cross: W_5(2,3)*ta_2(p)*ta_3(q) + W_5(3,2)*ta_3(p)*ta_2(q) - beta*(ta_2(p)*ta_3(q)+ta_3(p)*ta_2(q)+ta_2(p)*ta_3(p)+... )
        # Wait: kappa_5 = ta_5 - beta*ta_2*ta_3
        # kappa_5(p⊞q) = ta_5(r) - beta*ta_2(r)*ta_3(r)
        #               = [ta_5(p)+ta_5(q)+W52_cross] - beta*(ta_2(p)+ta_2(q))*(ta_3(p)+ta_3(q))
        # where W52_cross = W_5(2,3)*ta_2(p)*ta_3(q) + W_5(3,2)*ta_3(p)*ta_2(q)
        # = ta_5(p)-beta*ta_2(p)*ta_3(p) + ta_5(q)-beta*ta_2(q)*ta_3(q)
        #   + W_5(2,3)*ta_2(p)*ta_3(q) + W_5(3,2)*ta_3(p)*ta_2(q)
        #   - beta*ta_2(p)*ta_3(q) - beta*ta_3(p)*ta_2(q)
        # = kappa_5(p) + kappa_5(q) + (W_5(2,3)-beta)*ta_2(p)*ta_3(q) + (W_5(3,2)-beta)*ta_3(p)*ta_2(q)
        # For additivity: W_5(2,3) = beta AND W_5(3,2) = beta
        # So beta = W_5(2,3) (need W_5(2,3) = W_5(3,2))
        assert abs(W523 - W532) < 1e-10, f"W_5(2,3)={W523} ≠ W_5(3,2)={W532}"
        print(f"  kappa_5 = ta_5 - {Fraction(W523).limit_denominator(1000)} * ta_2 * ta_3")

    if n >= 6:
        W624 = compute_W(n, 6, 2, 4)
        W642 = compute_W(n, 6, 4, 2)
        W633 = compute_W(n, 6, 3, 3)
        # ta_6(p⊞q) cross terms (only i,j>=2):
        # W_6(2,4)*ta_2(p)*ta_4(q) + W_6(4,2)*ta_4(p)*ta_2(q) + W_6(3,3)*ta_3(p)*ta_3(q)
        # But ta_4 = kappa_4 + alpha_4*ta_2^2, so:
        # ta_4(p)*ta_2(q) = (kappa_4(p)+alpha_4*ta_2(p)^2)*ta_2(q)
        # = kappa_4(p)*ta_2(q) + alpha_4*ta_2(p)^2*ta_2(q)
        #
        # This mixes kappa_4 and ta_2, making the formula for kappa_6 more complex.
        # The general pattern requires the full cumulant recursion.
        #
        # Let me just compute kappa_6 numerically by fitting.

        print(f"  W_6(2,4)={Fraction(W624).limit_denominator(1000)}, W_6(4,2)={Fraction(W642).limit_denominator(1000)}, W_6(3,3)={Fraction(W633).limit_denominator(1000)}")

        # For kappa_6, the complete formula involves:
        # kappa_6 = ta_6 - gamma_24*ta_2*ta_4 - gamma_33*ta_3^2 + delta_222*ta_2^3
        # where gamma_24 = W_6(2,4), gamma_33 = W_6(3,3)/2
        # and delta_222 accounts for the triple correction.
        # The triple correction: since ta_4 = kappa_4 + alpha*ta_2^2,
        # the cross term W_6(2,4)*ta_2*ta_4 = W_6(2,4)*ta_2*(kappa_4+alpha*ta_2^2)
        # introduces a ta_2^3 term that needs to be cancelled.
        #
        # Full recursion:
        # kappa_6 = ta_6 - W_6(2,4)*ta_2*kappa_4 - (W_6(3,3)/2)*kappa_3^2
        #           - W_6(2,4)*alpha_4*ta_2^3 ... no, this isn't right.
        #
        # Actually: the correct formula comes from the LOGARITHM of the convolution algebra.
        # Let me just compute it numerically.

# ============================================================
# NUMERICAL FITTING FOR kappa_6
# ============================================================
print("\n" + "=" * 70)
print("FITTING kappa_6 for n=6")
print("=" * 70)

n = 6
# We need kappa_6 = ta_6 + a*ta_2*ta_4 + b*ta_3^2 + c*ta_2^3
# to be additive. (Higher-order terms in ta involve products of 3+ ta's.)
# Collect data from MSS convolutions and fit.

ta_data_p = []
ta_data_q = []
ta_data_r = []

for trial in range(500):
    rp = np.random.randn(n) * 1.5; rp -= np.mean(rp); rp = np.sort(rp)
    rq = np.random.randn(n) * 1.5; rq -= np.mean(rq); rq = np.sort(rq)
    if min(np.diff(rp)) < 0.15 or min(np.diff(rq)) < 0.15:
        continue
    rr = mss_convolve(rp, rq)
    ta_data_p.append(get_ta(rp))
    ta_data_q.append(get_ta(rq))
    ta_data_r.append(get_ta(rr))

print(f"  Generated {len(ta_data_p)} valid samples")

# For kappa_6 = ta_6 + a*ta_2*ta_4 + b*ta_3^2 + c*ta_2^3:
# Additivity: kappa_6(r) = kappa_6(p) + kappa_6(q)
# i.e., ta_6(r) + a*ta_2(r)*ta_4(r) + b*ta_3(r)^2 + c*ta_2(r)^3
#     = [ta_6(p) + a*ta_2(p)*ta_4(p) + b*ta_3(p)^2 + c*ta_2(p)^3]
#     + [ta_6(q) + a*ta_2(q)*ta_4(q) + b*ta_3(q)^2 + c*ta_2(q)^3]
#
# This means: the CROSS TERMS must cancel:
# ta_6(r) - ta_6(p) - ta_6(q) = -(a,b,c)*(cross terms of ta_2*ta_4, ta_3^2, ta_2^3)
#
# Let me define: delta_6 = ta_6(r) - ta_6(p) - ta_6(q)
# cross_24 = ta_2(r)*ta_4(r) - ta_2(p)*ta_4(p) - ta_2(q)*ta_4(q)
# cross_33 = ta_3(r)^2 - ta_3(p)^2 - ta_3(q)^2
# cross_222 = ta_2(r)^3 - ta_2(p)^3 - ta_2(q)^3
# Then: delta_6 + a*cross_24 + b*cross_33 + c*cross_222 = 0

delta_6_list = []
cross_24_list = []
cross_33_list = []
cross_222_list = []

for i in range(len(ta_data_p)):
    tp, tq, tr = ta_data_p[i], ta_data_q[i], ta_data_r[i]
    delta_6_list.append(tr[6] - tp[6] - tq[6])
    cross_24_list.append(tr[2]*tr[4] - tp[2]*tp[4] - tq[2]*tq[4])
    cross_33_list.append(tr[3]**2 - tp[3]**2 - tq[3]**2)
    cross_222_list.append(tr[2]**3 - tp[2]**3 - tq[2]**3)

# Solve: delta_6 = -a*cross_24 - b*cross_33 - c*cross_222
X_k6 = np.column_stack([-np.array(cross_24_list), -np.array(cross_33_list), -np.array(cross_222_list)])
y_k6 = np.array(delta_6_list)
abc, residual_k6, _, _ = np.linalg.lstsq(X_k6, y_k6, rcond=None)
a6, b6, c6 = abc
print(f"  kappa_6 = ta_6 + ({a6:.6f})*ta_2*ta_4 + ({b6:.6f})*ta_3^2 + ({c6:.6f})*ta_2^3")
print(f"  As fractions: a={Fraction(a6).limit_denominator(1000)}, b={Fraction(b6).limit_denominator(1000)}, c={Fraction(c6).limit_denominator(1000)}")

# Verify
errors_k6 = []
for i in range(len(ta_data_p)):
    tp, tq, tr = ta_data_p[i], ta_data_q[i], ta_data_r[i]
    k6_p = tp[6] + a6*tp[2]*tp[4] + b6*tp[3]**2 + c6*tp[2]**3
    k6_q = tq[6] + a6*tq[2]*tq[4] + b6*tq[3]**2 + c6*tq[2]**3
    k6_r = tr[6] + a6*tr[2]*tr[4] + b6*tr[3]**2 + c6*tr[2]**3
    errors_k6.append(k6_r - k6_p - k6_q)

print(f"  Max |additivity error|: {max(abs(e) for e in errors_k6):.2e}")

# ============================================================
# PART 3: Express 1/Phi_n in terms of CORRECT additive cumulants
# ============================================================
print("\n" + "=" * 70)
print("PART 3: 1/Phi_n in additive cumulants")
print("=" * 70)

# Summary of cumulant formulas (centered, for any n):
# kappa_2 = ta_2
# kappa_3 = ta_3
# kappa_4 = ta_4 - 3*ta_2^2  (with W_4(2,2)/2 = 3 for n>=4)
# kappa_5 = ta_5 - W_5(2,3)*ta_2*ta_3  (W_5(2,3) depends on n)
# kappa_6 = ta_6 + a6*ta_2*ta_4 + b6*ta_3^2 + c6*ta_2^3

# Wait, I need to check: does W_4(2,2) depend on n?
print("Checking if W_k(i,j) depends on n:")
for n_val in [4, 5, 6, 7, 8, 10]:
    W = compute_W(n_val, 4, 2, 2)
    print(f"  n={n_val}: W_4(2,2) = {Fraction(W).limit_denominator(10000)}")

for n_val in [5, 6, 7, 8, 10]:
    W = compute_W(n_val, 5, 2, 3)
    print(f"  n={n_val}: W_5(2,3) = {Fraction(W).limit_denominator(10000)}")

# Oh! W_4(2,2) might depend on n!
# Let me compute: W_4(2,2) = w(n,2,2)*C(n,2)^2/C(n,4)
# w(n,2,2) = (n-2)!*(n-2)!/(n!*(n-4)!)
# = [(n-2)!]^2 / [n!*(n-4)!]
# = [(n-2)(n-3)(n-4)!]^2 / [n(n-1)(n-2)!*(n-4)!]
# Hmm let me just compute directly.

print("\nDetailed W_4(2,2) computation:")
for n_val in [4, 5, 6, 7, 8, 10, 20, 100]:
    w = factorial(n_val-2)**2 / (factorial(n_val) * factorial(n_val-4))
    C_n_2 = comb(n_val, 2)
    C_n_4 = comb(n_val, 4)
    W = w * C_n_2**2 / C_n_4
    print(f"  n={n_val}: w={w:.6f}, C(n,2)^2={C_n_2**2}, C(n,4)={C_n_4}, W_4(2,2)={W:.6f} = {Fraction(W).limit_denominator(100000)}")

# ============================================================
# KEY REALIZATION: W_4(2,2) = 6/(n-1) for degree n polynomials!
# Let me verify this.
# ============================================================
print("\nIs W_4(2,2) = 6/(n-1)?")
for n_val in [4, 5, 6, 7, 8, 10, 20]:
    W = compute_W(n_val, 4, 2, 2)
    conj = 6.0/(n_val - 1)
    print(f"  n={n_val}: W_4(2,2)={W:.8f}, 6/(n-1)={conj:.8f}, match={abs(W-conj)<1e-10}")

print("\nIs W_5(2,3) = some formula?")
for n_val in [5, 6, 7, 8, 10, 20]:
    W = compute_W(n_val, 5, 2, 3)
    # Try: 10/(n-1), 10/(n-1)^2, 12/(n-1), etc.
    for formula_name, formula_val in [
        ("10/(n-1)", 10.0/(n_val-1)),
        ("20/((n-1)(n-2))", 20.0/((n_val-1)*(n_val-2))),
        ("10*(n-2)/((n-1)*(n-3))", 10.0*(n_val-2)/((n_val-1)*(n_val-3)) if n_val > 3 else 0),
        ("30/((n-1)*(n-2))", 30.0/((n_val-1)*(n_val-2))),
    ]:
        if abs(W - formula_val) < 1e-8:
            print(f"  n={n_val}: W_5(2,3)={W:.8f} = {formula_name}")
            break
    else:
        print(f"  n={n_val}: W_5(2,3)={W:.8f}, no simple formula found")
        print(f"    1/W = {1/W:.6f}")
        print(f"    (n-1)*(n-2)/W = {(n_val-1)*(n_val-2)/W:.6f}")

# ============================================================
# IMPORTANT: The cumulants DEPEND ON n!
# This means the decomposition 1/Phi_n = C_n*kappa_2 + R_n
# has kappa_2 = ta_2 (which is n-independent for centered polys)
# but kappa_4, kappa_5, etc. are n-dependent.
# ============================================================

# Since ta_2, ta_3 are the first two additive cumulants (same for all n>=3),
# and 1/Phi_n is a FIXED function of ta_2,...,ta_n for each n,
# the decomposition makes sense in terms of ta_k.

# Let me now express 1/Phi_n for each n as a rational function of ta_k
# and extract the "C_n" coefficient of ta_2.

# But we saw that at ta_3=...=ta_n=0, Phi becomes undefined (repeated roots).
# So "C_n" is not the value at ta_3=...=0.
# Instead, C_n should be extracted from a SERIES EXPANSION or as the
# dominant term in some regime.

# For n=3: 1/Phi_3 = (-2/3)*ta_2 - (1/6)*ta_3^2/ta_2^2
# Here "C_3" = -2/3 is the part linear in ta_2 (with no ta_3).
# R_3 = -(1/6)*ta_3^2/ta_2^2 depends on BOTH ta_2 and ta_3.

# For n=4: 1/Phi_4 = (-81*t2^4*t4 + 54*t2^3*t3^2 + 18*t2^2*t4^2 - 54*t2*t3^2*t4 + 27*t3^4 - t4^3) /
#                     (9*(3*t2^2+t4)*(9*t2^3-t2*t4+3*t3^2))
# This doesn't have a simple "linear in t2" term since the whole thing is a ratio.

# Let me try a DIFFERENT decomposition:
# Write 1/Phi_n = A_n(ta) / B_n(ta) and compute A_n/B_n.
# Then C_n is the leading coefficient when ta_3,...,ta_n are "infinitesimally small"
# compared to ta_2.

# For n=4, at t4=0:
# 1/Phi_4 = t3^2*(2*t2^3+t3^2) / (3*t2^2*(3*t2^3+t3^2))
# = (2*t2^3*t3^2 + t3^4) / (9*t2^5 + 3*t2^2*t3^2)
# As t3 -> 0: both num and den -> 0 (order t3^2)
# Divide: (2*t2^3 + t3^2) / (9*t2^5/t3^2 + 3*t2^2) -> hmm
# Actually: at leading order in t3:
# Num ~ 2*t2^3*t3^2, Den ~ 9*t2^5
# So 1/Phi_4 ~ (2/9)*t3^2/t2^2 as t3 -> 0 (with t4=0)
# This has NO linear-in-t2 term!

# This suggests that for n=4, 1/Phi_4 does NOT have a "C_4 * ta_2" term
# in the way n=3 does.

# Wait, for n=3 at t3=0: 1/Phi_3 = (-2/3)*ta_2 is well-defined.
# For n=4 at t3=t4=0: 1/Phi_4 = 0/0 (undefined).
# This is because t3=t4=0 gives REPEATED ROOTS for n=4 but not for n=3.
# For n=3 at t3=0: x^3 + 3*t2*x = x*(x^2+3*t2). With t2<0, roots are 0, ±sqrt(-3*t2). Distinct! ✓
# For n=4 at t3=t4=0: x^4 + 6*t2*x^2 = x^2*(x^2+6*t2). Root 0 is REPEATED. ✗

# So the decomposition 1/Phi_n = C_n*kappa_2 + R_n ONLY makes sense when
# the polynomial space is restricted to SIMPLE ROOTS.
# The "C_n" coefficient needs to be defined more carefully.

print("\n" + "=" * 70)
print("RE-EXAMINING THE DECOMPOSITION")
print("=" * 70)

# The original statement: for monic real-rooted polynomials with SIMPLE roots,
# 1/Phi_n(p ⊞ q) >= 1/Phi_n(p) + 1/Phi_n(q)
#
# The decomposition 1/Phi_n = C_n*k2 + R_n is an identity as a function of
# the cumulants, valid wherever Phi_n is defined (i.e., simple roots).
#
# For n=3: 1/Phi_3 = (2/9)*k2 + R_3(k2,k3) where R_3 = -(2/27)*k3^2/k2^2
# C_3*k2 is the "linear part" and R_3 is the "nonlinear correction".
# Since k2 is additive and C_3>0, the C_3*k2 part is automatically superadditive.
# The inequality then requires R_3 to be "sufficiently superadditive".
#
# For n=4, we need to similarly split 1/Phi_4 into "additive part" + "correction".
# But the function is a rational function of (t2, t3, t4).
#
# Let me think about this differently.
# Write f(t2,t3,t4) = 1/Phi_4.
# We want: f(t2+s2, t3+s3, t4+s4) >= f(t2,t3,t4) + f(s2,s3,s4)
# Note: under MSS, t2 and t3 are additive (t2(r)=t2(p)+t2(q), t3(r)=t3(p)+t3(q)).
# But t4 is NOT additive: t4(r) = t4(p)+t4(q)+W_4(2,2)*t2(p)*t2(q).
# And for n>=5, t5 also has a cross term.
#
# So the correct formulation is:
# f(t2(p)+t2(q), t3(p)+t3(q), t4(p)+t4(q)+W*t2(p)*t2(q))
# >= f(t2(p),t3(p),t4(p)) + f(t2(q),t3(q),t4(q))
#
# This is MORE COMPLEX than simple superadditivity in additive variables!

# ============================================================
# APPROACH: Work with the CORRECT additive cumulants
# ============================================================
# For n=4, the additive cumulants are:
# kappa_2 = t2, kappa_3 = t3, kappa_4 = t4 - (W_4(2,2)/2)*t2^2

# Express 1/Phi_4 in terms of (kappa_2, kappa_3, kappa_4):
# t2 = kappa_2, t3 = kappa_3, t4 = kappa_4 + 3*kappa_2^2 (for n=4 where W/2=3)

t2, t3, t4 = symbols('t2 t3 t4')
K2, K3, K4 = symbols('K2 K3 K4')

# For n=4:
# 1/Phi_4 = (-81*t2^4*t4 + 54*t2^3*t3^2 + 18*t2^2*t4^2 - 54*t2*t3^2*t4 + 27*t3^4 - t4^3) /
#           (9*(3*t2^2+t4)*(9*t2^3-t2*t4+3*t3^2))

num_4 = -81*t2**4*t4 + 54*t2**3*t3**2 + 18*t2**2*t4**2 - 54*t2*t3**2*t4 + 27*t3**4 - t4**3
den_4 = 9*(3*t2**2+t4)*(9*t2**3-t2*t4+3*t3**2)

# Substitute t2=K2, t3=K3, t4=K4+3*K2^2
num_4_K = num_4.subs([(t2, K2), (t3, K3), (t4, K4 + 3*K2**2)])
den_4_K = den_4.subs([(t2, K2), (t3, K3), (t4, K4 + 3*K2**2)])

num_4_K = expand(num_4_K)
den_4_K = expand(den_4_K)

inv_phi4_K = cancel(num_4_K / den_4_K)
num_4_K_s, den_4_K_s = sp.fraction(inv_phi4_K)
num_4_K_s = expand(num_4_K_s)
den_4_K_s = expand(den_4_K_s)

print(f"\n1/Phi_4 in additive cumulants (K2,K3,K4):")
print(f"  Num: {num_4_K_s}")
print(f"  Den: {den_4_K_s}")
print(f"  Num factored: {factor(num_4_K_s)}")
print(f"  Den factored: {factor(den_4_K_s)}")

# Now check: at K3=K4=0, what is 1/Phi_4?
inv_phi4_K0 = inv_phi4_K.subs([(K3, 0), (K4, 0)])
inv_phi4_K0 = cancel(inv_phi4_K0)
print(f"\n  At K3=K4=0: 1/Phi_4 = {inv_phi4_K0}")
# If this is well-defined, then C_4 = coefficient of K2!

# ============================================================
# GRAND COMPUTATION: C_n for n=2,3,4,5,6
# ============================================================
print("\n" + "=" * 70)
print("GRAND COMPUTATION: C_n for n=2,3,4,5,6")
print("=" * 70)

# For n=2: 1/Phi_2 = -2*t2, C_2 = -2 (in ta-variables)
print("\nn=2: 1/Phi_2 = -2*ta_2")
print("  C_2(ta) = -2")

# For n=3: 1/Phi_3 = (-2/3)*t2 + R_3
print("\nn=3: 1/Phi_3 = (-2/3)*ta_2 - (1/6)*ta_3^2/ta_2^2")
print("  C_3(ta) = -2/3")

# For n=4: compute 1/Phi_4(K2, 0, 0) where K2 = ta_2, K3=K4=0 means ta_3=0, ta_4=3*ta_2^2
# So we're looking at polys with ta_3=0 and ta_4=3*ta_2^2.
# These are x^4 + 6*t2*x^2 + 3*t2^2 = x^4 + 6t2*x^2 + 3t2^2
# Hmm, does this have simple roots?
# Discriminant: at t3=0, t4=3*t2^2:
# disc = 256*(3t2^2)^3 - 128*(6t2)^2*(3t2^2)^2 + 0 - 0 + 16*(6t2)^4*(3t2^2) - 0
# = 256*27*t2^6 - 128*36*9*t2^6 + 16*1296*3*t2^6
# = 6912*t2^6 - 41472*t2^6 + 62208*t2^6 = 27648*t2^6
# This is > 0 for t2 ≠ 0! So these are valid polynomials with simple roots!

print(f"\nn=4: At K3=K4=0 (i.e., ta_3=0, ta_4=3*ta_2^2):")
print(f"  1/Phi_4 = {inv_phi4_K0}")
# This should give C_4 * K2 = C_4 * ta_2

# Verify numerically
# x^4 + 6t2*x^2 + 3t2^2 = 0 => x^2 = (-6t2 ± sqrt(36t2^2 - 12t2^2))/2 = -3t2 ± sqrt(6)*t2
# For t2 < 0 (which is the valid case):
# x^2 = -3*t2 ± sqrt(6)*|t2|
# Both positive since -3*t2 > 0 and sqrt(6)*|t2| < 3*|t2|.
# So roots are ±sqrt(-3*t2 + sqrt(6)*|t2|), ±sqrt(-3*t2 - sqrt(6)*|t2|)
t2_test = -1.0
roots_test = []
disc_quad = 36*t2_test**2 - 12*t2_test**2
x2_plus = (-6*t2_test + np.sqrt(disc_quad))/2
x2_minus = (-6*t2_test - np.sqrt(disc_quad))/2
roots_test = np.array([np.sqrt(x2_plus), -np.sqrt(x2_plus), np.sqrt(x2_minus), -np.sqrt(x2_minus)])
roots_test = np.sort(roots_test)
print(f"  Test roots (t2={t2_test}): {roots_test}")
print(f"  Sum: {sum(roots_test):.10f}")
phi_test = phi_n_num(roots_test)
print(f"  1/Phi_4 actual: {1/phi_test:.10f}")
C4_val = float(inv_phi4_K0.subs(K2, t2_test))
print(f"  1/Phi_4 formula: {C4_val:.10f}")

# So C_4(ta) is the rational function inv_phi4_K0 evaluated at K2=ta_2.
# If inv_phi4_K0 = alpha*K2 for some constant alpha, then C_4 = alpha.
print(f"\n  Checking if 1/Phi_4(K2,0,0) is proportional to K2:")
for t2_val in [-0.5, -1.0, -2.0, -3.0, -5.0]:
    val = float(inv_phi4_K0.subs(K2, t2_val))
    ratio = val / t2_val
    print(f"    t2={t2_val}: 1/Phi_4 = {val:.10f}, ratio = {ratio:.10f}")
