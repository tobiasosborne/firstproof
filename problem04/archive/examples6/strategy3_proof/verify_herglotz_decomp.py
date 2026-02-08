#!/usr/bin/env python3
"""
Verify the Herglotz decomposition approach to <h,alpha> >= 0.

Key identity: omega_1(z) = z + c + sum_j m_j/(z - p_j) with m_j > 0
=> omega_1''(z)/2 = sum_j m_j/(z-p_j)^3
=> alpha_k = sum_j m_j / (nu_k - p_j)^3

So <h,alpha> = sum_k h_k * alpha_k = sum_j m_j * [sum_k h_k/(nu_k - p_j)^3]
            = sum_j m_j * A_j

where A_j = sum_k H_r(nu_k) / (nu_k - p_j)^3.

If A_j >= 0 for all j, then <h,alpha> >= 0 since m_j > 0.

We test whether A_j >= 0. The pole p_j lies between consecutive roots of r
(interlacing), so we can study the sign of A_j.
"""

import numpy as np
from numpy.polynomial import polynomial as P
from math import factorial

np.random.seed(42)

def boxplus_n(p_coeffs, q_coeffs, n):
    a, b = p_coeffs, q_coeffs
    r = np.zeros(n + 1)
    r[0] = 1.0
    for k in range(1, n + 1):
        c_k = 0.0
        for i in range(0, k + 1):
            j = k - i
            if i <= n and j <= n:
                coeff = (factorial(n - i) * factorial(n - j)) / (factorial(n) * factorial(n - k))
                c_k += coeff * a[i] * b[j]
        r[k] = c_k
    return r

def H_values(roots):
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H

def generate_random_poly(n, spread=5.0):
    roots = np.sort(np.random.uniform(-spread, spread, n))
    for i in range(1, n):
        if roots[i] - roots[i-1] < 0.5:
            roots[i] = roots[i-1] + 0.5 + np.random.uniform(0, 0.5)
    return roots

def compute_A_j(r_roots, p_j):
    """Compute A_j = sum_k H_r(nu_k) / (nu_k - p_j)^3"""
    H_r = H_values(r_roots)
    return np.sum(H_r / (r_roots - p_j)**3)

def compute_B_j(r_roots, p_j):
    """Compute B_j = sum_k H_r(nu_k) / (nu_k - p_j)^2 (for comparison)"""
    H_r = H_values(r_roots)
    return np.sum(H_r / (r_roots - p_j)**2)

print("=" * 70)
print("TESTING: Does A_j = sum_k h_k/(nu_k - p_j)^3 have definite sign")
print("for poles p_j of omega_1 that interlace with roots nu_k of r?")
print("=" * 70)

# First: study A_j as a function of p for general root configurations
for n in [3, 4, 5]:
    print(f"\n{'='*60}")
    print(f"n = {n}")
    print(f"{'='*60}")

    np.random.seed(42 + n)
    count_Aj_neg = 0
    count_Aj_total = 0

    for trial in range(200):
        r_roots = generate_random_poly(n)
        H_r = H_values(r_roots)

        # Test with p_j between consecutive roots (where poles of omega_1 typically are)
        for gap_idx in range(n - 1):
            # p_j between nu_{gap_idx} and nu_{gap_idx+1}
            for frac in [0.1, 0.3, 0.5, 0.7, 0.9]:
                p_j = r_roots[gap_idx] + frac * (r_roots[gap_idx + 1] - r_roots[gap_idx])
                A_j = compute_A_j(r_roots, p_j)
                count_Aj_total += 1
                if A_j < -1e-12:
                    count_Aj_neg += 1

        # Also test p_j outside all roots
        for p_j in [r_roots[0] - 1.0, r_roots[0] - 0.5, r_roots[-1] + 0.5, r_roots[-1] + 1.0]:
            A_j = compute_A_j(r_roots, p_j)
            count_Aj_total += 1
            if A_j < -1e-12:
                count_Aj_neg += 1

    print(f"  A_j < 0 in {count_Aj_neg}/{count_Aj_total} cases")

# Detailed analysis for small n
print("\n" + "=" * 70)
print("DETAILED ANALYSIS FOR n=3")
print("=" * 70)

n = 3
np.random.seed(42)
for trial in range(5):
    r_roots = generate_random_poly(n)
    H_r = H_values(r_roots)

    print(f"\n  Trial {trial}: r_roots = {r_roots}")
    print(f"  H_r = {H_r}")

    # Sweep p through the gaps
    for gap in range(n - 1):
        p_vals = np.linspace(r_roots[gap] + 0.001, r_roots[gap+1] - 0.001, 20)
        A_vals = [compute_A_j(r_roots, p) for p in p_vals]
        signs = ['+' if a > 0 else '-' if a < 0 else '0' for a in A_vals]
        print(f"  Gap [{r_roots[gap]:.2f}, {r_roots[gap+1]:.2f}]: "
              f"A_j signs = {''.join(signs)}")
        print(f"    A_j range: [{min(A_vals):.4f}, {max(A_vals):.4f}]")

    # Outside
    for p_j in [r_roots[0] - 2, r_roots[-1] + 2]:
        A_j = compute_A_j(r_roots, p_j)
        print(f"  p_j = {p_j:.2f} (outside): A_j = {A_j:.6f}")


print("\n" + "=" * 70)
print("PARTIAL FRACTION ANALYSIS: What does A_j actually look like?")
print("=" * 70)

# For n=3 with roots nu_1 < nu_2 < nu_3:
# H_r(nu_1) = 1/(nu_1-nu_2) + 1/(nu_1-nu_3) < 0  (both terms negative)
# H_r(nu_2) = 1/(nu_2-nu_1) + 1/(nu_2-nu_3)       (positive + negative, sign varies)
# H_r(nu_3) = 1/(nu_3-nu_1) + 1/(nu_3-nu_2) > 0  (both terms positive)
#
# A_j(p) = H_r(nu_1)/(nu_1-p)^3 + H_r(nu_2)/(nu_2-p)^3 + H_r(nu_3)/(nu_3-p)^3
#
# For p between nu_1 and nu_2:
# (nu_1-p)^3 < 0 (since nu_1 < p), H_r(nu_1) < 0, so first term > 0
# (nu_2-p)^3 > 0 (since nu_2 > p), H_r(nu_2) = ??, sign varies
# (nu_3-p)^3 > 0 (since nu_3 > p), H_r(nu_3) > 0, so third term > 0

# This shows A_j >= 0 is NOT obvious from sign analysis alone.
# There must be cancellations.

print("\nSign structure of individual terms:")
n = 3
np.random.seed(42)
r_roots = generate_random_poly(n)
H_r = H_values(r_roots)
print(f"  r_roots = {r_roots}, H_r = {H_r}")

for gap in range(n-1):
    p = (r_roots[gap] + r_roots[gap+1]) / 2
    print(f"\n  p = {p:.4f} (midpoint of gap {gap}):")
    for k in range(n):
        term = H_r[k] / (r_roots[k] - p)**3
        print(f"    k={k}: H_r[k]={H_r[k]:.4f}, (nu_k-p)^3={(r_roots[k]-p)**3:.4f}, "
              f"term={term:.4f}")
    A_j = compute_A_j(r_roots, p)
    print(f"    A_j = {A_j:.6f}")


print("\n" + "=" * 70)
print("ALTERNATIVE: Consider A_j as derivative of a known-sign function")
print("=" * 70)

# Note that A_j(p) = sum_k h_k / (nu_k - p)^3 = -(1/2) d/dp [sum_k h_k / (nu_k-p)^2]
# = -(1/2) d/dp [B_j(p)]
# where B_j(p) = sum_k h_k / (nu_k - p)^2
#
# Now B_j(p) = sum_k h_k / (nu_k-p)^2
# This is related to: d/dp [sum_k h_k / (nu_k-p)] = sum_k h_k / (nu_k-p)^2 = B_j(p)
# And sum_k h_k / (nu_k-p) = sum_k [sum_{l!=k} 1/(nu_k-nu_l)] / (nu_k-p)
#
# Consider the identity:
# sum_k 1/((nu_k-p)(nu_k-nu_l)) = [1/(nu_l-p)] * sum_k 1/(nu_k-p) - correction?
# Actually by partial fractions:
# 1/((x-a)(x-b)) = [1/(a-b)] * [1/(x-a) - 1/(x-b)]
# So 1/((nu_k-p)(nu_k-nu_l)) = [1/(p-nu_l)] * [1/(nu_k-p) - 1/(nu_k-nu_l)]
# when l != k.
#
# sum_{l!=k} 1/((nu_k-p)(nu_k-nu_l)) = 1/(nu_k-p) * h_k
# So h_k/(nu_k-p) = sum_{l!=k} 1/((nu_k-nu_l)(nu_k-p))  ... this is just definition.
#
# sum_k h_k/(nu_k-p) = sum_k sum_{l!=k} 1/((nu_k-nu_l)(nu_k-p))
#
# Let's use partial fractions:
# 1/((nu_k-nu_l)(nu_k-p)) = [1/(nu_l-p)] * [1/(nu_k-nu_l) - 1/(nu_k-p)] ... no
# Actually: 1/((a-b)(a-c)) = 1/((a-b)(a-c))
# partial fraction in a: not useful since a = nu_k varies.
#
# Partial fraction identity:
# h_k / (nu_k - p) = [sum_{l!=k} 1/(nu_k-nu_l)] / (nu_k-p)
#
# sum_k h_k/(nu_k-p) = sum_k sum_{l!=k} 1/((nu_k-nu_l)*(nu_k-p))
#
# Swap sums: = sum_{k,l: k!=l} 1/((nu_k-nu_l)*(nu_k-p))
# = sum_{k<l} [1/((nu_k-nu_l)*(nu_k-p)) + 1/((nu_l-nu_k)*(nu_l-p))]
# = sum_{k<l} 1/(nu_k-nu_l) * [1/(nu_k-p) - 1/(nu_l-p)]
# = sum_{k<l} 1/(nu_k-nu_l) * (nu_l-nu_k)/((nu_k-p)(nu_l-p))
# = -sum_{k<l} 1/((nu_k-p)(nu_l-p))
# = -(1/2) * [sum_k 1/(nu_k-p)]^2 + (1/2)*sum_k 1/(nu_k-p)^2
# = -(1/2)*[nG_r(p)]^2*n^2 + (1/2)*sum_k 1/(nu_k-p)^2
# Wait let me redo this:
# [sum_k 1/(nu_k-p)]^2 = sum_k 1/(nu_k-p)^2 + 2*sum_{k<l} 1/((nu_k-p)(nu_l-p))
# So sum_{k<l} 1/((nu_k-p)(nu_l-p)) = (1/2)*{[sum 1/(nu_k-p)]^2 - sum 1/(nu_k-p)^2}
# => sum_k h_k/(nu_k-p) = -sum_{k<l} 1/((nu_k-p)(nu_l-p))
#    = -(1/2)*{[nG_r(p)]^2*n^2 - sum_k 1/(nu_k-p)^2}
#    Hmm wait, nG_r(p) = (1/n)*sum_k 1/(p-nu_k) = -(1/n)*sum_k 1/(nu_k-p)
#    So sum_k 1/(nu_k-p) = -n*G_r(p)*n = ... no, nG_r(p) = sum_k 1/(p-nu_k)/1
#    Actually nG_r(z) = r'(z)/r(z) = sum_k 1/(z-nu_k)
#    So sum_k 1/(nu_k-p) = -sum_k 1/(p-nu_k) = -nG_r(p)  (where nG_r = r'/r = sum 1/(z-nu_k))
#
# OK: sum_k 1/(nu_k-p) = -nG_r(p)
# [sum_k 1/(nu_k-p)]^2 = [nG_r(p)]^2
# sum_{k<l} 1/((nu_k-p)(nu_l-p)) = (1/2)*{[nG_r(p)]^2 - sum_k 1/(nu_k-p)^2}
#
# And sum_k h_k/(nu_k-p) = -sum_{k<l} 1/((nu_k-p)(nu_l-p))
#                         = -(1/2)*{[nG_r(p)]^2 - sum_k 1/(nu_k-p)^2}
#                         = (1/2)*{sum_k 1/(nu_k-p)^2 - [nG_r(p)]^2}
#                         = (1/2)*{-[nG_r]'(p) - [nG_r(p)]^2}
#                         (using [nG_r]'(z) = -sum 1/(z-nu_k)^2)
#                         = -(1/2)*{[nG_r]'(p) + [nG_r(p)]^2}
#                         = -(1/2) * d/dp [nG_r(p)] + ??? hmm
#
# Actually [nG_r]' = -sum 1/(z-nu_k)^2 and sum 1/(nu_k-p)^2 = sum 1/(p-nu_k)^2 = -[nG_r]'(p)
# So: sum_k h_k/(nu_k-p) = (1/2)*{-[nG_r]'(p) - [nG_r(p)]^2}
# = -(1/2)*{[nG_r(p)]^2 + [nG_r]'(p)}
# = -(1/2) * d/dp [nG_r(p)] ... no, d/dp [nG_r(p)] = [nG_r]'(p)
# (nG_r)^2 + (nG_r)' is not a simple derivative.
#
# But note: d/dp [nG_r(p)] = [nG_r]'(p) = -sum 1/(p-nu_k)^2
# and (nG_r(p))^2 = [sum 1/(p-nu_k)]^2
#
# KEY IDENTITY:
# sum_k h_k/(nu_k-p) = -(1/2)*[(nG_r(p))^2 + (nG_r)'(p)]
#                     = -(1/2)*[r'(p)/r(p)]^2 - (1/2)*[r'(p)/r(p)]'
#                     = -(1/2)*[r'(p)/r(p)]^2 - (1/2)*[r''(p)*r(p) - r'(p)^2]/r(p)^2
#                     = -(1/2)*r'(p)^2/r(p)^2 - (1/2)*r''(p)/r(p) + (1/2)*r'(p)^2/r(p)^2
#                     = -(1/2)*r''(p)/r(p)
#
# BEAUTIFUL!
# sum_k h_k/(nu_k-p) = -(1/2)*r''(p)/r(p)

# Let me verify this:
n = 4
np.random.seed(42)
r_roots = generate_random_poly(n)
H_r = H_values(r_roots)
r_poly = np.poly(r_roots)
r_deriv = np.polyder(r_poly)
r_deriv2 = np.polyder(r_deriv)

p_test = r_roots[0] + 0.3 * (r_roots[1] - r_roots[0])
lhs = np.sum(H_r / (r_roots - p_test))
rhs = -0.5 * np.polyval(r_deriv2, p_test) / np.polyval(r_poly, p_test)
print(f"\nVerifying sum_k h_k/(nu_k-p) = -(1/2)*r''(p)/r(p)")
print(f"  LHS = {lhs:.10f}")
print(f"  RHS = {rhs:.10f}")
print(f"  Match: {abs(lhs - rhs) < 1e-8}")

# More tests
for n in [3, 4, 5]:
    np.random.seed(42+n)
    r_roots = generate_random_poly(n)
    H_r = H_values(r_roots)
    r_poly = np.poly(r_roots)
    r_deriv2 = np.polyder(np.polyder(r_poly))

    all_match = True
    for _ in range(20):
        p_test = np.random.uniform(r_roots[0] - 3, r_roots[-1] + 3)
        # Avoid being too close to a root
        if min(abs(p_test - r_roots)) < 0.01:
            continue
        lhs = np.sum(H_r / (r_roots - p_test))
        rhs = -0.5 * np.polyval(r_deriv2, p_test) / np.polyval(r_poly, p_test)
        if abs(lhs - rhs) > 1e-6:
            all_match = False
            print(f"  MISMATCH at n={n}, p={p_test:.4f}: LHS={lhs:.8f}, RHS={rhs:.8f}")
    if all_match:
        print(f"  n={n}: Identity verified for 20 random p values")

print("\n" + "=" * 70)
print("KEY IDENTITY PROVED:")
print("sum_k H_r(nu_k) / (nu_k - p) = -(1/2) * r''(p) / r(p)")
print("=" * 70)

# Now differentiate with respect to p:
# d/dp [sum_k h_k/(nu_k-p)] = sum_k h_k/(nu_k-p)^2 = B_j(p)
# d/dp [-(1/2)*r''(p)/r(p)] = -(1/2)*[r'''(p)*r(p) - r''(p)*r'(p)] / r(p)^2
#                             = -(1/2)*r'''(p)/r(p) + (1/2)*r''(p)*r'(p)/r(p)^2
#
# So B_j(p) = sum_k h_k/(nu_k-p)^2 = -(1/2)*r'''(p)/r(p) + (1/2)*r''(p)*r'(p)/r(p)^2

# Differentiate once more:
# A_j(p) = sum_k h_k/(nu_k-p)^3 = (1/2) * d/dp [sum_k h_k/(nu_k-p)^2]
# Wait: d/dp [1/(nu_k-p)^2] = 2/(nu_k-p)^3, so d/dp sum h_k/(nu_k-p)^2 = 2*sum h_k/(nu_k-p)^3 = 2*A_j
# So A_j = (1/2) * d/dp [B_j]

# But B_j = d/dp [-(1/2)*r''/r]
# So A_j = (1/2) * d^2/dp^2 [-(1/2)*r''/r] = -(1/4)*d^2/dp^2 [r''/r]

# Let F(p) = r''(p)/r(p). Then:
# sum_k h_k/(nu_k-p) = -(1/2)*F(p)
# sum_k h_k/(nu_k-p)^2 = (1/2)*F'(p)    (differentiate, noting d/dp(nu_k-p)^{-1} = (nu_k-p)^{-2})
# sum_k h_k/(nu_k-p)^3 = -(1/4)*F''(p)   (differentiate again)

# So A_j(p) = -(1/4)*F''(p) where F(p) = r''(p)/r(p)

# Verify:
print("\nVerifying A_j(p) = -(1/4)*F''(p) where F(p) = r''(p)/r(p)")

for n in [3, 4, 5]:
    np.random.seed(42+n)
    r_roots = generate_random_poly(n)
    H_r = H_values(r_roots)
    r_poly = np.poly(r_roots)
    r_d1 = np.polyder(r_poly)
    r_d2 = np.polyder(r_d1)
    r_d3 = np.polyder(r_d2)
    r_d4 = np.polyder(r_d3)

    all_match = True
    for _ in range(20):
        p_test = np.random.uniform(r_roots[0] - 3, r_roots[-1] + 3)
        if min(abs(p_test - r_roots)) < 0.1:
            continue

        A_j_direct = np.sum(H_r / (r_roots - p_test)**3)

        # Compute F''(p) numerically
        h = 1e-5
        F_plus = np.polyval(r_d2, p_test + h) / np.polyval(r_poly, p_test + h)
        F_0 = np.polyval(r_d2, p_test) / np.polyval(r_poly, p_test)
        F_minus = np.polyval(r_d2, p_test - h) / np.polyval(r_poly, p_test - h)
        F_pp = (F_plus - 2*F_0 + F_minus) / h**2

        A_j_formula = -0.25 * F_pp

        if abs(A_j_direct - A_j_formula) > abs(A_j_direct) * 0.01 + 1e-4:
            all_match = False
            print(f"  MISMATCH at n={n}, p={p_test:.4f}: "
                  f"direct={A_j_direct:.6f}, formula={A_j_formula:.6f}")

    if all_match:
        print(f"  n={n}: A_j = -(1/4)*F'' verified for 20 random p values")

print("\n" + "=" * 70)
print("KEY CONSEQUENCE: <h,alpha> = sum_j m_j * A_j = -(1/4)*sum_j m_j * F''(p_j)")
print("where F(p) = r''(p)/r(p) and p_j are poles of omega_1 with residues m_j > 0")
print("=" * 70)

# Now: F(p) = r''(p)/r(p)
# Near a root nu_k of r: r(nu_k) = 0, r'(nu_k) != 0
# F(p) ~ r''(nu_k) / [r'(nu_k) * (p - nu_k)] = 2*H_r(nu_k) / (p - nu_k)
# So F has simple poles at nu_k with residue 2*H_r(nu_k) = 2*h_k
#
# F is a rational function: degree (n-2) / degree n = O(1/p^2) at infinity.
# F(p) = 2*sum_k h_k/(p - nu_k) + [polynomial part?]
# Wait: r''(p) has degree n-2, r(p) has degree n, so F = O(1/p^2) at infinity.
# The partial fraction of F is EXACTLY:
# F(p) = sum_k [r''(nu_k)/r'(nu_k)] / (p - nu_k) = sum_k 2*h_k / (p - nu_k)
# (since the degree of numerator < degree of denominator, there's no polynomial part)
#
# This is consistent with our identity: sum_k h_k/(nu_k-p) = -(1/2)*F(p)
# => F(p) = -2*sum_k h_k/(nu_k-p) = 2*sum_k h_k/(p-nu_k) ✓

# Now F''(p) = 2*sum_k 2*h_k / (p-nu_k)^3 = 4*sum_k h_k/(p-nu_k)^3
# And A_j = -(1/4)*F''(p_j) = -sum_k h_k/(p_j-nu_k)^3 = sum_k h_k/(nu_k-p_j)^3 ✓

# So A_j >= 0 iff sum_k h_k/(nu_k-p_j)^3 >= 0
# or equivalently F''(p_j) <= 0.

# This means: <h,alpha> >= 0 iff F(p) = r''(p)/r(p) is CONCAVE at the poles p_j of omega_1.

print("\nQuestion reduced to: Is F(p) = r''(p)/r(p) concave at the poles of omega_1?")
print("F''(p_j) <= 0 for all poles p_j of omega_1?")

# Check this numerically
for n in [3, 4, 5]:
    np.random.seed(42+n)
    count_concave = 0
    count_total = 0

    for trial in range(200):
        r_roots = generate_random_poly(n)
        H_r = H_values(r_roots)
        r_poly = np.poly(r_roots)
        r_d2 = np.polyder(np.polyder(r_poly))

        # Check F''(p) at midpoints of gaps (proxy for poles of omega_1)
        for gap in range(n-1):
            p_mid = (r_roots[gap] + r_roots[gap+1]) / 2
            h = 1e-5
            F_plus = np.polyval(r_d2, p_mid+h) / np.polyval(r_poly, p_mid+h)
            F_0 = np.polyval(r_d2, p_mid) / np.polyval(r_poly, p_mid)
            F_minus = np.polyval(r_d2, p_mid-h) / np.polyval(r_poly, p_mid-h)
            F_pp = (F_plus - 2*F_0 + F_minus) / h**2

            count_total += 1
            if F_pp <= 1e-6:  # concave
                count_concave += 1

    print(f"  n={n}: F''(p_mid) <= 0 in {count_concave}/{count_total} cases")

# IMPORTANT: The midpoints of gaps are NOT the actual poles of omega_1.
# The poles depend on p and q, not just r. But let's see the general picture.

# Actually, let's check for ALL p in the gaps, not just midpoints
print("\nChecking sign of A_j = sum_k h_k/(nu_k-p)^3 for ALL p between roots:")
for n in [3, 4, 5]:
    np.random.seed(42+n)
    count_neg = 0
    count_total = 0

    for trial in range(100):
        r_roots = generate_random_poly(n)
        H_r = H_values(r_roots)

        for gap in range(n-1):
            p_vals = np.linspace(r_roots[gap]+0.001, r_roots[gap+1]-0.001, 50)
            for p in p_vals:
                A = np.sum(H_r / (r_roots - p)**3)
                count_total += 1
                if A < -1e-10:
                    count_neg += 1

    print(f"  n={n}: A_j < 0 in {count_neg}/{count_total} cases")

print("\n" + "=" * 70)
print("Checking A_j for p OUTSIDE the root range too:")
for n in [3, 4, 5]:
    np.random.seed(42+n)
    count_neg_in = 0
    count_neg_out = 0
    count_in = 0
    count_out = 0

    for trial in range(100):
        r_roots = generate_random_poly(n)
        H_r = H_values(r_roots)

        # Outside roots: p < nu_1 or p > nu_n
        for p in np.linspace(r_roots[0] - 5, r_roots[0] - 0.01, 20):
            A = np.sum(H_r / (r_roots - p)**3)
            count_out += 1
            if A < -1e-10:
                count_neg_out += 1

        for p in np.linspace(r_roots[-1] + 0.01, r_roots[-1] + 5, 20):
            A = np.sum(H_r / (r_roots - p)**3)
            count_out += 1
            if A < -1e-10:
                count_neg_out += 1

    print(f"  n={n}: A_j < 0 outside roots: {count_neg_out}/{count_out}")
