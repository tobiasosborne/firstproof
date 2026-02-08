#!/usr/bin/env python3
"""
Verify properties of D_{kl} = (lambda_k - lambda_l)/(nu_k - nu_l)
and the divided difference identity for <h,alpha>.
"""
import numpy as np
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


print("=" * 70)
print("COMPREHENSIVE CHECK OF D_{kl} PROPERTIES")
print("=" * 70)

for n in [3, 4, 5]:
    np.random.seed(42 + n)

    count_D_lt_1 = 0
    count_D_gt_1 = 0
    count_total = 0
    min_D = float('inf')
    max_D = float('-inf')

    # Adjacent vs non-adjacent
    count_adj_lt_1 = 0
    count_adj_gt_1 = 0
    count_adj_total = 0
    count_nonadj_lt_1 = 0
    count_nonadj_gt_1 = 0
    count_nonadj_total = 0

    for trial in range(500):
        p_roots = generate_random_poly(n)
        q_roots = generate_random_poly(n)
        p_coeffs = np.poly(p_roots)
        q_coeffs = np.poly(q_roots)
        r_coeffs = boxplus_n(p_coeffs, q_coeffs, n)
        r_roots_raw = np.roots(r_coeffs)
        if np.max(np.abs(np.imag(r_roots_raw))) > 1e-6:
            continue
        r_roots = np.sort(np.real(r_roots_raw))
        if np.min(np.diff(r_roots)) < 1e-8:
            continue

        for k in range(n):
            for l in range(k+1, n):
                D = (p_roots[k] - p_roots[l]) / (r_roots[k] - r_roots[l])
                count_total += 1
                if D < 1 - 1e-10:
                    count_D_lt_1 += 1
                elif D > 1 + 1e-10:
                    count_D_gt_1 += 1
                min_D = min(min_D, D)
                max_D = max(max_D, D)

                if l == k + 1:  # adjacent
                    count_adj_total += 1
                    if D < 1 - 1e-10:
                        count_adj_lt_1 += 1
                    elif D > 1 + 1e-10:
                        count_adj_gt_1 += 1
                else:  # non-adjacent
                    count_nonadj_total += 1
                    if D < 1 - 1e-10:
                        count_nonadj_lt_1 += 1
                    elif D > 1 + 1e-10:
                        count_nonadj_gt_1 += 1

    print(f"\nn = {n}:")
    print(f"  ALL pairs: D < 1 in {count_D_lt_1}/{count_total}, D > 1 in {count_D_gt_1}/{count_total}")
    print(f"  D range: [{min_D:.6f}, {max_D:.6f}]")
    print(f"  Adjacent: D < 1 in {count_adj_lt_1}/{count_adj_total}, D > 1 in {count_adj_gt_1}/{count_adj_total}")
    print(f"  Non-adj:  D < 1 in {count_nonadj_lt_1}/{count_nonadj_total}, D > 1 in {count_nonadj_gt_1}/{count_nonadj_total}")


print("\n" + "=" * 70)
print("UNDERSTANDING WHY D_{kl} < 1:")
print("The roots of r = p boxplus_n q are MORE SPREAD than roots of p.")
print("So |nu_k - nu_l| > |lambda_k - lambda_l| for most pairs.")
print("This means D_{kl} = |lambda_k - lambda_l|/|nu_k - nu_l| < 1.")
print("=" * 70)

# Verify: are roots of r always more spread?
for n in [3, 4, 5]:
    np.random.seed(42 + n)
    count_more_spread = 0
    count_total = 0
    for trial in range(500):
        p_roots = generate_random_poly(n)
        q_roots = generate_random_poly(n)
        p_coeffs = np.poly(p_roots)
        q_coeffs = np.poly(q_roots)
        r_coeffs = boxplus_n(p_coeffs, q_coeffs, n)
        r_roots_raw = np.roots(r_coeffs)
        if np.max(np.abs(np.imag(r_roots_raw))) > 1e-6:
            continue
        r_roots = np.sort(np.real(r_roots_raw))
        if np.min(np.diff(r_roots)) < 1e-8:
            continue

        count_total += 1
        if (r_roots[-1] - r_roots[0]) >= (p_roots[-1] - p_roots[0]):
            count_more_spread += 1

    print(f"  n={n}: r more spread than p in {count_more_spread}/{count_total}")

# The key insight: even though D_{kl} < 1, the H_r divided differences
# have the right sign pattern to make the weighted sum positive.

print("\n" + "=" * 70)
print("STRUCTURE OF H_r DIVIDED DIFFERENCES")
print("=" * 70)

# H_r[k,l] = (h_k - h_l)/(nu_k - nu_l)
# For k < l: nu_k < nu_l
# h_k tends to go from negative (k=1) to positive (k=n)
# So h_k - h_l < 0 and nu_k - nu_l < 0, giving H_r[k,l] > 0... usually.
#
# But H_r is NOT always monotonically increasing! For n >= 4, middle roots
# can have H values that are not monotone.

for n in [3, 4, 5]:
    np.random.seed(42 + n)
    count_monotone = 0
    count_total = 0
    for trial in range(500):
        r_roots = generate_random_poly(n)
        H_r = H_values(r_roots)
        count_total += 1
        if all(H_r[i] <= H_r[i+1] for i in range(n-1)):
            count_monotone += 1
    print(f"  n={n}: H_r monotone increasing in {count_monotone}/{count_total}")

# For n=3: H_r is always h_1 < h_2 < h_3 (since there are only 3 roots
# and the Hilbert transform is monotonically increasing for 3 points).
# For n >= 4: not always monotone.

# KEY STRUCTURAL RESULT:
# For n=3, H_r is always monotone, so ALL divided differences H_r[k,l] > 0,
# and since 1/D_{kl} - 1 > 0, we get <h,alpha> = sum of positive terms >= 0.
# This PROVES <h,alpha> >= 0 for n=3!

print("\n" + "=" * 70)
print("PROOF FOR n=3: H_r monotone + D < 1 implies <h,alpha> >= 0")
print("=" * 70)

n = 3
np.random.seed(42)
count_ha_neg = 0
count_total = 0
for trial in range(2000):
    p_roots = generate_random_poly(n)
    q_roots = generate_random_poly(n)
    p_coeffs = np.poly(p_roots)
    q_coeffs = np.poly(q_roots)
    r_coeffs = boxplus_n(p_coeffs, q_coeffs, n)
    r_roots_raw = np.roots(r_coeffs)
    if np.max(np.abs(np.imag(r_roots_raw))) > 1e-6:
        continue
    r_roots = np.sort(np.real(r_roots_raw))
    if np.min(np.diff(r_roots)) < 1e-8:
        continue

    H_r = H_values(r_roots)
    H_p = H_values(p_roots)
    ha = np.dot(H_r, H_p - H_r)
    count_total += 1
    if ha < -1e-10:
        count_ha_neg += 1

    # Also verify H_r monotone
    assert H_r[0] <= H_r[1] <= H_r[2], f"H_r not monotone: {H_r}"

print(f"n=3: <h,alpha> < 0 in {count_ha_neg}/{count_total} trials")
print("AND H_r is always monotone for n=3.")
print()

# For n=3: sum_{k<l} H_r[k,l] * (1/D_{kl} - 1) with all factors >= 0
# gives <h,alpha> >= 0. QED for n=3.

# Now let's check: for which n does H_r fail to be monotone?
# And when it fails, does the identity still give >= 0?

print("Checking n=4:")
n = 4
np.random.seed(42)
count_monotone = 0
count_total = 0
count_neg_dd = 0
count_neg_dd_term_dominates = 0

for trial in range(2000):
    p_roots = generate_random_poly(n)
    q_roots = generate_random_poly(n)
    p_coeffs = np.poly(p_roots)
    q_coeffs = np.poly(q_roots)
    r_coeffs = boxplus_n(p_coeffs, q_coeffs, n)
    r_roots_raw = np.roots(r_coeffs)
    if np.max(np.abs(np.imag(r_roots_raw))) > 1e-6:
        continue
    r_roots = np.sort(np.real(r_roots_raw))
    if np.min(np.diff(r_roots)) < 1e-8:
        continue

    H_r = H_values(r_roots)
    count_total += 1

    if all(H_r[i] <= H_r[i+1] for i in range(n-1)):
        count_monotone += 1
    else:
        count_neg_dd += 1
        # Check which divided differences are negative
        neg_terms = []
        pos_terms = []
        for k in range(n):
            for l in range(k+1, n):
                D = (p_roots[k] - p_roots[l]) / (r_roots[k] - r_roots[l])
                Hr_dd = (H_r[k] - H_r[l]) / (r_roots[k] - r_roots[l])
                term = Hr_dd * (1.0/D - 1.0)
                if term < 0:
                    neg_terms.append(term)
                else:
                    pos_terms.append(term)

        if sum(neg_terms) + sum(pos_terms) < 0:
            count_neg_dd_term_dominates += 1

print(f"n=4: H_r monotone in {count_monotone}/{count_total}")
print(f"     H_r non-monotone: {count_neg_dd}")
print(f"     Non-monotone AND <h,alpha> < 0: {count_neg_dd_term_dominates}")

# Prove H_r monotone for n=3 analytically:
# For n=3: H_r(nu_k) = sum_{j!=k} 1/(nu_k - nu_j)
# h_1 = 1/(nu_1-nu_2) + 1/(nu_1-nu_3) < 0 (both negative since nu_1 < nu_2, nu_3)
# h_2 = 1/(nu_2-nu_1) + 1/(nu_2-nu_3) = positive + negative
# h_3 = 1/(nu_3-nu_1) + 1/(nu_3-nu_2) > 0 (both positive)
# h_1 < 0 < h_3 always.
# h_2 = 1/(nu_2-nu_1) + 1/(nu_2-nu_3) = 1/(nu_2-nu_1) - 1/(nu_3-nu_2)
# h_2 > 0 iff 1/(nu_2-nu_1) > 1/(nu_3-nu_2) iff nu_3-nu_2 > nu_2-nu_1
# i.e., nu_2 is closer to nu_1 than to nu_3. Not always true.
# So h_2 can be positive or negative.
#
# h_2 - h_1 = [1/(nu_2-nu_1) + 1/(nu_2-nu_3)] - [1/(nu_1-nu_2) + 1/(nu_1-nu_3)]
# = 1/(nu_2-nu_1) + 1/(nu_2-nu_3) + 1/(nu_2-nu_1) - 1/(nu_1-nu_3)
# (using 1/(nu_1-nu_2) = -1/(nu_2-nu_1))
# = 2/(nu_2-nu_1) + 1/(nu_2-nu_3) + 1/(nu_3-nu_1)
# (using 1/(nu_1-nu_3) = -1/(nu_3-nu_1))
# = 2/(nu_2-nu_1) - 1/(nu_3-nu_2) + 1/(nu_3-nu_1)
#
# Let a = nu_2 - nu_1 > 0, b = nu_3 - nu_2 > 0. Then nu_3 - nu_1 = a + b.
# h_2 - h_1 = 2/a - 1/b + 1/(a+b)
# = 2/a + (1/(a+b) - 1/b) = 2/a - a/(b(a+b))
# = 2/a - a/(b(a+b))
# = (2*b*(a+b) - a^2) / (a*b*(a+b))
# = (2*a*b + 2*b^2 - a^2) / (a*b*(a+b))
#
# For this to be >= 0: 2ab + 2b^2 >= a^2, i.e., 2b(a+b) >= a^2.
# This is NOT always true. E.g., a = 10, b = 1: 2*1*11 = 22, a^2 = 100. FAILS.
#
# WAIT: But h_2 - h_1 should ALWAYS be >= 0 for n=3 based on our numerical check!
# Let me verify.

print("\n\nAnalytical check for n=3:")
print("h_2 - h_1 = (2ab + 2b^2 - a^2) / (ab(a+b))")
print("This is >= 0 iff a^2 <= 2b(a+b) = 2ab + 2b^2")
print("i.e., a^2 - 2ab - 2b^2 <= 0")
print("Discriminant: 4b^2 + 8b^2 = 12b^2, so a <= (2b + 2b*sqrt(3))/2 = b(1+sqrt(3))")
print(f"Critical ratio: a/b <= 1 + sqrt(3) = {1+np.sqrt(3):.6f}")
print()

# So h_2 - h_1 >= 0 iff a/b <= 1 + sqrt(3) ~ 2.732
# This CAN fail! E.g., a/b = 3: h_2 - h_1 = (6b^2 + 2b^2 - 9b^2)/(3b^3(4)) = -b^2/(12b^3) < 0

# Let me check:
a, b = 3.0, 1.0
nu = np.array([0.0, a, a+b])
H = H_values(nu)
print(f"a=3, b=1: nu = {nu}, H = {H}")
print(f"h_2 - h_1 = {H[1] - H[0]:.6f}")
print(f"Expected: (2*3*1 + 2*1 - 9)/(3*1*4) = {(6+2-9)/(3*4):.6f}")

# So H_r is NOT always monotone for n=3!
# My earlier numerical check was wrong because I only used generate_random_poly
# which constrains minimum spacing.

# Let me re-check with wider spacing ratios:
print("\n\nRe-checking H_r monotonicity for n=3 with extreme spacing:")
count_non_mono = 0
count_total = 0
for trial in range(10000):
    # Generate with varying spacing ratios
    a = np.random.exponential(2.0)
    b = np.random.exponential(2.0)
    base = np.random.uniform(-10, 10)
    nu = np.array([base, base+a, base+a+b])
    H = H_values(nu)
    count_total += 1
    if not (H[0] <= H[1] <= H[2]):
        count_non_mono += 1

print(f"n=3: H_r non-monotone in {count_non_mono}/{count_total}")

# AH HA! So H_r is NOT always monotone for n=3.
# The divided difference proof for n=3 does NOT work as stated!

# But <h,alpha> is STILL always >= 0. Let me understand why.

# For n=3 with a/b > 1+sqrt(3), we have h_1 > h_2 (h goes UP then DOWN or similar).
# But the divided difference identity still gives >= 0.
# This means the positive terms outweigh the negative.

print("\n\nChecking divided difference sign pattern for non-monotone H_r:")
n = 3
count_all_pos = 0
count_some_neg = 0
count_sum_neg = 0
for trial in range(10000):
    a = np.random.exponential(2.0) + 0.01
    b = np.random.exponential(2.0) + 0.01
    base = np.random.uniform(-10, 10)
    nu = np.array([base, base+a, base+a+b])
    H = H_values(nu)

    # Use random lambda with D < 1
    # lambda must be order-preserving with smaller gaps
    scale = np.random.uniform(0.3, 0.99)
    center = np.mean(nu)
    lam = center + scale * (nu - center) + np.random.uniform(-0.1, 0.1, 3)
    lam = np.sort(lam)

    # Check D < 1 for all pairs
    all_D_ok = True
    for k in range(3):
        for l in range(k+1, 3):
            D = (lam[k]-lam[l])/(nu[k]-nu[l])
            if D >= 1 or D <= 0:
                all_D_ok = False
    if not all_D_ok:
        continue

    H_lam = H_values(lam)  # stand-in for H_p
    alpha = H_lam - H

    ha = np.dot(H, alpha)
    if ha < -1e-8:
        count_sum_neg += 1

    # Check sign of terms
    some_neg = False
    for k in range(3):
        for l in range(k+1, 3):
            D = (lam[k]-lam[l])/(nu[k]-nu[l])
            Hr_dd = (H[k]-H[l])/(nu[k]-nu[l])
            term = Hr_dd * (1.0/D - 1)
            if term < -1e-10:
                some_neg = True
    if some_neg:
        count_some_neg += 1
    else:
        count_all_pos += 1

print(f"Random (nu, lambda) pairs with 0 < D < 1:")
print(f"  All terms positive: {count_all_pos}")
print(f"  Some terms negative: {count_some_neg}")
print(f"  Sum negative: {count_sum_neg}")

# CRITICAL TEST: For actual MSS convolution (not random lambda),
# is H_r monotone?

print("\n\nFinal check: Is H_r monotone for ACTUAL MSS convolution at n=3?")
n = 3
count_non_mono = 0
count_total = 0
for trial in range(5000):
    p_roots = np.sort(np.random.uniform(-10, 10, n))
    while np.min(np.diff(p_roots)) < 0.01:
        p_roots = np.sort(np.random.uniform(-10, 10, n))
    q_roots = np.sort(np.random.uniform(-10, 10, n))
    while np.min(np.diff(q_roots)) < 0.01:
        q_roots = np.sort(np.random.uniform(-10, 10, n))

    p_coeffs = np.poly(p_roots)
    q_coeffs = np.poly(q_roots)
    r_coeffs = boxplus_n(p_coeffs, q_coeffs, n)
    r_roots_raw = np.roots(r_coeffs)
    if np.max(np.abs(np.imag(r_roots_raw))) > 1e-6:
        continue
    r_roots = np.sort(np.real(r_roots_raw))
    if np.min(np.diff(r_roots)) < 1e-8:
        continue

    H_r = H_values(r_roots)
    count_total += 1
    if not (H_r[0] <= H_r[1] + 1e-10 and H_r[1] <= H_r[2] + 1e-10):
        count_non_mono += 1
        if count_non_mono <= 5:
            print(f"  Non-monotone: r_roots={r_roots}, H_r={H_r}")
            a = r_roots[1]-r_roots[0]
            b = r_roots[2]-r_roots[1]
            print(f"    a/b = {a/b:.4f} (threshold: {1+np.sqrt(3):.4f})")

print(f"\nn=3 MSS: H_r non-monotone in {count_non_mono}/{count_total}")
print(f"Ratio: {count_non_mono/max(count_total,1):.4f}")
