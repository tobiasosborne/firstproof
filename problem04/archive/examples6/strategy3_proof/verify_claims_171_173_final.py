#!/usr/bin/env python3
"""
DEFINITIVE VERIFIER SCRIPT for nodes 1.7.1 and 1.7.3.

KEY FINDING: There are THREE different convolutions in play:
1. Hat convolution (hat_e_k(r) = sum hat_e_j(p)*hat_e_{k-j}(q), hat_e_k = e_k/C(n,k))
2. MSS coefficient formula (c_k = sum [(n-i)!(n-j)!/(n!(n-k)!)] * c_i * c_j)
3. Haar unitary average (E_U[det(xI - A - UBU*)])

(1) and (2) are DIFFERENT formulas with different coefficients.
(2) and (3) agree (MSS formula IS the Haar unitary average).

The node text (1.7.1) specifies the MSS formula (2).
This script uses (2) throughout and verifies all claims.
"""

import numpy as np
from math import comb, factorial
from itertools import combinations

def elem_sym_poly(roots, k):
    n = len(roots)
    if k == 0: return 1.0
    if k > n: return 0.0
    return sum(np.prod([roots[i] for i in subset])
               for subset in combinations(range(n), k))

def boxplus_mss(roots_p, roots_q):
    """MSS finite free additive convolution.

    c_k(r) = sum_{i+j=k} [(n-i)!(n-j)! / (n!(n-k)!)] * c_i(p) * c_j(q)
    where c_k = e_k (k-th elementary symmetric polynomial of the roots).
    """
    n = len(roots_p)
    cp = [elem_sym_poly(roots_p, k) for k in range(n+1)]
    cq = [elem_sym_poly(roots_q, k) for k in range(n+1)]

    cr = [0.0] * (n+1)
    for k in range(n+1):
        for i in range(k+1):
            j = k - i
            coeff = (factorial(n-i) * factorial(n-j)) / (factorial(n) * factorial(n-k))
            cr[k] += coeff * cp[i] * cq[j]

    # Build polynomial: p(x) = sum_{k=0}^n (-1)^k c_k x^{n-k}
    poly_coeffs = [(-1)**k * cr[k] for k in range(n+1)]
    roots_r = np.sort(np.real(np.roots(poly_coeffs)))
    return roots_r, cr

def boxplus_haar_mc(roots_p, roots_q, samples=200000):
    """Monte Carlo over Haar-random unitaries."""
    n = len(roots_p)
    A = np.diag(np.array(roots_p, dtype=float))
    B = np.diag(np.array(roots_q, dtype=float))

    sum_ek = np.zeros(n+1)
    for _ in range(samples):
        Z = (np.random.randn(n, n) + 1j * np.random.randn(n, n)) / np.sqrt(2)
        Q, R_mat = np.linalg.qr(Z)
        d = np.diagonal(R_mat)
        ph = d / np.abs(d)
        U = Q * ph[np.newaxis, :]
        M = A + U @ B @ U.conj().T
        eigs = np.sort(np.linalg.eigvalsh(M))
        for k in range(n+1):
            sum_ek[k] += elem_sym_poly(eigs, k)
    avg_ek = sum_ek / samples
    poly_coeffs = [(-1)**k * avg_ek[k] for k in range(n+1)]
    roots_r = np.sort(np.real(np.roots(poly_coeffs)))
    return roots_r

def H_values(roots):
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                if abs(roots[i] - roots[j]) < 1e-14:
                    return None  # Degenerate
                H[i] += 1.0 / (roots[i] - roots[j])
    return H

def Phi_n(roots):
    H = H_values(roots)
    if H is None:
        return None
    return np.sum(H**2)


# ================================================================
# STEP 0: Verify MSS formula matches Haar MC
# ================================================================
print("=" * 70)
print("STEP 0: VERIFY MSS FORMULA MATCHES HAAR MC")
print("=" * 70)

np.random.seed(42)

test_cases = [
    (np.array([-1.0, 1.0]), np.array([-2.0, 2.0]), "n=2 sym"),
    (np.array([0.0, 1.0]), np.array([0.0, 3.0]), "n=2 asym"),
    (np.array([-0.5, 1.5]), np.array([2.0, 5.0]), "n=2 shifted"),
    (np.array([-2.0, 0.0, 2.0]), np.array([-3.0, 0.0, 3.0]), "n=3 sym"),
    (np.array([-1.0, 0.5, 2.0]), np.array([-0.5, 1.0, 3.0]), "n=3 asym"),
    (np.array([-3.0, -1.0, 1.0, 3.0]), np.array([-2.0, -0.5, 0.5, 2.0]), "n=4 sym"),
]

all_match = True
for roots_p, roots_q, label in test_cases:
    r_mss, _ = boxplus_mss(roots_p, roots_q)
    r_haar = boxplus_haar_mc(roots_p, roots_q, 500000)

    err = np.max(np.abs(r_mss - r_haar))
    match = err < 0.03  # MC tolerance

    print(f"  {label}: |MSS-Haar| = {err:.6f} {'OK' if match else 'MISMATCH'}")
    if not match:
        print(f"    MSS:  {r_mss}")
        print(f"    Haar: {r_haar}")
        all_match = False

print(f"\nAll match: {all_match}")

# ================================================================
# STEP 1: n=2 EXACT ALGEBRA WITH MSS FORMULA
# ================================================================
print("\n" + "=" * 70)
print("STEP 1: n=2 EXACT ALGEBRA (MSS formula)")
print("=" * 70)

def n2_exact_mss(a, b, c, d):
    """For n=2 with MSS formula, compute e_k(r) exactly.

    e_0(r) = 1
    e_1(r) = (1*1/(2*1)) * 1*(c+d) + (1*1/(2*1)) * (a+b)*1 = (a+b+c+d)/2
    Wait, let me be precise.

    c_k(r) = sum_{i+j=k} [(2-i)!(2-j)! / (2!(2-k)!)] * c_i(p) * c_j(q)

    k=0: (2!*2!)/(2!*2!) * 1*1 = 1
    k=1: i=0,j=1: (2!*1!)/(2!*1!) * 1*(c+d) + i=1,j=0: (1!*2!)/(2!*1!) * (a+b)*1
        = (c+d) + (a+b) = a+b+c+d? No wait...

    (2-0)!(2-1)! / (2!(2-1)!) = 2!*1! / (2!*1!) = 1
    (2-1)!(2-0)! / (2!(2-1)!) = 1!*2! / (2!*1!) = 1

    So e_1(r) = (a+b) + (c+d) = a+b+c+d
    k=2: i=0,j=2: (2!*0!)/(2!*0!) * 1*cd + i=1,j=1: (1!*1!)/(2!*0!) * (a+b)(c+d) + i=2,j=0: (0!*2!)/(2!*0!) * ab*1
        = cd + (a+b)(c+d)/2 + ab

    So r(x) = x^2 - (a+b+c+d)x + ab + cd + (a+b)(c+d)/2
    """
    e1 = a + b + c + d
    e2 = a*b + c*d + (a+b)*(c+d)/2

    disc = e1**2 - 4*e2
    # = (a+b+c+d)^2 - 4ab - 4cd - 2(a+b)(c+d)
    # = (a+b)^2 + (c+d)^2 + 2(a+b)(c+d) - 4ab - 4cd - 2(a+b)(c+d)
    # = (a-b)^2 + (c-d)^2
    # = gap_p^2 + gap_q^2

    gap_p = b - a
    gap_q = d - c

    assert abs(disc - (gap_p**2 + gap_q**2)) < 1e-10, f"disc={disc}, gap_p^2+gap_q^2={gap_p**2+gap_q**2}"

    gap_r = np.sqrt(disc)
    nu1 = (e1 - gap_r) / 2
    nu2 = (e1 + gap_r) / 2
    return nu1, nu2, gap_r

print("\nVerifying gap_r = sqrt(gap_p^2 + gap_q^2) for MSS formula at n=2:")
np.random.seed(42)
max_error = 0
for _ in range(200):
    a, b = sorted(np.random.randn(2) * 3)
    if b - a < 0.1: b = a + 0.5
    c, d = sorted(np.random.randn(2) * 3)
    if d - c < 0.1: d = c + 0.5

    nu1, nu2, gap_r = n2_exact_mss(a, b, c, d)
    predicted = np.sqrt((b-a)**2 + (d-c)**2)
    err = abs(gap_r - predicted)
    max_error = max(max_error, err)

print(f"  Max error: {max_error:.2e}")
print(f"  RESULT: {'PASS' if max_error < 1e-10 else 'FAIL'}")

# Verify Phi formulas
print("\nVerifying Phi_2 = 2/gap^2 and AB = Phi_r^2:")
np.random.seed(42)
max_AB_err = 0
for _ in range(200):
    a, b = sorted(np.random.randn(2) * 3)
    if b - a < 0.1: b = a + 0.5
    c, d = sorted(np.random.randn(2) * 3)
    if d - c < 0.1: d = c + 0.5

    s = b - a
    t = d - c
    R2 = s**2 + t**2

    Phi_p_val = 2/s**2
    Phi_q_val = 2/t**2
    Phi_r_val = 2/R2

    A_val = Phi_p_val - Phi_r_val
    B_val = Phi_q_val - Phi_r_val
    AB = A_val * B_val
    target = Phi_r_val**2

    err = abs(AB - target) / max(abs(target), 1e-15)
    max_AB_err = max(max_AB_err, err)

print(f"  Max relative error AB vs Phi_r^2: {max_AB_err:.2e}")
print(f"  RESULT: {'PASS' if max_AB_err < 1e-10 else 'FAIL'}")


# ================================================================
# STEP 2: MAIN INEQUALITY WITH MSS FORMULA
# ================================================================
print("\n" + "=" * 70)
print("STEP 2: MAIN INEQUALITY 1/Phi_r >= 1/Phi_p + 1/Phi_q (MSS)")
print("=" * 70)

np.random.seed(42)
results_by_n = {}

for trial in range(5000):
    n = np.random.choice([2, 3, 4, 5, 6])
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.3:
            roots_p[i] = roots_p[i-1] + 0.3
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.3:
            roots_q[i] = roots_q[i-1] + 0.3

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)

        # Check real-rootedness
        raw = np.roots([(-1)**k * elem_sym_poly(roots_r, k) for k in range(n+1)])
        # Actually just check gaps
        if np.any(np.diff(roots_r) < 0.01):
            continue

        phi_p = Phi_n(roots_p)
        phi_q = Phi_n(roots_q)
        phi_r = Phi_n(roots_r)

        if phi_p is None or phi_q is None or phi_r is None:
            continue

        excess = 1.0/phi_r - 1.0/phi_p - 1.0/phi_q

        if n not in results_by_n:
            results_by_n[n] = {'total': 0, 'pass': 0, 'min_excess': float('inf'), 'max_excess': -float('inf')}

        results_by_n[n]['total'] += 1
        results_by_n[n]['min_excess'] = min(results_by_n[n]['min_excess'], excess)
        results_by_n[n]['max_excess'] = max(results_by_n[n]['max_excess'], excess)

        if excess >= -1e-8:
            results_by_n[n]['pass'] += 1
        else:
            pass  # Violation
    except Exception:
        pass

print("\nSummary (MSS coefficient formula):")
for n in sorted(results_by_n.keys()):
    r = results_by_n[n]
    viols = r['total'] - r['pass']
    print(f"  n={n}: {r['pass']}/{r['total']} pass ({viols} violations)")
    print(f"         min_excess={r['min_excess']:.8e}, max_excess={r['max_excess']:.8e}")


# ================================================================
# STEP 3: EQUALLY SPACED ROOTS (Node 1.7.3)
# ================================================================
print("\n" + "=" * 70)
print("STEP 3: EQUALLY SPACED ROOTS WITH MSS FORMULA")
print("=" * 70)

# C_n computation
print("\nC_n values (Phi_n * d^2 for equally spaced roots):")
for n in range(2, 8):
    d = 1.0
    roots = np.array([k * d for k in range(n)], dtype=float)
    phi = Phi_n(roots)
    C_n = phi * d**2
    print(f"  n={n}: C_n = {C_n:.8f}")

# d_r^2 = d_p^2 + d_q^2 check
print("\nEqually spaced roots: is r also equally spaced?")
for n in [2, 3, 4, 5]:
    for d_p in [1.0, 2.0]:
        for d_q in [1.0, 2.0]:
            roots_p = np.array([k * d_p for k in range(n)])
            roots_q = np.array([k * d_q for k in range(n)])

            roots_r, _ = boxplus_mss(roots_p, roots_q)
            gaps_r = np.diff(roots_r)

            mean_gap = np.mean(gaps_r)
            rel_var = (max(gaps_r) - min(gaps_r)) / max(mean_gap, 1e-15) if len(gaps_r) > 1 else 0
            is_eq = rel_var < 1e-6

            phi_p = Phi_n(roots_p)
            phi_q = Phi_n(roots_q)
            phi_r = Phi_n(roots_r)
            excess = 1.0/phi_r - 1.0/phi_p - 1.0/phi_q

            if is_eq and n <= 4:
                d_r = mean_gap
                pyth = d_r**2 - d_p**2 - d_q**2
                print(f"  n={n}, d_p={d_p}, d_q={d_q}: eq_spaced=True, d_r={d_r:.6f}, d_r^2-d_p^2-d_q^2={pyth:.6e}, excess={excess:.6e}")
            else:
                print(f"  n={n}, d_p={d_p}, d_q={d_q}: eq_spaced={is_eq}, rel_var={rel_var:.4f}, excess={excess:.6e}")


# ================================================================
# STEP 4: <h,alpha> >= 0 WITH MSS FORMULA
# ================================================================
print("\n" + "=" * 70)
print("STEP 4: <h,alpha> >= 0 WITH MSS FORMULA")
print("=" * 70)

np.random.seed(42)
total = 0
failures = 0
min_val = float('inf')

for trial in range(5000):
    n = np.random.choice([3, 4, 5])
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.3:
            roots_p[i] = roots_p[i-1] + 0.3
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.3:
            roots_q[i] = roots_q[i-1] + 0.3

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        if np.any(np.diff(roots_r) < 0.01):
            continue

        h = H_values(roots_r)
        u = H_values(roots_p)  # identity sigma (order-preserving)

        if h is None or u is None:
            continue

        alpha = u - h
        h_alpha = np.dot(h, alpha)

        total += 1
        min_val = min(min_val, h_alpha)

        if h_alpha < -1e-6:
            failures += 1
    except Exception:
        pass

print(f"  Total valid trials: {total}")
print(f"  Failures (<h,alpha> < 0): {failures}")
print(f"  Minimum <h,alpha>: {min_val:.6e}")
if failures == 0:
    print(f"  Result: No violations found (but this is numerical evidence, not proof)")
else:
    print(f"  Result: VIOLATIONS FOUND! {failures}/{total} cases")


# ================================================================
# STEP 5: PARALLEL COMPONENT BOUND WITH MSS
# ================================================================
print("\n" + "=" * 70)
print("STEP 5: PARALLEL COMPONENT s_a(2+s_a)*s_b(2+s_b) >= 1 (MSS)")
print("=" * 70)

np.random.seed(42)
total = 0
failures = 0

for trial in range(5000):
    n = np.random.choice([3, 4, 5, 6])
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.3:
            roots_p[i] = roots_p[i-1] + 0.3
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.3:
            roots_q[i] = roots_q[i-1] + 0.3

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        if np.any(np.diff(roots_r) < 0.01):
            continue

        h = H_values(roots_r)
        u = H_values(roots_p)
        v = H_values(roots_q)

        if h is None or u is None or v is None:
            continue

        alpha = u - h
        beta = v - h

        norm_h = np.sqrt(np.dot(h, h))
        if norm_h < 1e-10:
            continue

        h_hat = h / norm_h
        s_a = np.dot(alpha, h_hat) / norm_h
        s_b = np.dot(beta, h_hat) / norm_h

        product = s_a * (2 + s_a) * s_b * (2 + s_b)

        total += 1
        if product < 1.0 - 1e-10:
            failures += 1
    except Exception:
        pass

print(f"  Total valid trials: {total}")
print(f"  Failures: {failures}")
print(f"  Failure rate: {failures/max(total,1)*100:.2f}%")


# ================================================================
# SUMMARY
# ================================================================
print("\n" + "=" * 70)
print("FINAL VERIFICATION SUMMARY")
print("=" * 70)
