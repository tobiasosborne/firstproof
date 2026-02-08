#!/usr/bin/env python3
"""
Investigate the subordination pairing sigma, tau more carefully.
The main script showed non-bijective sigma/tau (e.g., sigma=[0,0,2]),
which should not happen if sigma is a bijection. This may indicate
a bug in the pairing algorithm.

We verify the pairing by checking the consistency condition:
If sigma is correct, then alpha = u - h should satisfy
A = 2<h,alpha> + ||alpha||^2 = Phi(p) - Phi(r)
"""

import numpy as np
from math import factorial

np.random.seed(42)

def roots_to_monic_coeffs(roots):
    return np.poly(roots)

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

def Phi_n(roots):
    H = H_values(roots)
    return np.sum(H**2)

def generate_random_poly(n, spread=5.0):
    roots = np.sort(np.random.uniform(-spread, spread, n))
    for i in range(1, n):
        if roots[i] - roots[i-1] < 0.01:
            roots[i] = roots[i-1] + 0.01 + np.random.uniform(0, 0.1)
    return roots

def find_pairing_by_consistency(p_roots, q_roots, r_roots, n):
    """
    Since alpha = u - h and A = ||u||^2 - ||h||^2 = Phi(p) - Phi(r),
    the pairing sigma is irrelevant for the inequality verification!

    Phi(p) = sum of all H_p(lambda_i)^2 regardless of ordering.
    Phi(r) = sum of all H_r(nu_k)^2 regardless of ordering.

    So A = Phi(p) - Phi(r) doesn't depend on sigma at all.
    Similarly B = Phi(q) - Phi(r) doesn't depend on tau.

    The pairing only matters for understanding the STRUCTURE of alpha, beta.
    But for the main inequality verification, we don't need the pairing.

    However, for <h, alpha> we DO need the pairing.
    Let's try all permutations for small n and find the correct one.
    """
    from itertools import permutations

    H_p = H_values(p_roots)
    H_r = H_values(r_roots)

    # For each permutation sigma, compute A_check = 2<h, u_sigma - h> + ||u_sigma - h||^2
    # This should equal Phi(p) - Phi(r) for ANY bijection sigma.
    # Because ||u_sigma||^2 = sum H_p(lambda_{sigma(k)})^2 = sum H_p(lambda_i)^2 = Phi(p)
    # regardless of sigma (since sigma is a bijection).

    A_target = Phi_n(p_roots) - Phi_n(r_roots)

    best_sigma = None
    best_err = np.inf

    for perm in permutations(range(n)):
        sigma = np.array(perm)
        u = H_p[sigma]
        h = H_r
        alpha = u - h
        A_check = 2 * np.dot(h, alpha) + np.dot(alpha, alpha)
        err = abs(A_check - A_target)
        if err < best_err:
            best_err = err
            best_sigma = sigma

    return best_sigma


def find_pairing_continuation(p_roots, r_roots, n):
    """
    Find the subordination pairing by numerical continuation.
    Start from z = large real number, where omega_1(z) ~ z.
    Track omega_1(z) as z decreases towards the roots of r.

    The key equation: G_r(z) = G_p(omega_1(z))
    i.e., r'(z) * p(omega_1(z)) = p'(omega_1(z)) * r(z)

    For real z between roots of r, this has a unique real solution
    omega_1(z) between the corresponding roots of p (by the interlacing property).
    """
    p_poly = np.poly(p_roots)
    r_poly = np.poly(r_roots)
    dp = np.polyder(p_poly)
    dr = np.polyder(r_poly)

    # Start from z = max root + 100
    z_start = max(np.max(p_roots), np.max(r_roots)) + 100.0

    # At z_start, omega_1(z) ~ z, so omega_1(z_start) ~ z_start
    # Solve the polynomial equation to find omega_1(z_start) near z_start

    def solve_omega1(z):
        """Solve r'(z)*p(w) - p'(w)*r(z) = 0 for w."""
        rp_z = np.polyval(dr, z)
        r_z = np.polyval(r_poly, z)
        # rp_z * p(w) - p'(w) * r_z = 0
        # This is a polynomial in w:
        # rp_z * [w^n + a1*w^{n-1} + ... ] - r_z * [n*w^{n-1} + ...] = 0
        # Build combined polynomial
        combined = rp_z * p_poly.copy()
        padded_dp = np.zeros(n + 1)
        padded_dp[1:] = dp
        combined = combined - r_z * padded_dp
        roots = np.roots(combined)
        return roots

    # At z_start, find the root closest to z_start
    w_roots = solve_omega1(z_start)
    # Filter real roots
    real_mask = np.abs(np.imag(w_roots)) < 1e-6
    real_roots = np.real(w_roots[real_mask])
    # The one closest to z_start is our omega_1
    if len(real_roots) == 0:
        return np.arange(n)  # fallback

    idx = np.argmin(np.abs(real_roots - z_start))
    omega1_current = real_roots[idx]

    # Now trace: decrease z towards each root of r
    # Track which root of p omega_1 approaches
    sigma = np.zeros(n, dtype=int)

    for k in range(n - 1, -1, -1):
        # Trace omega_1 from current z to just above r_roots[k]
        nu_k = r_roots[k]
        # Evaluate at z = nu_k + small positive imaginary part
        z_eval = nu_k + 1e-8j

        w_roots = solve_omega1(z_eval)

        # Find the root closest to the last tracked omega_1
        # For the first step, use the one closest to nu_k + something
        if k == n - 1:
            # First root, omega_1 should be near the largest root of p
            target = p_roots[-1] + (nu_k - r_roots[-1]) if abs(nu_k - r_roots[-1]) < 10 else nu_k
        else:
            target = omega1_current

        # Among all solutions, find the one that is:
        # 1. Nearly real (small imaginary part)
        # 2. Close to a root of p
        best_p_idx = -1
        best_dist = np.inf

        for w in w_roots:
            # Check which root of p this is close to
            for m, lam in enumerate(p_roots):
                d = abs(w - lam)
                if d < best_dist:
                    best_dist = d
                    best_p_idx = m

        sigma[k] = best_p_idx

    return sigma


print("=" * 60)
print("INVESTIGATING SUBORDINATION PAIRING")
print("=" * 60)

print("\nKEY INSIGHT: The pairing sigma doesn't affect the main inequality!")
print("A = Phi(p) - Phi(r) and B = Phi(q) - Phi(r) are independent of sigma, tau.")
print("This is because ||u_sigma||^2 = sum H_p(lambda_{sigma(k)})^2 = Phi(p)")
print("for ANY bijection sigma (just reordering the sum).")
print()
print("The pairing DOES matter for decomposing A and B into <h,alpha> and ||alpha||^2.")
print("But the MAIN inequality (A*B >= ||h||^4) is pairing-independent.")

print("\n--- Checking consistency for n=3 ---")
n = 3
np.random.seed(12345)

for trial in range(10):
    p_roots = generate_random_poly(n)
    q_roots = generate_random_poly(n)
    r_coeffs = boxplus_n(roots_to_monic_coeffs(p_roots),
                         roots_to_monic_coeffs(q_roots), n)
    r_roots_raw = np.roots(r_coeffs)
    if np.max(np.abs(np.imag(r_roots_raw))) > 1e-8:
        continue
    r_roots = np.sort(np.real(r_roots_raw))
    if np.min(np.diff(r_roots)) < 1e-10:
        continue

    H_p = H_values(p_roots)
    H_r = H_values(r_roots)

    phi_p = np.sum(H_p**2)
    phi_r = np.sum(H_r**2)
    A_target = phi_p - phi_r

    # Check: for ANY bijection sigma, A_check should equal A_target
    from itertools import permutations
    all_consistent = True
    for perm in permutations(range(n)):
        sigma = np.array(perm)
        u = H_p[sigma]
        h = H_r
        alpha = u - h
        A_check = 2 * np.dot(h, alpha) + np.dot(alpha, alpha)
        if abs(A_check - A_target) > 1e-10:
            all_consistent = False
            break

    print(f"  Trial {trial}: A_target={A_target:.6e}, "
          f"all permutations give same A: {all_consistent}")

print("\nVERIFICATION: A = ||u||^2 - ||h||^2 = Phi(p) - Phi(r) for ANY bijection.")
print("This confirms the main inequality is pairing-independent.")

print("\n--- But <h,alpha> DOES depend on pairing ---")
print("Checking whether <h,alpha> >= 0 for the CORRECT (identity) pairing:")

np.random.seed(42)
n = 3
count_neg = 0
total = 0
for trial in range(1000):
    p_roots = generate_random_poly(n)
    q_roots = generate_random_poly(n)
    r_coeffs = boxplus_n(roots_to_monic_coeffs(p_roots),
                         roots_to_monic_coeffs(q_roots), n)
    r_roots_raw = np.roots(r_coeffs)
    if np.max(np.abs(np.imag(r_roots_raw))) > 1e-8:
        continue
    r_roots = np.sort(np.real(r_roots_raw))
    if np.min(np.diff(r_roots)) < 1e-10:
        continue

    H_p = H_values(p_roots)
    H_r = H_values(r_roots)

    total += 1

    # Try identity pairing sigma = id (order-preserving)
    u_id = H_p  # sigma = identity
    alpha_id = u_id - H_r
    h_alpha_id = np.dot(H_r, alpha_id)

    if h_alpha_id < -1e-10:
        count_neg += 1

print(f"\n  n=3, identity pairing sigma=id:")
print(f"    <h,alpha> < 0: {count_neg}/{total}")

# Now try the order-preserving pairing for all n
for n in [2, 3, 4, 5]:
    np.random.seed(42)
    count_neg_ha = 0
    count_neg_hb = 0
    total = 0
    min_ha = np.inf
    min_hb = np.inf

    for trial in range(1000):
        p_roots = generate_random_poly(n)
        q_roots = generate_random_poly(n)
        r_coeffs = boxplus_n(roots_to_monic_coeffs(p_roots),
                             roots_to_monic_coeffs(q_roots), n)
        r_roots_raw = np.roots(r_coeffs)
        if np.max(np.abs(np.imag(r_roots_raw))) > 1e-8:
            continue
        r_roots = np.sort(np.real(r_roots_raw))
        if np.min(np.diff(r_roots)) < 1e-10:
            continue

        H_p = H_values(p_roots)
        H_q = H_values(q_roots)
        H_r = H_values(r_roots)

        total += 1

        # Identity pairing (order-preserving)
        ha = np.dot(H_r, H_p - H_r)
        hb = np.dot(H_r, H_q - H_r)

        min_ha = min(min_ha, ha)
        min_hb = min(min_hb, hb)

        if ha < -1e-10:
            count_neg_ha += 1
        if hb < -1e-10:
            count_neg_hb += 1

    print(f"\n  n={n}, identity pairing:")
    print(f"    <h,alpha> < 0: {count_neg_ha}/{total}, min={min_ha:.6e}")
    print(f"    <h,beta>  < 0: {count_neg_hb}/{total}, min={min_hb:.6e}")

print("\n\nCRITICAL OBSERVATION:")
print("The subordination pairing sigma is NOT necessarily the identity.")
print("However, the main inequality A*B >= Phi(r)^2 does NOT depend on the pairing.")
print("So our numerical verification of the inequality is VALID regardless of pairing bugs.")
