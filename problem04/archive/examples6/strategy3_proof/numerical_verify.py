#!/usr/bin/env python3
"""
Numerical verification of the reformulated Fisher information inequality:

  (Phi_n(p) - Phi_n(r)) * (Phi_n(q) - Phi_n(r)) >= Phi_n(r)^2

where r = p boxplus_n q (MSS finite free additive convolution).

Also verifies:
(a) Vectors alpha, beta from subordination
(b) Inner products <h, alpha>, <h, beta>
(c) Whether A = 2<h,alpha> + ||alpha||^2 > 0 always (Phi_n(p) > Phi_n(r))
"""

import numpy as np
from math import factorial

np.random.seed(42)

# ============================================================
# Core Functions
# ============================================================

def roots_to_monic_coeffs(roots):
    """Return coefficients [1, a_1, ..., a_n] in descending order."""
    return np.poly(roots)

def boxplus_n(p_coeffs, q_coeffs, n):
    """
    Compute r = p boxplus_n q using MSS coefficient formula.
    c_k = sum_{i+j=k} [(n-i)!(n-j)! / (n!(n-k)!)] * a_i * b_j
    """
    a = p_coeffs
    b = q_coeffs
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
    """H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j)"""
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H

def Phi_n(roots):
    """Phi_n(p) = sum_i H_p(lambda_i)^2"""
    H = H_values(roots)
    return np.sum(H**2)

def generate_random_poly(n, spread=5.0):
    """Generate random monic polynomial of degree n with n distinct real roots."""
    roots = np.sort(np.random.uniform(-spread, spread, n))
    for i in range(1, n):
        if roots[i] - roots[i-1] < 0.01:
            roots[i] = roots[i-1] + 0.01 + np.random.uniform(0, 0.1)
    return roots

def find_subordination_pairing(p_roots, q_roots, r_roots, n):
    """
    Find the bijections sigma, tau such that
    omega_1(nu_k) = lambda_{sigma(k)}, omega_2(nu_k) = mu_{tau(k)}.

    Uses the equation: at z = nu_k (root of r), solve
    p'(w) / p(w) = r'(nu_k) / r(nu_k) ... but r(nu_k) = 0.

    Instead, from F(z,w) = r'(z)*p(w) - p'(w)*r(z) = 0:
    At z = nu_k: r'(nu_k)*p(w) = 0, so w is a root of p.
    omega_1(nu_k) is a root of p.

    To find the pairing: use numerical continuation.
    Evaluate omega_1 at z = nu_k + eps*i (small imaginary perturbation)
    by solving G_p(w) = G_r(z) as a polynomial equation.
    """
    p_poly = np.poly(p_roots)
    r_poly = np.poly(r_roots)
    dp_coeffs = np.polyder(p_poly)
    dr_coeffs = np.polyder(r_poly)

    sigma = np.zeros(n, dtype=int)
    for k in range(n):
        nu_k = r_roots[k]
        z_test = nu_k + 1e-6j

        # G_r(z_test) = r'(z_test) / (n * r(z_test))
        target_G = np.polyval(dr_coeffs, z_test) / (n * np.polyval(r_poly, z_test))

        # Solve: p'(w) - target_G * n * p(w) = 0
        # p'(w) has degree n-1, p(w) has degree n
        # So target_G * n * p(w) has degree n, p'(w) has degree n-1
        # Combined polynomial has degree n

        # Build the combined polynomial coefficients
        combined = -target_G * n * p_poly.copy()
        # Add p'(w) coefficients (need to pad to length n+1)
        padded_dp = np.zeros(n + 1, dtype=complex)
        padded_dp[1:] = dp_coeffs  # dp has n coefficients, pad with 0 at front
        combined = combined + padded_dp

        w_solutions = np.roots(combined)

        # Find solution closest to a root of p
        best_idx = -1
        best_dist = np.inf
        for w in w_solutions:
            if abs(w.imag) < 0.5:
                for m, lam in enumerate(p_roots):
                    dist = abs(w - lam)
                    if dist < best_dist:
                        best_dist = dist
                        best_idx = m
        sigma[k] = best_idx if best_idx >= 0 else k

    # Similarly for tau
    q_poly = np.poly(q_roots)
    dq_coeffs = np.polyder(q_poly)

    tau = np.zeros(n, dtype=int)
    for k in range(n):
        nu_k = r_roots[k]
        z_test = nu_k + 1e-6j
        target_G = np.polyval(dr_coeffs, z_test) / (n * np.polyval(r_poly, z_test))

        combined = -target_G * n * q_poly.copy()
        padded_dq = np.zeros(n + 1, dtype=complex)
        padded_dq[1:] = dq_coeffs
        combined = combined + padded_dq

        w_solutions = np.roots(combined)

        best_idx = -1
        best_dist = np.inf
        for w in w_solutions:
            if abs(w.imag) < 0.5:
                for m, mu in enumerate(q_roots):
                    dist = abs(w - mu)
                    if dist < best_dist:
                        best_dist = dist
                        best_idx = m
        tau[k] = best_idx if best_idx >= 0 else k

    return sigma, tau

def compute_vectors(p_roots, q_roots, r_roots, sigma, tau):
    """
    Compute u, v, h, alpha, beta vectors.
    u_k = H_p(lambda_{sigma(k)}), v_k = H_q(mu_{tau(k)}), h_k = H_r(nu_k).
    alpha = u - h, beta = v - h.
    """
    H_p = H_values(p_roots)
    H_q = H_values(q_roots)
    H_r = H_values(r_roots)
    u = H_p[sigma]
    v = H_q[tau]
    h = H_r
    alpha = u - h
    beta = v - h
    return u, v, h, alpha, beta

# ============================================================
# Verification Functions
# ============================================================

def verify_superadditivity(n, num_trials=1000):
    """Directly verify 1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q)"""
    print(f"\n{'='*60}")
    print(f"DIRECT SUPERADDITIVITY CHECK FOR n={n}")
    print(f"{'='*60}")

    violations = 0
    min_excess = np.inf
    valid = 0

    for _ in range(num_trials):
        p_roots = generate_random_poly(n)
        q_roots = generate_random_poly(n)
        r_coeffs = boxplus_n(roots_to_monic_coeffs(p_roots),
                             roots_to_monic_coeffs(q_roots), n)
        r_roots_raw = np.roots(r_coeffs)
        if np.max(np.abs(np.imag(r_roots_raw))) > 1e-8:
            continue
        r_roots = np.sort(np.real(r_roots_raw))
        if len(r_roots) < 2 or np.min(np.diff(r_roots)) < 1e-10:
            continue

        phi_p, phi_q, phi_r = Phi_n(p_roots), Phi_n(q_roots), Phi_n(r_roots)
        if phi_p <= 0 or phi_q <= 0 or phi_r <= 0:
            continue

        valid += 1
        excess = 1.0/phi_r - 1.0/phi_p - 1.0/phi_q
        min_excess = min(min_excess, excess)
        if excess < -1e-10:
            violations += 1
            print(f"  *** VIOLATION: excess = {excess:.8e}")

    print(f"  Valid trials: {valid}")
    print(f"  Violations: {violations}")
    print(f"  Min excess: {min_excess:.12e}")
    print(f"  Superadditivity {'HOLDS' if violations == 0 else 'FAILS'}")
    return violations

def verify_reformulated(n, num_trials=1000, verbose=False):
    """Verify (Phi(p)-Phi(r))*(Phi(q)-Phi(r)) >= Phi(r)^2"""
    print(f"\n{'='*60}")
    print(f"REFORMULATED INEQUALITY FOR n={n}, {num_trials} trials")
    print(f"{'='*60}")

    violations = 0
    valid = 0
    ratios = []
    h_alpha_list = []
    h_beta_list = []
    A_list = []
    B_list = []

    for trial in range(num_trials):
        p_roots = generate_random_poly(n)
        q_roots = generate_random_poly(n)
        r_coeffs = boxplus_n(roots_to_monic_coeffs(p_roots),
                             roots_to_monic_coeffs(q_roots), n)
        r_roots_raw = np.roots(r_coeffs)
        if np.max(np.abs(np.imag(r_roots_raw))) > 1e-8:
            continue
        r_roots = np.sort(np.real(r_roots_raw))
        if len(r_roots) < 2 or np.min(np.diff(r_roots)) < 1e-10:
            continue

        phi_p, phi_q, phi_r = Phi_n(p_roots), Phi_n(q_roots), Phi_n(r_roots)
        A = phi_p - phi_r
        B = phi_q - phi_r
        LHS = A * B
        RHS = phi_r**2
        ratio = LHS / RHS if RHS > 0 else np.inf

        valid += 1
        ratios.append(ratio)
        A_list.append(A)
        B_list.append(B)

        if LHS < RHS - 1e-10 * max(abs(LHS), abs(RHS), 1.0):
            violations += 1
            print(f"  *** VIOLATION trial {trial}: ratio={ratio:.10e}")
            print(f"      p={p_roots}, q={q_roots}, r={r_roots}")

        # Compute subordination
        sigma, tau = find_subordination_pairing(p_roots, q_roots, r_roots, n)
        u, v, h, alpha, beta = compute_vectors(p_roots, q_roots, r_roots, sigma, tau)
        h_alpha_list.append(np.dot(h, alpha))
        h_beta_list.append(np.dot(h, beta))

        if verbose and trial < 3:
            print(f"\n  Trial {trial}:")
            print(f"    sigma={sigma}, tau={tau}")
            print(f"    u={u}, v={v}, h={h}")
            print(f"    alpha={alpha}, beta={beta}")
            print(f"    <h,alpha>={np.dot(h,alpha):.6e}, <h,beta>={np.dot(h,beta):.6e}")
            print(f"    A={A:.6e}, B={B:.6e}")
            print(f"    A*B={LHS:.6e}, Phi(r)^2={RHS:.6e}, ratio={ratio:.6e}")

    ratios = np.array(ratios)
    h_alpha_arr = np.array(h_alpha_list)
    h_beta_arr = np.array(h_beta_list)
    A_arr = np.array(A_list)
    B_arr = np.array(B_list)

    print(f"\n--- Summary for n={n} ---")
    print(f"  Valid trials: {valid}/{num_trials}")
    print(f"  VIOLATIONS: {violations}")
    print(f"  Ratio LHS/RHS: min={np.min(ratios):.12e}, max={np.max(ratios):.6e}, "
          f"mean={np.mean(ratios):.6e}")
    print(f"  A (Phi(p)-Phi(r)): min={np.min(A_arr):.8e}, all>=0: {np.all(A_arr >= -1e-10)}")
    print(f"  B (Phi(q)-Phi(r)): min={np.min(B_arr):.8e}, all>=0: {np.all(B_arr >= -1e-10)}")
    print(f"  <h,alpha>: min={np.min(h_alpha_arr):.8e}, all>=0: {np.all(h_alpha_arr >= -1e-10)}")
    print(f"  <h,beta>:  min={np.min(h_beta_arr):.8e}, all>=0: {np.all(h_beta_arr >= -1e-10)}")
    neg_ha = np.sum(h_alpha_arr < -1e-10)
    neg_hb = np.sum(h_beta_arr < -1e-10)
    print(f"  <h,alpha> < 0 count: {neg_ha}/{valid}")
    print(f"  <h,beta>  < 0 count: {neg_hb}/{valid}")

    if n == 2:
        print(f"\n  n=2 EQUALITY TEST:")
        print(f"    max |ratio - 1| = {np.max(np.abs(ratios - 1)):.6e}")
        print(f"    Equality {'CONFIRMED' if np.max(np.abs(ratios - 1)) < 1e-6 else 'NOT confirmed'}")

    return ratios, violations

# ============================================================
# Detailed n=2 exact analysis
# ============================================================

def n2_exact():
    """
    For n=2, p(x) = (x-a)(x-b), q(x) = (x-c)(x-d).
    Phi_2(p) = 2/(b-a)^2.
    boxplus_2: r(x) = x^2 + (a+b+c+d)/2 * x + (a+b)(c+d)/2 ... no.
    Let's compute via the formula.
    """
    print(f"\n{'='*60}")
    print("EXACT ANALYSIS FOR n=2")
    print(f"{'='*60}")

    n = 2
    ratios = []
    for trial in range(1000):
        a, b = np.sort(np.random.uniform(-5, 5, 2))
        if b - a < 0.1: b = a + 0.1
        c, d = np.sort(np.random.uniform(-5, 5, 2))
        if d - c < 0.1: d = c + 0.1

        p_roots = np.array([a, b])
        q_roots = np.array([c, d])
        r_coeffs = boxplus_n(roots_to_monic_coeffs(p_roots),
                             roots_to_monic_coeffs(q_roots), n)
        r_roots = np.sort(np.real(np.roots(r_coeffs)))

        phi_p = 2.0 / (b - a)**2
        phi_q = 2.0 / (d - c)**2
        phi_r = Phi_n(r_roots)

        A, B = phi_p - phi_r, phi_q - phi_r
        ratio = A * B / phi_r**2 if phi_r > 0 else np.inf
        ratios.append(ratio)

        if trial < 5:
            print(f"  Trial {trial}: a={a:.3f}, b={b:.3f}, c={c:.3f}, d={d:.3f}")
            print(f"    Phi(p)={phi_p:.6e}, Phi(q)={phi_q:.6e}, Phi(r)={phi_r:.6e}")
            print(f"    A={A:.6e}, B={B:.6e}, ratio={ratio:.12e}")

    ratios = np.array(ratios)
    print(f"\n  1000 trials:")
    print(f"    min ratio = {np.min(ratios):.15e}")
    print(f"    max ratio = {np.max(ratios):.15e}")
    print(f"    max |ratio-1| = {np.max(np.abs(ratios - 1)):.6e}")
    print(f"    EQUALITY {'CONFIRMED' if np.max(np.abs(ratios - 1)) < 1e-6 else 'NOT confirmed'}")

# ============================================================
# Subordination detail for small n
# ============================================================

def subordination_detail(n, num_trials=200):
    """Detailed subordination analysis."""
    print(f"\n{'='*60}")
    print(f"SUBORDINATION DETAIL FOR n={n}, {num_trials} trials")
    print(f"{'='*60}")

    h_alpha_all = []
    h_beta_all = []
    alpha_norm_all = []
    beta_norm_all = []

    for trial in range(num_trials):
        p_roots = generate_random_poly(n)
        q_roots = generate_random_poly(n)
        r_coeffs = boxplus_n(roots_to_monic_coeffs(p_roots),
                             roots_to_monic_coeffs(q_roots), n)
        r_roots_raw = np.roots(r_coeffs)
        if np.max(np.abs(np.imag(r_roots_raw))) > 1e-8:
            continue
        r_roots = np.sort(np.real(r_roots_raw))
        if len(r_roots) < 2 or np.min(np.diff(r_roots)) < 1e-10:
            continue

        sigma, tau = find_subordination_pairing(p_roots, q_roots, r_roots, n)
        u, v, h, alpha, beta = compute_vectors(p_roots, q_roots, r_roots, sigma, tau)

        ha = np.dot(h, alpha)
        hb = np.dot(h, beta)
        h_alpha_all.append(ha)
        h_beta_all.append(hb)
        alpha_norm_all.append(np.linalg.norm(alpha))
        beta_norm_all.append(np.linalg.norm(beta))

        if trial < 5:
            print(f"\n  Trial {trial}:")
            print(f"    p_roots = {np.array2string(p_roots, precision=4)}")
            print(f"    q_roots = {np.array2string(q_roots, precision=4)}")
            print(f"    r_roots = {np.array2string(r_roots, precision=4)}")
            print(f"    sigma={sigma}, tau={tau}")
            print(f"    h     = {np.array2string(h, precision=6)}")
            print(f"    alpha = {np.array2string(alpha, precision=6)}")
            print(f"    beta  = {np.array2string(beta, precision=6)}")
            print(f"    <h,alpha>={ha:.8e}, <h,beta>={hb:.8e}")
            print(f"    ||alpha||={np.linalg.norm(alpha):.6e}, ||beta||={np.linalg.norm(beta):.6e}")
            A_check = 2*ha + np.dot(alpha, alpha)
            B_check = 2*hb + np.dot(beta, beta)
            print(f"    A=2<h,a>+|a|^2 = {A_check:.6e}")
            print(f"    B=2<h,b>+|b|^2 = {B_check:.6e}")

    h_alpha_all = np.array(h_alpha_all)
    h_beta_all = np.array(h_beta_all)

    print(f"\n  --- Summary for n={n} ---")
    print(f"  Valid: {len(h_alpha_all)}")
    print(f"  <h,alpha>: min={np.min(h_alpha_all):.8e}, max={np.max(h_alpha_all):.6e}, "
          f"mean={np.mean(h_alpha_all):.6e}")
    print(f"  <h,beta>:  min={np.min(h_beta_all):.8e}, max={np.max(h_beta_all):.6e}, "
          f"mean={np.mean(h_beta_all):.6e}")
    print(f"  <h,alpha> >= 0 always: {np.all(h_alpha_all >= -1e-10)}")
    print(f"  <h,beta>  >= 0 always: {np.all(h_beta_all >= -1e-10)}")
    neg_a = np.sum(h_alpha_all < -1e-10)
    neg_b = np.sum(h_beta_all < -1e-10)
    print(f"  <h,alpha> < 0: {neg_a}/{len(h_alpha_all)}")
    print(f"  <h,beta>  < 0: {neg_b}/{len(h_beta_all)}")

    return h_alpha_all, h_beta_all

# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    print("=" * 60)
    print("NUMERICAL VERIFICATION OF FISHER INFORMATION INEQUALITY")
    print("(Phi(p)-Phi(r))*(Phi(q)-Phi(r)) >= Phi(r)^2")
    print("where r = p boxplus_n q")
    print("=" * 60)

    # Part 0: Direct superadditivity
    for n in [2, 3, 4, 5]:
        verify_superadditivity(n, num_trials=1000)

    # Part 1: Reformulated inequality
    all_ratios = {}
    total_violations = 0
    for n in [2, 3, 4, 5]:
        ratios, viol = verify_reformulated(n, num_trials=1000, verbose=(n <= 3))
        all_ratios[n] = ratios
        total_violations += viol

    # Part 2: n=2 exact analysis
    n2_exact()

    # Part 3: Subordination detail
    for n in [2, 3]:
        subordination_detail(n, num_trials=200)

    # Final summary
    print(f"\n{'='*60}")
    print("FINAL SUMMARY")
    print(f"{'='*60}")
    print(f"Total violations of reformulated inequality: {total_violations}")
    for n in [2, 3, 4, 5]:
        r = all_ratios[n]
        print(f"  n={n}: min ratio={np.min(r):.12e}, trials={len(r)}")
