"""
VERIFIER-12: Adversarial verification of PROVER-13 claims.

Claims under test:
1. Finite free de Bruijn identity: dS/dt = Phi_n(p)
2. Entropy power inequality: N(p ⊞ q) >= N(p) + N(q)
3. Monotone gap: F(t) - G(t) is decreasing
4. Gaussian splitting: (p ⊞ G_{t/2}) ⊞ (q ⊞ G_{t/2}) = (p ⊞ q) ⊞ G_t
"""

import numpy as np
from math import factorial, comb
import warnings
import sys
warnings.filterwarnings('ignore')

# ============================================================
# Core infrastructure (from PROVER-13, verified independently)
# ============================================================

def finite_free_convolution(p_coeffs, q_coeffs):
    """MSS finite free convolution."""
    n = len(p_coeffs) - 1
    assert len(q_coeffs) - 1 == n
    c = np.zeros(n + 1)
    for k in range(n + 1):
        s = 0.0
        for i in range(k + 1):
            j = k - i
            if i <= n and j <= n:
                coeff = (factorial(n - i) * factorial(n - j)) / (factorial(n) * factorial(n - k))
                s += coeff * p_coeffs[i] * q_coeffs[j]
        c[k] = s
    return c

def coeffs_from_roots(roots):
    coeffs = np.polynomial.polynomial.polyfromroots(roots)[::-1]
    return coeffs / coeffs[0]

def roots_from_coeffs(coeffs):
    return np.sort(np.roots(coeffs)).real

def H_values(roots):
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                diff = roots[i] - roots[j]
                if abs(diff) < 1e-15:
                    return None
                H[i] += 1.0 / diff
    return H

def Phi_n(roots):
    H = H_values(roots)
    if H is None:
        return float('inf')
    return np.sum(H**2)

def hermite_prob_coeffs(n):
    coeffs = np.zeros(n + 1)
    for m in range(n // 2 + 1):
        k = 2 * m
        coeff = factorial(n) * ((-1)**m) / (factorial(m) * (2**m) * factorial(n - 2*m))
        coeffs[k] = coeff
    return coeffs

def gaussian_poly_coeffs(n, s):
    he = hermite_prob_coeffs(n)
    gs = np.zeros(n + 1)
    for k in range(n + 1):
        gs[k] = he[k] * s**(k / 2.0)
    return gs

def S_entropy(roots):
    n = len(roots)
    val = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            diff = abs(roots[j] - roots[i])
            if diff < 1e-15:
                return float('-inf')
            val += np.log(diff)
    return val

def N_power(roots):
    n = len(roots)
    m = n * (n - 1) // 2
    S = S_entropy(roots)
    if S == float('-inf'):
        return 0.0
    return np.exp(2 * S / m)


# ============================================================
# INDEPENDENT VERIFICATION: MSS convolution formula
# ============================================================

def verify_mss_convolution():
    """
    Verify the MSS convolution formula by checking against the known
    random matrix interpretation: eigenvalues of A + UBU*.
    We check a known case: G_s ⊞ G_t = G_{s+t}.
    """
    print("=" * 70)
    print("VERIFICATION 0: MSS Convolution Formula Correctness")
    print("=" * 70)

    # The MSS formula for finite free convolution is:
    # (p ⊞ q)(x) = sum_k c_k x^{n-k}
    # c_k = sum_{i+j=k} [C(n-i,n-k)*C(n-j,n-k)/C(n,n-k)] * a_i * b_j
    # which simplifies to the formula used.
    # Let's verify by an independent computation.

    # Cross-check: for n=2, p(x)=(x-a)(x-b), q(x)=(x-c)(x-d)
    # Coefficients: p = [1, -(a+b), ab], q = [1, -(c+d), cd]
    # MSS result should have:
    # c_0 = 1
    # c_1 = a_1 + b_1 = -(a+b+c+d)  [first moments add]
    # c_2 = ...

    for n in [2, 3, 4, 5]:
        # Verify with two Gaussians
        for s, t in [(1.0, 1.0), (0.5, 2.0), (1.0, 3.0)]:
            gs = gaussian_poly_coeffs(n, s)
            gt = gaussian_poly_coeffs(n, t)
            conv = finite_free_convolution(gs, gt)
            expected = gaussian_poly_coeffs(n, s + t)
            err = np.max(np.abs(conv - expected))
            if err > 1e-10:
                print(f"  FAIL: n={n}, G_{s} ⊞ G_{t} != G_{s+t}, err={err:.2e}")
                return False

    # Independent check: MSS coefficients match the Marcus-Spielman-Srivastava formula
    # For n=2: c_k = sum_{i+j=k} (2-i)!(2-j)! / (2!(2-k)!) * a_i * b_j
    n = 2
    a = np.array([1.0, -3.0, 2.0])  # (x-1)(x-2)
    b = np.array([1.0, -7.0, 10.0])  # (x-2)(x-5)

    c_manual = np.zeros(3)
    for k in range(3):
        for i in range(k+1):
            j = k - i
            if i <= 2 and j <= 2:
                coeff = (factorial(2-i) * factorial(2-j)) / (factorial(2) * factorial(2-k))
                c_manual[k] += coeff * a[i] * b[j]

    c_code = finite_free_convolution(a, b)
    err = np.max(np.abs(c_manual - c_code))
    print(f"  MSS manual vs code for n=2: err={err:.2e}")
    if err > 1e-14:
        print("  FAIL: MSS implementation differs from manual calculation")
        return False

    print("  MSS convolution formula: VERIFIED")
    return True


# ============================================================
# CLAIM 1: De Bruijn Identity
# ============================================================

def verify_debruijn():
    """
    PROVER-13 claims: dS/dt|_{t=0} = Phi_n(p)
    where S(p) = sum_{i<j} log|r_i - r_j| and p_t = p ⊞ G_t.

    Check:
    1. Sign of the identity
    2. Adversarial test cases
    3. Correctness of the finite difference approximation
    """
    print("\n" + "=" * 70)
    print("CLAIM 1: DE BRUIJN IDENTITY: dS/dt = Phi_n(p)")
    print("=" * 70)

    results = {
        'sign_errors': 0,
        'ratio_errors': 0,
        'total_tests': 0,
        'max_rel_error': 0.0,
        'adversarial_failures': [],
    }

    # --- Test 1a: Sign check ---
    # Classical de Bruijn: dH/dt = J(X) > 0 (entropy increases under heat flow)
    # Here S is like entropy (sum of log-gaps). It should INCREASE under heat flow
    # (roots spread apart). So dS/dt should be POSITIVE.
    # Phi_n = sum H_i^2 >= 0 always.
    # So dS/dt = +Phi_n is consistent with S increasing. SIGN IS CORRECT.
    print("\n--- Sign check ---")
    for n in [3, 4, 5]:
        roots_p = np.arange(1, n+1, dtype=float)
        coeffs_p = coeffs_from_roots(roots_p)
        S0 = S_entropy(roots_p)

        dt = 1e-6
        gt = gaussian_poly_coeffs(n, dt)
        ct = finite_free_convolution(coeffs_p, gt)
        rt = np.sort(roots_from_coeffs(ct).real)
        St = S_entropy(rt)

        dSdt = (St - S0) / dt
        phi = Phi_n(roots_p)

        print(f"  n={n}: S(0)={S0:.6f}, S(dt)={St:.6f}, dS/dt={dSdt:.6f}, Phi={phi:.6f}")
        print(f"    S increases: {St > S0}, dS/dt > 0: {dSdt > 0}, Phi > 0: {phi > 0}")

        if dSdt < 0:
            results['sign_errors'] += 1
            print("    *** SIGN ERROR: S should increase under heat flow ***")

    # --- Test 1b: Adversarial cases ---
    print("\n--- Adversarial: nearly-coincident roots ---")
    for n in [3, 4, 5]:
        for gap in [0.01, 0.001, 0.0001]:
            # Roots with one very small gap
            roots_p = np.arange(1, n+1, dtype=float)
            roots_p[1] = roots_p[0] + gap  # make gap between first two roots tiny
            roots_p = np.sort(roots_p)

            coeffs_p = coeffs_from_roots(roots_p)
            phi = Phi_n(roots_p)
            S0 = S_entropy(roots_p)

            # Use multiple step sizes for Richardson extrapolation
            dSdt_vals = []
            for dt in [1e-5, 5e-6, 2.5e-6]:
                gt = gaussian_poly_coeffs(n, dt)
                ct = finite_free_convolution(coeffs_p, gt)
                rt = np.sort(roots_from_coeffs(ct).real)
                St = S_entropy(rt)
                dSdt_vals.append((St - S0) / dt)

            # Richardson extrapolation
            rich = (4 * dSdt_vals[1] - dSdt_vals[0]) / 3
            rel_err = abs(rich - phi) / abs(phi) if phi > 0 else float('nan')
            results['total_tests'] += 1
            results['max_rel_error'] = max(results['max_rel_error'], rel_err)

            if rel_err > 0.01:
                results['adversarial_failures'].append(
                    f"n={n}, gap={gap}: rel_err={rel_err:.2e}")
                print(f"  n={n}, gap={gap}: Phi={phi:.4f}, dS/dt(Rich)={rich:.4f}, "
                      f"rel_err={rel_err:.2e} *** POTENTIAL ISSUE")
            else:
                print(f"  n={n}, gap={gap}: Phi={phi:.4f}, dS/dt(Rich)={rich:.4f}, "
                      f"rel_err={rel_err:.2e} OK")

    # --- Test 1c: Large n ---
    print("\n--- Adversarial: large n ---")
    for n in [10, 15, 20]:
        roots_p = np.arange(1, n+1, dtype=float) * 1.0
        coeffs_p = coeffs_from_roots(roots_p)
        phi = Phi_n(roots_p)
        S0 = S_entropy(roots_p)

        dSdt_vals = []
        for dt in [1e-4, 5e-5, 2.5e-5]:
            gt = gaussian_poly_coeffs(n, dt)
            ct = finite_free_convolution(coeffs_p, gt)
            rt = roots_from_coeffs(ct)
            # For large n, roots might become complex
            if np.any(np.abs(np.imag(rt)) > 1e-6):
                dSdt_vals.append(float('nan'))
                continue
            rt = np.sort(rt.real)
            St = S_entropy(rt)
            dSdt_vals.append((St - S0) / dt)

        if all(np.isnan(dSdt_vals)):
            print(f"  n={n}: roots became complex -- convolution unstable for large n")
            results['adversarial_failures'].append(
                f"n={n}: numerical instability (complex roots)")
        else:
            valid = [v for v in dSdt_vals if not np.isnan(v)]
            if len(valid) >= 2:
                rich = (4 * valid[1] - valid[0]) / 3
            else:
                rich = valid[0] if valid else float('nan')

            rel_err = abs(rich - phi) / abs(phi) if phi > 0 and not np.isnan(rich) else float('nan')
            results['total_tests'] += 1
            if not np.isnan(rel_err):
                results['max_rel_error'] = max(results['max_rel_error'], rel_err)

            status = "OK" if rel_err < 0.01 else "*** ISSUE ***"
            print(f"  n={n}: Phi={phi:.4f}, dS/dt~{rich:.4f}, rel_err={rel_err:.2e} {status}")

    # --- Test 1d: Verify dr_i/dt = H_i exactly (root dynamics) ---
    print("\n--- Root dynamics verification: dr_i/dt = H_i ---")
    # PROVER-13 claims dr_i/dt = H_i, but in de Bruijn Part 2, ratios showed
    # dr/dt / H is close to 1.0. Let me check more carefully.
    for n in [3, 4, 5, 6]:
        roots_p = np.arange(1, n+1, dtype=float) * 2.0
        coeffs_p = coeffs_from_roots(roots_p)
        H0 = H_values(roots_p)

        # Multiple dt values for consistency check
        ratios_per_dt = []
        for dt in [1e-5, 1e-6, 1e-7]:
            gt = gaussian_poly_coeffs(n, dt)
            ct = finite_free_convolution(coeffs_p, gt)
            rt = np.sort(roots_from_coeffs(ct).real)
            dr = (rt - roots_p) / dt
            ratios = dr / H0
            ratios_per_dt.append(np.mean(ratios))

        # Extrapolate to dt=0
        print(f"  n={n}: dr/dt / H at dt=1e-5: {ratios_per_dt[0]:.8f}, "
              f"dt=1e-6: {ratios_per_dt[1]:.8f}, dt=1e-7: {ratios_per_dt[2]:.8f}")
        # If ratio converges to 1.0, dr/dt = H_i exactly.
        # If to some other value c, dr/dt = c*H_i.

    # --- Test 1e: The algebraic identity sum_{i<j}(H_i-H_j)/(r_i-r_j) = sum_i H_i^2 ---
    print("\n--- Algebraic identity verification ---")
    print("  Testing: sum_{i<j} (H_i - H_j)/(r_i - r_j) = sum_i H_i^2 = Phi_n")
    for n in [3, 4, 5, 6, 7, 8]:
        for trial in range(5):
            if trial == 0:
                roots = np.arange(1, n+1, dtype=float)
            else:
                roots = np.sort(np.random.RandomState(trial*100+n).randn(n)*2 + np.arange(n)*3)
            if np.min(np.diff(roots)) < 0.2:
                continue

            H = H_values(roots)
            if H is None:
                continue
            phi = Phi_n(roots)

            lhs = 0.0
            for i in range(n):
                for j in range(i+1, n):
                    lhs += (H[i] - H[j]) / (roots[i] - roots[j])

            rel_err = abs(lhs - phi) / abs(phi) if phi > 0 else 0
            if rel_err > 1e-10:
                print(f"  n={n}, trial {trial}: LHS={lhs:.8f}, Phi={phi:.8f}, "
                      f"rel_err={rel_err:.2e} *** MISMATCH")
            elif trial == 0:
                print(f"  n={n}: sum (H_i-H_j)/(r_i-r_j) = {lhs:.8f}, Phi = {phi:.8f}, match: YES")

    print(f"\n  SUMMARY: {results['sign_errors']} sign errors, "
          f"{len(results['adversarial_failures'])} adversarial failures, "
          f"max rel error = {results['max_rel_error']:.2e}")
    return results


# ============================================================
# CLAIM 2: Entropy Power Inequality
# ============================================================

def verify_epi():
    """
    PROVER-13 claims: N(p ⊞ q) >= N(p) + N(q)
    where N(p) = exp(2*S(p)/m), m = n(n-1)/2, S = sum_{i<j} log|r_i - r_j|.

    Adversarial tests to find counterexamples.
    """
    print("\n" + "=" * 70)
    print("CLAIM 2: ENTROPY POWER INEQUALITY")
    print("=" * 70)

    total_violations = 0
    total_tests = 0
    min_ratios = {}

    # --- Test 2a: n=2 (simplest case) ---
    print("\n--- n=2 (simplest case) ---")
    n = 2
    violations_n2 = 0
    total_n2 = 0
    min_ratio_n2 = float('inf')
    np.random.seed(12345)

    for trial in range(5000):
        a = np.random.randn() * 3
        b = a + abs(np.random.randn()) * 0.5 + 0.01
        c = np.random.randn() * 3
        d = c + abs(np.random.randn()) * 0.5 + 0.01
        roots_p = np.sort([a, b])
        roots_q = np.sort([c, d])

        coeffs_p = coeffs_from_roots(roots_p)
        coeffs_q = coeffs_from_roots(roots_q)
        coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)
        roots_pq = roots_from_coeffs(coeffs_pq)

        if np.any(np.abs(np.imag(roots_pq)) > 1e-8):
            continue
        roots_pq = np.sort(roots_pq.real)
        if np.min(np.diff(roots_pq)) < 1e-12:
            continue

        N_p = N_power(roots_p)
        N_q = N_power(roots_q)
        N_pq = N_power(roots_pq)

        total_n2 += 1
        ratio = N_pq / (N_p + N_q)
        min_ratio_n2 = min(min_ratio_n2, ratio)

        if N_pq < N_p + N_q - 1e-10:
            violations_n2 += 1
            if violations_n2 <= 3:
                print(f"  VIOLATION: roots_p={roots_p}, roots_q={roots_q}")
                print(f"    N(p)={N_p:.8f}, N(q)={N_q:.8f}, N(p+q)={N_pq:.8f}, sum={N_p+N_q:.8f}")

    print(f"  n=2: {violations_n2}/{total_n2} violations, min ratio = {min_ratio_n2:.6f}")
    total_violations += violations_n2
    total_tests += total_n2
    min_ratios[2] = min_ratio_n2

    # For n=2: analytic check.
    # p(x) = (x-a)(x-b), q(x) = (x-c)(x-d)
    # S(p) = log|a-b|, m = 1, N(p) = exp(2*log|a-b|) = (a-b)^2
    # S(q) = log|c-d|, N(q) = (c-d)^2
    # For n=2 MSS: p ⊞ q has coefficients:
    #   c_0 = 1, c_1 = a_1 + b_1, c_2 = a_1*b_1/(n) + a_2 + b_2 [need to compute]
    # Actually let me just compute for specific n=2.

    print("\n--- n=2 analytic check ---")
    # For n=2, MSS convolution:
    # c_0 = 1
    # c_1 = (2-0)!(2-1)!/(2!(2-1)!) * a_0*b_1 + (2-1)!(2-0)!/(2!(2-1)!) * a_1*b_0
    #     = (2*1)/(2*1) * b_1 + (1*2)/(2*1) * a_1 = b_1 + a_1
    # c_2 = sum over i+j=2: (2-i)!(2-j)!/(2!*0!) * a_i*b_j
    #     = (2!*0!)/(2!*1) * a_0*b_2 + (1!*1!)/(2!*1) * a_1*b_1 + (0!*2!)/(2!*1) * a_2*b_0
    #     = b_2 + a_1*b_1/2 + a_2
    # So for p = [1, -(a+b), ab], q = [1, -(c+d), cd]:
    #   c_1 = -(a+b+c+d)
    #   c_2 = cd + (a+b)(c+d)/2 + ab
    # Roots of (x^2 + c_1*x + c_2):
    #   mu_{1,2} = (-c_1 +/- sqrt(c_1^2 - 4c_2)) / 2
    # Gap: mu_2 - mu_1 = sqrt(c_1^2 - 4c_2)
    # N(p⊞q) = (mu_2 - mu_1)^2 = c_1^2 - 4c_2

    a_val, b_val = 0.0, 2.0
    c_val, d_val = 0.0, 3.0

    N_p_analytic = (b_val - a_val)**2  # = 4
    N_q_analytic = (d_val - c_val)**2  # = 9

    c1 = -(a_val + b_val + c_val + d_val)
    c2 = c_val*d_val + (a_val+b_val)*(c_val+d_val)/2 + a_val*b_val
    gap_sq = c1**2 - 4*c2
    N_pq_analytic = gap_sq

    print(f"  p = (x-{a_val})(x-{b_val}), q = (x-{c_val})(x-{d_val})")
    print(f"  N(p) = {N_p_analytic}, N(q) = {N_q_analytic}")
    print(f"  N(p⊞q) = {N_pq_analytic}")
    print(f"  N(p⊞q) >= N(p) + N(q)? {N_pq_analytic} >= {N_p_analytic + N_q_analytic}? "
          f"{N_pq_analytic >= N_p_analytic + N_q_analytic}")

    # Verify with code
    cp = coeffs_from_roots([a_val, b_val])
    cq = coeffs_from_roots([c_val, d_val])
    cpq = finite_free_convolution(cp, cq)
    rpq = roots_from_coeffs(cpq)
    N_pq_code = N_power(np.sort(rpq.real))
    print(f"  Code N(p⊞q) = {N_pq_code:.6f}, analytic = {N_pq_analytic:.6f}")

    # --- Test 2b: ADVERSARIAL - very unequal spacings ---
    print("\n--- Adversarial: very unequal root spacings ---")
    np.random.seed(54321)
    for n in [3, 4, 5]:
        n_viol = 0
        n_tot = 0
        min_r = float('inf')
        for trial in range(2000):
            # p: tightly clustered roots, q: widely spread roots
            spread_p = np.random.exponential(0.1)
            spread_q = np.random.exponential(10.0)
            center_p = np.random.randn() * 5
            center_q = np.random.randn() * 5
            roots_p = np.sort(center_p + np.arange(n) * spread_p)
            roots_q = np.sort(center_q + np.arange(n) * spread_q)
            if np.min(np.diff(roots_p)) < 1e-8 or np.min(np.diff(roots_q)) < 1e-8:
                continue

            coeffs_p = coeffs_from_roots(roots_p)
            coeffs_q = coeffs_from_roots(roots_q)
            coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)
            roots_pq = roots_from_coeffs(coeffs_pq)
            if np.any(np.abs(np.imag(roots_pq)) > 1e-6):
                continue
            roots_pq = np.sort(roots_pq.real)
            if np.min(np.diff(roots_pq)) < 1e-12:
                continue

            N_p = N_power(roots_p)
            N_q = N_power(roots_q)
            N_pq = N_power(roots_pq)

            if N_p <= 0 or N_q <= 0 or N_pq <= 0:
                continue

            n_tot += 1
            ratio = N_pq / (N_p + N_q)
            min_r = min(min_r, ratio)
            if N_pq < N_p + N_q - 1e-8:
                n_viol += 1
                if n_viol <= 2:
                    print(f"  VIOLATION n={n}: spread_p={spread_p:.4f}, spread_q={spread_q:.4f}")
                    print(f"    N(p)={N_p:.6f}, N(q)={N_q:.6f}, N(p⊞q)={N_pq:.6f}")

        print(f"  n={n}: {n_viol}/{n_tot} violations, min ratio = {min_r:.6f}")
        total_violations += n_viol
        total_tests += n_tot
        min_ratios[n] = min_r

    # --- Test 2c: Near-degenerate polynomials ---
    print("\n--- Adversarial: near-degenerate (roots very close) ---")
    np.random.seed(99999)
    for n in [3, 4, 5]:
        n_viol = 0
        n_tot = 0
        min_r = float('inf')
        for trial in range(1000):
            # Make roots very close together for p
            base = np.random.randn() * 5
            eps = np.random.exponential(0.001)
            roots_p = base + np.arange(n) * eps
            roots_q = np.sort(np.random.randn(n) * 2 + np.arange(n) * 1.5)

            if np.min(np.diff(roots_p)) < 1e-10 or np.min(np.diff(roots_q)) < 0.1:
                continue

            coeffs_p = coeffs_from_roots(roots_p)
            coeffs_q = coeffs_from_roots(roots_q)
            coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)
            roots_pq = roots_from_coeffs(coeffs_pq)
            if np.any(np.abs(np.imag(roots_pq)) > 1e-6):
                continue
            roots_pq = np.sort(roots_pq.real)
            if np.min(np.diff(roots_pq)) < 1e-12:
                continue

            N_p = N_power(roots_p)
            N_q = N_power(roots_q)
            N_pq = N_power(roots_pq)

            if N_p <= 0 or N_q <= 0 or N_pq <= 0:
                continue

            n_tot += 1
            ratio = N_pq / (N_p + N_q)
            min_r = min(min_r, ratio)
            if N_pq < N_p + N_q - 1e-8:
                n_viol += 1

        print(f"  n={n}: {n_viol}/{n_tot} violations, min ratio = {min_r:.6f}")
        total_violations += n_viol
        total_tests += n_tot

    # --- Test 2d: Large n ---
    print("\n--- Adversarial: large n ---")
    np.random.seed(7777)
    for n in [10, 15, 20]:
        n_viol = 0
        n_tot = 0
        min_r = float('inf')
        for trial in range(100):
            roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n) * 2.0)
            roots_q = np.sort(np.random.randn(n) * 2 + np.arange(n) * 1.5)
            if np.min(np.diff(roots_p)) < 0.1 or np.min(np.diff(roots_q)) < 0.1:
                continue

            coeffs_p = coeffs_from_roots(roots_p)
            coeffs_q = coeffs_from_roots(roots_q)
            try:
                coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)
            except:
                continue
            roots_pq = roots_from_coeffs(coeffs_pq)
            if np.any(np.abs(np.imag(roots_pq)) > 1e-6):
                continue
            roots_pq = np.sort(roots_pq.real)
            if len(roots_pq) < n or np.min(np.diff(roots_pq)) < 1e-10:
                continue

            N_p = N_power(roots_p)
            N_q = N_power(roots_q)
            N_pq = N_power(roots_pq)

            if N_p <= 0 or N_q <= 0 or N_pq <= 0 or np.isinf(N_pq):
                continue

            n_tot += 1
            ratio = N_pq / (N_p + N_q)
            min_r = min(min_r, ratio)
            if N_pq < N_p + N_q - 1e-6:
                n_viol += 1
                if n_viol <= 2:
                    print(f"  VIOLATION n={n}: N(p)={N_p:.4f}, N(q)={N_q:.4f}, N(p⊞q)={N_pq:.4f}")

        print(f"  n={n}: {n_viol}/{n_tot} violations, min ratio = {min_r:.6f}" if n_tot > 0
              else f"  n={n}: no valid tests (numerical instability)")
        total_violations += n_viol
        total_tests += n_tot
        if n_tot > 0:
            min_ratios[n] = min_r

    # --- Test 2e: Check normalization m = n(n-1)/2 ---
    print("\n--- Normalization check: is m = n(n-1)/2 the right choice? ---")
    # For Gaussian G_s: S = (m/2)*log(s) + S_0, so N = exp(2S/m) = s * exp(2S_0/m)
    # This means N(G_s) = C_n * s, which is linear in s.
    # Then EPI for G_s ⊞ G_t: N(G_{s+t}) = C_n*(s+t) >= C_n*s + C_n*t = N(G_s) + N(G_t)
    # This holds with EQUALITY! Good -- Gaussians are the equality case.
    for n in [3, 4, 5, 6]:
        N_vals = []
        s_vals = [0.5, 1.0, 2.0, 5.0]
        for s in s_vals:
            roots_g = np.sort(roots_from_coeffs(gaussian_poly_coeffs(n, s)).real)
            N_vals.append(N_power(roots_g))
        # Check linearity: N(G_s) = C_n * s
        ratios = [N_vals[i] / s_vals[i] for i in range(len(s_vals))]
        C_n = np.mean(ratios)
        max_dev = max(abs(r/C_n - 1) for r in ratios)
        print(f"  n={n}: N(G_s)/s ratios = {[f'{r:.6f}' for r in ratios]}, max deviation = {max_dev:.2e}")
        if max_dev > 1e-8:
            print(f"    *** N(G_s) is NOT linear in s! Normalization may be wrong. ***")

        # Verify Gaussian equality: N(G_s ⊞ G_t) = N(G_s) + N(G_t)
        s, t = 1.0, 2.0
        gs = gaussian_poly_coeffs(n, s)
        gt = gaussian_poly_coeffs(n, t)
        gst = finite_free_convolution(gs, gt)
        roots_s = np.sort(roots_from_coeffs(gs).real)
        roots_t = np.sort(roots_from_coeffs(gt).real)
        roots_st = np.sort(roots_from_coeffs(gst).real)
        N_s = N_power(roots_s)
        N_t = N_power(roots_t)
        N_st = N_power(roots_st)
        print(f"    G equality: N(G_1⊞G_2) = {N_st:.6f}, N(G_1)+N(G_2) = {N_s+N_t:.6f}, "
              f"err = {abs(N_st - N_s - N_t):.2e}")

    print(f"\n  TOTAL VIOLATIONS: {total_violations}/{total_tests}")
    print(f"  Min ratios by n: {min_ratios}")
    return total_violations, total_tests, min_ratios


# ============================================================
# CLAIM 3: Monotone Gap
# ============================================================

def verify_monotone_gap():
    """
    PROVER-13 claims: F(t) - G(t) is decreasing where
    F(t) = 1/Phi((p⊞q)⊞G_{2t}), G(t) = 1/Phi(p⊞G_t) + 1/Phi(q⊞G_t).

    Critical checks:
    1. Does the gap actually decrease monotonically?
    2. Does gap -> 0 as t -> inf (from ABOVE)?
    3. Does this actually prove the conjecture?
    """
    print("\n" + "=" * 70)
    print("CLAIM 3: MONOTONE GAP")
    print("=" * 70)

    # --- Test 3a: Basic monotonicity ---
    print("\n--- Basic gap monotonicity test ---")
    np.random.seed(2024)
    gap_increase_count = 0
    total_pairs = 0

    for n in [3, 4, 5, 6]:
        for trial in range(50):
            roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n) * 1.5)
            roots_q = np.sort(np.random.randn(n) * 1.5 + np.arange(n) * 1.0)
            if np.min(np.diff(roots_p)) < 0.1 or np.min(np.diff(roots_q)) < 0.1:
                continue

            coeffs_p = coeffs_from_roots(roots_p)
            coeffs_q = coeffs_from_roots(roots_q)
            coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)

            t_vals = np.linspace(0.001, 5.0, 300)
            gaps = []
            valid_ts = []

            for t in t_vals:
                g2t = gaussian_poly_coeffs(n, 2*t)
                gt = gaussian_poly_coeffs(n, t)

                c_pq_2t = finite_free_convolution(coeffs_pq, g2t)
                c_p_t = finite_free_convolution(coeffs_p, gt)
                c_q_t = finite_free_convolution(coeffs_q, gt)

                r_pq = roots_from_coeffs(c_pq_2t)
                r_p = roots_from_coeffs(c_p_t)
                r_q = roots_from_coeffs(c_q_t)

                if (np.any(np.abs(np.imag(r_pq)) > 1e-8) or
                    np.any(np.abs(np.imag(r_p)) > 1e-8) or
                    np.any(np.abs(np.imag(r_q)) > 1e-8)):
                    continue

                r_pq = np.sort(r_pq.real)
                r_p = np.sort(r_p.real)
                r_q = np.sort(r_q.real)

                if (np.min(np.diff(r_pq)) < 1e-10 or
                    np.min(np.diff(r_p)) < 1e-10 or
                    np.min(np.diff(r_q)) < 1e-10):
                    continue

                F = 1.0 / Phi_n(r_pq)
                G = 1.0 / Phi_n(r_p) + 1.0 / Phi_n(r_q)
                gaps.append(F - G)
                valid_ts.append(t)

            if len(gaps) < 20:
                continue

            total_pairs += 1
            gaps = np.array(gaps)

            # Check monotonicity
            dgaps = np.diff(gaps)
            n_increase = np.sum(dgaps > 1e-8)

            if n_increase > 0:
                gap_increase_count += 1
                # Find the worst increase
                worst_idx = np.argmax(dgaps)
                print(f"  n={n}, trial {trial}: GAP INCREASES at t~{valid_ts[worst_idx]:.3f}, "
                      f"delta={dgaps[worst_idx]:.2e}, total increases: {n_increase}")
                print(f"    gap[0]={gaps[0]:.6f}, gap[-1]={gaps[-1]:.6f}")

    print(f"\n  Gap increases in {gap_increase_count}/{total_pairs} polynomial pairs")

    # --- Test 3b: CRITICAL LOGICAL CHECK ---
    print("\n--- CRITICAL: Logical structure of the argument ---")
    print("  PROVER-13's argument:")
    print("    (a) gap(t) is decreasing")
    print("    (b) gap(t) -> 0 as t -> inf")
    print("    (c) Therefore gap(0) >= 0")
    print()
    print("  This argument requires gap(t) -> 0 from ABOVE.")
    print("  If gap(t) -> 0 from BELOW (gap is negative for large t),")
    print("  then decreasing + limit 0 would mean gap(0) <= 0, WRONG direction!")
    print()
    print("  Let's check the SIGN of the gap for large t:")

    np.random.seed(42)
    for n in [3, 4, 5]:
        for trial in range(5):
            roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n) * 1.5)
            roots_q = np.sort(np.random.randn(n) * 1.5 + np.arange(n) * 1.0)
            if np.min(np.diff(roots_p)) < 0.2 or np.min(np.diff(roots_q)) < 0.2:
                continue

            coeffs_p = coeffs_from_roots(roots_p)
            coeffs_q = coeffs_from_roots(roots_q)
            coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)

            # Check gap at large t
            for t in [10.0, 50.0, 100.0]:
                g2t = gaussian_poly_coeffs(n, 2*t)
                gt = gaussian_poly_coeffs(n, t)

                c_pq = finite_free_convolution(coeffs_pq, g2t)
                c_p = finite_free_convolution(coeffs_p, gt)
                c_q = finite_free_convolution(coeffs_q, gt)

                r_pq = roots_from_coeffs(c_pq)
                r_p = roots_from_coeffs(c_p)
                r_q = roots_from_coeffs(c_q)

                if (np.any(np.abs(np.imag(r_pq)) > 1e-6) or
                    np.any(np.abs(np.imag(r_p)) > 1e-6) or
                    np.any(np.abs(np.imag(r_q)) > 1e-6)):
                    continue

                F = 1.0 / Phi_n(np.sort(r_pq.real))
                G = 1.0 / Phi_n(np.sort(r_p.real)) + 1.0 / Phi_n(np.sort(r_q.real))
                gap = F - G

                # Expected asymptotics: F ~ 8t/(n(n-1)), G ~ 8t/(n(n-1))
                F_asymp = 8*t / (n*(n-1))
                G_asymp = 8*t / (n*(n-1))

                if trial == 0:
                    print(f"  n={n}, t={t}: gap={gap:.8f}, F={F:.6f}, G={G:.6f}, "
                          f"F_asymp={F_asymp:.6f}")

            break

    # --- Test 3c: Asymptotic analysis ---
    print("\n--- Asymptotic correction terms ---")
    print("  As t -> inf:")
    print("  1/Phi(p ⊞ G_t) = 4t/(n(n-1)) + correction(p)/t + O(1/t^2)")
    print("  F(t) = 1/Phi((p⊞q) ⊞ G_{2t}) = 8t/(n(n-1)) + correction(p⊞q)/(2t) + ...")
    print("  G(t) = 1/Phi(p ⊞ G_t) + 1/Phi(q ⊞ G_t)")
    print("       = 8t/(n(n-1)) + [correction(p) + correction(q)]/t + ...")
    print("  gap(t) = [correction(p⊞q)/2 - correction(p) - correction(q)] / t + ...")
    print("  The sign of the correction determines whether gap -> 0 from above or below.")
    print()

    # Numerically extract the correction
    np.random.seed(42)
    for n in [3, 4, 5]:
        roots_p = np.arange(1, n+1, dtype=float)
        roots_q = np.arange(1, n+1, dtype=float) * 0.5
        coeffs_p = coeffs_from_roots(roots_p)
        coeffs_q = coeffs_from_roots(roots_q)
        coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)

        large_t_vals = [50.0, 100.0, 200.0, 500.0]
        for t in large_t_vals:
            g2t = gaussian_poly_coeffs(n, 2*t)
            gt = gaussian_poly_coeffs(n, t)

            c_pq = finite_free_convolution(coeffs_pq, g2t)
            c_p = finite_free_convolution(coeffs_p, gt)
            c_q = finite_free_convolution(coeffs_q, gt)

            r_pq = roots_from_coeffs(c_pq)
            r_p = roots_from_coeffs(c_p)
            r_q = roots_from_coeffs(c_q)

            if (np.any(np.abs(np.imag(r_pq)) > 1e-6) or
                np.any(np.abs(np.imag(r_p)) > 1e-6) or
                np.any(np.abs(np.imag(r_q)) > 1e-6)):
                continue

            F = 1.0 / Phi_n(np.sort(r_pq.real))
            G = 1.0 / Phi_n(np.sort(r_p.real)) + 1.0 / Phi_n(np.sort(r_q.real))
            gap = F - G
            gap_times_t = gap * t  # should converge to a constant

            if t == large_t_vals[0]:
                print(f"  n={n}: gap*t for large t:", end="")
            print(f" {gap_times_t:.6f}", end="")
        print(f"  (sign: {'POSITIVE' if gap_times_t > 0 else 'NEGATIVE'})")

    return gap_increase_count, total_pairs


# ============================================================
# CLAIM 4: Gaussian Splitting
# ============================================================

def verify_gaussian_splitting():
    """
    PROVER-13 claims: (p ⊞ G_{t/2}) ⊞ (q ⊞ G_{t/2}) = (p ⊞ q) ⊞ G_t

    Verify:
    1. Numerically
    2. Via cumulant additivity argument
    3. Check if this is exact or approximate
    """
    print("\n" + "=" * 70)
    print("CLAIM 4: GAUSSIAN SPLITTING")
    print("=" * 70)

    # --- Test 4a: Extensive numerical verification ---
    print("\n--- Numerical verification ---")
    np.random.seed(2025)
    max_err = 0.0
    total = 0

    for n in [2, 3, 4, 5, 6, 7, 8, 10]:
        n_total = 0
        n_max_err = 0.0
        for trial in range(50):
            roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n) * 1.5)
            roots_q = np.sort(np.random.randn(n) * 1.5 + np.arange(n) * 1.0)
            coeffs_p = coeffs_from_roots(roots_p)
            coeffs_q = coeffs_from_roots(roots_q)

            for t in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
                # LHS: (p ⊞ q) ⊞ G_t
                pq = finite_free_convolution(coeffs_p, coeffs_q)
                gt = gaussian_poly_coeffs(n, t)
                lhs = finite_free_convolution(pq, gt)

                # RHS: (p ⊞ G_{t/2}) ⊞ (q ⊞ G_{t/2})
                gth = gaussian_poly_coeffs(n, t / 2)
                p_h = finite_free_convolution(coeffs_p, gth)
                q_h = finite_free_convolution(coeffs_q, gth)
                rhs = finite_free_convolution(p_h, q_h)

                err = np.max(np.abs(lhs - rhs))
                n_max_err = max(n_max_err, err)
                n_total += 1
                total += 1
                max_err = max(max_err, err)

        print(f"  n={n}: max error over {n_total} tests = {n_max_err:.2e}")

    print(f"  Overall max error: {max_err:.2e}")

    # --- Test 4b: The cumulant argument ---
    print("\n--- Cumulant additivity argument ---")
    print("  PROVER-13 claims this follows from cumulant additivity.")
    print("  Let's verify that finite free cumulants ARE additive under MSS convolution.")

    # Compute finite free cumulants
    def finite_free_cumulants(coeffs):
        """Compute finite free cumulants: kappa_k = (-1)^k * a_k * n^{k-1} / C(n,k)"""
        n = len(coeffs) - 1
        kappas = [0.0]  # kappa_0 placeholder
        for k in range(1, n + 1):
            kappas.append((-1)**k * coeffs[k] * n**(k-1) / comb(n, k))
        return kappas

    np.random.seed(42)
    cumulant_additive = True
    for n in [3, 4, 5, 6]:
        roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n))
        roots_q = np.sort(np.random.randn(n) * 1.5 + np.arange(n) * 0.5)
        coeffs_p = coeffs_from_roots(roots_p)
        coeffs_q = coeffs_from_roots(roots_q)
        coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)

        kp = finite_free_cumulants(coeffs_p)
        kq = finite_free_cumulants(coeffs_q)
        kpq = finite_free_cumulants(coeffs_pq)

        max_err_k = 0.0
        for k in range(1, n+1):
            err = abs(kpq[k] - kp[k] - kq[k])
            max_err_k = max(max_err_k, err)

        if max_err_k > 1e-8:
            cumulant_additive = False
            print(f"  n={n}: cumulants NOT additive! max error = {max_err_k:.2e}")
        else:
            print(f"  n={n}: cumulant additivity confirmed, max error = {max_err_k:.2e}")

    if cumulant_additive:
        print("  Cumulants are additive, so Gaussian splitting follows EXACTLY.")
        print("  (Because G_s has kappa_2 = s and kappa_k = 0 for k >= 3)")
    else:
        print("  *** Cumulant additivity FAILED -- splitting may be approximate! ***")

    # Verify Gaussian cumulants
    print("\n--- Gaussian cumulant structure ---")
    for n in [3, 4, 5, 6, 8]:
        for s in [1.0, 2.0]:
            gs = gaussian_poly_coeffs(n, s)
            kg = finite_free_cumulants(gs)
            print(f"  n={n}, s={s}: kappa_1={kg[1]:.8f}, kappa_2={kg[2]:.8f}", end="")
            if len(kg) > 3:
                print(f", kappa_3={kg[3]:.8f}", end="")
            if len(kg) > 4:
                print(f", kappa_4={kg[4]:.8f}", end="")
            print()
            # Check: kappa_1 = 0, kappa_2 = s, kappa_k = 0 for k >= 3
            if abs(kg[1]) > 1e-10:
                print(f"    *** kappa_1 != 0! ***")
            if abs(kg[2] - s) > 1e-8:
                print(f"    *** kappa_2 != s! Expected {s}, got {kg[2]} ***")
            for k in range(3, min(len(kg), n+1)):
                if abs(kg[k]) > 1e-8:
                    print(f"    *** kappa_{k} != 0! Got {kg[k]} ***")

    # --- Test 4c: More general splitting ---
    print("\n--- General splitting: (p ⊞ G_s) ⊞ (q ⊞ G_t) = (p ⊞ q) ⊞ G_{s+t} ---")
    np.random.seed(123)
    gen_max_err = 0.0
    for n in [3, 4, 5]:
        for trial in range(20):
            roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n) * 1.5)
            roots_q = np.sort(np.random.randn(n) * 1.5 + np.arange(n))
            coeffs_p = coeffs_from_roots(roots_p)
            coeffs_q = coeffs_from_roots(roots_q)

            for s_val in [0.3, 1.0, 2.0]:
                for t_val in [0.5, 1.5, 3.0]:
                    gs = gaussian_poly_coeffs(n, s_val)
                    gt = gaussian_poly_coeffs(n, t_val)
                    gst = gaussian_poly_coeffs(n, s_val + t_val)

                    lhs = finite_free_convolution(
                        finite_free_convolution(coeffs_p, gs),
                        finite_free_convolution(coeffs_q, gt))
                    rhs = finite_free_convolution(
                        finite_free_convolution(coeffs_p, coeffs_q),
                        gst)

                    err = np.max(np.abs(lhs - rhs))
                    gen_max_err = max(gen_max_err, err)

    print(f"  Max error: {gen_max_err:.2e}")
    if gen_max_err < 1e-10:
        print("  CONFIRMED: General splitting holds exactly.")
    else:
        print(f"  *** General splitting has error {gen_max_err:.2e} ***")

    return max_err, cumulant_additive


# ============================================================
# ADDITIONAL CHECKS
# ============================================================

def check_debruijn_sign_carefully():
    """
    The classical de Bruijn identity says:
      d/dt H(X + sqrt(t)Z) = (1/2) J(X + sqrt(t)Z)
    where H is entropy and J is Fisher information.

    PROVER-13 claims: dS/dt = Phi_n(p)  (POSITIVE)
    where S is the log-Vandermonde entropy.

    But the classical identity has a factor of 1/2!
    Let's check if the constant is really 1 or 1/2 or something else.
    """
    print("\n" + "=" * 70)
    print("DEEP CHECK: Exact constant in de Bruijn identity")
    print("=" * 70)

    # Use very precise computation
    for n in [3, 4, 5, 6]:
        roots_p = np.arange(1, n+1, dtype=float) * 3  # well-separated
        coeffs_p = coeffs_from_roots(roots_p)
        phi = Phi_n(roots_p)
        S0 = S_entropy(roots_p)

        # Use very small dt and high-order finite differences
        dt_vals = [1e-3, 5e-4, 2.5e-4, 1.25e-4]
        dSdt_vals = []
        for dt in dt_vals:
            gt = gaussian_poly_coeffs(n, dt)
            ct = finite_free_convolution(coeffs_p, gt)
            rt = np.sort(roots_from_coeffs(ct).real)
            St = S_entropy(rt)
            dSdt_vals.append((St - S0) / dt)

        # Richardson extrapolation (order 2)
        rich1 = (4*dSdt_vals[1] - dSdt_vals[0]) / 3
        rich2 = (4*dSdt_vals[2] - dSdt_vals[1]) / 3
        rich3 = (16*rich2 - rich1) / 15  # order 4

        ratio = rich3 / phi
        print(f"  n={n}: dS/dt = {rich3:.8f}, Phi = {phi:.8f}, ratio dS/dt / Phi = {ratio:.10f}")
        # If ratio = 1, then dS/dt = Phi (PROVER-13's claim)
        # If ratio = 1/2, it's the classical factor

    # Also check using root dynamics
    print("\n  Root dynamics: if dr_i/dt = c * H_i, and dS/dt = Phi, what is c?")
    print("  dS/dt = sum_{i<j} (dr_i - dr_j)/(r_i - r_j) = c * sum_{i<j} (H_i-H_j)/(r_i-r_j)")
    print("  By the algebraic identity: sum_{i<j} (H_i-H_j)/(r_i-r_j) = Phi")
    print("  So dS/dt = c * Phi")
    print("  Combined with dS/dt = Phi (the claim), we get c = 1.")
    print("  Let's verify c = 1 from the root dynamics:")

    for n in [3, 4, 5, 6]:
        roots_p = np.arange(1, n+1, dtype=float) * 3
        coeffs_p = coeffs_from_roots(roots_p)
        H0 = H_values(roots_p)

        dt_vals = [1e-5, 1e-6, 1e-7]
        ratios = []
        for dt in dt_vals:
            gt = gaussian_poly_coeffs(n, dt)
            ct = finite_free_convolution(coeffs_p, gt)
            rt = np.sort(roots_from_coeffs(ct).real)
            dr = (rt - roots_p) / dt
            ratios.append(np.mean(dr / H0))

        # Extrapolate
        print(f"  n={n}: c at dt=1e-5: {ratios[0]:.8f}, dt=1e-6: {ratios[1]:.8f}, "
              f"dt=1e-7: {ratios[2]:.8f}")


def check_concavity_of_inv_phi():
    """
    PROVER-13 claims 1/Phi(p ⊞ G_t) is concave in t.
    This is crucial for the monotone gap argument.
    """
    print("\n" + "=" * 70)
    print("DEEP CHECK: Concavity of 1/Phi along heat flow")
    print("=" * 70)

    np.random.seed(2026)
    violations = 0
    total = 0

    # Use fine grid for precise concavity check
    for n in [3, 4, 5, 6, 7]:
        n_viol = 0
        n_tot = 0
        for trial in range(100):
            roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n) * 2)
            if np.min(np.diff(roots_p)) < 0.1:
                continue

            coeffs_p = coeffs_from_roots(roots_p)
            t_vals = np.linspace(0.001, 5.0, 500)
            inv_phi = []

            for t in t_vals:
                gt = gaussian_poly_coeffs(n, t)
                ct = finite_free_convolution(coeffs_p, gt)
                rt = roots_from_coeffs(ct)
                if np.any(np.abs(np.imag(rt)) > 1e-8):
                    inv_phi.append(float('nan'))
                    continue
                rt = np.sort(rt.real)
                if np.min(np.diff(rt)) < 1e-12:
                    inv_phi.append(float('nan'))
                    continue
                phi = Phi_n(rt)
                inv_phi.append(1.0 / phi if phi > 0 else float('nan'))

            inv_phi = np.array(inv_phi)
            # Concavity: f(t_mid) >= (f(t_left) + f(t_right)) / 2
            # Use second differences as proxy
            d2 = np.diff(inv_phi, 2)
            valid = ~np.isnan(d2)
            if np.sum(valid) < 20:
                continue

            n_tot += 1
            total += 1
            max_d2 = np.max(d2[valid])
            if max_d2 > 1e-5:  # strict threshold
                n_viol += 1
                violations += 1

        print(f"  n={n}: {n_viol} concavity violations / {n_tot}")

    print(f"  Total: {violations} / {total}")
    return violations, total


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("=" * 70)
    print("VERIFIER-12: ADVERSARIAL VERIFICATION OF PROVER-13 CLAIMS")
    print("=" * 70)
    print()

    # Verify infrastructure first
    mss_ok = verify_mss_convolution()
    if not mss_ok:
        print("*** MSS convolution verification FAILED -- all results suspect ***")

    # Claim 1: De Bruijn
    debruijn_results = verify_debruijn()

    # Deep check on de Bruijn constant
    check_debruijn_sign_carefully()

    # Claim 2: EPI
    epi_violations, epi_total, epi_min_ratios = verify_epi()

    # Claim 3: Monotone gap
    gap_violations, gap_total = verify_monotone_gap()

    # Concavity check
    concavity_viols, concavity_total = check_concavity_of_inv_phi()

    # Claim 4: Gaussian splitting
    splitting_err, cumulant_ok = verify_gaussian_splitting()

    # ============================================================
    # FINAL SUMMARY
    # ============================================================
    print("\n" + "=" * 70)
    print("VERIFIER-12 FINAL SUMMARY")
    print("=" * 70)

    print(f"\nCLAIM 1 (de Bruijn): dS/dt = Phi_n(p)")
    print(f"  Sign errors: {debruijn_results['sign_errors']}")
    print(f"  Adversarial failures: {len(debruijn_results['adversarial_failures'])}")
    print(f"  Max relative error: {debruijn_results['max_rel_error']:.2e}")

    print(f"\nCLAIM 2 (EPI): N(p ⊞ q) >= N(p) + N(q)")
    print(f"  Violations: {epi_violations}/{epi_total}")
    print(f"  Min ratios: {epi_min_ratios}")

    print(f"\nCLAIM 3 (Monotone Gap):")
    print(f"  Gap increases detected: {gap_violations}/{gap_total}")

    print(f"\nCLAIM 3 support (Concavity of 1/Phi):")
    print(f"  Concavity violations: {concavity_viols}/{concavity_total}")

    print(f"\nCLAIM 4 (Gaussian Splitting):")
    print(f"  Max numerical error: {splitting_err:.2e}")
    print(f"  Cumulant additivity: {'CONFIRMED' if cumulant_ok else 'FAILED'}")
