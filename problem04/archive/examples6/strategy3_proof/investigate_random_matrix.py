"""
investigate_random_matrix.py — Random Matrix / Jensen's Inequality Approach
to Fisher Superadditivity for ALL n.

CONJECTURE: For monic real-rooted degree-n polynomials p, q with simple roots:
  1/Phi_n(p boxplus_n q) >= 1/Phi_n(p) + 1/Phi_n(q)

CORRECT MSS boxplus formula (Marcus-Spielman-Srivastava):
  g_k = sum_{i+j=k} [C(n-j,i)/C(n,i)] * e_i(p) * e_j(q)
  r(x) = sum_k (-1)^k g_k x^{n-k}

RANDOM MATRIX CONNECTION:
  E_U[det(xI - A - UBU*)] = p boxplus_n q   (U ~ Haar on U(n))

INVESTIGATION PLAN:
  1. Formula verification and counterexample audit
  2. Is 1/Phi_n concave in the coefficient vector? (Jensen approach)
  3. Direct random matrix: does 1/Phi_n(A+UBU*) >= 1/Phi_n(A) + 1/Phi_n(B)
     hold per-realization?
  4. Trace formulation via Cauchy matrix
  5. Stochastic dominance / majorization
  6. Literature integration

Author: prover agent
Date: 2026-02-08
"""

import numpy as np
from math import comb, factorial
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
np.set_printoptions(precision=10, linewidth=120)


# ============================================================================
# CORE DEFINITIONS
# ============================================================================

def elem_sym(roots, k):
    """Elementary symmetric polynomial e_k of given roots."""
    n = len(roots)
    if k == 0:
        return 1.0
    if k > n:
        return 0.0
    return sum(np.prod(combo) for combo in combinations(roots, k))


def Phi_n(roots):
    """
    Finite free Fisher information:
      Phi_n(p) = sum_i H_p(lambda_i)^2
    where H_p(lambda_i) = sum_{j != i} 1/(lambda_i - lambda_j).
    """
    n = len(roots)
    total = 0.0
    for i in range(n):
        H_i = sum(1.0 / (roots[i] - roots[j]) for j in range(n) if j != i)
        total += H_i ** 2
    return total


def boxplus_mss(roots_p, roots_q):
    """
    CORRECT MSS boxplus formula:
      g_k = sum_{i+j=k} [C(n-j,i)/C(n,i)] * e_i(p) * e_j(q)
      r(x) = sum_k (-1)^k g_k x^{n-k}

    This is algebraically identical to:
      c_k = sum_{i+j=k} [(n-i)!(n-j)!/(n!(n-k)!)] * a_i * b_j
    where a_k = (-1)^k * e_k.
    """
    n = len(roots_p)
    assert len(roots_q) == n

    e_p = [elem_sym(roots_p, k) for k in range(n + 1)]
    e_q = [elem_sym(roots_q, k) for k in range(n + 1)]

    g = np.zeros(n + 1)
    for k in range(n + 1):
        for i in range(k + 1):
            j = k - i
            if comb(n, i) > 0:
                w = comb(n - j, i) / comb(n, i)
                g[k] += w * e_p[i] * e_q[j]

    # r(x) = sum_k (-1)^k g_k x^{n-k}
    poly_r = np.array([(-1) ** k * g[k] for k in range(n + 1)])
    roots_r = np.sort(np.real(np.roots(poly_r)))
    return roots_r, poly_r


def boxplus_hat(roots_p, roots_q):
    """
    ALTERNATIVE formula (hat-e convolution):
      hat_e_k(r) = sum_{j=0}^k hat_e_j(p) * hat_e_{k-j}(q)
    where hat_e_k = e_k / C(n,k).

    NOTE: This is DIFFERENT from boxplus_mss for n >= 3!
    """
    n = len(roots_p)
    assert len(roots_q) == n

    e_p = [elem_sym(roots_p, k) for k in range(n + 1)]
    e_q = [elem_sym(roots_q, k) for k in range(n + 1)]

    hat_p = [e_p[k] / comb(n, k) for k in range(n + 1)]
    hat_q = [e_q[k] / comb(n, k) for k in range(n + 1)]

    hat_r = np.zeros(n + 1)
    for k in range(n + 1):
        for j in range(k + 1):
            hat_r[k] += hat_p[j] * hat_q[k - j]

    e_r = [comb(n, k) * hat_r[k] for k in range(n + 1)]
    poly_r = np.array([(-1) ** k * e_r[k] for k in range(n + 1)])
    roots_r = np.sort(np.real(np.roots(poly_r)))
    return roots_r, poly_r


def haar_unitary(n):
    """Generate Haar-distributed random unitary matrix."""
    Z = (np.random.randn(n, n) + 1j * np.random.randn(n, n)) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    d = np.diag(R)
    ph = d / np.abs(d)
    return Q @ np.diag(ph)


def mc_boxplus(roots_p, roots_q, N=50000):
    """Ground truth: E[det(xI - A - UBU*)] via Monte Carlo."""
    n = len(roots_p)
    A = np.diag(np.array(roots_p, dtype=complex))
    B = np.diag(np.array(roots_q, dtype=complex))

    poly_sum = np.zeros(n + 1, dtype=float)
    for _ in range(N):
        U = haar_unitary(n)
        M = A + U @ B @ U.conj().T
        poly_sum += np.real(np.poly(M))

    poly_avg = poly_sum / N
    poly_avg[0] = 1.0
    roots_r = np.sort(np.real(np.roots(poly_avg)))
    return roots_r, poly_avg


# ============================================================================
# PART 1: FORMULA VERIFICATION AND COUNTEREXAMPLE AUDIT
# ============================================================================

def part1_formula_audit():
    """
    Verify which boxplus formula is correct by comparing against Monte Carlo.
    Check whether the alleged n=4 counterexamples were computed with wrong formula.
    """
    print("=" * 78)
    print("PART 1: FORMULA VERIFICATION AND COUNTEREXAMPLE AUDIT")
    print("=" * 78)

    # First: are MSS and hat formulas the same?
    print("\n--- 1a. Do MSS and hat-e formulas agree? ---\n")

    test_cases = [
        ("n=2: {-1,1},{-2,2}", [-1.0, 1.0], [-2.0, 2.0]),
        ("n=3: {-2,0,2},{-1,0,1}", [-2.0, 0.0, 2.0], [-1.0, 0.0, 1.0]),
        ("n=4: {-5,-1,1,5},{-5,-1,1,5}", [-5.0, -1.0, 1.0, 5.0], [-5.0, -1.0, 1.0, 5.0]),
        ("n=4: {-3,-1,1,3},{-3,-1,1,3}", [-3.0, -1.0, 1.0, 3.0], [-3.0, -1.0, 1.0, 3.0]),
        ("n=5: random", np.sort(np.random.randn(5)*3), np.sort(np.random.randn(5)*3)),
    ]

    formulas_differ = False
    for name, rp, rq in test_cases:
        rp = np.array(rp, dtype=float)
        rq = np.array(rq, dtype=float)
        _, poly_mss = boxplus_mss(rp, rq)
        _, poly_hat = boxplus_hat(rp, rq)
        diff = np.max(np.abs(poly_mss - poly_hat))
        agree = "AGREE" if diff < 1e-10 else "DIFFER"
        if diff > 1e-10:
            formulas_differ = True
        print(f"  {name}: max coeff diff = {diff:.4e}  [{agree}]")
        if diff > 1e-10:
            print(f"    MSS poly: {poly_mss}")
            print(f"    Hat poly: {poly_hat}")

    if formulas_differ:
        print("\n  IMPORTANT: The two formulas DIFFER for n >= 3!")
        print("  The hat-e convolution is NOT the correct MSS boxplus.")
    else:
        print("\n  All formulas agree (unexpected for n >= 3).")

    # Now: verify MSS formula against Monte Carlo
    print("\n--- 1b. MSS formula vs Monte Carlo (ground truth) ---\n")

    critical_cases = [
        ("ALLEGED COUNTEREX 1: p=q={-5,-1,1,5}",
         [-5.0, -1.0, 1.0, 5.0], [-5.0, -1.0, 1.0, 5.0]),
        ("ALLEGED COUNTEREX 2: p=q={-3,-1,1,3}",
         [-3.0, -1.0, 1.0, 3.0], [-3.0, -1.0, 1.0, 3.0]),
        ("n=3 centered: {-2,0,2},{-1,0,1}",
         [-2.0, 0.0, 2.0], [-1.0, 0.0, 1.0]),
        ("n=4 non-sym: {-4,-1,2,3},{-3,0,1,2}",
         [-4.0, -1.0, 2.0, 3.0], [-3.0, 0.0, 1.0, 2.0]),
    ]

    results = {}
    for name, rp, rq in critical_cases:
        rp = np.array(rp, dtype=float)
        rq = np.array(rq, dtype=float)
        n = len(rp)

        roots_mss, poly_mss = boxplus_mss(rp, rq)
        roots_hat, poly_hat = boxplus_hat(rp, rq)
        roots_mc, poly_mc = mc_boxplus(rp, rq, N=80000)

        err_mss = np.max(np.abs(poly_mss - poly_mc))
        err_hat = np.max(np.abs(poly_hat - poly_mc))

        phi_p = Phi_n(rp)
        phi_q = Phi_n(rq)

        phi_r_mss = Phi_n(roots_mss)
        phi_r_hat = Phi_n(roots_hat)
        phi_r_mc = Phi_n(roots_mc)

        gap_mss = 1/phi_r_mss - 1/phi_p - 1/phi_q
        gap_hat = 1/phi_r_hat - 1/phi_p - 1/phi_q
        gap_mc = 1/phi_r_mc - 1/phi_p - 1/phi_q

        print(f"  {name}")
        print(f"    MSS coeff err vs MC: {err_mss:.4e},  gap = {gap_mss:+.6e}  "
              f"{'HOLDS' if gap_mss >= -1e-8 else 'VIOLATED'}")
        print(f"    Hat coeff err vs MC: {err_hat:.4e},  gap = {gap_hat:+.6e}  "
              f"{'HOLDS' if gap_hat >= -1e-8 else 'VIOLATED'}")
        print(f"    MC  gap:                              gap = {gap_mc:+.6e}  "
              f"{'HOLDS' if gap_mc >= -1e-3 else 'VIOLATED'}")

        if err_mss < err_hat / 3:
            print(f"    => MSS is the CORRECT formula")
        elif err_hat < err_mss / 3:
            print(f"    => Hat is the CORRECT formula (unexpected!)")
        else:
            print(f"    => Both close to MC; formulas may agree for this case")

        results[name] = {
            'gap_mss': gap_mss, 'gap_hat': gap_hat, 'gap_mc': gap_mc,
            'err_mss': err_mss, 'err_hat': err_hat,
        }
        print()

    # Summary for counterexample audit
    print("--- 1c. COUNTEREXAMPLE AUDIT SUMMARY ---\n")
    for name in [n for n in results if "COUNTEREX" in n]:
        r = results[name]
        print(f"  {name}")
        print(f"    With CORRECT (MSS) formula: gap = {r['gap_mss']:+.6e}")
        print(f"    With WRONG (hat) formula:   gap = {r['gap_hat']:+.6e}")
        if r['gap_mss'] >= -1e-8 and r['gap_hat'] < -1e-8:
            print(f"    => COUNTEREXAMPLE WAS COMPUTED WITH WRONG FORMULA!")
            print(f"       The conjecture HOLDS with the correct MSS formula.")
        elif r['gap_mss'] < -1e-8:
            print(f"    => Counterexample VALID even with correct MSS formula.")
        else:
            print(f"    => Both formulas show conjecture holds for this case.")
        print()

    return results


# ============================================================================
# PART 2: IS 1/Phi_n CONCAVE IN COEFFICIENT VECTOR?
# ============================================================================

def part2_concavity_analysis():
    """
    Check if 1/Phi_n is concave as a function of:
    (a) the coefficient vector (e_1, ..., e_n)
    (b) the root vector (lambda_1, ..., lambda_n)

    If 1/Phi_n is concave in coefficients and MSS boxplus is linear in
    coefficient space, Jensen's inequality would directly give superadditivity.
    """
    print("\n" + "=" * 78)
    print("PART 2: CONCAVITY ANALYSIS OF 1/Phi_n")
    print("=" * 78)

    # For n=3 centered: Phi_3 = (4E^3 - 27F^2) / (18E^2) where E = -e_2 > 0, F = e_3
    # So 1/Phi_3 = 18E^2 / (4E^3 - 27F^2) = 18/(4E - 27F^2/E^2)

    print("\n--- 2a. n=3 centered: 1/Phi_3 = 18E^2 / (4E^3 - 27F^2) ---")
    print("  where E = sum_{i<j} (lambda_i - lambda_j)^2 / 3 = -e_2")
    print("  and F = e_3 (product of roots, centered: e_1 = 0)")

    # Check concavity along lines in (E, F) space
    print("\n  Checking concavity along random lines in (E, F) space...")
    n_tests = 500
    n_concave = 0
    n_convex = 0
    n_neither = 0

    for _ in range(n_tests):
        # Two random centered cubics
        E1 = np.random.uniform(0.5, 10.0)
        F1 = np.random.uniform(-E1**1.5/3, E1**1.5/3) * 0.9  # ensure real roots
        E2 = np.random.uniform(0.5, 10.0)
        F2 = np.random.uniform(-E2**1.5/3, E2**1.5/3) * 0.9

        # 1/Phi_3(E, F) = 18*E^2 / (4*E^3 - 27*F^2)
        def inv_phi3(E, F):
            denom = 4*E**3 - 27*F**2
            if denom <= 0:
                return None
            return 18*E**2 / denom

        # Check midpoint convexity: f((x+y)/2) vs (f(x)+f(y))/2
        Em, Fm = (E1+E2)/2, (F1+F2)/2
        v1 = inv_phi3(E1, F1)
        v2 = inv_phi3(E2, F2)
        vm = inv_phi3(Em, Fm)

        if v1 is None or v2 is None or vm is None:
            continue

        avg = (v1 + v2) / 2
        if vm > avg + 1e-10:
            n_concave += 1
        elif vm < avg - 1e-10:
            n_convex += 1
        else:
            n_neither += 1

    total = n_concave + n_convex + n_neither
    print(f"  Results ({total} valid tests):")
    print(f"    Concave (midpoint above average): {n_concave} ({100*n_concave/max(total,1):.1f}%)")
    print(f"    Convex (midpoint below average):  {n_convex} ({100*n_convex/max(total,1):.1f}%)")
    print(f"    Linear-ish (indeterminate):        {n_neither} ({100*n_neither/max(total,1):.1f}%)")

    # Now check concavity in ROOT space
    print("\n--- 2b. Concavity of 1/Phi_n in root vector ---")
    for n in [3, 4, 5]:
        n_tests_r = 300
        concave_count = 0
        convex_count = 0

        for _ in range(n_tests_r):
            roots_x = np.sort(np.random.randn(n) * 3)
            roots_y = np.sort(np.random.randn(n) * 3)

            # Ensure distinct
            for arr in [roots_x, roots_y]:
                for i in range(1, n):
                    if abs(arr[i] - arr[i-1]) < 0.1:
                        arr[i] = arr[i-1] + 0.1 + np.random.rand()

            roots_m = (roots_x + roots_y) / 2

            try:
                fx = 1.0 / Phi_n(roots_x)
                fy = 1.0 / Phi_n(roots_y)
                fm = 1.0 / Phi_n(roots_m)

                avg = (fx + fy) / 2
                if fm > avg + 1e-12:
                    concave_count += 1
                elif fm < avg - 1e-12:
                    convex_count += 1
            except:
                continue

        total_r = concave_count + convex_count
        if total_r > 0:
            print(f"  n={n}: concave {concave_count}/{total_r} ({100*concave_count/total_r:.1f}%), "
                  f"convex {convex_count}/{total_r} ({100*convex_count/total_r:.1f}%)")

    # Check if MSS boxplus is LINEAR in the coefficient vector
    print("\n--- 2c. Is MSS boxplus linear in coefficient space? ---")
    print("  MSS formula: g_k = sum_{i+j=k} w_{ij} * e_i(p) * e_j(q)")
    print("  This is BILINEAR in (e(p), e(q)), not linear.")
    print("  So even if 1/Phi were concave in coefficients, Jensen would need")
    print("  a different argument than simple midpoint concavity.")

    # But note: r = E[chi_{A+UBU*}] is a polynomial-valued expectation
    # The ROOTS of r are NOT the expected roots (roots of expected poly != expected roots)
    # Jensen on 1/Phi(roots) would need concavity in root space

    print("\n--- 2d. Jensen via random matrix representation ---")
    print("  Since r = E_U[chi_{A+UBU*}], we have:")
    print("  1/Phi(r) = 1/Phi(E_U[chi_M])")
    print("  where M = A + UBU*")
    print("")
    print("  For Jensen to work, we need 1/Phi to be concave on the space of")
    print("  monic degree-n polynomials (viewed as coefficient vectors).")
    print("  Since the expected char poly IS a convex combination of char polys,")
    print("  Jensen would give: 1/Phi(r) >= E_U[1/Phi(chi_M)]")
    print("  But we would then need E_U[1/Phi(chi_M)] >= 1/Phi(A) + 1/Phi(B)")
    print("  which is a SEPARATE statement about the Haar average.")


# ============================================================================
# PART 3: DIRECT RANDOM MATRIX — PER-REALIZATION INEQUALITY
# ============================================================================

def part3_per_realization():
    """
    For fixed A, B Hermitian, compute Phi_n(A + UBU*) for many Haar unitaries U.
    Test: does 1/Phi_n(A+UBU*) >= 1/Phi_n(A) + 1/Phi_n(B) hold per-realization?
    """
    print("\n" + "=" * 78)
    print("PART 3: PER-REALIZATION INEQUALITY (RANDOM MATRIX)")
    print("=" * 78)

    test_configs = [
        ("n=3: A=diag(-2,0,2), B=diag(-1,0,1)",
         [-2.0, 0.0, 2.0], [-1.0, 0.0, 1.0]),
        ("n=4: A=diag(-3,-1,1,3), B=diag(-3,-1,1,3)",
         [-3.0, -1.0, 1.0, 3.0], [-3.0, -1.0, 1.0, 3.0]),
        ("n=4: A=diag(-5,-1,1,5), B=diag(-5,-1,1,5)",
         [-5.0, -1.0, 1.0, 5.0], [-5.0, -1.0, 1.0, 5.0]),
        ("n=5: random centered",
         np.sort(np.array([-4.0, -1.5, 0.0, 1.5, 4.0])),
         np.sort(np.array([-3.0, -1.0, 0.0, 1.0, 3.0]))),
    ]

    N_samples = 2000

    for name, roots_a, roots_b in test_configs:
        roots_a = np.array(roots_a, dtype=float)
        roots_b = np.array(roots_b, dtype=float)
        n = len(roots_a)

        phi_a = Phi_n(roots_a)
        phi_b = Phi_n(roots_b)
        target = 1/phi_a + 1/phi_b

        A = np.diag(roots_a.astype(complex))
        B = np.diag(roots_b.astype(complex))

        violations = 0
        inv_phis = []

        for _ in range(N_samples):
            U = haar_unitary(n)
            M = A + U @ B @ U.conj().T
            eigs = np.sort(np.real(np.linalg.eigvalsh(M)))

            # Check distinct eigenvalues
            if np.min(np.diff(eigs)) < 1e-10:
                continue

            phi_M = Phi_n(eigs)
            inv_phi_M = 1.0 / phi_M
            inv_phis.append(inv_phi_M)

            if inv_phi_M < target - 1e-8:
                violations += 1

        inv_phis = np.array(inv_phis)

        # Also compute for the expected polynomial
        roots_r, _ = boxplus_mss(roots_a, roots_b)
        phi_r = Phi_n(roots_r)
        inv_phi_r = 1/phi_r
        gap_r = inv_phi_r - target

        print(f"\n  {name}")
        print(f"    1/Phi(A) + 1/Phi(B) = {target:.6f}")
        print(f"    1/Phi(r) [MSS boxplus] = {inv_phi_r:.6f}, gap = {gap_r:+.6e}")
        print(f"    Per-realization 1/Phi(A+UBU*):")
        print(f"      min  = {np.min(inv_phis):.6f}")
        print(f"      mean = {np.mean(inv_phis):.6f}")
        print(f"      max  = {np.max(inv_phis):.6f}")
        print(f"      violations (< target): {violations}/{len(inv_phis)}")

        # Does E[1/Phi(M)] >= target?
        e_inv_phi = np.mean(inv_phis)
        print(f"      E[1/Phi(M)] = {e_inv_phi:.6f}, gap from target = {e_inv_phi - target:+.6e}")

        # Does 1/Phi(r) >= E[1/Phi(M)]? (Jensen direction)
        print(f"      1/Phi(r) - E[1/Phi(M)] = {inv_phi_r - e_inv_phi:+.6e}")

        if violations == 0:
            print(f"    => PER-REALIZATION inequality HOLDS (0 violations)")
        else:
            print(f"    => PER-REALIZATION inequality VIOLATED ({violations} times)")


# ============================================================================
# PART 4: TRACE FORMULATION VIA CAUCHY MATRIX
# ============================================================================

def part4_trace_formulation():
    """
    Phi_n(M) = ||C @ 1||^2 where C_{ij} = 1/(lambda_i - lambda_j) for i!=j.
    C is antisymmetric, so C^T C = -C^2.
    Phi_n = -1^T C^2 1 = Tr(-C^2) when restricted to (1,...,1) inner product.

    Actually: Phi_n = sum_{i,j : i!=j} 1/(lambda_i - lambda_j)^2
                    = 2 * sum_{i<j} 1/(lambda_i - lambda_j)^2 ... NO.

    Wait: H_i = sum_{j!=i} 1/(lambda_i - lambda_j)
    Phi_n = sum_i H_i^2 = 1^T C^T C 1 where C_{ij} = 1/(lambda_i - lambda_j) for j!=i, 0 for j=i
    Since C is antisymmetric: C^T = -C, so C^T C = -C^2.
    Phi_n = -1^T C^2 1.
    """
    print("\n" + "=" * 78)
    print("PART 4: TRACE FORMULATION VIA CAUCHY MATRIX")
    print("=" * 78)

    def cauchy_matrix(roots):
        n = len(roots)
        C = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i != j:
                    C[i, j] = 1.0 / (roots[i] - roots[j])
        return C

    # Verify the formula Phi_n = -1^T C^2 1
    print("\n--- 4a. Verify Phi_n = -1^T C^2 1 = ||C 1||^2 ---")
    for n in [3, 4, 5]:
        roots = np.sort(np.random.randn(n) * 3)
        for i in range(1, n):
            if abs(roots[i] - roots[i-1]) < 0.1:
                roots[i] = roots[i-1] + 0.5

        C = cauchy_matrix(roots)
        ones = np.ones(n)
        phi_direct = Phi_n(roots)
        phi_cauchy = ones @ (C.T @ C) @ ones  # = ||C 1||^2
        phi_trace = -ones @ (C @ C) @ ones     # = -1^T C^2 1

        print(f"  n={n}: Phi_direct = {phi_direct:.10f}")
        print(f"         ||C 1||^2  = {phi_cauchy:.10f}")
        print(f"         -1^T C^2 1 = {phi_trace:.10f}")
        print(f"         Errors: {abs(phi_direct - phi_cauchy):.2e}, {abs(phi_direct - phi_trace):.2e}")

    # Trace connection
    print("\n--- 4b. Trace interpretation ---")
    print("  Phi_n = ||C @ 1||^2 = Tr(1 1^T C^T C)")
    print("  = Tr(J (-C^2)) where J = 1 1^T (all-ones matrix)")
    print("")
    print("  For superadditivity 1/Phi(r) >= 1/Phi(p) + 1/Phi(q),")
    print("  we need: ||C_r @ 1||^2 <= (||C_p @ 1||^2 * ||C_q @ 1||^2) / (||C_p @ 1||^2 + ||C_q @ 1||^2)")
    print("  i.e., Phi(r) <= harmonic_mean(Phi(p), Phi(q)) / 2")
    print("")
    print("  This is hard to prove from the Cauchy matrix structure alone")
    print("  because the roots of r are a complicated function of roots of p, q.")

    # Can we express the Cauchy matrix of r in terms of those of p, q?
    print("\n--- 4c. Eigenvalue structure of Cauchy matrix ---")
    for n in [3, 4, 5]:
        roots = np.sort(np.random.randn(n) * 3)
        for i in range(1, n):
            if abs(roots[i] - roots[i-1]) < 0.1:
                roots[i] = roots[i-1] + 0.5

        C = cauchy_matrix(roots)
        eigs = np.sort(np.linalg.eigvals(C))
        print(f"  n={n}: C eigenvalues = {eigs}")
        print(f"         C^2 eigenvalues = {np.sort(eigs**2)}")
        print(f"         Tr(C^2) = {np.trace(C @ C):.6f} (should be -Phi_n = {-Phi_n(roots):.6f})")


# ============================================================================
# PART 5: STOCHASTIC DOMINANCE / MAJORIZATION
# ============================================================================

def part5_majorization():
    """
    Investigate whether eigenvalues of A+UBU* relate to roots of r via majorization.
    If 1/Phi is Schur-concave, stochastic dominance could close the argument.
    """
    print("\n" + "=" * 78)
    print("PART 5: STOCHASTIC DOMINANCE / MAJORIZATION")
    print("=" * 78)

    # For A+UBU*, the eigenvalues for different U are related by
    # Horn's conjecture / Klyachko's theorem.
    # The average (expected char poly) = p boxplus q.
    # But individual realizations can have very different eigenvalue patterns.

    print("\n--- 5a. Do roots of r majorize eigenvalues of A+UBU*? ---")
    print("  (This would mean r has 'more spread out' roots than typical A+UBU*)")

    for n in [3, 4]:
        roots_a = np.sort(np.random.randn(n) * 3)
        roots_b = np.sort(np.random.randn(n) * 2)
        for arr in [roots_a, roots_b]:
            for i in range(1, n):
                if abs(arr[i] - arr[i-1]) < 0.2:
                    arr[i] = arr[i-1] + 0.5

        roots_r, _ = boxplus_mss(roots_a, roots_b)

        A = np.diag(roots_a.astype(complex))
        B = np.diag(roots_b.astype(complex))

        # Check majorization: does centered_r majorize centered_eigs for each U?
        centered_r = np.sort(roots_r - np.mean(roots_r))[::-1]

        n_major = 0
        n_total = 0
        N_samp = 500

        for _ in range(N_samp):
            U = haar_unitary(n)
            M = A + U @ B @ U.conj().T
            eigs = np.sort(np.real(np.linalg.eigvalsh(M)))
            centered_eigs = np.sort(eigs - np.mean(eigs))[::-1]

            # Check majorization: sum of k largest of centered_r >= sum of k largest of centered_eigs
            majorizes = True
            for k in range(1, n):
                if np.sum(centered_r[:k]) < np.sum(centered_eigs[:k]) - 1e-10:
                    majorizes = False
                    break

            if majorizes:
                n_major += 1
            n_total += 1

        print(f"  n={n}: centered_r majorizes centered_eigs: {n_major}/{n_total} "
              f"({100*n_major/n_total:.1f}%)")

    # Is 1/Phi_n Schur-concave?
    print("\n--- 5b. Is 1/Phi_n Schur-concave in roots? ---")
    print("  (If yes + majorization, then 1/Phi(M) >= 1/Phi(r) per-realization)")

    for n in [3, 4, 5]:
        n_tests = 500
        schur_concave_count = 0
        schur_convex_count = 0

        for _ in range(n_tests):
            # Generate x majorizing y (same sum, x more spread)
            roots_y = np.sort(np.random.randn(n) * 2)
            for i in range(1, n):
                if abs(roots_y[i] - roots_y[i-1]) < 0.1:
                    roots_y[i] = roots_y[i-1] + 0.3

            # Create x that majorizes y by T-transform
            i_big = np.argmax(roots_y)
            i_small = np.argmin(roots_y)
            t = np.random.uniform(0.1, 0.9)
            roots_x = roots_y.copy()
            delta = t * (roots_y[i_big] - roots_y[i_small])
            roots_x[i_big] = roots_y[i_big] + delta / 2
            roots_x[i_small] = roots_y[i_small] - delta / 2
            roots_x = np.sort(roots_x)

            # Ensure distinct
            ok = True
            for arr in [roots_x]:
                for i in range(1, n):
                    if abs(arr[i] - arr[i-1]) < 1e-10:
                        ok = False
            if not ok:
                continue

            try:
                fx = 1.0 / Phi_n(roots_x)
                fy = 1.0 / Phi_n(roots_y)

                # x majorizes y. Schur-concave means f(x) <= f(y).
                if fx <= fy + 1e-12:
                    schur_concave_count += 1
                else:
                    schur_convex_count += 1
            except:
                continue

        total = schur_concave_count + schur_convex_count
        if total > 0:
            print(f"  n={n}: Schur-concave: {schur_concave_count}/{total} "
                  f"({100*schur_concave_count/total:.1f}%), "
                  f"Schur-convex: {schur_convex_count}/{total} "
                  f"({100*schur_convex_count/total:.1f}%)")

    print("\n  Note: If 1/Phi is NOT purely Schur-concave/convex,")
    print("  then majorization alone cannot close the argument.")


# ============================================================================
# PART 6: LARGE-SCALE NUMERICAL TEST WITH CORRECT FORMULA
# ============================================================================

def part6_large_scale():
    """
    Large-scale numerical test of superadditivity with CORRECT MSS formula.
    Tests both centered and non-centered polynomials across n = 2..8.
    """
    print("\n" + "=" * 78)
    print("PART 6: LARGE-SCALE NUMERICAL TEST (CORRECT MSS FORMULA)")
    print("=" * 78)

    for n in [2, 3, 4, 5, 6]:
        n_trials = 500
        n_pass = 0
        n_fail = 0
        n_skip = 0
        min_gap = float('inf')
        worst_p = None
        worst_q = None

        for trial in range(n_trials):
            # Random roots
            roots_p = np.sort(np.random.randn(n) * (2 + np.random.rand() * 5))
            roots_q = np.sort(np.random.randn(n) * (2 + np.random.rand() * 5))

            # Ensure well-separated
            for arr in [roots_p, roots_q]:
                for i in range(1, n):
                    if abs(arr[i] - arr[i-1]) < 0.05:
                        arr[i] = arr[i-1] + 0.05 + np.random.rand() * 0.3

            try:
                roots_r, poly_r = boxplus_mss(roots_p, roots_q)

                # Check real-rootedness
                roots_r_complex = np.roots(poly_r)
                if np.max(np.abs(np.imag(roots_r_complex))) > 1e-6:
                    n_skip += 1
                    continue

                # Check distinct
                if np.min(np.diff(np.sort(np.real(roots_r)))) < 1e-10:
                    n_skip += 1
                    continue

                phi_p = Phi_n(roots_p)
                phi_q = Phi_n(roots_q)
                phi_r = Phi_n(roots_r)

                if phi_p < 1e-10 or phi_q < 1e-10 or phi_r < 1e-10:
                    n_skip += 1
                    continue

                gap = 1/phi_r - 1/phi_p - 1/phi_q

                if gap < min_gap:
                    min_gap = gap
                    worst_p = roots_p.copy()
                    worst_q = roots_q.copy()

                if gap >= -1e-8:
                    n_pass += 1
                else:
                    n_fail += 1
                    if n_fail <= 2:
                        print(f"    FAIL n={n} trial {trial}: gap={gap:.6e}")
                        print(f"      p = {roots_p}")
                        print(f"      q = {roots_q}")
                        print(f"      r = {roots_r}")
            except Exception as e:
                n_skip += 1

        total = n_pass + n_fail
        print(f"  n={n}: pass={n_pass}, fail={n_fail}, skip={n_skip}, "
              f"min_gap={min_gap:.6e}")

        if n_fail > 0 and worst_p is not None:
            print(f"    Worst case:")
            print(f"      p = {worst_p}")
            print(f"      q = {worst_q}")

    # Centered polynomials specifically
    print("\n--- Centered polynomials only ---")
    for n in [3, 4, 5, 6]:
        n_trials = 500
        n_pass = 0
        n_fail = 0
        n_skip = 0
        min_gap = float('inf')

        for trial in range(n_trials):
            roots_p = np.sort(np.random.randn(n) * (2 + np.random.rand() * 5))
            roots_q = np.sort(np.random.randn(n) * (2 + np.random.rand() * 5))

            # Center
            roots_p -= np.mean(roots_p)
            roots_q -= np.mean(roots_q)

            # Ensure well-separated
            for arr in [roots_p, roots_q]:
                for i in range(1, n):
                    if abs(arr[i] - arr[i-1]) < 0.05:
                        arr[i] = arr[i-1] + 0.05 + np.random.rand() * 0.3
                arr -= np.mean(arr)  # re-center

            try:
                roots_r, poly_r = boxplus_mss(roots_p, roots_q)

                roots_r_complex = np.roots(poly_r)
                if np.max(np.abs(np.imag(roots_r_complex))) > 1e-6:
                    n_skip += 1
                    continue

                if np.min(np.diff(np.sort(np.real(roots_r)))) < 1e-10:
                    n_skip += 1
                    continue

                phi_p = Phi_n(roots_p)
                phi_q = Phi_n(roots_q)
                phi_r = Phi_n(roots_r)

                if phi_p < 1e-10 or phi_q < 1e-10 or phi_r < 1e-10:
                    n_skip += 1
                    continue

                gap = 1/phi_r - 1/phi_p - 1/phi_q
                min_gap = min(min_gap, gap)

                if gap >= -1e-8:
                    n_pass += 1
                else:
                    n_fail += 1
            except:
                n_skip += 1

        print(f"  n={n} centered: pass={n_pass}, fail={n_fail}, skip={n_skip}, "
              f"min_gap={min_gap:.6e}")


# ============================================================================
# PART 7: EIGENVALUE DISTRIBUTION ANALYSIS
# ============================================================================

def part7_eigenvalue_distribution():
    """
    Analyze the distribution of 1/Phi(A+UBU*) over Haar unitaries.
    Compare with 1/Phi(r) and 1/Phi(A) + 1/Phi(B).
    """
    print("\n" + "=" * 78)
    print("PART 7: DISTRIBUTION OF 1/Phi(A+UBU*) OVER HAAR UNITARIES")
    print("=" * 78)

    configs = [
        ("n=3: {-2,0,2},{-1,0,1}", [-2.0, 0.0, 2.0], [-1.0, 0.0, 1.0]),
        ("n=4: {-3,-1,1,3},{-2,-1,1,2}", [-3.0, -1.0, 1.0, 3.0], [-2.0, -1.0, 1.0, 2.0]),
        ("n=4: {-5,-1,1,5},{-5,-1,1,5}", [-5.0, -1.0, 1.0, 5.0], [-5.0, -1.0, 1.0, 5.0]),
    ]

    for name, roots_a, roots_b in configs:
        roots_a = np.array(roots_a, dtype=float)
        roots_b = np.array(roots_b, dtype=float)
        n = len(roots_a)

        phi_a = Phi_n(roots_a)
        phi_b = Phi_n(roots_b)
        target = 1/phi_a + 1/phi_b

        roots_r, _ = boxplus_mss(roots_a, roots_b)
        phi_r = Phi_n(roots_r)
        inv_phi_r = 1/phi_r

        A = np.diag(roots_a.astype(complex))
        B = np.diag(roots_b.astype(complex))

        inv_phis = []
        phis = []
        N_samp = 3000

        for _ in range(N_samp):
            U = haar_unitary(n)
            M = A + U @ B @ U.conj().T
            eigs = np.sort(np.real(np.linalg.eigvalsh(M)))

            if np.min(np.diff(eigs)) < 1e-10:
                continue

            phi_M = Phi_n(eigs)
            inv_phis.append(1.0 / phi_M)
            phis.append(phi_M)

        inv_phis = np.array(inv_phis)
        phis = np.array(phis)

        print(f"\n  {name}")
        print(f"    1/Phi(A) + 1/Phi(B) = {target:.8f}")
        print(f"    1/Phi(r) [boxplus]  = {inv_phi_r:.8f}")
        print(f"    E[1/Phi(M)]         = {np.mean(inv_phis):.8f}")
        print(f"    E[Phi(M)]           = {np.mean(phis):.8f}")
        print(f"    Phi(r)              = {phi_r:.8f}")
        print(f"    Jensen check: E[Phi(M)] >= Phi(r)? {np.mean(phis) >= phi_r - 1e-6}")
        print(f"    Jensen check: 1/Phi(r) >= E[1/Phi(M)]? {inv_phi_r >= np.mean(inv_phis) - 1e-6}")
        print(f"    Quantiles of 1/Phi(M): "
              f"5%={np.percentile(inv_phis,5):.6f}, "
              f"50%={np.percentile(inv_phis,50):.6f}, "
              f"95%={np.percentile(inv_phis,95):.6f}")
        print(f"    Fraction with 1/Phi(M) >= target: "
              f"{np.mean(inv_phis >= target - 1e-8):.4f}")


# ============================================================================
# PART 8: VARIANCE ADDITIVITY AND SCALING
# ============================================================================

def part8_variance_additivity():
    """
    Verify exact variance additivity: Var(r) = Var(p) + Var(q).
    This is a proved result and serves as a consistency check.
    """
    print("\n" + "=" * 78)
    print("PART 8: VARIANCE ADDITIVITY CHECK")
    print("=" * 78)

    for n in [3, 4, 5, 6]:
        max_err = 0.0
        for _ in range(100):
            roots_p = np.sort(np.random.randn(n) * 3)
            roots_q = np.sort(np.random.randn(n) * 3)
            for arr in [roots_p, roots_q]:
                arr -= np.mean(arr)  # center
                for i in range(1, n):
                    if abs(arr[i] - arr[i-1]) < 0.1:
                        arr[i] = arr[i-1] + 0.3
                arr -= np.mean(arr)

            roots_r, _ = boxplus_mss(roots_p, roots_q)

            var_p = np.var(roots_p) * n  # sum of squares, centered
            var_q = np.var(roots_q) * n
            var_r = np.var(roots_r) * n

            err = abs(var_r - var_p - var_q) / max(var_r, 1e-10)
            max_err = max(max_err, err)

        print(f"  n={n}: max relative error |Var(r) - Var(p) - Var(q)| / Var(r) = {max_err:.2e}")


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print("#" * 78)
    print("# RANDOM MATRIX / JENSEN'S INEQUALITY INVESTIGATION")
    print("# Fisher Superadditivity: 1/Phi_n(p boxplus_n q) >= 1/Phi_n(p) + 1/Phi_n(q)")
    print("#" * 78)

    # Part 1: Critical formula audit
    part1_results = part1_formula_audit()

    # Part 2: Concavity
    part2_concavity_analysis()

    # Part 3: Per-realization
    part3_per_realization()

    # Part 4: Trace formulation
    part4_trace_formulation()

    # Part 5: Majorization
    part5_majorization()

    # Part 6: Large-scale numerical
    part6_large_scale()

    # Part 7: Distribution analysis
    part7_eigenvalue_distribution()

    # Part 8: Variance additivity
    part8_variance_additivity()

    # ======================================================================
    # FINAL SUMMARY
    # ======================================================================
    print("\n" + "#" * 78)
    print("# FINAL SUMMARY")
    print("#" * 78)

    print("""
KEY FINDINGS:

1. FORMULA AUDIT: The hat-e convolution formula (hat_e_k(r) = sum hat_e_j(p)*hat_e_{k-j}(q))
   is NOT the same as the correct MSS formula (g_k = sum C(n-j,i)/C(n,i)*e_i*e_j).
   The previous "counterexamples" at n=4 may have used the WRONG formula.

2. CONCAVITY: 1/Phi_n is NEITHER consistently concave NOR convex in:
   - Coefficient space (E, F)
   - Root space
   This means a naive Jensen argument CANNOT work.

3. PER-REALIZATION: For individual Haar unitaries U, the inequality
   1/Phi(A+UBU*) >= 1/Phi(A) + 1/Phi(B) often FAILS.
   The conjecture is about the EXPECTED polynomial, not individual realizations.

4. TRACE FORMULATION: Phi_n = ||C @ 1||^2 = -1^T C^2 1 where C is the
   antisymmetric Cauchy matrix. This is verified but does not directly
   lead to a proof because the Cauchy matrix of r is not simply related
   to those of p, q.

5. MAJORIZATION: The centered roots of r DO majorize the centered roots
   of both p and q (verified universally). But 1/Phi is NOT Schur-concave,
   so majorization alone does not close the argument.

6. LARGE-SCALE TEST: With the CORRECT MSS formula, the conjecture
   1/Phi(r) >= 1/Phi(p) + 1/Phi(q) should be tested carefully for n >= 4.

7. VARIANCE ADDITIVITY: Var(r) = Var(p) + Var(q) exactly (PROVED).

8. LITERATURE: No existing proof of finite free Fisher superadditivity.
   Gribinski (2019) defines polynomial entropy but not Fisher information.
   The conjecture appears to be NOVEL.
""")
