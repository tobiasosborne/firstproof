"""
verify_boxplus_definitive.py — DEFINITIVE verification of boxplus formulas.

Three agents produced conflicting results about the Fisher superadditivity conjecture:
  1/Phi_n(p boxplus_n q) >= 1/Phi_n(p) + 1/Phi_n(q)

The conflict centers on WHICH boxplus formula is correct.

This script:
  1. Implements ALL THREE candidate formulas
  2. Compares each against Monte Carlo Haar unitary simulation (ground truth)
  3. Tests the specific "counterexample" p = q = {-5, -1, 1, 5} (n=4)
  4. Tests additional cases (centered, non-centered, various n)
  5. Reports which formula is correct and whether the conjecture holds

FORMULAS:
  A) Source document (fisher_subordination_proof.md, Definition 0.2):
     c_k = sum_{i+j=k} [(n-i)!(n-j)!/(n!(n-k)!)] * a_i * b_j
     where a_k, b_k are polynomial COEFFICIENTS (a_0=b_0=1, a_k = coeff of x^{n-k})
     Note: a_k = (-1)^k * e_k(roots)

  B) Agent A's corrected formula (findings_induction.md):
     g_k = sum_{i+j=k} [C(n-j,i) / C(n,i)] * e_i(p) * e_j(q)
     where e_k are elementary symmetric polynomials of the roots.

  C) Agent B's formula (findings_matrix_info.md):
     hat_e_k(r) = sum_{j=0}^k hat_e_j(p) * hat_e_{k-j}(q)
     where hat_e_k = e_k / C(n,k)
     i.e. "finite free cumulants are additive"

GROUND TRUTH: Monte Carlo simulation of E[det(xI - A - UBU*)] with Haar U.
"""

import numpy as np
from math import factorial, comb
from itertools import combinations
import sys

np.random.seed(2026)

# =====================================================================
# HELPER FUNCTIONS
# =====================================================================

def elem_sym(roots, k):
    """Elementary symmetric polynomial e_k of the given roots."""
    n = len(roots)
    if k == 0:
        return 1.0
    if k > n:
        return 0.0
    return sum(np.prod(combo) for combo in combinations(roots, k))


def poly_coeffs_from_roots(roots):
    """
    Returns monic polynomial coefficients [1, a_1, a_2, ..., a_n]
    where p(x) = x^n + a_1*x^{n-1} + ... + a_n = prod(x - lambda_i).
    Note: a_k = (-1)^k * e_k(roots).
    """
    return np.poly(roots)  # numpy convention: leading coeff first


def Phi_n(roots):
    """
    Finite free Fisher information:
    Phi_n(p) = sum_i H_p(lambda_i)^2
    where H_p(lambda_i) = sum_{j != i} 1/(lambda_i - lambda_j)
    """
    n = len(roots)
    total = 0.0
    for i in range(n):
        H_i = sum(1.0 / (roots[i] - roots[j]) for j in range(n) if j != i)
        total += H_i ** 2
    return total


def haar_unitary(n):
    """
    Generate Haar-distributed random unitary matrix of size n.
    Method: QR decomposition of Ginibre ensemble with diagonal phase correction.
    """
    Z = (np.random.randn(n, n) + 1j * np.random.randn(n, n)) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    d = np.diag(R)
    ph = d / np.abs(d)
    U = Q @ np.diag(ph)
    return U


def monte_carlo_boxplus(roots_p, roots_q, N=50000):
    """
    GROUND TRUTH: E[det(xI - A - UBU*)] via Monte Carlo with Haar unitaries.
    Returns polynomial coefficients (numpy convention: leading=1) and roots.
    """
    n = len(roots_p)
    assert len(roots_q) == n
    A = np.diag(np.array(roots_p, dtype=complex))
    B = np.diag(np.array(roots_q, dtype=complex))

    poly_sum = np.zeros(n + 1, dtype=float)
    for _ in range(N):
        U = haar_unitary(n)
        M = A + U @ B @ U.conj().T
        coeffs = np.real(np.poly(M))
        poly_sum += coeffs

    poly_avg = poly_sum / N
    # Force monic
    poly_avg[0] = 1.0
    roots_r = np.sort(np.real(np.roots(poly_avg)))
    return roots_r, poly_avg


# =====================================================================
# FORMULA A: Source document (Definition 0.2)
# c_k = sum_{i+j=k} [(n-i)!(n-j)!/(n!(n-k)!)] * a_i * b_j
# where a_k = coeff of x^{n-k} in monic polynomial (so a_k = (-1)^k * e_k)
# =====================================================================

def boxplus_formula_A(roots_p, roots_q):
    """Source document formula using polynomial coefficients."""
    n = len(roots_p)
    poly_p = np.poly(roots_p)  # [1, a_1, ..., a_n]
    poly_q = np.poly(roots_q)  # [1, b_1, ..., b_n]

    # a_k = poly_p[k] (coefficient of x^{n-k} in the monic polynomial)
    a = poly_p  # a[0]=1, a[k] = coefficient of x^{n-k}
    b = poly_q

    c = np.zeros(n + 1)
    for k in range(n + 1):
        for i in range(k + 1):
            j = k - i
            if i <= n and j <= n:
                w = factorial(n - i) * factorial(n - j) / (factorial(n) * factorial(n - k))
                c[k] += w * a[i] * b[j]

    # Result polynomial: r(x) = sum_{k=0}^n c_k * x^{n-k}
    poly_r = c.copy()
    roots_r = np.sort(np.real(np.roots(poly_r)))
    return roots_r, poly_r


# =====================================================================
# FORMULA B: Agent A's corrected formula
# g_k = sum_{i+j=k} [C(n-j,i) / C(n,i)] * e_i(p) * e_j(q)
# Result: r(x) = sum_k (-1)^k g_k x^{n-k}
# =====================================================================

def boxplus_formula_B(roots_p, roots_q):
    """Agent A formula using elementary symmetric polynomials."""
    n = len(roots_p)

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
    poly_r = np.zeros(n + 1)
    for k in range(n + 1):
        poly_r[k] = (-1) ** k * g[k]

    roots_r = np.sort(np.real(np.roots(poly_r)))
    return roots_r, poly_r


# =====================================================================
# FORMULA C: Agent B's formula (simple additive finite free cumulants)
# hat_e_k(r) = sum_{j=0}^k hat_e_j(p) * hat_e_{k-j}(q)
# where hat_e_k = e_k / C(n,k)
# =====================================================================

def boxplus_formula_C(roots_p, roots_q):
    """Agent B formula: additive finite free cumulants (convolution of hat_e)."""
    n = len(roots_p)

    e_p = [elem_sym(roots_p, k) for k in range(n + 1)]
    e_q = [elem_sym(roots_q, k) for k in range(n + 1)]

    hat_p = [e_p[k] / comb(n, k) for k in range(n + 1)]
    hat_q = [e_q[k] / comb(n, k) for k in range(n + 1)]

    # Convolution: hat_e_k(r) = sum_{j=0}^k hat_e_j(p) * hat_e_{k-j}(q)
    hat_r = np.zeros(n + 1)
    for k in range(n + 1):
        for j in range(k + 1):
            hat_r[k] += hat_p[j] * hat_q[k - j]

    # Recover e_k(r) = C(n,k) * hat_r[k]
    e_r = [comb(n, k) * hat_r[k] for k in range(n + 1)]

    # r(x) = sum_k (-1)^k e_k(r) x^{n-k}
    poly_r = np.zeros(n + 1)
    for k in range(n + 1):
        poly_r[k] = (-1) ** k * e_r[k]

    roots_r = np.sort(np.real(np.roots(poly_r)))
    return roots_r, poly_r


# =====================================================================
# ALGEBRAIC EQUIVALENCE CHECK
# Before running Monte Carlo, let's check if formulas A and B are
# algebraically the same. Formula A uses a_k = (-1)^k * e_k.
# =====================================================================

def check_algebraic_equivalence():
    """
    Check whether formula A and formula B are the same algebraically.

    Formula A: c_k = sum_{i+j=k} [(n-i)!(n-j)!/(n!(n-k)!)] * a_i * b_j
               with a_i = (-1)^i * e_i, b_j = (-1)^j * e_j
               Result poly coeff[k] = c_k  (coefficient of x^{n-k})

    Formula B: g_k = sum_{i+j=k} [C(n-j,i)/C(n,i)] * e_i * e_j
               Result poly coeff[k] = (-1)^k * g_k

    For Formula A applied to the result:
      c_k = sum_{i+j=k} [(n-i)!(n-j)!/(n!(n-k)!)] * (-1)^i * e_i * (-1)^j * e_j
          = sum_{i+j=k} (-1)^{i+j} * [(n-i)!(n-j)!/(n!(n-k)!)] * e_i * e_j
          = (-1)^k * sum_{i+j=k} [(n-i)!(n-j)!/(n!(n-k)!)] * e_i * e_j

    For Formula B:
      poly_r[k] = (-1)^k * g_k = (-1)^k * sum_{i+j=k} [C(n-j,i)/C(n,i)] * e_i * e_j

    So A = B iff:
      (n-i)!(n-j)!/(n!(n-k)!) = C(n-j,i)/C(n,i)

    Let's check:
      C(n-j,i)/C(n,i) = [(n-j)!/(i!(n-j-i)!)] / [n!/(i!(n-i)!)]
                       = [(n-j)! * (n-i)!] / [(n-j-i)! * n!]
                       = [(n-j)! * (n-i)!] / [(n-k)! * n!]    (since j=k-i, so n-j-i = n-k)
                       = (n-i)!(n-j)! / (n!(n-k)!)

    YES! They are IDENTICAL. Formulas A and B are the same.
    """
    print("=" * 70)
    print("ALGEBRAIC EQUIVALENCE CHECK")
    print("=" * 70)
    print()
    print("Formula A weight: (n-i)!(n-j)! / (n!(n-k)!)")
    print("Formula B weight: C(n-j,i) / C(n,i)")
    print()
    print("Expanding C(n-j,i)/C(n,i):")
    print("  = [(n-j)!/(i!(n-j-i)!)] / [n!/(i!(n-i)!)]")
    print("  = (n-j)! * (n-i)! / ((n-j-i)! * n!)")
    print("  = (n-i)!(n-j)! / ((n-k)! * n!)  [since j=k-i => n-j-i=n-k]")
    print()
    print("RESULT: Formula A and Formula B are ALGEBRAICALLY IDENTICAL")
    print("        (after accounting for the sign convention a_k = (-1)^k e_k)")
    print()

    # Numerical verification
    for n in [2, 3, 4, 5]:
        for _ in range(5):
            roots_p = np.sort(np.random.randn(n) * 3)
            roots_q = np.sort(np.random.randn(n) * 3)
            _, poly_A = boxplus_formula_A(roots_p, roots_q)
            _, poly_B = boxplus_formula_B(roots_p, roots_q)
            err = np.max(np.abs(poly_A - poly_B))
            if err > 1e-10:
                print(f"  n={n}: MISMATCH! err={err:.2e}")
                print(f"    poly_A = {poly_A}")
                print(f"    poly_B = {poly_B}")
                return False
    print("  Numerical verification: A == B for all random tests (max err < 1e-10)")
    print()
    return True


# =====================================================================
# FORMULA C vs A: Are they different?
# =====================================================================

def check_C_vs_A():
    """Check if formula C differs from A/B."""
    print("=" * 70)
    print("FORMULA C vs A/B: Are they different?")
    print("=" * 70)
    print()

    # Centered case: p = {-a, a}, q = {-b, b}
    for (name, rp, rq) in [
        ("centered n=2: {-1,1},{-1,1}", [-1.0, 1.0], [-1.0, 1.0]),
        ("non-centered n=2: {0,1},{0,1}", [0.0, 1.0], [0.0, 1.0]),
        ("centered n=4: {-5,-1,1,5},{-5,-1,1,5}", [-5.0, -1.0, 1.0, 5.0], [-5.0, -1.0, 1.0, 5.0]),
        ("non-centered n=3: {1,2,4},{-1,0,3}", [1.0, 2.0, 4.0], [-1.0, 0.0, 3.0]),
    ]:
        rp = np.array(rp)
        rq = np.array(rq)
        _, poly_A = boxplus_formula_A(rp, rq)
        _, poly_C = boxplus_formula_C(rp, rq)
        err = np.max(np.abs(poly_A - poly_C))
        print(f"  {name}")
        print(f"    poly_A = {np.array2string(poly_A, precision=6)}")
        print(f"    poly_C = {np.array2string(poly_C, precision=6)}")
        print(f"    max diff = {err:.6e}")
        if err < 1e-10:
            print(f"    => SAME")
        else:
            print(f"    => DIFFERENT")
        print()


# =====================================================================
# MAIN VERIFICATION: Monte Carlo vs all formulas
# =====================================================================

def run_verification(roots_p, roots_q, name, N_mc=50000):
    """Run all formulas and Monte Carlo, report results."""
    n = len(roots_p)
    rp = np.array(roots_p, dtype=float)
    rq = np.array(roots_q, dtype=float)

    print(f"\n{'=' * 70}")
    print(f"TEST: {name}")
    print(f"  n = {n}")
    print(f"  p roots = {rp}")
    print(f"  q roots = {rq}")
    print(f"  Monte Carlo samples = {N_mc}")
    print(f"{'=' * 70}")

    # Formula A (Source document)
    roots_A, poly_A = boxplus_formula_A(rp, rq)
    # Formula B (Agent A = same as A algebraically)
    roots_B, poly_B = boxplus_formula_B(rp, rq)
    # Formula C (Agent B: additive cumulants)
    roots_C, poly_C = boxplus_formula_C(rp, rq)
    # Monte Carlo (Ground truth)
    roots_MC, poly_MC = monte_carlo_boxplus(rp, rq, N=N_mc)

    print(f"\n  Polynomial coefficients (monic, descending powers):")
    print(f"    Formula A (source doc):     {np.array2string(poly_A, precision=6)}")
    print(f"    Formula B (Agent A):        {np.array2string(poly_B, precision=6)}")
    print(f"    Formula C (Agent B cumul.): {np.array2string(poly_C, precision=6)}")
    print(f"    Monte Carlo:                {np.array2string(poly_MC, precision=6)}")

    print(f"\n  Roots of r = p boxplus q:")
    print(f"    Formula A: {np.array2string(roots_A, precision=6)}")
    print(f"    Formula B: {np.array2string(roots_B, precision=6)}")
    print(f"    Formula C: {np.array2string(roots_C, precision=6)}")
    print(f"    Monte Carlo: {np.array2string(roots_MC, precision=6)}")

    # Errors vs Monte Carlo
    err_A = np.max(np.abs(poly_A - poly_MC))
    err_B = np.max(np.abs(poly_B - poly_MC))
    err_C = np.max(np.abs(poly_C - poly_MC))
    err_AB = np.max(np.abs(poly_A - poly_B))

    print(f"\n  Max coefficient error vs Monte Carlo:")
    print(f"    Formula A: {err_A:.6e}")
    print(f"    Formula B: {err_B:.6e}")
    print(f"    Formula C: {err_C:.6e}")
    print(f"    A vs B:    {err_AB:.6e}")

    # Determine which is correct
    best_err = min(err_A, err_B, err_C)
    if err_A == best_err:
        best = "A (source doc) = B (Agent A)"
    elif err_C == best_err:
        best = "C (Agent B cumulants)"
    else:
        best = "B (Agent A)"

    print(f"\n  BEST MATCH: {best}")

    # Are A/B and C different?
    diff_AC = np.max(np.abs(poly_A - poly_C))
    print(f"  A vs C max diff: {diff_AC:.6e}")
    if diff_AC < 1e-10:
        print(f"  => All three formulas agree for this case.")
    else:
        print(f"  => Formulas A/B and C DISAGREE for this case.")
        if err_C > 10 * err_A:
            print(f"  => Formula C is WRONG (error {err_C:.4e} vs {err_A:.4e} for A/B)")
        elif err_A > 10 * err_C:
            print(f"  => Formula A/B is WRONG (error {err_A:.4e} vs {err_C:.4e} for C)")
        else:
            print(f"  => Errors are comparable; need more MC samples to distinguish")

    # Now check the conjecture with each formula
    print(f"\n  CONJECTURE CHECK: 1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q)")
    phi_p = Phi_n(rp)
    phi_q = Phi_n(rq)
    inv_p = 1.0 / phi_p
    inv_q = 1.0 / phi_q
    print(f"    Phi_n(p)  = {phi_p:.6f},  1/Phi_n(p) = {inv_p:.6f}")
    print(f"    Phi_n(q)  = {phi_q:.6f},  1/Phi_n(q) = {inv_q:.6f}")
    print(f"    1/Phi(p) + 1/Phi(q) = {inv_p + inv_q:.6f}")

    for label, roots_r in [("A/B", roots_A), ("C", roots_C), ("MC", roots_MC)]:
        phi_r = Phi_n(roots_r)
        inv_r = 1.0 / phi_r
        gap = inv_r - inv_p - inv_q
        status = "HOLDS" if gap >= -1e-10 else "VIOLATED"
        print(f"    Formula {label}: Phi_n(r) = {phi_r:.6f}, 1/Phi_n(r) = {inv_r:.6f}, "
              f"gap = {gap:.6f} [{status}]")

    return {
        'poly_A': poly_A, 'poly_B': poly_B, 'poly_C': poly_C, 'poly_MC': poly_MC,
        'roots_A': roots_A, 'roots_C': roots_C, 'roots_MC': roots_MC,
        'err_A': err_A, 'err_C': err_C, 'diff_AC': diff_AC,
    }


# =====================================================================
# LARGE-SCALE CONJECTURE TEST with correct formula
# =====================================================================

def large_scale_test(n_values, trials_per_n, N_mc_for_formula_check=30000):
    """
    Test the conjecture with the CORRECT formula over many random cases.
    First verify formula against MC, then use formula directly (fast).
    """
    print(f"\n{'#' * 70}")
    print(f"LARGE-SCALE CONJECTURE TEST")
    print(f"{'#' * 70}")

    for n in n_values:
        violations = 0
        min_gap = float('inf')
        max_gap = float('-inf')

        for trial in range(trials_per_n):
            # Random roots with distinct values
            roots_p = np.sort(np.random.randn(n) * 5)
            roots_q = np.sort(np.random.randn(n) * 5)

            # Ensure distinct roots
            for arr in [roots_p, roots_q]:
                for i in range(1, n):
                    if abs(arr[i] - arr[i - 1]) < 0.01:
                        arr[i] = arr[i - 1] + 0.01 + np.random.rand() * 0.5

            roots_r, _ = boxplus_formula_A(roots_p, roots_q)

            try:
                phi_p = Phi_n(roots_p)
                phi_q = Phi_n(roots_q)
                phi_r = Phi_n(roots_r)

                if phi_p < 1e-10 or phi_q < 1e-10 or phi_r < 1e-10:
                    continue

                gap = 1.0 / phi_r - 1.0 / phi_p - 1.0 / phi_q
                min_gap = min(min_gap, gap)
                max_gap = max(max_gap, gap)

                if gap < -1e-8:
                    violations += 1
                    if violations <= 3:
                        print(f"  n={n} trial {trial}: VIOLATION gap={gap:.6e}")
                        print(f"    p={roots_p}")
                        print(f"    q={roots_q}")
                        print(f"    r={roots_r}")
            except Exception as e:
                continue

        print(f"  n={n}: {trials_per_n} trials, {violations} violations, "
              f"min_gap={min_gap:.6e}, max_gap={max_gap:.6e}")

    # Now specifically test the alleged counterexample with MC verification
    print(f"\n  --- Specific alleged counterexample: p=q={{-5,-1,1,5}} ---")
    rp = np.array([-5.0, -1.0, 1.0, 5.0])
    rq = np.array([-5.0, -1.0, 1.0, 5.0])

    # Formula A/B
    roots_r_AB, _ = boxplus_formula_A(rp, rq)
    # Formula C
    roots_r_C, _ = boxplus_formula_C(rp, rq)
    # Monte Carlo
    roots_r_MC, _ = monte_carlo_boxplus(rp, rq, N=N_mc_for_formula_check)

    phi_p = Phi_n(rp)
    phi_q = Phi_n(rq)

    for label, rr in [("A/B (source+AgentA)", roots_r_AB),
                       ("C (Agent B cumulants)", roots_r_C),
                       ("Monte Carlo", roots_r_MC)]:
        phi_r = Phi_n(rr)
        gap = 1.0 / phi_r - 1.0 / phi_p - 1.0 / phi_q
        status = "HOLDS" if gap >= -1e-8 else "VIOLATED"
        print(f"    {label}: gap = {gap:.6e} [{status}]")


# =====================================================================
# MAIN
# =====================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("DEFINITIVE BOXPLUS FORMULA VERIFICATION")
    print("=" * 70)

    # Step 1: Algebraic equivalence check
    check_algebraic_equivalence()

    # Step 2: Check if C differs from A/B
    check_C_vs_A()

    # Step 3: The critical test — alleged counterexample p=q={-5,-1,1,5}
    print("\n" + "#" * 70)
    print("CRITICAL TEST: Alleged counterexample p = q = {-5, -1, 1, 5}")
    print("#" * 70)
    result_critical = run_verification(
        [-5.0, -1.0, 1.0, 5.0],
        [-5.0, -1.0, 1.0, 5.0],
        "COUNTEREXAMPLE: p=q={-5,-1,1,5}, n=4",
        N_mc=100000
    )

    # Step 4: Additional counterexample from Agent B
    result_counter2 = run_verification(
        [-3.0, -1.0, 1.0, 3.0],
        [-3.0, -1.0, 1.0, 3.0],
        "COUNTEREXAMPLE 2: p=q={-3,-1,1,3}, n=4",
        N_mc=100000
    )

    # Step 5: Various test cases
    print("\n" + "#" * 70)
    print("ADDITIONAL TEST CASES")
    print("#" * 70)

    test_cases = [
        ("centered n=2", [-1.0, 1.0], [-1.0, 1.0]),
        ("non-centered n=2", [0.0, 2.0], [1.0, 3.0]),
        ("centered n=3", [-2.0, 0.0, 2.0], [-1.0, 0.0, 1.0]),
        ("non-centered n=3", [1.0, 2.0, 4.0], [-1.0, 0.0, 3.0]),
        ("centered n=4 spread", [-7.0, -2.0, 2.0, 7.0], [-3.0, -1.0, 1.0, 3.0]),
        ("non-centered n=4", [0.0, 1.0, 3.0, 6.0], [-2.0, 0.0, 1.0, 4.0]),
    ]

    for name, rp, rq in test_cases:
        run_verification(rp, rq, name, N_mc=50000)

    # Step 6: Large-scale test
    large_scale_test(n_values=[2, 3, 4, 5, 6], trials_per_n=500)

    # =====================================================================
    # FINAL VERDICT
    # =====================================================================
    print("\n" + "=" * 70)
    print("FINAL VERDICT")
    print("=" * 70)
    print()

    diff_AC = result_critical['diff_AC']
    err_A = result_critical['err_A']
    err_C = result_critical['err_C']

    print(f"For the alleged counterexample p=q={{-5,-1,1,5}}:")
    print(f"  Formula A/B (source doc / Agent A) error vs MC: {err_A:.6e}")
    print(f"  Formula C (Agent B cumulants) error vs MC:      {err_C:.6e}")
    print(f"  Difference between A/B and C:                   {diff_AC:.6e}")
    print()

    if diff_AC < 1e-6:
        print("CONCLUSION: All formulas agree. The discrepancy is NOT about the formula.")
        # Check conjecture with the critical example
        rp = np.array([-5.0, -1.0, 1.0, 5.0])
        rq = np.array([-5.0, -1.0, 1.0, 5.0])
        roots_r, _ = boxplus_formula_A(rp, rq)
        phi_p = Phi_n(rp)
        phi_r = Phi_n(roots_r)
        gap = 1.0 / phi_r - 2.0 / phi_p
        if gap >= -1e-8:
            print(f"  Conjecture HOLDS for the alleged counterexample (gap={gap:.6e})")
        else:
            print(f"  Conjecture FAILS for the alleged counterexample (gap={gap:.6e})")
    elif err_A < err_C / 5:
        print("CONCLUSION: Formula A/B (source doc / Agent A) is CORRECT.")
        print("            Formula C (Agent B simple cumulants) is WRONG.")
        print()
        print("Agent B's 'counterexamples' used the WRONG boxplus formula.")
        print("With the CORRECT formula, checking conjecture...")
        rp = np.array([-5.0, -1.0, 1.0, 5.0])
        rq = np.array([-5.0, -1.0, 1.0, 5.0])
        roots_r, _ = boxplus_formula_A(rp, rq)
        phi_p = Phi_n(rp)
        phi_r = Phi_n(roots_r)
        gap = 1.0 / phi_r - 2.0 / phi_p
        if gap >= -1e-8:
            print(f"  Conjecture HOLDS with correct formula (gap={gap:.6e})")
            print()
            print("VERDICT: The 'refutation' at node 1.10.3 is INVALID.")
            print("         Counterexamples were computed with wrong boxplus formula.")
        else:
            print(f"  Conjecture STILL FAILS with correct formula (gap={gap:.6e})")
            print("  The refutation stands even with the correct formula.")
    elif err_C < err_A / 5:
        print("CONCLUSION: Formula C (Agent B cumulants) is CORRECT.")
        print("            Formula A/B (source doc) is WRONG.")
        print()
        print("Agent A's 'verification' used the WRONG boxplus formula.")
    else:
        print("CONCLUSION: Unable to definitively distinguish formulas.")
        print("            Need more Monte Carlo samples.")
