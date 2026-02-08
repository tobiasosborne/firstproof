"""
investigate_schur_convexity.py -- Schur convexity / majorization approach to Fisher superadditivity.

Investigates whether majorization theory and Schur convexity/concavity can prove:
    1/Phi_n(p boxplus_n q) >= 1/Phi_n(p) + 1/Phi_n(q)

Key questions:
  Q1: Does the gap vector of r = p boxplus q majorize that of p (or vice versa)?
  Q2: Is 1/Phi_n Schur-concave in the gap vector?
  Q3: Is Phi_n Schur-convex in the gap vector?
  Q4: Is there a variable change making Phi_n or 1/Phi_n Schur-convex/concave?
  Q5: Does the MSS convolution produce root vectors related by majorization?
  Q6: Can the derivative identity r' = n*(p^{(1)} boxplus_{n-1} q^{(1)}) be used?

Uses the CORRECT MSS boxplus formula (Formula A = Formula B, verified):
  e_k(r) = sum_{i+j=k} C(n-j,i)/C(n,i) * e_i(p) * e_j(q)

Author: prover agent (Schur convexity investigation)
Date: 2026-02-08
"""

import numpy as np
from math import comb, factorial
from itertools import combinations
import sys

np.random.seed(2026)

# =====================================================================
# CORE FUNCTIONS
# =====================================================================

def elem_sym(roots, k):
    """Elementary symmetric polynomial e_k of the given roots."""
    n = len(roots)
    if k == 0:
        return 1.0
    if k > n:
        return 0.0
    return sum(np.prod(combo) for combo in combinations(roots, k))


def boxplus(roots_p, roots_q):
    """
    Correct MSS finite free additive convolution.
    Formula: e_k(r) = sum_{i+j=k} C(n-j,i)/C(n,i) * e_i(p) * e_j(q)
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
    poly_r = np.zeros(n + 1)
    for k in range(n + 1):
        poly_r[k] = (-1) ** k * g[k]

    roots_r = np.sort(np.real(np.roots(poly_r)))
    return roots_r


def H_values(roots):
    """H_p(lambda_i) = sum_{j != i} 1/(lambda_i - lambda_j)"""
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H


def Phi_n(roots):
    """Phi_n(p) = sum_i H_p(lambda_i)^2"""
    H = H_values(roots)
    return np.sum(H**2)


def gap_vector(roots):
    """Sorted gap vector: (lambda_2 - lambda_1, ..., lambda_n - lambda_{n-1})"""
    r = np.sort(roots)
    return np.diff(r)


def majorizes(x, y, tol=1e-10):
    """
    Check if x majorizes y (x >- y).
    Requires: sum of k largest of x >= sum of k largest of y, for all k.
    And: total sums are equal.
    """
    if len(x) != len(y):
        return False
    n = len(x)
    xs = np.sort(x)[::-1]  # descending
    ys = np.sort(y)[::-1]
    # Check total sum
    if abs(np.sum(xs) - np.sum(ys)) > tol * max(1, abs(np.sum(xs))):
        return False
    # Check partial sums
    for k in range(1, n + 1):
        if np.sum(xs[:k]) < np.sum(ys[:k]) - tol:
            return False
    return True


def weakly_majorizes(x, y, tol=1e-10):
    """
    Check if x weakly majorizes y (x >-_w y).
    Only requires: sum of k largest of x >= sum of k largest of y, for all k.
    (No total sum equality needed.)
    """
    if len(x) != len(y):
        return False
    n = len(x)
    xs = np.sort(x)[::-1]
    ys = np.sort(y)[::-1]
    for k in range(1, n + 1):
        if np.sum(xs[:k]) < np.sum(ys[:k]) - tol:
            return False
    return True


def random_roots(n, spread=3.0, min_gap=0.3):
    """Generate n random sorted roots with minimum gap."""
    roots = np.sort(np.random.randn(n) * spread)
    for i in range(1, n):
        if roots[i] - roots[i - 1] < min_gap:
            roots[i] = roots[i - 1] + min_gap + np.random.rand() * 0.3
    return roots


# =====================================================================
# Q1: MAJORIZATION OF GAP VECTORS UNDER MSS CONVOLUTION
# =====================================================================

def investigate_gap_majorization():
    """
    Does the gap vector of r = p boxplus q majorize that of p or q (or vice versa)?
    If g_r >- g_p, combined with Schur concavity of 1/Phi, this could work.
    """
    print("=" * 70)
    print("Q1: MAJORIZATION OF GAP VECTORS UNDER MSS CONVOLUTION")
    print("=" * 70)
    print()

    results = {
        'gr_maj_gp': 0,  # g_r majorizes g_p
        'gp_maj_gr': 0,  # g_p majorizes g_r
        'gr_wmaj_gp': 0,  # g_r weakly majorizes g_p
        'gp_wmaj_gr': 0,  # g_p weakly majorizes g_r
        'gr_maj_gq': 0,
        'gq_maj_gr': 0,
        'neither': 0,
        'total': 0,
    }

    for trial in range(500):
        for n in [3, 4, 5, 6]:
            roots_p = random_roots(n)
            roots_q = random_roots(n)

            try:
                roots_r = boxplus(roots_p, roots_q)
                gaps_r = gap_vector(roots_r)
                gaps_p = gap_vector(roots_p)
                gaps_q = gap_vector(roots_q)

                if np.any(gaps_r < 0.01) or np.any(gaps_p < 0.01) or np.any(gaps_q < 0.01):
                    continue

                results['total'] += 1

                if majorizes(gaps_r, gaps_p):
                    results['gr_maj_gp'] += 1
                if majorizes(gaps_p, gaps_r):
                    results['gp_maj_gr'] += 1
                if majorizes(gaps_r, gaps_q):
                    results['gr_maj_gq'] += 1
                if majorizes(gaps_q, gaps_r):
                    results['gq_maj_gr'] += 1
                if weakly_majorizes(gaps_r, gaps_p):
                    results['gr_wmaj_gp'] += 1
                if weakly_majorizes(gaps_p, gaps_r):
                    results['gp_wmaj_gr'] += 1

                if (not majorizes(gaps_r, gaps_p) and not majorizes(gaps_p, gaps_r) and
                    not majorizes(gaps_r, gaps_q) and not majorizes(gaps_q, gaps_r)):
                    results['neither'] += 1

            except Exception:
                continue

    total = results['total']
    print(f"Total valid trials: {total}")
    print(f"  g_r majorizes g_p:   {results['gr_maj_gp']}/{total} = {results['gr_maj_gp']/total:.3f}")
    print(f"  g_p majorizes g_r:   {results['gp_maj_gr']}/{total} = {results['gp_maj_gr']/total:.3f}")
    print(f"  g_r majorizes g_q:   {results['gr_maj_gq']}/{total} = {results['gr_maj_gq']/total:.3f}")
    print(f"  g_q majorizes g_r:   {results['gq_maj_gr']}/{total} = {results['gq_maj_gr']/total:.3f}")
    print(f"  g_r w-majorizes g_p: {results['gr_wmaj_gp']}/{total} = {results['gr_wmaj_gp']/total:.3f}")
    print(f"  g_p w-majorizes g_r: {results['gp_wmaj_gr']}/{total} = {results['gp_wmaj_gr']/total:.3f}")
    print(f"  Neither direction:   {results['neither']}/{total} = {results['neither']/total:.3f}")
    print()

    # Note: gap vectors have different sums in general, so exact majorization
    # requires equal sums. Let's also check NORMALIZED gap vectors.
    print("--- NORMALIZED gap vectors (sum to 1) ---")
    norm_results = {
        'gr_maj_gp': 0, 'gp_maj_gr': 0,
        'gr_maj_gq': 0, 'gq_maj_gr': 0,
        'neither': 0, 'total': 0,
    }

    for trial in range(500):
        for n in [3, 4, 5, 6]:
            roots_p = random_roots(n)
            roots_q = random_roots(n)

            try:
                roots_r = boxplus(roots_p, roots_q)
                gaps_r = gap_vector(roots_r)
                gaps_p = gap_vector(roots_p)
                gaps_q = gap_vector(roots_q)

                if np.any(gaps_r < 0.01) or np.any(gaps_p < 0.01) or np.any(gaps_q < 0.01):
                    continue

                # Normalize to sum 1
                ng_r = gaps_r / np.sum(gaps_r)
                ng_p = gaps_p / np.sum(gaps_p)
                ng_q = gaps_q / np.sum(gaps_q)

                norm_results['total'] += 1

                if majorizes(ng_r, ng_p):
                    norm_results['gr_maj_gp'] += 1
                if majorizes(ng_p, ng_r):
                    norm_results['gp_maj_gr'] += 1
                if majorizes(ng_r, ng_q):
                    norm_results['gr_maj_gq'] += 1
                if majorizes(ng_q, ng_r):
                    norm_results['gq_maj_gr'] += 1

                if (not majorizes(ng_r, ng_p) and not majorizes(ng_p, ng_r) and
                    not majorizes(ng_r, ng_q) and not majorizes(ng_q, ng_r)):
                    norm_results['neither'] += 1

            except Exception:
                continue

    total = norm_results['total']
    print(f"Total valid trials: {total}")
    print(f"  ng_r majorizes ng_p: {norm_results['gr_maj_gp']}/{total} = {norm_results['gr_maj_gp']/total:.3f}")
    print(f"  ng_p majorizes ng_r: {norm_results['gp_maj_gr']}/{total} = {norm_results['gp_maj_gr']/total:.3f}")
    print(f"  ng_r majorizes ng_q: {norm_results['gr_maj_gq']}/{total} = {norm_results['gr_maj_gq']/total:.3f}")
    print(f"  ng_q majorizes ng_r: {norm_results['gq_maj_gr']}/{total} = {norm_results['gq_maj_gr']/total:.3f}")
    print(f"  Neither:             {norm_results['neither']}/{total} = {norm_results['neither']/total:.3f}")
    print()

    return results, norm_results


# =====================================================================
# Q2: IS 1/Phi_n SCHUR-CONCAVE IN THE GAP VECTOR?
# =====================================================================

def investigate_schur_concavity_of_inv_phi():
    """
    Fix the mean and total span. Vary the gaps subject to majorization.
    If x >- y (x majorizes y), check whether 1/Phi(x) <= 1/Phi(y) (Schur-concave)
    or 1/Phi(x) >= 1/Phi(y) (Schur-convex).

    We parametrize: roots from gap vector g = (g_1, ..., g_{n-1}) with sum(g) = S fixed.
    Construct roots as: lambda_k = start + sum_{j<k} g_j
    """
    print("=" * 70)
    print("Q2: IS 1/Phi_n SCHUR-CONCAVE IN THE GAP VECTOR?")
    print("=" * 70)
    print()

    def roots_from_gaps(gaps, start=0.0):
        """Construct sorted roots from gap vector, starting at `start`."""
        n = len(gaps) + 1
        roots = np.zeros(n)
        roots[0] = start
        for i in range(len(gaps)):
            roots[i + 1] = roots[i] + gaps[i]
        return roots

    def test_schur_property(n, trials=2000):
        """For degree n polynomials, test Schur property of 1/Phi as function of gaps."""
        schur_concave_violations = 0
        schur_convex_violations = 0
        total_comparable = 0

        for _ in range(trials):
            # Generate two gap vectors where one majorizes the other
            # Method: start with uniform gaps, apply a T-transform to create majorized pair
            S = 5.0 + np.random.rand() * 5.0  # total span
            m = n - 1  # number of gaps

            # Start with random gap vector
            g1 = np.random.dirichlet(np.ones(m)) * S
            g1 = np.sort(g1)[::-1]  # sort descending for majorization check

            # Create g2 by a Robin Hood (T-transform) operation:
            # move mass from a large gap to a small gap
            g2 = g1.copy()
            if m >= 2:
                i, j = 0, m - 1  # largest and smallest
                transfer = np.random.rand() * min(g2[i] - g2[j], g2[i]) * 0.4
                g2[i] -= transfer
                g2[j] += transfer
                g2 = np.sort(g2)[::-1]

            # Now g1 should majorize g2 (since we transferred from larger to smaller)
            if not majorizes(g1, g2, tol=1e-8):
                continue

            total_comparable += 1

            # Build roots and compute Phi
            roots1 = roots_from_gaps(np.sort(g1), start=0.0)
            roots2 = roots_from_gaps(np.sort(g2), start=0.0)

            try:
                phi1 = Phi_n(roots1)
                phi2 = Phi_n(roots2)
                inv_phi1 = 1.0 / phi1
                inv_phi2 = 1.0 / phi2

                # g1 >- g2 (g1 majorizes g2, g1 is "more spread")
                # Schur-concave: g1 >- g2 => f(g1) <= f(g2)
                # Schur-convex: g1 >- g2 => f(g1) >= f(g2)

                if inv_phi1 > inv_phi2 + 1e-12:
                    schur_concave_violations += 1
                if inv_phi1 < inv_phi2 - 1e-12:
                    schur_convex_violations += 1

            except Exception:
                total_comparable -= 1
                continue

        return total_comparable, schur_concave_violations, schur_convex_violations

    for n in [3, 4, 5, 6, 7]:
        total, sc_viol, sx_viol = test_schur_property(n, trials=3000)
        print(f"  n={n}: {total} comparable pairs")
        print(f"    Schur-concave violations of 1/Phi: {sc_viol}/{total}", end="")
        if total > 0:
            print(f" ({sc_viol/total:.3f})")
        else:
            print()
        print(f"    Schur-convex violations of 1/Phi:  {sx_viol}/{total}", end="")
        if total > 0:
            print(f" ({sx_viol/total:.3f})")
        else:
            print()

        if sc_viol == 0 and total > 0:
            print(f"    => 1/Phi_n appears SCHUR-CONCAVE in gaps!")
        elif sx_viol == 0 and total > 0:
            print(f"    => 1/Phi_n appears SCHUR-CONVEX in gaps!")
        else:
            print(f"    => 1/Phi_n is NEITHER Schur-convex nor Schur-concave in gaps")
        print()


# =====================================================================
# Q3: IS Phi_n SCHUR-CONVEX IN THE GAP VECTOR?
# =====================================================================

def investigate_schur_convexity_of_phi():
    """
    Check if Phi_n is Schur-convex in the gap vector.
    If g1 majorizes g2 (g1 is more spread/unequal), is Phi(g1) >= Phi(g2)?
    """
    print("=" * 70)
    print("Q3: IS Phi_n SCHUR-CONVEX IN THE GAP VECTOR?")
    print("=" * 70)
    print()

    def roots_from_gaps(gaps, start=0.0):
        n = len(gaps) + 1
        roots = np.zeros(n)
        roots[0] = start
        for i in range(len(gaps)):
            roots[i + 1] = roots[i] + gaps[i]
        return roots

    for n in [3, 4, 5, 6]:
        sx_viol = 0
        sc_viol = 0
        total = 0

        for _ in range(3000):
            S = 5.0 + np.random.rand() * 5.0
            m = n - 1
            g1 = np.random.dirichlet(np.ones(m)) * S
            g1 = np.sort(g1)[::-1]

            g2 = g1.copy()
            if m >= 2:
                i, j = 0, m - 1
                transfer = np.random.rand() * min(g2[i] - g2[j], g2[i]) * 0.4
                g2[i] -= transfer
                g2[j] += transfer
                g2 = np.sort(g2)[::-1]

            if not majorizes(g1, g2, tol=1e-8):
                continue

            total += 1
            roots1 = roots_from_gaps(np.sort(g1))
            roots2 = roots_from_gaps(np.sort(g2))

            try:
                phi1 = Phi_n(roots1)
                phi2 = Phi_n(roots2)

                # g1 >- g2
                # Schur-convex: phi1 >= phi2
                if phi1 < phi2 - 1e-12:
                    sx_viol += 1
                if phi1 > phi2 + 1e-12:
                    sc_viol += 1
            except:
                total -= 1

        print(f"  n={n}: {total} comparable pairs")
        print(f"    Phi Schur-convex violations:  {sx_viol}/{total}")
        print(f"    Phi Schur-concave violations: {sc_viol}/{total}")
        if sx_viol == 0 and total > 0:
            print(f"    => Phi_n appears SCHUR-CONVEX in gaps (more unequal gaps => larger Phi)")
        elif sc_viol == 0 and total > 0:
            print(f"    => Phi_n appears SCHUR-CONCAVE in gaps")
        else:
            print(f"    => NEITHER")
        print()


# =====================================================================
# Q4: ALTERNATIVE COORDINATES
# =====================================================================

def investigate_alternative_coordinates():
    """
    Try different parameterizations:
    - Log gaps: log(g_k)
    - Squared gaps: g_k^2
    - Reciprocal gaps: 1/g_k
    - Gap ratios
    Check Schur properties in these alternative coordinates.
    """
    print("=" * 70)
    print("Q4: SCHUR PROPERTIES IN ALTERNATIVE COORDINATES")
    print("=" * 70)
    print()

    def roots_from_gaps(gaps, start=0.0):
        n = len(gaps) + 1
        roots = np.zeros(n)
        roots[0] = start
        for i in range(len(gaps)):
            roots[i + 1] = roots[i] + gaps[i]
        return roots

    n = 4
    m = n - 1

    # Test: is Phi_n Schur-convex in log(gaps)?
    # i.e., is Phi_n(g) a Schur-convex function of (log g_1, ..., log g_m)?
    # This means: log-majorization of gaps implies ordering of Phi.
    print(f"Testing n={n}, m={m} gaps")
    print()

    transforms = {
        'gaps': lambda g: g,
        'log_gaps': lambda g: np.log(g),
        'sq_gaps': lambda g: g**2,
        'inv_gaps': lambda g: 1.0/g,
    }

    for name, transform in transforms.items():
        sx_viol = 0
        sc_viol = 0
        total = 0

        for _ in range(3000):
            S = 5.0 + np.random.rand() * 5.0

            # Create gap pair where transformed vectors have majorization relation
            g1 = np.random.dirichlet(np.ones(m)) * S
            g1 = np.maximum(g1, 0.1)  # ensure positivity for log
            g1 = np.sort(g1)[::-1]

            g2 = g1.copy()
            if m >= 2:
                i, j = 0, m - 1
                transfer = np.random.rand() * min(g2[i] - g2[j], g2[i]) * 0.3
                g2[i] -= transfer
                g2[j] += transfer
                g2 = np.maximum(g2, 0.1)
                g2 = np.sort(g2)[::-1]

            t1 = transform(g1)
            t2 = transform(g2)

            if not majorizes(t1, t2, tol=1e-8):
                continue

            total += 1
            roots1 = roots_from_gaps(np.sort(g1))
            roots2 = roots_from_gaps(np.sort(g2))

            try:
                inv_phi1 = 1.0 / Phi_n(roots1)
                inv_phi2 = 1.0 / Phi_n(roots2)

                if inv_phi1 > inv_phi2 + 1e-12:
                    sc_viol += 1  # Schur-concave violation
                if inv_phi1 < inv_phi2 - 1e-12:
                    sx_viol += 1  # Schur-convex violation
            except:
                total -= 1

        print(f"  Coordinate: {name}")
        print(f"    {total} comparable pairs")
        if total > 0:
            print(f"    1/Phi Schur-concave violations: {sc_viol}/{total} ({sc_viol/total:.3f})")
            print(f"    1/Phi Schur-convex violations:  {sx_viol}/{total} ({sx_viol/total:.3f})")
            if sc_viol == 0:
                print(f"    => 1/Phi appears SCHUR-CONCAVE in {name}")
            elif sx_viol == 0:
                print(f"    => 1/Phi appears SCHUR-CONVEX in {name}")
            else:
                print(f"    => NEITHER")
        print()


# =====================================================================
# Q5: MAJORIZATION OF ROOT VECTORS UNDER MSS CONVOLUTION
# =====================================================================

def investigate_root_majorization():
    """
    Check whether the root vector of r = p boxplus q is related to those of p, q
    by majorization (after centering/normalizing).

    In the random matrix model A + UBU*, Horn's inequalities constrain the
    eigenvalues of the sum. Do they imply majorization?
    """
    print("=" * 70)
    print("Q5: MAJORIZATION OF ROOT VECTORS (not gaps)")
    print("=" * 70)
    print()

    results = {
        'r_maj_p': 0, 'p_maj_r': 0,
        'r_wmaj_p': 0, 'p_wmaj_r': 0,
        'total': 0,
    }

    for trial in range(500):
        for n in [3, 4, 5]:
            roots_p = random_roots(n)
            roots_q = random_roots(n)

            try:
                roots_r = boxplus(roots_p, roots_q)

                # Center: subtract means so sums are 0
                centered_p = roots_p - np.mean(roots_p)
                centered_q = roots_q - np.mean(roots_q)
                centered_r = roots_r - np.mean(roots_r)

                results['total'] += 1

                if majorizes(centered_r, centered_p):
                    results['r_maj_p'] += 1
                if majorizes(centered_p, centered_r):
                    results['p_maj_r'] += 1
                if weakly_majorizes(centered_r, centered_p):
                    results['r_wmaj_p'] += 1
                if weakly_majorizes(centered_p, centered_r):
                    results['p_wmaj_r'] += 1

            except Exception:
                continue

    total = results['total']
    print(f"Total valid trials: {total}")
    print(f"  centered_r majorizes centered_p:   {results['r_maj_p']}/{total} ({results['r_maj_p']/total:.3f})")
    print(f"  centered_p majorizes centered_r:   {results['p_maj_r']}/{total} ({results['p_maj_r']/total:.3f})")
    print(f"  centered_r w-majorizes centered_p: {results['r_wmaj_p']}/{total} ({results['r_wmaj_p']/total:.3f})")
    print(f"  centered_p w-majorizes centered_r: {results['p_wmaj_r']}/{total} ({results['p_wmaj_r']/total:.3f})")
    print()

    # Also check absolute values (|lambda_k|): does |roots_r| majorize |roots_p|?
    print("--- Absolute value majorization ---")
    abs_results = {'r_maj_p': 0, 'p_maj_r': 0, 'total': 0}

    for trial in range(500):
        for n in [3, 4, 5]:
            roots_p = random_roots(n) - 2.0  # shift to have mix of signs
            roots_q = random_roots(n) - 2.0

            try:
                roots_r = boxplus(roots_p, roots_q)

                abs_p = np.sort(np.abs(roots_p - np.mean(roots_p)))[::-1]
                abs_r = np.sort(np.abs(roots_r - np.mean(roots_r)))[::-1]

                abs_results['total'] += 1
                if weakly_majorizes(abs_r, abs_p):
                    abs_results['r_maj_p'] += 1
                if weakly_majorizes(abs_p, abs_r):
                    abs_results['p_maj_r'] += 1

            except Exception:
                continue

    total = abs_results['total']
    print(f"  |r| w-majorizes |p|: {abs_results['r_maj_p']}/{total} ({abs_results['r_maj_p']/total:.3f})")
    print(f"  |p| w-majorizes |r|: {abs_results['p_maj_r']}/{total} ({abs_results['p_maj_r']/total:.3f})")
    print()


# =====================================================================
# Q6: SQUARED DIFFERENCES MAJORIZATION
# =====================================================================

def investigate_squared_diff_majorization():
    """
    Check whether the vector of ALL pairwise squared differences
    D_ij = (lambda_i - lambda_j)^2 is related by majorization.
    This is motivated by the fact that Phi_n = sum_i (sum_j 1/D_ij)^2.
    """
    print("=" * 70)
    print("Q6: PAIRWISE SQUARED DIFFERENCES MAJORIZATION")
    print("=" * 70)
    print()

    def squared_diffs(roots):
        n = len(roots)
        diffs = []
        for i in range(n):
            for j in range(i + 1, n):
                diffs.append((roots[i] - roots[j]) ** 2)
        return np.sort(diffs)[::-1]

    results = {'r_wmaj_p': 0, 'p_wmaj_r': 0, 'total': 0}

    for trial in range(500):
        for n in [3, 4, 5]:
            roots_p = random_roots(n)
            roots_q = random_roots(n)

            try:
                roots_r = boxplus(roots_p, roots_q)

                sd_r = squared_diffs(roots_r)
                sd_p = squared_diffs(roots_p)
                sd_q = squared_diffs(roots_q)

                results['total'] += 1

                if weakly_majorizes(sd_r, sd_p):
                    results['r_wmaj_p'] += 1
                if weakly_majorizes(sd_p, sd_r):
                    results['p_wmaj_r'] += 1

            except Exception:
                continue

    total = results['total']
    print(f"Total valid trials: {total}")
    print(f"  sq_diffs(r) w-majorizes sq_diffs(p): {results['r_wmaj_p']}/{total} ({results['r_wmaj_p']/total:.3f})")
    print(f"  sq_diffs(p) w-majorizes sq_diffs(r): {results['p_wmaj_r']}/{total} ({results['p_wmaj_r']/total:.3f})")
    print()

    # Also check sum of squared diffs (= n * power sum p_2 of centered roots)
    print("--- Sum of squared differences ---")
    sum_larger = 0
    total2 = 0
    for trial in range(500):
        for n in [3, 4, 5]:
            roots_p = random_roots(n)
            roots_q = random_roots(n)
            try:
                roots_r = boxplus(roots_p, roots_q)
                s_p = np.sum(squared_diffs(roots_p))
                s_q = np.sum(squared_diffs(roots_q))
                s_r = np.sum(squared_diffs(roots_r))
                total2 += 1
                if s_r >= s_p + s_q - 1e-8:
                    sum_larger += 1
            except:
                continue

    print(f"  sum_sq_diffs(r) >= sum_sq_diffs(p) + sum_sq_diffs(q): "
          f"{sum_larger}/{total2} ({sum_larger/total2:.3f})")
    print("  (This is the Pythagorean property generalization)")
    print()


# =====================================================================
# Q7: DIRECT SCHUR ANALYSIS OF THE SUPERADDITIVITY
# =====================================================================

def investigate_direct_schur_bound():
    """
    For the specific case p = q (which is the hardest case for the inequality),
    check:
      1/Phi(p boxplus p) >= 2/Phi(p)

    In terms of gaps, if g is the gap vector of p, what is the gap vector of
    p boxplus p? And is 1/Phi(2*gap) >= 2/Phi(gap) by Schur arguments?
    """
    print("=" * 70)
    print("Q7: p=q CASE AND SCALING ANALYSIS")
    print("=" * 70)
    print()

    # For p = q with gaps g, what are the gaps of p boxplus p?
    print("Gap analysis for p boxplus p vs p:")

    for n in [3, 4, 5, 6]:
        print(f"\n  n = {n}:")
        ratios_data = []

        for trial in range(200):
            roots_p = random_roots(n)
            try:
                roots_r = boxplus(roots_p, roots_p)
                gaps_p = gap_vector(roots_p)
                gaps_r = gap_vector(roots_r)

                # Gap ratios
                ratios = gaps_r / gaps_p
                ratios_data.append(ratios)

            except Exception:
                continue

        if ratios_data:
            ratios_arr = np.array(ratios_data)
            print(f"    gap_r / gap_p ratios:")
            print(f"      mean: {np.mean(ratios_arr, axis=0)}")
            print(f"      min:  {np.min(ratios_arr, axis=0)}")
            print(f"      max:  {np.max(ratios_arr, axis=0)}")

            # For equally spaced roots, gap_r should be sqrt(2)*gap_p
            roots_eq = np.linspace(-3, 3, n)
            roots_eq_r = boxplus(roots_eq, roots_eq)
            gap_eq = gap_vector(roots_eq)
            gap_eq_r = gap_vector(roots_eq_r)
            print(f"    Equally spaced: gap_r/gap_p = {gap_eq_r / gap_eq}")
            print(f"    Expected sqrt(2) = {np.sqrt(2):.6f}")

    # Now test the superadditivity ratio for p = q
    print("\n\n  Superadditivity ratio for p = q:")
    print("  R = 1/Phi(p boxplus p) / (2/Phi(p)) = Phi(p) / (2*Phi(p boxplus p))")
    print()

    for n in [3, 4, 5, 6]:
        min_R = float('inf')
        max_R = 0
        total = 0
        for trial in range(500):
            roots_p = random_roots(n)
            try:
                roots_r = boxplus(roots_p, roots_p)
                phi_p = Phi_n(roots_p)
                phi_r = Phi_n(roots_r)
                R = phi_p / (2 * phi_r)
                min_R = min(min_R, R)
                max_R = max(max_R, R)
                total += 1
            except:
                continue

        print(f"  n={n}: R in [{min_R:.6f}, {max_R:.6f}] over {total} trials")
        if min_R >= 1 - 1e-8:
            print(f"    => 1/Phi(p boxplus p) >= 2/Phi(p) HOLDS (min_R = {min_R:.6f})")
        else:
            print(f"    => VIOLATED!")


# =====================================================================
# Q8: PARTIAL DERIVATIVE SIGNS (SCHUR CONDITION)
# =====================================================================

def investigate_schur_condition():
    """
    The Schur-Ostrowski criterion: f is Schur-convex iff
      (x_i - x_j)(df/dx_i - df/dx_j) >= 0  for all i,j.

    Test this numerically for Phi_n as a function of the gaps.
    """
    print("=" * 70)
    print("Q8: SCHUR-OSTROWSKI CRITERION FOR Phi_n(gaps)")
    print("=" * 70)
    print()

    def roots_from_gaps(gaps, start=0.0):
        n = len(gaps) + 1
        roots = np.zeros(n)
        roots[0] = start
        for i in range(len(gaps)):
            roots[i + 1] = roots[i] + gaps[i]
        return roots

    def phi_of_gaps(gaps):
        return Phi_n(roots_from_gaps(gaps))

    eps = 1e-6

    for n in [3, 4, 5]:
        m = n - 1
        violations = 0
        total = 0

        for trial in range(2000):
            gaps = np.random.rand(m) * 3 + 0.5

            # Compute partial derivatives numerically
            partials = np.zeros(m)
            phi0 = phi_of_gaps(gaps)
            for k in range(m):
                gaps_plus = gaps.copy()
                gaps_plus[k] += eps
                partials[k] = (phi_of_gaps(gaps_plus) - phi0) / eps

            # Check Schur-Ostrowski condition for Schur-CONVEX
            for i in range(m):
                for j in range(i + 1, m):
                    total += 1
                    val = (gaps[i] - gaps[j]) * (partials[i] - partials[j])
                    if val < -1e-4:  # violation of Schur-convexity
                        violations += 1

        print(f"  n={n}: Schur-Ostrowski condition for Phi_n Schur-convex:")
        print(f"    Violations: {violations}/{total} ({violations/total:.4f})")
        if violations == 0:
            print(f"    => Phi_n IS Schur-convex in gaps!")
        print()

        # Also check for 1/Phi_n Schur-concave
        violations_inv = 0
        total_inv = 0

        for trial in range(2000):
            gaps = np.random.rand(m) * 3 + 0.5

            partials_inv = np.zeros(m)
            inv_phi0 = 1.0 / phi_of_gaps(gaps)
            for k in range(m):
                gaps_plus = gaps.copy()
                gaps_plus[k] += eps
                partials_inv[k] = (1.0 / phi_of_gaps(gaps_plus) - inv_phi0) / eps

            for i in range(m):
                for j in range(i + 1, m):
                    total_inv += 1
                    val = (gaps[i] - gaps[j]) * (partials_inv[i] - partials_inv[j])
                    if val > 1e-4:  # violation of Schur-concavity
                        violations_inv += 1

        print(f"  n={n}: Schur-Ostrowski condition for 1/Phi_n Schur-concave:")
        print(f"    Violations: {violations_inv}/{total_inv} ({violations_inv/total_inv:.4f})")
        if violations_inv == 0:
            print(f"    => 1/Phi_n IS Schur-concave in gaps!")
        print()


# =====================================================================
# Q9: CAN WE COMBINE MAJORIZATION + SCHUR TO GET SUPERADDITIVITY?
# =====================================================================

def investigate_combined_approach():
    """
    The ideal proof path:
    1. Show g_r majorizes some combination of g_p, g_q
    2. Show 1/Phi is Schur-concave in the relevant coordinates
    3. Conclude 1/Phi_r >= f(1/Phi_p, 1/Phi_q)

    But superadditivity 1/Phi_r >= 1/Phi_p + 1/Phi_q involves TWO DIFFERENT
    polynomials p, q with DIFFERENT gap structures. We need to combine them.

    Let's check: if g1, g2 are gap vectors of p, q, and g_r is the gap vector
    of p boxplus q, is there a relation like:
      1/Phi(g_r) >= 1/Phi(g1) + 1/Phi(g2)?

    This is the original conjecture in gap-vector language.
    """
    print("=" * 70)
    print("Q9: COMBINING MAJORIZATION WITH SCHUR CONCAVITY")
    print("=" * 70)
    print()

    print("The challenge: even if 1/Phi_n is Schur-concave in gaps,")
    print("we need majorization BETWEEN g_r and some function of g_p, g_q.")
    print("But g_r, g_p, g_q have different total sums in general,")
    print("so ordinary majorization does not apply.")
    print()

    # Test: is there a function f(g_p, g_q) such that g_r majorizes f(g_p, g_q)?
    # Natural candidate: f(g_p, g_q) = sqrt(g_p^2 + g_q^2) (Pythagorean)
    print("Test: g_r vs element-wise sqrt(g_p^2 + g_q^2)")

    pyth_match = 0
    pyth_close = 0
    pyth_total = 0

    for trial in range(500):
        for n in [3, 4, 5]:
            roots_p = random_roots(n)
            roots_q = random_roots(n)

            try:
                roots_r = boxplus(roots_p, roots_q)
                gaps_p = gap_vector(roots_p)
                gaps_q = gap_vector(roots_q)
                gaps_r = gap_vector(roots_r)

                pyth = np.sqrt(gaps_p**2 + gaps_q**2)
                rel_err = np.max(np.abs(gaps_r - pyth) / gaps_r)

                pyth_total += 1
                if rel_err < 0.01:
                    pyth_match += 1
                if rel_err < 0.1:
                    pyth_close += 1

            except Exception:
                continue

    print(f"  Exact (1% rel err): {pyth_match}/{pyth_total} ({pyth_match/pyth_total:.3f})")
    print(f"  Close (10% rel err): {pyth_close}/{pyth_total} ({pyth_close/pyth_total:.3f})")
    print()

    # For equally spaced roots, the Pythagorean property holds exactly.
    # For non-equally spaced, it's approximate.
    # So the Pythagorean approach only works for equally spaced.

    # Alternative: check if Phi_r <= harmonic_mean(Phi_p, Phi_q)
    # i.e., 1/Phi_r >= 1/Phi_p + 1/Phi_q ... this is the original conjecture!
    print("Verifying the conjecture directly one more time:")
    viol = 0
    total = 0
    min_ratio = float('inf')

    for trial in range(1000):
        for n in [3, 4, 5, 6]:
            roots_p = random_roots(n)
            roots_q = random_roots(n)

            try:
                roots_r = boxplus(roots_p, roots_q)
                phi_p = Phi_n(roots_p)
                phi_q = Phi_n(roots_q)
                phi_r = Phi_n(roots_r)

                ratio = (1.0/phi_r) / (1.0/phi_p + 1.0/phi_q)
                total += 1
                min_ratio = min(min_ratio, ratio)
                if ratio < 1 - 1e-8:
                    viol += 1
            except:
                continue

    print(f"  {total} trials, {viol} violations, min ratio = {min_ratio:.8f}")
    if viol == 0:
        print(f"  => Conjecture holds with ratio >= {min_ratio:.8f}")
    print()


# =====================================================================
# Q10: DERIVATIVE IDENTITY AND INDUCTION
# =====================================================================

def investigate_derivative_identity():
    """
    The MSS derivative identity: (p boxplus_n q)' = n * (p^{(1)} boxplus_{n-1} q^{(1)})
    where p^{(1)} = p'/n is the "normalized derivative" with roots = critical points.

    If we can relate Phi_n to Phi_{n-1} via this identity, an inductive proof
    might work.

    Key: the roots of r' are the critical points of r, and they interlace with
    the roots of r.
    """
    print("=" * 70)
    print("Q10: DERIVATIVE IDENTITY AND INDUCTION")
    print("=" * 70)
    print()

    for n in [3, 4, 5]:
        print(f"\n  n = {n}:")

        for trial in range(5):
            roots_p = random_roots(n)
            roots_q = random_roots(n)

            try:
                roots_r = boxplus(roots_p, roots_q)

                # Derivative roots (critical points) of r
                r_poly = np.poly1d(np.poly(roots_r))
                r_deriv = r_poly.deriv()
                crit_r = np.sort(np.real(np.roots(r_deriv.coeffs)))

                # Normalized derivative: roots of p^{(1)} = critical points of p
                p_poly = np.poly1d(np.poly(roots_p))
                q_poly = np.poly1d(np.poly(roots_q))
                crit_p = np.sort(np.real(np.roots(p_poly.deriv().coeffs)))
                crit_q = np.sort(np.real(np.roots(q_poly.deriv().coeffs)))

                # p^{(1)} boxplus_{n-1} q^{(1)}
                roots_deriv_conv = boxplus(crit_p, crit_q)

                print(f"    Trial {trial}:")
                print(f"      Critical pts of r:  {np.array2string(crit_r, precision=4)}")
                print(f"      p'^{(1)} boxplus q'^{(1)}: {np.array2string(roots_deriv_conv, precision=4)}")
                print(f"      Match: {np.allclose(crit_r, roots_deriv_conv, atol=1e-4)}")

                # Check Phi relationship
                phi_r = Phi_n(roots_r)
                phi_crit_r = Phi_n(crit_r) if len(crit_r) > 1 else 0
                phi_p = Phi_n(roots_p)
                phi_crit_p = Phi_n(crit_p) if len(crit_p) > 1 else 0

                print(f"      Phi(roots_r) = {phi_r:.4f}")
                print(f"      Phi(crit_r)  = {phi_crit_r:.4f}")
                print(f"      Phi(crit_r)/Phi(roots_r) = {phi_crit_r/phi_r:.4f}")

            except Exception as e:
                print(f"    Trial {trial}: error - {e}")


# =====================================================================
# Q11: CONVEXITY OF Phi AS FUNCTION OF ROOTS DIRECTLY
# =====================================================================

def investigate_phi_convexity_in_roots():
    """
    Is Phi_n convex or concave as a function of the root vector?
    This would give different information than Schur properties.
    """
    print("=" * 70)
    print("Q11: CONVEXITY/CONCAVITY OF Phi_n IN ROOTS")
    print("=" * 70)
    print()

    for n in [3, 4, 5]:
        # Test convexity: Phi(t*x + (1-t)*y) <= t*Phi(x) + (1-t)*Phi(y)?
        convex_viol = 0
        concave_viol = 0
        total = 0

        for trial in range(2000):
            roots_x = random_roots(n)
            roots_y = random_roots(n)
            t = np.random.rand()

            roots_mid = t * roots_x + (1 - t) * roots_y
            # Check if mid roots are distinct
            gaps_mid = np.diff(np.sort(roots_mid))
            if np.any(gaps_mid < 0.01):
                continue

            try:
                phi_x = Phi_n(roots_x)
                phi_y = Phi_n(roots_y)
                phi_mid = Phi_n(np.sort(roots_mid))

                total += 1
                if phi_mid > t * phi_x + (1 - t) * phi_y + 1e-10:
                    convex_viol += 1
                if phi_mid < t * phi_x + (1 - t) * phi_y - 1e-10:
                    concave_viol += 1

            except:
                continue

        print(f"  n={n}: {total} trials")
        print(f"    Phi convex violations: {convex_viol}/{total}")
        print(f"    Phi concave violations: {concave_viol}/{total}")
        if convex_viol == 0:
            print(f"    => Phi_n appears CONVEX in roots!")
        elif concave_viol == 0:
            print(f"    => Phi_n appears CONCAVE in roots!")
        else:
            print(f"    => NEITHER convex nor concave")
        print()


# =====================================================================
# MAIN
# =====================================================================

if __name__ == "__main__":
    print("#" * 70)
    print("# SCHUR CONVEXITY / MAJORIZATION INVESTIGATION")
    print("# For Fisher superadditivity: 1/Phi(r) >= 1/Phi(p) + 1/Phi(q)")
    print("#" * 70)
    print()

    investigate_gap_majorization()
    investigate_schur_concavity_of_inv_phi()
    investigate_schur_convexity_of_phi()
    investigate_alternative_coordinates()
    investigate_root_majorization()
    investigate_squared_diff_majorization()
    investigate_direct_schur_bound()
    investigate_schur_condition()
    investigate_combined_approach()
    investigate_derivative_identity()
    investigate_phi_convexity_in_roots()

    print("\n" + "#" * 70)
    print("# INVESTIGATION COMPLETE")
    print("#" * 70)
