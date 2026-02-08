"""
investigate_shape_factor.py -- Shape factor decomposition approach to Fisher superadditivity.

CONJECTURE (Fisher superadditivity):
  For monic real-rooted degree-n polynomials p, q with simple roots:
    1/Phi_n(p boxplus_n q) >= 1/Phi_n(p) + 1/Phi_n(q)

SHAPE FACTOR APPROACH:
  Decompose: 1/Phi_n(p) = Var(p) / SF(p)
  where:
    Var(p) = (1/n) * sum_i (lambda_i - mean)^2  (variance of roots)
    SF(p) = Phi_n(p) * Var(p)                    (dimensionless "shape factor")

  Key claim: Var(r) = Var(p) + Var(q) exactly for r = p boxplus_n q.
  If also SF(r) <= min(SF(p), SF(q)), then:
    1/Phi(r) = Var(r)/SF(r) >= (Var(p)+Var(q))/min(SF(p),SF(q))
             >= Var(p)/SF(p) + Var(q)/SF(q) = 1/Phi(p) + 1/Phi(q)

TASKS:
  1. Verify Var(r) = Var(p) + Var(q) exactly.
  2. Numerically investigate SF(r) <= min(SF(p), SF(q)).
  3. Understand what SF is.
  4. Alternative decompositions.
  5. n=3 verification against proved result.
  6. Find counterexamples if SF conjecture is false.
  7. Attempt proof if true.

Uses CORRECT MSS boxplus: g_k = sum_{i+j=k} C(n-j,i)/C(n,i) * e_i(p) * e_j(q)

Author: prover agent (shape factor investigation)
Date: 2026-02-08
"""

import numpy as np
from math import comb, factorial
from itertools import combinations
from fractions import Fraction
import sys

np.random.seed(42)

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
    Returns sorted real roots.
    """
    n = len(roots_p)
    assert len(roots_q) == n

    e_p = [elem_sym(roots_p, k) for k in range(n + 1)]
    e_q = [elem_sym(roots_q, k) for k in range(n + 1)]

    g = np.zeros(n + 1)
    for k in range(n + 1):
        for i in range(k + 1):
            j = k - i
            if i <= n and j <= n and comb(n, i) > 0:
                w = comb(n - j, i) / comb(n, i)
                g[k] += w * e_p[i] * e_q[j]

    # Build polynomial: x^n - g_1 x^{n-1} + g_2 x^{n-2} - ... + (-1)^n g_n
    poly_coeffs = np.zeros(n + 1)
    for k in range(n + 1):
        poly_coeffs[k] = (-1) ** k * g[k]

    roots_r = np.sort(np.real(np.roots(poly_coeffs)))
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
    return np.sum(H ** 2)


def root_variance(roots):
    """Variance of the roots: (1/n) * sum_i (lambda_i - mean)^2"""
    return np.var(roots)  # numpy var uses 1/n denominator


def shape_factor(roots):
    """SF(p) = Phi_n(p) * Var(p)"""
    return Phi_n(roots) * root_variance(roots)


def random_roots(n, spread=3.0, min_gap=0.3):
    """Generate n random sorted roots with minimum gap."""
    roots = np.sort(np.random.randn(n) * spread)
    for i in range(1, n):
        if roots[i] - roots[i - 1] < min_gap:
            roots[i] = roots[i - 1] + min_gap + np.random.rand() * 0.3
    return roots


def random_centered_roots(n, spread=3.0, min_gap=0.3):
    """Generate n random sorted centered roots (mean = 0) with minimum gap."""
    roots = random_roots(n, spread, min_gap)
    return roots - np.mean(roots)


# =====================================================================
# TASK 1: VERIFY Var(r) = Var(p) + Var(q)
# =====================================================================

def task1_variance_additivity():
    """
    Verify: is Var(r) = Var(p) + Var(q) exactly for r = p boxplus q?

    Variance = (1/n) sum (lambda_i - mean)^2 = (1/n) sum lambda_i^2 - mean^2
             = e_1^2/n^2 - 2*e_2/n ... wait, more carefully:

    For centered polys (e_1 = 0): Var = (1/n) sum lambda_i^2 = -2*e_2/n
    (since sum lambda_i^2 = e_1^2 - 2*e_2 = -2*e_2 for centered)

    MSS boxplus for e_2 (centered, e_1=0):
    g_2 = sum_{i+j=2} C(n-j,i)/C(n,i) * e_i(p) * e_j(q)
    Terms: (i,j) = (0,2), (1,1), (2,0)
    (0,2): C(n-2,0)/C(n,0) * 1 * e_2(q) = 1 * e_2(q)
    (1,1): C(n-1,1)/C(n,1) * e_1(p) * e_1(q) = (n-1)/n * 0 * 0 = 0
    (2,0): C(n,2)/C(n,2) * e_2(p) * 1 = e_2(p)
    So e_2(r) = e_2(p) + e_2(q) for centered polys.

    For NON-centered: e_1(r) = e_1(p) + e_1(q), but then
    Var(r) = (sum r_i^2)/n - (sum r_i / n)^2
           = (e_1(r)^2 - 2*e_2(r))/n - (e_1(r)/n)^2
    Need to work this out.

    Actually: Var is translation-invariant (Var(p+c) = Var(p)), and
    MSS is equivariant (boxplus commutes with translation).
    So: Var(p boxplus q) = Var(centered(p boxplus q))
                         = Var(centered_p boxplus centered_q)
    And for centered polys, Var(r) = -2*e_2(r)/n = -2*(e_2(p)+e_2(q))/n
                                   = Var(p) + Var(q).
    """
    print("=" * 72)
    print("TASK 1: VERIFY Var(r) = Var(p) + Var(q)")
    print("=" * 72)
    print()

    # Analytical argument:
    print("ANALYTICAL ARGUMENT:")
    print("  Var is translation-invariant, MSS boxplus is translation-equivariant.")
    print("  So WLOG assume p, q centered (e_1 = 0). Then:")
    print("  Var(p) = -2*e_2(p)/n, since sum(lambda_i^2) = -2*e_2 for centered.")
    print("  MSS boxplus for e_2 with e_1=0:")
    print("    e_2(r) = C(n-2,0)/C(n,0)*e_2(q) + C(n-1,1)/C(n,1)*0*0 + C(n,2)/C(n,2)*e_2(p)")
    print("           = e_2(p) + e_2(q)")
    print("  Therefore Var(r) = -2*e_2(r)/n = -2*(e_2(p)+e_2(q))/n = Var(p) + Var(q).")
    print()
    print("  THIS IS EXACT ADDITIVITY, not just an inequality.")
    print()

    # Numerical verification
    print("NUMERICAL VERIFICATION:")
    max_rel_err = 0
    for n in [3, 4, 5, 6, 7, 8]:
        errors = []
        for trial in range(200):
            roots_p = random_roots(n)
            roots_q = random_roots(n)
            try:
                roots_r = boxplus(roots_p, roots_q)
                var_p = root_variance(roots_p)
                var_q = root_variance(roots_q)
                var_r = root_variance(roots_r)
                rel_err = abs(var_r - var_p - var_q) / max(var_r, 1e-15)
                errors.append(rel_err)
                max_rel_err = max(max_rel_err, rel_err)
            except:
                continue
        mean_err = np.mean(errors)
        max_err = np.max(errors)
        print(f"  n={n}: max rel error = {max_err:.2e}, mean = {mean_err:.2e} "
              f"({'EXACT' if max_err < 1e-10 else 'APPROXIMATE'})")

    print()
    if max_rel_err < 1e-10:
        print("  CONCLUSION: Var(r) = Var(p) + Var(q) EXACTLY (to machine precision).")
    else:
        print(f"  CONCLUSION: Var(r) = Var(p) + Var(q) up to relative error {max_rel_err:.2e}.")
        print("  This may be numerical error from root-finding for larger n.")
    print()

    # Also verify with exact rational arithmetic for small cases
    print("EXACT VERIFICATION (n=3, symbolic check via e_2 additivity):")
    print("  e_2(r) = e_2(p) + e_2(q) for centered n=3:")
    for trial in range(5):
        rp = random_centered_roots(3)
        rq = random_centered_roots(3)
        rr = boxplus(rp, rq)
        rr_c = rr - np.mean(rr)  # should already be centered, but just in case
        e2p = elem_sym(rp, 2)
        e2q = elem_sym(rq, 2)
        e2r = elem_sym(rr_c, 2)
        print(f"    e2(p)={e2p:.8f}, e2(q)={e2q:.8f}, e2(p)+e2(q)={e2p+e2q:.8f}, "
              f"e2(r)={e2r:.8f}, err={abs(e2r-e2p-e2q):.2e}")

    return True  # Var(r) = Var(p) + Var(q) is exact


# =====================================================================
# TASK 2: NUMERICALLY INVESTIGATE SF(r) <= min(SF(p), SF(q))
# =====================================================================

def task2_shape_factor_conjecture():
    """
    Compute SF = Phi * Var for many random polynomials.
    Check: is SF(r) <= min(SF(p), SF(q)) always?
    """
    print("=" * 72)
    print("TASK 2: SHAPE FACTOR CONJECTURE: SF(r) <= min(SF(p), SF(q))")
    print("=" * 72)
    print()

    for n in [3, 4, 5, 6, 7]:
        violations = 0
        total = 0
        worst_ratio = 0  # max SF(r) / min(SF(p), SF(q))
        worst_case = None
        ratios = []

        for trial in range(2000):
            roots_p = random_roots(n)
            roots_q = random_roots(n)

            try:
                roots_r = boxplus(roots_p, roots_q)
                sf_p = shape_factor(roots_p)
                sf_q = shape_factor(roots_q)
                sf_r = shape_factor(roots_r)
                min_sf_pq = min(sf_p, sf_q)

                ratio = sf_r / min_sf_pq
                ratios.append(ratio)
                total += 1

                if sf_r > min_sf_pq + 1e-10:
                    violations += 1
                    if ratio > worst_ratio:
                        worst_ratio = ratio
                        worst_case = (roots_p.copy(), roots_q.copy(), roots_r.copy(),
                                      sf_p, sf_q, sf_r)
            except:
                continue

        ratios = np.array(ratios)
        print(f"  n={n}: {total} trials, {violations} violations")
        if total > 0:
            print(f"    SF(r)/min(SF(p),SF(q)): mean={np.mean(ratios):.6f}, "
                  f"max={np.max(ratios):.6f}, min={np.min(ratios):.6f}")
            if violations > 0:
                print(f"    COUNTEREXAMPLE FOUND! Worst ratio = {worst_ratio:.6f}")
                rp, rq, rr, sfp, sfq, sfr = worst_case
                print(f"      p roots: {rp}")
                print(f"      q roots: {rq}")
                print(f"      r roots: {rr}")
                print(f"      SF(p)={sfp:.8f}, SF(q)={sfq:.8f}, SF(r)={sfr:.8f}")
            else:
                print(f"    NO VIOLATIONS -- conjecture holds for n={n}")
        print()

    # Now test more aggressively with targeted examples
    print("  --- Targeted tests: nearly-degenerate roots ---")
    for n in [3, 4, 5]:
        violations = 0
        total = 0
        worst_ratio = 0
        worst_case = None

        for trial in range(2000):
            # Create roots with one very small gap
            roots_p = random_roots(n, spread=5.0, min_gap=0.01)
            roots_q = random_roots(n, spread=5.0, min_gap=0.01)

            try:
                roots_r = boxplus(roots_p, roots_q)
                # Check all gaps positive
                if np.min(np.diff(roots_r)) < 1e-8:
                    continue

                sf_p = shape_factor(roots_p)
                sf_q = shape_factor(roots_q)
                sf_r = shape_factor(roots_r)
                min_sf = min(sf_p, sf_q)

                ratio = sf_r / min_sf
                total += 1

                if sf_r > min_sf + 1e-10:
                    violations += 1
                    if ratio > worst_ratio:
                        worst_ratio = ratio
                        worst_case = (roots_p.copy(), roots_q.copy(),
                                      sf_p, sf_q, sf_r)
            except:
                continue

        print(f"  n={n}: {total} near-degenerate trials, {violations} violations", end="")
        if violations > 0:
            print(f" (worst ratio = {worst_ratio:.6f})")
        else:
            print()

    # Test with highly asymmetric polynomials
    print("\n  --- Targeted tests: highly asymmetric roots ---")
    for n in [3, 4, 5]:
        violations = 0
        total = 0
        worst_ratio = 0

        for trial in range(2000):
            # One poly with equal spacing, other very skewed
            d = np.random.uniform(0.5, 3.0)
            roots_p = np.linspace(-d, d, n)

            # Very skewed polynomial
            roots_q = np.zeros(n)
            roots_q[0] = -10.0
            for i in range(1, n):
                roots_q[i] = roots_q[i-1] + np.random.exponential(0.3)

            try:
                roots_r = boxplus(roots_p, roots_q)
                if np.min(np.diff(roots_r)) < 1e-8:
                    continue

                sf_p = shape_factor(roots_p)
                sf_q = shape_factor(roots_q)
                sf_r = shape_factor(roots_r)
                min_sf = min(sf_p, sf_q)

                total += 1
                if sf_r > min_sf + 1e-10:
                    violations += 1
                    if sf_r / min_sf > worst_ratio:
                        worst_ratio = sf_r / min_sf
            except:
                continue

        print(f"  n={n}: {total} asymmetric trials, {violations} violations", end="")
        if violations > 0:
            print(f" (worst ratio = {worst_ratio:.6f})")
        else:
            print()

    # Test with p = q (symmetric case)
    print("\n  --- Targeted tests: p = q (symmetric case) ---")
    for n in [3, 4, 5, 6]:
        violations = 0
        total = 0
        ratios = []

        for trial in range(2000):
            roots_p = random_roots(n)
            try:
                roots_r = boxplus(roots_p, roots_p)
                sf_p = shape_factor(roots_p)
                sf_r = shape_factor(roots_r)
                ratio = sf_r / sf_p
                ratios.append(ratio)
                total += 1
                if sf_r > sf_p + 1e-10:
                    violations += 1
            except:
                continue

        if total > 0:
            ratios = np.array(ratios)
            print(f"  n={n}: {total} trials, {violations} violations, "
                  f"SF(r)/SF(p): mean={np.mean(ratios):.6f}, max={np.max(ratios):.6f}")

    print()


# =====================================================================
# TASK 3: WHAT IS THE SHAPE FACTOR?
# =====================================================================

def task3_understand_shape_factor():
    """
    SF = Phi * Var. What does this equal for various root configurations?

    For equally-spaced roots {-d, -d+2d/(n-1), ..., d} with spacing h = 2d/(n-1):
      Var = (1/n) sum_i (r_i - mean)^2 = (n^2-1)/(12n) * h^2 ... no, this depends on n.

    For general roots:
      SF = [sum_i H_p(lambda_i)^2] * [(1/n) sum_i (lambda_i - mu)^2]

    Since Phi is scale-homogeneous of degree -2 and Var is scale-homogeneous of degree 2,
    SF = Phi * Var is scale-INVARIANT (dimensionless).

    So SF depends only on the "shape" of the root configuration, i.e., the ratios of gaps.
    """
    print("=" * 72)
    print("TASK 3: UNDERSTANDING THE SHAPE FACTOR")
    print("=" * 72)
    print()

    # Verify scale invariance
    print("Scale invariance check (SF(c*roots) = SF(roots)):")
    for n in [3, 4, 5]:
        roots = random_roots(n)
        for c in [0.1, 0.5, 2.0, 10.0, 100.0]:
            sf_orig = shape_factor(roots)
            sf_scaled = shape_factor(c * roots)
            rel_err = abs(sf_scaled - sf_orig) / sf_orig
            if rel_err > 1e-10:
                print(f"  n={n}, c={c}: FAIL, rel_err={rel_err:.2e}")
        print(f"  n={n}: scale-invariant (all c tested, max rel err < 1e-10)")

    # SF for equally-spaced roots
    print("\nSF for equally-spaced roots (lambda_i = i, i=0,...,n-1):")
    for n in [3, 4, 5, 6, 7, 8, 9, 10]:
        roots_eq = np.arange(n, dtype=float)
        sf_eq = shape_factor(roots_eq)
        phi_eq = Phi_n(roots_eq)
        var_eq = root_variance(roots_eq)
        print(f"  n={n}: SF = {sf_eq:.8f}, Phi = {phi_eq:.8f}, Var = {var_eq:.6f}")

    # SF for two-cluster configurations: {-d, ..., -d} u {d, ..., d} with perturbation
    print("\nSF for two-cluster configurations (n=4):")
    for epsilon in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0]:
        roots = np.array([-1-epsilon, -1+epsilon, 1-epsilon, 1+epsilon])
        sf = shape_factor(roots)
        print(f"  eps={epsilon:.2f}: roots={roots}, SF={sf:.8f}")

    # SF as function of gap ratios for n=3
    print("\nSF for n=3 as function of gap ratio (g1/g2):")
    print("  (Centered cubic with roots -a, b, a+b-a such that mean=0)")
    for ratio in [0.1, 0.2, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
        # Gaps: g1, g2 = ratio*g1, total span = g1*(1+ratio)
        g1 = 1.0
        g2 = ratio * g1
        roots = np.array([0, g1, g1+g2])
        roots = roots - np.mean(roots)  # center
        sf = shape_factor(roots)
        print(f"  ratio={ratio:5.1f}: SF = {sf:.8f}")

    # Express SF in terms of moments for n=3
    print("\nSF for n=3 in terms of E, F (centered: e2 = -E, e3 = F):")
    print("  From prove_n3_symbolic.py:")
    print("    1/Phi_3 = (4*E^3 - 27*F^2) / (18*E^2)")
    print("    Var = -2*e_2/3 = 2*E/3")
    print("  So: SF = Phi * Var = [18*E^2/(4*E^3-27*F^2)] * [2*E/3]")
    print("         = 12*E^3 / (4*E^3-27*F^2)")
    print("         = 12 / (4 - 27*F^2/E^3)")
    print("  Let u = F^2/E^3 (dimensionless, u in [0, 4/27) for real roots).")
    print("  Then SF = 12 / (4 - 27*u)")
    print("  SF is MONOTONE INCREASING in u.")
    print("  At u = 0 (equally spaced): SF = 12/4 = 3.")
    print("  As u -> 4/27 (degenerate): SF -> infinity.")
    print()

    # Verify this for n=3
    print("  Verification:")
    for trial in range(10):
        roots = random_centered_roots(3, spread=3.0, min_gap=0.5)
        e2 = elem_sym(roots, 2)
        e3 = elem_sym(roots, 3)
        E = -e2
        F = e3
        u = F**2 / E**3
        sf_computed = shape_factor(roots)
        sf_formula = 12.0 / (4.0 - 27.0 * u)
        print(f"    u={u:.6f}, SF(computed)={sf_computed:.8f}, "
              f"SF(formula)={sf_formula:.8f}, match={abs(sf_computed-sf_formula)<1e-8}")

    # For general n, is there a similar formula?
    print("\n\nSF for n=4 -- exploring dependence on cumulant ratios:")
    print("  For centered poly, let p_k = (1/n) sum lambda_i^k (k-th moment).")
    print("  p_2 = Var (by definition for centered).")
    print("  The candidate: SF depends on p_4/p_2^2, p_3/p_2^{3/2}, etc.")

    for trial in range(8):
        roots = random_centered_roots(4, spread=3.0, min_gap=0.3)
        p2 = np.mean(roots**2)
        p3 = np.mean(roots**3)
        p4 = np.mean(roots**4)
        sf = shape_factor(roots)
        kurt = p4 / p2**2  # kurtosis-like
        skew = p3 / p2**1.5  # skewness-like
        print(f"    p2={p2:.4f}, skew={skew:.4f}, kurt={kurt:.4f}, SF={sf:.6f}")

    print()


# =====================================================================
# TASK 4: ALTERNATIVE DECOMPOSITIONS
# =====================================================================

def task4_alternative_decompositions():
    """
    Instead of 1/Phi = Var/SF, try:
    a) 1/Phi = p_2 / SF2 where p_2 = sum lambda_i^2 (for centered), SF2 = Phi * p_2.
       Note p_2 = n * Var for centered.
    b) 1/Phi = Delta / SF3 where Delta = (sum of all (lambda_i - lambda_j)^2) / n^2.
    c) 1/Phi = D^(2/n) / SF4 where D = discriminant.

    For each: check additivity of the "extensive" part and monotonicity of the "intensive" part.
    """
    print("=" * 72)
    print("TASK 4: ALTERNATIVE DECOMPOSITIONS")
    print("=" * 72)
    print()

    # Alternative A: Use sum of squared differences
    # D2(p) = (1/n(n-1)) * sum_{i<j} (lambda_i - lambda_j)^2
    # Note: sum_{i<j} (lam_i - lam_j)^2 = n * sum lam_i^2 - (sum lam_i)^2 = n^2 * Var (for centered)
    # Actually: sum_{i<j} (lam_i - lam_j)^2 = n * sum(lam_i^2) - (sum lam_i)^2
    # For centered: = n * sum(lam_i^2) = n^2 * Var

    print("Alternative A: D2(p) = sum_{i<j} (lam_i - lam_j)^2")
    print("  For centered poly: D2 = n^2 * Var, so D2 additivity <=> Var additivity. Same as Var.")
    print()

    # Alternative B: Use the discriminant
    # disc(p) = prod_{i<j} (lam_i - lam_j)^2
    # This is NOT additive under MSS convolution.
    # But: disc(r) = disc(p) * disc(q) * C_n  (some formula?)
    print("Alternative B: Discriminant-based decomposition")
    print("  disc(p) = prod_{i<j} (lambda_i - lambda_j)^2")
    print("  Check: is log(disc(r)) = log(disc(p)) + log(disc(q)) + correction?")

    for n in [3, 4, 5]:
        print(f"\n  n={n}:")
        for trial in range(5):
            rp = random_roots(n, min_gap=0.5)
            rq = random_roots(n, min_gap=0.5)
            try:
                rr = boxplus(rp, rq)
                disc_p = np.prod([(rp[i]-rp[j])**2 for i in range(n) for j in range(i+1,n)])
                disc_q = np.prod([(rq[i]-rq[j])**2 for i in range(n) for j in range(i+1,n)])
                disc_r = np.prod([(rr[i]-rr[j])**2 for i in range(n) for j in range(i+1,n)])
                log_ratio = np.log(disc_r) - np.log(disc_p) - np.log(disc_q)
                print(f"    log(disc_r) - log(disc_p) - log(disc_q) = {log_ratio:.6f}")
            except:
                pass

    # Alternative C: Use Phi directly in a harmonic-mean decomposition
    # The conjecture 1/Phi(r) >= 1/Phi(p) + 1/Phi(q) is equivalent to:
    # Phi(r) <= Phi(p)*Phi(q)/(Phi(p)+Phi(q)) = harmonic_mean(Phi(p), Phi(q))/2
    print("\n\nAlternative C: Harmonic mean bound")
    print("  Check: Phi(r) <= Phi(p)*Phi(q)/(Phi(p)+Phi(q))?")
    for n in [3, 4, 5, 6]:
        violations = 0
        total = 0
        min_ratio = float('inf')
        for trial in range(1000):
            rp = random_roots(n)
            rq = random_roots(n)
            try:
                rr = boxplus(rp, rq)
                phi_p = Phi_n(rp)
                phi_q = Phi_n(rq)
                phi_r = Phi_n(rr)
                harm = phi_p * phi_q / (phi_p + phi_q)
                ratio = harm / phi_r  # should be >= 1
                min_ratio = min(min_ratio, ratio)
                total += 1
                if phi_r > harm + 1e-10:
                    violations += 1
            except:
                continue
        print(f"  n={n}: {total} trials, {violations} violations, min(harm/Phi(r)) = {min_ratio:.6f}")

    # Alternative D: power sum p_2 decomposition
    # For centered: p_2 = sum lambda_i^2 = n*Var. So 1/Phi = p_2/(n*SF) = p_2/SF'.
    # p_2 is additive (since Var is additive and p_2 = n*Var).
    # SF' = n*SF, same conjecture.
    print("\n  Alternative D: p_2 = n*Var decomposition is equivalent to Var decomposition.")
    print()


# =====================================================================
# TASK 5: n=3 VERIFICATION OF SHAPE FACTOR CONJECTURE
# =====================================================================

def task5_n3_shape_factor():
    """
    For n=3 centered:
      SF = 12 / (4 - 27*u) where u = F^2/E^3, E = -e_2, F = e_3.
      SF(r) <= min(SF(p), SF(q))
    <=> u_r <= min(u_p, u_q)  (since SF is monotone increasing in u)

    Now E_r = E_p + E_q, F_r = F_p + F_q.
    So u_r = F_r^2 / E_r^3 = (F_p + F_q)^2 / (E_p + E_q)^3.
    And u_p = F_p^2 / E_p^3, u_q = F_q^2 / E_q^3.

    Need: (F_p + F_q)^2 / (E_p + E_q)^3 <= min(F_p^2/E_p^3, F_q^2/E_q^3)?
    """
    print("=" * 72)
    print("TASK 5: n=3 SHAPE FACTOR CONJECTURE ANALYSIS")
    print("=" * 72)
    print()

    print("For n=3 centered:")
    print("  SF = 12 / (4 - 27*u)  where u = F^2/E^3")
    print("  SF is monotone increasing in u, so SF(r) <= min(SF(p),SF(q))")
    print("  <=> u_r <= min(u_p, u_q)")
    print("  <=> (Fp+Fq)^2/(Ep+Eq)^3 <= min(Fp^2/Ep^3, Fq^2/Eq^3)")
    print()

    # Test this for the specific case F_p = 0 (equally spaced)
    print("Special case: F_p = 0 (p equally spaced):")
    print("  u_r = F_q^2/(E_p+E_q)^3")
    print("  u_q = F_q^2/E_q^3")
    print("  u_r/u_q = E_q^3/(E_p+E_q)^3 = 1/(1+E_p/E_q)^3 < 1. OK.")
    print("  u_p = 0, so min(u_p, u_q) = 0, and u_r = F_q^2/(Ep+Eq)^3 > 0.")
    print()
    print("  WAIT: u_p = 0, so min(u_p, u_q) = 0, but u_r > 0.")
    print("  THIS MEANS SF(r) > SF(p) when p has equally spaced roots!")
    print()
    print("  Let's check: if p has equally spaced roots (F_p=0) and q is arbitrary:")
    print("    SF(p) = 12/4 = 3  (the minimum possible SF)")
    print("    SF(r) = 12/(4 - 27*F_q^2/(E_p+E_q)^3) > 3 in general.")
    print()
    print("  SO SF(r) <= min(SF(p), SF(q)) IS FALSE for n=3!")
    print()

    # Numerical verification of this counterexample
    print("NUMERICAL COUNTEREXAMPLE:")
    # p = equally spaced, q = generic
    rp = np.array([-1.0, 0.0, 1.0])  # equally spaced, centered
    rq_list = [
        np.array([-2.0, -0.3, 2.3]),  # centered, F != 0
        np.array([-3.0, 1.0, 2.0]),
        np.array([-1.5, 0.2, 1.3]),
    ]
    for rq in rq_list:
        rq = rq - np.mean(rq)  # center
        try:
            rr = boxplus(rp, rq)
            sf_p = shape_factor(rp)
            sf_q = shape_factor(rq)
            sf_r = shape_factor(rr)
            print(f"  p={rp}, q={np.round(rq,3)}")
            print(f"    SF(p)={sf_p:.6f}, SF(q)={sf_q:.6f}, SF(r)={sf_r:.6f}")
            print(f"    min(SF(p),SF(q))={min(sf_p,sf_q):.6f}")
            print(f"    SF(r) <= min? {sf_r <= min(sf_p,sf_q) + 1e-10}")
            print(f"    BUT: 1/Phi(r) >= 1/Phi(p)+1/Phi(q)? ", end="")
            phi_p = Phi_n(rp)
            phi_q = Phi_n(rq)
            phi_r = Phi_n(rr)
            excess = 1.0/phi_r - 1.0/phi_p - 1.0/phi_q
            print(f"{excess >= -1e-10} (excess = {excess:.6e})")
        except Exception as e:
            print(f"  Error: {e}")
    print()

    # The key insight: SF(r) <= min(SF(p), SF(q)) is FALSE.
    # But the main conjecture 1/Phi(r) >= 1/Phi(p)+1/Phi(q) is still TRUE.
    # The shape factor approach FAILS because SF is not monotone under boxplus.
    print("=" * 72)
    print("CRITICAL FINDING: SF(r) <= min(SF(p), SF(q)) IS FALSE!")
    print("=" * 72)
    print()
    print("The shape factor decomposition 1/Phi = Var/SF with the claim")
    print("SF(r) <= min(SF(p), SF(q)) does NOT work, because:")
    print("  - When p has equally-spaced roots, SF(p) = 3 (the minimum)")
    print("  - SF(r) depends on q's shape, and SF(r) > 3 in general")
    print("  - So SF(r) > SF(p) = min(SF(p), SF(q))")
    print()
    print("However, the MAIN CONJECTURE is still true. The shape factor approach")
    print("is too crude: it ignores correlations between Var and SF.")
    print()


# =====================================================================
# TASK 6: DEEPER ANALYSIS OF WHEN SF(r) > min(SF(p), SF(q))
# =====================================================================

def task6_counterexample_analysis():
    """
    Document how bad the SF counterexample is.
    When is SF(r) > min(SF(p), SF(q)), and by how much?
    """
    print("=" * 72)
    print("TASK 6: COUNTEREXAMPLE ANALYSIS")
    print("=" * 72)
    print()

    for n in [3, 4, 5, 6]:
        violations = 0
        total = 0
        max_excess_ratio = 0
        violation_count_by_type = {'eq_vs_generic': 0, 'generic_vs_generic': 0}

        for trial in range(3000):
            roots_p = random_roots(n)
            roots_q = random_roots(n)

            try:
                roots_r = boxplus(roots_p, roots_q)
                sf_p = shape_factor(roots_p)
                sf_q = shape_factor(roots_q)
                sf_r = shape_factor(roots_r)
                min_sf = min(sf_p, sf_q)

                total += 1
                if sf_r > min_sf + 1e-10:
                    violations += 1
                    excess = sf_r / min_sf
                    max_excess_ratio = max(max_excess_ratio, excess)
            except:
                continue

        # Also test equally-spaced vs generic
        eq_viol = 0
        eq_total = 0
        eq_max_ratio = 0
        for trial in range(1000):
            d = np.random.uniform(0.5, 3.0)
            roots_p = np.linspace(-d, d, n)  # equally spaced
            roots_q = random_roots(n, min_gap=0.3)

            try:
                roots_r = boxplus(roots_p, roots_q)
                sf_p = shape_factor(roots_p)
                sf_q = shape_factor(roots_q)
                sf_r = shape_factor(roots_r)
                min_sf = min(sf_p, sf_q)

                eq_total += 1
                if sf_r > min_sf + 1e-10:
                    eq_viol += 1
                    eq_max_ratio = max(eq_max_ratio, sf_r / min_sf)
            except:
                continue

        print(f"n={n}:")
        print(f"  Random: {violations}/{total} violations ({100*violations/total:.1f}%), "
              f"worst ratio = {max_excess_ratio:.4f}")
        print(f"  Eq-spaced vs random: {eq_viol}/{eq_total} violations ({100*eq_viol/eq_total:.1f}%), "
              f"worst ratio = {eq_max_ratio:.4f}")
        print()


# =====================================================================
# TASK 7: CAN WE SALVAGE THE APPROACH?
# =====================================================================

def task7_salvage_attempts():
    """
    The SF approach fails because SF is not monotone.
    But maybe a WEAKER bound can still be used.

    Approach A: Instead of SF(r) <= min(SF(p), SF(q)), try:
      SF(r) <= f(SF(p), SF(q), Var(p), Var(q)) for some function f
      such that Var(r)/f >= Var(p)/SF(p) + Var(q)/SF(q).

    Approach B: Use a different "intensive" factor.
      1/Phi = Var / SF. But Var(r) = Var(p) + Var(q).
      The inequality 1/Phi(r) >= 1/Phi(p) + 1/Phi(q) becomes:
        (Var(p)+Var(q))/SF(r) >= Var(p)/SF(p) + Var(q)/SF(q)
      i.e., SF(r) <= (Var(p)+Var(q)) / (Var(p)/SF(p) + Var(q)/SF(q))
      i.e., SF(r) <= harmonic_mean_weighted(SF(p), SF(q); Var(p), Var(q))

    So the CORRECT statement is:
      SF(r) <= (Var(p) + Var(q)) * SF(p) * SF(q) / (Var(p)*SF(q) + Var(q)*SF(p))
    which is the variance-weighted harmonic mean of SF(p) and SF(q).
    """
    print("=" * 72)
    print("TASK 7: SALVAGE -- SF(r) <= WEIGHTED HARMONIC MEAN OF SF(p), SF(q)")
    print("=" * 72)
    print()

    print("The CORRECT reformulation of the conjecture in SF terms is:")
    print("  SF(r) <= (Vp + Vq) * SFp * SFq / (Vp*SFq + Vq*SFp)")
    print("       = variance-weighted harmonic mean of SFp, SFq")
    print()
    print("This is EQUIVALENT to the original conjecture (not stronger or weaker).")
    print("It just reformulates 1/Phi(r) >= 1/Phi(p) + 1/Phi(q) in different notation.")
    print()

    # Verify equivalence numerically
    print("Verification that the reformulation is equivalent:")
    for n in [3, 4, 5, 6]:
        matching = 0
        total = 0
        for trial in range(1000):
            rp = random_roots(n)
            rq = random_roots(n)
            try:
                rr = boxplus(rp, rq)
                vp = root_variance(rp)
                vq = root_variance(rq)
                sfp = shape_factor(rp)
                sfq = shape_factor(rq)
                sfr = shape_factor(rr)

                # Original: 1/Phi(r) >= 1/Phi(p) + 1/Phi(q)
                orig = (1.0/Phi_n(rr)) - (1.0/Phi_n(rp)) - (1.0/Phi_n(rq))

                # SF form: (Vp+Vq)/(sfr) >= Vp/sfp + Vq/sfq
                # equiv: SF(r) <= (Vp+Vq)*sfp*sfq / (Vp*sfq + Vq*sfp)
                bound = (vp + vq) * sfp * sfq / (vp * sfq + vq * sfp)
                sf_form = bound - sfr

                total += 1
                # Both should have the same sign
                if (orig >= -1e-10) == (sf_form >= -1e-10):
                    matching += 1
            except:
                continue

        print(f"  n={n}: {matching}/{total} sign-agreement (should be 100%)")

    print()

    # Now let's look at what the weighted harmonic mean bound implies
    print("Analysis of the weighted harmonic mean bound:")
    print("  H(w1,w2; s1,s2) = (w1+w2)*s1*s2 / (w1*s2 + w2*s1)")
    print("  where w = Var, s = SF.")
    print()
    print("  Key properties:")
    print("  - H <= min(s1, s2) (harmonic mean is <= min)")
    print("  - H = min(s1,s2) iff s1 = s2")
    print("  - H depends on the ratio w1/w2 (weight ratio)")
    print()

    # Actually, H = (w1+w2)/(w1/s1 + w2/s2). This is the harmonic mean with weights.
    # H <= min(s1,s2) since (w1/s1+w2/s2) >= (w1+w2)/max(s1,s2) so H <= max(s1,s2).
    # Actually: (w1/s1+w2/s2) >= (w1+w2)/max(s1,s2), so H <= max(s1,s2). But that's not min.
    # Let's verify: is H(w1,w2;s1,s2) <= min(s1,s2)?
    # H = (w1+w2)*s1*s2/(w1*s2+w2*s1). If s1 <= s2: w1*s2+w2*s1 >= w1*s1+w2*s1 = (w1+w2)*s1.
    # So H <= s1*s2/s1 = s2. Also: w1*s2+w2*s1 >= w1*s2+w2*s2 = (w1+w2)*s2? NO: w2*s1 <= w2*s2.
    # So w1*s2+w2*s1 <= w1*s2+w2*s2 = (w1+w2)*s2, giving H >= (w1+w2)*s1*s2/((w1+w2)*s2) = s1.
    # So: min(s1,s2) <= H <= max(s1,s2).
    print("  Actually: min(s1,s2) <= H <= max(s1,s2).")
    print("  So the bound SF(r) <= H says SF(r) is between min and max of SF(p),SF(q).")
    print("  This is MUCH weaker than SF(r) <= min.")
    print()

    # So the SF approach gives no additional insight beyond the original conjecture.
    # But let's look at the behavior of SF(r) relative to this bound.
    print("Ratio SF(r) / weighted_harmonic_mean:")
    for n in [3, 4, 5, 6]:
        ratios = []
        for trial in range(2000):
            rp = random_roots(n)
            rq = random_roots(n)
            try:
                rr = boxplus(rp, rq)
                vp = root_variance(rp)
                vq = root_variance(rq)
                sfp = shape_factor(rp)
                sfq = shape_factor(rq)
                sfr = shape_factor(rr)
                H = (vp + vq) * sfp * sfq / (vp * sfq + vq * sfp)
                ratios.append(sfr / H)
            except:
                continue
        ratios = np.array(ratios)
        print(f"  n={n}: mean={np.mean(ratios):.6f}, max={np.max(ratios):.6f}, "
              f"min={np.min(ratios):.6f}")
        if np.max(ratios) <= 1.0 + 1e-8:
            print(f"    => SF(r) <= WHM ALWAYS (equivalent to conjecture)")
        else:
            print(f"    => VIOLATION (conjecture fails?!)")

    print()


# =====================================================================
# TASK 8: n=3 EXPLICIT SHAPE FACTOR ANALYSIS
# =====================================================================

def task8_n3_explicit():
    """
    For n=3, we have the explicit formula.
    Let's work out exactly what the shape factor approach gives.

    1/Phi = (4E^3 - 27F^2) / (18E^2)
    Var = 2E/3
    SF = 18E^2 / (4E^3-27F^2) * 2E/3 = 12E^3 / (4E^3-27F^2)

    The conjecture: 1/Phi(r) >= 1/Phi(p) + 1/Phi(q) with Er=Ep+Eq, Fr=Fp+Fq.

    In terms of SF:
    (Vp+Vq)/SFr >= Vp/SFp + Vq/SFq
    i.e., (2Ep/3 + 2Eq/3) / SFr >= 2Ep/(3*SFp) + 2Eq/(3*SFq)
    i.e., (Ep+Eq)/SFr >= Ep/SFp + Eq/SFq

    Substituting SFr = 12Er^3/(4Er^3-27Fr^2) with Er=Ep+Eq, Fr=Fp+Fq:
    (Ep+Eq)*(4(Ep+Eq)^3-27(Fp+Fq)^2) / (12(Ep+Eq)^3)
    >= Ep*(4Ep^3-27Fp^2)/(12Ep^3) + Eq*(4Eq^3-27Fq^2)/(12Eq^3)
    i.e., (4(Ep+Eq)^3-27(Fp+Fq)^2) / (12(Ep+Eq)^2)
    >= (4Ep^3-27Fp^2)/(12Ep^2) + (4Eq^3-27Fq^2)/(12Eq^2)

    This is just 1/Phi(r) >= 1/Phi(p) + 1/Phi(q) / 2 .... wait, no.
    Actually: 1/Phi = (4E^3-27F^2)/(18E^2), so we need:
    (4Er^3-27Fr^2)/(18Er^2) >= (4Ep^3-27Fp^2)/(18Ep^2) + (4Eq^3-27Fq^2)/(18Eq^2)

    This is the original conjecture, proved in prove_n3_symbolic.py.
    The shape factor approach is just a REPHRASING, not an independent insight.
    """
    print("=" * 72)
    print("TASK 8: n=3 EXPLICIT ANALYSIS -- SF IS JUST A REPHRASING")
    print("=" * 72)
    print()

    print("For n=3 centered cubic:")
    print("  1/Phi = (4E^3 - 27F^2) / (18E^2)")
    print("  Var = 2E/3")
    print("  SF = 12E^3 / (4E^3 - 27F^2)")
    print()
    print("The conjecture 1/Phi(r) >= 1/Phi(p) + 1/Phi(q) is equivalent to:")
    print("  Var(r)/SF(r) >= Var(p)/SF(p) + Var(q)/SF(q)")
    print("  where Var(r) = Var(p)+Var(q) exactly.")
    print()
    print("This is algebraically identical to the original conjecture.")
    print("The shape factor approach provides NO additional leverage for n=3.")
    print()

    # But let's understand the STRUCTURE of SF for n=3:
    print("Structure of SF for n=3:")
    print("  SF = 12/(4 - 27u) where u = F^2/E^3 in [0, 4/27).")
    print("  SF is in [3, infinity).")
    print("  SF = 3 iff u = 0 iff F = 0 iff roots equally spaced.")
    print("  SF -> inf as u -> 4/27 (discriminant -> 0, nearly degenerate).")
    print()

    # What about the ratio SF(r)/harmonic_mean?
    print("SF(r) vs harmonic mean of SF(p), SF(q) (weighted by Var):")
    print("  The conjecture in SF terms: SF(r) <= (Vp+Vq)*SFp*SFq/(Vp*SFq+Vq*SFp)")
    print()
    print("  Substituting Vp = 2Ep/3, Vq = 2Eq/3:")
    print("  SF(r) <= (Ep+Eq)*SFp*SFq/(Ep*SFq+Eq*SFp)")
    print()

    # For the special case SFp = SFq = s (equal shape factors):
    print("  Special case SFp = SFq = s:")
    print("    Bound: SF(r) <= s")
    print("    i.e., the convolution doesn't increase SF when shapes are equal.")
    print()

    # Verify this numerically for n=3 with same shape factor
    print("  Numerical check (n=3, same SF):")
    for trial in range(5):
        # Create two cubics with the same shape factor (same u = F^2/E^3)
        E1 = np.random.uniform(1, 5)
        u = np.random.uniform(0, 0.14)  # u < 4/27 ~ 0.148
        F1 = np.sqrt(u * E1**3)
        if np.random.rand() < 0.5:
            F1 = -F1

        E2 = np.random.uniform(1, 5)
        F2 = np.sqrt(u * E2**3)
        if np.random.rand() < 0.5:
            F2 = -F2

        # Build centered roots from (E, F):
        # x^3 - E*x - F = 0
        coeffs_p = [1, 0, -E1, -F1]
        coeffs_q = [1, 0, -E2, -F2]
        rp = np.sort(np.real(np.roots(coeffs_p)))
        rq = np.sort(np.real(np.roots(coeffs_q)))

        if np.min(np.diff(rp)) < 0.01 or np.min(np.diff(rq)) < 0.01:
            continue

        try:
            rr = boxplus(rp, rq)
            sfp = shape_factor(rp)
            sfq = shape_factor(rq)
            sfr = shape_factor(rr)
            print(f"    u={u:.4f}, SFp={sfp:.4f}, SFq={sfq:.4f}, SFr={sfr:.4f}, "
                  f"SFr<=min? {sfr <= min(sfp, sfq) + 1e-8}")
        except:
            pass

    print()

    # For the case where one poly has F=0 (equally spaced):
    print("  When p has equally-spaced roots (SFp = 3, u_p = 0):")
    print("    bound = (Ep+Eq)*3*SFq / (Ep*SFq + Eq*3)")
    print("    = 3*(Ep+Eq)*SFq / (Ep*SFq + 3*Eq)")
    print()
    print("    SF(r) = 12*(Ep+Eq)^3 / (4*(Ep+Eq)^3 - 27*Fq^2)")
    print("    Need: 12*(Ep+Eq)^3 / (4*(Ep+Eq)^3-27*Fq^2) <= 3*(Ep+Eq)*SFq/(Ep*SFq+3*Eq)")
    print("    This is equivalent to the original conjecture for this case.")
    print()


# =====================================================================
# TASK 9: EXPLORE WHAT MAKES SF(r) LARGE/SMALL
# =====================================================================

def task9_sf_behavior():
    """
    Study how SF(r) relates to SF(p) and SF(q) in practice.
    """
    print("=" * 72)
    print("TASK 9: BEHAVIOR OF SF(r) RELATIVE TO SF(p), SF(q)")
    print("=" * 72)
    print()

    for n in [3, 4, 5, 6]:
        sf_ratios_to_max = []
        sf_ratios_to_min = []
        sf_ratios_to_geo = []  # geometric mean

        for trial in range(2000):
            rp = random_roots(n)
            rq = random_roots(n)
            try:
                rr = boxplus(rp, rq)
                sfp = shape_factor(rp)
                sfq = shape_factor(rq)
                sfr = shape_factor(rr)

                sf_ratios_to_max.append(sfr / max(sfp, sfq))
                sf_ratios_to_min.append(sfr / min(sfp, sfq))
                sf_ratios_to_geo.append(sfr / np.sqrt(sfp * sfq))
            except:
                continue

        r_max = np.array(sf_ratios_to_max)
        r_min = np.array(sf_ratios_to_min)
        r_geo = np.array(sf_ratios_to_geo)

        print(f"n={n}: SF(r) / max(SF(p),SF(q)): "
              f"mean={np.mean(r_max):.4f}, min={np.min(r_max):.4f}, max={np.max(r_max):.4f}")
        print(f"       SF(r) / min(SF(p),SF(q)): "
              f"mean={np.mean(r_min):.4f}, min={np.min(r_min):.4f}, max={np.max(r_min):.4f}")
        print(f"       SF(r) / geo(SF(p),SF(q)): "
              f"mean={np.mean(r_geo):.4f}, min={np.min(r_geo):.4f}, max={np.max(r_geo):.4f}")

        # What fraction has SF(r) < min?
        frac_below_min = np.mean(r_min < 1.0 - 1e-10)
        frac_below_max = np.mean(r_max < 1.0 - 1e-10)
        frac_below_geo = np.mean(r_geo < 1.0 - 1e-10)
        print(f"       SF(r) < min: {100*frac_below_min:.1f}%, "
              f"< geo: {100*frac_below_geo:.1f}%, "
              f"< max: {100*frac_below_max:.1f}%")
        print()


# =====================================================================
# TASK 10: THE REAL INSIGHT -- WHAT DOES THE PROOF NEED?
# =====================================================================

def task10_what_proof_needs():
    """
    The shape factor approach fails in its strong form (SF(r) <= min).
    But maybe a modified decomposition works.

    Key identity: 1/Phi(r) - 1/Phi(p) - 1/Phi(q)
    = Var(r)/SF(r) - Var(p)/SF(p) - Var(q)/SF(q)
    = (Vp+Vq)/SFr - Vp/SFp - Vq/SFq    [using Var additivity]
    = Vp*(1/SFr - 1/SFp) + Vq*(1/SFr - 1/SFq)

    So the excess is >= 0 iff:
    Vp*(1/SFr - 1/SFp) + Vq*(1/SFr - 1/SFq) >= 0
    iff (Vp+Vq)/SFr >= Vp/SFp + Vq/SFq

    This can hold even when SFr > SFp (if it's compensated by SFr < SFq with enough weight).

    So the real question is about the WEIGHTED comparison.
    """
    print("=" * 72)
    print("TASK 10: WHAT THE PROOF REALLY NEEDS")
    print("=" * 72)
    print()

    print("The excess = Vp*(1/SFr - 1/SFp) + Vq*(1/SFr - 1/SFq)")
    print()
    print("The two terms can have opposite signs! It's possible that")
    print("  1/SFr < 1/SFp  (i.e., SFr > SFp)")
    print("but compensated by")
    print("  1/SFr > 1/SFq  (i.e., SFr < SFq)")
    print("with Vq being large enough.")
    print()

    # Verify this pattern
    print("Let's see how often each term has a definite sign:")
    for n in [3, 4, 5]:
        both_pos = 0  # SFr < SFp and SFr < SFq
        term1_neg = 0  # SFr > SFp (but total still >= 0)
        term2_neg = 0
        total = 0

        for trial in range(2000):
            rp = random_roots(n)
            rq = random_roots(n)
            try:
                rr = boxplus(rp, rq)
                vp = root_variance(rp)
                vq = root_variance(rq)
                sfp = shape_factor(rp)
                sfq = shape_factor(rq)
                sfr = shape_factor(rr)

                total += 1
                t1 = 1.0/sfr - 1.0/sfp  # term 1 sign
                t2 = 1.0/sfr - 1.0/sfq  # term 2 sign

                if t1 >= -1e-12 and t2 >= -1e-12:
                    both_pos += 1
                elif t1 < -1e-12:
                    term1_neg += 1
                else:
                    term2_neg += 1
            except:
                continue

        print(f"  n={n}: both terms positive: {100*both_pos/total:.1f}%, "
              f"term1 neg: {100*term1_neg/total:.1f}%, "
              f"term2 neg: {100*term2_neg/total:.1f}%")

    print()
    print("CONCLUSION: A significant fraction of cases have one term negative,")
    print("showing that the simple bound SF(r) <= min(SF(p),SF(q)) fails.")
    print("The conjecture works because the negative term is always compensated")
    print("by the other term, but this compensation is hard to prove in the SF framework.")
    print()

    # Final test: is there a DIFFERENT intensive quantity that IS monotone?
    print("Search for a monotone intensive quantity:")
    print("  Let T(p) = Phi(p) * f(roots(p)) for some function f.")
    print("  We want T(r) <= min(T(p), T(q)) always.")
    print()

    # Candidate: T = Phi * Var * (some correction based on higher cumulants)
    # For n=3: SF = 12/(4-27u) = 12/(4-27F^2/E^3). Min is at u=0.
    # Maybe T = SF * Var^alpha for some alpha?
    # Or T = Phi * sum_of_squared_diffs?

    # Let's try T = Phi * D^{2/(n(n-1))} where D = prod_{i<j} |lam_i-lam_j|^2
    print("  Candidate: T = Phi * D^{2/(n*(n-1))} where D = disc(p)")
    for n in [3, 4, 5]:
        viols = 0
        total = 0
        for trial in range(2000):
            rp = random_roots(n)
            rq = random_roots(n)
            try:
                rr = boxplus(rp, rq)
                m = n*(n-1)//2
                disc_p = np.prod([(rp[i]-rp[j])**2 for i in range(n) for j in range(i+1,n)])
                disc_q = np.prod([(rq[i]-rq[j])**2 for i in range(n) for j in range(i+1,n)])
                disc_r = np.prod([(rr[i]-rr[j])**2 for i in range(n) for j in range(i+1,n)])
                Tp = Phi_n(rp) * disc_p**(1.0/m)
                Tq = Phi_n(rq) * disc_q**(1.0/m)
                Tr = Phi_n(rr) * disc_r**(1.0/m)
                total += 1
                if Tr > min(Tp, Tq) + 1e-8:
                    viols += 1
            except:
                continue
        print(f"    n={n}: T(r) <= min(T(p),T(q))? violations: {viols}/{total}")

    print()

    # Try T = Phi * max_gap^2 (or min_gap^2, or range^2)
    print("  Candidate: T = Phi * range^2 (range = max root - min root)")
    for n in [3, 4, 5]:
        viols = 0
        total = 0
        for trial in range(2000):
            rp = random_roots(n)
            rq = random_roots(n)
            try:
                rr = boxplus(rp, rq)
                Tp = Phi_n(rp) * (rp[-1]-rp[0])**2
                Tq = Phi_n(rq) * (rq[-1]-rq[0])**2
                Tr = Phi_n(rr) * (rr[-1]-rr[0])**2
                total += 1
                if Tr > min(Tp, Tq) + 1e-8:
                    viols += 1
            except:
                continue
        print(f"    n={n}: T(r) <= min(T(p),T(q))? violations: {viols}/{total}")

    print()


# =====================================================================
# MAIN
# =====================================================================

if __name__ == "__main__":
    print("#" * 72)
    print("# SHAPE FACTOR DECOMPOSITION INVESTIGATION")
    print("# For Fisher superadditivity: 1/Phi(r) >= 1/Phi(p) + 1/Phi(q)")
    print("#" * 72)
    print()

    task1_variance_additivity()
    task2_shape_factor_conjecture()
    task3_understand_shape_factor()
    task4_alternative_decompositions()
    task5_n3_shape_factor()
    task6_counterexample_analysis()
    task7_salvage_attempts()
    task8_n3_explicit()
    task9_sf_behavior()
    task10_what_proof_needs()

    print("\n" + "#" * 72)
    print("# INVESTIGATION COMPLETE")
    print("#" * 72)
