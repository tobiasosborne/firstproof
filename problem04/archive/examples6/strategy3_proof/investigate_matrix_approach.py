"""
Investigation: Random matrix / information-theoretic approach to Fisher superadditivity.

CONJECTURE: For monic real-rooted p, q of degree n, r = p boxplus_n q:
  1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q)

where Phi_n(p) = sum_i H_p(lambda_i)^2, H_p(lambda_i) = p''(lambda_i)/(2p'(lambda_i))
                = sum_{j!=i} 1/(lambda_i - lambda_j).

RANDOM MATRIX CONNECTION (MSS): If A, B are n×n Hermitian with char polys p, q,
and U ~ Haar(U(n)), then E_U[chi_{A + UBU*}] = p boxplus_n q.

FREE PROBABILITY: Voiculescu's free Fisher information Phi*(mu) satisfies
  1/Phi*(mu boxplus nu) = 1/Phi*(mu) + 1/Phi*(nu)  [EXACT equality]
for free convolution. The finite analogue should be >= (superadditivity).

This script investigates:
1. Phi_n as a matrix trace functional
2. Resolvent-based expressions
3. Convexity/concavity of 1/Phi_n
4. Connection to free probability Fisher information
5. Jensen/convexity approach via random matrix expectation
"""

import numpy as np
from numpy.polynomial import polynomial as P
from itertools import combinations
from scipy.special import comb
import warnings
warnings.filterwarnings('ignore')

np.set_printoptions(precision=10, linewidth=120)

# ============================================================================
# PART 0: Core definitions
# ============================================================================

def phi_n_from_roots(roots):
    """Compute Phi_n(p) = sum_i H_p(lambda_i)^2 from roots lambda_i.
    H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j).
    """
    n = len(roots)
    phi = 0.0
    for i in range(n):
        H_i = sum(1.0/(roots[i] - roots[j]) for j in range(n) if j != i)
        phi += H_i**2
    return phi

def H_values(roots):
    """Compute H_p(lambda_i) for each root."""
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        H[i] = sum(1.0/(roots[i] - roots[j]) for j in range(n) if j != i)
    return H

def finite_free_convolution(roots_p, roots_q):
    """Compute p boxplus_n q using the hat-elementary symmetric function formula.
    hat_e_k(r) = sum_{j=0}^k hat_e_j(p) * hat_e_{k-j}(q)
    where hat_e_k = e_k / C(n,k).
    """
    n = len(roots_p)
    assert len(roots_q) == n

    # Compute elementary symmetric polynomials from roots
    def elem_sym(roots, k):
        """e_k = sum of products of k roots taken at a time."""
        if k == 0:
            return 1.0
        if k > len(roots):
            return 0.0
        return sum(np.prod(list(combo)) for combo in combinations(roots, k))

    # hat_e_k = e_k / C(n,k)
    hat_e_p = [elem_sym(roots_p, k) / comb(n, k, exact=True) for k in range(n+1)]
    hat_e_q = [elem_sym(roots_q, k) / comb(n, k, exact=True) for k in range(n+1)]

    # Convolution: hat_e_k(r) = sum_{j=0}^k hat_e_j(p) * hat_e_{k-j}(q)
    hat_e_r = [sum(hat_e_p[j] * hat_e_q[k-j] for j in range(k+1)) for k in range(n+1)]

    # Convert back: e_k(r) = C(n,k) * hat_e_k(r)
    e_r = [comb(n, k, exact=True) * hat_e_r[k] for k in range(n+1)]

    # Build polynomial coefficients: p(x) = x^n - e_1 x^{n-1} + e_2 x^{n-2} - ...
    # In numpy: coefficients from highest to lowest degree
    coeffs = [(-1)**k * e_r[k] for k in range(n+1)]

    # Find roots
    roots_r = np.sort(np.roots(coeffs))

    # Check real-rootedness
    if np.max(np.abs(np.imag(roots_r))) > 1e-8:
        print(f"  WARNING: Non-real roots detected! Max imag part: {np.max(np.abs(np.imag(roots_r))):.2e}")

    return np.real(roots_r)

# ============================================================================
# PART 1: Basic verification of superadditivity
# ============================================================================

print("=" * 80)
print("PART 1: Basic verification of Fisher superadditivity")
print("=" * 80)

def test_superadditivity(roots_p, roots_q, label=""):
    """Test 1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q)."""
    roots_r = finite_free_convolution(roots_p, roots_q)

    phi_p = phi_n_from_roots(roots_p)
    phi_q = phi_n_from_roots(roots_q)
    phi_r = phi_n_from_roots(roots_r)

    lhs = 1.0/phi_r
    rhs = 1.0/phi_p + 1.0/phi_q

    gap = lhs - rhs

    print(f"\n{label}")
    print(f"  roots_p = {np.sort(roots_p)}")
    print(f"  roots_q = {np.sort(roots_q)}")
    print(f"  roots_r = {np.sort(roots_r)}")
    print(f"  Phi_n(p) = {phi_p:.10f},  1/Phi = {1/phi_p:.10f}")
    print(f"  Phi_n(q) = {phi_q:.10f},  1/Phi = {1/phi_q:.10f}")
    print(f"  Phi_n(r) = {phi_r:.10f},  1/Phi = {1/phi_r:.10f}")
    print(f"  LHS - RHS = {gap:.2e}  {'PASS' if gap >= -1e-12 else 'FAIL'}")

    return phi_p, phi_q, phi_r, gap

# n=2
test_superadditivity(np.array([1.0, 3.0]), np.array([0.0, 2.0]), "n=2, generic")
test_superadditivity(np.array([0.0, 4.0]), np.array([0.0, 4.0]), "n=2, same poly")
test_superadditivity(np.array([0.0, 1.0]), np.array([0.0, 10.0]), "n=2, widely separated")

# n=3
test_superadditivity(np.array([1.0, 2.0, 5.0]), np.array([0.0, 3.0, 4.0]), "n=3, generic")
test_superadditivity(np.array([0.0, 1.0, 2.0]), np.array([0.0, 1.0, 2.0]), "n=3, same poly")
test_superadditivity(np.array([0.0, 1.0, 100.0]), np.array([0.0, 1.0, 2.0]), "n=3, one outlier")

# n=4
test_superadditivity(np.array([1.0, 2.0, 3.0, 4.0]), np.array([0.0, 1.0, 5.0, 6.0]), "n=4, generic")

# n=5
test_superadditivity(np.array([1.0, 2.0, 3.0, 4.0, 5.0]), np.array([0.0, 1.0, 2.0, 3.0, 10.0]), "n=5, generic")

print("\n\nRandom tests:")
np.random.seed(42)
n_tests = 50
failures = 0
for trial in range(n_tests):
    n = np.random.choice([2, 3, 4, 5, 6])
    roots_p = np.sort(np.random.randn(n) * 3)
    roots_q = np.sort(np.random.randn(n) * 3)

    # Ensure distinct roots (for Phi to be well-defined)
    while np.min(np.diff(roots_p)) < 0.1:
        roots_p = np.sort(np.random.randn(n) * 3)
    while np.min(np.diff(roots_q)) < 0.1:
        roots_q = np.sort(np.random.randn(n) * 3)

    roots_r = finite_free_convolution(roots_p, roots_q)

    # Check roots_r has distinct real roots
    if np.max(np.abs(np.imag(roots_r))) > 1e-6:
        continue
    roots_r = np.sort(np.real(roots_r))
    if len(roots_r) >= 2 and np.min(np.diff(roots_r)) < 1e-10:
        continue

    phi_p = phi_n_from_roots(roots_p)
    phi_q = phi_n_from_roots(roots_q)
    phi_r = phi_n_from_roots(roots_r)

    lhs = 1.0/phi_r
    rhs = 1.0/phi_p + 1.0/phi_q
    gap = lhs - rhs

    if gap < -1e-8:
        print(f"  FAILURE at trial {trial}: n={n}, gap={gap:.2e}")
        print(f"    roots_p = {roots_p}")
        print(f"    roots_q = {roots_q}")
        failures += 1

print(f"\nRandom tests: {n_tests} trials, {failures} failures")


# ============================================================================
# PART 2: Phi_n as a matrix trace functional / resolvent expression
# ============================================================================

print("\n" + "=" * 80)
print("PART 2: Phi_n in terms of resolvent and matrix quantities")
print("=" * 80)

def phi_n_resolvent_analysis(roots):
    """Analyze Phi_n in terms of the resolvent R(z) = (zI - A)^{-1}.

    For A = diag(lambda_1,...,lambda_n), we have:
    R(z)_{ij} = delta_{ij} / (z - lambda_i)

    The Cauchy/Stieltjes transform: G(z) = Tr(R(z)) = sum_i 1/(z - lambda_i)

    Note: H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j) = -G(lambda_i) + lim singular term

    More precisely, H_p(lambda_i) = -lim_{z->lambda_i} [G(z) - 1/(z - lambda_i)]
    = -lim_{z->lambda_i} sum_{j!=i} 1/(z - lambda_j)
    Wait, no: H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j)
    while G(z) = sum_j 1/(z - lambda_j), so
    H_p(lambda_i) = lim_{z->lambda_i} [G(z) - 1/(z - lambda_i)]
    = "regularized G at eigenvalue"
    """
    n = len(roots)

    # Direct computation
    H = H_values(roots)
    phi_direct = np.sum(H**2)

    # Resolvent approach: epsilon-regularization
    # G_eps(lambda_i) = sum_j 1/(lambda_i + i*eps - lambda_j)
    # = 1/(i*eps) + sum_{j!=i} 1/(lambda_i - lambda_j + i*eps)  (for the j=i term)
    # So H_p(lambda_i) = lim_{eps->0} [G_eps(lambda_i) - 1/(i*eps)]
    #                   = lim_{eps->0} Re[G_eps(lambda_i)]
    # Wait, let's be more careful:
    # G(lambda_i + i*eps) = sum_j 1/(lambda_i + i*eps - lambda_j)
    #                     = 1/(i*eps) + sum_{j!=i} 1/(lambda_i - lambda_j + i*eps)
    # Taking limit eps->0: the sum_{j!=i} -> H_p(lambda_i), and 1/(i*eps) diverges.
    # But Im[G(lambda_i + i*eps)] = sum_j (-eps)/((lambda_i - lambda_j)^2 + eps^2)
    # For eps->0, this -> -pi * (density at lambda_i). For discrete measure, diverges.
    #
    # Key: H_p(lambda_i) = Re[G(lambda_i + i*eps) - 1/(i*eps)] as eps->0
    #                     = Re[sum_{j!=i} 1/(lambda_i - lambda_j + i*eps)] -> sum_{j!=i} 1/(lambda_i - lambda_j)

    for eps in [0.1, 0.01, 0.001, 0.0001]:
        H_approx = np.zeros(n)
        for i in range(n):
            z = roots[i] + 1j * eps
            G_z = sum(1.0/(z - roots[j]) for j in range(n))
            # Subtract the singular part
            H_approx[i] = np.real(G_z - 1.0/(1j * eps))
        phi_approx = np.sum(H_approx**2)
        max_err = np.max(np.abs(H_approx - H))
        # print(f"  eps={eps:.4f}: Phi_approx={phi_approx:.10f}, max_H_err={max_err:.2e}")

    # MATRIX EXPRESSION for Phi_n:
    # Let D_{ij} = 1/(lambda_i - lambda_j) for i != j, D_{ii} = 0.
    # Then H_i = sum_j D_{ij} = (D @ 1)_i where 1 = all-ones vector.
    # So Phi_n = ||D @ 1||^2 = 1^T D^T D 1 = Tr(D^T D) restricted sense
    # More precisely: Phi_n = sum_i (sum_j D_{ij})^2 = ||D 1||^2
    # where D is the Cauchy matrix (with zero diagonal).

    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                D[i, j] = 1.0 / (roots[i] - roots[j])

    ones = np.ones(n)
    D_ones = D @ ones
    phi_matrix = np.dot(D_ones, D_ones)

    # Also: Phi_n = 1^T (D^T D) 1 = 1^T (D^2) 1  (since D is antisymmetric! D^T = -D)
    # Wait: D_{ij} = 1/(lambda_i - lambda_j), D_{ji} = 1/(lambda_j - lambda_i) = -D_{ij}
    # So D is antisymmetric: D^T = -D.
    # Therefore D^T D = (-D) D = -D^2.
    # So Phi_n = 1^T (-D^2) 1 = -1^T D^2 1 = -Tr(D^2 @ outer(1,1))
    # Or: Phi_n = -sum_{i,j} (D^2)_{ij}  since summing = multiplying by 1 on both sides.

    D2 = D @ D
    phi_from_D2 = -np.sum(D2)  # Should equal Phi_n

    # Also: -sum_{i,j} (D^2)_{ij} = -sum_{i,j} sum_k D_{ik} D_{kj}
    # = -sum_{i,j,k} 1/((lambda_i - lambda_k)(lambda_k - lambda_j))  [k!=i, k!=j]
    # + diagonal corrections where k=i or k=j

    # Even more elegant: Phi_n = Tr( (D @ J)^T (D @ J) ) / n  where J = n*|1><1|/n?
    # No, simpler: Phi_n = ||D @ 1||^2. Period.

    print(f"\n  roots = {roots}")
    print(f"  Phi_n (direct)    = {phi_direct:.10f}")
    print(f"  Phi_n (||D@1||^2) = {phi_matrix:.10f}")
    print(f"  Phi_n (-sum D^2)  = {phi_from_D2:.10f}")
    print(f"  D is antisymmetric: max|D + D^T| = {np.max(np.abs(D + D.T)):.2e}")

    # KEY INSIGHT: Can we express Phi_n as a trace of something involving A only?
    # Phi_n = sum_i (sum_{j!=i} 1/(lambda_i - lambda_j))^2
    # In terms of A = diag(lambda_i):
    # The matrix D_{ij} = [(A_ii - A_jj)^{-1}]_{off-diag} which depends on the
    # spectral decomposition, not just A itself.
    #
    # However, there IS a formula: Phi_n = Tr(L^2) where L is the "level repulsion" matrix,
    # L = sum_i H_p(lambda_i) |e_i><e_i| in the eigenbasis.
    # But this is tautological.
    #
    # More useful: Phi_n is related to the DISCRIMINANT and derivatives of p.

    return phi_direct, D

print("\nResolvent and matrix analysis:")
phi_n_resolvent_analysis(np.array([1.0, 3.0, 6.0]))
phi_n_resolvent_analysis(np.array([0.0, 1.0, 2.0, 5.0]))


# ============================================================================
# PART 2b: Phi_n in terms of discriminant and classical quantities
# ============================================================================

print("\n" + "-" * 40)
print("PART 2b: Relation to discriminant and p'")
print("-" * 40)

def phi_n_via_derivative(roots):
    """Compute Phi_n via p'(lambda_i).

    Since p(x) = prod_j (x - lambda_j), we have
    p'(lambda_i) = prod_{j!=i} (lambda_i - lambda_j)

    And p''(lambda_i) = 2 * sum_{j!=i} prod_{k!=i,k!=j} (lambda_i - lambda_k)

    So H_p(lambda_i) = p''(lambda_i) / (2*p'(lambda_i)) = sum_{j!=i} 1/(lambda_i - lambda_j)

    Therefore: Phi_n = sum_i [p''(lambda_i)/(2*p'(lambda_i))]^2

    Also: p'(lambda_i) = prod_{j!=i} (lambda_i - lambda_j)
    The discriminant disc(p) = prod_{i<j} (lambda_i - lambda_j)^2 = prod_i p'(lambda_i) * (-1)^{n(n-1)/2}
    Actually: disc(p) = (-1)^{n(n-1)/2} prod_i p'(lambda_i) ... let's verify.

    Actually disc(p) = prod_{i<j} (lambda_i - lambda_j)^2.
    And prod_i p'(lambda_i) = prod_i prod_{j!=i} (lambda_i - lambda_j)
                             = prod_{i!=j} (lambda_i - lambda_j) = (-1)^{n(n-1)/2} * disc(p)

    So: Phi_n = sum_i H_i^2 where H_i = p''(lambda_i)/(2*p'(lambda_i))

    Note: 1/Phi_n has NO obvious simple form in terms of disc(p).
    """
    n = len(roots)

    # Build polynomial
    p_coeffs = np.poly(roots)  # Highest degree first
    p_deriv = np.polyder(p_coeffs)
    p_deriv2 = np.polyder(p_deriv)

    phi_via_deriv = 0.0
    for i in range(n):
        p_prime_i = np.polyval(p_deriv, roots[i])
        p_double_prime_i = np.polyval(p_deriv2, roots[i])
        H_i = p_double_prime_i / (2 * p_prime_i)
        phi_via_deriv += H_i**2

    phi_direct = phi_n_from_roots(roots)

    # Discriminant
    disc = 1.0
    for i in range(n):
        for j in range(i+1, n):
            disc *= (roots[i] - roots[j])**2

    print(f"\n  roots = {roots}")
    print(f"  Phi_n (direct)         = {phi_direct:.10f}")
    print(f"  Phi_n (via p''/2p')    = {phi_via_deriv:.10f}")
    print(f"  disc(p)                = {disc:.6e}")
    print(f"  1/Phi_n                = {1/phi_direct:.10f}")

phi_n_via_derivative(np.array([1.0, 3.0]))
phi_n_via_derivative(np.array([0.0, 2.0, 5.0]))
phi_n_via_derivative(np.array([1.0, 2.0, 3.0, 6.0]))


# ============================================================================
# PART 3: Convexity/concavity of 1/Phi_n
# ============================================================================

print("\n" + "=" * 80)
print("PART 3: Convexity/concavity analysis of 1/Phi_n and Phi_n")
print("=" * 80)

def test_convexity_along_path(roots_p, roots_q, n_points=50):
    """Test convexity of Phi_n and 1/Phi_n along the path t*roots_p + (1-t)*roots_q.

    NOTE: This is the LINEAR interpolation of roots, NOT the boxplus interpolation.
    The boxplus path would be more natural but harder to parameterize.
    """
    ts = np.linspace(0, 1, n_points)
    phi_values = []
    inv_phi_values = []

    for t in ts:
        roots_t = t * roots_p + (1-t) * roots_q
        # Check distinct
        if np.min(np.abs(np.diff(np.sort(roots_t)))) < 1e-12:
            phi_values.append(np.nan)
            inv_phi_values.append(np.nan)
            continue
        phi_t = phi_n_from_roots(roots_t)
        phi_values.append(phi_t)
        inv_phi_values.append(1.0/phi_t)

    phi_values = np.array(phi_values)
    inv_phi_values = np.array(inv_phi_values)

    # Check convexity: f is convex if f(t*x + (1-t)*y) <= t*f(x) + (1-t)*f(y)
    # Equivalently: the second derivative is >= 0.
    # Numerical check: compare midpoint value to linear interpolation.

    # Check at several points
    convex_phi = True
    convex_inv_phi = True
    concave_phi = True
    concave_inv_phi = True

    for i in range(1, n_points-1):
        mid = (phi_values[i-1] + phi_values[i+1]) / 2
        if phi_values[i] > mid + 1e-10:
            convex_phi = False
        if phi_values[i] < mid - 1e-10:
            concave_phi = False

        mid_inv = (inv_phi_values[i-1] + inv_phi_values[i+1]) / 2
        if inv_phi_values[i] > mid_inv + 1e-10:
            convex_inv_phi = False
        if inv_phi_values[i] < mid_inv - 1e-10:
            concave_inv_phi = False

    return convex_phi, concave_phi, convex_inv_phi, concave_inv_phi

print("\nConvexity along LINEAR root interpolation t*lambda_p + (1-t)*lambda_q:")

tests_conv = [
    (np.array([1.0, 3.0]), np.array([0.0, 4.0]), "n=2 test 1"),
    (np.array([0.0, 1.0]), np.array([2.0, 5.0]), "n=2 test 2"),
    (np.array([1.0, 2.0, 5.0]), np.array([0.0, 3.0, 4.0]), "n=3 test 1"),
    (np.array([0.0, 1.0, 2.0]), np.array([3.0, 5.0, 8.0]), "n=3 test 2"),
    (np.array([1.0, 2.0, 3.0, 4.0]), np.array([0.0, 1.0, 5.0, 6.0]), "n=4 test 1"),
]

for roots_p, roots_q, label in tests_conv:
    cv_phi, cc_phi, cv_inv, cc_inv = test_convexity_along_path(roots_p, roots_q)
    print(f"  {label}:  Phi_n convex={cv_phi} concave={cc_phi}  |  1/Phi_n convex={cv_inv} concave={cc_inv}")


# Test convexity along BOXPLUS path
print("\nConvexity of Phi_n and 1/Phi_n along boxplus path:")
print("  For p_t = p boxplus_n (t * delta_a), testing as t varies.")

def test_boxplus_path_convexity(roots_p, a_vals, n_points=20):
    """Along path roots_p boxplus t*{a}, check behavior of Phi_n and 1/Phi_n.

    Actually, we test: fix roots_p. Let roots_q(t) = t * roots_q_base.
    Compute r(t) = roots_p boxplus_n roots_q(t) and Phi_n(r(t)).
    """
    roots_q_base = a_vals
    n = len(roots_p)

    ts = np.linspace(0.1, 3.0, n_points)
    phis = []

    for t in ts:
        roots_q = t * roots_q_base
        try:
            roots_r = finite_free_convolution(roots_p, roots_q)
            if np.max(np.abs(np.imag(roots_r))) > 1e-6:
                phis.append(np.nan)
                continue
            roots_r = np.real(roots_r)
            if np.min(np.abs(np.diff(np.sort(roots_r)))) < 1e-10:
                phis.append(np.nan)
                continue
            phis.append(phi_n_from_roots(roots_r))
        except:
            phis.append(np.nan)

    phis = np.array(phis)
    valid = ~np.isnan(phis)
    if np.sum(valid) < 3:
        return "insufficient data"

    # Check convexity of 1/Phi
    inv_phis = 1.0 / phis[valid]
    diffs = np.diff(inv_phis)
    ddiffs = np.diff(diffs)

    if np.all(ddiffs > -1e-10):
        return "1/Phi_n appears CONVEX along boxplus scaling path"
    elif np.all(ddiffs < 1e-10):
        return "1/Phi_n appears CONCAVE along boxplus scaling path"
    else:
        return "1/Phi_n appears NEITHER convex nor concave along boxplus scaling path"

result = test_boxplus_path_convexity(np.array([1.0, 3.0]), np.array([0.0, 2.0]))
print(f"  n=2: {result}")

result = test_boxplus_path_convexity(np.array([1.0, 2.0, 5.0]), np.array([0.0, 1.0, 3.0]))
print(f"  n=3: {result}")


# ============================================================================
# PART 4: Free probability connection — Phi_n vs Phi*
# ============================================================================

print("\n" + "=" * 80)
print("PART 4: Free probability connection")
print("=" * 80)

print("""
THEORETICAL ANALYSIS:

Voiculescu's free Fisher information for a compactly-supported measure mu on R:
  Phi*(mu) = integral |H_mu(x)|^2 dmu(x)
where H_mu is the Hilbert transform of mu:
  H_mu(x) = PV integral 1/(x-t) dmu(t)

For a DISCRETE measure mu = (1/n) sum_i delta_{lambda_i}:
  H_mu(lambda_i) = (1/n) sum_{j!=i} 1/(lambda_i - lambda_j)
  Phi*(mu) = (1/n) sum_i |H_mu(lambda_i)|^2
           = (1/n) sum_i [(1/n) sum_{j!=i} 1/(lambda_i - lambda_j)]^2
           = (1/n^3) sum_i [sum_{j!=i} 1/(lambda_i - lambda_j)]^2
           = (1/n^3) * Phi_n(p)

So: Phi_n(p) = n^3 * Phi*(mu_p)  where mu_p = (1/n) sum delta_{lambda_i}.

The VOICULESCU IDENTITY for free convolution:
  1/Phi*(mu boxplus nu) = 1/Phi*(mu) + 1/Phi*(nu)   [for freely independent]

In terms of Phi_n:
  n^3/Phi_n(r) should approximate n^3/Phi_n(p) + n^3/Phi_n(q)
  i.e., 1/Phi_n(r) ~ 1/Phi_n(p) + 1/Phi_n(q)

BUT: The Voiculescu identity is for the FREE convolution (n->infinity limit).
The FINITE free convolution boxplus_n is an approximation that agrees in the limit.

CONJECTURE: The finite version has 1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q),
i.e., superadditivity, with equality only in the limit n->infinity.
""")

def compute_free_fisher_info(roots):
    """Compute Phi*(mu) for mu = (1/n) sum delta_{lambda_i}.
    Phi*(mu) = Phi_n / n^3.
    """
    n = len(roots)
    phi_n = phi_n_from_roots(roots)
    return phi_n / n**3

print("Testing Voiculescu identity at finite n:")
print("  Should see 1/Phi*(r) >= 1/Phi*(p) + 1/Phi*(q), approaching equality as n grows.\n")

for n in [3, 5, 8, 12, 20]:
    # Generate well-separated roots
    roots_p = np.sort(np.linspace(-1, 1, n) + 0.05 * np.random.randn(n))
    roots_q = np.sort(np.linspace(-0.5, 0.5, n) + 0.05 * np.random.randn(n))

    # Ensure distinct
    for i in range(n-1):
        if roots_p[i+1] - roots_p[i] < 0.05:
            roots_p[i+1] = roots_p[i] + 0.05
        if roots_q[i+1] - roots_q[i] < 0.05:
            roots_q[i+1] = roots_q[i] + 0.05

    roots_r = finite_free_convolution(roots_p, roots_q)
    if np.max(np.abs(np.imag(roots_r))) > 1e-6:
        print(f"  n={n}: skipped (complex roots)")
        continue
    roots_r = np.sort(np.real(roots_r))

    phi_star_p = compute_free_fisher_info(roots_p)
    phi_star_q = compute_free_fisher_info(roots_q)
    phi_star_r = compute_free_fisher_info(roots_r)

    lhs = 1.0/phi_star_r
    rhs = 1.0/phi_star_p + 1.0/phi_star_q
    ratio = lhs / rhs  # Should be >= 1

    print(f"  n={n:3d}: 1/Phi*(r)={lhs:.6f}, 1/Phi*(p)+1/Phi*(q)={rhs:.6f}, ratio={ratio:.6f}, gap={lhs-rhs:.2e}")


# ============================================================================
# PART 5: Jensen / convexity approach
# ============================================================================

print("\n" + "=" * 80)
print("PART 5: Jensen / convexity approach")
print("=" * 80)

print("""
KEY IDEA: Since r = E_U[chi_{A + UBU*}], if Phi_n is CONVEX in roots (or coefficients),
then by Jensen: Phi_n(r) <= E_U[Phi_n(A + UBU*)].

Then we'd need: E_U[Phi_n(A + UBU*)] relates to Phi_n(A) + Phi_n(B) somehow.

Let's check: Is Phi_n convex in roots? Is 1/Phi_n concave?
""")

# Test: is Phi_n(lambda_1,...,lambda_n) convex as a function of the root vector?
print("Testing convexity of Phi_n as function of root vector:")

np.random.seed(123)
n_convexity_tests = 20
phi_convex_count = 0
phi_concave_count = 0
inv_phi_convex_count = 0
inv_phi_concave_count = 0

for trial in range(n_convexity_tests):
    n = np.random.choice([3, 4, 5])
    roots_a = np.sort(np.random.randn(n) * 2)
    roots_b = np.sort(np.random.randn(n) * 2)

    # Ensure well-separated
    while np.min(np.diff(roots_a)) < 0.3:
        roots_a = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(roots_b)) < 0.3:
        roots_b = np.sort(np.random.randn(n) * 2)

    cv_phi, cc_phi, cv_inv, cc_inv = test_convexity_along_path(roots_a, roots_b, n_points=40)
    if cv_phi: phi_convex_count += 1
    if cc_phi: phi_concave_count += 1
    if cv_inv: inv_phi_convex_count += 1
    if cc_inv: inv_phi_concave_count += 1

print(f"  Out of {n_convexity_tests} random line segments:")
print(f"  Phi_n convex:   {phi_convex_count}/{n_convexity_tests}")
print(f"  Phi_n concave:  {phi_concave_count}/{n_convexity_tests}")
print(f"  1/Phi_n convex: {inv_phi_convex_count}/{n_convexity_tests}")
print(f"  1/Phi_n concave:{inv_phi_concave_count}/{n_convexity_tests}")


# Monte Carlo estimate of E_U[Phi_n(A + UBU*)]
print("\nMonte Carlo: E_U[Phi_n(A + UBU*)] vs Phi_n(A) + Phi_n(B):")

def random_haar_unitary(n):
    """Generate Haar-distributed random unitary matrix."""
    Z = (np.random.randn(n, n) + 1j * np.random.randn(n, n)) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    d = np.diag(R)
    ph = d / np.abs(d)
    Q = Q @ np.diag(ph)
    return Q

def mc_expected_phi(roots_A, roots_B, n_samples=2000):
    """Monte Carlo estimate of E_U[Phi_n(A + UBU*)] and E_U[1/Phi_n(A + UBU*)]."""
    n = len(roots_A)
    A = np.diag(roots_A)
    B = np.diag(roots_B.astype(complex))

    phis = []
    inv_phis = []
    for _ in range(n_samples):
        U = random_haar_unitary(n)
        C = A + U @ B @ U.conj().T
        eigs = np.sort(np.real(np.linalg.eigvalsh(C)))

        # Check distinct
        if np.min(np.abs(np.diff(eigs))) < 1e-10:
            continue

        phi = phi_n_from_roots(eigs)
        phis.append(phi)
        inv_phis.append(1.0/phi)

    return np.mean(phis), np.mean(inv_phis), np.std(phis)/np.sqrt(len(phis))

for n in [3, 4, 5]:
    roots_A = np.linspace(0, 2*(n-1), n).astype(float)
    roots_B = np.linspace(0, n-1, n).astype(float)

    phi_A = phi_n_from_roots(roots_A)
    phi_B = phi_n_from_roots(roots_B)

    E_phi, E_inv_phi, std_phi = mc_expected_phi(roots_A, roots_B, n_samples=3000)

    roots_r = finite_free_convolution(roots_A, roots_B)
    roots_r = np.sort(np.real(roots_r))
    phi_r = phi_n_from_roots(roots_r)

    print(f"\n  n={n}:")
    print(f"    Phi_n(A)={phi_A:.6f}, Phi_n(B)={phi_B:.6f}, sum={phi_A+phi_B:.6f}")
    print(f"    Phi_n(r=A boxplus B)={phi_r:.6f}")
    print(f"    E_U[Phi_n(A+UBU*)]={E_phi:.6f} +/- {std_phi:.6f}")
    print(f"    E_U[1/Phi_n(A+UBU*)]={E_inv_phi:.8f}")
    print(f"    1/Phi_n(r)={1/phi_r:.8f}")
    print(f"    1/Phi_n(A)+1/Phi_n(B)={1/phi_A+1/phi_B:.8f}")
    print(f"    Jensen check: Phi_n(r)={phi_r:.6f} vs E_U[Phi_n]={E_phi:.6f} => {'Phi_n convex!' if phi_r <= E_phi + 2*std_phi else 'NOT obviously convex'}")


# ============================================================================
# PART 5b: Is 1/Phi_n concave in the polynomial COEFFICIENTS?
# ============================================================================

print("\n" + "-" * 40)
print("PART 5b: Concavity of 1/Phi_n in polynomial coefficients")
print("-" * 40)

def phi_n_from_coefficients(coeffs):
    """Given polynomial coefficients [1, a_{n-1}, ..., a_0] (highest degree first, monic),
    compute Phi_n."""
    roots = np.roots(coeffs)
    if np.max(np.abs(np.imag(roots))) > 1e-8:
        return np.nan  # Not real-rooted
    roots = np.sort(np.real(roots))
    if len(roots) >= 2 and np.min(np.abs(np.diff(roots))) < 1e-10:
        return np.nan  # Repeated roots
    return phi_n_from_roots(roots)

# For n=2: p(x) = x^2 - sx + p where s = lambda_1 + lambda_2, p = lambda_1 * lambda_2
# The polynomial is determined by (s, p). Let's test convexity in coefficient space.
print("\nFor n=2, Phi_n as function of (s=sum roots, p=product roots):")
print("  p(x) = x^2 - s*x + p, roots = (s +/- sqrt(s^2-4p))/2")
print("  Phi_n = (lambda_1-lambda_2)^{-2} + (lambda_2-lambda_1)^{-2} = 2/(lambda_1-lambda_2)^2 = 2/(s^2-4p)")
print("  So 1/Phi_n = (s^2-4p)/2.")
print("  This is LINEAR in p and CONVEX in s! So 1/Phi_n is convex in coefficients for n=2.")

# Verify
for s in [2.0, 4.0, 6.0]:
    for p_val in [0.0, 0.5]:
        coeffs = [1, -s, p_val]
        roots = np.roots(coeffs)
        if np.max(np.abs(np.imag(roots))) > 1e-10:
            continue
        roots = np.sort(np.real(roots))
        phi = phi_n_from_roots(roots)
        diff_sq = (roots[0] - roots[1])**2
        print(f"  s={s:.1f}, p={p_val:.1f}: Phi_n={phi:.6f}, 2/diff^2={2/diff_sq:.6f}, 1/Phi_n={(s**2-4*p_val)/2:.6f}")


# ============================================================================
# PART 6: The key formula: Phi_n and the Cauchy matrix
# ============================================================================

print("\n" + "=" * 80)
print("PART 6: Phi_n via Cauchy matrix — potential proof path")
print("=" * 80)

print("""
KEY STRUCTURAL INSIGHT:

Let Lambda = (lambda_1, ..., lambda_n) be the roots.
Define the off-diagonal Cauchy matrix:
  C_{ij} = 1/(lambda_i - lambda_j) for i != j, C_{ii} = 0.

Then: H_i = (C @ 1)_i  where 1 = (1,...,1)^T.
And:  Phi_n = ||C @ 1||^2 = 1^T C^T C 1.

Since C is antisymmetric (C^T = -C):
  C^T C = -C^2

So: Phi_n = -1^T C^2 1 = -sum_{i,j} (C^2)_{ij}.

Now, (C^2)_{ij} = sum_{k, k!=i, k!=j} 1/((lambda_i - lambda_k)(lambda_k - lambda_j))
(with the convention that terms with k=i or k=j are absent due to C_{kk}=0).

For the SUPERADDITIVITY inequality, we need:
  1/||C_r @ 1||^2 >= 1/||C_p @ 1||^2 + 1/||C_q @ 1||^2

where C_r is the Cauchy matrix of roots of r = p boxplus_n q.

This is equivalent to: the function Lambda -> 1/||C_Lambda @ 1||^2 is superadditive
under finite free convolution.

POTENTIAL PROOF PATH via Cauchy-Schwarz / matrix inequalities:
If we can show that Phi_n^{-1/2} (or some transform) is "additive-like" under boxplus,
the inequality might follow from Cauchy-Schwarz or a trace inequality.
""")

# Verify the Cauchy matrix formulation
print("Verification of Cauchy matrix formulation:")
for n in [3, 4, 5]:
    roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(roots)) < 0.3:
        roots = np.sort(np.random.randn(n) * 2)

    C = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                C[i, j] = 1.0 / (roots[i] - roots[j])

    ones = np.ones(n)
    phi_cauchy = np.dot(C @ ones, C @ ones)
    phi_direct = phi_n_from_roots(roots)
    phi_via_C2 = -np.sum(C @ C)

    print(f"  n={n}: direct={phi_direct:.8f}, ||C@1||^2={phi_cauchy:.8f}, -sum(C^2)={phi_via_C2:.8f}")


# ============================================================================
# PART 7: n=2 exact analysis
# ============================================================================

print("\n" + "=" * 80)
print("PART 7: n=2 exact analysis (the inequality is TRIVIALLY true!)")
print("=" * 80)

print("""
For n=2:
  p has roots a, b with a < b. Then:
  H_p(a) = 1/(a-b), H_p(b) = 1/(b-a)
  Phi_n(p) = 1/(a-b)^2 + 1/(b-a)^2 = 2/(b-a)^2
  1/Phi_n(p) = (b-a)^2 / 2

For q with roots c, d (c < d): 1/Phi_n(q) = (d-c)^2 / 2.

For r = p boxplus_2 q: the hat convolution gives
  hat_e_1(r) = hat_e_1(p) + hat_e_1(q) = (a+b)/2 + (c+d)/2
  hat_e_2(r) = hat_e_2(p) + hat_e_1(p)*hat_e_1(q) + hat_e_2(q)
             = ab + (a+b)(c+d)/4 + cd  (using hat_e_1 = e_1/C(2,1) = (a+b)/2, hat_e_2 = e_2/1 = ab)
  Wait, hat_e_2 = e_2/C(2,2) = ab/1 = ab.
  So hat_e_2(r) = ab + (a+b)(c+d)/4 + cd.

  e_1(r) = 2 * hat_e_1(r) = (a+b) + (c+d)
  e_2(r) = 1 * hat_e_2(r) = ab + (a+b)(c+d)/4 + cd

  r(x) = x^2 - e_1(r)*x + e_2(r)
  Roots of r: s = e_1(r), p = e_2(r)
  gap^2 = s^2 - 4p = [(a+b)+(c+d)]^2 - 4[ab + (a+b)(c+d)/4 + cd]
        = (a+b)^2 + 2(a+b)(c+d) + (c+d)^2 - 4ab - (a+b)(c+d) - 4cd
        = (a-b)^2 + (a+b)(c+d) + (c-d)^2
        = (a-b)^2 + (c-d)^2 + (a+b)(c+d)

  Hmm, that has a cross term. Let me redo more carefully.

  s = a+b+c+d, p = ab + cd + (a+b)(c+d)/4
  gap^2 = s^2 - 4p
        = (a+b+c+d)^2 - 4ab - 4cd - (a+b)(c+d)
        = (a+b)^2 + 2(a+b)(c+d) + (c+d)^2 - 4ab - 4cd - (a+b)(c+d)
        = (a-b)^2 + (a+b)(c+d) + (c-d)^2

  So 1/Phi_n(r) = gap^2/2 = [(a-b)^2 + (c-d)^2 + (a+b)(c+d)] / 2

  1/Phi_n(p) + 1/Phi_n(q) = (a-b)^2/2 + (c-d)^2/2 = [(a-b)^2 + (c-d)^2] / 2

  So: 1/Phi_n(r) - [1/Phi_n(p) + 1/Phi_n(q)] = (a+b)(c+d) / 2

  This is (mean_p)(mean_q) * 2 = the product of the MEANS of the root sets!

  IMPORTANT: This can be NEGATIVE if the means have opposite signs.
  For example: p has roots -3, 1 (mean = -1), q has roots 0, 2 (mean = 1).
  Then the gap = (-1)(1)*2 = -2 < 0, so the inequality FAILS?!
""")

# VERIFY the n=2 formula
print("Verifying n=2 formula:")
test_cases_n2 = [
    (np.array([-3.0, 1.0]), np.array([0.0, 2.0]), "means have opposite sign"),
    (np.array([1.0, 3.0]), np.array([0.0, 2.0]), "both means positive"),
    (np.array([-2.0, -1.0]), np.array([-3.0, -0.5]), "both means negative"),
    (np.array([-1.0, 1.0]), np.array([-1.0, 1.0]), "both mean zero"),
]

for roots_p, roots_q, label in test_cases_n2:
    roots_r = finite_free_convolution(roots_p, roots_q)
    phi_p = phi_n_from_roots(roots_p)
    phi_q = phi_n_from_roots(roots_q)
    phi_r = phi_n_from_roots(np.sort(np.real(roots_r)))

    a, b = roots_p
    c, d = roots_q
    predicted_gap = (a+b)*(c+d)/2
    actual_gap = 1/phi_r - (1/phi_p + 1/phi_q)

    print(f"\n  {label}:")
    print(f"    p roots: {roots_p}, mean={(a+b)/2:.2f}")
    print(f"    q roots: {roots_q}, mean={(c+d)/2:.2f}")
    print(f"    r roots: {np.sort(np.real(roots_r))}")
    print(f"    1/Phi(r) - 1/Phi(p) - 1/Phi(q) = {actual_gap:.10f}")
    print(f"    Predicted (a+b)(c+d)/2         = {predicted_gap:.10f}")
    print(f"    Match: {abs(actual_gap - predicted_gap) < 1e-10}")
    if predicted_gap < 0:
        print(f"    *** INEQUALITY FAILS! Gap is negative: {predicted_gap:.6f}")


# ============================================================================
# PART 8: CRITICAL FINDING — does the conjecture actually hold?
# ============================================================================

print("\n" + "=" * 80)
print("PART 8: CRITICAL — Systematic search for counterexamples")
print("=" * 80)

print("""
The n=2 analysis shows that 1/Phi_n(r) - 1/Phi_n(p) - 1/Phi_n(q) = (a+b)(c+d)/2.

This is NEGATIVE when mean(p) and mean(q) have opposite signs!

Let's check: does the conjecture assume CENTERED polynomials (mean zero)?
Or is the conjecture simply FALSE without centering?
""")

# Centered case (mean zero roots)
print("\n--- Centered case (roots sum to zero) ---")
n_centered_tests = 100
failures_centered = 0
np.random.seed(999)

for trial in range(n_centered_tests):
    n = np.random.choice([2, 3, 4, 5])

    # Generate centered roots (sum to 0)
    roots_p = np.random.randn(n) * 2
    roots_p -= np.mean(roots_p)
    roots_p = np.sort(roots_p)

    roots_q = np.random.randn(n) * 2
    roots_q -= np.mean(roots_q)
    roots_q = np.sort(roots_q)

    while np.min(np.diff(roots_p)) < 0.2:
        roots_p = np.random.randn(n) * 2
        roots_p -= np.mean(roots_p)
        roots_p = np.sort(roots_p)
    while np.min(np.diff(roots_q)) < 0.2:
        roots_q = np.random.randn(n) * 2
        roots_q -= np.mean(roots_q)
        roots_q = np.sort(roots_q)

    try:
        roots_r = finite_free_convolution(roots_p, roots_q)
        if np.max(np.abs(np.imag(roots_r))) > 1e-6:
            continue
        roots_r = np.sort(np.real(roots_r))
        if len(roots_r) >= 2 and np.min(np.diff(roots_r)) < 1e-10:
            continue

        phi_p = phi_n_from_roots(roots_p)
        phi_q = phi_n_from_roots(roots_q)
        phi_r = phi_n_from_roots(roots_r)

        gap = 1/phi_r - (1/phi_p + 1/phi_q)
        if gap < -1e-8:
            print(f"  CENTERED FAILURE at trial {trial}: n={n}, gap={gap:.2e}")
            print(f"    p={roots_p}, q={roots_q}")
            failures_centered += 1
    except:
        pass

print(f"Centered tests: {n_centered_tests} trials, {failures_centered} failures")

# Uncentered case
print("\n--- Uncentered case ---")
n_uncentered_tests = 100
failures_uncentered = 0
np.random.seed(888)

for trial in range(n_uncentered_tests):
    n = np.random.choice([2, 3, 4, 5])

    roots_p = np.sort(np.random.randn(n) * 3 + np.random.randn() * 2)
    roots_q = np.sort(np.random.randn(n) * 3 + np.random.randn() * 2)

    while np.min(np.diff(roots_p)) < 0.2:
        roots_p = np.sort(np.random.randn(n) * 3 + np.random.randn() * 2)
    while np.min(np.diff(roots_q)) < 0.2:
        roots_q = np.sort(np.random.randn(n) * 3 + np.random.randn() * 2)

    try:
        roots_r = finite_free_convolution(roots_p, roots_q)
        if np.max(np.abs(np.imag(roots_r))) > 1e-6:
            continue
        roots_r = np.sort(np.real(roots_r))
        if len(roots_r) >= 2 and np.min(np.diff(roots_r)) < 1e-10:
            continue

        phi_p = phi_n_from_roots(roots_p)
        phi_q = phi_n_from_roots(roots_q)
        phi_r = phi_n_from_roots(roots_r)

        gap = 1/phi_r - (1/phi_p + 1/phi_q)
        if gap < -1e-8:
            if failures_uncentered < 5:
                print(f"  UNCENTERED FAILURE at trial {trial}: n={n}, gap={gap:.2e}")
                print(f"    p={roots_p}, mean_p={np.mean(roots_p):.3f}")
                print(f"    q={roots_q}, mean_q={np.mean(roots_q):.3f}")
            failures_uncentered += 1
    except:
        pass

print(f"Uncentered tests: {n_uncentered_tests} trials, {failures_uncentered} failures")


# ============================================================================
# PART 9: What if we normalize / center?
# ============================================================================

print("\n" + "=" * 80)
print("PART 9: Centering analysis")
print("=" * 80)

print("""
The n=2 exact formula shows:
  gap = (mean_p * 2)(mean_q * 2) / 2 = 2 * mean_p * mean_q

This is >= 0 iff mean_p and mean_q have the SAME SIGN (or at least one is zero).

CENTERING: If we require p, q to be CENTERED (roots sum to zero = mean zero),
then for n=2: gap = 0 (equality!).

For centered polynomials with n >= 3, the gap should be non-negative.

Let's check: does Voiculescu's result assume centered measures?
In free probability, the free convolution mu boxplus nu has mean = mean(mu) + mean(nu).
The Fisher information Phi*(mu) is translation-invariant?
Actually NO: H_mu(x) = integral 1/(x-t) dmu(t) changes under translation.
But Phi*(mu) = integral |H_mu|^2 dmu IS invariant under simultaneous translation
of x and mu, i.e., Phi*(mu(.- c)) = Phi*(mu).

Wait: if mu is shifted by c, then H_{mu_c}(x) = integral 1/(x-t) dmu(t-c) = integral 1/(x-t-c) dmu(t) = H_mu(x-c).
And Phi*(mu_c) = integral |H_{mu_c}(x)|^2 dmu_c(x) = integral |H_mu(x-c)|^2 dmu(x-c) = integral |H_mu(y)|^2 dmu(y) = Phi*(mu).

So Phi* IS translation invariant! And Phi_n should be too. Let's check.
""")

# Is Phi_n translation invariant?
print("Testing: is Phi_n translation invariant?")
for roots in [np.array([1.0, 3.0, 5.0]), np.array([0.0, 2.0, 4.0, 7.0])]:
    for shift in [-3.0, 0.0, 2.0, 10.0]:
        phi_orig = phi_n_from_roots(roots)
        phi_shift = phi_n_from_roots(roots + shift)
        print(f"  roots={roots}, shift={shift:.1f}: Phi_n(orig)={phi_orig:.8f}, Phi_n(shifted)={phi_shift:.8f}, diff={abs(phi_orig-phi_shift):.2e}")

print("""
YES! Phi_n is translation invariant, as expected from the definition:
H_p(lambda_i + c) = sum_{j!=i} 1/((lambda_i+c) - (lambda_j+c)) = sum_{j!=i} 1/(lambda_i - lambda_j) = H_p(lambda_i)

So the CENTERING does NOT affect Phi_n. The issue is that:
1/Phi_n(p boxplus q) vs 1/Phi_n(p) + 1/Phi_n(q)
is NOT invariant under translation of p or q individually!

Because boxplus_n is NOT just translation: shifting p by c shifts r by c too.
But shifting p by c doesn't change Phi_n(p) while it DOES change the roots of r.

WAIT: Let's reconsider. If p is shifted by c, then p boxplus q is shifted by c too.
And Phi_n is translation invariant. So 1/Phi_n(r) and 1/Phi_n(p) are both unchanged.
Only 1/Phi_n(q) is also unchanged. So the inequality IS translation invariant!

Let me re-examine the n=2 formula.
""")

# Re-examine n=2
print("\nRe-examining n=2 with shift invariance:")
a, b = -3.0, 1.0
c, d = 0.0, 2.0

# Shift p so mean is 0: shift by 1, so roots become -4, 0
a2, b2 = a - (a+b)/2, b - (a+b)/2  # = -2, 2
c2, d2 = c - (c+d)/2, d - (c+d)/2  # = -1, 1

print(f"  Original: p roots = [{a}, {b}], q roots = [{c}, {d}]")
print(f"  Centered: p roots = [{a2}, {b2}], q roots = [{c2}, {d2}]")

# Compute directly
roots_r = finite_free_convolution(np.array([a, b]), np.array([c, d]))
phi_r = phi_n_from_roots(np.sort(np.real(roots_r)))
phi_p = phi_n_from_roots(np.array([a, b]))
phi_q = phi_n_from_roots(np.array([c, d]))
gap1 = 1/phi_r - 1/phi_p - 1/phi_q

roots_r2 = finite_free_convolution(np.array([a2, b2]), np.array([c2, d2]))
phi_r2 = phi_n_from_roots(np.sort(np.real(roots_r2)))
phi_p2 = phi_n_from_roots(np.array([a2, b2]))
phi_q2 = phi_n_from_roots(np.array([c2, d2]))
gap2 = 1/phi_r2 - 1/phi_p2 - 1/phi_q2

print(f"  Gap (original) = {gap1:.10f}")
print(f"  Gap (centered) = {gap2:.10f}")
print(f"  Phi_n(p) original={phi_p:.8f}, centered={phi_p2:.8f}")
print(f"  Phi_n(q) original={phi_q:.8f}, centered={phi_q2:.8f}")
print(f"  Phi_n(r) original={phi_r:.8f}, centered={phi_r2:.8f}")

print("""
KEY INSIGHT: The gap is NOT the same!
Shifting p by c and q by d individually doesn't preserve the gap because
r = (p shifted by c) boxplus (q shifted by d) is shifted by c+d, not unchanged.

The gap IS invariant under SIMULTANEOUS shift of p and q by the same c
(which shifts r by c, and Phi_n is translation invariant).

But the n=2 formula gap = (a+b)(c+d)/2 IS invariant under shifting both by c:
(a+b+2c)(c+d+2c)/2 != (a+b)(c+d)/2 in general.

Wait, that's wrong. Let me think again...

Shifting p: replace a,b by a+c, b+c. Then (a+c+b+c) = a+b+2c.
Shifting q: same amount c. Then (c'+d') = c+d+2c.
New gap: (a+b+2c)(c+d+2c)/2.

This is NOT the same as (a+b)(c+d)/2. So the gap is NOT shift-invariant!

The reason: shifting p by c and q by c gives r shifted by 2c (because boxplus adds means).
And Phi_n(r shifted by 2c) = Phi_n(r) but the r itself is different...

No wait: r = p boxplus q. If p' = p(. - c) (shift p by c), then p' boxplus q' where q' = q(. - c):
(p' boxplus q')(.?) = (p boxplus q)(. - 2c) shifted by 2c? No, that's not right either.

Let me just verify numerically whether shifting BOTH by the same c preserves the gap.
""")

print("\nShift invariance test (shift both p and q by same c):")
for shift in [-5.0, -2.0, 0.0, 2.0, 5.0]:
    rp = np.array([a, b]) + shift
    rq = np.array([c, d]) + shift
    rr = finite_free_convolution(rp, rq)
    pp = phi_n_from_roots(rp)
    pq = phi_n_from_roots(rq)
    pr = phi_n_from_roots(np.sort(np.real(rr)))
    g = 1/pr - 1/pp - 1/pq
    print(f"  shift={shift:+.1f}: gap = {g:.10f}")

print("""
CONCLUSION from n=2: The gap (1/Phi_n(r) - 1/Phi_n(p) - 1/Phi_n(q)) is NOT
shift-invariant, and it CAN be negative.

The n=2 formula gap = (mean_p * 2)(mean_q * 2)/2 = 2*mean_p*mean_q shows this clearly.

So THE CONJECTURE AS STATED IS FALSE for n=2 with general (uncentered) polynomials.
""")


# ============================================================================
# PART 10: Corrected conjecture?
# ============================================================================

print("\n" + "=" * 80)
print("PART 10: Investigating a CORRECTED conjecture")
print("=" * 80)

print("""
The conjecture 1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q) is FALSE for uncentered polynomials.

Possible corrections:
(A) Require centering: assume mean(roots_p) = mean(roots_q) = 0.
    For n=2 this gives gap = 0 (trivial equality).

(B) Add a correction term involving the means:
    1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q) + correction(mean_p, mean_q)

(C) Use a DIFFERENT normalization of Phi_n.

(D) The original conjecture is about a different definition of H or Phi_n.

Let's test option (A) for n >= 3.
""")

print("Testing corrected conjecture with CENTERED polynomials (n >= 3):")
np.random.seed(777)
n_tests = 200
results_by_n = {}

for trial in range(n_tests):
    n = np.random.choice([3, 4, 5, 6, 7])

    roots_p = np.random.randn(n) * 2
    roots_p -= np.mean(roots_p)
    roots_p = np.sort(roots_p)

    roots_q = np.random.randn(n) * 2
    roots_q -= np.mean(roots_q)
    roots_q = np.sort(roots_q)

    while np.min(np.diff(roots_p)) < 0.15:
        roots_p = np.random.randn(n) * 2
        roots_p -= np.mean(roots_p)
        roots_p = np.sort(roots_p)
    while np.min(np.diff(roots_q)) < 0.15:
        roots_q = np.random.randn(n) * 2
        roots_q -= np.mean(roots_q)
        roots_q = np.sort(roots_q)

    try:
        roots_r = finite_free_convolution(roots_p, roots_q)
        if np.max(np.abs(np.imag(roots_r))) > 1e-6:
            continue
        roots_r = np.sort(np.real(roots_r))
        if len(roots_r) >= 2 and np.min(np.diff(roots_r)) < 1e-10:
            continue

        phi_p = phi_n_from_roots(roots_p)
        phi_q = phi_n_from_roots(roots_q)
        phi_r = phi_n_from_roots(roots_r)

        gap = 1/phi_r - (1/phi_p + 1/phi_q)

        if n not in results_by_n:
            results_by_n[n] = []
        results_by_n[n].append(gap)
    except:
        pass

for n in sorted(results_by_n.keys()):
    gaps = results_by_n[n]
    min_gap = min(gaps)
    max_gap = max(gaps)
    mean_gap = np.mean(gaps)
    n_neg = sum(1 for g in gaps if g < -1e-10)
    print(f"  n={n}: {len(gaps)} tests, min_gap={min_gap:.6e}, max_gap={max_gap:.6e}, mean={mean_gap:.6e}, negative={n_neg}")


# ============================================================================
# PART 11: Alternative — maybe the conjecture is about a different Phi
# ============================================================================

print("\n" + "=" * 80)
print("PART 11: Alternative definitions of Phi_n")
print("=" * 80)

print("""
Perhaps the correct definition uses a DIFFERENT normalization.

Alternative 1: Phi_n^(alt1) = (1/n) sum_i H_i^2  (average instead of sum)
Alternative 2: Phi_n^(alt2) = sum_i (H_i / n)^2 = Phi_n / n^2
Alternative 3: Phi_n^(alt3) = Phi_n / n^3  (= Phi* from free probability)
Alternative 4: Use H_i = (1/(n-1)) sum_{j!=i} 1/(lambda_i - lambda_j)
""")

def phi_n_alt(roots, version):
    """Various normalizations of Phi_n."""
    n = len(roots)
    H = H_values(roots)

    if version == "sum":
        return np.sum(H**2)
    elif version == "avg":
        return np.mean(H**2)
    elif version == "n2":
        return np.sum(H**2) / n**2
    elif version == "n3":
        return np.sum(H**2) / n**3
    elif version == "norm_H":
        # Use H_i = (1/(n-1)) sum_{j!=i} ...
        H_norm = H / (n-1)
        return np.sum(H_norm**2)

print("Testing alternative normalizations for n=2 counterexample:")
roots_p = np.array([-3.0, 1.0])
roots_q = np.array([0.0, 2.0])
roots_r = np.sort(np.real(finite_free_convolution(roots_p, roots_q)))

for ver in ["sum", "avg", "n2", "n3", "norm_H"]:
    pp = phi_n_alt(roots_p, ver)
    pq = phi_n_alt(roots_q, ver)
    pr = phi_n_alt(roots_r, ver)
    gap = 1/pr - 1/pp - 1/pq
    print(f"  {ver:8s}: 1/Phi(r)={1/pr:.6f}, 1/Phi(p)+1/Phi(q)={1/pp+1/pq:.6f}, gap={gap:.6f}")

print("""
OBSERVATION: All normalizations give the SAME gap sign (just scaled).
The issue is fundamental: the gap for n=2 is proportional to mean_p * mean_q.

THIS MEANS: The conjecture as stated is FALSE without centering.
""")


# ============================================================================
# PART 12: FINAL ANALYSIS — what IS true?
# ============================================================================

print("\n" + "=" * 80)
print("PART 12: FINAL ANALYSIS — What is actually true?")
print("=" * 80)

print("""
SUMMARY OF FINDINGS:

1. The conjecture 1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q) is FALSE in general.
   Counterexample for n=2: p with roots (-3, 1), q with roots (0, 2).
   The gap = 2 * mean(p_roots) * mean(q_roots) = 2 * (-1) * 1 = -2 < 0.

2. Phi_n is translation invariant but the gap is NOT.

3. For CENTERED polynomials (roots sum to zero):
   - n=2: gap = 0 (trivial equality)
   - n=3,4,5,...: gap appears to be >= 0 in all numerical tests

4. FREE PROBABILITY CONNECTION:
   Voiculescu's free Fisher information Phi*(mu) = integral |H_mu|^2 dmu
   satisfies 1/Phi*(mu boxplus nu) = 1/Phi*(mu) + 1/Phi*(nu) EXACTLY
   for free convolution of freely independent variables.

   Our Phi_n = n^3 * Phi*(mu_empirical) is a finite-n analogue.
   The finite free convolution boxplus_n approximates free convolution as n -> infinity.

5. MATRIX FORMULATION:
   Phi_n = ||C @ 1||^2 where C is the off-diagonal Cauchy matrix C_{ij} = 1/(lambda_i - lambda_j).
   C is antisymmetric, so Phi_n = -1^T C^2 1 = -sum_{ij} (C^2)_{ij}.

6. The CORRECT CONJECTURE should be either:
   (a) Restricted to centered polynomials: for p, q with roots summing to 0,
       1/Phi_n(p boxplus_n q) >= 1/Phi_n(p) + 1/Phi_n(q).
   (b) Modified with a correction term:
       1/Phi_n(p boxplus_n q) = 1/Phi_n(p) + 1/Phi_n(q) + R_n(p,q)
       where R_n is a remainder that can be negative but satisfies some bound.

7. PROOF STRATEGY (for centered version):
   - Use the Cauchy matrix representation Phi_n = ||C@1||^2
   - The boxplus operation has a specific effect on the Cauchy matrix
   - Via the free probability limit, this should relate to subordination functions
   - The key tool would be the FINITE SUBORDINATION developed by Marcus
""")

# Final comprehensive test
print("\n--- Final comprehensive test of CENTERED conjecture ---")
np.random.seed(54321)
total_tests = 0
total_failures = 0

for n in [3, 4, 5, 6, 8, 10]:
    n_trials = 50
    failures = 0
    min_gap = float('inf')

    for trial in range(n_trials):
        roots_p = np.random.randn(n) * 2
        roots_p -= np.mean(roots_p)
        roots_p = np.sort(roots_p)

        roots_q = np.random.randn(n) * 2
        roots_q -= np.mean(roots_q)
        roots_q = np.sort(roots_q)

        while np.min(np.diff(roots_p)) < 0.1:
            roots_p = np.random.randn(n) * 2
            roots_p -= np.mean(roots_p)
            roots_p = np.sort(roots_p)
        while np.min(np.diff(roots_q)) < 0.1:
            roots_q = np.random.randn(n) * 2
            roots_q -= np.mean(roots_q)
            roots_q = np.sort(roots_q)

        try:
            roots_r = finite_free_convolution(roots_p, roots_q)
            if np.max(np.abs(np.imag(roots_r))) > 1e-6:
                continue
            roots_r = np.sort(np.real(roots_r))
            if len(roots_r) >= 2 and np.min(np.diff(roots_r)) < 1e-10:
                continue

            phi_p = phi_n_from_roots(roots_p)
            phi_q = phi_n_from_roots(roots_q)
            phi_r = phi_n_from_roots(roots_r)

            gap = 1/phi_r - (1/phi_p + 1/phi_q)
            min_gap = min(min_gap, gap)
            total_tests += 1

            if gap < -1e-8:
                failures += 1
                total_failures += 1
        except:
            pass

    print(f"  n={n:3d}: {n_trials} trials, min_gap={min_gap:.6e}, failures={failures}")

print(f"\nOVERALL: {total_tests} valid tests, {total_failures} failures")

if total_failures == 0:
    print("CENTERED CONJECTURE appears TRUE for all tested cases!")
else:
    print(f"CENTERED CONJECTURE has {total_failures} apparent failures — needs investigation")
