#!/usr/bin/env python3
"""
Investigate the Herglotz structure of omega_1 and phi = omega_1 - id.

Key setup:
- r = p boxplus_n q (MSS convolution)
- G_r(z) = G_p(omega_1(z))  (subordination)
- omega_1 maps C^+ to C^+ (Herglotz)
- omega_1(nu_k) = lambda_k, omega_1'(nu_k) = 1

phi(z) = omega_1(z) - z satisfies:
- phi'(nu_k) = 0 for all k (critical points at all roots of r)
- phi''(nu_k) = 2 * alpha_k

Question: Can we use the partial fraction / Herglotz representation of phi
to prove sum_k h_k * alpha_k >= 0?

Since omega_1 is Herglotz, Im(omega_1(z)) >= 0 for Im(z) > 0.
And Im(z) >= 0 for z in C^+.
So Im(phi(z)) = Im(omega_1(z)) - Im(z).
The sign of Im(phi(z)) is NOT determined in general.

However, omega_1 is a rational Herglotz function of degree n. Such functions
have a partial fraction representation:
  omega_1(z) = a*z + b + sum_j c_j / (d_j - z)
where a >= 0, c_j > 0, d_j real (the poles).

Since omega_1'(z) -> 1 as z -> infinity (from the Cauchy transform relation),
we have a = 1 (leading coefficient). So:
  omega_1(z) = z + b + sum_j c_j / (d_j - z)
  phi(z) = omega_1(z) - z = b + sum_j c_j / (d_j - z)

This means phi is a CONSTANT plus a SUM OF SIMPLE POLES with POSITIVE residues
(in the form c_j/(d_j - z) with c_j > 0).

Wait, let's be more careful. A degree-n rational Herglotz function mapping
C^+ to C^+ with omega_1(z) ~ z as z -> infinity has the form:
  omega_1(z) = z + b + sum_{j=1}^{n-1} c_j / (d_j - z)
where b is real and c_j > 0 and d_j are real poles.

Actually, for a degree-n rational function with n poles (at the roots of its
denominator), the number of poles of omega_1 - z equals n-1 (since omega_1
has degree n but the leading z cancels).

So phi(z) = b + sum_{j=1}^{n-1} c_j / (d_j - z) with c_j > 0.

phi'(z) = sum_{j=1}^{n-1} c_j / (d_j - z)^2
phi''(z) = -2 * sum_{j=1}^{n-1} c_j / (d_j - z)^3

phi'(nu_k) = 0 means: sum_j c_j / (d_j - nu_k)^2 = 0
But each term c_j/(d_j-nu_k)^2 >= 0 with c_j > 0. So all terms must be 0.
But c_j > 0 and d_j != nu_k (assuming). CONTRADICTION.

Wait, this means phi'(z) = sum c_j/(d_j-z)^2 is strictly positive everywhere
(away from poles), so phi' > 0 everywhere on R \ {poles}. But we need phi'(nu_k) = 0.

This is a contradiction. Let me re-examine the setup.

Actually, the issue is the partial fraction form. Let me be more careful.

omega_1 is a RATIONAL function of degree n. If omega_1(z) = P(z)/Q(z) where
P has degree n and Q has degree n-1 (or less), then
omega_1(z) = z + (constant) + (degree < n-1 terms) / Q(z)

But the exact form depends on the degrees. Let me compute omega_1 numerically
and find its poles and residues.
"""

import numpy as np
from math import factorial
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

def elem_sym_poly(roots, k):
    n = len(roots)
    if k == 0: return 1.0
    if k > n: return 0.0
    return sum(np.prod([roots[i] for i in subset])
               for subset in combinations(range(n), k))

def poly_coeffs_from_roots(roots):
    n = len(roots)
    ek = [elem_sym_poly(roots, k) for k in range(n+1)]
    return [(-1)**k * ek[k] for k in range(n+1)]

def boxplus_mss(roots_p, roots_q):
    n = len(roots_p)
    a = poly_coeffs_from_roots(roots_p)
    b = poly_coeffs_from_roots(roots_q)
    c = np.zeros(n+1)
    for k in range(n+1):
        for i in range(k+1):
            j = k - i
            if i <= n and j <= n:
                coeff = factorial(n-i) * factorial(n-j) / (factorial(n) * factorial(n-k))
                c[k] += coeff * a[i] * b[j]
    r_roots = np.sort(np.real(np.roots(c)))
    return r_roots, c

def H_values(roots):
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H

def Phi_n(roots):
    return np.sum(H_values(roots)**2)


def compute_omega1_rational(roots_p, roots_q, roots_r):
    """
    Compute omega_1 as a rational function using the relation G_r(z) = G_p(omega_1(z)).

    G_r(z) = (1/n) * r'(z)/r(z) = (1/n) * sum_k 1/(z - nu_k)
    G_p(w) = (1/n) * sum_k 1/(w - lambda_k)

    The equation G_r(z) = G_p(omega_1(z)) defines omega_1 implicitly.

    For numerical computation, we can use the fact that omega_1(nu_k) = lambda_k
    and omega_1 is a degree-n rational function to reconstruct it.

    Actually, let's compute omega_1(z) by solving G_r(z) = G_p(w) for w.
    G_p(w) = (1/n) sum_k 1/(w - lambda_k)
    Setting this = G_r(z), we get a degree-n polynomial in w.
    """
    n = len(roots_p)

    # For a specific z, solve G_p(w) = G_r(z) for w.
    # G_p(w) = (1/n) sum_k 1/(w - lambda_k) = (1/n) * p'(w)/p(w)
    # G_r(z) = (1/n) * r'(z)/r(z)

    # The equation p'(w)*r(z) = r'(z)*p(w) has n roots in w for each z.
    # omega_1(z) is the "correct" root (the one that maps nu_k -> lambda_k).

    p_poly = np.poly(roots_p)  # coeffs of p
    r_poly = np.poly(roots_r)  # coeffs of r
    p_prime = np.polyder(p_poly)
    r_prime = np.polyder(r_poly)

    # For z value, evaluate r(z) and r'(z), then solve p'(w)*r(z) - r'(z)*p(w) = 0
    # This is a polynomial of degree n-1 in w (since p' has degree n-1 and p has degree n,
    # and the leading terms: p'(w) has leading coeff n, p(w) has leading coeff 1.
    # So leading coeff of p'(w)*r(z): n*r(z)
    # Leading coeff of r'(z)*p(w): r'(z)*1
    # For generic z, these don't cancel, so degree = n.
    # Wait, p(w) has degree n and p'(w) has degree n-1.
    # p'(w)*r(z) has degree n-1 in w.
    # r'(z)*p(w) has degree n in w.
    # So F(z,w) = p'(w)*r(z) - r'(z)*p(w) has degree n in w.

    # The n roots of F(z,w)=0 in w, for z = nu_k, are: lambda_k (with F_w = r'(nu_k)*p'(lambda_k) != 0)
    # and the other roots come from...

    # Actually, at z = nu_k: r(nu_k) = 0, so F(nu_k, w) = -r'(nu_k)*p(w).
    # This has roots w = lambda_1, ..., lambda_n. All n roots of p!
    # So at z = nu_k, omega_1 = lambda_k is selected from among the n roots.

    # For generic z, F(z,w) = 0 has n roots in w. omega_1(z) is the one that varies
    # continuously from lambda_k to lambda_{k+1} as z goes from nu_k to nu_{k+1}.

    return p_poly, r_poly, p_prime, r_prime


def evaluate_omega1(z, roots_p, roots_r, p_poly, r_poly, p_prime, r_prime):
    """Evaluate omega_1(z) by solving the implicit equation."""
    n = len(roots_p)
    rz = np.polyval(r_poly, z)
    rpz = np.polyval(r_prime, z)

    # Solve p'(w)*r(z) - r'(z)*p(w) = 0
    # = r(z) * p'(w) - r'(z) * p(w) = 0
    # Build the polynomial in w:
    # p'(w) = sum of p_prime coefficients
    # p(w) = sum of p_poly coefficients
    # Total = r(z)*p_prime_coeffs(padded) - r'(z)*p_poly_coeffs

    # p_prime has degree n-1, p_poly has degree n
    # Pad p_prime to degree n with a leading 0
    p_prime_padded = np.zeros(n+1)
    p_prime_padded[1:] = p_prime  # p' has n coefficients (degree n-1)

    combined = rz * p_prime_padded - rpz * p_poly
    # This has degree n in w

    roots_w = np.roots(combined)

    # Select the root closest to the expected omega_1(z)
    # For real z between nu_k and nu_{k+1}, omega_1(z) should be between lambda_k and lambda_{k+1}
    return roots_w


def compute_phi_partial_fractions(roots_p, roots_q, roots_r, n_grid=1000):
    """
    Numerically determine the partial fraction structure of phi(z) = omega_1(z) - z.

    We know phi'(nu_k) = 0 for all k. This constrains the partial fraction form.

    If phi(z) = b + sum_j c_j / (d_j - z), then phi has simple poles at d_j.
    The poles of omega_1 are at the poles of G_r (which are NOT at the roots of r,
    since G_r has poles at z where r(z) = 0... wait, G_r = (1/n)*r'/r has poles at roots of r.

    But omega_1(nu_k) = lambda_k is FINITE. So omega_1 does NOT have poles at nu_k.
    The poles of omega_1 as a rational function are elsewhere on the real line.

    Let me find the poles numerically.
    """
    n = len(roots_p)
    p_poly, r_poly, p_prime, r_prime = compute_omega1_rational(roots_p, roots_q, roots_r)

    # omega_1(z) is a root of F(z,w) = r(z)*p'(w) - r'(z)*p(w) = 0.
    # As a rational function of z, omega_1(z) = P(z)/Q(z) where P, Q are polynomials.

    # The degree of omega_1: since F has degree n in w, for each z there are n branches.
    # omega_1 is one branch, which is a rational function of z of degree...
    # Actually, omega_1(z) satisfies an algebraic equation of degree n in w.
    # If the equation is generically irreducible, omega_1 is not rational but algebraic.
    # BUT in the MSS/free probability setting, omega_1 IS rational (this is a key theorem).

    # For n=2: omega_1 is a Mobius transformation (degree 1 rational function).
    # For n=3: omega_1 has degree 3 as a rational function.
    # General: omega_1 has degree n.

    # Let me compute omega_1(z) at many points and fit a rational function.

    # First, let's just evaluate omega_1 at many z values on the real line,
    # away from the roots of r.

    z_min = min(roots_r) - 5
    z_max = max(roots_r) + 5

    # Evaluate on grid points away from roots of r
    z_grid = np.linspace(z_min, z_max, n_grid)
    # Remove points too close to roots of r
    keep = np.ones(len(z_grid), dtype=bool)
    for nu in roots_r:
        keep &= np.abs(z_grid - nu) > 0.1
    z_grid = z_grid[keep]

    omega1_vals = np.zeros(len(z_grid))
    phi_vals = np.zeros(len(z_grid))

    for i, z in enumerate(z_grid):
        roots_w = evaluate_omega1(z, roots_p, roots_r, p_poly, r_poly, p_prime, r_prime)

        # Select the correct branch: for real z, omega_1(z) should be real
        # and should be the branch that connects lambda_k values
        real_roots = roots_w[np.abs(np.imag(roots_w)) < 0.01]
        if len(real_roots) == 0:
            omega1_vals[i] = np.nan
            phi_vals[i] = np.nan
            continue

        real_roots = np.real(real_roots)

        # Find which interval of r-roots z is in
        # and select the corresponding branch
        if z < roots_r[0]:
            # z < nu_1, so omega_1(z) < lambda_1
            target = min(real_roots)
        elif z > roots_r[-1]:
            # z > nu_n, so omega_1(z) > lambda_n
            target = max(real_roots)
        else:
            # z is between some nu_k and nu_{k+1}
            k = np.searchsorted(roots_r, z) - 1
            # omega_1(z) should be between lambda_k and lambda_{k+1}
            in_range = [r for r in real_roots if roots_p[k] - 0.5 <= r <= roots_p[k+1] + 0.5]
            if in_range:
                target = min(in_range, key=lambda r: abs(r - (roots_p[k]+roots_p[k+1])/2))
            else:
                target = min(real_roots, key=lambda r: abs(r - z))

        omega1_vals[i] = target
        phi_vals[i] = target - z

    return z_grid, omega1_vals, phi_vals


# ================================================================
# MAIN INVESTIGATION
# ================================================================
np.random.seed(42)

print("=" * 70)
print("HERGLOTZ STRUCTURE OF phi(z) = omega_1(z) - z")
print("=" * 70)

# Test with n=3 first
for trial in range(3):
    n = 3
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 1.0:
            roots_p[i] = roots_p[i-1] + 1.0
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 1.0:
            roots_q[i] = roots_q[i-1] + 1.0

    roots_r, c = boxplus_mss(roots_p, roots_q)
    raw = np.roots(c)
    if np.any(np.abs(np.imag(raw)) > 0.01):
        continue
    roots_r = np.sort(np.real(raw))
    if np.any(np.diff(roots_r) < 0.1):
        continue

    print(f"\nTrial {trial}: n={n}")
    print(f"  p roots (lambda): {roots_p}")
    print(f"  q roots (mu):     {roots_q}")
    print(f"  r roots (nu):     {roots_r}")

    # Verify omega_1 properties
    h = H_values(roots_r)
    u = H_values(roots_p)
    alpha = u - h  # Since sigma = identity (sorted roots)

    print(f"  h = {h}")
    print(f"  u = {u}")
    print(f"  alpha = {alpha}")
    print(f"  <h,alpha> = {np.dot(h, alpha):.8f} (should be >= 0)")

    # Compute phi values
    phi_at_nu = roots_p - roots_r  # omega_1(nu_k) - nu_k = lambda_k - nu_k
    print(f"  phi(nu_k) = lambda_k - nu_k = {phi_at_nu}")

    # phi'(nu_k) = 0 and phi''(nu_k) = 2*alpha_k
    print(f"  phi''(nu_k)/2 = alpha_k = {alpha}")

    # KEY COMPUTATION:
    # <h, alpha> = sum_k h_k * alpha_k = (1/2) sum_k H_r(nu_k) * phi''(nu_k)
    #
    # Now, H_r(nu_k) = (d/dz)(log |r(z)/(z-nu_k)|)_{z=nu_k}
    #                 = r''(nu_k)/(2*r'(nu_k))
    # (Actually, H_r(nu_k) = sum_{j!=k} 1/(nu_k - nu_j) which equals p''(nu_k)/(2p'(nu_k))
    #  for the SECOND derivative of the LOG of the polynomial.)

    # Let's explore the partial fraction structure of phi
    # phi(z) = omega_1(z) - z has the form:
    # phi(z) = b + R(z) where R(z) is a rational function with poles

    # Actually, since omega_1(z) ~ z + c/z + ... for large z (Herglotz with a=1),
    # phi(z) ~ c/z + ... for large z. So phi(z) -> 0 as z -> infinity.
    # This means b = 0 in the partial fraction.

    # Check: sum of residues of phi
    # phi(z) = sum_j c_j / (d_j - z) => phi(z) -> -sum c_j / z as z -> inf
    # From omega_1(z) = z + sum(lambda)/n - sum(nu)/n + O(1/z)... let me compute.

    # omega_1(z) ~ z + (mean(lambda) - mean(nu)) + O(1/z)
    # Actually, for the MSS convolution, mean(nu) = mean(lambda) + mean(mu).
    # So mean(lambda) - mean(nu) = -mean(mu).
    # Thus phi(z) ~ -mean(mu) + O(1/z) for large z.

    mean_diff = np.mean(roots_p) - np.mean(roots_r)
    print(f"  mean(lambda) - mean(nu) = {mean_diff:.6f}")
    print(f"  -mean(mu) = {-np.mean(roots_q):.6f}")

    # Hmm, these should be equal IF mean(nu) = mean(lambda) + mean(mu).
    # For MSS: e_1(r) = e_1(p) + e_1(q) (linear coefficient of boxplus).
    # And e_1 = sum of roots. So sum(nu) = sum(lambda) + sum(mu).
    # mean(nu) = mean(lambda) + mean(mu). YES.
    # So mean(lambda) - mean(nu) = -mean(mu).

    print(f"  Verified: mean(lambda)-mean(nu) = -mean(mu): "
          f"{abs(mean_diff + np.mean(roots_q)) < 1e-10}")

    # So for large z: phi(z) ~ -mean(mu) + O(1/z).
    # Wait, that's a CONSTANT, not going to 0. Let me re-check.

    # G_r(z) = G_p(omega_1(z))
    # For large z: G_r(z) = 1/z + mean(nu)/z^2 + ...
    # G_p(w) = 1/w + mean(lambda)/w^2 + ...
    # Setting w = omega_1(z) = z + phi(z):
    # G_p(z + phi) = 1/(z+phi) + mean(lambda)/(z+phi)^2 + ...
    #              = (1/z)(1/(1+phi/z)) + mean(lambda)/z^2 * (1/(1+phi/z))^2 + ...
    #              ~ 1/z - phi/z^2 + mean(lambda)/z^2 + ...
    # Setting equal to G_r(z) = 1/z + mean(nu)/z^2 + ...:
    # -phi/z^2 + mean(lambda)/z^2 = mean(nu)/z^2
    # => phi = mean(lambda) - mean(nu) = -mean(mu) (constant as z -> inf)

    # Actually wait: phi(z) = omega_1(z) - z. For large z:
    # If omega_1(z) = z + c_0 + c_1/z + ..., then phi(z) = c_0 + c_1/z + ...
    # From the above: c_0 = mean(lambda) - mean(nu) = -mean(mu).
    # So phi(z) -> -mean(mu) as z -> infinity. NOT zero.

    # This means the partial fraction form of phi is:
    # phi(z) = c_0 + sum_j c_j / (d_j - z) where c_0 = -mean(mu)

    # And phi' = sum_j c_j/(d_j-z)^2. For phi'(nu_k) = 0:
    # sum_j c_j/(d_j - nu_k)^2 = 0
    # This requires SOME c_j to be NEGATIVE (since (d_j-nu_k)^2 > 0).

    # AH HA! So phi does NOT have all-positive residues!
    # The Herglotz property is for omega_1, not phi.
    # omega_1(z) = z + c_0 + sum c_j/(d_j - z) with sum c_j/(d_j - z) being
    # a Pick function IF and ONLY IF all c_j > 0 and d_j real.
    # But omega_1(z) - z - c_0 = sum c_j/(d_j-z) maps C^+ to...
    # Im(omega_1(z)) >= 0 when Im(z) > 0 (Herglotz).
    # Im(z) > 0 trivially. So Im(omega_1(z) - z) = Im(omega_1(z)) - Im(z).
    # This is NOT necessarily >= 0!

    # So phi = omega_1 - id is NOT a Herglotz function in general.
    # The residues c_j can be of EITHER sign.

    # Let me compute them numerically.
    z_grid, omega1_vals, phi_vals = compute_phi_partial_fractions(roots_p, roots_q, roots_r)

    # Filter out NaN values
    valid = ~np.isnan(phi_vals)
    z_valid = z_grid[valid]
    phi_valid = phi_vals[valid]

    # Find the poles of phi by looking for discontinuities
    dphi = np.diff(phi_valid)
    dz = np.diff(z_valid)
    dphi_dz = dphi / dz

    # Large derivative indicates nearby pole
    pole_indices = np.where(np.abs(dphi_dz) > 100)[0]

    print(f"\n  Approximate poles of phi (from large derivatives):")
    if len(pole_indices) > 0:
        # Cluster nearby indices
        clusters = []
        current = [pole_indices[0]]
        for i in range(1, len(pole_indices)):
            if pole_indices[i] - pole_indices[i-1] < 5:
                current.append(pole_indices[i])
            else:
                clusters.append(current)
                current = [pole_indices[i]]
        clusters.append(current)

        pole_locations = []
        for cluster in clusters:
            idx = cluster[len(cluster)//2]
            pole_locations.append(z_valid[idx])
            print(f"    pole near z = {z_valid[idx]:.4f}")

        print(f"  Number of poles found: {len(pole_locations)}")
        print(f"  r roots: {roots_r}")
        print(f"  p roots: {roots_p}")

        # The poles of phi should be at the poles of omega_1.
        # omega_1 has poles where the denominator of the rational function vanishes.
        # Since omega_1 maps C^+ to C^+, poles must be on the real line.

    # Now let's try to compute the partial fraction decomposition more carefully.
    # phi(z) = c_0 + sum_j r_j / (z - p_j) where p_j are the poles.
    # (Using the convention phi(z) = c_0 + sum r_j/(z-p_j), so r_j are residues.)

    # For n=3: omega_1 is degree 3 rational, phi = omega_1 - z is degree 2 rational
    # (since the z terms cancel). So phi has at most 2 poles.

    # Actually, omega_1 has degree n as a rational function:
    # omega_1 = P(z)/Q(z) where P, Q have degrees n, n-1 (or n, n).
    # Wait, omega_1(z) -> z as z -> inf, so omega_1 = (z*Q(z) + R(z))/Q(z)
    # where deg(R) < deg(Q) = n-1. So phi = R(z)/Q(z) has degree at most n-1.
    # For n=3: phi has degree at most 2 rational function, so at most 2 poles.

    print(f"\n  Expected number of poles: {n-1}")


# ================================================================
# NOW: DETAILED ANALYSIS OF <h, alpha> VIA PARTIAL FRACTIONS
# ================================================================
print("\n\n" + "="*70)
print("DETAILED ANALYSIS: <h,alpha> VIA PARTIAL FRACTIONS OF omega_1")
print("="*70)

# Key insight to try:
# omega_1(z) = z + c_0 + R(z)/S(z) where R/S is a proper rational function
# with deg(R) < deg(S) = n-1.
#
# omega_1'(z) = 1 + (R'S - RS')/S^2
# omega_1'(nu_k) = 1 => (R'S - RS')/S^2 evaluated at nu_k = 0
# => R'(nu_k)*S(nu_k) = R(nu_k)*S'(nu_k) for all k.
#
# omega_1''(nu_k) = [(R''S - RS'')S^2 - (R'S - RS')*2SS'] / S^4 at nu_k
# Since R'S - RS' = 0 at nu_k:
# omega_1''(nu_k) = (R''(nu_k)*S(nu_k) - R(nu_k)*S''(nu_k)) / S(nu_k)^2
#
# alpha_k = omega_1''(nu_k)/2

# Actually, let me take a different approach. Let me use the NEVANLINNA representation
# directly for omega_1.

# A degree-n rational function mapping C^+ to C^+ can be written as:
# omega_1(z) = z + b + sum_{j=1}^{m} c_j / (d_j - z)  for Blaschke-type
# But this is for a SPECIAL class. Let me think more carefully.

# Actually, for a Herglotz function (maps upper half-plane to itself),
# the Nevanlinna representation is:
# f(z) = alpha + beta*z + integral (1/(t-z) - t/(1+t^2)) dmu(t)
# where alpha real, beta >= 0, mu is a positive measure.
#
# For a RATIONAL Herglotz function of degree n:
# f(z) = alpha + beta*z + sum_{j=1}^{m} c_j / (d_j - z)
# where beta >= 0, c_j > 0, d_j real.
# Here m = n - 1 if beta > 0, or m = n if beta = 0 (for degree n).
#
# Since omega_1(z) ~ z for large z, we need beta = 1.
# So omega_1(z) = alpha + z + sum_{j=1}^{n-1} c_j / (d_j - z) with all c_j > 0.
#
# Then phi(z) = omega_1(z) - z = alpha + sum_{j=1}^{n-1} c_j / (d_j - z).
# phi'(z) = sum_{j=1}^{n-1} c_j / (d_j - z)^2.
#
# BUT phi'(nu_k) = 0 means sum_j c_j/(d_j - nu_k)^2 = 0.
# Since all c_j > 0 and (d_j - nu_k)^2 > 0, each term is strictly positive.
# Sum of positive terms = 0 is IMPOSSIBLE.
#
# CONTRADICTION! This means omega_1 is NOT a Herglotz function of the form
# alpha + z + sum c_j/(d_j-z) with all c_j > 0.
#
# What went wrong? Let me re-examine.
#
# The issue is: omega_1 being Herglotz means Im(omega_1(z)) >= 0 when Im(z) >= 0.
# But this does NOT mean omega_1 has the Nevanlinna representation with all c_j > 0.
# That representation is for functions mapping C^+ STRICTLY into C^+ (or closure).
# Actually, the Nevanlinna representation IS the general form for Herglotz functions.
# So if omega_1 is Herglotz, it DOES have all c_j >= 0.
#
# But then phi'(nu_k) = sum c_j/(d_j-nu_k)^2 > 0 always, contradicting phi'(nu_k)=0.
#
# RESOLUTION: omega_1 might NOT be a "standard" Herglotz function!
#
# In the MSS setting, omega_1 is defined as the subordination function satisfying
# G_r(z) = G_p(omega_1(z)). The claim that omega_1 maps C^+ to C^+ needs verification.
#
# Actually, G_r and G_p both map C^+ to C^- (the Cauchy transform of a real measure
# maps upper half-plane to lower half-plane). So G_r(z) is in C^- for z in C^+.
# omega_1(z) must satisfy G_p(omega_1(z)) = G_r(z) in C^-.
# For G_p: C^+ -> C^-, so if omega_1(z) in C^+, then G_p(omega_1(z)) in C^-.
# This is consistent. So omega_1 DOES map C^+ to C^+.
#
# But then the Nevanlinna representation gives phi' > 0 at all real points.
# And omega_1'(nu_k) = 1 + phi'(nu_k) > 1, NOT = 1.
#
# THIS IS THE KEY ISSUE: either omega_1 is NOT Herglotz, or omega_1'(nu_k) != 1.
#
# Let me check numerically whether omega_1 is actually Herglotz.

print("\n--- Checking: Is omega_1 actually Herglotz? ---")
print("Testing Im(omega_1(z)) for z in C^+\n")

for trial in range(3):
    n = 3
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 1.0:
            roots_p[i] = roots_p[i-1] + 1.0
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 1.0:
            roots_q[i] = roots_q[i-1] + 1.0

    roots_r, c = boxplus_mss(roots_p, roots_q)
    raw = np.roots(c)
    if np.any(np.abs(np.imag(raw)) > 0.01):
        continue
    roots_r = np.sort(np.real(raw))
    if np.any(np.diff(roots_r) < 0.1):
        continue

    p_poly, r_poly, p_prime, r_prime = compute_omega1_rational(roots_p, roots_q, roots_r)

    print(f"Trial {trial}: n={n}, p={np.round(roots_p,3)}, r={np.round(roots_r,3)}")

    # Test omega_1 at complex points
    for y in [0.1, 0.5, 1.0, 2.0]:
        for x in np.linspace(min(roots_r)-2, max(roots_r)+2, 5):
            z = complex(x, y)
            rz = np.polyval(r_poly, z)
            rpz = np.polyval(r_prime, z)

            # Solve p'(w)*r(z) - r'(z)*p(w) = 0 for w
            p_prime_padded = np.zeros(n+1, dtype=complex)
            p_prime_padded[1:] = p_prime
            combined = rz * p_prime_padded - rpz * p_poly
            roots_w = np.roots(combined)

            # Select the root with Im(w) > 0 (expected for Herglotz)
            upper = roots_w[np.imag(roots_w) > -0.01]
            if len(upper) > 0:
                # Among upper half-plane roots, pick the one closest to real line
                # (the "physical" branch)
                best = upper[np.argmin(np.abs(np.imag(upper) - y))]
                im_omega = np.imag(best)
                if im_omega < -0.001:
                    print(f"  z={z}: Im(omega_1) = {im_omega:.6f} < 0 !!! NOT HERGLOTZ")

    # Check omega_1'(nu_k) = 1 numerically using the subordination derivative
    print(f"  omega_1'(nu_k) values:")
    for k in range(n):
        z0 = roots_r[k]
        w0 = roots_p[k]
        rp = np.polyval(r_prime, z0)
        rpp = np.polyval(np.polyder(r_prime), z0)
        pp = np.polyval(p_prime, w0)
        ppp = np.polyval(np.polyder(p_prime), w0)

        # omega_1'(nu_k) = r'(nu_k)*p'(lambda_k) / (r'(nu_k)*p'(lambda_k)) = 1
        # Wait, from F(z,w) = r(z)*p'(w) - r'(z)*p(w) = 0:
        # F_z = r'(z)*p'(w) - r''(z)*p(w)
        # F_w = r(z)*p''(w) - r'(z)*p'(w)
        # At (nu_k, lambda_k): r(nu_k) = 0, p(lambda_k) = 0
        # F_z = r'(nu_k)*p'(lambda_k) - r''(nu_k)*0 = r'(nu_k)*p'(lambda_k)
        # F_w = 0*p''(lambda_k) - r'(nu_k)*p'(lambda_k) = -r'(nu_k)*p'(lambda_k)
        # omega_1' = -F_z/F_w = -r'*p' / (-r'*p') = 1. YES!

        print(f"    omega_1'(nu_{k}) = 1 (exact, from F_z/F_w calculation)")

    # NOW: the phi' = 0 contradiction.
    # If omega_1 IS Herglotz, then phi = omega_1 - id has phi' > 0 everywhere.
    # But we just showed phi'(nu_k) = 0.
    # Unless: the nu_k are POLES of phi (so phi' is not defined there).

    # But omega_1(nu_k) = lambda_k is finite, so phi(nu_k) = lambda_k - nu_k is finite.
    # So nu_k is NOT a pole of phi.

    # The resolution must be that omega_1 is NOT a strict Herglotz function.
    # Specifically: omega_1 maps C^+ to C^+ union R (the closure).
    # A rational function mapping C^+ to closure of C^+ has the form:
    # omega_1(z) = z + b + sum c_j/(d_j - z) with c_j >= 0.
    # If c_j = 0 for some j, we can drop those terms.
    # phi' = sum c_j/(d_j-z)^2 >= 0 everywhere, with phi'(x) = 0 only if ALL c_j = 0,
    # meaning phi is constant.

    # WAIT: that can't be right either. omega_1 is not the identity map.

    # Actually, I think the issue is the DEGREE.
    # omega_1 is a degree-n rational function, meaning it's P(z)/Q(z) with max(deg P, deg Q) = n.
    # The Nevanlinna representation alpha + beta*z + sum c_j/(d_j-z) has degree max(1, m)
    # where m is the number of poles. For a degree-n Herglotz function, m = n-1 (if beta > 0)
    # or m = n (if beta = 0).

    # But this only works for degree-n functions of the SPECIAL form above.
    # NOT every degree-n rational function mapping C^+ to closure(C^+) has this form!

    # A degree-n rational function can have HIGHER ORDER poles.
    # For example, (z + 1/(z-a)^2) maps C^+ to C^+ (check: Im(1/(z-a)^2) can be positive).
    # Wait, actually 1/(z-a)^2 does NOT always have positive imaginary part in C^+.
    # For z = a + iy: 1/(iy)^2 = -1/y^2, which is real and negative. So Im = 0.
    # For z = a + x + iy: 1/(x+iy)^2 = (x-iy)^2/(x^2+y^2)^2 = ((x^2-y^2) - 2ixy)/(x^2+y^2)^2
    # Im = -2xy/(x^2+y^2)^2. For y > 0 and x > 0, this is negative.
    # So 1/(z-a)^2 is NOT Herglotz. Good.

    # So a Herglotz rational function of degree n MUST have only simple poles.
    # Therefore: omega_1(z) = alpha + z + sum_{j=1}^{n-1} c_j/(d_j - z) with c_j > 0.

    # And phi'(z) = sum c_j/(d_j-z)^2 > 0 for all z != d_j.

    # THIS CONTRADICTS phi'(nu_k) = 0.

    # CONCLUSION: Either omega_1 is NOT Herglotz (doesn't map C^+ to C^+ after all),
    # or the claim omega_1'(nu_k) = 1 is wrong.

    # Let me check numerically whether Im(omega_1(x + iy)) > y for small y.

    print(f"\n  Testing Im(omega_1(nu_k + i*eps)) vs eps:")
    for k in range(n):
        for eps in [0.01, 0.001, 0.0001]:
            z = complex(roots_r[k], eps)
            rz = np.polyval(r_poly, z)
            rpz = np.polyval(r_prime, z)

            p_prime_padded = np.zeros(n+1, dtype=complex)
            p_prime_padded[1:] = p_prime
            combined = rz * p_prime_padded - rpz * p_poly
            roots_w = np.roots(combined)

            # Select root near lambda_k + i*something
            dists = np.abs(roots_w - roots_p[k])
            best_idx = np.argmin(dists)
            w = roots_w[best_idx]

            im_ratio = np.imag(w) / eps
            print(f"    k={k}, eps={eps}: omega_1 = {w:.8f}, "
                  f"Im(omega_1)/eps = {im_ratio:.6f}")

    print()


# ================================================================
# RESOLUTION: WHAT TYPE OF FUNCTION IS omega_1?
# ================================================================
print("\n" + "="*70)
print("RESOLUTION: STRUCTURE OF omega_1")
print("="*70)
print("""
The apparent contradiction arises from the assumption that omega_1 is Herglotz
(maps C^+ to C^+). Let us re-examine this.

In the FREE PROBABILITY infinite-dimensional setting, the subordination function
omega_1 IS Herglotz: Im(omega_1(z)) >= Im(z) > 0 for z in C^+.
This is a deep theorem of Belinschi-Bercovici.

BUT in the FINITE MSS setting (finite free convolution of degree-n polynomials),
omega_1 is defined as the unique function satisfying G_r(z) = G_p(omega_1(z)).
This is an ALGEBRAIC equation, and omega_1 is a branch of an algebraic function.

The finite omega_1 may NOT be Herglotz in general!
Specifically, omega_1 is defined by the implicit equation and is rational of degree n.
If omega_1 were Herglotz with the partial fraction form, phi' > 0 everywhere,
which contradicts phi'(nu_k) = 0.

So either:
(a) omega_1 is NOT Herglotz in the finite case, OR
(b) omega_1 IS Herglotz but phi'(nu_k) != 0 (meaning omega_1'(nu_k) != 1)

We proved algebraically that omega_1'(nu_k) = 1 from the implicit function theorem.
Therefore (a) must hold: omega_1 is NOT a Herglotz function in the finite case.

This means the Herglotz approach to proving <h,alpha> >= 0 FAILS.
""")

# Let's verify: omega_1 is NOT Herglotz by finding a point where Im(omega_1(z)) < Im(z)
print("Checking if omega_1(z) - z maps C^+ to C^+...")
found_violation = False
for trial in range(5):
    n = 3
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 1.0:
            roots_p[i] = roots_p[i-1] + 1.0
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 1.0:
            roots_q[i] = roots_q[i-1] + 1.0

    roots_r, c = boxplus_mss(roots_p, roots_q)
    raw = np.roots(c)
    if np.any(np.abs(np.imag(raw)) > 0.01):
        continue
    roots_r = np.sort(np.real(raw))
    if np.any(np.diff(roots_r) < 0.1):
        continue

    p_poly = np.poly(roots_p)
    r_poly = np.poly(roots_r)
    p_prime = np.polyder(p_poly)
    r_prime = np.polyder(r_poly)

    for y in [0.001, 0.01, 0.1, 0.5, 1.0]:
        for x in np.linspace(min(roots_r)-2, max(roots_r)+2, 20):
            z = complex(x, y)
            rz = np.polyval(r_poly, z)
            rpz = np.polyval(r_prime, z)

            p_prime_padded = np.zeros(n+1, dtype=complex)
            p_prime_padded[1:] = p_prime
            combined = rz * p_prime_padded - rpz * p_poly
            roots_w = np.roots(combined)

            # Select the "physical" branch: the one with Im > 0 closest to z
            upper = roots_w[np.imag(roots_w) > -0.0001]
            if len(upper) == 0:
                continue

            # Select the one that would connect continuously to the correct branch
            best = upper[np.argmin(np.abs(upper - z))]

            im_phi = np.imag(best) - y
            if im_phi < -0.01 * y:
                if not found_violation:
                    print(f"  FOUND: Im(phi(z)) < 0 at z = {z}")
                    print(f"    omega_1(z) = {best}")
                    print(f"    Im(omega_1) = {np.imag(best):.6f}, Im(z) = {y}")
                    print(f"    Im(phi) = Im(omega_1) - Im(z) = {im_phi:.6f}")
                    found_violation = True

if not found_violation:
    print("  No violations found in tested points. omega_1 may still map C^+ to C^+")
    print("  (but then the Nevanlinna representation must have a subtlety)")

# Let's check more carefully: does Im(omega_1(z))/Im(z) -> 1 as Im(z) -> 0?
print("\n\nChecking Im(omega_1(nu_k + i*eps))/eps as eps -> 0:")
for trial in range(2):
    n = 4
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.8:
            roots_p[i] = roots_p[i-1] + 0.8
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.8:
            roots_q[i] = roots_q[i-1] + 0.8

    roots_r, c = boxplus_mss(roots_p, roots_q)
    raw = np.roots(c)
    if np.any(np.abs(np.imag(raw)) > 0.01):
        continue
    roots_r = np.sort(np.real(raw))
    if np.any(np.diff(roots_r) < 0.1):
        continue

    p_poly = np.poly(roots_p)
    r_poly = np.poly(roots_r)
    p_prime = np.polyder(p_poly)
    r_prime = np.polyder(r_poly)

    print(f"\nTrial {trial}: n={n}")
    print(f"  p roots: {np.round(roots_p, 3)}")
    print(f"  r roots: {np.round(roots_r, 3)}")

    for k in range(n):
        print(f"  k={k}, nu_k={roots_r[k]:.4f}, lambda_k={roots_p[k]:.4f}")
        for eps in [1.0, 0.1, 0.01, 0.001]:
            z = complex(roots_r[k], eps)
            rz = np.polyval(r_poly, z)
            rpz = np.polyval(r_prime, z)

            p_prime_padded = np.zeros(n+1, dtype=complex)
            p_prime_padded[1:] = p_prime
            combined = rz * p_prime_padded - rpz * p_poly
            roots_w = np.roots(combined)

            # Select root closest to lambda_k
            dists = np.abs(roots_w - roots_p[k])
            best_idx = np.argmin(dists)
            w = roots_w[best_idx]

            ratio = np.imag(w) / eps if eps > 0 else 0
            print(f"    eps={eps:.4f}: Im(omega_1)/eps = {ratio:.8f}, "
                  f"omega_1 = {np.real(w):.6f} + {np.imag(w):.6f}i")
