#!/usr/bin/env python3
"""
Numerical verification of contour integral / residue identities
for sum_k h_k * alpha_k.

We verify various candidate functions whose residues at nu_k
give h_k * alpha_k or related quantities.
"""

import numpy as np
from math import factorial
from itertools import permutations

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
        if roots[i] - roots[i-1] < 0.5:
            roots[i] = roots[i-1] + 0.5 + np.random.uniform(0, 0.5)
    return roots

def G_p(z, roots):
    """Cauchy transform: (1/n) * sum 1/(z - lambda_i)"""
    n = len(roots)
    return np.sum(1.0 / (z - roots)) / n

def nG_p(z, roots):
    """n * Cauchy transform: sum 1/(z - lambda_i)"""
    return np.sum(1.0 / (z - roots))

def solve_omega1(z, p_roots, r_roots):
    """
    Solve G_r(z) = G_p(w) for w numerically.
    Return the solution in C^+ (or closest to expected).
    """
    n = len(p_roots)
    p_poly = np.poly(p_roots)
    r_poly = np.poly(r_roots)
    dp = np.polyder(p_poly)
    dr = np.polyder(r_poly)

    # G_r(z) = G_p(w) means r'(z)*p(w)/(n*r(z)) = p'(w)*1/(n) ... actually
    # nG_r(z) = nG_p(w)
    # r'(z)/r(z) = p'(w)/p(w)
    # r'(z)*p(w) = p'(w)*r(z)

    rp_z = np.polyval(dr, z)
    r_z = np.polyval(r_poly, z)

    # rp_z * p(w) - r_z * p'(w) = 0
    # Build as polynomial in w
    combined = rp_z * p_poly.copy()
    padded_dp = np.zeros(n + 1)
    padded_dp[1:] = dp
    combined = combined - r_z * padded_dp

    w_roots = np.roots(combined)
    return w_roots

def get_omega1_branch(z, p_roots, r_roots):
    """Get the Herglotz branch of omega_1(z)."""
    w_roots = solve_omega1(z, p_roots, r_roots)
    # Select the one in C^+ (positive imaginary part)
    if np.imag(z) > 0:
        candidates = [w for w in w_roots if np.imag(w) > 0]
        if candidates:
            # Pick the one closest to z (at large |z|, omega_1(z) ~ z)
            return min(candidates, key=lambda w: abs(w - z))
    # For real z, pick by continuation
    return min(w_roots, key=lambda w: abs(np.imag(w)))

def omega1_at_root(nu_k, p_roots, r_roots):
    """
    Compute omega_1(nu_k) = which lambda_j it maps to.
    Use z = nu_k + i*epsilon and take limit.
    """
    eps_vals = [1e-4, 1e-5, 1e-6, 1e-7]
    results = []
    for eps in eps_vals:
        z = nu_k + 1j * eps
        w = get_omega1_branch(z, p_roots, r_roots)
        results.append(np.real(w))
    # Extrapolate
    return results[-1]

def omega1_derivative_at_root(nu_k, p_roots, r_roots, order=1):
    """
    Compute omega_1^(m)(nu_k) numerically via finite differences on C.
    Use z = nu_k + i*eps.
    """
    eps = 1e-5
    h = 1e-6
    z0 = nu_k + 1j * eps

    if order == 1:
        w_plus = get_omega1_branch(z0 + h, p_roots, r_roots)
        w_minus = get_omega1_branch(z0 - h, p_roots, r_roots)
        deriv = (w_plus - w_minus) / (2 * h)
    elif order == 2:
        w_plus = get_omega1_branch(z0 + h, p_roots, r_roots)
        w0 = get_omega1_branch(z0, p_roots, r_roots)
        w_minus = get_omega1_branch(z0 - h, p_roots, r_roots)
        deriv = (w_plus - 2*w0 + w_minus) / h**2
    return deriv


print("=" * 70)
print("VERIFICATION OF RESIDUE IDENTITIES FOR <h, alpha>")
print("=" * 70)

# Test for several n values
for n in [2, 3, 4]:
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")

    np.random.seed(100 + n)

    for trial in range(3):
        p_roots = generate_random_poly(n)
        q_roots = generate_random_poly(n)

        p_coeffs = roots_to_monic_coeffs(p_roots)
        q_coeffs = roots_to_monic_coeffs(q_roots)
        r_coeffs = boxplus_n(p_coeffs, q_coeffs, n)
        r_roots_raw = np.roots(r_coeffs)

        if np.max(np.abs(np.imag(r_roots_raw))) > 1e-6:
            print(f"  Trial {trial}: complex roots, skipping")
            continue
        r_roots = np.sort(np.real(r_roots_raw))
        if np.min(np.diff(r_roots)) < 1e-8:
            print(f"  Trial {trial}: repeated roots, skipping")
            continue

        H_p = H_values(p_roots)
        H_r = H_values(r_roots)

        h = H_r
        u = H_p  # identity pairing
        alpha = u - h

        ha_dot = np.dot(h, alpha)
        hu_dot = np.dot(h, u)
        phi_r = np.dot(h, h)
        phi_p = np.dot(u, u)

        print(f"\n  Trial {trial}:")
        print(f"    p_roots = {p_roots}")
        print(f"    r_roots = {r_roots}")
        print(f"    H_p = {H_p}")
        print(f"    H_r = {H_r}")
        print(f"    <h, alpha> = {ha_dot:.10f}")
        print(f"    <h, u>     = {hu_dot:.10f}")
        print(f"    ||h||^2    = {phi_r:.10f}")
        print(f"    <h,u> >= ||h||^2: {hu_dot >= phi_r - 1e-10}")

        # ============================================================
        # RESIDUE IDENTITY 1: f(z) = [nG_r(z)]^2 * omega_1'(z)
        #
        # Res_{z=nu_k} f(z) = 2*h_k + 2*alpha_k = 2*u_k
        # So sum of residues = 2*sum u_k = 2*sum H_p(lambda_k)
        #
        # By residue theorem on large contour:
        # sum_k Res f = -(1/2pi i) * integral on large circle
        # For large z: nG_r(z) ~ n/(nz) = 1/z, omega_1'(z) ~ 1
        # So f(z) ~ 1/z^2, integral on large circle -> 0
        # Hence sum Res = 0
        # This gives sum u_k = 0, which we know (sum of H_p = 0)
        # ============================================================

        sum_u = np.sum(u)
        sum_h = np.sum(H_r)
        print(f"\n    Identity check 1: sum(u_k) = {sum_u:.10f} (should be 0)")
        print(f"    Identity check 1: sum(h_k) = {sum_h:.10f} (should be 0)")

        # ============================================================
        # RESIDUE IDENTITY 2: g(z) = [nG_r(z)]^3 * omega_1'(z)
        #
        # Near z = nu_k:
        # [nG_r(z)]^3 = 1/(z-nu_k)^3 + 3*h_k/(z-nu_k)^2 + (3h_k^2 + S_k)/(z-nu_k) + ...
        # where S_k is the coefficient of (z-nu_k)^0 in the expansion of [nG_r]^2
        #
        # Actually [nG_r(z)]^2 = 1/(z-nu_k)^2 + 2h_k/(z-nu_k) + (h_k^2 + T_k) + ...
        # where T_k = sum_{j!=k} 1/(nu_k - nu_j)^2 - H_r(nu_k)^2 + H_r(nu_k)^2
        # Wait, let's be more careful.
        #
        # nG_r(z) = sum_j 1/(z - nu_j) = 1/(z-nu_k) + sum_{j!=k} 1/(z-nu_j)
        # Let f_k(z) = sum_{j!=k} 1/(z-nu_j)
        # f_k(nu_k) = H_r(nu_k) = h_k
        # f_k'(nu_k) = -sum_{j!=k} 1/(nu_k - nu_j)^2 =: -S_k  (negative definite!)
        #
        # [nG_r(z)]^2 = 1/(z-nu_k)^2 + 2*h_k/(z-nu_k) + (h_k^2 + 2*f_k'(nu_k)/1! * something)
        # Actually:
        # [1/(z-nu_k) + f_k(z)]^2 = 1/(z-nu_k)^2 + 2*f_k(z)/(z-nu_k) + f_k(z)^2
        # f_k(z) = h_k + f_k'(nu_k)*(z-nu_k) + ...
        # 2*f_k(z)/(z-nu_k) = 2*h_k/(z-nu_k) + 2*f_k'(nu_k) + ...
        # f_k(z)^2 = h_k^2 + ...
        # So [nG_r]^2 = 1/(z-nu_k)^2 + 2*h_k/(z-nu_k) + (h_k^2 + 2*f_k'(nu_k)) + ...
        # where f_k'(nu_k) = -sum_{j!=k} 1/(nu_k-nu_j)^2
        #
        # Now multiply by omega_1'(z) = 1 + 2*alpha_k*(z-nu_k) + ...
        # Res of [nG_r]^2 * omega_1' at nu_k = 2*h_k + 2*alpha_k (as computed in prompt)
        # = 2*u_k

        # For [nG_r]^3:
        # [nG_r]^3 = [nG_r]^2 * nG_r
        # = [1/(z-nu_k)^2 + 2h_k/(z-nu_k) + c0 + ...] * [1/(z-nu_k) + h_k + ...]
        # = 1/(z-nu_k)^3 + h_k/(z-nu_k)^2 + 2h_k/(z-nu_k)^2 + 2h_k^2/(z-nu_k)
        #   + c0/(z-nu_k) + h_k*[...]
        # = 1/(z-nu_k)^3 + 3h_k/(z-nu_k)^2 + (3h_k^2 + c0 - something...)/(z-nu_k) + ...
        #
        # Hmm, let me just compute the residues numerically instead.
        # ============================================================

        # Compute S_k = sum_{j!=k} 1/(nu_k - nu_j)^2
        S = np.zeros(n)
        for k in range(n):
            for j in range(n):
                if j != k:
                    S[k] += 1.0 / (r_roots[k] - r_roots[j])**2

        print(f"\n    S_k (sum 1/(nu_k-nu_j)^2) = {S}")

        # ============================================================
        # KEY APPROACH: Consider F(z) = d/dz [nG_r(z) * omega_1'(z)]
        # or equivalently look at what contour integral gives <h,alpha>
        #
        # We want sum_k h_k * alpha_k.
        #
        # Consider the function:
        #   psi(z) = [nG_r(z)]^2 * [omega_1'(z) - 1] / 2
        #
        # Near z = nu_k:
        # omega_1'(z) - 1 = 2*alpha_k*(z-nu_k) + (omega_1'''(nu_k)/2)*(z-nu_k)^2 + ...
        # [nG_r]^2 = 1/(z-nu_k)^2 + 2h_k/(z-nu_k) + (h_k^2 - S_k) + ...
        #    Wait: c0 = h_k^2 + 2*f_k'(nu_k) = h_k^2 - 2*S_k
        #
        # psi(z) = (1/2) * [1/(z-nu_k)^2 + 2h_k/(z-nu_k) + ...] * [2*alpha_k*(z-nu_k) + ...]
        # = (1/2) * [2*alpha_k/(z-nu_k) + 4*h_k*alpha_k + ...]
        # = alpha_k/(z-nu_k) + 2*h_k*alpha_k + ...
        #
        # So Res_{z=nu_k} psi(z) = alpha_k
        # And sum_k Res psi = sum_k alpha_k = sum_k (u_k - h_k) = 0 - 0 = 0
        # (since both sum u_k = 0 and sum h_k = 0)
        # This gives 0 = 0, not useful directly.
        # ============================================================

        sum_alpha = np.sum(alpha)
        print(f"\n    sum(alpha_k) = {sum_alpha:.10f} (should be 0)")

        # ============================================================
        # IDENTITY EXPLORATION: What function has Res = h_k * alpha_k ?
        #
        # We need a function whose residue at nu_k is h_k * alpha_k.
        #
        # Attempt: phi(z) = nG_r(z) * [omega_1'(z) - 1]
        # Near nu_k:
        # nG_r(z) = 1/(z-nu_k) + h_k + ...
        # omega_1'(z)-1 = 2*alpha_k*(z-nu_k) + ...
        # Product = 2*alpha_k + (h_k * 2*alpha_k + ...)(z-nu_k) + ...
        # This is REGULAR at nu_k! No residue.
        #
        # Attempt: phi(z) = [nG_r(z)]^2 * nG_r(z) * [omega_1'(z) - 1] / 2
        # = (1/2) * [nG_r]^3 * [omega_1'(z) - 1]
        # Near nu_k:
        # [nG_r]^3 = 1/(z-nu_k)^3 + 3h_k/(z-nu_k)^2 + c1/(z-nu_k) + ...
        # [omega_1'-1] = 2*alpha_k*(z-nu_k) + ...
        # Product = 2*alpha_k/(z-nu_k)^2 + 6*h_k*alpha_k/(z-nu_k) + ...
        # Res = (1/2) * 6*h_k*alpha_k = 3*h_k*alpha_k
        #
        # So sum Res phi = 3 * sum h_k*alpha_k = 3*<h,alpha>
        #
        # But what is the integral of phi on a large contour?
        # [nG_r]^3 ~ 1/z^3, omega_1'-1 ~ O(1/z^2) (since omega_1(z) = z + c + O(1/z))
        # So phi ~ 1/z^5, integral on large contour -> 0
        #
        # WAIT: This would give <h,alpha> = 0 which is WRONG!
        # Let me recalculate...
        # ============================================================

        # Let me recalculate more carefully.
        # omega_1(z) = z + a_0 + a_{-1}/z + a_{-2}/z^2 + ...
        # omega_1'(z) = 1 - a_{-1}/z^2 - 2*a_{-2}/z^3 - ...
        # omega_1'(z) - 1 = -a_{-1}/z^2 + O(1/z^3)
        #
        # [nG_r(z)]^3 = [1/z + O(1/z^2)]^3 = 1/z^3 + O(1/z^4)
        #
        # phi(z) = (1/2) * [1/z^3 + ...] * [-a_{-1}/z^2 + ...] = -a_{-1}/(2z^5) + ...
        # Integral on large circle: ~ 1/R^4 -> 0
        # So indeed sum Res = 0 on large contour, meaning 3*<h,alpha> = ???
        #
        # BUT: phi might have OTHER poles besides the nu_k!
        # omega_1'(z) could have poles (from the poles of omega_1).
        # omega_1(z) = z + c + sum_j m_j/(z - p_j) with m_j > 0, p_j real
        # omega_1'(z) = 1 - sum_j m_j/(z-p_j)^2
        # omega_1'(z) has double poles at the p_j.
        # So omega_1'(z) - 1 has double poles at p_j.
        #
        # phi(z) = (1/2)[nG_r(z)]^3 * [omega_1'(z)-1] has poles at:
        # 1. nu_k (from [nG_r]^3, order 3) * (omega_1'-1 has a zero of order 1
        #    at nu_k since omega_1'(nu_k)=1), net order 2, Res computed above
        # 2. p_j (poles of omega_1', double poles)
        #
        # So sum_{k} Res_{nu_k} phi + sum_{j} Res_{p_j} phi = 0
        # => 3*<h,alpha> = -sum_j Res_{p_j} phi
        #
        # This is INTERESTING! The sign of <h,alpha> depends on the
        # residues at the poles of omega_1.
        # ============================================================

        print(f"\n    Exploring phi(z) = (1/2)[nG_r(z)]^3 * [omega_1'(z)-1]")
        print(f"    Res at nu_k should be 3*h_k*alpha_k")

        # Verify numerically by computing residue via contour integral
        for k in range(n):
            nu_k = r_roots[k]
            # Small circle around nu_k
            R_circ = 1e-3
            N_pts = 1000
            theta = np.linspace(0, 2*np.pi, N_pts, endpoint=False)
            z_pts = nu_k + R_circ * np.exp(1j * theta)

            integrand_vals = np.zeros(N_pts, dtype=complex)
            for idx, z in enumerate(z_pts):
                nGr = nG_p(z, r_roots)
                # omega_1'(z) numerically
                h_step = 1e-7
                w_plus = get_omega1_branch(z + h_step, p_roots, r_roots)
                w_minus = get_omega1_branch(z - h_step, p_roots, r_roots)
                om1_prime = (w_plus - w_minus) / (2 * h_step)

                integrand_vals[idx] = 0.5 * nGr**3 * (om1_prime - 1)

            # Residue = (1/2pi i) * integral
            dz = 1j * R_circ * np.exp(1j * theta)
            residue = np.sum(integrand_vals * dz * (2*np.pi/N_pts)) / (2*np.pi*1j)
            expected = 3 * h[k] * alpha[k]

            print(f"    k={k}: Res_numerical = {np.real(residue):.6f}, "
                  f"3*h_k*alpha_k = {expected:.6f}, "
                  f"match = {abs(np.real(residue) - expected) < 0.01}")

print("\n" + "="*70)
print("APPROACH 2: Use [nG_r(z)]^2 * phi'(z) where phi(z) = omega_1(z) - z")
print("="*70)

# phi(z) = omega_1(z) - z
# phi(nu_k) = lambda_k - nu_k
# phi'(nu_k) = omega_1'(nu_k) - 1 = 0
# phi''(nu_k) = omega_1''(nu_k) = 2*alpha_k

# Consider: Psi(z) = phi'(z) / (phi(z) - c) for some parameter c?
# Or: Psi(z) = [nG_r(z)]' * phi(z) = -(sum 1/(z-nu_k)^2) * phi(z)

# Actually, let's think about this differently.
# <h,alpha> = (1/2) sum_k h_k * omega_1''(nu_k)
# = (1/2) sum_k h_k * phi''(nu_k)   where phi = omega_1 - id
#
# If we could write this as a contour integral involving phi...

# KEY IDEA:
# sum_k h_k * phi''(nu_k) = sum_k Res_{z=nu_k} [nG_r(z) * phi''(z)]
# because Res_{z=nu_k} nG_r(z) = 1 and phi''(z) is regular at nu_k
# (since phi has no poles at nu_k if omega_1 has no poles there)
#
# Wait: omega_1(z) does have poles! phi(z) = omega_1(z) - z has the same poles.
# But phi''(z) = omega_1''(z) is regular at nu_k if omega_1 is analytic there.
# Is omega_1 analytic at nu_k? YES - omega_1 is a rational function with poles
# at the p_j (distinct from the nu_k).
#
# So sum_k h_k * phi''(nu_k) = sum_k Res_{z=nu_k} [nG_r(z) * phi''(z)]
# = (1/2pi i) oint nG_r(z) * phi''(z) dz - sum_j Res_{p_j} [nG_r(z)*phi''(z)]
#
# On large contour: nG_r(z) ~ 1/z, phi''(z) = omega_1''(z) ~ O(1/z^3)
# (since omega_1(z) = z + c + a/z + ..., omega_1''(z) = 2a/z^3 + ...)
# So nG_r * phi'' ~ 1/z^4, integral -> 0
#
# => sum_k h_k*phi''(nu_k) = -sum_j Res_{p_j} [nG_r(z)*phi''(z)]
#
# At a pole p_j of omega_1 (with residue m_j > 0):
# phi(z) = m_j/(z-p_j) + (regular)
# phi'(z) = -m_j/(z-p_j)^2 + (regular)
# phi''(z) = 2*m_j/(z-p_j)^3 + (regular)
#
# nG_r(z) is regular at p_j (since p_j are not roots of r).
# So Res_{p_j} [nG_r(z)*phi''(z)] = 0 (phi'' has pole of order 3, nG_r regular,
# product has pole of order 3 but no 1/(z-p_j) term unless...)
#
# Wait: nG_r(z) = nG_r(p_j) + nG_r'(p_j)*(z-p_j) + (1/2)nG_r''(p_j)*(z-p_j)^2 + ...
# phi''(z) = 2m_j/(z-p_j)^3 + phi'''(p_j??)...
# Actually phi''(z) = sum_j 2m_j/(z-p_j)^3
# At pole p_j: phi''(z) = 2m_j/(z-p_j)^3 + (terms from other poles, regular at p_j)
# Let R_j = sum_{l!=j} 2m_l/(p_j-p_l)^3 + correction...
# phi''(z) = 2m_j/(z-p_j)^3 + R_j''(p_j) + ... where R_j is regular
#
# Product: nG_r(z)*phi''(z) near p_j:
# = [nG_r(p_j) + nG_r'(p_j)(z-p_j) + (1/2)nG_r''(p_j)(z-p_j)^2 + ...]
#   * [2m_j/(z-p_j)^3 + C_j/(z-p_j)^2 + D_j/(z-p_j) + ...]
#
# Wait, phi''(z) has the Laurent expansion at p_j:
# phi''(z) = 2m_j/(z-p_j)^3 + [regular part]
# since only the j-th pole contributes a singularity.
#
# So: nG_r*phi'' = 2m_j*nG_r(p_j)/(z-p_j)^3 + [2m_j*nG_r'(p_j) + ...]/(z-p_j)^2
#                  + [m_j*nG_r''(p_j) + ...]/(z-p_j) + ...
#
# Res_{p_j} = m_j * nG_r''(p_j) + ... (need to be more careful)
# Actually just the 1/(z-p_j) coefficient:
# from 2m_j/(z-p_j)^3 * (1/2)*nG_r''(p_j)*(z-p_j)^2 = m_j*nG_r''(p_j)/(z-p_j)
# PLUS from the regular part of phi'' at p_j:
# We need the full Laurent expansion. But the key point is:

# Res_{p_j}[nG_r(z)*phi''(z)] = m_j * nG_r''(p_j)  (leading contribution)
# where nG_r''(p_j) = 2*sum_k 1/(p_j - nu_k)^3

# So <h,alpha> = (1/2)*sum_k h_k*phi''(nu_k)
#              = -(1/2)*sum_j Res_{p_j}[nG_r*phi'']
#              = -(1/2)*sum_j m_j * nG_r''(p_j) * (leading order)

# nG_r''(p_j) = sum_k 2/(p_j-nu_k)^3
# This has MIXED SIGNS depending on whether p_j is above or below nu_k.

# Hmm, this doesn't immediately give a sign. Let me try a different function.

print("\n" + "="*70)
print("APPROACH 3: Integration by parts")
print("Consider sum_k h_k * alpha_k = (1/2) sum_k [nG_r residue] * omega_1''(nu_k)")
print("="*70)

# Use the identity: for f(z) with simple poles at nu_k with residue 1,
# and g(z) analytic at nu_k:
# sum_k g(nu_k) = sum_k Res_{z=nu_k} [f(z)*g(z)]
# = (1/2pi i) oint f(z)*g(z) dz - (residues at other poles)
#
# Here f(z) = nG_r(z), g(z) = omega_1''(z)/2
#
# Actually this doesn't directly give h_k * alpha_k because we'd get
# 1 * alpha_k, not h_k * alpha_k.
#
# For h_k * alpha_k we need f(z) with Res = h_k at z = nu_k.
# What function has Res = h_k at nu_k?
#
# Consider: [nG_r(z)]^2 / 2
# Res_{nu_k} [nG_r(z)]^2 / 2 = h_k (from the Laurent expansion)
# [nG_r(z)]^2 = 1/(z-nu_k)^2 + 2h_k/(z-nu_k) + ...
# Res = 2*h_k, so [nG_r]^2/2 has residue h_k. YES!
#
# So: sum_k h_k * alpha_k = sum_k Res_{nu_k} { [nG_r(z)]^2/2 * alpha_k }
# But alpha_k is a constant, not a function of z!
# We need: sum_k h_k * f(nu_k) where f(nu_k) = alpha_k = omega_1''(nu_k)/2
#
# So: sum_k h_k * omega_1''(nu_k)/2
#   = sum_k Res_{nu_k} [nG_r(z)]^2/2 * omega_1''(z)/2
#   + "cross terms"
#
# No! The residue of [nG_r]^2/2 * omega_1''/2 at nu_k is NOT h_k * alpha_k.
# Because omega_1''(z) varies, and the residue picks up the value at the pole.
# Res_{nu_k} [nG_r(z)]^2/2 * omega_1''(z)/2
# = (1/2)*(1/2)* Res_{nu_k} [1/(z-nu_k)^2 + 2h_k/(z-nu_k) + ...]*[omega_1''(nu_k)+...]
# = (1/4) * [2*h_k * omega_1''(nu_k) + d/dz|_{nu_k} omega_1''(z) * 1]
# Hmm wait. Let me be careful:
# [nG_r(z)]^2 = 1/(z-nu_k)^2 + 2h_k/(z-nu_k) + c_0 + ...
# omega_1''(z) = 2*alpha_k + omega_1'''(nu_k)*(z-nu_k) + ...
# Product = [1/(z-nu_k)^2]*2alpha_k + [2h_k/(z-nu_k)]*2alpha_k
#           + [1/(z-nu_k)^2]*omega_1'''(nu_k)*(z-nu_k) + ...
# = 2alpha_k/(z-nu_k)^2 + [4h_k*alpha_k + omega_1'''(nu_k)]/(z-nu_k) + ...
# Res = 4*h_k*alpha_k + omega_1'''(nu_k)
# So (1/4) * Res = h_k*alpha_k + omega_1'''(nu_k)/4
# Not quite what we want -- there's an extra omega_1'''(nu_k)/4 term.

print("\nKey finding: direct residue approach picks up omega_1''' correction terms.")
print("These extra terms complicate the contour integral approach.")

print("\n" + "="*70)
print("APPROACH 4: The 'right' function for <h,alpha>")
print("="*70)

# Let me try: Theta(z) = [nG_r(z)]^2 * [nG_p(omega_1(z)) - nG_r(z)]
#
# Near z = nu_k:
# nG_p(omega_1(z)) = nG_r(z) (by subordination!)
# So Theta(z) = 0 identically!
# This is useless.
#
# Try: Theta(z) = nG_r(z) * [nG_p(omega_1(z))]'
# = nG_r(z) * [nG_p]'(omega_1(z)) * omega_1'(z)
# = nG_r(z) * (-sum_i 1/(omega_1(z)-lambda_i)^2) * omega_1'(z)
# Hmm.
#
# Better approach: Consider the DERIVATIVE of the subordination identity.
# G_r(z) = G_p(omega_1(z))
# => G_r'(z) = G_p'(omega_1(z)) * omega_1'(z)
# => n*G_r'(z) = n*G_p'(omega_1(z)) * omega_1'(z)
#
# Now: n*G_r'(z) = -sum_k 1/(z-nu_k)^2
# At z = nu_k this diverges (double pole).
#
# Consider ratio: G_r'(z) / G_r(z)^2
# = G_p'(omega_1(z)) * omega_1'(z) / G_p(omega_1(z))^2
# = omega_1'(z) * [G_p'/G_p^2](omega_1(z))
#
# At z = nu_k (pole of G_r):
# G_r(z) = 1/(n(z-nu_k)) + h_k/n + ...
# G_r'(z) = -1/(n(z-nu_k)^2) + ...
# G_r'/G_r^2 = [-1/(n(z-nu_k)^2)] / [1/(n^2(z-nu_k)^2)] * (1 + ...)
# = -n * (1 + ...)
# So the ratio is regular at nu_k.

# Let me try yet another approach:
#
# We know that sum_k h_k * alpha_k = <h,u> - ||h||^2
# = sum_k H_r(nu_k) * H_p(lambda_k) - sum_k H_r(nu_k)^2
#
# The second term is Phi_n(r).
# The first term sum_k H_r(nu_k) * H_p(lambda_k) is what we need to express.
#
# Using residues of nG_r at nu_k:
# sum_k H_r(nu_k) * H_p(lambda_k)
# = sum_k Res_{nu_k}[nG_r(z)] * H_p(omega_1(nu_k))
# but the residue of nG_r at nu_k is 1, not H_r(nu_k)!
#
# For H_r(nu_k) as residue:
# Res_{nu_k} [nG_r(z)]^2 / 2 ... no, [nG_r]^2 has Res = 2*h_k
#
# Or: define F(z) = (1/2) * d/dz [(nG_r(z))^2] = -nG_r(z) * sum_k 1/(z-nu_k)^2
# Hmm, d/dz [nG_r] = -sum 1/(z-nu_k)^2, so d/dz [nG_r^2] = 2*nG_r * d/dz[nG_r]
# The residue of [nG_r]^2 at nu_k is 2*h_k, but this is the integral quantity.

# CLEAN APPROACH:
# We want sum_k h_k * g(nu_k) where g(nu_k) = alpha_k.
#
# Lagrange interpolation:
# g(z) = sum_k g(nu_k) * prod_{j!=k} (z-nu_j)/(nu_k-nu_j)
# = sum_k alpha_k * r(z)/((z-nu_k)*r'(nu_k))
#
# Then: sum_k h_k * g(nu_k) = sum_k h_k * alpha_k = <h,alpha>
#
# Can we express g(z) = omega_1''(z)/2 in terms of r and p?
# If so, we can evaluate sum_k h_k * alpha_k via a contour integral
# = (1/2pi i) oint [nG_r(z)]^2 / 2 * [omega_1''(z)/2] * [some correction] dz

# This is getting circular. Let me try the ENERGY functional approach.

print("\n" + "="*70)
print("APPROACH 5: Energy / generating function")
print("Consider E(t) = sum_k [H_{r_t}(nu_k(t))]^2 = Phi(r_t)")
print("where r_t = p boxplus_n (q scaled by t)")
print("="*70)

# The idea: if we can show d/dt Phi(r_t) <= 0 at t=1, that relates to <h,alpha>.
# But this is a different kind of derivative (parametric).

# Actually, let's try the simplest possible contour integral identity:
#
# CLAIM: <h,alpha> = (1/4pi i) oint [r''(z)/(r'(z))]' * [omega_1(z) - z] dz
# where the contour encloses all roots of r.
#
# r''(z)/(2r'(z)) restricted to the k-th root gives H_r(nu_k).
# omega_1(z) - z evaluated at nu_k gives lambda_k - nu_k.
# But alpha_k = H_p(lambda_k) - H_r(nu_k), NOT lambda_k - nu_k.
# So this won't directly give <h,alpha>.

# OK, let me try to verify the key identities numerically.
# Focus on: what is the behavior of omega_1 at its poles?

for n in [3]:
    np.random.seed(200)
    p_roots = generate_random_poly(n)
    q_roots = generate_random_poly(n)

    p_coeffs = roots_to_monic_coeffs(p_roots)
    q_coeffs = roots_to_monic_coeffs(q_roots)
    r_coeffs = boxplus_n(p_coeffs, q_coeffs, n)
    r_roots_raw = np.roots(r_coeffs)
    r_roots = np.sort(np.real(r_roots_raw))

    print(f"\n  n={n}")
    print(f"  p_roots = {p_roots}")
    print(f"  r_roots = {r_roots}")

    # Map out omega_1 on the real line to find its poles
    x_vals = np.linspace(min(r_roots) - 5, max(r_roots) + 5, 1000)
    omega1_vals = []
    for x in x_vals:
        z = x + 1j*1e-6
        try:
            w = get_omega1_branch(z, p_roots, r_roots)
            omega1_vals.append(np.real(w))
        except:
            omega1_vals.append(np.nan)
    omega1_vals = np.array(omega1_vals)

    # Find where omega_1 jumps (indicates poles)
    diffs = np.abs(np.diff(omega1_vals))
    pole_indices = np.where(diffs > 1.0)[0]
    if len(pole_indices) > 0:
        print(f"  Approximate poles of omega_1 at x ~= {x_vals[pole_indices]}")

    # Verify omega_1(nu_k) = lambda_k
    for k in range(n):
        w_k = omega1_at_root(r_roots[k], p_roots, r_roots)
        print(f"  omega_1(nu_{k}) = {w_k:.6f}, lambda_{k} = {p_roots[k]:.6f}, "
              f"match = {abs(w_k - p_roots[k]) < 0.01}")

print("\n" + "="*70)
print("APPROACH 6: Partial fraction of omega_1 and sign analysis")
print("="*70)

# omega_1(z) = z + c + sum_j m_j / (z - p_j)  with m_j > 0, p_j real
# omega_1'(z) = 1 - sum_j m_j / (z-p_j)^2
# omega_1''(z) = 2*sum_j m_j / (z-p_j)^3
# alpha_k = omega_1''(nu_k)/2 = sum_j m_j / (nu_k - p_j)^3
#
# h_k = H_r(nu_k) = sum_{l!=k} 1/(nu_k - nu_l)
#
# <h,alpha> = sum_k h_k * alpha_k
# = sum_k [sum_{l!=k} 1/(nu_k-nu_l)] * [sum_j m_j/(nu_k-p_j)^3]
# = sum_j m_j * sum_k [sum_{l!=k} 1/(nu_k-nu_l)] * 1/(nu_k-p_j)^3
# = sum_j m_j * A_j
# where A_j = sum_k H_r(nu_k) / (nu_k - p_j)^3
#
# If we can show A_j >= 0 for all j, then <h,alpha> >= 0 since m_j > 0.
#
# A_j = sum_k h_k / (nu_k - p_j)^3
# h_k = sum_{l!=k} 1/(nu_k - nu_l) has alternating sign: h_1 < 0, h_n > 0 (for sorted roots)
# 1/(nu_k - p_j)^3 also has a sign depending on whether nu_k > p_j or nu_k < p_j.

# Let's verify this decomposition and check if A_j >= 0

for n in [3, 4]:
    np.random.seed(300 + n)

    for trial in range(5):
        p_roots = generate_random_poly(n)
        q_roots = generate_random_poly(n)

        p_coeffs = roots_to_monic_coeffs(p_roots)
        q_coeffs = roots_to_monic_coeffs(q_roots)
        r_coeffs = boxplus_n(p_coeffs, q_coeffs, n)
        r_roots_raw = np.roots(r_coeffs)

        if np.max(np.abs(np.imag(r_roots_raw))) > 1e-6:
            continue
        r_roots = np.sort(np.real(r_roots_raw))
        if np.min(np.diff(r_roots)) < 1e-8:
            continue

        H_r = H_values(r_roots)
        H_p = H_values(p_roots)

        # Compute omega_1 data at each root
        alpha = H_p - H_r
        ha = np.dot(H_r, alpha)

        # Try to find the poles of omega_1 by solving omega_1'(z) = 0
        # For the Herglotz function omega_1(z) = z + c + sum m_j/(z-p_j),
        # omega_1' = 0 gives sum m_j/(z-p_j)^2 = 1
        # The poles p_j of omega_1 interlace with the nu_k (roots of r)
        # Actually, omega_1 has n-1 poles (it's degree n rational function
        # since it comes from solving a degree-n equation)

        # Find poles by looking at omega_1 on real line
        x_test = np.linspace(min(r_roots) - 3, max(r_roots) + 3, 5000)
        om_vals = []
        for x in x_test:
            z = x + 1j*1e-8
            try:
                w = get_omega1_branch(z, p_roots, r_roots)
                om_vals.append(np.real(w))
            except:
                om_vals.append(np.nan)
        om_vals = np.array(om_vals)

        # Detect jumps (poles)
        diffs_abs = np.abs(np.diff(om_vals))
        threshold = np.nanmedian(diffs_abs) * 100
        pole_idx = np.where(diffs_abs > max(threshold, 0.5))[0]

        if len(pole_idx) > 0:
            pole_locs = x_test[pole_idx]
            print(f"\n  n={n}, trial {trial}: <h,alpha> = {ha:.6f}")
            print(f"    Poles of omega_1 near: {pole_locs}")
            print(f"    r_roots = {r_roots}")

            # For each detected pole, compute A_j
            for pj in pole_locs:
                A_j = np.sum(H_r / (r_roots - pj)**3)
                print(f"    A_j at p={pj:.4f}: {A_j:.6f}")

print("\n" + "="*70)
print("FINAL: Testing the decomposition <h,alpha> = sum_j m_j * A_j")
print("="*70)
