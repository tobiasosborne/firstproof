#!/usr/bin/env python3
"""
Clean verification of contour integral identities.

KEY PROVEN IDENTITY: sum_k h_k/(nu_k - p) = -(1/2)*r''(p)/r(p)

We use this to study <h,alpha> = sum_k h_k * alpha_k
via the contour integral of [nG_r(z)]^2 * Psi(z) dz
where Psi(z) = z - omega_1(z).

PROVEN: (1/2pi i) oint [nG_r(z)]^2 * Psi(z) dz = 0
gives: sum_k 2*h_k*(nu_k - lambda_k) = sum_j m_j * [nG_r(p_j)]^2 >= 0

QUESTION: What contour integral gives <h,alpha>?

IDEA: Instead of Psi(z) = z - omega_1(z) (position displacement),
use the "Hilbert transform displacement" which is the regular part
of nG_p(omega_1(z)) minus nG_r(z) at nu_k.

But nG_p(omega_1(z)) = nG_r(z) identically! So this is zero.

The key insight is that alpha_k comes from the SECOND derivative of omega_1,
not from the first-order displacement.

Let me try: what if we integrate [nG_r(z)]^3 * Psi(z)?
"""

import numpy as np
from math import factorial

np.random.seed(42)

def boxplus_n(p_coeffs, q_coeffs, n):
    a, b = p_coeffs, q_coeffs
    r_arr = np.zeros(n + 1)
    r_arr[0] = 1.0
    for k in range(1, n + 1):
        c_k = 0.0
        for i in range(0, k + 1):
            j = k - i
            if i <= n and j <= n:
                coeff = (factorial(n - i) * factorial(n - j)) / (factorial(n) * factorial(n - k))
                c_k += coeff * a[i] * b[j]
        r_arr[k] = c_k
    return r_arr

def H_values(roots):
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H

def generate_random_poly(n, spread=5.0):
    roots = np.sort(np.random.uniform(-spread, spread, n))
    for i in range(1, n):
        if roots[i] - roots[i-1] < 0.5:
            roots[i] = roots[i-1] + 0.5 + np.random.uniform(0, 0.5)
    return roots


# Let me think about what <h,alpha> really is in terms of known quantities.
#
# <h,alpha> = sum_k h_k*alpha_k = sum_k h_k*u_k - sum_k h_k^2
#
# The SECOND sum is Phi(r).
# The FIRST sum is the "cross Fisher information" between r and p.
#
# KEY: Can we express sum_k h_k*u_k as a contour integral?
#
# h_k = H_r(nu_k) = regular part of nG_r at nu_k
# u_k = H_p(lambda_k) = regular part of nG_p at lambda_k
#
# h_k*u_k: this is a product of two "regular parts" at two DIFFERENT points
# (nu_k and lambda_k = omega_1(nu_k)).
#
# APPROACH: Express both h_k and u_k as contour integrals.
# h_k = (1/2pi i) oint_small [nG_r(z) - 1/(z-nu_k)] dz   ... no, this is trivially 0.
# Actually: h_k = (1/2pi i) oint_small [nG_r(z) - 1/(z-nu_k)] * (z-nu_k)^{-1} ... no.
# h_k = regular part at nu_k = coefficient of (z-nu_k)^0 in the Laurent expansion.
# = (1/2pi i) oint nG_r(z) * 1 dz / (some correction)
#
# More useful: Phi_n(r) = sum h_k^2.
# The standard formula: Phi_n(r) = sum_k H_r(nu_k)^2
# This can be written as:
# Phi_n(r) = sum_k Res_{nu_k} [((nG_r)^2 - sum_k 1/(z-nu_k)^2) * nG_r / ???]
# This is getting nowhere useful.

# Let me try a COMPLETELY DIFFERENT approach.
# Instead of contour integrals, try to prove <h,alpha> >= 0 directly from properties
# of the polynomial r and the subordination function omega_1.

# KEY OBSERVATION: <h,alpha> = sum_k H_r(nu_k) * [H_p(lambda_k) - H_r(nu_k)]
# This is the inner product of H_r with the "discrete derivative of H" along omega_1.

# CAUCHY INTERLACING IDEA:
# The roots of r and p satisfy lambda_k = omega_1(nu_k), and omega_1 is a
# monotone-increasing rational function on each interval between its poles.
# The poles of omega_1 interlace with the roots of r (they are the critical
# points of r). So omega_1 maps (nu_k, nu_{k+1}) to an interval that
# contains lambda_k (or lambda_{k+1}).

# Actually, the ORDER-PRESERVING property omega_1(nu_k) = lambda_k with
# lambda_1 < ... < lambda_n and nu_1 < ... < nu_n means omega_1 preserves
# the ordering. Combined with omega_1'(nu_k) = 1, the function omega_1
# has slope 1 at each root of r, and its poles (where it jumps from +inf to -inf)
# are at the critical points of r, which lie between consecutive roots.

# The SIGN of alpha_k = H_p(lambda_k) - H_r(nu_k) depends on whether the
# "spreading of roots" from p to r increases or decreases H at that index.

# For the EXTREMAL roots (k=1 and k=n):
# H_r(nu_1) < 0 (most negative), H_r(nu_n) > 0 (most positive)
# If the roots of r are more spread than those of p (which they are, since
# r = p boxplus q), then lambda_1 > nu_1 and lambda_n < nu_n.
# This means H_p(lambda_1) is less negative than H_r(nu_1) in some cases...
# Actually this isn't clear without more detailed analysis.

# Let me instead verify whether <h,alpha> has a nice expression in terms of
# the polynomials themselves.

print("=" * 70)
print("SEARCH FOR POLYNOMIAL IDENTITY")
print("=" * 70)

# <h,alpha> = sum_k H_r(nu_k) * H_p(lambda_k) - Phi_n(r)
#
# Phi_n(r) = sum_k [r''(nu_k)/(2*r'(nu_k))]^2
# = (1/4) * sum_k [r''(nu_k)]^2 / [r'(nu_k)]^2
#
# sum_k H_r(nu_k)*H_p(lambda_k) = sum_k [r''(nu_k)/(2*r'(nu_k))] * [p''(lambda_k)/(2*p'(lambda_k))]
# = (1/4) * sum_k r''(nu_k)*p''(lambda_k) / (r'(nu_k)*p'(lambda_k))

# Using omega_1(nu_k) = lambda_k and omega_1'(nu_k) = 1:
# The chain rule nG_r = nG_p(omega_1) differentiated gives:
# [nG_r]' = [nG_p]'(omega_1)*omega_1'
# At z = nu_k (AWAY from the pole, i.e., considering the REGULAR part):
# We can't just evaluate at nu_k because these have poles.

# Actually from the subordination identity r'(z)/r(z) = p'(omega_1(z))/p(omega_1(z)):
# Differentiate: [r''*r - (r')^2]/r^2 = [p''(omega)*omega'*p(omega) - p'(omega)^2*omega']/p(omega)^2
# At z = nu_k: r(nu_k) = 0, so the LHS blows up. Not directly useful.

# Let me try: what is the polynomial that gives <h,u> = sum H_r(nu_k)*H_p(lambda_k)?
# Using the fact that p has roots lambda_k:
# p'(lambda_k) = prod_{j!=k} (lambda_k - lambda_j)
# p''(lambda_k) = 2 * sum_{j!=k} prod_{l!=k,l!=j} (lambda_k - lambda_l)
# So H_p(lambda_k) = p''(lambda_k)/(2*p'(lambda_k))

# Consider the polynomial Q(z) = sum_k [r''(nu_k)*p''(lambda_k)] / (r'(nu_k)*p'(lambda_k)) * L_k(z)
# where L_k is the Lagrange basis polynomial for the points nu_1,...,nu_n.
# Q(nu_k) = r''(nu_k)*p''(lambda_k) / (r'(nu_k)*p'(lambda_k)) = 4*h_k*u_k

# Hmm, this is circular. Let me focus on numerical verification of various candidate
# expressions for <h,alpha>.

for n in [3, 4, 5]:
    np.random.seed(42 + n)

    for trial in range(3):
        p_roots = generate_random_poly(n)
        q_roots = generate_random_poly(n)
        p_coeffs = np.poly(p_roots)
        q_coeffs = np.poly(q_roots)
        r_coeffs = boxplus_n(p_coeffs, q_coeffs, n)
        r_roots_raw = np.roots(r_coeffs)
        if np.max(np.abs(np.imag(r_roots_raw))) > 1e-6:
            continue
        r_roots = np.sort(np.real(r_roots_raw))
        if np.min(np.diff(r_roots)) < 1e-8:
            continue

        H_r = H_values(r_roots)
        H_p = H_values(p_roots)

        ha = np.dot(H_r, H_p - H_r)
        hu = np.dot(H_r, H_p)
        phi_r = np.dot(H_r, H_r)
        phi_p = np.dot(H_p, H_p)

        # Test candidate expressions:
        r_poly = np.poly(r_roots)
        p_poly = np.poly(p_roots)
        r_d1 = np.polyder(r_poly)
        r_d2 = np.polyder(r_d1)
        p_d1 = np.polyder(p_poly)
        p_d2 = np.polyder(p_d1)

        # Candidate 1: sum r''(nu_k)*p''(lambda_k)/(r'(nu_k)*p'(lambda_k))
        cand1 = 0.0
        for k in range(n):
            cand1 += (np.polyval(r_d2, r_roots[k]) * np.polyval(p_d2, p_roots[k]) /
                     (np.polyval(r_d1, r_roots[k]) * np.polyval(p_d1, p_roots[k])))
        cand1 /= 4  # since H = r''/(2r'), product has factor 4

        # Candidate 2: (1/2)*sum [r''(nu_k)/r'(nu_k)] * [p''(lambda_k)/p'(lambda_k)]
        cand2 = 0.0
        for k in range(n):
            cand2 += ((np.polyval(r_d2, r_roots[k]) / np.polyval(r_d1, r_roots[k])) *
                     (np.polyval(p_d2, p_roots[k]) / np.polyval(p_d1, p_roots[k])))
        cand2 /= 4

        # Candidate 3: - (1/2) sum_k r''(nu_k) * p(nu_k) / (r'(nu_k) * p'(nu_k)) ... nah

        # Candidate 4: Residue of r''/(2r) * p''/(2p) at... doesn't make sense

        # Candidate 5: sum_k 1/r'(nu_k) * p''(omega(nu_k))/(2*p'(omega(nu_k)))
        # = sum_k [1/r'(nu_k)] * H_p(lambda_k)
        # But sum 1/r'(nu_k) * g(nu_k) = sum Res_{nu_k} [g(z)/r(z)] * ... hmm

        # A clean identity: for any function g regular at all nu_k,
        # sum_k g(nu_k)/r'(nu_k) = (1/2pi i) oint g(z)/r(z) dz
        # This is just the Lagrange interpolation / residue sum.

        # So sum_k H_p(lambda_k)/r'(nu_k) = (1/2pi i) oint H_p(omega_1(z))/r(z) dz?
        # No, H_p(omega_1(z)) is not a simple function.

        # But: sum_k H_r(nu_k) * H_p(lambda_k) is NOT of the form sum g(nu_k)/r'(nu_k).
        # It's sum_k g(nu_k)*h(nu_k) where g(nu_k) = H_r(nu_k), h(nu_k) = H_p(lambda_k).
        # This is a bilinear sum.

        print(f"\nn={n}, trial {trial}:")
        print(f"  <h,u> = {hu:.8f}")
        print(f"  cand1 (from r'',p'') = {cand1:.8f}")
        print(f"  cand2 (same) = {cand2:.8f}")
        print(f"  ||h||^2 = {phi_r:.8f}")
        print(f"  <h,alpha> = {ha:.8f}")

        # KEY: Express <h,u> via r'' and the Cauchy transform.
        # H_r(nu_k) = nG_r(z)|_{reg at nu_k} = sum_{j!=k} 1/(nu_k - nu_j)
        # This equals (d/dz log r(z))|_{z=nu_k} minus the singular part 1/(z-nu_k).
        # So H_r(nu_k) = [d/dz log r(z) - 1/(z-nu_k)]|_{z=nu_k}
        # = [r'(z)/r(z) - 1/(z-nu_k)]|_{z=nu_k}
        # = lim_{z->nu_k} [r'(z)*(z-nu_k) - r(z)] / [r(z)*(z-nu_k)]
        # Using L'Hopital or Taylor: r(z) = r'(nu_k)*(z-nu_k) + (1/2)*r''(nu_k)*(z-nu_k)^2 + ...
        # r'(z) = r'(nu_k) + r''(nu_k)*(z-nu_k) + ...
        # r'(z)*(z-nu_k) - r(z) = [r'(nu_k)+r''(nu_k)*eps+...]*eps - [r'(nu_k)*eps+(1/2)*r''(nu_k)*eps^2+...]
        # = (1/2)*r''(nu_k)*eps^2 + ...
        # Denominator: r(z)*(z-nu_k) = r'(nu_k)*eps^2 + ...
        # So H_r(nu_k) = (1/2)*r''(nu_k)/r'(nu_k) = r''(nu_k)/(2*r'(nu_k))

        # Now: <h,u> = sum_k H_r(nu_k) * H_p(lambda_k)
        # = sum_k [r''(nu_k)/(2*r'(nu_k))] * [p''(lambda_k)/(2*p'(lambda_k))]
        # = (1/4) * sum_k [r''(nu_k)/r'(nu_k)] * [p''(lambda_k)/p'(lambda_k)]

        # And Phi_r = sum_k [r''(nu_k)/(2*r'(nu_k))]^2
        # = (1/4) * sum_k [r''(nu_k)/r'(nu_k)]^2

        # So <h,alpha> = (1/4) * sum_k [r''(nu_k)/r'(nu_k)] *
        #                       {[p''(lambda_k)/p'(lambda_k)] - [r''(nu_k)/r'(nu_k)]}

        # Define a_k = r''(nu_k)/r'(nu_k) = 2*H_r(nu_k) = 2*h_k
        # Define b_k = p''(lambda_k)/p'(lambda_k) = 2*H_p(lambda_k) = 2*u_k
        # Then <h,alpha> = (1/4) * sum_k a_k * (b_k - a_k)

        # We want to show sum a_k*(b_k - a_k) >= 0, i.e., sum a_k*b_k >= sum a_k^2.

        # This is <a,b> >= ||a||^2 where a_k = 2*h_k, b_k = 2*u_k.
        # Which is just <h,u> >= ||h||^2 restated.

        # So the question reduces to: for the two real-rooted polynomials r and p
        # related by omega_1 with omega_1'(nu_k) = 1, show:
        # sum_k [r''(nu_k)/r'(nu_k)] * [p''(lambda_k)/p'(lambda_k)] >=
        # sum_k [r''(nu_k)/r'(nu_k)]^2

        # Note: this has a Cauchy-Schwarz "flavor" if a and b were related by a contraction.
        # If |b_k| <= |a_k| for all k, then <a,b> <= ||a||^2 (WRONG direction).
        # If b_k = a_k + c_k with <a,c> >= 0, that's what we want.

        # WAIT: there's a simpler form. <h,u> >= ||h||^2 is equivalent to
        # <h, u-h> >= 0, i.e., <h, alpha> >= 0. So we're going in circles.

        # Let me try to find an energy functional or generating function.

print("\n" + "=" * 70)
print("APPROACH: Energy functional E(t) = sum [nG_{r_t}(nu_k(t))]_{reg}^2")
print("where r_t interpolates between p and r")
print("=" * 70)

# Consider the parametric family: omega_t(z) such that
# omega_0(z) = z (so r_0 = r) and omega_1(z) = omega_1(z) (so the roots are lambda_k).
# Define omega_t(z) = (1-t)*z + t*omega_1(z) for t in [0,1].
# Then omega_t(nu_k) = (1-t)*nu_k + t*lambda_k
# omega_t'(nu_k) = (1-t) + t*1 = 1 for all t  (since omega_1'(nu_k) = 1).

# The roots "at time t" are: xi_k(t) = omega_t^{-1}(mu_k(t)) ... hmm, this isn't right.
# Let me think differently.

# Actually, consider the family of root vectors:
# x_k(t) = (1-t)*nu_k + t*lambda_k
# H(x(t))_k = sum_{j!=k} 1/(x_k(t) - x_j(t))

# Then <H(x(t)), x'(t)> = sum_k H(x(t))_k * (lambda_k - nu_k)
# and d/dt [sum H(x(t))_k^2] = 2*sum H(x(t))_k * d/dt[H(x(t))_k]

# This is getting complicated. Let me try the simplest possible thing:
# just verify that <h,alpha> can be written as a sum of explicitly positive terms.

# IDEA: <h,alpha> = sum_k h_k*(u_k - h_k)
# = sum_k h_k*u_k - sum_k h_k^2
# = sum_{k,j:j!=k} [1/(nu_k-nu_j)] * [1/(lambda_k-lambda_j)] * delta_sign
#   - sum_{k,j:j!=k} [1/(nu_k-nu_j)] * [1/(nu_k-nu_j)]
# where the first sum is the "cross" term.

# Let me expand:
# h_k*u_k = [sum_{j!=k} 1/(nu_k-nu_j)] * [sum_{l!=k} 1/(lambda_k-lambda_l)]
# = sum_{j!=k, l!=k} 1/((nu_k-nu_j)(lambda_k-lambda_l))

# h_k^2 = [sum_{j!=k} 1/(nu_k-nu_j)]^2
# = sum_{j!=k, l!=k} 1/((nu_k-nu_j)(nu_k-nu_l))

# <h,alpha> = sum_k sum_{j!=k, l!=k} [1/((nu_k-nu_j)(lambda_k-lambda_l)) - 1/((nu_k-nu_j)(nu_k-nu_l))]
# = sum_k sum_{j!=k, l!=k} 1/(nu_k-nu_j) * [1/(lambda_k-lambda_l) - 1/(nu_k-nu_l)]
# = sum_k sum_{j!=k, l!=k} 1/(nu_k-nu_j) * [(nu_k-nu_l) - (lambda_k-lambda_l)] / [(lambda_k-lambda_l)(nu_k-nu_l)]
# = sum_k sum_{j!=k, l!=k} [(delta_k - delta_l)] / [(nu_k-nu_j)(lambda_k-lambda_l)(nu_k-nu_l)]
# where delta_i = nu_i - lambda_i

# This triple sum decomposition might help if we can show the summand has a definite sign.
# But it's a triple sum over (k,j,l) with j!=k, l!=k, and the sign of each term
# depends on the relative positions.

# Let's try something even more explicit: use the PROVEN identity
# sum_k h_k/(nu_k-p) = -(1/2)*r''(p)/r(p) = -(1/2)*F(p)
# and the analogous identity for p:
# sum_k u_k/(lambda_k-q) = -(1/2)*p''(q)/p(q)

# <h,u> = sum_k h_k*u_k
# = sum_k h_k * sum_{l!=k} 1/(lambda_k - lambda_l)
# = sum_{k,l: l!=k} h_k / (lambda_k - lambda_l)

# For fixed l: sum_{k!=l} h_k / (lambda_k - lambda_l)
# This is NOT the same as sum_k h_k/(nu_k - p) evaluated at p = lambda_l
# because the sum excludes k=l.

# sum_{k!=l} h_k/(lambda_k - lambda_l)
# = [sum_k h_k/(lambda_k - lambda_l)] - h_l/(lambda_l - lambda_l)
# The second term is undefined (0/0)! So this decomposition doesn't work directly.

# Let me try the substitution lambda_k = omega_1(nu_k):
# <h,u> = sum_{k,l: l!=k} h_k / (omega_1(nu_k) - omega_1(nu_l))
#
# Now: omega_1(nu_k) - omega_1(nu_l) = (nu_k - nu_l) * omega_1[nu_k,nu_l]
# where omega_1[nu_k,nu_l] is the divided difference.
# Since omega_1'(nu_k) = 1 and omega_1 is monotone, omega_1[nu_k,nu_l] > 0.
#
# <h,u> = sum_{k,l:l!=k} h_k / ((nu_k-nu_l) * omega_1[nu_k,nu_l])
#
# Compare with ||h||^2 = sum_{k,l:l!=k} h_k / (nu_k - nu_l) * sum_{m!=k} 1/(nu_k-nu_m) ... no
# ||h||^2 = sum_k h_k^2 = sum_k h_k * sum_{l!=k} 1/(nu_k-nu_l)
# = sum_{k,l:l!=k} h_k / (nu_k - nu_l)
#
# So <h,u> - ||h||^2 = sum_{k,l:l!=k} h_k * [1/((nu_k-nu_l)*omega_1[nu_k,nu_l]) - 1/(nu_k-nu_l)]
# = sum_{k,l:l!=k} h_k / (nu_k-nu_l) * [1/omega_1[nu_k,nu_l] - 1]
# = sum_{k,l:l!=k} h_k / (nu_k-nu_l) * (1 - omega_1[nu_k,nu_l]) / omega_1[nu_k,nu_l]

print("\nTesting divided difference decomposition:")
for n in [3, 4]:
    np.random.seed(42 + n)
    for trial in range(3):
        p_roots = generate_random_poly(n)
        q_roots = generate_random_poly(n)
        p_coeffs = np.poly(p_roots)
        q_coeffs = np.poly(q_roots)
        r_coeffs = boxplus_n(p_coeffs, q_coeffs, n)
        r_roots_raw = np.roots(r_coeffs)
        if np.max(np.abs(np.imag(r_roots_raw))) > 1e-6:
            continue
        r_roots = np.sort(np.real(r_roots_raw))
        if np.min(np.diff(r_roots)) < 1e-8:
            continue

        H_r = H_values(r_roots)
        H_p = H_values(p_roots)
        ha = np.dot(H_r, H_p - H_r)

        # Compute omega_1 divided differences
        # omega_1[k,l] = (lambda_k - lambda_l) / (nu_k - nu_l)
        ha_decomp = 0.0
        for k in range(n):
            for l in range(n):
                if l == k:
                    continue
                omega_dd = (p_roots[k] - p_roots[l]) / (r_roots[k] - r_roots[l])
                ha_decomp += H_r[k] / (r_roots[k] - r_roots[l]) * (1 - omega_dd) / omega_dd

        print(f"\n  n={n}, trial {trial}:")
        print(f"    <h,alpha> direct = {ha:.8f}")
        print(f"    <h,alpha> decomp = {ha_decomp:.8f}")
        print(f"    Match: {abs(ha - ha_decomp) < 1e-6}")

        # Check sign of individual terms
        pos_count = 0
        neg_count = 0
        for k in range(n):
            for l in range(n):
                if l == k:
                    continue
                omega_dd = (p_roots[k] - p_roots[l]) / (r_roots[k] - r_roots[l])
                term = H_r[k] / (r_roots[k] - r_roots[l]) * (1 - omega_dd) / omega_dd
                if term > 0:
                    pos_count += 1
                else:
                    neg_count += 1

        print(f"    Terms: {pos_count} positive, {neg_count} negative")

        # Check divided differences
        print(f"    omega_1 divided differences:")
        for k in range(n):
            for l in range(k+1, n):
                omega_dd = (p_roots[k] - p_roots[l]) / (r_roots[k] - r_roots[l])
                print(f"      omega_1[{k},{l}] = {omega_dd:.6f}")


print("\n" + "=" * 70)
print("KEY IDENTITY: <h,alpha> = sum_{k!=l} h_k/(nu_k-nu_l) * (1/D_{kl} - 1)")
print("where D_{kl} = omega_1[nu_k,nu_l] = (lambda_k-lambda_l)/(nu_k-nu_l)")
print("=" * 70)

# Note: h_k/(nu_k-nu_l) = 1/((nu_k-nu_l)^2) * sum_{m!=k} 1/... no.
# h_k = sum_{m!=k} 1/(nu_k-nu_m), so h_k/(nu_k-nu_l) is a specific term
# in the partial fraction expansion.

# The divided difference D_{kl} = (lambda_k - lambda_l)/(nu_k - nu_l) has:
# - D_{kl} > 0 (since lambda and nu are both sorted, order-preserving)
# - D_{kk} = omega_1'(nu_k) = 1 (by the chain rule identity)
# - Interlacing: each D_{kl} > 0

# Now: 1/D_{kl} - 1 = (1 - D_{kl})/D_{kl}
# If D_{kl} < 1: 1/D_{kl} > 1, so 1/D_{kl} - 1 > 0
# If D_{kl} > 1: 1/D_{kl} < 1, so 1/D_{kl} - 1 < 0

# So the sign of each term depends on h_k/(nu_k-nu_l) and 1/D_{kl} - 1.

# MEAN VALUE THEOREM: For omega_1 with omega_1'(nu_k) = omega_1'(nu_l) = 1,
# by the mean value theorem, D_{kl} = omega_1'(xi) for some xi between nu_k and nu_l.
# Since omega_1'(z) = 1 + sum m_j/(z-p_j)^2 >= 1 + ... wait:
# omega_1(z) = z + b + sum m_j/(p_j-z), omega_1'(z) = 1 + sum m_j/(p_j-z)^2
# Since m_j > 0 and (p_j-z)^2 > 0 for real z != p_j:
# omega_1'(z) = 1 + sum m_j/(p_j-z)^2 >= 1

# BEAUTIFUL! omega_1'(z) >= 1 for all real z (not at poles).
# Therefore by MVT: D_{kl} = omega_1'(xi) >= 1 for all k != l.
# So D_{kl} >= 1, which means 1/D_{kl} <= 1, which means 1/D_{kl} - 1 <= 0.

# So each factor (1/D_{kl} - 1) <= 0.
# And h_k/(nu_k - nu_l): when k < l, nu_k < nu_l, so nu_k - nu_l < 0,
# and h_k has sign... varying.

# Hmm, the individual terms have indefinite sign. Let me check.

print("\nVerification that omega_1'(z) >= 1 everywhere (since omega_1' = 1 + positive terms):")
print("This follows from omega_1(z) = z + c + sum m_j/(p_j - z) with m_j > 0")
print("omega_1'(z) = 1 + sum m_j/(p_j - z)^2 >= 1 for real z != p_j")
print()
print("By MVT: D_{kl} = (lambda_k - lambda_l)/(nu_k - nu_l) >= 1")
print("This means 1/D_{kl} - 1 <= 0 for all k != l")
print()

# So <h,alpha> = sum_{k!=l} h_k/(nu_k-nu_l) * (1/D_{kl} - 1)
# where 1/D_{kl} - 1 <= 0 for all k,l.
# The sign of h_k/(nu_k-nu_l) depends on the indices.

# Let's check whether we can symmetrize.
# <h,alpha> = sum_{k!=l} h_k/(nu_k-nu_l) * (1/D_{kl} - 1)
# Swap k,l:
# = sum_{k!=l} h_l/(nu_l-nu_k) * (1/D_{lk} - 1)
# Note D_{lk} = (lambda_l - lambda_k)/(nu_l - nu_k) = D_{kl} (the same!)
# And h_l/(nu_l-nu_k) = -h_l/(nu_k-nu_l)
# So the "swapped" sum = -sum_{k!=l} h_l/(nu_k-nu_l) * (1/D_{kl} - 1)

# Average: <h,alpha> = (1/2)*sum_{k!=l} (h_k - h_l)/(nu_k-nu_l) * (1/D_{kl} - 1)

# Now: (h_k - h_l)/(nu_k - nu_l) is the divided difference of H_r!
# And (1/D_{kl} - 1) <= 0.

# So <h,alpha> = (1/2)*sum_{k!=l} H_r[nu_k,nu_l] * (1/D_{kl} - 1)
# where H_r[nu_k,nu_l] = (H_r(nu_k) - H_r(nu_l))/(nu_k - nu_l)

# For <h,alpha> >= 0 with all (1/D_{kl} - 1) <= 0, we need:
# sum_{k!=l} H_r[nu_k,nu_l] * (1/D_{kl} - 1) >= 0
# Since 1/D - 1 <= 0, this requires that the "weighted average" of H_r[k,l]
# is <= 0.

# But H_r[nu_k,nu_l] = (h_k - h_l)/(nu_k - nu_l).
# For k < l: nu_k < nu_l. What's the sign of h_k - h_l?
# Since h_k = H_r(nu_k) goes from negative (small k) to positive (large k),
# h_k - h_l < 0 for k < l.
# And nu_k - nu_l < 0 for k < l.
# So H_r[k,l] = (h_k - h_l)/(nu_k - nu_l) = negative/negative > 0 ?
# Hmm, not necessarily. H_r is not monotone increasing in general for n >= 4.

print("\nChecking sign of H_r divided differences:")
for n in [3, 4, 5]:
    np.random.seed(42 + n)
    count_neg = 0
    count_total = 0
    for trial in range(100):
        r_roots = generate_random_poly(n)
        H_r = H_values(r_roots)
        for k in range(n):
            for l in range(k+1, n):
                dd = (H_r[k] - H_r[l]) / (r_roots[k] - r_roots[l])
                count_total += 1
                if dd < -1e-10:
                    count_neg += 1
    print(f"  n={n}: H_r[k,l] < 0 in {count_neg}/{count_total}")

# KEY IDENTITY (proven):
# <h,alpha> = (1/2)*sum_{k<l} [(h_k-h_l)/(nu_k-nu_l)] * 2*(1/D_{kl}-1)
# Wait, I need to be careful about the factor of 2 from the symmetrization.
# sum_{k!=l} = 2*sum_{k<l}
# <h,alpha> = (1/2)*2*sum_{k<l} [(h_k-h_l)/(nu_k-nu_l)] * (1/D_{kl}-1)
# = sum_{k<l} [(h_k-h_l)/(nu_k-nu_l)] * (1/D_{kl}-1)

# Verify this:
print("\nVerifying: <h,alpha> = sum_{k<l} H_r[k,l] * (1/D_{kl} - 1)")
for n in [3, 4]:
    np.random.seed(42 + n)
    for trial in range(3):
        p_roots = generate_random_poly(n)
        q_roots = generate_random_poly(n)
        p_coeffs = np.poly(p_roots)
        q_coeffs = np.poly(q_roots)
        r_coeffs = boxplus_n(p_coeffs, q_coeffs, n)
        r_roots_raw = np.roots(r_coeffs)
        if np.max(np.abs(np.imag(r_roots_raw))) > 1e-6:
            continue
        r_roots = np.sort(np.real(r_roots_raw))
        if np.min(np.diff(r_roots)) < 1e-8:
            continue

        H_r = H_values(r_roots)
        H_p = H_values(p_roots)
        ha = np.dot(H_r, H_p - H_r)

        ha_sym = 0.0
        for k in range(n):
            for l in range(k+1, n):
                D_kl = (p_roots[k] - p_roots[l]) / (r_roots[k] - r_roots[l])
                Hr_dd = (H_r[k] - H_r[l]) / (r_roots[k] - r_roots[l])
                ha_sym += Hr_dd * (1.0/D_kl - 1.0)

        print(f"  n={n}: <h,alpha> = {ha:.8f}, symmetric = {ha_sym:.8f}, "
              f"match = {abs(ha - ha_sym) < 1e-6}")

print("\n" + "=" * 70)
print("FINAL KEY IDENTITY:")
print("<h,alpha> = -sum_{k<l} [(h_k-h_l)/(nu_k-nu_l)] * [1 - 1/D_{kl}]")
print("where D_{kl} >= 1 (proven) and (h_k-h_l)/(nu_k-nu_l) has definite sign?")
print("=" * 70)
