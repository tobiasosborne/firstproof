#!/usr/bin/env python3
"""
CAREFUL sign verification for the Herglotz decomposition.

The Herglotz representation for a function f: C+ -> C+ is:
f(z) = a*z + b + integral dmu(t)/(t-z)  with a >= 0, mu positive measure

For rational functions with simple real poles:
f(z) = a*z + b + sum_j m_j/(p_j - z)  with a >= 0, m_j > 0

Note: 1/(p_j - z) maps C+ to C+ (Im(1/(p_j-z)) > 0 for Im(z) > 0 when p_j real).

For omega_1: maps C+ to C+, omega_1(z) ~ z at infinity
omega_1(z) = z + b + sum_j m_j/(p_j - z)   with m_j > 0

Derivatives:
omega_1'(z) = 1 + sum_j m_j/(p_j - z)^2   = 1 + sum_j m_j/(z - p_j)^2  (same!)
omega_1''(z) = -2*sum_j m_j/(p_j - z)^3   = 2*sum_j m_j/(z - p_j)^3 ... WAIT

d/dz [m/(p-z)] = m/(p-z)^2
d^2/dz^2 [m/(p-z)] = -2m/(p-z)^3 = 2m/(z-p)^3

So omega_1''(z) = -2*sum_j m_j/(p_j - z)^3 = 2*sum_j m_j/(z-p_j)^3

Hmm wait: 1/(p-z)^3 and 1/(z-p)^3 = -1/(p-z)^3, so:
-2*m/(p-z)^3 = -2*m*(-1)/(z-p)^3 = 2*m/(z-p)^3

OK so omega_1''(z) = 2*sum_j m_j/(z-p_j)^3 regardless of notation.
"""

import numpy as np
from math import factorial

np.random.seed(42)

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

def generate_random_poly(n, spread=5.0):
    roots = np.sort(np.random.uniform(-spread, spread, n))
    for i in range(1, n):
        if roots[i] - roots[i-1] < 0.5:
            roots[i] = roots[i-1] + 0.5 + np.random.uniform(0, 0.5)
    return roots


# Step 1: Verify the Herglotz representation and pole locations for omega_1
# For n=2, we can compute omega_1 explicitly.

print("=" * 70)
print("EXPLICIT omega_1 FOR n=2")
print("=" * 70)

# For n=2: p(x) = (x-a)(x-b), r(x) = (x-c)(x-d) with a < b, c < d
# nG_r(z) = 1/(z-c) + 1/(z-d)
# nG_p(w) = 1/(w-a) + 1/(w-b)
# Equation: 1/(z-c) + 1/(z-d) = 1/(w-a) + 1/(w-b)
# i.e., (2z-c-d)/((z-c)(z-d)) = (2w-a-b)/((w-a)(w-b))
# Let s_r = c+d, s_p = a+b, g_r = cd, g_p = ab (where p(x) = x^2 - s_p*x + g_p, etc.)
# Then: (2z-s_r)/(z^2-s_r*z+g_r) = (2w-s_p)/(w^2-s_p*w+g_p)
# Cross-multiply: (2z-s_r)(w^2-s_p*w+g_p) = (2w-s_p)(z^2-s_r*z+g_r)
# This is a quadratic in w. Solve for w.
# (2z-s_r)*w^2 - (2z-s_r)*s_p*w + (2z-s_r)*g_p = (z^2-s_r*z+g_r)*2*w - (z^2-s_r*z+g_r)*s_p
# (2z-s_r)*w^2 - [(2z-s_r)*s_p + 2*(z^2-s_r*z+g_r)]*w + [(2z-s_r)*g_p + (z^2-s_r*z+g_r)*s_p] = 0
# Coefficient of w^2: 2z - s_r
# Coefficient of w: -(2z-s_r)*s_p - 2z^2 + 2*s_r*z - 2*g_r
#                  = -2z*s_p + s_r*s_p - 2z^2 + 2*s_r*z - 2*g_r
#                  = -2z^2 + (2*s_r - 2*s_p)*z + (s_r*s_p - 2*g_r)
# Actually this is getting messy. Let me just compute numerically for a specific example.

a, b = 0.0, 2.0   # p roots
n = 2
q_roots_test = np.array([0.0, 1.0])
p_roots_test = np.array([a, b])
p_coeffs = np.poly(p_roots_test)
q_coeffs = np.poly(q_roots_test)
r_coeffs = boxplus_n(p_coeffs, q_coeffs, n)
r_roots_raw = np.roots(r_coeffs)
r_roots_test = np.sort(np.real(r_roots_raw))

print(f"p = (x-0)(x-2), q = x(x-1)")
print(f"r_roots = {r_roots_test}")
print(f"r has one critical point at: {np.mean(r_roots_test):.6f}")

# omega_1 has n-1 = 1 pole, which should be at the critical point of r.
crit = np.mean(r_roots_test)
r_poly = np.poly(r_roots_test)
r_d2 = np.polyder(np.polyder(r_poly))
m1 = -n * np.polyval(r_poly, crit) / np.polyval(r_d2, crit)
print(f"m_1 = -n*r(c)/r''(c) = {m1:.6f}")
print(f"omega_1(z) = z + const + m_1/(c - z) where m_1 = {m1:.6f}, c = {crit:.6f}")

# Verify: omega_1(nu_k) = lambda_k
for k in range(n):
    nu_k = r_roots_test[k]
    # omega_1(nu_k) = nu_k + const + m_1/(crit - nu_k)
    # We need to find const so that omega_1(nu_1) = a = 0 and omega_1(nu_2) = b = 2
    # Two equations:
    # nu_1 + const + m1/(crit - nu_1) = a
    # nu_2 + const + m1/(crit - nu_2) = b

nu_1, nu_2 = r_roots_test
const_from_1 = a - nu_1 - m1/(crit - nu_1)
const_from_2 = b - nu_2 - m1/(crit - nu_2)
print(f"\nConst from nu_1: {const_from_1:.6f}")
print(f"Const from nu_2: {const_from_2:.6f}")
print(f"Match: {abs(const_from_1 - const_from_2) < 1e-6}")

const = const_from_1
print(f"\nomega_1(z) = z + {const:.6f} + {m1:.6f}/({crit:.6f} - z)")

# Now compute omega_1'' at nu_k:
# omega_1(z) = z + const + m1/(crit - z)
# omega_1'(z) = 1 + m1/(crit-z)^2
# omega_1''(z) = -2*m1/(crit-z)^3 = 2*m1/(z-crit)^3

for k in range(n):
    nu_k = r_roots_test[k]
    om_pp = 2*m1/(nu_k - crit)**3
    alpha_k_computed = om_pp / 2
    alpha_k_from_H = H_values(p_roots_test)[k] - H_values(r_roots_test)[k]
    print(f"\nk={k}: nu_k = {nu_k:.6f}")
    print(f"  omega_1''(nu_k) = {om_pp:.6f}")
    print(f"  alpha_k from Herglotz = {alpha_k_computed:.6f}")
    print(f"  alpha_k from H_p - H_r = {alpha_k_from_H:.6f}")
    print(f"  Match: {abs(alpha_k_computed - alpha_k_from_H) < 0.01}")

# Compute <h,alpha> both ways
H_r = H_values(r_roots_test)
H_p = H_values(p_roots_test)
alpha_true = H_p - H_r
h_alpha_true = np.dot(H_r, alpha_true)

# From Herglotz decomposition:
# <h,alpha> = sum_k h_k * alpha_k where alpha_k = m1/(nu_k - crit)^3
# (using omega_1(z) = z + const + m1/(crit - z), so omega_1'' = 2*m1/(z-crit)^3)
# Wait: omega_1(z) = z + const + m1/(crit-z)
# d/dz[m1/(crit-z)] = m1/(crit-z)^2
# d^2/dz^2[m1/(crit-z)] = -2*m1/(crit-z)^3

# Let me be EXTREMELY careful.
# f(z) = m/(c-z) where c = crit, m = m1
# f'(z) = m/(c-z)^2   (since d/dz[1/(c-z)] = 1/(c-z)^2)
# f''(z) = -2m/(c-z)^3 (since d/dz[1/(c-z)^2] = 2/(c-z)^3... wait
# d/dz[(c-z)^{-2}] = -2*(c-z)^{-3}*(-1) = 2/(c-z)^3
# So f''(z) = 2m/(c-z)^3

# omega_1''(z) = 0 + 0 + f''(z) = 2*m1/(crit-z)^3
# alpha_k = omega_1''(nu_k)/2 = m1/(crit - nu_k)^3

print(f"\n--- Sign check ---")
for k in range(n):
    nu_k = r_roots_test[k]
    alpha_from_herglotz = m1 / (crit - nu_k)**3
    alpha_from_H = H_p[k] - H_r[k]
    print(f"k={k}: alpha_herglotz = {alpha_from_herglotz:.6f}, alpha_H = {alpha_from_H:.6f}")

print(f"\nm1 = {m1:.6f}")
print(f"crit = {crit:.6f}")

# Now <h,alpha> = sum_k h_k * m1/(crit-nu_k)^3 = m1 * sum_k h_k/(crit-nu_k)^3
# = m1 * (-1)^3 * sum_k h_k/(nu_k-crit)^3  ... no
# (crit-nu_k)^3 and (nu_k-crit)^3 = -(crit-nu_k)^3
# So 1/(crit-nu_k)^3 = -1/(nu_k-crit)^3

# <h,alpha> = m1 * sum_k h_k/(crit-nu_k)^3 = -m1 * sum_k h_k/(nu_k-crit)^3

sum_term = np.sum(H_r / (crit - r_roots_test)**3)
h_alpha_herglotz = m1 * sum_term
print(f"\n<h,alpha> from direct: {h_alpha_true:.6f}")
print(f"<h,alpha> from Herglotz: {h_alpha_herglotz:.6f}")
print(f"Match: {abs(h_alpha_true - h_alpha_herglotz) < 0.01}")

# Now relate to F'' at the critical point:
# F(p) = r''(p)/r(p), F'' = ?
# We showed: sum_k h_k/(nu_k-p) = -(1/2)*r''(p)/r(p)
# Note: 1/(nu_k-p) = -1/(p-nu_k)
# And r''(p)/r(p) = F(p)
# So: -sum_k h_k/(p-nu_k) = -(1/2)*F(p), i.e., sum_k h_k/(p-nu_k) = (1/2)*F(p)

# Differentiate with respect to p:
# -sum_k h_k/(p-nu_k)^2 = (1/2)*F'(p)
# d/dp again: 2*sum_k h_k/(p-nu_k)^3 = (1/2)*F''(p)
# So sum_k h_k/(p-nu_k)^3 = (1/4)*F''(p)

# Now 1/(crit-nu_k)^3 = -1/(nu_k-crit)^3 = 1/(crit-nu_k)^3
# And (crit-nu_k) = -(nu_k-crit), so (crit-nu_k)^3 = -(nu_k-crit)^3
# Hmm: 1/(crit-nu_k)^3 = 1/(-(nu_k-crit))^3 = -1/(nu_k-crit)^3

# sum_k h_k/(crit-nu_k)^3 = -sum_k h_k/(nu_k-crit)^3

# Using our identity with p = crit:
# sum_k h_k/(crit-nu_k)^3 = (1/4)*F''(crit)   (? need to check)

# Our identity: sum_k h_k/(p-nu_k)^3 = (1/4)*F''(p)
# At p = crit: sum_k h_k/(crit-nu_k)^3 = (1/4)*F''(crit)

h_term = np.sum(H_r / (crit - r_roots_test)**3)
# F''(crit):
h_step = 1e-6
F_plus = np.polyval(r_d2, crit+h_step) / np.polyval(r_poly, crit+h_step)
F_0 = np.polyval(r_d2, crit) / np.polyval(r_poly, crit)
F_minus = np.polyval(r_d2, crit-h_step) / np.polyval(r_poly, crit-h_step)
Fpp_crit = (F_plus - 2*F_0 + F_minus) / h_step**2

print(f"\nsum h_k/(c-nu_k)^3 = {h_term:.6f}")
print(f"(1/4)*F''(c) = {0.25*Fpp_crit:.6f}")
print(f"Match: {abs(h_term - 0.25*Fpp_crit) < 0.1}")

# So <h,alpha> = m1 * sum h_k/(crit-nu_k)^3 = m1 * (1/4)*F''(crit)
print(f"\n<h,alpha> = m1 * (1/4)*F''(crit) = {m1 * 0.25 * Fpp_crit:.6f}")
print(f"<h,alpha> direct = {h_alpha_true:.6f}")

print(f"\nm1 = {m1:.6f}")
print(f"F''(crit) = {Fpp_crit:.6f}")
print(f"m1 * F''(crit) / 4 = {m1 * Fpp_crit / 4:.6f}")

# If F''(crit) < 0 and m1 > 0, then m1*F''(crit)/4 < 0
# But <h,alpha> should be >= 0. CHECK!

# Hmm: let's also check the sign of the identity for n=2 which we know exactly.

print("\n" + "=" * 70)
print("EXACT CHECK FOR n=2")
print("=" * 70)

# For n=2: p = (x-a)(x-b), q = (x-c)(x-d), r has gap_r = sqrt(gap_p^2 + gap_q^2)
# <h,alpha> = 2*(gap_r - gap_p)/(gap_r^2 * gap_p) > 0
# This is POSITIVE because gap_r > gap_p.

gap_p = b - a  # = 2
gap_q = q_roots_test[1] - q_roots_test[0]  # = 1
gap_r = r_roots_test[1] - r_roots_test[0]
print(f"gap_p = {gap_p}, gap_q = {gap_q}, gap_r = {gap_r}")
print(f"Expected gap_r = sqrt(gap_p^2 + gap_q^2) = {np.sqrt(gap_p**2+gap_q**2):.6f}")
print(f"<h,alpha> = 2*(gap_r-gap_p)/(gap_r^2*gap_p) = {2*(gap_r-gap_p)/(gap_r**2*gap_p):.6f}")
print(f"<h,alpha> direct = {h_alpha_true:.6f}")

# The issue is that I have: <h,alpha> = m1 * (1/4)*F''(crit) with m1 > 0
# but F''(crit) might be positive for n=2.
print(f"\nFor this n=2 case:")
print(f"  m1 = {m1:.6f}")
print(f"  F''(crit) = {Fpp_crit:.6f}")
print(f"  m1 * F''(crit) / 4 = {m1 * Fpp_crit / 4:.6f}")

# OHHH! For n=2, F''(crit) > 0 and m1 > 0, so the product is POSITIVE!
# The earlier test was for ARBITRARY roots of r, not for roots coming from MSS convolution.
# But wait, the formula <h,alpha> = m1*F''(crit)/4 uses the poles of omega_1 AT the
# critical points of r, and this F'' is evaluated on r which is the SPECIFIC polynomial
# from the MSS convolution. There's no constraint on r here -- the formula should hold
# for ANY real-rooted polynomial r with simple roots and ANY Herglotz omega_1.

# Let me re-run the test but with EXPLICIT n=3 examples from MSS convolution.

print("\n" + "=" * 70)
print("F''(c_j) for SPECIFIC MSS convolution examples")
print("=" * 70)

for n in [3, 4, 5]:
    np.random.seed(42 + n)

    for trial in range(5):
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
        alpha = H_p - H_r
        ha = np.dot(H_r, alpha)

        r_poly = np.poly(r_roots)
        r_d1 = np.polyder(r_poly)
        r_d2 = np.polyder(r_d1)

        crit_pts = np.sort(np.real(np.roots(r_d1)))

        # Compute m_j and F''(c_j) at each critical point
        print(f"\n  n={n}, trial {trial}: <h,alpha> = {ha:.6f}")
        total_from_formula = 0.0
        for j, c_j in enumerate(crit_pts):
            m_j = -n * np.polyval(r_poly, c_j) / np.polyval(r_d2, c_j)

            h_step = 1e-6
            F_plus = np.polyval(r_d2, c_j+h_step) / np.polyval(r_poly, c_j+h_step)
            F_0 = np.polyval(r_d2, c_j) / np.polyval(r_poly, c_j)
            F_minus = np.polyval(r_d2, c_j-h_step) / np.polyval(r_poly, c_j-h_step)
            Fpp = (F_plus - 2*F_0 + F_minus) / h_step**2

            contribution = m_j * Fpp / 4
            total_from_formula += contribution
            print(f"    c_{j} = {c_j:.4f}: m_j = {m_j:.6f}, F''(c_j) = {Fpp:.4f}, "
                  f"contribution = {contribution:.6f}")

        print(f"    Total from formula = {total_from_formula:.6f}")
        print(f"    Direct <h,alpha> = {ha:.6f}")
        print(f"    Match: {abs(total_from_formula - ha) < abs(ha)*0.1 + 0.01}")

        # Also check: these are NOT the actual poles of omega_1!
        # The poles of omega_1 are at the critical points of r, but the RESIDUES
        # m_j should come from the subordination equation, not from the simple formula.
        # Let me verify whether the critical points ARE the actual poles.

        # For the actual poles: we need to solve the subordination equation and find
        # where omega_1 diverges. We showed analytically that omega_1 diverges where
        # r'(z) = 0, which IS the critical points. But the residue formula
        # m_j = n*r(c_j)/r''(c_j) came from the leading-order expansion.
        # Wait: we had m_j = -n*r(c_j)/r''(c_j) (with the sign for the Herglotz representation
        # omega_1(z) = z + const + sum m_j/(p_j - z)), so
        # omega_1 near c_j: omega_1 ~ m_j/(c_j - z)
        # But our analytical computation gave omega_1 ~ n*r(c_j)/(r''(c_j)*(z-c_j))
        # = -n*r(c_j)/(r''(c_j)*(c_j-z))
        # So m_j = -n*r(c_j)/r''(c_j) in the convention omega_1 ~ m_j/(c_j-z).
        # Since r(c_j)/r''(c_j) < 0 (extremum), -r(c_j)/r''(c_j) > 0, so m_j > 0. Good.
        #
        # But ACTUALLY, this is the correct m_j only if omega_1 has NO other poles.
        # For degree-n omega_1, it has n-1 poles (the n-1 critical points of r).
        # Plus the "regular" part z + const. So the Herglotz decomposition is:
        # omega_1(z) = z + const + sum_{j=1}^{n-1} m_j/(c_j - z)
        # This has n-1 poles at c_j. The residues m_j are uniquely determined by
        # the polynomial equation for omega_1. The leading-order asymptotic
        # omega_1 ~ m_j/(c_j - z) near c_j gives the EXACT residue (since the other
        # terms are regular at c_j).
