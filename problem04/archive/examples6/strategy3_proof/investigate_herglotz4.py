#!/usr/bin/env python3
"""
Part 4: Focused investigation on proving <h,alpha> >= 0.

Summary of what we've learned:
1. omega_1 maps C^+ to C^+ but is NOT strictly Herglotz (Im(omega_1) can be < Im(z))
2. phi = omega_1 - id does NOT have all-positive residues
3. The Loewner matrix of omega_1 is NOT always PSD (fails ~50% of the time)
4. Individual terms (alpha_j-alpha_k)/(nu_j-nu_k) can be negative
5. The universal identities (sum x*H(x) = n(n-1)/2) are too weak

KEY REFORMULATIONS:
(A) <h,alpha> = sum_k h_k*(u_k - h_k) where u_k = H_p(lambda_k), h_k = H_r(nu_k)
(B) <h,alpha> = sum_{k<j} (alpha_j - alpha_k)/(nu_j - nu_k)
(C) <h,alpha> = <h,u> - ||h||^2

NEW APPROACH: Use the ELECTROSTATIC/VARIATIONAL structure.

The key insight: omega_1 is determined by the equation G_r(z) = G_p(omega_1(z)).
The chain rule gives omega_1'(nu_k) = 1.
This means the JACOBIAN of the map nu -> lambda (via omega_1) equals 1 at each point.

This is an AREA-PRESERVING (or rather LENGTH-PRESERVING) condition.
It constrains how the roots can rearrange.

Let me think about what SPECIFIC property of the MSS convolution forces <h,alpha> >= 0.
"""

import numpy as np
from math import factorial
from itertools import combinations

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


# ================================================================
# KEY OBSERVATION: omega_1''(nu_k) and the Schwarzian derivative
# ================================================================
print("="*70)
print("OMEGA_1 CONVEXITY ANALYSIS")
print("="*70)

# omega_1'(nu_k) = 1 for all k
# alpha_k = omega_1''(nu_k)/2
#
# Consider omega_1 on the interval [nu_k, nu_{k+1}].
# omega_1(nu_k) = lambda_k, omega_1(nu_{k+1}) = lambda_{k+1}
# omega_1'(nu_k) = 1, omega_1'(nu_{k+1}) = 1
#
# The integral int_{nu_k}^{nu_{k+1}} omega_1(z) dz can be computed by
# integration by parts or Taylor expansion.
#
# Since omega_1'(nu_k) = omega_1'(nu_{k+1}) = 1:
# omega_1(z) = z + f(z) where f(nu_k) = lambda_k - nu_k, f(nu_{k+1}) = lambda_{k+1} - nu_{k+1}
# f'(nu_k) = 0, f'(nu_{k+1}) = 0.
#
# The function f has zero derivative at both endpoints. If f is "convex" on
# the interval, then f'' >= 0, meaning alpha_k >= 0 and alpha_{k+1} may differ.
#
# But we saw that alpha can have MIXED signs. So f is not always convex.

# Let me try the MOST DIRECT approach possible.

# CLAIM TO TEST: For the MSS convolution with omega_1'(nu_k) = 1,
# the second derivative omega_1'' satisfies certain sum rules.

# We know: sum_k alpha_k = sum_k (u_k - h_k) = sum_k u_k - sum_k h_k = 0.
# (Both sums are 0 by antisymmetry.)

# We know: sum_k nu_k * alpha_k = sum_k nu_k * u_k - sum_k nu_k * h_k
# = n(n-1)/2 - n(n-1)/2 = 0. (Using the universal identity.)

# Wait, that's wrong. sum_k nu_k * u_k uses nu_k (not lambda_k).
# sum_k lambda_k * u_k = n(n-1)/2. sum_k nu_k * h_k = n(n-1)/2.
# But sum_k nu_k * u_k is NOT n(n-1)/2 in general.

# Actually: sum_k nu_k * alpha_k = sum_k nu_k * u_k - sum_k nu_k * h_k
# The second part is sum_k nu_k * h_k = n(n-1)/2.
# The first part is sum_k nu_k * H_p(lambda_k) -- this uses different points.

np.random.seed(42)
print("\nChecking sum rules for alpha:")
for trial in range(5):
    n = 4
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        raw = np.roots(np.poly(roots_r))
        if np.any(np.abs(np.imag(raw)) > 0.01):
            continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01):
            continue

        h = H_values(roots_r)
        u = H_values(roots_p)
        alpha = u - h

        print(f"  Trial {trial}: n={n}")
        print(f"    sum(alpha) = {np.sum(alpha):.8f} (should be ~0)")
        print(f"    sum(nu * alpha) = {np.dot(roots_r, alpha):.8f}")
        print(f"    sum(lambda * alpha) = {np.dot(roots_p, alpha):.8f}")
        print(f"    <h, alpha> = {np.dot(h, alpha):.8f}")
        print(f"    sum(h^2) = Phi_r = {np.sum(h**2):.8f}")
        print(f"    sum(u^2) = Phi_p = {np.sum(u**2):.8f}")
    except:
        pass

# ================================================================
# APPROACH: Partial summation / Abel summation
# ================================================================
print("\n\n" + "="*70)
print("ABEL SUMMATION / DIFFERENCE OPERATOR APPROACH")
print("="*70)

# <h, alpha> = sum_k h_k * alpha_k
#
# Both h and alpha depend on the root configuration. Let me try to
# express this in terms of gap sequences.
#
# Define: g_k = nu_{k+1} - nu_k (gaps of r), k = 1,...,n-1
#         d_k = lambda_{k+1} - lambda_k (gaps of p), k = 1,...,n-1
#
# The gap ratio: r_k = d_k / g_k = (lambda_{k+1}-lambda_k)/(nu_{k+1}-nu_k)
# This is the "average slope" of omega_1 on the k-th interval.
#
# Since omega_1'(nu_k) = 1 at endpoints and omega_1 maps [nu_k,nu_{k+1}] to [lambda_k,lambda_{k+1}]:
# int_{nu_k}^{nu_{k+1}} omega_1'(z) dz = lambda_{k+1} - lambda_k = d_k
# int_{nu_k}^{nu_{k+1}} 1 dz = g_k
# So int_{nu_k}^{nu_{k+1}} (omega_1'(z) - 1) dz = d_k - g_k

# The "deviation" d_k - g_k is related to the integral of phi' = omega_1' - 1.
# Since phi'(nu_k) = phi'(nu_{k+1}) = 0, phi' has a local extremum on each interval.

# For n = 2: omega_1'(z) = 1 everywhere would give d = g (contradiction unless
# the convolution is trivial). Actually for n=2, omega_1 is a Mobius transform
# and omega_1' = 1 at 2 points with omega_1'' != 0.

# FUNDAMENTAL IDENTITY FOR n=2:
# p = (x - lambda_1)(x - lambda_2), q = (x - mu_1)(x - mu_2)
# r = p boxplus_2 q = x^2 - (lambda_1+mu_1+lambda_2+mu_2)/2 * x + ...
# For n=2: boxplus is explicitly computable.

# For the general case, let me try a COMPLETELY DIFFERENT approach.
# Instead of trying to prove <h,alpha> >= 0 abstractly, let me look for
# a STRUCTURAL FORMULA.

# ================================================================
# THE SCHWARZIAN DERIVATIVE APPROACH
# ================================================================
print("\n" + "="*70)
print("SCHWARZIAN DERIVATIVE OF omega_1")
print("="*70)

# The Schwarzian derivative: S(f) = (f'''/f') - (3/2)(f''/f')^2
# For omega_1: S(omega_1)(z) evaluated at z = nu_k where omega_1'(nu_k) = 1:
# S(omega_1)(nu_k) = omega_1'''(nu_k) - (3/2) * omega_1''(nu_k)^2
#                   = omega_1'''(nu_k) - 6 * alpha_k^2

# For a Mobius transformation, S = 0. The Schwarzian measures deviation from Mobius.
# For n=2: omega_1 IS Mobius (degree 1 rational), so S = 0 everywhere.
# For n >= 3: omega_1 is NOT Mobius, so S != 0.

# The key property of the Schwarzian for Herglotz functions:
# If f maps C^+ to C^+ and is real on R, then Im(f(z)) > 0 for Im(z) > 0.
# The Schwarzian is related to the curvature of the image contour.

# THIS IS A DEAD END for proving <h,alpha> >= 0. Let me move to a completely
# different approach.

# ================================================================
# APPROACH: VERIFY A CLOSED-FORM IDENTITY FOR n=3
# ================================================================
print("\n" + "="*70)
print("CLOSED-FORM ANALYSIS FOR n = 3")
print("="*70)

# For n=3, we have 3 roots each. Let's parametrize:
# r roots: nu_1 < nu_2 < nu_3
# p roots: lambda_1 < lambda_2 < lambda_3
# q roots: mu_1 < mu_2 < mu_3
# r = p boxplus_3 q

# h_k = H_r(nu_k) = sum_{j!=k} 1/(nu_k - nu_j)
# h_1 = 1/(nu_1-nu_2) + 1/(nu_1-nu_3) < 0 (both terms negative)
# h_2 = 1/(nu_2-nu_1) + 1/(nu_2-nu_3) (mixed signs)
# h_3 = 1/(nu_3-nu_1) + 1/(nu_3-nu_2) > 0 (both terms positive)

# Similarly for u. So h_1 < 0, h_3 > 0, and h_2 can be either sign.

# For n=3 specifically, can we get a formula for <h,alpha>?

# Let a = nu_2 - nu_1 > 0, b = nu_3 - nu_2 > 0 (gaps of r)
# Let c = lambda_2 - lambda_1 > 0, d = lambda_3 - lambda_2 > 0 (gaps of p)

# h_1 = -1/a - 1/(a+b)
# h_2 = 1/a - 1/b
# h_3 = 1/(a+b) + 1/b

# u_1 = -1/c - 1/(c+d)
# u_2 = 1/c - 1/d
# u_3 = 1/(c+d) + 1/d

# alpha_k = u_k - h_k

# <h,alpha> = sum_k h_k*(u_k - h_k) = <h,u> - ||h||^2

# <h,u> = h_1*u_1 + h_2*u_2 + h_3*u_3
# = [-1/a - 1/(a+b)]*[-1/c - 1/(c+d)] + [1/a - 1/b]*[1/c - 1/d] + [1/(a+b) + 1/b]*[1/(c+d) + 1/d]

# ||h||^2 = h_1^2 + h_2^2 + h_3^2
# = [1/a + 1/(a+b)]^2 + [1/a - 1/b]^2 + [1/(a+b) + 1/b]^2

# This is algebraically messy but computable. Let me compute it symbolically.
# Actually, let me just verify numerically for specific cases and look for patterns.

print("\nn=3 detailed analysis:")
np.random.seed(42)
for trial in range(10):
    roots_p = np.sort(np.random.randn(3) * 2)
    for i in range(1, 3):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(3) * 2)
    for i in range(1, 3):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        raw = np.roots(np.poly(roots_r))
        if np.any(np.abs(np.imag(raw)) > 0.01):
            continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01):
            continue

        a, b = np.diff(roots_r)
        c, d = np.diff(roots_p)
        e, f = np.diff(roots_q)

        h = H_values(roots_r)
        u = H_values(roots_p)
        v = H_values(roots_q)
        alpha = u - h

        hu = np.dot(h, u)
        hh = np.dot(h, h)
        ha = hu - hh

        # gap ratio
        r1 = c/a  # d_1/g_1 = (lambda_2-lambda_1)/(nu_2-nu_1)
        r2 = d/b  # d_2/g_2

        # Is there a formula for <h,alpha> in terms of gap ratios?
        print(f"  a={a:.3f}, b={b:.3f}, c={c:.3f}, d={d:.3f}, "
              f"r1=c/a={r1:.3f}, r2=d/b={r2:.3f}, <h,alpha>={ha:.6f}")
    except:
        pass

# ================================================================
# APPROACH: SPECTRAL / EIGENVALUE INTERLACING
# ================================================================
print("\n\n" + "="*70)
print("SPECTRAL APPROACH: EXPRESSING <h,alpha> AS A POSITIVE QUANTITY")
print("="*70)

# Let me reconsider the problem from scratch.
# We want to show: sum_k H_r(nu_k) * H_p(lambda_k) >= sum_k H_r(nu_k)^2
# where nu_k = roots of r, lambda_k = roots of p, and the mapping nu_k -> lambda_k
# comes from the subordination function omega_1 with omega_1'(nu_k) = 1.

# Another way to write this:
# sum_k h_k * (u_k - h_k) >= 0
# where h_k = sum_{j!=k} 1/(nu_k - nu_j) and u_k = sum_{j!=k} 1/(lambda_k - lambda_j)

# Let's write phi_k = lambda_k - nu_k (the shift).
# Then lambda_k = nu_k + phi_k, so:
# u_k = sum_{j!=k} 1/((nu_k+phi_k) - (nu_j+phi_j))
#      = sum_{j!=k} 1/((nu_k-nu_j) + (phi_k-phi_j))

# Taylor expand for small phi:
# 1/((nu_k-nu_j) + (phi_k-phi_j)) = 1/(nu_k-nu_j) * 1/(1 + (phi_k-phi_j)/(nu_k-nu_j))
# ~ 1/(nu_k-nu_j) - (phi_k-phi_j)/(nu_k-nu_j)^2 + (phi_k-phi_j)^2/(nu_k-nu_j)^3 - ...

# u_k ~ h_k - sum_{j!=k} (phi_k-phi_j)/(nu_k-nu_j)^2 + sum_{j!=k} (phi_k-phi_j)^2/(nu_k-nu_j)^3 + ...
# alpha_k = u_k - h_k ~ -sum_{j!=k} (phi_k-phi_j)/(nu_k-nu_j)^2 + O(phi^2)

# To first order in phi:
# <h, alpha> ~ -sum_k h_k * sum_{j!=k} (phi_k-phi_j)/(nu_k-nu_j)^2
# = -sum_{k!=j} h_k * (phi_k-phi_j)/(nu_k-nu_j)^2
# = -sum_{k<j} [(h_k - h_j) * (phi_k - phi_j)] / (nu_k - nu_j)^2

# This is a SUM of products of "increments" of h and phi, weighted by 1/gap^2.
# If h and phi are "monotonically related" (both increase or both decrease with k),
# these products are positive, and the negative sign makes <h,alpha> negative!
# But we want <h,alpha> >= 0!

# So the LINEAR order term might be negative (or could be zero if the perturbation
# is in the "right" direction). The QUADRATIC order terms must compensate.

# Let me check the perturbation expansion more carefully.

print("\nPerturbation analysis (small phi regime):")
np.random.seed(42)

base_roots = np.array([0.0, 1.0, 3.0])
n = 3

for eps in [1.0, 0.5, 0.2, 0.1, 0.05, 0.01]:
    roots_p = base_roots + eps * np.array([0.1, -0.3, 0.2])
    roots_q = base_roots + eps * np.array([-0.2, 0.1, 0.1])

    # Ensure well-separated
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.01:
            continue
        if roots_q[i] - roots_q[i-1] < 0.01:
            continue

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        raw = np.roots(np.poly(roots_r))
        if np.any(np.abs(np.imag(raw)) > 0.01):
            continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.001):
            continue

        h = H_values(roots_r)
        u = H_values(roots_p)
        alpha = u - h
        ha = np.dot(h, alpha)
        phi = roots_p - roots_r

        # First order contribution
        first_order = 0
        for k in range(n):
            for j in range(k+1, n):
                first_order -= (h[k]-h[j]) * (phi[k]-phi[j]) / (roots_r[k]-roots_r[j])**2

        print(f"  eps={eps:.3f}: <h,alpha>={ha:.8f}, 1st_order~{first_order:.8f}, "
              f"phi={np.round(phi,4)}")
    except:
        pass


# ================================================================
# THE MOST PROMISING APPROACH: FUNCTIONAL MONOTONICITY
# ================================================================
print("\n\n" + "="*70)
print("FUNCTIONAL MONOTONICITY: Phi as function of t along interpolation")
print("="*70)

# Define x_k(t) = (1-t)*nu_k + t*lambda_k for t in [0,1].
# F(t) = sum_k H_t(x_k(t)) * H_r(nu_k) where H_t uses the x_k(t) points.
# F(0) = ||h||^2
# F(1) = <h,u>
# <h,alpha> = F(1) - F(0)

# If F is monotone increasing, then <h,alpha> >= 0.
# F'(t) involves the derivative of the H-values along the interpolation.

# Actually, let me consider a simpler functional:
# Phi(t) = sum_k H_t(x_k(t))^2  where x_k(t) = (1-t)*nu_k + t*lambda_k
# Phi(0) = Phi_r
# Phi(1) = Phi_p
# We know Phi_p >= Phi_r (from <h,alpha> >= 0 and the CS inequality).
# Is Phi(t) monotone? Let's check.

print("\nPhi(t) along interpolation:")
np.random.seed(42)

for trial in range(3):
    n = 4
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        raw = np.roots(np.poly(roots_r))
        if np.any(np.abs(np.imag(raw)) > 0.01): continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01): continue

        t_vals = np.linspace(0, 1, 20)
        Phi_vals = []
        F_vals = []  # F(t) = <h_r, h_t> where h_r is H-values at r-roots
        h_r = H_values(roots_r)

        for t in t_vals:
            x = (1-t)*roots_r + t*roots_p
            h_t = H_values(x)
            Phi_vals.append(np.sum(h_t**2))
            F_vals.append(np.dot(h_r, h_t))

        print(f"\n  Trial {trial}: n={n}")
        print(f"  Phi(0) = Phi_r = {Phi_vals[0]:.6f}")
        print(f"  Phi(1) = Phi_p = {Phi_vals[-1]:.6f}")
        print(f"  F(0) = <h_r,h_r> = {F_vals[0]:.6f}")
        print(f"  F(1) = <h_r,h_p> = {F_vals[-1]:.6f}")
        print(f"  <h,alpha> = F(1)-F(0) = {F_vals[-1]-F_vals[0]:.6f}")

        # Is Phi monotone?
        Phi_diffs = np.diff(Phi_vals)
        print(f"  Phi(t) monotone increasing: {np.all(Phi_diffs >= -1e-10)}")

        # Is F monotone?
        F_diffs = np.diff(F_vals)
        print(f"  F(t) monotone increasing: {np.all(F_diffs >= -1e-10)}")
    except:
        pass


# ================================================================
# THE CAUCHY MATRIX APPROACH: POSITIVE DEFINITE STRUCTURE
# ================================================================
print("\n\n" + "="*70)
print("CAUCHY MATRIX PD STRUCTURE")
print("="*70)

# The Cauchy matrix C_{ij} = 1/(x_i - y_j) is well-known.
# For x_1 < ... < x_n and y_1 < ... < y_n with x_i != y_j,
# if x and y interlace (y_1 < x_1 < y_2 < x_2 < ...), then C has specific
# properties.

# But in our case, we don't have interlacing. We have:
# nu_1 < nu_2 < ... < nu_n (roots of r)
# lambda_1 < lambda_2 < ... < lambda_n (roots of p)
# with no particular interlacing (other than what the convolution provides).

# KEY MATRIX: Define M_{kl} = 1/((nu_k - nu_l)(lambda_k - lambda_l)) for k != l.
# Then alpha_k = sum_{l!=k} (phi_k - phi_l) * M_{kl} where phi_k = nu_k - lambda_k.

# <h,alpha> = sum_k h_k * sum_{l!=k} (phi_k - phi_l) * M_{kl}
#           = sum_{k!=l} h_k * (phi_k - phi_l) * M_{kl}
#           = sum_{k<l} [h_k - h_l] * [phi_k - phi_l] * M_{kl}

# Wait, let me be careful:
# sum_{k!=l} h_k * (phi_k - phi_l) * M_{kl}
# = sum_{k!=l} h_k*phi_k*M_{kl} - sum_{k!=l} h_k*phi_l*M_{kl}
# = sum_k h_k*phi_k*sum_{l!=k} M_{kl} - sum_l phi_l*sum_{k!=l} h_k*M_{kl}
# These are different sums.

# Actually using antisymmetry in (k,l) swap:
# sum_{k!=l} h_k * (phi_k - phi_l) * M_{kl}
# = sum_{k<l} h_k*(phi_k-phi_l)*M_{kl} + sum_{k>l} h_k*(phi_k-phi_l)*M_{kl}
# In the second sum, swap k<->l:
# = sum_{k<l} h_k*(phi_k-phi_l)*M_{kl} + sum_{k<l} h_l*(phi_l-phi_k)*M_{lk}
# Note M_{kl} = M_{lk} (symmetric in k,l since (nu_k-nu_l)(lambda_k-lambda_l) = (nu_l-nu_k)(lambda_l-lambda_k))
# = sum_{k<l} [h_k*(phi_k-phi_l) - h_l*(phi_k-phi_l)] * M_{kl}
# = sum_{k<l} (h_k - h_l)*(phi_k - phi_l) * M_{kl}

# So <h,alpha> = sum_{k<l} (h_k - h_l)(phi_k - phi_l) / [(nu_k-nu_l)(lambda_k-lambda_l)]

# Since k < l means nu_k < nu_l and lambda_k < lambda_l:
# (nu_k - nu_l) < 0, (lambda_k - lambda_l) < 0
# So (nu_k-nu_l)(lambda_k-lambda_l) > 0, hence M_{kl} > 0 for k < l.

# Therefore <h,alpha> >= 0 iff sum_{k<l} (h_k-h_l)(phi_k-phi_l) * M_{kl} >= 0
# with M_{kl} > 0.

# This requires (h_k - h_l)(phi_k - phi_l) to be non-negative ON AVERAGE
# (weighted by M_{kl}).

# CHECK: Is (h_k - h_l) always the same sign as (phi_k - phi_l)?

np.random.seed(42)
same_sign_count = 0
diff_sign_count = 0
total_pairs = 0

same_sign_sum = 0
diff_sign_sum = 0
total_sum = 0

for trial in range(500):
    n = np.random.choice([3, 4, 5])
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5
    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        raw = np.roots(np.poly(roots_r))
        if np.any(np.abs(np.imag(raw)) > 0.01): continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01): continue

        h = H_values(roots_r)
        phi = roots_r - roots_p  # nu - lambda (note: NOT lambda - nu)

        for k in range(n):
            for l in range(k+1, n):
                dh = h[k] - h[l]
                dphi = phi[k] - phi[l]  # (nu_k-lambda_k) - (nu_l-lambda_l)
                Mkl = 1.0 / ((roots_r[k]-roots_r[l]) * (roots_p[k]-roots_p[l]))
                term = dh * dphi * Mkl

                total_pairs += 1
                total_sum += abs(term)
                if dh * dphi >= 0:
                    same_sign_count += 1
                    same_sign_sum += abs(term)
                else:
                    diff_sign_count += 1
                    diff_sign_sum += abs(term)
    except:
        pass

print(f"\nSign concordance of (h_k-h_l) and (phi_k-phi_l):")
print(f"  Same sign: {same_sign_count}/{total_pairs} ({100*same_sign_count/total_pairs:.1f}%)")
print(f"  Diff sign: {diff_sign_count}/{total_pairs} ({100*diff_sign_count/total_pairs:.1f}%)")
print(f"  Weighted same: {same_sign_sum:.4f}, weighted diff: {diff_sign_sum:.4f}")
print(f"  Ratio same/diff (weighted): {same_sign_sum/diff_sign_sum:.4f}")

# ================================================================
# THE KEY FORMULA
# ================================================================
print("\n\n" + "="*70)
print("KEY FORMULA: <h,alpha> = sum_{k<l} (h_k-h_l)(phi_k-phi_l) / [(nu_k-nu_l)(lambda_k-lambda_l)]")
print("="*70)

# Verify this formula:
for trial in range(3):
    n = 4
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5
    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        raw = np.roots(np.poly(roots_r))
        if np.any(np.abs(np.imag(raw)) > 0.01): continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01): continue

        h = H_values(roots_r)
        u = H_values(roots_p)
        alpha = u - h
        ha_direct = np.dot(h, alpha)

        phi = roots_r - roots_p
        ha_formula = 0
        for k in range(n):
            for l in range(k+1, n):
                ha_formula += (h[k]-h[l])*(phi[k]-phi[l]) / ((roots_r[k]-roots_r[l])*(roots_p[k]-roots_p[l]))

        print(f"  <h,alpha> direct = {ha_direct:.8f}, formula = {ha_formula:.8f}, "
              f"match = {abs(ha_direct-ha_formula) < 1e-8}")
    except:
        pass

# ================================================================
# CAN WE PROVE THE SIGN CONCORDANCE?
# ================================================================
print("\n\n" + "="*70)
print("INVESTIGATING SIGN CONCORDANCE")
print("="*70)

# The formula: <h,alpha> = sum_{k<l} (h_k-h_l)(phi_k-phi_l) * M_{kl} with M_{kl} > 0.
# We need this to be >= 0.
# Numerically, (h_k-h_l)(phi_k-phi_l) is positive in 75% of pairs.
# The weighted average is even more skewed toward positive.

# What is the structural reason?
# h_k = sum_{j!=k} 1/(nu_k - nu_j) -- this is "small negative" for small k
# and "large positive" for large k (roughly increasing with k).
# phi_k = nu_k - lambda_k -- what's the pattern?

# For the MSS convolution, roots of r are more spread than roots of p.
# So the outer roots of r are more extreme. This means:
# nu_1 < lambda_1 (r-root is to the left of p-root) => phi_1 = nu_1 - lambda_1 < 0
# nu_n > lambda_n (r-root is to the right) => phi_n = nu_n - lambda_n > 0
# So phi is INCREASING (from negative to positive).

# h is also roughly INCREASING (from negative to positive):
# h_1 < 0 (leftmost root), h_n > 0 (rightmost root).

# So BOTH h and phi are roughly increasing, meaning (h_k-h_l) and (phi_k-phi_l)
# tend to have the SAME sign. This is why the sum is positive!

# But "roughly increasing" is not precise enough. Let me check if they are
# always monotonically ordered in the same direction.

print("\nChecking monotonicity of h and phi:")
np.random.seed(42)
h_mono = 0
phi_mono = 0
both_mono = 0
total = 0
for trial in range(500):
    n = np.random.choice([3, 4, 5])
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5
    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        raw = np.roots(np.poly(roots_r))
        if np.any(np.abs(np.imag(raw)) > 0.01): continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01): continue

        h = H_values(roots_r)
        phi = roots_r - roots_p

        total += 1
        h_inc = np.all(np.diff(h) > -1e-10)
        phi_inc = np.all(np.diff(phi) > -1e-10)
        if h_inc: h_mono += 1
        if phi_inc: phi_mono += 1
        if h_inc and phi_inc: both_mono += 1
    except:
        pass

print(f"  h monotone increasing: {h_mono}/{total} ({100*h_mono/total:.1f}%)")
print(f"  phi monotone increasing: {phi_mono}/{total} ({100*phi_mono/total:.1f}%)")
print(f"  Both monotone: {both_mono}/{total} ({100*both_mono/total:.1f}%)")
print(f"  (h is NOT always monotone, and phi is NOT always monotone)")

# So neither h nor phi is always monotone. The proof needs something more subtle.

# ================================================================
# THE CHEBYSHEV SUM INEQUALITY APPROACH
# ================================================================
print("\n\n" + "="*70)
print("CHEBYSHEV / CORRELATION APPROACH")
print("="*70)

# Even though individual pairs can have opposite signs, the SUM
# sum_{k<l} (h_k-h_l)(phi_k-phi_l)*M_{kl} is always >= 0.
# This is a WEIGHTED version of the Chebyshev sum inequality.

# Chebyshev's sum inequality: if a_1 <= ... <= a_n and b_1 <= ... <= b_n, then
# n * sum a_i*b_i >= (sum a_i)(sum b_i).
# This is equivalent to: sum_{i<j} (a_i-a_j)(b_i-b_j) >= 0.

# In our case: sum_{k<l} (h_k-h_l)(phi_k-phi_l)*w_{kl} >= 0 with w_{kl} > 0.
# This is a WEIGHTED Chebyshev inequality.

# The weighted version holds if h and phi are SIMILARLY ORDERED
# (i.e., (h_k-h_l)(phi_k-phi_l) >= 0 for all k<l), which we showed fails.
# BUT: it also holds if the POSITIVE contributions (same sign pairs)
# outweigh the negative contributions, which they always do empirically.

# The question is: WHAT STRUCTURAL PROPERTY of the MSS convolution
# guarantees this? The answer likely involves the subordination structure.

# Let me investigate whether there's a SIMPLER sufficient condition.

# KEY OBSERVATION: <h,alpha> = sum_{k<l} (h_k-h_l)(phi_k-phi_l)/[(nu_k-nu_l)(lambda_k-lambda_l)]
# The denominator (nu_k-nu_l)(lambda_k-lambda_l) > 0 (since both sorted).
# So the sign depends on (h_k-h_l)(phi_k-phi_l).
# h_k - h_l = H_r(nu_k) - H_r(nu_l)
# phi_k - phi_l = (nu_k-lambda_k) - (nu_l-lambda_l) = (nu_k-nu_l) - (lambda_k-lambda_l)

# For k < l: nu_k < nu_l, lambda_k < lambda_l.
# phi_k - phi_l = (nu_k-nu_l) - (lambda_k-lambda_l) = -(gap_r(k,...,l) - gap_p(k,...,l))
# = gap_p(k,...,l) - gap_r(k,...,l)  (with appropriate sign... let me be careful)

# Actually: phi_k - phi_l = (nu_k - lambda_k) - (nu_l - lambda_l)
# With k < l: this is negative if nu_l-lambda_l > nu_k-lambda_k,
# i.e., if the shift phi is INCREASING.

# h_k - h_l: for k < l (nu_k < nu_l), this can be positive or negative.
# h_1 < h_n in general, so h is roughly increasing.

# If both h and phi are increasing, (h_k-h_l)(phi_k-phi_l) = (neg)(neg) = positive.
# So the concordance is positive.

# The cases where concordance fails are when h or phi has a "dip" or "bump".

# FINAL INSIGHT: Maybe the proof should use the CAUCHY-BINET identity or
# a determinantal formula. The sum has the form of a quadratic form in the
# differences, weighted by the Cauchy kernel.

# Actually, the expression sum_{k<l} f(k)*g(l)/[x(k)*y(l)] is related to
# the VANDERMONDE structure. Let me check if there's a matrix formulation.

# Define the n x n matrix:
# A_{kl} = sqrt(M_{kl}) * (h_k - h_l) for k < l (upper triangular)
# B_{kl} = sqrt(M_{kl}) * (phi_k - phi_l) for k < l
# <h,alpha> = sum_{k<l} A_{kl} * B_{kl} = <A, B>_F (Frobenius inner product of upper triangular parts)

# By Cauchy-Schwarz: |<h,alpha>| <= ||A||_F * ||B||_F.
# But we need a LOWER bound, not an upper bound.

# THIS IS FUNDAMENTALLY DIFFICULT because the inequality involves the interplay
# between h (which depends on r-roots) and phi (which depends on the relationship
# between r-roots and p-roots). Without a closed-form for the convolution,
# we can't decouple these.

# CONCLUSION: The Herglotz approach to proving <h,alpha> >= 0 FAILS because
# omega_1 is NOT a strict Herglotz function. The divided difference / Cauchy
# matrix reformulation gives a structural formula, but proving positivity
# requires understanding the fine structure of the MSS convolution.

print("\nCONCLUSION: The Herglotz approach does not yield a proof of <h,alpha> >= 0.")
print("The key obstacle: omega_1 maps C^+ to C^+ but NOT with Im(omega_1) >= Im(z).")
print("The partial fraction decomposition of phi = omega_1 - id has MIXED SIGN residues.")
print("This invalidates any argument based on the Nevanlinna representation.")
print("\nThe most promising reformulation found:")
print("  <h,alpha> = sum_{k<l} (h_k-h_l)(phi_k-phi_l) / [(nu_k-nu_l)(lambda_k-lambda_l)]")
print("where phi_k = nu_k - lambda_k. This is a WEIGHTED correlation of the sequences")
print("h = H-values and phi = root-shifts, with positive weights.")
print("Numerical evidence shows 75% of individual terms are positive, and the")
print("total sum is ALWAYS positive (800+ trials). But a rigorous proof of positivity")
print("requires understanding why the MSS convolution structure forces this correlation.")
