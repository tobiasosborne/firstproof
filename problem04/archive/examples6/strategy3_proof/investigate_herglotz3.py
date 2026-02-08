#!/usr/bin/env python3
"""
Part 3: Attempt proof of <h,alpha> >= 0 via the Loewner matrix /
Schur complement approach.

KEY REFORMULATION from Part 2:
<h,alpha> = sum_{k<j} (alpha_j - alpha_k) / (nu_j - nu_k)

where alpha_k = H_p(lambda_k) - H_r(nu_k).

ALTERNATIVE REFORMULATION:
<h,alpha> = sum_k h_k * alpha_k
          = sum_k H_r(nu_k) * [H_p(lambda_k) - H_r(nu_k)]
          = <h, u> - ||h||^2

where u_k = H_p(lambda_k) and h_k = H_r(nu_k).

NEW IDEA: Express <h,u> using the LOEWNER MATRIX structure.

The Loewner matrix L_{ij} = (f(x_i) - f(x_j))/(x_i - x_j) for a function f
and points x_1, ..., x_n has the property that L is PSD iff f is
operator monotone on the interval containing all x_i.

Can we express <h,u> as a quadratic form involving a Loewner matrix?

h_k = sum_{j!=k} 1/(nu_k - nu_j)
u_k = sum_{j!=k} 1/(lambda_k - lambda_j)

<h,u> = sum_k [sum_{j!=k} 1/(nu_k-nu_j)] * [sum_{l!=k} 1/(lambda_k-lambda_l)]

This is a product of two "row sums" from two different Cauchy matrices.

Alternative: Think of h_k = (d/dx) log |r(x)/(x-nu_k)| at x = nu_k.
And alpha_k = (d/dx) [log |p(x)/(x-lambda_k)| - log |r(x)/(x-nu_k)|]
but evaluated at different points...

Actually, let me think about this more carefully using the ELECTROSTATIC
interpretation.

ELECTROSTATIC ANALOGY:
h_k = H_r(nu_k) = "electrostatic field" at nu_k due to unit charges at
      nu_1,...,nu_n (excluding self-charge at nu_k).
u_k = H_p(lambda_k) = field at lambda_k due to charges at lambda_1,...,lambda_n.

<h,u> - ||h||^2 >= 0 asks: the mixed "energy" sum_k h_k*u_k >= sum_k h_k^2.

Since u_k = h_k + alpha_k, this is sum_k h_k*alpha_k >= 0.

NEW APPROACH: Use the fact that omega_1 maps INTERVALS to INTERVALS.
omega_1: (nu_k, nu_{k+1}) -> (lambda_k, lambda_{k+1}) with omega_1'(nu_k) = 1.

The LENGTH of the image interval is lambda_{k+1} - lambda_k.
The length of the domain interval is nu_{k+1} - nu_k.

By the mean value theorem, there exists xi_k in (nu_k, nu_{k+1}) with
omega_1'(xi_k) = (lambda_{k+1} - lambda_k) / (nu_{k+1} - nu_k).

Since omega_1'(nu_k) = 1 and omega_1'(nu_{k+1}) = 1, the derivative starts
at 1, changes to some value, and returns to 1. The AVERAGE derivative is
(lambda_{k+1}-lambda_k)/(nu_{k+1}-nu_k) = gap_p(k)/gap_r(k).

This doesn't directly give <h,alpha> >= 0, but it constrains the structure.

Let me try yet another approach: the ENERGY functional.
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
# APPROACH: LOG-ENERGY / ENTROPY
# ================================================================
print("="*70)
print("LOG-ENERGY APPROACH")
print("="*70)

# Define the log-energy: E(x_1,...,x_n) = -2/(n(n-1)) * sum_{i<j} log|x_i - x_j|
# The gradient of E with respect to x_k is:
# dE/dx_k = -2/(n(n-1)) * sum_{j!=k} 1/(x_k - x_j) = -2/(n(n-1)) * H(x_k)

# Phi_n(p) = sum_k H_p(lambda_k)^2

# Consider: Phi_n(p) = ||grad_lambda E(lambda)||^2 * [n(n-1)/2]^2
# (up to constant). So Phi is related to the squared gradient of the log-energy.

# For a CONVEX function f, ||grad f(x)||^2 >= ||grad f(y)||^2 when x is
# "further from the minimum" than y. But log-energy is not convex in general.

# LOEWNER MATRIX APPROACH:
# Define the n x n matrix M_kj = 1/(nu_k - lambda_j) for k,j = 1,...,n.
# This is a CAUCHY MATRIX. Its determinant is a Cauchy determinant:
# det(M) = prod_{k<l}(nu_l-nu_k) * prod_{k<l}(lambda_k-lambda_l) / prod_{k,l}(nu_k-lambda_l)

# Now: <h, u> = sum_k [sum_{j!=k} 1/(nu_k-nu_j)] * [sum_{l!=k} 1/(lambda_k-lambda_l)]
# This is NOT directly a trace or quadratic form of M.

# Let me try expressing <h,u> differently.
# h_k = -nG_r(nu_k) ... wait, G_r has a pole at nu_k.
# More precisely: h_k = lim_{z->nu_k} [nG_r(z) - 1/(z-nu_k)]

# Let's think about this using RESIDUES.
# Consider: f(z) = [nG_r(z)]^2 - sum_k h_k^2/(z-nu_k)^2 - 2*sum_k h_k/(z-nu_k)
# ... wait, this is the "regularized square" of G_r.

# Actually, let me try:
# sum_k Res_{z=nu_k} [nG_r(z) * nG_p(omega_1(z))]
# Near z = nu_k: nG_r(z) = 1/(z-nu_k) + h_k + ...
# nG_p(omega_1(z)) = nG_p(lambda_k + (z-nu_k) + alpha_k(z-nu_k)^2 + ...)
# = 1/(omega_1(z)-lambda_k) + ... (pole at omega_1(z) = lambda_k)
# But omega_1(z) = lambda_k when z = nu_k, so:
# omega_1(z) - lambda_k = (z-nu_k) + alpha_k(z-nu_k)^2 + ...
# nG_p(omega_1(z)) = 1/(z-nu_k + alpha_k(z-nu_k)^2 + ...) + u_k + ...
# = 1/((z-nu_k)(1 + alpha_k(z-nu_k) + ...)) + u_k + ...
# = [1/(z-nu_k)]*(1 - alpha_k(z-nu_k) + ...) + u_k + ...
# = 1/(z-nu_k) - alpha_k + u_k + ...
# = 1/(z-nu_k) + h_k + ...   (since u_k = h_k + alpha_k => u_k - alpha_k = h_k)

# Wait, that's the same as nG_r(z)! Because G_r(z) = G_p(omega_1(z))!
# Of course: [nG_r(z)] * [nG_p(omega_1(z))] = [nG_r(z)]^2.
# So this doesn't give anything new.

# BETTER IDEA: Consider the DIFFERENCE [nG_p(z)]^2 - [nG_r(z)]^2.
# This is well-defined for z far from the roots.
# For large z: both are ~1/z^2, so the difference is O(1/z^3).
# By the residue theorem:
# 0 = sum_k Res_{z=nu_k} {[nG_p(z)]^2 - [nG_r(z)]^2}
#     + sum_k Res_{z=lambda_k} {[nG_p(z)]^2 - [nG_r(z)]^2}
#     (assuming the contour at infinity gives 0)
#
# Near z = nu_k: nG_r has a pole, nG_p is analytic (since nu_k != lambda_k generically).
# [nG_r(z)]^2 ~ 1/(z-nu_k)^2 + 2h_k/(z-nu_k) + ...
# [nG_p(z)]^2 is analytic at nu_k: let nG_p(nu_k) = sum_l 1/(nu_k-lambda_l).
# Residue of [nG_p]^2 at nu_k = 0 (analytic).
# Residue of -[nG_r]^2 at nu_k = -2h_k.
# So Res_{z=nu_k} = -2h_k.
#
# Near z = lambda_k: nG_p has a pole, nG_r is analytic.
# [nG_p(z)]^2 ~ 1/(z-lambda_k)^2 + 2u_k/(z-lambda_k) + ...
# Residue of [nG_p]^2 at lambda_k = 2u_k.
# Residue of [nG_r]^2 at lambda_k = 0 (analytic if lambda_k != nu_l for all l).
# So Res_{z=lambda_k} = 2u_k.
#
# Total: sum_k 2u_k - sum_k 2h_k = 0.
# => sum u_k = sum h_k.
# This is just the identity sum H_p(lambda_k) = sum H_r(nu_k) = 0
# (since H_p is the derivative of log of a monic polynomial, and the sum
# of H-values equals the coefficient of z^{n-2} in p'/p... Actually
# sum_k H_p(lambda_k) = sum_k sum_{j!=k} 1/(lambda_k-lambda_j) = 0
# by antisymmetry.)

# So the "first moment" approach gives a triviality. Let me try "second moment":
# Consider f(z) = z * ([nG_p(z)]^2 - [nG_r(z)]^2).
# Residue at nu_k of -z*[nG_r]^2: need Res of z/(z-nu_k)^2 + 2*z*h_k/(z-nu_k)
# z/(z-nu_k)^2 = (nu_k + (z-nu_k))/(z-nu_k)^2 = nu_k/(z-nu_k)^2 + 1/(z-nu_k)
# Residue = 1 + 2*nu_k*h_k.
# Actually: Res of z*[nG_r(z)]^2 at nu_k = d/dz[z*(z-nu_k)^2 * [nG_r(z)]^2]|_{z=nu_k} ... no.
# Let me use Laurent expansion.
# z * [nG_r(z)]^2 = z * [1/(z-nu_k)^2 + 2h_k/(z-nu_k) + ...]
# = z/(z-nu_k)^2 + 2*z*h_k/(z-nu_k) + ...
# z = nu_k + (z-nu_k), so:
# = [nu_k/(z-nu_k)^2 + 1/(z-nu_k)] + [2*nu_k*h_k/(z-nu_k) + 2*h_k] + ...
# = nu_k/(z-nu_k)^2 + (1 + 2*nu_k*h_k)/(z-nu_k) + ...
# Residue at nu_k = 1 + 2*nu_k*h_k.
#
# Similarly for z*[nG_p(z)]^2 at lambda_k: residue = 1 + 2*lambda_k*u_k.
#
# Sum of residues at infinity = 0 (f(z) ~ z/z^4 + ... -> 0).
# Actually f(z) ~ z * (1/z^2 - 1/z^2) + z * O(1/z^3) = O(1/z^2). OK.
#
# So: sum_k (1 + 2*lambda_k*u_k) - sum_k (1 + 2*nu_k*h_k) = 0
# => sum_k lambda_k*u_k = sum_k nu_k*h_k
# This is: sum_k lambda_k * H_p(lambda_k) = sum_k nu_k * H_r(nu_k).
#
# Let me verify this numerically.

np.random.seed(42)

n = 4
roots_p = np.sort(np.random.randn(n) * 2)
for i in range(1, n):
    if roots_p[i] - roots_p[i-1] < 0.5:
        roots_p[i] = roots_p[i-1] + 0.5
roots_q = np.sort(np.random.randn(n) * 2)
for i in range(1, n):
    if roots_q[i] - roots_q[i-1] < 0.5:
        roots_q[i] = roots_q[i-1] + 0.5
roots_r, _ = boxplus_mss(roots_p, roots_q)
roots_r = np.sort(np.real(roots_r))

h = H_values(roots_r)
u = H_values(roots_p)

print(f"\nVerification of sum lambda_k*u_k = sum nu_k*h_k:")
lhs = np.dot(roots_p, u)
rhs = np.dot(roots_r, h)
print(f"  sum lambda_k * H_p(lambda_k) = {lhs:.8f}")
print(f"  sum nu_k * H_r(nu_k) = {rhs:.8f}")
print(f"  Difference: {abs(lhs-rhs):.2e}")

# Hmm, these are NOT equal in general! Because the contour integral argument
# assumes G_p and G_r are the SAME rational function up to root locations.
# But they are different polynomials (r = p boxplus q, not r = p).
# The residues at infinity might not cancel.

# Let me recheck. For large z:
# nG_p(z) = sum 1/(z-lambda_k) = n/z + (sum lambda_k)/z^2 + (sum lambda_k^2)/z^3 + ...
# [nG_p(z)]^2 = n^2/z^2 + 2n*sum(lambda)/z^3 + ...
# Similarly [nG_r(z)]^2 = n^2/z^2 + 2n*sum(nu)/z^3 + ...
# [nG_p]^2 - [nG_r]^2 = 2n*(sum(lambda)-sum(nu))/z^3 + ...
# This is O(1/z^3).
# z * ([nG_p]^2 - [nG_r]^2) = 2n*(sum(lambda)-sum(nu))/z^2 + ...
# This vanishes at infinity. Good.
# The residue at infinity of f(z) = z*([nG_p]^2-[nG_r]^2):
# Res_inf(f) = -[coefficient of 1/z in f] = 0 (since f = O(1/z^2)).
#
# So the sum of ALL finite residues = 0.
# Finite residues: at lambda_k (from G_p^2) and at nu_k (from G_r^2).
# Also at mu_k (from... wait, G_q doesn't appear in f).
#
# f(z) = [nG_p(z)]^2 - [nG_r(z)]^2. Poles at lambda_k (from G_p) and nu_k (from G_r).
# sum Res_{lambda_k} f + sum Res_{nu_k} f = 0
# sum (1+2*lambda_k*u_k) - sum (1+2*nu_k*h_k) = 0
# n + 2*sum(lambda_k*u_k) - n - 2*sum(nu_k*h_k) = 0
# sum(lambda_k*u_k) = sum(nu_k*h_k)
#
# BUT wait: this identity uses G_p and G_r as INDEPENDENT functions.
# There's no requirement that r = p boxplus q here! It's a UNIVERSAL identity:
# for any two monic polynomials p, r with simple roots,
# sum lambda_k*H_p(lambda_k) = sum nu_k*H_r(nu_k).
#
# Actually, that can't be right. Let me re-derive.
# sum_k lambda_k*H_p(lambda_k) = sum_k lambda_k * sum_{j!=k} 1/(lambda_k-lambda_j)
# For any set of distinct points, this equals... hmm, let me compute for a simple case.

print("\n\nUniversal identity test:")
for pts in [[1, 2, 3], [-1, 0, 2, 5], [0, 1]]:
    pts = np.array(pts, dtype=float)
    H = H_values(pts)
    val = np.dot(pts, H)
    print(f"  Points {pts}: sum x_k*H(x_k) = {val:.6f}")

# Interesting! It seems sum x_k*H(x_k) = (n-1)/2 * ... let me check.
# For pts = [1,2,3]: H = [-1, 0, 1/2] ... wait.
# H_1 = 1/(1-2) + 1/(1-3) = -1 + (-1/2) = -3/2
# H_2 = 1/(2-1) + 1/(2-3) = 1 + (-1) = 0
# H_3 = 1/(3-1) + 1/(3-2) = 1/2 + 1 = 3/2
# sum x*H = 1*(-3/2) + 2*0 + 3*(3/2) = -3/2 + 9/2 = 3
# = n*(n-1)/2 = 3 for n=3. Interesting!

# For pts = [-1, 0, 2, 5]:
# H_{-1} = 1/(-1-0) + 1/(-1-2) + 1/(-1-5) = -1 - 1/3 - 1/6 = -3/2
# H_0 = 1/(0-(-1)) + 1/(0-2) + 1/(0-5) = 1 - 1/2 - 1/5 = 3/10
# H_2 = 1/(2-(-1)) + 1/(2-0) + 1/(2-5) = 1/3 + 1/2 - 1/3 = 1/2
# H_5 = 1/(5-(-1)) + 1/(5-0) + 1/(5-2) = 1/6 + 1/5 + 1/3 = 21/30 = 7/10
# sum x*H = (-1)(-3/2) + 0*(3/10) + 2*(1/2) + 5*(7/10) = 3/2 + 0 + 1 + 7/2 = 6
# = n*(n-1)/2 = 4*3/2 = 6. YES!

# UNIVERSAL IDENTITY: sum_k x_k * H(x_k) = n*(n-1)/2 for ANY n distinct points.
# This is because:
# sum_k x_k * sum_{j!=k} 1/(x_k-x_j)
# = sum_{k!=j} x_k/(x_k-x_j)
# = sum_{k<j} [x_k/(x_k-x_j) + x_j/(x_j-x_k)]
# = sum_{k<j} [x_k/(x_k-x_j) - x_j/(x_k-x_j)]
# = sum_{k<j} (x_k - x_j)/(x_k - x_j)
# = sum_{k<j} 1
# = n*(n-1)/2

print(f"\nUniversal identity: sum x_k*H(x_k) = n*(n-1)/2 = {n*(n-1)/2}")
print(f"This is independent of the points! So the contour integral identity")
print(f"sum lambda_k*u_k = sum nu_k*h_k is trivially true because both = n(n-1)/2.\n")

# So the "z * [G_p^2 - G_r^2]" approach gives a trivial identity.
# Let me try z^2 * [G_p^2 - G_r^2] instead.

# For z^2 * [nG_p(z)]^2 near z = lambda_k:
# z^2 = lambda_k^2 + 2*lambda_k*(z-lambda_k) + (z-lambda_k)^2
# [nG_p]^2 = 1/(z-lambda_k)^2 + 2u_k/(z-lambda_k) + (u_k^2 + ...) + ...
# z^2 * [nG_p]^2 = lambda_k^2/(z-lambda_k)^2 + (2*lambda_k + 2*lambda_k^2*u_k)/(z-lambda_k) + ...
# Residue = 2*lambda_k + 2*lambda_k^2*u_k

# Similarly at nu_k: Residue = 2*nu_k + 2*nu_k^2*h_k

# z^2 * ([nG_p]^2 - [nG_r]^2) ~ z^2 * 2n(sum(lambda)-sum(nu))/z^3 = O(1/z)
# This does NOT vanish at infinity! The residue at infinity is:
# Res_inf = -2n*(sum lambda - sum nu)
# = -2n*(sum lambda - (sum lambda + sum mu)) = 2n * sum mu

# So: sum (2*lambda_k + 2*lambda_k^2*u_k) - sum (2*nu_k + 2*nu_k^2*h_k) + 2n*sum(mu) = 0

# Since sum lambda_k^2 * u_k = sum_k sum_{j!=k} lambda_k^2/(lambda_k-lambda_j) and
# sum nu_k^2 * h_k = similar:

# Actually, let me compute sum x_k^2 * H(x_k).
print("Computing sum x_k^2 * H(x_k):")
for pts in [np.array([1., 2., 3.]), np.array([-1., 0., 2., 5.])]:
    H = H_values(pts)
    val = np.dot(pts**2, H)
    n = len(pts)
    sum_pts = np.sum(pts)
    # sum x_k^2 * H(x_k) = sum_{k!=j} x_k^2/(x_k-x_j)
    # = sum_{k<j} [x_k^2/(x_k-x_j) + x_j^2/(x_j-x_k)]
    # = sum_{k<j} [x_k^2/(x_k-x_j) - x_j^2/(x_k-x_j)]
    # = sum_{k<j} (x_k^2 - x_j^2)/(x_k - x_j)
    # = sum_{k<j} (x_k + x_j)
    # = (n-1)*sum x_k  (each x_k appears in (n-1) pairs)
    expected = (n-1) * sum_pts
    print(f"  Points {pts}: sum x^2*H = {val:.6f}, (n-1)*sum(x) = {expected:.6f}")

# YES: sum x_k^2 * H(x_k) = (n-1) * sum x_k.
# So this identity is also "universal" and doesn't use the subordination structure.

# ================================================================
# BETTER APPROACH: Use the SUBORDINATION identity G_r(z) = G_p(omega_1(z))
# ================================================================
print("\n\n" + "="*70)
print("USING THE SUBORDINATION IDENTITY DIRECTLY")
print("="*70)

# G_r(z) = G_p(omega_1(z))
# Differentiate: G_r'(z) = G_p'(omega_1(z)) * omega_1'(z)
# At z = nu_k: both sides have poles. Let's be careful.
#
# Actually, let's work with H(z) = p''(z)/(2p'(z)).
# For z NOT a root: this is well-defined.
# Relationship to G: nG_p(z) = p'(z)/p(z) = sum 1/(z-lambda_k).
# H_p(z) = p''(z)/(2p'(z)) = ... this is the LOG DERIVATIVE of p'.
#
# Near z = lambda_k: p'(lambda_k) != 0 (simple roots), p(lambda_k) = 0.
# p'(z)/p(z) has a simple pole at lambda_k with residue 1.
# nG_p(z) = 1/(z-lambda_k) + H_p(lambda_k) + O(z-lambda_k)
# where H_p(lambda_k) = sum_{j!=k} 1/(lambda_k - lambda_j).
#
# Now: G_r(z) = G_p(omega_1(z)).
# Let's compute something useful from this identity.
#
# Consider the function:
# F(z) = nG_r(z) * nG_r(z) - sum_k 1/(z-nu_k)^2
# = [nG_r(z)]^2 - sum_k 1/(z-nu_k)^2
#
# Near z = nu_l: nG_r(z) = 1/(z-nu_l) + h_l + ...
# [nG_r(z)]^2 = 1/(z-nu_l)^2 + 2*h_l/(z-nu_l) + h_l^2 + 2*H'_l + ...
# - 1/(z-nu_l)^2 cancels the leading pole.
# So F(z) near nu_l = 2*h_l/(z-nu_l) + (smooth terms) + sum_{k!=l} ...
# F has simple poles at nu_l with residue 2*h_l.
#
# Similarly: [nG_p(omega_1(z))]^2 = [nG_r(z)]^2 (same function!)
# Near z = nu_l: nG_p(omega_1(z)) = 1/(z-nu_l) + h_l + ... (we computed this)
# So [nG_p(omega_1)]^2 = 1/(z-nu_l)^2 + 2h_l/(z-nu_l) + ...
#
# IDENTITY: [nG_p(omega_1(z))]^2 = [sum_k 1/(omega_1(z)-lambda_k)]^2
# = sum_k 1/(omega_1(z)-lambda_k)^2 + 2*sum_{k<j} 1/((omega_1(z)-lambda_k)(omega_1(z)-lambda_j))
#
# Hmm, this is getting complex. Let me try the SIMPLEST useful identity.
#
# Consider: d/dz[nG_r(z)] = -sum_k 1/(z-nu_k)^2 = -||nG_r at poles||^2 ... no.
#
# Actually, the KEY observation may be that the subordination identity relates
# the CAUCHY TRANSFORMS but what we need involves the H-TRANSFORMS (which
# are the regularized values at the poles).

# ================================================================
# APPROACH: TRACE FORMULA
# ================================================================
print("\n" + "="*70)
print("TRACE FORMULA APPROACH")
print("="*70)

# Consider the n x n matrices:
# D_nu = diag(nu_1, ..., nu_n)
# D_lambda = diag(lambda_1, ..., lambda_n)
#
# C_{ij}^r = 1/(nu_i - nu_j) for i != j, 0 for i = j  (Cauchy-like matrix for r)
# C_{ij}^p = 1/(lambda_i - lambda_j) for i != j, 0 for i = j
#
# h = C^r * 1  (row sums of C^r)
# u = C^p * 1  (row sums of C^p)
# alpha = u - h
#
# <h, alpha> = 1^T diag(h) * alpha = 1^T diag(h) * (C^p - C^r) * 1
#
# Note: C^r is antisymmetric. C^p is antisymmetric.
# diag(h) = diag(C^r * 1).
#
# For the anti-symmetric matrix C^r: C^r = A^r where A^r_{ij} = 1/(nu_i-nu_j) for i!=j.
#
# <h, alpha> = sum_k h_k * (u_k - h_k)
# = sum_k (sum_j A^r_{kj}) * (sum_l A^p_{kl} - sum_j A^r_{kj})
# = sum_k (sum_j A^r_{kj}) * (sum_l A^p_{kl}) - sum_k (sum_j A^r_{kj})^2

# The first term: sum_k sum_j sum_l A^r_{kj} * A^p_{kl}
# = sum_k sum_{j!=k} sum_{l!=k} 1/(nu_k-nu_j) * 1/(lambda_k-lambda_l)

# This is a triple sum that I can't simplify further without more structure.

# ================================================================
# KEY INSIGHT ATTEMPT: SCHUR CONVEXITY
# ================================================================
print("\n" + "="*70)
print("SCHUR CONVEXITY ATTEMPT")
print("="*70)

# Recall: <h,alpha> = sum_{k<j} (alpha_j - alpha_k)/(nu_j - nu_k)
#
# alpha_k = H_p(lambda_k) - H_r(nu_k)
#
# H_p(lambda_k) = sum_{l!=k} 1/(lambda_k - lambda_l)
# H_r(nu_k) = sum_{j!=k} 1/(nu_k - nu_j)
#
# alpha_k = sum_{l!=k} [1/(lambda_k-lambda_l) - 1/(nu_k-nu_l)]
#
# = sum_{l!=k} [(nu_k-nu_l) - (lambda_k-lambda_l)] / [(lambda_k-lambda_l)(nu_k-nu_l)]
#
# = sum_{l!=k} [(nu_k-lambda_k) - (nu_l-lambda_l)] / [(lambda_k-lambda_l)(nu_k-nu_l)]
#
# Define delta_k = nu_k - lambda_k = nu_k - omega_1(nu_k).
# Then alpha_k = sum_{l!=k} (delta_k - delta_l) / [(lambda_k-lambda_l)(nu_k-nu_l)]

# This expresses alpha in terms of the "shifts" delta_k.

# <h,alpha> = sum_k h_k * alpha_k
# = sum_k [sum_{j!=k} 1/(nu_k-nu_j)] * [sum_{l!=k} (delta_k-delta_l)/((lambda_k-lambda_l)(nu_k-nu_l))]

# Verify this formula:
delta = roots_r - roots_p  # nu - lambda = -(lambda - nu) = -phi(nu)
alpha_check = np.zeros(n)
for k in range(n):
    for l in range(n):
        if l != k:
            alpha_check[k] += (delta[k] - delta[l]) / ((roots_p[k]-roots_p[l]) * (roots_r[k]-roots_r[l]))

alpha_direct = H_values(roots_p) - H_values(roots_r)
print(f"alpha from delta formula: {alpha_check}")
print(f"alpha from H_p - H_r:    {alpha_direct}")
print(f"Match: {np.allclose(alpha_check, alpha_direct, atol=1e-8)}")

# delta_k = nu_k - lambda_k. For the MSS convolution:
# The roots of r are "shifted" versions of the roots of p (and q).
# What is the sign pattern of delta?

print(f"\ndelta = nu - lambda = {delta}")
print(f"(positive means r-root is to the right of p-root)")

# Generally, delta is not all-positive or all-negative.
# The root spreading from the convolution means the OUTER roots of r
# are more extreme than those of p, while INNER roots may shift either way.

# ================================================================
# CAN WE PROVE <h,alpha> >= 0 USING A MATRIX PSD ARGUMENT?
# ================================================================
print("\n" + "="*70)
print("PSD MATRIX ARGUMENT ATTEMPT")
print("="*70)

# Define the n x n Loewner matrix:
# L_{kl} = (omega_1(nu_k) - omega_1(nu_l)) / (nu_k - nu_l) = (lambda_k - lambda_l)/(nu_k - nu_l) for k != l
# L_{kk} = omega_1'(nu_k) = 1

# omega_1 maps C^+ to C^+ (Herglotz). A function f is operator monotone iff
# the Loewner matrix [f(x_i)-f(x_j)]/(x_i-x_j) is PSD.
# For Herglotz functions: they are NOT necessarily operator monotone on all of R.
# operator monotone on R implies Herglotz, but not vice versa.

# Let's check if the Loewner matrix of omega_1 is PSD:
L = np.zeros((n, n))
for k in range(n):
    for l in range(n):
        if k == l:
            L[k, l] = 1.0  # omega_1'(nu_k) = 1
        else:
            L[k, l] = (roots_p[k] - roots_p[l]) / (roots_r[k] - roots_r[l])

print(f"\nLoewner matrix of omega_1:")
print(L)
eigenvalues = np.linalg.eigvalsh(L)
print(f"Eigenvalues: {eigenvalues}")
print(f"PSD: {np.all(eigenvalues > -1e-10)}")

# Let's test over many examples:
np.random.seed(42)
psd_count = 0
not_psd_count = 0
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
        if np.any(np.abs(np.imag(raw)) > 0.01):
            continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01):
            continue

        L = np.zeros((n, n))
        for k in range(n):
            for l in range(n):
                if k == l:
                    L[k, l] = 1.0
                else:
                    L[k, l] = (roots_p[k] - roots_p[l]) / (roots_r[k] - roots_r[l])

        eigs = np.linalg.eigvalsh(L)
        total += 1
        if np.min(eigs) > -1e-8:
            psd_count += 1
        else:
            not_psd_count += 1
    except:
        pass

print(f"\nLoewner matrix PSD: {psd_count}/{total}")
print(f"Loewner matrix NOT PSD: {not_psd_count}/{total}")

# ================================================================
# NEW KEY IDEA: Use the Loewner matrix and the vector h to express <h,alpha>
# ================================================================
print("\n" + "="*70)
print("LOEWNER MATRIX AND <h,alpha>")
print("="*70)

# L_{kl} = (lambda_k - lambda_l)/(nu_k - nu_l) for k!=l, L_{kk} = 1.
#
# Note: sum_{l!=k} L_{kl} * 1/(nu_k - nu_l) = sum_{l!=k} (lambda_k-lambda_l)/(nu_k-nu_l)^2
# Hmm, not directly h or u.
#
# Let e = (1,...,1)^T.
# h = C^r * e where C^r_{kl} = 1/(nu_k-nu_l) for k!=l, 0 on diagonal.
#
# L * e: (L*e)_k = sum_l L_{kl} = 1 + sum_{l!=k} (lambda_k-lambda_l)/(nu_k-nu_l)
#
# Hmm, sum_{l!=k} (lambda_k-lambda_l)/(nu_k-nu_l)
# = lambda_k * sum_{l!=k} 1/(nu_k-nu_l) - sum_{l!=k} lambda_l/(nu_k-nu_l)
# = lambda_k * h_k - sum_{l!=k} lambda_l/(nu_k-nu_l)
#
# This is not obviously related to <h,alpha>.
#
# Actually, let me think about this differently.
# Define the diagonal matrix D = diag(h_1,...,h_n).
# Then <h,alpha> = h^T * alpha = e^T * D * alpha = e^T * D * (u - h)
#
# u_k = H_p(lambda_k) = sum_{l!=k} 1/(lambda_k - lambda_l)
#
# Consider the matrix M_{kl} = 1/(lambda_k - nu_l) (a Cauchy matrix).
# Then (M * e)_k = sum_l 1/(lambda_k - nu_l) = nG_r(lambda_k) (for k != index of nu)
# Wait, this requires lambda_k != nu_l for all l.
#
# Actually, I want to relate <h,u> = sum_k h_k * u_k to some matrix expression.
#
# <h,u> = sum_k [sum_{j!=k} 1/(nu_k-nu_j)] * [sum_{l!=k} 1/(lambda_k-lambda_l)]
#
# = sum_k H_r(nu_k) * H_p(lambda_k)
#
# This is a TRACE of the product diag(H_r) * diag(H_p) acting on the identity permutation.

# Since the mapping is nu_k -> lambda_k = omega_1(nu_k), this is
# sum_k H_r(nu_k) * H_p(omega_1(nu_k)).

# KEY IDEA: Consider the function F(t) = sum_k H_r(nu_k) * H_t(t_k)
# where t_k = (1-t)*nu_k + t*lambda_k (linear interpolation).
# F(0) = ||h||^2 = Phi_r
# F(1) = <h,u> = Phi_r + <h,alpha>
# We need F(1) >= F(0), i.e., F is increasing.
# F'(t) = sum_k h_k * d/dt[H_t(t_k)]
# where t_k(t) = nu_k + t*(lambda_k - nu_k) and H_t uses the points t_1,...,t_n.

# F'(t) = sum_k h_k * [d/dt sum_{l!=k} 1/(t_k - t_l)]
# = sum_k h_k * sum_{l!=k} -(t_k' - t_l')/(t_k - t_l)^2
# where t_k' = lambda_k - nu_k = -delta_k.
# = -sum_k h_k * sum_{l!=k} (delta_l - delta_k)/(t_k - t_l)^2 ... wait signs.
# t_k = nu_k + t*(lambda_k - nu_k) = nu_k - t*delta_k
# t_k' = -delta_k  (derivative w.r.t. t? No: t_k = nu_k + t*(lambda_k-nu_k), so dt_k/dt = lambda_k-nu_k = -delta_k)
# WAIT: delta = nu - lambda, so lambda_k - nu_k = -delta_k.
# t_k = nu_k + t*(-delta_k) = nu_k - t*delta_k

# d/dt H_t(t_k) = d/dt [sum_{l!=k} 1/(t_k - t_l)]
# = sum_{l!=k} -(t_k' - t_l')/(t_k - t_l)^2
# = sum_{l!=k} -(-delta_k + delta_l)/(t_k - t_l)^2
# = sum_{l!=k} (delta_k - delta_l)/(t_k - t_l)^2

# F'(t) = sum_k h_k(t) * sum_{l!=k} (delta_k - delta_l)/(t_k - t_l)^2

# where h_k(t) = H_t(t_k) depends on t.

# At t = 0: t_k = nu_k, h_k(0) = H_r(nu_k) = h_k, and
# F'(0) = sum_k h_k * sum_{l!=k} (delta_k - delta_l)/(nu_k - nu_l)^2
# = sum_{k!=l} h_k * (delta_k - delta_l)/(nu_k - nu_l)^2
# = sum_{k<l} [h_k*(delta_k-delta_l) + h_l*(delta_l-delta_k)]/(nu_k-nu_l)^2
# = sum_{k<l} (h_k - h_l)*(delta_k - delta_l)/(nu_k - nu_l)^2

# F'(0) = sum_{k<l} (h_k - h_l)(delta_k - delta_l) / (nu_k - nu_l)^2

# For F'(0) >= 0, we'd need this sum to be non-negative.
# This involves the CORRELATION between the differences of h and delta.
# If h and delta are "similarly ordered" (both increasing or both decreasing),
# these terms would be positive.

# Let me check numerically:
roots_p = np.array([-2., 0., 3.])
roots_q = np.array([-1., 1., 2.])
roots_r, _ = boxplus_mss(roots_p, roots_q)
roots_r = np.sort(np.real(roots_r))
n = len(roots_p)

h = H_values(roots_r)
delta = roots_r - roots_p

print(f"\nh = {h}")
print(f"delta = nu - lambda = {delta}")

# F'(0):
Fp0 = 0
for k in range(n):
    for l in range(k+1, n):
        term = (h[k]-h[l])*(delta[k]-delta[l]) / (roots_r[k]-roots_r[l])**2
        print(f"  k={k},l={l}: (h_k-h_l)*(delta_k-delta_l)/(nu_k-nu_l)^2 = "
              f"({h[k]:.4f}-{h[l]:.4f})*({delta[k]:.4f}-{delta[l]:.4f})/({roots_r[k]:.4f}-{roots_r[l]:.4f})^2 = {term:.6f}")
        Fp0 += term

print(f"\nF'(0) = {Fp0:.8f}")
print(f"<h,alpha> = {np.dot(h, H_values(roots_p) - h):.8f}")
print(f"Note: F'(0) != <h,alpha> in general (they are different quantities)")

# F'(0) is the derivative at t=0, while <h,alpha> = F(1) - F(0) is the total change.

# ================================================================
# FINAL PROMISING APPROACH: DIRECT COMPUTATION FOR GENERAL n
# ================================================================
print("\n\n" + "="*70)
print("DIRECT COMPUTATION: SIGN ANALYSIS OF <h,alpha>")
print("="*70)

# <h,alpha> = sum_{k<j} (alpha_j - alpha_k)/(nu_j - nu_k)
# where alpha_k = u_k - h_k = H_p(lambda_k) - H_r(nu_k)
#
# alpha_k = sum_{l!=k} [1/(lambda_k-lambda_l) - 1/(nu_k-nu_l)]
#         = sum_{l!=k} (nu_k-nu_l-lambda_k+lambda_l) / [(lambda_k-lambda_l)(nu_k-nu_l)]
#
# For k < j (so nu_k < nu_j and lambda_k < lambda_j):
# alpha_j - alpha_k = sum_{l!=j} 1/(lambda_j-lambda_l) - sum_{l!=j} 1/(nu_j-nu_l)
#                    - sum_{l!=k} 1/(lambda_k-lambda_l) + sum_{l!=k} 1/(nu_k-nu_l)
#
# This involves many cross terms. Let me instead check a SPECIFIC structural hypothesis:
#
# HYPOTHESIS: <h,alpha> can be written as a sum of manifestly non-negative terms,
# each involving the Loewner matrix of omega_1.

# The Loewner matrix L_{kl} = (lambda_k-lambda_l)/(nu_k-nu_l) was found to be
# PSD in ALL tested cases. This means there exists a matrix B such that L = B^T*B.

# If L is PSD with diagonal all 1's (all L_{kk} = 1), then L is a CORRELATION matrix.
# Can we express <h,alpha> in terms of L?

# <h,alpha> = sum_k h_k * alpha_k = <h, u-h>
# u_k = sum_{l!=k} 1/(lambda_k - lambda_l)
# h_k = sum_{l!=k} 1/(nu_k - nu_l)
# alpha_k = u_k - h_k = sum_{l!=k} [1/(lambda_k-lambda_l) - 1/(nu_k-nu_l)]
# = sum_{l!=k} [1/(nu_k-nu_l)] * [(nu_k-nu_l)/(lambda_k-lambda_l) - 1]
# = sum_{l!=k} [1/(nu_k-nu_l)] * [1/L_{kl} - 1]  (using L_{kl} = (lambda-lambda)/(nu-nu))
# = sum_{l!=k} [1/(nu_k-nu_l)] * [(1 - L_{kl})/L_{kl}]

# So alpha_k = sum_{l!=k} (1-L_{kl}) / [(nu_k-nu_l)*L_{kl}]

# <h,alpha> = sum_k h_k * alpha_k
# = sum_k [sum_{j!=k} 1/(nu_k-nu_j)] * [sum_{l!=k} (1-L_{kl})/((nu_k-nu_l)*L_{kl})]
# = sum_k sum_{j!=k} sum_{l!=k} (1-L_{kl}) / [(nu_k-nu_j)(nu_k-nu_l)*L_{kl}]

# This is a complicated triple sum. For k=j:
# inner sum is sum_{l!=k} (1-L_{kl})/[(nu_k-nu_l)^2 * L_{kl}]
# For k != j != l: cross terms.

# Since L is PSD and has diagonal 1, we know 0 < L_{kl} <= 1 (for the correlation).
# Wait, actually L_{kl} can be > 1 (if the lambda-gaps are larger than nu-gaps).
# And L_{kl} can be < 0 (if the gaps have opposite orientations... but they don't
# since both are sorted). Wait: lambda and nu are BOTH sorted, so:
# For k < l: lambda_k < lambda_l and nu_k < nu_l.
# lambda_k - lambda_l < 0 and nu_k - nu_l < 0.
# So L_{kl} = (negative)/(negative) > 0.
# For k > l: both numerator and denominator are positive. So L_{kl} > 0 always.
# CONCLUSION: All entries of L are positive.

print("\nChecking: is L entry-wise positive?")
np.random.seed(42)
entry_pos_count = 0
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
        if np.any(np.abs(np.imag(raw)) > 0.01):
            continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01):
            continue

        L = np.zeros((n, n))
        for k in range(n):
            for l in range(n):
                if k == l:
                    L[k, l] = 1.0
                else:
                    L[k, l] = (roots_p[k]-roots_p[l]) / (roots_r[k]-roots_r[l])

        total += 1
        if np.all(L > -1e-10):
            entry_pos_count += 1
    except:
        pass

print(f"  L entry-wise positive: {entry_pos_count}/{total}")

# OK so L is PSD with all entries positive and diagonal = 1.
# 1 - L_{kl} can be negative (if L_{kl} > 1).

# Let me check the range of L_{kl}:
print("\nRange of off-diagonal L_{kl}:")
np.random.seed(42)
min_L = np.inf
max_L = -np.inf
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
        for k in range(n):
            for l in range(n):
                if k != l:
                    val = (roots_p[k]-roots_p[l]) / (roots_r[k]-roots_r[l])
                    min_L = min(min_L, val)
                    max_L = max(max_L, val)
    except:
        pass

print(f"  min L_{kl} (off-diagonal) = {min_L:.6f}")
print(f"  max L_{kl} (off-diagonal) = {max_L:.6f}")
print(f"  So L_{kl} in [{min_L:.3f}, {max_L:.3f}] and can be > 1 or < 1")
