#!/usr/bin/env python3
"""
Deeper analysis of the <h,alpha> >= 0 problem.

KEY FINDINGS FROM PREVIOUS SCRIPT:
- A_j = sum_k h_k/(nu_k - p_j)^3 is NOT always >= 0 for arbitrary p between roots
- A_j is ALWAYS < 0 for p outside the root range
- So the pointwise decomposition <h,alpha> = sum_j m_j * A_j cannot be used directly.

NEW APPROACHES:
1. Direct contour integral representation
2. The identity <h,u> = sum_k H_r(nu_k) * H_p(lambda_k) as a bilinear form
3. Using the implicit function F(z,w) = r'(z)*p(w) - p'(w)*r(z) = 0
"""

import numpy as np
from math import factorial
from scipy.optimize import brentq

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

print("=" * 70)
print("APPROACH A: Direct bilinear form analysis")
print("=" * 70)

# <h,u> = sum_k H_r(nu_k) * H_p(lambda_k) where lambda_k = omega_1(nu_k)
#
# Using the identity sum_i h_i/(nu_i - p) = -(1/2)*r''/r evaluated at
# p = lambda_k, we get... but that gives sum over i, not a single term.
#
# Let's think about this differently.
# <h,u> - ||h||^2 = sum_k h_k*(u_k - h_k) = sum_k h_k*alpha_k
#
# Now alpha_k = H_p(lambda_k) - H_r(nu_k)
#
# H_p(lambda_k) = sum_{j!=k} 1/(lambda_k - lambda_j)
# H_r(nu_k) = sum_{j!=k} 1/(nu_k - nu_j)
#
# alpha_k = sum_{j!=k} [1/(lambda_k - lambda_j) - 1/(nu_k - nu_j)]
#         = sum_{j!=k} [(nu_k - nu_j) - (lambda_k - lambda_j)] / [(lambda_k - lambda_j)(nu_k - nu_j)]
#         = sum_{j!=k} [(nu_k - lambda_k) - (nu_j - lambda_j)] / [(lambda_k - lambda_j)(nu_k - nu_j)]
#
# Let delta_k = nu_k - lambda_k = nu_k - omega_1(nu_k).
# Then: alpha_k = sum_{j!=k} (delta_k - delta_j) / [(lambda_k - lambda_j)(nu_k - nu_j)]
#
# This is interesting! If delta_k = delta for all k (uniform shift), then alpha_k = 0.
# The non-uniformity of delta drives alpha.

print("\nComputing delta_k = nu_k - lambda_k:")

for n in [3, 4]:
    np.random.seed(100 + n)
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

        delta = r_roots - p_roots  # nu_k - lambda_k
        H_r = H_values(r_roots)
        H_p = H_values(p_roots)
        alpha = H_p - H_r

        print(f"\n  n={n}, trial {trial}:")
        print(f"    lambda = {p_roots}")
        print(f"    nu     = {r_roots}")
        print(f"    delta  = {delta}")
        print(f"    H_p    = {H_p}")
        print(f"    H_r    = {H_r}")
        print(f"    alpha  = {alpha}")
        print(f"    <h,alpha> = {np.dot(H_r, alpha):.8f}")

        # Check: are deltas monotone?
        delta_diffs = np.diff(delta)
        print(f"    delta differences = {delta_diffs}")
        print(f"    delta monotone increasing: {np.all(delta_diffs > -1e-10)}")
        print(f"    delta monotone decreasing: {np.all(delta_diffs < 1e-10)}")


print("\n" + "=" * 70)
print("APPROACH B: Use the relation [nG_r]' = [nG_p(omega_1)]' = [nG_p]'(omega_1)*omega_1'")
print("to express sum_k h_k*alpha_k in terms of nG_r and nG_p")
print("=" * 70)

# Let's use the identity:
# nG_r(z) = nG_p(omega_1(z))
# Differentiate:
# [nG_r]'(z) = [nG_p]'(omega_1(z)) * omega_1'(z)
# At z = nu_k: [nG_r]'(nu_k) = lim_{z->nu_k} -sum_j 1/(z-nu_j)^2
# This diverges! The pole at nu_k contributes -1/(z-nu_k)^2.
#
# More carefully, near z = nu_k:
# nG_r(z) = 1/(z-nu_k) + h_k + h_k^{(1)}*(z-nu_k) + ...
# where h_k^{(1)} = -S_k + (higher order terms of f_k)
# Actually: f_k(z) = sum_{j!=k} 1/(z-nu_j) so f_k'(z) = -sum_{j!=k} 1/(z-nu_j)^2
# h_k^{(1)} = f_k'(nu_k) = -sum_{j!=k} 1/(nu_k-nu_j)^2 = -S_k
#
# [nG_r]'(z) = -1/(z-nu_k)^2 + f_k'(z) = -1/(z-nu_k)^2 - S_k + ...
#
# omega_1(z) = lambda_k + (z-nu_k) + alpha_k*(z-nu_k)^2 + ...
# omega_1'(z) = 1 + 2*alpha_k*(z-nu_k) + ...
# [nG_p]'(w) = -sum_j 1/(w-lambda_j)^2
# [nG_p]'(omega_1(z)) near z=nu_k:
#   omega_1 ~ lambda_k + (z-nu_k)
#   [nG_p]'(lambda_k + (z-nu_k)) = -1/(z-nu_k)^2 + ... (from the pole at lambda_k)
#   More precisely: [nG_p]'(w) = -1/(w-lambda_k)^2 + sum_{j!=k} -1/(w-lambda_j)^2
#   At w = lambda_k + epsilon: [nG_p]'(w) = -1/epsilon^2 - S_k^{(p)} + ...
#   where S_k^{(p)} = sum_{j!=k} 1/(lambda_k - lambda_j)^2
#
# So [nG_p]'(omega_1(z)) * omega_1'(z)
#   = [-1/(z-nu_k)^2 - S_k^{(p)} + ...] * [1 + 2alpha_k*(z-nu_k) + ...]
#   = -1/(z-nu_k)^2 - 2alpha_k/(z-nu_k) - S_k^{(p)} + ...
#
# But this must equal [nG_r]'(z) = -1/(z-nu_k)^2 - S_k + ...
#
# Matching the 1/(z-nu_k) coefficient (there shouldn't be one in [nG_r]'):
# [nG_r]'(z) = -1/(z-nu_k)^2 + f_k'(nu_k) + f_k''(nu_k)*(z-nu_k) + ...
# There is no 1/(z-nu_k) term! So the 1/(z-nu_k) coefficient is 0.
#
# From the RHS: the 1/(z-nu_k) coefficient is -2*alpha_k
# (from the -1/(z-nu_k)^2 * 2*alpha_k*(z-nu_k) term... wait that gives -2*alpha_k, constant)
# Let me redo: -1/(z-nu_k)^2 * [1 + 2*alpha_k*(z-nu_k)] = -1/(z-nu_k)^2 - 2*alpha_k/(z-nu_k)
# So the 1/(z-nu_k) coefficient on RHS is -2*alpha_k
# On LHS: 0 (since [nG_r]' has a double pole with no simple pole)
# => 0 = -2*alpha_k + ??? No, there must be a contribution from the other terms.
#
# From [nG_p]'(omega_1(z)):
# -1/(omega_1-lambda_k)^2 = -1/(z-nu_k)^2 * 1/(1 + alpha_k*(z-nu_k) + ...)^2
# = -1/(z-nu_k)^2 * [1 - 2*alpha_k*(z-nu_k) + ...] = -1/(z-nu_k)^2 + 2*alpha_k/(z-nu_k) + ...
#
# Hmm I got the wrong sign. Let me redo.
# omega_1(z) = lambda_k + (z-nu_k) + alpha_k*(z-nu_k)^2 + ...
# omega_1(z) - lambda_k = (z-nu_k)[1 + alpha_k*(z-nu_k) + ...]
# 1/(omega_1 - lambda_k) = 1/((z-nu_k)[1+alpha_k*(z-nu_k)+...])
#                         = 1/(z-nu_k) * [1 - alpha_k*(z-nu_k) + ...]
#                         = 1/(z-nu_k) - alpha_k + ...
# 1/(omega_1 - lambda_k)^2 = [1/(z-nu_k) - alpha_k + ...]^2
#                            = 1/(z-nu_k)^2 - 2*alpha_k/(z-nu_k) + alpha_k^2 + ...
#
# So -1/(omega_1-lambda_k)^2 = -1/(z-nu_k)^2 + 2*alpha_k/(z-nu_k) - alpha_k^2 + ...
# And the other terms in [nG_p]'(omega_1): -sum_{j!=k} 1/(omega_1-lambda_j)^2
# = -sum_{j!=k} 1/(lambda_k-lambda_j+(z-nu_k)+...)^2
# = -sum_{j!=k} 1/(lambda_k-lambda_j)^2 * 1/[1+(z-nu_k)/(lambda_k-lambda_j)+...]^2
# = -S_k^{(p)} + 2*sum_{j!=k} (z-nu_k)/(lambda_k-lambda_j)^3 + ...
# = -S_k^{(p)} + O(z-nu_k)
#
# Total [nG_p]'(omega_1(z)):
# = [-1/(z-nu_k)^2 + 2*alpha_k/(z-nu_k) - alpha_k^2 + ...] + [-S_k^{(p)} + O(z-nu_k)]
#
# Multiply by omega_1'(z) = 1 + 2*alpha_k*(z-nu_k) + ...:
# = [-1/(z-nu_k)^2 + 2*alpha_k/(z-nu_k) + (-alpha_k^2 - S_k^{(p)}) + ...]
#   * [1 + 2*alpha_k*(z-nu_k) + ...]
# = -1/(z-nu_k)^2 + 2*alpha_k/(z-nu_k) + (-alpha_k^2 - S_k^{(p)})
#   - 2*alpha_k/(z-nu_k) + 4*alpha_k^2 + ...  [from cross terms]
# = -1/(z-nu_k)^2 + (3*alpha_k^2 - S_k^{(p)}) + ...
#
# Compare with [nG_r]'(z) = -1/(z-nu_k)^2 + (-S_k) + O(z-nu_k)
#
# Matching constant terms: -S_k = 3*alpha_k^2 - S_k^{(p)}
# => S_k^{(p)} - S_k = 3*alpha_k^2
# => sum_{j!=k} [1/(lambda_k-lambda_j)^2 - 1/(nu_k-nu_j)^2] = 3*alpha_k^2
#
# This is an identity relating alpha_k to the change in the "S" sums!

print("\nVerifying the identity: S_k^{(p)} - S_k = 3*alpha_k^2")
print("where S_k^{(p)} = sum_{j!=k} 1/(lambda_k-lambda_j)^2")
print("and   S_k = sum_{j!=k} 1/(nu_k-nu_j)^2")

for n in [3, 4, 5]:
    np.random.seed(100 + n)
    all_match = True
    for trial in range(20):
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

        H_p = H_values(p_roots)
        H_r = H_values(r_roots)
        alpha = H_p - H_r

        # Compute S_k and S_k^{(p)}
        for k in range(n):
            Sk = sum(1/(r_roots[k]-r_roots[j])**2 for j in range(n) if j != k)
            Skp = sum(1/(p_roots[k]-p_roots[j])**2 for j in range(n) if j != k)

            lhs = Skp - Sk
            rhs = 3 * alpha[k]**2

            if abs(lhs - rhs) > 1e-4 * (abs(lhs) + abs(rhs) + 1e-8):
                all_match = False
                print(f"  MISMATCH n={n}, trial {trial}, k={k}: "
                      f"S^(p)-S={lhs:.8f}, 3*alpha^2={rhs:.8f}")

    if all_match:
        print(f"  n={n}: Identity S_k^(p) - S_k = 3*alpha_k^2 verified in all trials")
    else:
        print(f"  n={n}: Some mismatches found")


print("\n" + "=" * 70)
print("Now: can we use S_k^{(p)} - S_k = 3*alpha_k^2 to prove <h,alpha> >= 0?")
print("=" * 70)

# We have: S_k^{(p)} = S_k + 3*alpha_k^2
# where S_k = sum_{j!=k} 1/(nu_k-nu_j)^2 and S_k^{(p)} = sum_{j!=k} 1/(lambda_k-lambda_j)^2
#
# Now consider:
# Phi_n(p) = sum_k H_p(lambda_k)^2 = sum_k u_k^2 = ||u||^2
# Phi_n(r) = sum_k H_r(nu_k)^2 = sum_k h_k^2 = ||h||^2
#
# Note: sum_k S_k = sum_k sum_{j!=k} 1/(nu_k-nu_j)^2 = sum_{k!=j} 1/(nu_k-nu_j)^2
# = 2*sum_{k<j} 1/(nu_k-nu_j)^2  (by symmetry (nu_k-nu_j)^2 = (nu_j-nu_k)^2)
# Similarly sum_k S_k^{(p)} = 2*sum_{k<j} 1/(lambda_k-lambda_j)^2
#
# From the identity: sum_k S_k^{(p)} - sum_k S_k = 3*sum_k alpha_k^2 = 3*||alpha||^2
# => 3*||alpha||^2 = sum_k S_k^{(p)} - sum_k S_k  ... (*)
#
# Also: ||u||^2 = ||h+alpha||^2 = ||h||^2 + 2<h,alpha> + ||alpha||^2
# => <h,alpha> = (||u||^2 - ||h||^2 - ||alpha||^2) / 2
# => <h,alpha> = (Phi_p - Phi_r - ||alpha||^2) / 2

# From (*): ||alpha||^2 = (1/3)*(sum_k S_k^{(p)} - sum_k S_k)
# => <h,alpha> = (Phi_p - Phi_r) / 2 - (1/6)*(sum_k S_k^{(p)} - sum_k S_k)

# Question: is there a relationship between Phi and sum S?
# Phi_n(r) = sum_k h_k^2 and sum_k S_k = sum_{k,j:k!=j} 1/(nu_k-nu_j)^2
#
# Note: H_r(nu_k)^2 = [sum_{j!=k} 1/(nu_k-nu_j)]^2 = sum_{j!=k} 1/(nu_k-nu_j)^2 + cross terms
# = S_k + 2*sum_{j<l, j!=k, l!=k} 1/((nu_k-nu_j)(nu_k-nu_l))
#
# So Phi_r = sum_k S_k + 2*sum_k sum_{j<l, j,l!=k} 1/((nu_k-nu_j)(nu_k-nu_l))
# The cross terms are complicated.

print("\nLet's verify: <h,alpha> = (Phi_p - Phi_r - ||alpha||^2)/2")
for n in [3, 4]:
    np.random.seed(100 + n)
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

        H_p = H_values(p_roots)
        H_r = H_values(r_roots)
        alpha = H_p - H_r

        lhs = np.dot(H_r, alpha)
        rhs = (np.dot(H_p, H_p) - np.dot(H_r, H_r) - np.dot(alpha, alpha)) / 2
        print(f"  n={n}: <h,alpha> = {lhs:.8f}, (Phi_p-Phi_r-||a||^2)/2 = {rhs:.8f}, "
              f"match = {abs(lhs-rhs) < 1e-8}")

# This is just the polarization identity, of course it holds.

print("\n" + "=" * 70)
print("APPROACH C: Consider the SECOND chain rule identity more carefully")
print("and try to express <h,alpha> in terms of spectral data")
print("=" * 70)

# Key identity derived: S_k^{(p)} - S_k = 3*alpha_k^2
# This means: alpha_k^2 = (1/3)*(sum_{j!=k} [1/(lambda_k-lambda_j)^2 - 1/(nu_k-nu_j)^2])
# The sign of alpha_k is determined, but alpha_k^2 is controlled by the gap structure.
#
# For <h,alpha> = sum_k h_k * alpha_k, we need to relate the sign of h_k*alpha_k
# to the gap structure.
#
# OBSERVATION from numerics: h_k*alpha_k can be negative for individual k,
# but the sum is always non-negative.
#
# IDEA: Use Cauchy-Schwarz in a WEIGHTED form:
# <h,alpha>^2 <= ||h||^2 * ||alpha||^2   (standard CS)
# This doesn't help since we need LOWER bound on <h,alpha>.
#
# IDEA: Express <h,alpha> as a sum of squares or use AM-GM.
#
# <h,alpha> = sum_k h_k*(u_k - h_k) = sum_k h_k*u_k - ||h||^2
# = <h,u> - ||h||^2
#
# So <h,alpha> >= 0 iff <h,u> >= ||h||^2, i.e., the "projection" of u onto h
# is at least ||h||.
#
# Since u = h + alpha, <h,u> = ||h||^2 + <h,alpha>, this is tautological.
# The key content is: sum_k H_r(nu_k) * H_p(lambda_k) >= sum_k H_r(nu_k)^2

# Let's write this as:
# sum_k H_r(nu_k) * [H_p(lambda_k) - H_r(nu_k)] >= 0
# sum_k H_r(nu_k) * alpha_k >= 0

# Using the NEW identity alpha_k = (omega_1''(nu_k))/2:
# omega_1''(nu_k)/2 = H_p(lambda_k) - H_r(nu_k)
# This holds because omega_1'(nu_k) = 1.

# Now consider: the function omega_1(z) - z is the "displacement" function.
# It's a rational function that maps C+ to C+ (since omega_1 does and z maps C+ to itself,
# but omega_1 - z is the DEVIATION; actually omega_1(z) - z maps C+ to C+ since omega_1 is Herglotz
# and omega_1(z) ~ z + c + ... at infinity, so omega_1(z) - z ~ c + a/z + ... maps C+ to C+).
# Wait: does phi(z) = omega_1(z) - z map C+ to C+?
# omega_1 maps C+ to C+ (Herglotz), but z is the identity which maps C+ to C+.
# Their difference has Im(phi(z)) = Im(omega_1(z)) - Im(z).
# Since omega_1 is Herglotz with omega_1(z) ~ z + c at infinity (c real),
# we have Im(omega_1(z)) >= Im(z) for z in C+ (by the Herglotz property and
# the specific normalization). Actually this is not automatic.
#
# omega_1(z) = z + c + sum_j m_j/(z-p_j) with m_j > 0
# Im(omega_1(z)) = Im(z) + sum_j m_j * Im(1/(z-p_j))
# For z in C+: Im(1/(z-p_j)) = -Im(z)/|z-p_j|^2 < 0
# So Im(omega_1(z)) = Im(z) * (1 - sum_j m_j/|z-p_j|^2)
# This can be < Im(z) or > Im(z)!

# Actually for omega_1(z) = z + c + sum m_j/(z-p_j):
# Im[omega_1(z) - z] = Im[c + sum m_j/(z-p_j)]
# = sum m_j * Im[1/(z-p_j)]
# = -sum m_j * Im(z)/|z-p_j|^2 < 0  for z in C+
# So phi(z) = omega_1(z) - z maps C+ to C- !!!
# i.e., -phi is Herglotz, or equivalently phi maps C+ to C-.

print("\nSo phi(z) = omega_1(z) - z maps C+ to C-")
print("This means -phi is a Nevanlinna/Herglotz function mapping C+ to C+")
print()

# NOW: define Psi(z) = -phi(z) = z - omega_1(z) = -c - sum_j m_j/(z-p_j)
# This is: Psi(z) = -c + sum_j m_j/(p_j - z)
# (careful with signs: -m_j/(z-p_j) = m_j/(p_j-z))
#
# Psi maps C+ to C+ (confirmed).
# At the roots of r: Psi(nu_k) = nu_k - lambda_k = -delta_k (recall delta_k = nu_k - lambda_k...
# Wait: delta_k = nu_k - lambda_k, so Psi(nu_k) = nu_k - omega_1(nu_k) = nu_k - lambda_k = delta_k
# Hmm: omega_1(nu_k) = lambda_k, so Psi(nu_k) = nu_k - lambda_k = delta_k.
#
# Psi'(z) = 1 - omega_1'(z) = sum_j m_j/(z-p_j)^2 > 0 for real z not at poles.
# At nu_k: Psi'(nu_k) = 1 - omega_1'(nu_k) = 1 - 1 = 0
# So the roots of r are critical points of Psi!
#
# Psi''(z) = -omega_1''(z) = -2*sum_j m_j/(z-p_j)^3
# alpha_k = omega_1''(nu_k)/2 = -Psi''(nu_k)/2

# Now: <h,alpha> = sum_k h_k * alpha_k = -sum_k h_k * Psi''(nu_k)/2
# = -(1/2)*sum_k h_k * Psi''(nu_k)

# Using the residue identity from the key identity:
# sum_k h_k * g(nu_k) = sum_k Res_{nu_k}[nG_r(z) * g(z)] (when g is regular at nu_k)
#   (since Res_{nu_k} nG_r(z) = 1, so Res_{nu_k}[nG_r(z)*g(z)] = g(nu_k) for regular g)
# Wait: this gives sum_k g(nu_k), not sum_k h_k*g(nu_k).
#
# For sum_k h_k*g(nu_k), we need:
# h_k = Res_{nu_k} [(nG_r)^2/2] ... no, Res is 2*h_k/2 = h_k, but only if g is constant.
# Actually Res_{nu_k}[(nG_r)^2 * g / 2] involves g and its derivatives.
#
# Let's use a different approach: integration by parts via contour integral.

print("=" * 70)
print("APPROACH D: Integration by parts via contour integral")
print("=" * 70)
print()
print("Consider I = (1/2pi i) oint nG_r(z) * Psi'(z) * nG_r(z) dz")
print("= (1/2pi i) oint [nG_r(z)]^2 * Psi'(z) dz")
print()
print("Near z = nu_k:")
print("  [nG_r]^2 = 1/(z-nu_k)^2 + 2h_k/(z-nu_k) + ...")
print("  Psi'(z) = Psi'(nu_k) + Psi''(nu_k)*(z-nu_k) + ...")
print("          = 0 + Psi''(nu_k)*(z-nu_k) + ...  (since Psi'(nu_k) = 0)")
print()
print("Product = Psi''(nu_k)/(z-nu_k) + 2*h_k*Psi''(nu_k) + ...")
print("Res = Psi''(nu_k) = -2*alpha_k")
print()
print("So sum of residues at nu_k = sum_k (-2*alpha_k) = 0  (since sum alpha_k = 0)")
print()
print("Psi' also has double poles at p_j (poles of omega_1):")
print("  Psi'(z) = sum m_j/(z-p_j)^2, so near p_j:")
print("  Psi'(z) = m_j/(z-p_j)^2 + [regular]")
print("  [nG_r]^2 is regular at p_j, with value [nG_r(p_j)]^2")
print("  Product ~ m_j*[nG_r(p_j)]^2/(z-p_j)^2 + ... ")
print("  This has a double pole but we need to be more careful about the residue.")
print()
print("Res_{p_j}[[nG_r(z)]^2 * Psi'(z)]")
print("= Res_{p_j}[[nG_r(z)]^2 * m_j/(z-p_j)^2 + ...]")
print("= m_j * d/dz|_{p_j} [nG_r(z)]^2  (since Res_{p_j} f(z)/(z-p_j)^2 = f'(p_j))")
print("  PLUS the contribution from the regular part of Psi' at p_j")
print("= m_j * 2*nG_r(p_j)*nG_r'(p_j) + [cross terms from other poles of Psi']")
print()

# This gives us:
# 0 = sum_k (-2*alpha_k) + sum_j [Res_{p_j}] = 0 + sum_j [Res_{p_j}]
# => sum_j Res_{p_j} = 0 (which is just consistency since sum alpha_k = 0)

# Let me try a DIFFERENT integrand that captures <h,alpha> directly.

print("=" * 70)
print("APPROACH E: The key integral oint [nG_r]^2 * [nG_r]' * [omega_1 - z] dz")
print("=" * 70)

# Consider I = (1/2pi i) oint nG_r(z) * nG_r(z) * Psi(z) dz
#            = (1/2pi i) oint [nG_r(z)]^2 * [z - omega_1(z)] dz
# (using Psi(z) = z - omega_1(z))
#
# Near z = nu_k:
# [nG_r]^2 = 1/(z-nu_k)^2 + 2h_k/(z-nu_k) + c_k + ...
# Psi(z) = Psi(nu_k) + Psi'(nu_k)*(z-nu_k) + (1/2)*Psi''(nu_k)*(z-nu_k)^2 + ...
#        = delta_k + 0*(z-nu_k) + (1/2)*Psi''(nu_k)*(z-nu_k)^2 + ...
#        = delta_k - alpha_k*(z-nu_k)^2 + ...
#
# Product = delta_k/(z-nu_k)^2 + 2*h_k*delta_k/(z-nu_k) + ...
#           + [-alpha_k + ...]  (from (1/(z-nu_k)^2)*(Psi''(nu_k)/2)*(z-nu_k)^2)
#
# Res at nu_k = 2*h_k*delta_k + [-alpha_k * coefficient... wait]
# Actually: the product of 1/(z-nu_k)^2 * [delta_k - alpha_k*(z-nu_k)^2 + ...] gives:
# delta_k/(z-nu_k)^2 - alpha_k + O(z-nu_k)
# And 2*h_k/(z-nu_k) * [delta_k + ...] gives:
# 2*h_k*delta_k/(z-nu_k) + O(1)
# And c_k * [delta_k + ...] gives: c_k*delta_k + O(z-nu_k)
#
# Res = 2*h_k*delta_k
#
# So sum Res at nu_k = 2*sum_k h_k*delta_k
#
# And the poles at p_j:
# [nG_r]^2 is regular at p_j; Psi(z) = sum m_l/(p_l-z) -c has simple poles at p_j.
# Psi(z) near p_j: m_j/(p_j-z) + [regular] = -m_j/(z-p_j) + [regular]
# [nG_r(z)]^2 * Psi(z) near p_j has a simple pole with residue:
# -m_j * [nG_r(p_j)]^2
#
# At infinity: [nG_r]^2 ~ 1/z^2, Psi ~ -c + O(1/z)
# Product ~ -c/z^2 + O(1/z^3), integral -> 0
#
# Residue theorem: 0 = 2*sum_k h_k*delta_k - sum_j m_j*[nG_r(p_j)]^2
# => sum_k h_k*delta_k = (1/2)*sum_j m_j*[nG_r(p_j)]^2

# Now: sum_j m_j*[nG_r(p_j)]^2 >= 0 since m_j > 0 and [nG_r]^2 >= 0 (it's real).
# Wait: p_j are real (poles of omega_1 are real since omega_1 is Herglotz with real poles),
# and nG_r(p_j) = sum_k 1/(p_j-nu_k) is real.
# So [nG_r(p_j)]^2 >= 0, and therefore sum_k h_k*delta_k >= 0.
#
# BUT: delta_k = nu_k - lambda_k, and we wanted alpha_k = H_p(lambda_k) - H_r(nu_k).
# These are DIFFERENT quantities! <h,delta> != <h,alpha>.

print("\n*** IMPORTANT FINDING ***")
print("sum_k h_k * delta_k = (1/2)*sum_j m_j*[nG_r(p_j)]^2 >= 0")
print("where delta_k = nu_k - lambda_k")
print("This proves sum_k H_r(nu_k) * (nu_k - lambda_k) >= 0 !!!")
print()
print("But we need sum_k H_r(nu_k) * (H_p(lambda_k) - H_r(nu_k)) >= 0")
print("which is a DIFFERENT statement.")

# Verify this identity:
print("\nVerifying: sum h_k*delta_k = (1/2)*sum m_j*[nG_r(p_j)]^2")
for n in [3, 4]:
    np.random.seed(100+n)
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
        delta = r_roots - p_roots  # nu_k - lambda_k
        lhs = np.dot(H_r, delta)

        # We don't know the poles of omega_1 easily, but we can verify <h,delta> >= 0
        print(f"  n={n}, trial {trial}: <h,delta> = {lhs:.8f} >= 0: {lhs >= -1e-10}")


print("\n" + "=" * 70)
print("APPROACH F: Try [nG_r]^2 * Psi''(z) / 2 directly")
print("but with INTEGRATION BY PARTS")
print("=" * 70)

# We want: <h,alpha> = -(1/2)*sum_k h_k * Psi''(nu_k)
#
# Consider: I = (1/2pi i) oint [nG_r(z)]^2 * Psi(z) dz
# We showed: I = sum of residues at nu_k + sum of residues at p_j
#           = 2*sum h_k*delta_k - sum m_j*[nG_r(p_j)]^2 = 0
#
# Now integrate by parts:
# [nG_r]^2 * Psi = -d/dz{something} * Psi + something * Psi'
# This is not straightforward for a contour integral.
#
# ALTERNATIVELY: Use a different kernel function.
# We need Res_{nu_k} f(z) = h_k * alpha_k.
#
# Recall: alpha_k = -Psi''(nu_k)/2 and h_k = H_r(nu_k).
#
# Consider: f(z) = (1/2) * d/dz{[nG_r(z)]^2 * Psi(z)}
# = (1/2)*{2*nG_r(z)*[nG_r]'(z)*Psi(z) + [nG_r(z)]^2*Psi'(z)}
# The residues of f = (1/2)*d/dz[...] at nu_k:
# Res_{nu_k} f = (1/2)*0 = 0 (derivative of a Laurent series has no residue)
# Wait, that's not right. The residue of f' is 0 only for exact derivatives.
# d/dz[sum a_n (z-z0)^n] = sum n*a_n*(z-z0)^{n-1}
# The coefficient of (z-z_0)^{-1} in f' is: 0*a_0 = 0 if f had a residue a_{-1},
# since d/dz[a_{-1}/(z-z_0)] = -a_{-1}/(z-z_0)^2 (no residue in f').
# Wait: the derivative of g(z) = sum_{n=-N}^{infty} a_n (z-z_0)^n is
# g'(z) = sum_{n=-N}^{infty} n*a_n*(z-z_0)^{n-1}
# The (z-z_0)^{-1} term in g' comes from n = 0: 0*a_0 = 0.
# So yes, Res(g') = 0 for any g.
#
# This means the integral of (d/dz[...]) dz around any closed contour = 0.
# So this approach gives 0 = 0. Not useful.

# Let me try a fundamentally different approach:
# Instead of contour integrals, use the QUADRATIC FORM interpretation.

print("\n" + "=" * 70)
print("APPROACH G: Quadratic form interpretation")
print("=" * 70)
print()
print("<h,alpha> = sum_k H_r(nu_k) * [H_p(lambda_k) - H_r(nu_k)]")
print("= sum_k H_r(nu_k) * H_p(lambda_k) - Phi(r)")
print()
print("Now use the identity from Approach D:")
print("sum_k H_r(nu_k) * (nu_k - lambda_k) = (1/2)*sum_j m_j*[nG_r(p_j)]^2")
print()
print("This gives a nice identity for <h,delta> but we need <h,alpha>.")
print("Can we relate alpha_k to delta_k?")
print()

# Relationship between alpha_k and delta_k:
# alpha_k = H_p(lambda_k) - H_r(nu_k)
#         = sum_{j!=k} 1/(lambda_k-lambda_j) - sum_{j!=k} 1/(nu_k-nu_j)
#
# delta_k = nu_k - lambda_k
#
# For small delta: alpha_k ~ sum_{j!=k} [(delta_j - delta_k)/(nu_k-nu_j)^2]
# (by Taylor expansion of 1/(lambda_k-lambda_j) around nu_k-nu_j)

for n in [3, 4]:
    np.random.seed(100+n)
    print(f"\nn={n}:")
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
        alpha = H_p - H_r
        delta = r_roots - p_roots

        h_alpha = np.dot(H_r, alpha)
        h_delta = np.dot(H_r, delta)

        print(f"  trial {trial}: <h,alpha> = {h_alpha:.6f}, <h,delta> = {h_delta:.6f}, "
              f"ratio = {h_alpha/h_delta if abs(h_delta) > 1e-10 else 'inf'}")

print("\n" + "=" * 70)
print("SUMMARY OF FINDINGS")
print("=" * 70)
