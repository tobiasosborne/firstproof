"""
Investigate the INDUCTIVE APPROACH to the Fisher inequality.

Key MSS identity: r'(x) = n * (p^{(1)} boxplus_{n-1} q^{(1)})(x)
where p^{(1)} = (1/n)*p' is the monic polynomial whose roots are the critical points of p.

QUESTION: Is there a relationship between Phi_n(p) and Phi_{n-1}(p^{(1)})?

If yes, and if the inequality holds at degree n-1, we might be able to "lift" it to degree n.

Also investigate: does Phi_n decompose in a nice way related to this derivative structure?
"""

import numpy as np
from math import factorial, comb
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
    return np.sort(np.real(np.roots(c))), c

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


def critical_points(roots):
    """Compute the critical points of p(x) = prod(x-r_i), i.e., roots of p'."""
    n = len(roots)
    poly = np.poly(roots)
    dpoly = np.polyder(poly)
    cps = np.sort(np.real(np.roots(dpoly)))
    return cps


np.random.seed(42)

print("="*70)
print("RELATIONSHIP BETWEEN Phi_n(p) AND Phi_{n-1}(p^(1))")
print("="*70)

for trial in range(20):
    n = np.random.choice([3, 4, 5, 6])
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5

    cps = critical_points(roots_p)

    if len(cps) == n - 1 and np.all(np.diff(cps) > 0.01):
        Ph_n = Phi_n(roots_p)
        Ph_n1 = Phi_n(cps)  # Phi_{n-1}(p^{(1)})

        # Is there a pattern?
        print(f"n={n}: Phi_n(p) = {Ph_n:.6f}, Phi_{n-1}(p^(1)) = {Ph_n1:.6f}, ratio = {Ph_n/Ph_n1:.6f}")


# Check the MSS derivative identity: r' = n * (p^{(1)} boxplus_{n-1} q^{(1)})
print("\n" + "="*70)
print("VERIFYING MSS DERIVATIVE IDENTITY")
print("="*70)

for trial in range(10):
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

        # Critical points of p, q, r
        cp_p = critical_points(roots_p)
        cp_q = critical_points(roots_q)
        cp_r = critical_points(roots_r)

        # Compute p^{(1)} boxplus_{n-1} q^{(1)} using MSS
        if len(cp_p) == n-1 and len(cp_q) == n-1:
            r1_roots, _ = boxplus_mss(cp_p, cp_q)
            r1_roots = np.sort(np.real(r1_roots))

            # r' = n * p^{(1)} boxplus_{n-1} q^{(1)}
            # So critical points of r should equal roots of p^{(1)} boxplus_{n-1} q^{(1)}
            error = np.max(np.abs(cp_r - r1_roots))
            print(f"n={n}: cp_r = {cp_r}")
            print(f"       r1_roots = {r1_roots}")
            print(f"       error = {error:.8f}")
            print()
    except Exception as e:
        print(f"  Error: {e}")


# Now check if the inequality is INDUCTIVE:
# If 1/Phi_{n-1}(p^{(1)} boxplus q^{(1)}) >= 1/Phi_{n-1}(p^{(1)}) + 1/Phi_{n-1}(q^{(1)}),
# does this help prove the degree-n inequality?
print("\n" + "="*70)
print("CHECKING INDUCTIVE STRUCTURE")
print("="*70)

for trial in range(20):
    n = np.random.choice([4, 5, 6])
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

        cp_p = critical_points(roots_p)
        cp_q = critical_points(roots_q)
        cp_r = critical_points(roots_r)

        if len(cp_p) != n-1 or len(cp_q) != n-1 or len(cp_r) != n-1:
            continue
        if np.any(np.diff(cp_p) < 0.01) or np.any(np.diff(cp_q) < 0.01) or np.any(np.diff(cp_r) < 0.01):
            continue

        # Degree n inequality
        Ph_p = Phi_n(roots_p)
        Ph_q = Phi_n(roots_q)
        Ph_r = Phi_n(roots_r)
        ineq_n = 1/Ph_r - 1/Ph_p - 1/Ph_q

        # Degree n-1 inequality (on critical points)
        Ph_cp = Phi_n(cp_p)
        Ph_cq = Phi_n(cp_q)
        Ph_cr = Phi_n(cp_r)
        ineq_n1 = 1/Ph_cr - 1/Ph_cp - 1/Ph_cq

        print(f"n={n}: ineq_n = {ineq_n:.8f}, ineq_{n-1} = {ineq_n1:.8f}")
        print(f"  Phi_n:   p={Ph_p:.4f}, q={Ph_q:.4f}, r={Ph_r:.4f}")
        print(f"  Phi_{n-1}: p'={Ph_cp:.4f}, q'={Ph_cq:.4f}, r'={Ph_cr:.4f}")
        print(f"  Phi_n(p)/Phi_{n-1}(p') = {Ph_p/Ph_cp:.6f}")
        print(f"  Phi_n(r)/Phi_{n-1}(r') = {Ph_r/Ph_cr:.6f}")
        print()
    except:
        pass


# ================================================================
# CRITICAL: Is there a NICE formula relating Phi_n(p) and Phi_{n-1}(p')?
# ================================================================
print("\n" + "="*70)
print("FORMULA FOR Phi_n(p) IN TERMS OF Phi_{n-1}(p') AND ROOT DATA")
print("="*70)

# Let lambda_1 < ... < lambda_n be roots of p.
# Let xi_1 < ... < xi_{n-1} be roots of p' (critical points, interlacing with lambda).
# H_p(lambda_i) = sum_{j != i} 1/(lambda_i - lambda_j)
# H_{p'}(xi_j) = sum_{k != j} 1/(xi_j - xi_k)
# Phi_n(p) = sum_i H_p(lambda_i)^2
# Phi_{n-1}(p') = sum_j H_{p'}(xi_j)^2

# Note: p''(lambda_i) = 2*p'(lambda_i)*H_p(lambda_i)
# And p'(lambda_i) = prod_{j != i} (lambda_i - lambda_j)

# From the relation p(x) = (x - lambda_k) * q_k(x) where q_k(x) = prod_{j!=k} (x - lambda_j),
# p'(x) = q_k(x) + (x-lambda_k)*q_k'(x)
# p'(lambda_k) = q_k(lambda_k) = prod_{j!=k} (lambda_k - lambda_j)
# H_p(lambda_k) = q_k'(lambda_k)/q_k(lambda_k) = sum_{j!=k} 1/(lambda_k - lambda_j)

# The relationship between Phi_n(p) and Phi_{n-1}(p') involves the interlacing structure.
# Interlacing: lambda_1 < xi_1 < lambda_2 < xi_2 < ... < xi_{n-1} < lambda_n

# A useful identity: sum_i H_p(lambda_i) = 0 (because sum_i sum_{j!=i} 1/(lambda_i-lambda_j) = 0
# by antisymmetry).

# Similarly sum_j H_{p'}(xi_j) = 0.

# Can we relate sum_i H_p(lambda_i)^2 to sum_j H_{p'}(xi_j)^2?

# One approach: use the residue calculus. Define F(z) = sum_i 1/(z - lambda_i).
# Then n*G_p(z) = F(z). And H_p(lambda_i) = Res[F(z)/(z-lambda_i)] evaluated at z=lambda_i.
# Actually H_p(lambda_i) = lim_{z->lambda_i} (z-lambda_i)*F(z) - 1 ... no.
# H_p(lambda_i) = lim_{z->lambda_i} [F(z) - 1/(z-lambda_i)] = sum_{j!=i} 1/(lambda_i - lambda_j).

# Integral representation:
# int_R (n*G_p(x))^2 dx / (2*pi*i) around contour encircling all roots
# = sum_i Res[(n*G_p)^2, lambda_i]
# The residue at lambda_i is:
# Res = lim_{z->lambda_i} (z-lambda_i) * (sum 1/(z-lambda_j))^2
# = lim (z-lambda_i) * [1/(z-lambda_i)^2 + 2*H_p(lambda_i)/(z-lambda_i) + ...]
# = 2*H_p(lambda_i) + 1/(... ) ... hmm, this doesn't simplify easily.

# Actually: n*G_p(z) = F(z) = sum 1/(z-lambda_j)
# F(z)^2 = sum_j 1/(z-lambda_j)^2 + 2*sum_{j<k} 1/((z-lambda_j)(z-lambda_k))
# Res[F^2, lambda_i] = d/dz[sum_{j!=i} 1/(z-lambda_j)]|_{z=lambda_i} + ...
# Wait, let me be more careful.
# F(z) = 1/(z-lambda_i) + sum_{j!=i} 1/(z-lambda_j)
# Near z = lambda_i: F(z) = 1/(z-lambda_i) + H_p(lambda_i) + O(z-lambda_i)
# F(z)^2 = 1/(z-lambda_i)^2 + 2*H_p(lambda_i)/(z-lambda_i) + [H_p(lambda_i)^2 + ...] + ...
# Res[F^2, lambda_i] = 2*H_p(lambda_i)

# So sum_i Res[F^2, lambda_i] = 2*sum_i H_p(lambda_i) = 0.

# And the integral of F^2 around a large contour encircling all roots:
# For |z| large: F(z) = n/z + O(1/z^2), so F^2 ~ n^2/z^2 + O(1/z^3).
# int F^2 dz/(2*pi*i) around large circle = 0 (residue at infinity is 0 for 1/z^2 term).

# This confirms sum H_p(lambda_i) = 0 but doesn't give Phi_n.

# For Phi_n, consider F(z)^3:
# Near lambda_i: F^3 = 1/(z-lambda_i)^3 + 3H/(z-lambda_i)^2 + (3H^2 + ...)/(z-lambda_i) + ...
# Res[F^3, lambda_i] = 3*H_p(lambda_i)^2 + d^2/dz^2 term...
# This is getting messy.

# Let me just check numerically: is there a constant c_n such that
# Phi_n(p) = c_n * Phi_{n-1}(p') + correction?

for trial in range(10):
    n = 5
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5

    cp = critical_points(roots_p)
    if len(cp) == n-1 and np.all(np.diff(cp) > 0.01):
        Ph_n = Phi_n(roots_p)
        Ph_n1 = Phi_n(cp)
        print(f"  Phi_5(p)={Ph_n:.6f}, Phi_4(p')={Ph_n1:.6f}, ratio={Ph_n/Ph_n1:.6f}")


# ================================================================
# ALTERNATIVE: Direct computation via the VARIANCE IDENTITY
# ================================================================
print("\n\n" + "="*70)
print("TESTING: 1/Phi_n as a VARIANCE-LIKE QUANTITY")
print("="*70)

# The reciprocal 1/Phi_n might have a simpler structure.
# 1/Phi_n(p) = 1/sum_i H_p(lambda_i)^2

# In free probability: 1/Phi*(mu) is related to the variance of the "conjugate variable".
# The free Fisher information satisfies: 1/Phi*(mu boxplus nu) >= 1/Phi*(mu) + 1/Phi*(nu).
# This is the free Cramer-Rao inequality.

# In the finite case: 1/Phi_n(p) can be thought of as a measure of "root spread" squared.
# For equally spaced roots lambda_k = k*d, H_p(lambda_k) = sum_{j!=k} 1/((k-j)*d) = H_k/d.
# Phi_n = (1/d^2) * sum H_k^2. So 1/Phi_n = d^2 / sum H_k^2.
# This is proportional to d^2 (gap squared), which IS additive for the free convolution!

# Let's test: for equally spaced roots, does 1/Phi_n = c * d^2?
print("\nEqually spaced roots test:")
for n in [3, 4, 5, 6]:
    for d in [1.0, 2.0, 0.5]:
        roots = np.array([i * d for i in range(n)])
        Ph = Phi_n(roots)
        print(f"  n={n}, d={d}: Phi={Ph:.6f}, 1/Phi={1/Ph:.6f}, d^2={d**2:.6f}, (1/Phi)/d^2={1/(Ph*d**2):.6f}")

# So 1/Phi = C_n * d^2 where C_n depends on n but not on d.
# C_n = 1/(Phi * d^2) = 1/(sum_k H_k^2) where H_k is the H-transform for unit spacing.

# For the MSS convolution of equally spaced roots with spacings d_p and d_q:
# What is the spacing of r = p boxplus_n q?
print("\n\nEqually spaced convolution test:")
for n in [3, 4, 5]:
    for dp, dq in [(1.0, 2.0), (1.0, 1.0), (0.5, 3.0)]:
        roots_p = np.array([i * dp for i in range(n)])
        roots_q = np.array([i * dq for i in range(n)])
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        roots_r = np.sort(np.real(roots_r))

        # Check if r is equally spaced
        gaps_r = np.diff(roots_r)
        is_equal = np.max(gaps_r) - np.min(gaps_r) < 0.001 * np.mean(gaps_r)
        dr = np.mean(gaps_r)

        Ph_p = Phi_n(roots_p)
        Ph_q = Phi_n(roots_q)
        Ph_r = Phi_n(roots_r)

        ineq = 1/Ph_r - 1/Ph_p - 1/Ph_q

        print(f"  n={n}, dp={dp}, dq={dq}: dr={dr:.4f}, equal_spaced={is_equal}")
        print(f"    Phi_p={Ph_p:.6f}, Phi_q={Ph_q:.6f}, Phi_r={Ph_r:.6f}")
        print(f"    1/Phi_r={1/Ph_r:.6f}, 1/Phi_p+1/Phi_q={1/Ph_p+1/Ph_q:.6f}")
        print(f"    ineq={ineq:.8f}")
        if is_equal:
            print(f"    dr^2={dr**2:.6f}, dp^2+dq^2={dp**2+dq**2:.6f}")
            print(f"    dr^2/(dp^2+dq^2) = {dr**2/(dp**2+dq**2):.6f}")
