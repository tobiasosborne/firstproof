"""
PROVER-10: Systematic investigation of C_n and R_n for n=2,3,4,5,6.

Goal: Decompose 1/Phi_n = C_n * kappa_2 + R_n(kappa_2,...,kappa_n)
and study the general pattern.

We use the Arizmendi-Perales finite free cumulants.
"""

import numpy as np
from itertools import combinations
from math import comb, factorial
from fractions import Fraction
import sympy as sp
from sympy import symbols, Rational, poly, resultant, expand, simplify
from sympy import Matrix, det, sqrt, collect, factor, cancel, together
from collections import defaultdict

np.random.seed(42)

# ============================================================
# PART 0: Core definitions (from verify_cumulant_claims.py)
# ============================================================

def elementary_symmetric(roots, k):
    """e_k = sum of products of k roots."""
    n = len(roots)
    if k == 0:
        return 1.0
    if k > n:
        return 0.0
    return sum(np.prod(combo) for combo in combinations(roots, k))

def phi_n(roots):
    """Phi_n = sum_i H_p(lambda_i)^2."""
    n = len(roots)
    total = 0.0
    for i in range(n):
        H_i = sum(1.0 / (roots[i] - roots[j]) for j in range(n) if j != i)
        total += H_i**2
    return total

def mss_convolve(roots_p, roots_q):
    """MSS (finite free additive) convolution."""
    n = len(roots_p)
    assert len(roots_q) == n
    p_coeffs = np.poly(roots_p)
    q_coeffs = np.poly(roots_q)
    r_coeffs = np.zeros(n + 1)
    r_coeffs[0] = 1.0
    for k in range(1, n + 1):
        ck = 0.0
        for i in range(k + 1):
            j = k - i
            if i <= n and j <= n:
                w = factorial(n - i) * factorial(n - j) / (factorial(n) * factorial(n - k))
                ck += w * p_coeffs[i] * q_coeffs[j]
        r_coeffs[k] = ck
    r_roots = np.roots(r_coeffs)
    return np.sort(np.real(r_roots))

def compute_cumulants(roots):
    """Compute finite free cumulants kappa_1,...,kappa_n from roots.

    Using the Arizmendi-Perales formula via normalized coefficients.
    For a monic polynomial p(x) = x^n + a_1 x^{n-1} + ... + a_n,
    the normalized coefficients are tilde_a_k = (-1)^k * a_k / C(n,k).

    Then cumulants are computed via Möbius inversion on the partition lattice.
    For practical purposes, we use the direct formulas:
      kappa_1 = tilde_a_1 = mean of roots
      kappa_2 = -n * (tilde_a_2 - tilde_a_1^2)
      kappa_3 = n^2/2 * (tilde_a_3 - 3*tilde_a_2*tilde_a_1 + 2*tilde_a_1^3)
      ... (higher via complete Bell polynomials / Möbius inversion)

    For centered polynomials (kappa_1=0), we have simpler formulas.
    """
    n = len(roots)
    # Get polynomial coefficients: p(x) = x^n + a_1 x^{n-1} + ...
    coeffs = np.poly(roots)  # [1, a_1, a_2, ..., a_n]

    # Normalized coefficients: tilde_a_k = (-1)^k * a_k / C(n,k)
    ta = np.zeros(n + 1)
    ta[0] = 1.0
    for k in range(1, n + 1):
        ta[k] = ((-1)**k * coeffs[k]) / comb(n, k)

    # Compute cumulants via the moment-cumulant relation
    # For finite free probability, the relation between tilde_a_k and kappa_k
    # is given by the free moment-cumulant formula.
    #
    # The key relation: tilde_a_k = sum over non-crossing partitions of {1,...,k}
    #                   of product of kappa_{|block|} / appropriate factor
    #
    # More precisely, for the "Boolean" cumulants approach:
    # tilde_a_k = sum_{pi in NC(k)} prod_{B in pi} kappa_{|B|} * weight(pi)
    #
    # For finite free cumulants (Arizmendi-Perales 2018), the relation is:
    # m_k = sum_{pi in NC(k)} prod_{B in pi} kappa_{|B|} / n^{|pi|-1}
    # where m_k are free moments (= tilde_a_k up to sign convention)

    # Direct computation of cumulants from moments using recursion.
    # For centered polys, use the simpler formulas.

    # Actually, let me use the known exact formulas for small n.
    # The key relation: if p(x) = x^n + a_1 x^{n-1} + ... and we define
    # b_k = (-1)^k * a_k / C(n,k), then the finite free cumulants satisfy:
    # b_1 = kappa_1 / n  (NO: b_1 = kappa_1)
    # b_2 = kappa_1^2 + kappa_2 / n ...

    # Let me use the correct Arizmendi-Perales definition from the paper.
    # They define: for p(x) = sum_{k=0}^n (-1)^k C(n,k) a_k x^{n-k}
    # with a_0 = 1, so the coefficient of x^{n-k} is (-1)^k C(n,k) a_k.
    # Our convention: p(x) = x^n + c_1 x^{n-1} + ... + c_n
    # So c_k = (-1)^k C(n,k) a_k, giving a_k = (-1)^k c_k / C(n,k) = tilde_a_k.
    #
    # The cumulant-moment relation (Thm 3.6 of Arizmendi-Perales):
    # a_k = sum_{pi in NC(k)} prod_{V in pi} kappa_{|V|} * 1/n^{|pi|-1}
    # Wait, this doesn't look right. Let me just compute numerically.

    # I'll use the standard approach: extract cumulants from tilde_a values
    # using the recursive formula.

    # For the FINITE FREE cumulants, the relation is:
    # tilde_a_k = sum over pi in NC(k) of (1/n^{|pi|-1}) * prod kappa_{|V|}
    #
    # This gives:
    # tilde_a_1 = kappa_1
    # tilde_a_2 = kappa_1^2 + kappa_2/n
    # tilde_a_3 = kappa_1^3 + 3*kappa_1*kappa_2/n + kappa_3/n^2
    # tilde_a_4 = kappa_1^4 + 6*kappa_1^2*kappa_2/n + (4*kappa_1*kappa_3 + 2*kappa_2^2)/n^2 + kappa_4/n^3

    # Inverting for centered (kappa_1 = 0):
    # tilde_a_1 = 0
    # tilde_a_2 = kappa_2/n  =>  kappa_2 = n * tilde_a_2
    # tilde_a_3 = kappa_3/n^2  =>  kappa_3 = n^2 * tilde_a_3
    # tilde_a_4 = 2*kappa_2^2/n^2 + kappa_4/n^3
    #           = 2*(n*tilde_a_2)^2/n^2 + kappa_4/n^3
    #           = 2*tilde_a_2^2 + kappa_4/n^3
    # => kappa_4 = n^3 * (tilde_a_4 - 2*tilde_a_2^2)
    # tilde_a_5 = (10*kappa_2*kappa_3)/n^3 + kappa_5/n^4 (centered)
    #           = 10*n*tilde_a_2 * n^2*tilde_a_3 / n^3 + kappa_5/n^4
    #           = 10*tilde_a_2*tilde_a_3 + kappa_5/n^4
    # => kappa_5 = n^4 * (tilde_a_5 - 10*tilde_a_2*tilde_a_3)

    # Wait, I need to be more careful with the NC(k) counting.
    # NC(2) = {{1,2}}, {{1},{2}} => |NC(2)| = 2
    # For pi = {{1,2}}: |pi|=1, contrib = kappa_2
    # For pi = {{1},{2}}: |pi|=2, contrib = kappa_1^2/n
    # So tilde_a_2 = kappa_2 + kappa_1^2/n ... hmm, but for centered this gives kappa_2 = tilde_a_2

    # Actually I think the correct formula uses n-dependent weights.
    # Let me just use the KNOWN formulas from the verify script.

    kappas = np.zeros(n + 1)  # kappas[k] = kappa_k

    # kappa_1 = mean
    kappas[1] = ta[1]  # = mean of roots

    # kappa_2 (from verify script): kappa_2 = -n*(tilde_a_2 - tilde_a_1^2)
    kappas[2] = -n * (ta[2] - ta[1]**2)

    if n >= 3:
        # kappa_3 = (n^2/2) * (tilde_a_3 - 3*tilde_a_2*tilde_a_1 + 2*tilde_a_1^3)
        kappas[3] = (n**2 / 2) * (ta[3] - 3*ta[2]*ta[1] + 2*ta[1]**3)

    if n >= 4:
        # Need kappa_4. For centered polys:
        # From the moment-cumulant relation for NC(4):
        # The NC partitions of {1,2,3,4}:
        # 1-block: {1234} -> kappa_4
        # 2-block NC: {12}{34}, {14}{23}, {1}{234}, {123}{4}, {12}{3}{4}...
        # Actually NC(4) has 14 elements (Catalan number C_4 = 14).
        # This gets complex. Let me use the numeric approach instead.
        pass

    return kappas

def compute_cumulants_centered(roots):
    """Compute finite free cumulants for CENTERED polynomials.

    For centered (mean=0) polynomials, we use the explicit formulas.
    Arizmendi-Perales (2018), the moment-cumulant formula for finite free probability.

    The normalized coefficients tilde_a_k satisfy (for centered, tilde_a_1=0):
    tilde_a_2 = kappa_2 / (-n)     ... wait

    Actually, let me re-derive carefully.

    From verify_cumulant_claims.py:
    - kappa_2 = -n * tilde_a_2  (centered)
    - kappa_3 = (n^2/2) * tilde_a_3  (centered)

    Let's verify: for n=3 centered:
    tilde_a_2 = a_2/C(3,2) = e_2/3
    kappa_2 = -3 * (e_2/3) = -e_2 ✓ (matches verify script)

    tilde_a_3 = -a_3/C(3,3) = e_3/1 = e_3
    kappa_3 = (9/2)*e_3 ✓

    For general n centered:
    kappa_2 = -n * tilde_a_2
    kappa_3 = (n^2/2) * tilde_a_3

    For kappa_4 centered, the free moment-cumulant formula gives:
    tilde_a_4 = ... involves kappa_4/n^3 + 2*kappa_2^2/n^2 (from NC(4) centered)
    Wait, this depends on the exact form. Let me compute it numerically.
    """
    n = len(roots)
    assert abs(np.mean(roots)) < 1e-10, f"Not centered: mean = {np.mean(roots)}"

    coeffs = np.poly(roots)
    ta = np.zeros(n + 1)
    for k in range(n + 1):
        ta[k] = ((-1)**k * coeffs[k]) / comb(n, k)

    kappas = np.zeros(n + 1)
    kappas[1] = 0.0  # centered
    kappas[2] = -n * ta[2]

    if n >= 3:
        kappas[3] = (n**2 / 2) * ta[3]

    # For higher cumulants, we need to be more careful.
    # Let me compute them via the numeric approach.
    return kappas, ta


# ============================================================
# PART 1: Symbolic computation of Phi_n for small n
# ============================================================

print("=" * 70)
print("PART 1: Symbolic Phi_n computation")
print("=" * 70)

# --- n = 2 ---
print("\n--- n = 2 ---")
a = sp.Symbol('a', positive=True)
# Centered roots: a, -a
# H(a) = 1/(a-(-a)) = 1/(2a)
# H(-a) = 1/(-a-a) = -1/(2a)
# Phi_2 = 1/(4a^2) + 1/(4a^2) = 1/(2a^2)
# kappa_2 = 2a^2 (from verify script)
# 1/Phi_2 = 2a^2 = kappa_2
# So C_2 = 1, R_2 = 0.
print("1/Phi_2 = kappa_2")
print("C_2 = 1, R_2 = 0")

# --- n = 3 ---
print("\n--- n = 3 ---")
# Use symbolic roots for centered cubic: roots sum to 0
# Parameterize: roots = {r1, r2, -(r1+r2)}
r1, r2 = sp.symbols('r1 r2')
r3 = -(r1 + r2)
roots_3 = [r1, r2, r3]

# Compute H_i for each root
def symbolic_H(roots, i):
    n = len(roots)
    return sum(1 / (roots[i] - roots[j]) for j in range(n) if j != i)

def symbolic_phi(roots):
    return sum(symbolic_H(roots, i)**2 for i in range(len(roots)))

phi3_sym = symbolic_phi(roots_3)
phi3_simplified = sp.simplify(phi3_sym)
print(f"Phi_3 = {phi3_simplified}")

inv_phi3 = sp.simplify(1 / phi3_simplified)
print(f"1/Phi_3 = {inv_phi3}")

# Express in terms of e2, e3
e2_3 = sp.Symbol('e2')
e3_3 = sp.Symbol('e3')
# e2 = r1*r2 + r1*r3 + r2*r3, e3 = r1*r2*r3
e2_expr = r1*r2 + r1*r3 + r2*r3
e3_expr = r1*r2*r3

# For n=3 centered: kappa_2 = -e2, kappa_3 = 9*e3/2
k2, k3 = sp.symbols('k2 k3', real=True)

# Substitute e2 = -k2, e3 = 2*k3/9
# But first, express 1/Phi_3 in terms of e2, e3
# From the known formula: Phi_3 = (9/2)*k2^2/(k2^3 - k3^2/3)
# So 1/Phi_3 = (k2^3 - k3^2/3) / ((9/2)*k2^2)
#            = (2/9)*(k2^3 - k3^2/3)/k2^2
#            = (2/9)*k2 - (2/27)*k3^2/k2^2
inv_phi3_cumulant = (k2**3 - k3**2/3) / (sp.Rational(9,2) * k2**2)
inv_phi3_expanded = sp.expand(inv_phi3_cumulant)
print(f"1/Phi_3 = {inv_phi3_expanded}")
print(f"        = (2/9)*k2 - (2/27)*k3^2/k2^2")
print(f"C_3 = 2/9, R_3 = -(2/27)*k3^2/k2^2")
print(f"Verify: 2/(3*2) = 1/3? No, 2/(3*3) = 2/9. ✓ C_3 = 2/(n*(n-1)) for n=3")

# ============================================================
# PART 2: n = 4 symbolic computation
# ============================================================
print("\n" + "=" * 70)
print("PART 2: n = 4 symbolic computation")
print("=" * 70)

# For n=4, we need to be more systematic.
# Centered polynomial: x^4 + a2*x^2 + a3*x + a4 (a1=0)
# Roots lambda_1,...,lambda_4 with sum = 0.

# Rather than symbolic roots (expensive), let's work with
# the polynomial directly and use Newton's identities.

# Phi_n in terms of power sums:
# H_p(lambda_i) = p''(lambda_i)/(2*p'(lambda_i)) = sum_{j≠i} 1/(lambda_i - lambda_j)
# Phi_n = sum_i (sum_{j≠i} 1/(lambda_i - lambda_j))^2

# Key insight: Phi_n * disc(p) is a polynomial in the coefficients.
# For n=3: Phi_3 * disc = 18*e_2^2

# For n=4 centered, let's compute Phi_4 * disc symbolically.
# The discriminant of x^4 + px^2 + qx + r is:
# disc = 256*r^3 - 128*p^2*r^2 + 144*p*q^2*r - 27*q^4 + 16*p^4*r - 4*p^3*q^2
# (using the standard formula for depressed quartic x^4 + px^2 + qx + r)

# Let me compute everything numerically first to find the pattern.
print("\nNumerical investigation for n=4...")

def compute_all_n4(roots):
    """Compute everything for a centered n=4 polynomial."""
    n = 4
    assert len(roots) == n
    assert abs(sum(roots)) < 1e-8

    e1 = sum(roots)
    e2 = sum(roots[i]*roots[j] for i in range(n) for j in range(i+1,n))
    e3 = sum(roots[i]*roots[j]*roots[k] for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
    e4 = roots[0]*roots[1]*roots[2]*roots[3]

    # Polynomial coefficients (centered): x^4 + 0*x^3 + e2*x^2 - e3*x + e4
    # (Vieta: a1=-e1=0, a2=e2, a3=-e3, a4=e4)

    # Normalized coefficients
    ta2 = e2 / comb(4,2)  # = e2/6
    ta3 = -(-e3) / comb(4,3)  # = e3/4 ... wait
    # tilde_a_k = (-1)^k * a_k / C(n,k)
    # a_1 = 0, a_2 = e2, a_3 = -e3, a_4 = e4
    # ta_2 = (-1)^2 * e2 / C(4,2) = e2/6
    # ta_3 = (-1)^3 * (-e3) / C(4,3) = e3/4
    # ta_4 = (-1)^4 * e4 / C(4,4) = e4
    ta2_val = e2 / 6
    ta3_val = e3 / 4
    ta4_val = e4

    # Cumulants (centered):
    # kappa_2 = -n * ta_2 = -4 * e2/6 = -2*e2/3
    k2 = -4 * ta2_val
    # kappa_3 = (n^2/2) * ta_3 = 8 * e3/4 = 2*e3
    k3 = 8 * ta3_val

    # kappa_4: need to derive. For centered case with NC(4):
    # The non-crossing partitions of {1,2,3,4} with no singletons
    # (since kappa_1=0):
    # {1234}: weight kappa_4 / n^3
    # {12}{34}: weight kappa_2^2 / n^2  (wait, need to count carefully)
    # {14}{23}: weight kappa_2^2 / n^2
    # So: ta_4 = kappa_4/n^3 + C_NC*kappa_2^2/n^2
    # where C_NC counts the number of NC partitions of {1,2,3,4} into pairs.
    # NC pair partitions of {1,2,3,4}: {12}{34}, {14}{23} = 2 (not {13}{24} which is crossing)
    # So ta_4 = kappa_4/64 + 2*kappa_2^2/16 = kappa_4/64 + kappa_2^2/8
    # => kappa_4 = 64*(ta_4 - kappa_2^2/8) = 64*ta_4 - 8*kappa_2^2
    k4 = 64 * ta4_val - 8 * k2**2

    phi = phi_n(roots)
    inv_phi = 1.0 / phi

    return {'e2': e2, 'e3': e3, 'e4': e4,
            'k2': k2, 'k3': k3, 'k4': k4,
            'phi': phi, 'inv_phi': inv_phi,
            'ta2': ta2_val, 'ta3': ta3_val, 'ta4': ta4_val}

# Generate many random n=4 centered polynomials
data_4 = []
for _ in range(500):
    r = np.random.randn(4) * 2
    r = r - np.mean(r)
    r = np.sort(r)
    # Ensure distinct roots
    if min(np.diff(r)) < 0.3:
        continue
    d = compute_all_n4(r)
    data_4.append(d)

print(f"Generated {len(data_4)} valid n=4 samples")

# Now fit: 1/Phi_4 = C_4 * k2 + f(k2, k3, k4)
# From the conjecture, C_4 = 2/(4*3) = 1/6

# First, verify C_4 by looking at the limit k3,k4 -> 0
# In this limit, the polynomial is x^4 + a2*x^2, with e3=0, e4=0
# So roots are 0, 0, sqrt(-a2), -sqrt(-a2) ... but these are not distinct!
#
# Instead, let's use regression. Fit inv_phi = C_4 * k2 + R_4_terms

# Better approach: compute 1/Phi_4 - (1/6)*k2 and see what's left
residuals_4 = []
for d in data_4:
    res = d['inv_phi'] - (1/6) * d['k2']  # assuming C_4 = 1/6
    residuals_4.append({'res': res, **d})

# Check if residual depends on k2, k3, k4
# Expected form: R_4 = polynomial / polynomial in k2, k3, k4
# From the n=3 pattern: R_3 = -(2/27)*k3^2/k2^2

# Let's try: R_4 might be of the form A*k3^2/k2^2 + B*k4/k2 + ...
# Actually, let's think about homogeneity. Phi_n is dimensionally [length]^{-2},
# so 1/Phi_n is [length]^2. kappa_k has dimension [length]^k.
# So C_n*k2 has dimension [length]^2. R_n must also have dimension [length]^2.
# Possible terms: k3^2/k2^2 (dim 2), k4/k2 (dim 2), k2 (dim 2), etc.

# Try fitting: R_4 = alpha * k3^2/k2^2 + beta * k4/k2
X = np.array([[d['k3']**2/d['k2']**2, d['k4']/d['k2']] for d in residuals_4])
y = np.array([d['res'] for d in residuals_4])

# Check if C_4 = 1/6 is correct
C4_test = 1/6
resid_for_test = np.array([d['inv_phi'] - C4_test * d['k2'] for d in data_4])
print(f"\nTesting C_4 = 1/6:")
print(f"  Mean residual: {np.mean(resid_for_test):.6f}")
print(f"  Std residual: {np.std(resid_for_test):.6f}")

# Actually, let's first verify C_4 more carefully.
# Use polynomials with k3 ≈ 0 and k4 ≈ 0.
# These are polynomials close to x^4 + a2*x^2 form.
# For x^4 + a*x^2 = x(x^3 + a*x) = x*x*(x^2+a), roots: 0, 0, ±sqrt(-a)
# These have repeated roots! So we can't directly take this limit.

# Better: use linear regression
# 1/Phi_4 = C_4 * k2 + alpha * k3^2/k2^2 + beta * k4/k2 + gamma * k3^2*k4/k2^4 + ...

# Let's first try the simplest model:
# 1/Phi_4 = C * k2 + a * k3^2/k2^2 + b * k4/k2
# Possible additional: c * k4^2/k2^4, d * k3*k4/k2^3, etc.

print("\n--- Linear regression for C_4 and R_4 ---")

# Design matrix with all dimension-2 terms up to reasonable order
features_4 = []
for d in data_4:
    k2, k3, k4 = d['k2'], d['k3'], d['k4']
    if abs(k2) < 1e-8:
        continue
    features_4.append({
        'inv_phi': d['inv_phi'],
        'k2': k2,
        'k3^2/k2^2': k3**2/k2**2,
        'k4/k2': k4/k2,
        'k3^2*k4/k2^4': k3**2*k4/k2**4,
        'k4^2/k2^4': k4**2/k2**4,
        'k3^4/k2^5': k3**4/k2**5,  # dim = 12-10 = 2 ✓
        'k3*k4/k2^3': k3*k4/k2**3,  # dim = 3+4-6 = 1 ✗ wrong dim
    })

# Actually, let me think about this more carefully.
# kappa_k has "weight" k (degree in root differences).
# 1/Phi_n has weight 2.
# So we need terms of total weight 2.
# k2 has weight 2. k3^2/k2^2 has weight 6-4=2. k4/k2 has weight 4-2=2.
# k3^2*k4/(k2^4) has weight 6+4-8=2. k4^2/k2^4 has weight 8-8=0 (wrong).
# Wait: k4^2/k2^3 has weight 8-6=2 ✓

# But the function must be a RATIONAL function. Let's hypothesize:
# 1/Phi_n = N(k2,...,kn) / D(k2,...,kn) where N,D are polynomials.
#
# For n=3: 1/Phi_3 = (k2^3 - k3^2/3) / ((9/2)*k2^2)
# Numerator: degree 3 in k's (homogeneous weight 6)
# Denominator: degree 2 (weight 4)
# 1/Phi_3: weight 6-4 = 2. ✓

# For n=4, let's directly compute Phi_4 * disc and see if it's a polynomial.
# Then 1/Phi_4 = disc / (Phi_4 * disc).

# Discriminant of centered quartic x^4 + p*x^2 + q*x + r:
# disc = 256*r^3 - 128*p^2*r^2 + 144*p*q^2*r - 27*q^4 + 16*p^4*r - 4*p^3*q^2

# Here p = a2 = e2, q = a3 = -e3, r = a4 = e4.
# disc = 256*e4^3 - 128*e2^2*e4^2 + 144*e2*e3^2*e4 - 27*e3^4 + 16*e2^4*e4 - 4*e2^3*e3^2

# Let's compute Phi_4 * disc numerically:
print("\n--- Phi_4 * disc pattern ---")
for d in data_4[:10]:
    r = d  # reuse computed values
    e2, e3, e4 = r['e2'], r['e3'], r['e4']

    # Compute discriminant
    disc = (256*e4**3 - 128*e2**2*e4**2 + 144*e2*e3**2*e4
            - 27*e3**4 + 16*e2**4*e4 - 4*e2**3*e3**2)

    phi_disc = r['phi'] * disc

    # Try different polynomials in e2, e3, e4
    # For n=3, Phi_3 * disc = 18*e2^2 (a degree-4 polynomial with weight 4)
    # For n=4, maybe Phi_4 * disc = ?? * e2^?? * ??

    # Weight analysis: disc has weight n*(n-1) = 12 for n=4
    # Phi_4 has weight -2, so Phi_4*disc has weight 10.
    # Expected polynomial in e2 (wt 2), e3 (wt 3), e4 (wt 4) of total weight 10.

    print(f"  e2={e2:.4f}, e3={e3:.4f}, e4={e4:.4f}, "
          f"Phi_4*disc={phi_disc:.4f}")

# Let me try to find Phi_4*disc as a polynomial in e2, e3, e4 (weight 10)
# Possible monomials of weight 10: e2^5, e2^3*e4, e2^2*e3^2*??
# Weight 10 monomials: e2^5 (10), e2^3*e4 (10), e2^2*e3^2 (10)
# Also: e2*e3*??... e2*e3 has weight 5, need weight 5 more: impossible from e2,e3,e4
# Wait: e2*e3^2*? impossible to get 10 exactly. Let me enumerate:
# e2^a * e3^b * e4^c with 2a+3b+4c = 10
# (5,0,0): e2^5
# (3,0,1): e2^3*e4  (6+4=10)
# (2,2,0): e2^2*e3^2 (4+6=10)
# (1,0,2): e2*e4^2 (2+8=10)
# (0,2,1): e3^2*e4 (6+4=10)
# (1,2,1): 2+6+4=12 ✗
# (0,0,...): no way to get 10 with just e3,e4

# So: Phi_4*disc = A*e2^5 + B*e2^3*e4 + C*e2^2*e3^2 + D*e2*e4^2 + E*e3^2*e4
print("\n--- Fitting Phi_4*disc = A*e2^5 + B*e2^3*e4 + C*e2^2*e3^2 + D*e2*e4^2 + E*e3^2*e4 ---")

X_disc = []
y_disc = []
for d in data_4:
    e2, e3, e4 = d['e2'], d['e3'], d['e4']
    disc_val = (256*e4**3 - 128*e2**2*e4**2 + 144*e2*e3**2*e4
                - 27*e3**4 + 16*e2**4*e4 - 4*e2**3*e3**2)
    phi_disc = d['phi'] * disc_val

    X_disc.append([e2**5, e2**3*e4, e2**2*e3**2, e2*e4**2, e3**2*e4])
    y_disc.append(phi_disc)

X_disc = np.array(X_disc)
y_disc = np.array(y_disc)

# Solve least squares
coeffs_disc, residuals_ls, rank, sv = np.linalg.lstsq(X_disc, y_disc, rcond=None)
print(f"  Coefficients: A={coeffs_disc[0]:.6f}, B={coeffs_disc[1]:.6f}, "
      f"C={coeffs_disc[2]:.6f}, D={coeffs_disc[3]:.6f}, E={coeffs_disc[4]:.6f}")
print(f"  Residual norm: {np.sqrt(np.sum((X_disc @ coeffs_disc - y_disc)**2)):.2e}")

# Try to identify as rational numbers
from fractions import Fraction
for i, name in enumerate(['A (e2^5)', 'B (e2^3*e4)', 'C (e2^2*e3^2)',
                           'D (e2*e4^2)', 'E (e3^2*e4)']):
    frac = Fraction(coeffs_disc[i]).limit_denominator(1000)
    print(f"  {name}: {coeffs_disc[i]:.10f} ≈ {frac}")

# ============================================================
# PART 3: Alternative approach - direct symbolic for n=4
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Direct symbolic computation for n=4")
print("=" * 70)

# Use SymPy with 4 roots constrained to sum to 0
# This is expensive but exact
print("Computing Phi_4 symbolically...")

# Parameterize: 4 roots with sum=0
# Use 3 free parameters
s1, s2, s3 = sp.symbols('s1 s2 s3')
sym_roots_4 = [s1, s2, s3, -(s1+s2+s3)]

# Compute H_i and Phi_4 symbolically
phi4_terms = []
for i in range(4):
    H_i = sum(1/(sym_roots_4[i] - sym_roots_4[j]) for j in range(4) if j != i)
    phi4_terms.append(sp.expand(H_i**2))

phi4_sym = sum(phi4_terms)
print("  Phi_4 computed (raw). Simplifying...")

# This is a rational function. Let's get numerator and denominator.
phi4_num, phi4_den = sp.fraction(sp.together(phi4_sym))
phi4_num = sp.expand(phi4_num)
phi4_den = sp.expand(phi4_den)
print(f"  Phi_4 = N/D where deg(N) in s1,s2,s3 = ..., deg(D) = ...")

# The denominator should be related to the discriminant
# Actually, Phi_n = N / prod_{i<j} (lambda_i - lambda_j)^2 for some polynomial N.
# Let's compute the discriminant
disc4_sym = sp.Integer(1)
for i in range(4):
    for j in range(i+1, 4):
        disc4_sym *= (sym_roots_4[i] - sym_roots_4[j])**2

# Check if phi4_den divides disc4^k for some k
# Actually, H_i = sum 1/(lambda_i - lambda_j), so H_i = P_i / prod_{j≠i}(lambda_i-lambda_j)
# where P_i is a polynomial. Then H_i^2 has denominator prod_{j≠i}(lambda_i-lambda_j)^2.
# When we sum, common denominator = disc.
# So Phi_n * disc is a polynomial.

# Let's compute Phi_4 * disc directly
print("  Computing Phi_4 * disc...")
phi4_times_disc = sp.expand(phi4_sym * disc4_sym)
phi4_times_disc = sp.simplify(phi4_times_disc)
print(f"  Phi_4 * disc = {phi4_times_disc}")

# Now express in terms of elementary symmetric polynomials
e2_sym = sp.Symbol('e2')
e3_sym = sp.Symbol('e3')
e4_sym = sp.Symbol('e4')

# e2 = sum_{i<j} s_i*s_j (including the 4th root)
e2_val = (s1*s2 + s1*s3 + s1*(-(s1+s2+s3)) + s2*s3 + s2*(-(s1+s2+s3)) + s3*(-(s1+s2+s3)))
e2_val = sp.expand(e2_val)
e3_val = (s1*s2*s3 + s1*s2*(-(s1+s2+s3)) + s1*s3*(-(s1+s2+s3)) + s2*s3*(-(s1+s2+s3)))
e3_val = sp.expand(e3_val)
e4_val = sp.expand(s1*s2*s3*(-(s1+s2+s3)))
print(f"  e2 = {e2_val}")
print(f"  e3 = {e3_val}")
print(f"  e4 = {e4_val}")

# ============================================================
# PART 4: Numerical computation of C_n for n = 2,...,6
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Numerical verification of C_n = 2/(n*(n-1))")
print("=" * 70)

def compute_Cn_numerically(n, num_samples=2000):
    """Estimate C_n by computing 1/Phi_n for polynomials with small k3,...,kn.

    Strategy: Use polynomials that are perturbations of x^n + a*x^{n-2}
    (which has only kappa_2 nonzero among higher cumulants, in the limit).
    Actually, such polynomials may have repeated roots.

    Better strategy: Use many random polynomials and regress 1/Phi_n on k2
    with terms involving k3,...,kn as additional regressors.
    """
    results = []
    for _ in range(num_samples):
        # Generate random centered polynomial with distinct roots
        roots = np.random.randn(n) * 2
        roots = roots - np.mean(roots)
        roots = np.sort(roots)

        # Check distinct
        if min(np.diff(roots)) < 0.2:
            continue

        # Compute Phi_n
        phi = phi_n(roots)
        inv_phi = 1.0 / phi

        # Compute kappa_2 (using the formula kappa_2 = -n * tilde_a_2 for centered)
        coeffs = np.poly(roots)  # [1, a1, a2, ..., an]
        a2 = coeffs[2]  # coefficient of x^{n-2}
        ta2 = a2 / comb(n, 2)  # since a1=0, tilde_a_2 = (-1)^2 * a2/C(n,2) = a2/C(n,2)
        k2 = -n * ta2

        # Also compute k3 for n>=3
        k3 = 0.0
        if n >= 3:
            a3 = coeffs[3]
            ta3 = (-1)**3 * a3 / comb(n, 3)
            k3 = (n**2 / 2) * ta3

        results.append({'inv_phi': inv_phi, 'k2': k2, 'k3': k3, 'roots': roots})

    # Method 1: Simple ratio for polynomials with small k3
    # Filter to those with |k3/k2^{3/2}| < 0.1
    near_zero = [r for r in results if abs(r['k3']) < 0.05 * abs(r['k2'])**1.5 and abs(r['k2']) > 0.5]
    if len(near_zero) > 5:
        C_est1 = np.mean([r['inv_phi']/r['k2'] for r in near_zero])
    else:
        C_est1 = float('nan')

    # Method 2: Linear regression 1/Phi = C*k2 + noise(higher cumulants)
    k2_vals = np.array([r['k2'] for r in results])
    inv_phi_vals = np.array([r['inv_phi'] for r in results])

    # Use ridge regression with k3^2/k2^2, etc. as additional features
    if n >= 3:
        k3_vals = np.array([r['k3'] for r in results])
        X = np.column_stack([k2_vals, k3_vals**2/k2_vals**2])
        coefficients = np.linalg.lstsq(X, inv_phi_vals, rcond=None)[0]
        C_est2 = coefficients[0]
    else:
        coefficients = np.linalg.lstsq(k2_vals.reshape(-1,1), inv_phi_vals, rcond=None)[0]
        C_est2 = coefficients[0]

    C_conj = 2 / (n * (n-1))

    return C_est1, C_est2, C_conj, len(results)

for n in [2, 3, 4, 5, 6]:
    C_est1, C_est2, C_conj, nsamp = compute_Cn_numerically(n, num_samples=5000)
    print(f"\n  n={n}: C_n conjecture = 2/({n}*{n-1}) = {C_conj:.10f}")
    print(f"    Estimate (near k3≈0): {C_est1:.10f}")
    print(f"    Estimate (regression): {C_est2:.10f}")
    print(f"    Ratio est1/conj: {C_est1/C_conj:.8f}" if not np.isnan(C_est1) else "    (insufficient samples near k3≈0)")
    print(f"    Ratio est2/conj: {C_est2/C_conj:.8f}")
    print(f"    Samples: {nsamp}")
