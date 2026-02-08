"""
VERIFIER-9: Rigorous verification of R_4 formula and superadditivity.
Adversarial checks on all claims in HANDOFF.md regarding cumulant decomposition.
"""

import numpy as np
from itertools import combinations
from math import comb, factorial
from fractions import Fraction
import sympy as sp
from sympy import symbols, expand, simplify, Rational, Poly, resultant, factor
from sympy import Symbol, poly, cancel, together, apart, collect, numer, denom

# ============================================================
# PART 0: Core numerical functions
# ============================================================

def elementary_symmetric(roots, k):
    """e_k of the roots."""
    n = len(roots)
    if k == 0:
        return 1.0
    if k > n:
        return 0.0
    return sum(np.prod(combo) for combo in combinations(roots, k))

def disc_from_roots(roots):
    """Discriminant = prod_{i<j} (r_i - r_j)^2."""
    n = len(roots)
    d = 1.0
    for i in range(n):
        for j in range(i + 1, n):
            d *= (roots[i] - roots[j]) ** 2
    return d

def phi_n(roots):
    """Phi_n = sum_i H_p(lambda_i)^2."""
    n = len(roots)
    total = 0.0
    for i in range(n):
        H = sum(1.0 / (roots[i] - roots[j]) for j in range(n) if j != i)
        total += H ** 2
    return total

def mss_convolve(roots_p, roots_q):
    """MSS finite free convolution."""
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

def finite_free_cumulants(roots, n_val):
    """Compute finite free cumulants kappa_1, ..., kappa_n from roots.

    Following Arizmendi-Perales convention:
    tilde_a_k = (-1)^k * a_k / C(n,k)  where a_k are coefficients of x^{n-k}

    Then kappa_k are defined recursively via moment-cumulant relations.
    For centered polynomials (kappa_1=0):
    kappa_2 = -n * tilde_a_2
    kappa_3 = n^2/2 * tilde_a_3  (when tilde_a_1=0)
    kappa_4 = -n^3/6 * (tilde_a_4 - 3*tilde_a_2^2)  (when tilde_a_1=0)  <-- NEED TO VERIFY
    """
    coeffs = np.poly(roots)  # [1, a_1, ..., a_n]

    # tilde_a_k = (-1)^k * a_k / C(n,k)
    tilde_a = [0.0] * (n_val + 1)
    tilde_a[0] = 1.0
    for k in range(1, n_val + 1):
        tilde_a[k] = ((-1) ** k * coeffs[k]) / comb(n_val, k)

    # Compute kappa from tilde_a using the recursive formula
    # kappa_1 = tilde_a_1
    # kappa_2 = -n * (tilde_a_2 - tilde_a_1^2)
    # kappa_3 = n^2/2 * (tilde_a_3 - 3*tilde_a_2*tilde_a_1 + 2*tilde_a_1^3)
    # kappa_4 = -n^3/6 * (tilde_a_4 - 4*tilde_a_3*tilde_a_1 - 3*tilde_a_2^2 + 12*tilde_a_2*tilde_a_1^2 - 6*tilde_a_1^4)

    kappa = [0.0] * (n_val + 1)
    kappa[1] = tilde_a[1]

    if n_val >= 2:
        kappa[2] = -n_val * (tilde_a[2] - tilde_a[1] ** 2)

    if n_val >= 3:
        kappa[3] = n_val ** 2 / 2 * (tilde_a[3] - 3 * tilde_a[2] * tilde_a[1] + 2 * tilde_a[1] ** 3)

    if n_val >= 4:
        kappa[4] = -n_val ** 3 / 6 * (
            tilde_a[4] - 4 * tilde_a[3] * tilde_a[1] - 3 * tilde_a[2] ** 2
            + 12 * tilde_a[2] * tilde_a[1] ** 2 - 6 * tilde_a[1] ** 4
        )

    if n_val >= 5:
        kappa[5] = n_val ** 4 / 24 * (
            tilde_a[5] - 5 * tilde_a[4] * tilde_a[1] - 10 * tilde_a[3] * tilde_a[2]
            + 20 * tilde_a[3] * tilde_a[1] ** 2 + 30 * tilde_a[2] ** 2 * tilde_a[1]
            - 60 * tilde_a[2] * tilde_a[1] ** 3 + 24 * tilde_a[1] ** 5
        )

    return kappa


# ============================================================
# PART A: VERIFY CUMULANT ADDITIVITY FIRST
# ============================================================
print("=" * 70)
print("PART A0: VERIFY CUMULANT ADDITIVITY UNDER MSS CONVOLUTION")
print("=" * 70)

np.random.seed(42)
for n_val in [3, 4]:
    print(f"\n  n = {n_val}:")
    max_err = 0.0
    for trial in range(20):
        roots_p = np.sort(np.random.randn(n_val) * 2)
        roots_p -= np.mean(roots_p)
        roots_q = np.sort(np.random.randn(n_val) * 2)
        roots_q -= np.mean(roots_q)

        if min(np.diff(roots_p)) < 0.2 or min(np.diff(roots_q)) < 0.2:
            continue

        try:
            roots_r = mss_convolve(roots_p, roots_q)

            kp = finite_free_cumulants(roots_p, n_val)
            kq = finite_free_cumulants(roots_q, n_val)
            kr = finite_free_cumulants(roots_r, n_val)

            for j in range(2, n_val + 1):
                err = abs(kr[j] - (kp[j] + kq[j]))
                max_err = max(max_err, err)
                if err > 1e-6:
                    print(f"    WARNING: kappa_{j} not additive! k_r={kr[j]:.8f}, k_p+k_q={kp[j]+kq[j]:.8f}, err={err:.2e}")
        except:
            pass

    print(f"  Max additivity error: {max_err:.2e}")

# ============================================================
# PART A1: SYMBOLIC DERIVATION OF 1/Phi_4
# ============================================================
print("\n" + "=" * 70)
print("PART A1: SYMBOLIC DERIVATION OF Phi_4 AND 1/Phi_4")
print("=" * 70)

# For n=4 centered polynomial p(x) = x^4 + e2*x^2 + e3*x + e4
# (with e1 = 0 since centered)
#
# We need Phi_4 = sum_i H(lambda_i)^2 where H(lambda_i) = p''(lambda_i)/(2*p'(lambda_i))
#
# Key identity: Phi_n * disc = N_n where N_n is a polynomial in e_k
#
# For general n: Phi_n = N_n / disc
# So 1/Phi_n = disc / N_n
#
# We need N_4 and disc_4 as polynomials in e2, e3, e4 (with e1=0)

# Let me compute N_4 = Phi_4 * disc symbolically
# Actually, let me verify the CLAIMED formulas numerically first

# The discriminant of x^4 + e2*x^2 + e3*x + e4 (monic, centered, e1=0)
# Standard formula for disc of quartic x^4 + px^2 + qx + r:
# disc = 256r^3 - 128p^2*r^2 + 144p*q^2*r - 27*q^4 + 16p^4*r - 4p^3*q^2
# Wait, this needs careful checking. Let me use sympy.

x, e2_s, e3_s, e4_s = symbols('x e2 e3 e4')
p4 = x**4 + e2_s*x**2 + e3_s*x + e4_s

# Compute discriminant symbolically
p4_poly = Poly(p4, x)
disc_sym = p4_poly.discriminant()
disc_sym_expanded = expand(disc_sym)
print(f"\nDiscriminant of x^4 + e2*x^2 + e3*x + e4:")
print(f"  disc = {disc_sym_expanded}")

# Compute N_4 = Phi_4 * disc
# Phi_4 = sum_i (p''(lambda_i)/(2p'(lambda_i)))^2
# = (1/4) * sum_i (p''(lambda_i))^2 / (p'(lambda_i))^2
#
# Note: p'(lambda_i) = prod_{j!=i} (lambda_i - lambda_j)
# So p'(lambda_i)^2 = prod_{j!=i} (lambda_i - lambda_j)^2
# And disc = prod_{i<j} (lambda_i - lambda_j)^2 = prod_i prod_{j>i} (lambda_i - lambda_j)^2
#
# Actually, disc = (-1)^{n(n-1)/2} * Res(p, p') / leading_coeff
# For monic, disc = (-1)^{n(n-1)/2} * prod_{i != j} (lambda_i - lambda_j)
# Actually disc = prod_{i < j} (lambda_i - lambda_j)^2
#
# And prod_i p'(lambda_i) = prod_i prod_{j!=i} (lambda_i - lambda_j)
# = prod_{i != j} (lambda_i - lambda_j) = (-1)^{n(n-1)/2} * disc (up to sign)
#
# The key identity is:
# Phi_n * disc = sum_i (p''(lambda_i))^2 / (4) * disc / p'(lambda_i)^2
#
# Let me instead compute this using resultants.
# Phi_4 * disc can be expressed as the resultant of certain polynomials.
#
# Better approach: compute numerically for many random polynomials and fit the polynomial.

print("\n  Computing N_4 = Phi_4 * disc numerically and fitting polynomial...")

# Strategy: sample many centered quartics, compute Phi_4 * disc, and identify
# the polynomial in (e2, e3, e4).
# We expect N_4 to be a polynomial of degree 8 in the e_k (since disc has degree 6
# for quartic, and we expect Phi_4 * disc to be polynomial).

# Actually, let's use a smarter approach. We know:
# Phi_n = Res(p'', p) / (4 * Res_related)... This is getting complicated.
# Let me just compute symbolically using sympy.

# Phi_4 = sum_i (p''(ri)/(2*p'(ri)))^2 where r_i are roots of p
# = (1/4) * sum_i (p''(ri)/p'(ri))^2
#
# p''(ri)/p'(ri) = d/dx[log p'(x)]|_{x=ri}
# But we need to express Phi_4 as a rational function of the coefficients.

# Use the identity: sum_i f(ri)/p'(ri) = [coefficient of x^{n-1} in the
# polynomial part of f(x)/p(x)] via partial fractions.
# But for f = (p''/p')^2 * p'^2 /4 = (p'')^2/4, this is more complex.

# Let me try a direct symbolic computation.
# For small n=4, we can compute with explicit roots... but that's circular.
#
# Instead, use Newton's identities to express sum_i g(lambda_i) in terms of e_k.

# Alternative: compute Phi_4 as a rational function of e_k by computing
# sum_i (p''(lambda_i))^2 * prod_{j!=i} (lambda_i - lambda_j)^2 / (4 * disc)
# = sum_i (p''(lambda_i))^2 * (disc / (lambda_i - ...)^2) / (4 * disc)
# Hmm, this is circular again.

# Best approach for n=4: use the fact that
# sum_i f(lambda_i) * g(lambda_i) = [extractable from resultants]
#
# Actually the cleanest is:
# N_4 = Phi_4 * disc = (1/4) * sum_i (p''(lambda_i))^2 * prod_{j!=i}(lambda_i - lambda_j)^2
# = (1/4) * sum_i (p''(lambda_i) * prod_{j!=i}(lambda_i - lambda_j))^2 / prod_{j!=i}(lambda_i-lambda_j)^2 * prod_{j!=i}(lambda_i-lambda_j)^2
# ... no, simpler:
# N_4 = (1/4) * sum_i (p''(lambda_i))^2 * (p'(lambda_i))^{-2} * disc
#
# Hmm, we need: sum_i (p''(lambda_i)/p'(lambda_i))^2.
# Let's define Q(x) = p''(x)/p'(x) = sum_j 1/(x - mu_j) where mu_j are roots of p'.
# Wait, that's not right. p'(x) = 4*prod(x - mu_j) for j=1..3 (critical points).
# So p'(x)/4 = prod(x - mu_j), and p''(x)/p'(x) = sum_j 1/(x - mu_j) + ... no.
# p'(x) = 4x^3 + 2e2*x + e3
# p''(x) = 12x^2 + 2e2
#
# H(lambda_i) = p''(lambda_i)/(2*p'(lambda_i))
#
# We can compute sum_i H(lambda_i)^2 using the identity:
# sum_i f(lambda_i)/p'(lambda_i) = sum of residues of f(x)/p(x)
# When deg(f) < n, this equals the coefficient extraction.
#
# Specifically: sum_i f(lambda_i)/p'(lambda_i) = [sum of residues of f(x)/p(x)]
# For f polynomial with deg(f) < deg(p) = n.
#
# We want: Phi_4 = sum_i (p''(lambda_i))^2 / (4 * p'(lambda_i)^2)
# = (1/4) * sum_i [(p''(lambda_i))^2 / p'(lambda_i)] * (1/p'(lambda_i))
# = (1/4) * sum_i g(lambda_i) / p'(lambda_i)  where g(x) = (p''(x))^2 / p'(x)
# But g(x) is a rational function, not a polynomial.
#
# Instead, let's just do it numerically. Generate many centered quartics,
# compute Phi_4 * disc, and determine the polynomial.

print("\n  Strategy: numerical fitting of N_4 = Phi_4 * disc as polynomial in (e2, e3, e4)")

from numpy.polynomial.polynomial import polyvander

def compute_phi_disc(e2_val, e3_val, e4_val):
    """Compute Phi_4 * disc for centered quartic x^4 + e2*x^2 + e3*x + e4."""
    coeffs = [1.0, 0.0, e2_val, e3_val, e4_val]
    roots = np.roots(coeffs)
    if np.max(np.abs(np.imag(roots))) > 1e-8:
        return None
    roots = np.sort(np.real(roots))
    if min(np.diff(roots)) < 1e-10:
        return None
    phi = phi_n(roots)
    d = disc_from_roots(roots)
    return phi * d

# Generate monomials in (e2, e3, e4) up to certain degree
# disc is degree 6 in coefficients, Phi_4 should be degree 2 in reciprocal,
# so N_4 = Phi_4 * disc should be degree ???
#
# disc(quartic) has degree 6 in coefficients (for x^4+...+e4)
# Phi_4 is a sum of rational functions of roots
# N_4 should be polynomial in e2, e3, e4
#
# For n=3: disc is degree 4 (-4e2^3 - 27e3^2), N_3 = 18*e2^2 (degree 2)
# so Phi_3 has "degree -2" in some sense.
#
# For n=4, let's figure out the degree of N_4.
# H(lambda_i) ~ 1/gap, so Phi_4 ~ 1/gap^2 ~ 1/(root spread)^2
# disc ~ (root spread)^{n(n-1)} for uniform spacing
# For n=4: disc ~ gap^12 (since 4*3/2 = 6 pairs, each contributing gap^2)
# Phi_4 ~ gap^{-2}
# So N_4 ~ gap^{10}, which in terms of coefficients...
# e2 ~ gap^2 (sum of products of pairs of roots), e3 ~ gap^3, e4 ~ gap^4
# N_4 ~ gap^{10} should be degree ... this is getting complicated.
#
# Let me just try fitting with monomials up to degree 8.

def monomial_basis(e2, e3, e4, max_deg):
    """Generate all monomials e2^a * e3^b * e4^c with a+b+c <= max_deg."""
    basis = []
    exponents = []
    for total in range(max_deg + 1):
        for a in range(total + 1):
            for b in range(total - a + 1):
                c = total - a - b
                basis.append(e2 ** a * e3 ** b * e4 ** c)
                exponents.append((a, b, c))
    return np.array(basis), exponents

# Sample many centered quartics with real distinct roots
np.random.seed(12345)
samples = []
values = []

for _ in range(500):
    # Generate random quartic with real roots
    roots = np.sort(np.random.randn(4) * 2)
    roots -= np.mean(roots)  # center

    if min(np.diff(roots)) < 0.3:
        continue

    e2 = elementary_symmetric(roots, 2)
    e3 = -elementary_symmetric(roots, 3)  # Note: for x^4 + 0*x^3 + a2*x^2 + a3*x + a4
    # p(x) = prod(x - ri) = x^4 - (sum ri)x^3 + (sum ri*rj)x^2 - (sum ri*rj*rk)x + prod ri
    # centered: sum ri = 0
    # a2 = sum_{i<j} ri*rj = e_2(roots)
    # a3 = -(sum_{i<j<k} ri*rj*rk) = -e_3(roots)  NO WAIT
    # Vieta: if p(x) = x^4 + c1*x^3 + c2*x^2 + c3*x + c4
    # c1 = -sum(ri), c2 = sum_{i<j}(ri*rj), c3 = -sum_{i<j<k}(ri*rj*rk), c4 = prod(ri)
    # So c2 = e2, c3 = -e3, c4 = e4
    # In our notation: p(x) = x^4 + e2_coeff*x^2 + e3_coeff*x + e4_coeff
    # where e2_coeff = c2 = e2(roots), e3_coeff = c3 = -e3(roots), e4_coeff = c4 = e4(roots)

    e2_coeff = sum(roots[i]*roots[j] for i in range(4) for j in range(i+1,4))
    e3_coeff = -sum(roots[i]*roots[j]*roots[k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4))
    e4_coeff = roots[0]*roots[1]*roots[2]*roots[3]

    phi = phi_n(roots)
    d = disc_from_roots(roots)
    N = phi * d

    samples.append([e2_coeff, e3_coeff, e4_coeff])
    values.append(N)

samples = np.array(samples)
values = np.array(values)
print(f"\n  Collected {len(samples)} samples")

# Fit polynomial: try degrees 4, 5, 6, 7, 8
for max_deg in range(4, 9):
    A_rows = []
    exp_list = []
    for total in range(max_deg + 1):
        for a in range(total + 1):
            for b in range(total - a + 1):
                c = total - a - b
                exp_list.append((a, b, c))

    A = np.zeros((len(samples), len(exp_list)))
    for idx, (e2v, e3v, e4v) in enumerate(samples):
        for j, (a, b, c) in enumerate(exp_list):
            A[idx, j] = e2v ** a * e3v ** b * e4v ** c

    # Least squares fit
    coeffs_fit, residuals, rank, sv = np.linalg.lstsq(A, values, rcond=None)

    predicted = A @ coeffs_fit
    max_rel_err = np.max(np.abs(predicted - values) / np.abs(values))

    if max_rel_err < 1e-6:
        print(f"\n  Degree {max_deg}: max relative error = {max_rel_err:.2e} -- GOOD FIT")

        # Print significant terms
        print(f"  N_4 = Phi_4 * disc =")
        terms = []
        for j, (a, b, c) in enumerate(exp_list):
            if abs(coeffs_fit[j]) > 1e-6:
                coeff = coeffs_fit[j]
                # Try to identify as simple fraction
                for denom_try in [1, 2, 3, 4, 6, 8, 9, 12, 16, 18, 24, 27, 36, 48, 54, 72, 108]:
                    val = coeff * denom_try
                    if abs(val - round(val)) < 1e-3:
                        num = round(val)
                        if denom_try == 1:
                            coeff_str = f"{num}"
                        else:
                            coeff_str = f"{num}/{denom_try}"
                        break
                else:
                    coeff_str = f"{coeff:.6f}"

                mon = ""
                if a > 0:
                    mon += f"e2^{a}" if a > 1 else "e2"
                if b > 0:
                    mon += f"*e3^{b}" if b > 1 else "*e3"
                if c > 0:
                    mon += f"*e4^{c}" if c > 1 else "*e4"
                if not mon:
                    mon = "1"

                terms.append(f"    {coeff_str} * {mon}")

        for t in terms:
            print(t)
        break
    else:
        print(f"  Degree {max_deg}: max relative error = {max_rel_err:.2e}")


# ============================================================
# PART A2: VERIFY N_4 SYMBOLICALLY WITH SYMPY
# ============================================================
print("\n" + "=" * 70)
print("PART A2: SYMBOLIC COMPUTATION OF N_4 = Phi_4 * disc")
print("=" * 70)

# Use sympy to compute Phi_4 symbolically for quartic
# p(x) = (x-a)(x-b)(x-c)(x-d) centered: a+b+c+d=0

a_s, b_s, c_s, d_s = symbols('a b c d')

# Centered: d = -(a+b+c)
d_expr = -(a_s + b_s + c_s)

roots_sym = [a_s, b_s, c_s, d_expr]

# Compute H(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j)
def H_sym(roots, i):
    return sum(1/(roots[i] - roots[j]) for j in range(4) if j != i)

# Compute Phi_4 = sum H^2
# This will be a rational function of a,b,c
# N_4 = Phi_4 * disc where disc = prod_{i<j}(ri-rj)^2

print("  Computing symbolically (this may take a moment)...")

# Compute disc symbolically
disc_pairs = 1
for i in range(4):
    for j in range(i+1, 4):
        disc_pairs *= (roots_sym[i] - roots_sym[j])**2

# Compute Phi_4 * disc
# Phi_4 = sum_i H_i^2
# H_i = sum_{j!=i} 1/(ri - rj)
# H_i * prod_{j!=i}(ri - rj) = sum_{j!=i} prod_{k!=i,k!=j}(ri - rk)
# So H_i^2 * prod_{j!=i}(ri - rj)^2 = (sum_{j!=i} prod_{k!=i,k!=j}(ri - rk))^2
# And Phi_4 * disc = sum_i (sum_{j!=i} prod_{k!=i,k!=j}(ri - rk))^2 * prod_{j>i}...
# Actually this is getting complicated. Let me just compute directly.

# H_i = (p''(ri)) / (2*p'(ri))
# p'(ri) = prod_{j!=i}(ri - rj)
# Phi_4 * disc = (1/4) * sum_i (p''(ri))^2 * disc / p'(ri)^2
# = (1/4) * sum_i (p''(ri))^2 * prod_{j<k, (j,k)!=(i,?)} ...
#
# Simpler: disc = prod_i p'(lambda_i) * (-1)^{n(n-1)/2} / leading^{n-2}
# For monic degree 4: disc = (-1)^6 * prod_i p'(lambda_i) / 1 = prod_i p'(lambda_i)...
# Actually no. disc = product of (lambda_i - lambda_j)^2 over i<j.
# prod_i p'(lambda_i) = prod_i prod_{j!=i}(lambda_i - lambda_j)
# = prod_{i!=j}(lambda_i - lambda_j) = prod_{i<j}(lambda_i-lambda_j) * prod_{i>j}(lambda_i-lambda_j)
# = prod_{i<j}(lambda_i-lambda_j) * prod_{j<i}(lambda_i-lambda_j)
# = prod_{i<j}(lambda_i-lambda_j) * prod_{i<j}(lambda_j-lambda_i)  [swapping indices]
# = prod_{i<j}(lambda_i-lambda_j) * prod_{i<j}(-(lambda_i-lambda_j))
# = (-1)^{n(n-1)/2} * (prod_{i<j}(lambda_i-lambda_j))^2
# = (-1)^6 * disc = disc  for n=4
#
# So prod_i p'(lambda_i) = disc for n=4.
#
# Therefore: Phi_4 * disc = (1/4) * sum_i (p''(lambda_i))^2 * disc / p'(lambda_i)^2
# = (1/4) * sum_i (p''(lambda_i))^2 * prod_{j!=i}(lambda_i-lambda_j)^2 * ... no
#
# disc / p'(lambda_i)^2 = disc / (prod_{j!=i}(lambda_i - lambda_j))^2
#
# disc = prod_{k<l}(lambda_k - lambda_l)^2
# p'(lambda_i)^2 = (prod_{j!=i}(lambda_i - lambda_j))^2
#
# disc / p'(lambda_i)^2 = prod_{k<l, k!=i and l!=i}(lambda_k - lambda_l)^2
# This is the discriminant of the polynomial with root lambda_i removed.
#
# So N_4 = (1/4) * sum_i (p''(lambda_i))^2 * disc_{p/(x-lambda_i)}

# For n=4 centered: p(x) = x^4 + e2*x^2 + e3*x + e4
# p'(x) = 4x^3 + 2*e2*x + e3
# p''(x) = 12x^2 + 2*e2

# Let me compute this using sympy with explicit roots
# Since we need a polynomial in e2, e3, e4, let me try a different approach.

# Use the identity: if p has roots lambda_i, then
# sum_i f(lambda_i) / p'(lambda_i)^(2k-1) can be expressed via resultants.
# But specifically, sum_i f(lambda_i) / p'(lambda_i) = integral formula.

# Let me try the residue approach:
# sum_i f(lambda_i)/p'(lambda_i) = sum of residues of f(x)/p(x) at roots
# = -residue at infinity of f(x)/p(x) when deg(f) < n.
#
# For Phi_4: sum_i (p''(lambda_i))^2 / (4*p'(lambda_i)^2)
# This involves 1/p'(lambda_i)^2, not 1/p'(lambda_i).
#
# Key identity: sum_i f(lambda_i)/p'(lambda_i)^2 = -sum_i d/dx[f(x)/p(x)]|_{x=lambda_i} residue
# = coefficient in the partial fraction of f(x)/(p(x))^2
# Hmm, this is getting complicated. Let me just use numerical data with rational reconstruction.

# Actually, let me try the direct sympy approach with 3 free variables
# (using centered quartic with 3 free roots, 4th determined by centering).

print("  Using direct symbolic computation with 3 free roots (centered)...")

# This is slow for sympy, so let me use a cleaner approach.
# Compute N_4 numerically with high precision and use rational reconstruction.

from decimal import Decimal, getcontext
getcontext().prec = 50

# Sample at specific rational points for e2, e3, e4
# Then N_4 should be a polynomial with rational coefficients
# Use many sample points and solve the linear system

# First, determine which monomials e2^a * e3^b * e4^c appear in N_4
# From homogeneity: if roots scale by t, then e_k scales by t^k (for centered poly)
# So e2 -> t^2*e2, e3 -> t^3*e3, e4 -> t^4*e4
# H scales by t^{-1}, Phi scales by t^{-2}
# disc scales by t^{n(n-1)} = t^12
# N_4 = Phi_4 * disc scales by t^{-2+12} = t^{10}
# So N_4 must be homogeneous of degree 10 in the e_k with weights (2,3,4)
# i.e., 2a + 3b + 4c = 10

# Possible monomials with 2a + 3b + 4c = 10:
# a=5, b=0, c=0: e2^5 (weight 10) ✓
# a=2, b=2, c=0: e2^2*e3^2 (weight 4+6=10) ✓
# a=0, b=0, c=2: Wait, need 4c <= 10 and b <= (10-4c)/3 etc.
# Let me enumerate:
print("\n  Monomials with weight 2a+3b+4c = 10:")
monomials_10 = []
for a in range(6):  # 2a <= 10
    for b in range((10-2*a)//3 + 1):
        rem = 10 - 2*a - 3*b
        if rem >= 0 and rem % 4 == 0:
            c = rem // 4
            monomials_10.append((a, b, c))
            print(f"    e2^{a} * e3^{b} * e4^{c}  (weight {2*a+3*b+4*c})")

print(f"  Total: {len(monomials_10)} monomials")

# Now solve for coefficients using numerical samples
# Need at least len(monomials_10) samples with good conditioning

np.random.seed(999)
n_monoms = len(monomials_10)
n_samples_needed = n_monoms + 20  # over-determined for robustness

A_mat = []
b_vec = []

for _ in range(200):
    # Generate centered quartic with real simple roots
    roots = np.sort(np.random.randn(4) * 2)
    roots -= np.mean(roots)

    if min(np.diff(roots)) < 0.3:
        continue

    e2v = sum(roots[i]*roots[j] for i in range(4) for j in range(i+1,4))
    e3_sym = sum(roots[i]*roots[j]*roots[k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4))
    e3v = -e3_sym  # coefficient of x in p(x) = x^4 + e2*x^2 + e3*x + e4
    e4v = roots[0]*roots[1]*roots[2]*roots[3]

    phi = phi_n(roots)
    d = disc_from_roots(roots)
    N = phi * d

    row = []
    for (a, b, c) in monomials_10:
        row.append(e2v**a * e3v**b * e4v**c)

    A_mat.append(row)
    b_vec.append(N)

    if len(A_mat) >= n_samples_needed:
        break

A_mat = np.array(A_mat)
b_vec = np.array(b_vec)

print(f"\n  Solving {A_mat.shape[0]}x{A_mat.shape[1]} system...")
coeffs_N4, residuals, rank, sv = np.linalg.lstsq(A_mat, b_vec, rcond=None)

# Check fit quality
predicted = A_mat @ coeffs_N4
max_rel_err = np.max(np.abs(predicted - b_vec) / (np.abs(b_vec) + 1e-30))
print(f"  Max relative error of fit: {max_rel_err:.2e}")

print(f"\n  N_4 (= Phi_4 * disc) =")
for j, (a, b, c) in enumerate(monomials_10):
    if abs(coeffs_N4[j]) > 1e-4:
        # Identify rational coefficient
        coeff = coeffs_N4[j]
        identified = False
        for denom_try in range(1, 200):
            val = coeff * denom_try
            if abs(val - round(val)) < 0.01:
                num = round(val)
                from math import gcd
                g = gcd(abs(num), denom_try)
                num //= g
                den = denom_try // g
                if den == 1:
                    coeff_str = f"{num}"
                else:
                    coeff_str = f"{num}/{den}"
                identified = True
                break
        if not identified:
            coeff_str = f"{coeff:.8f}"

        mon = f"e2^{a}" if a > 0 else ""
        if b > 0:
            mon += ("*" if mon else "") + f"e3^{b}"
        if c > 0:
            mon += ("*" if mon else "") + f"e4^{c}"
        if not mon:
            mon = "1"

        print(f"    {coeff_str:>10s} * {mon}")

# ============================================================
# PART A3: CONVERT TO CUMULANTS
# ============================================================
print("\n" + "=" * 70)
print("PART A3: CONVERT N_4 AND disc TO CUMULANTS")
print("=" * 70)

# For n=4 centered:
# tilde_a_k = (-1)^k * a_k / C(4,k) where p(x) = x^4 + a_2*x^2 + a_3*x + a_4
# tilde_a_1 = 0 (centered)
# tilde_a_2 = a_2/C(4,2) = a_2/6   (a_2 = e2_coeff)
# tilde_a_3 = -a_3/C(4,3) = -a_3/4  (a_3 = e3_coeff)
# tilde_a_4 = a_4/C(4,4) = a_4/1 = a_4  (a_4 = e4_coeff)
#
# kappa_1 = 0
# kappa_2 = -4*(tilde_a_2) = -4*a_2/6 = -2*a_2/3
# kappa_3 = 16/2 * tilde_a_3 = 8*(-a_3/4) = -2*a_3
# kappa_4 = -64/6 * (tilde_a_4 - 3*tilde_a_2^2) = (-32/3)*(a_4 - 3*(a_2/6)^2) = (-32/3)*(a_4 - a_2^2/12)

# Let me verify these formulas numerically
print("  Verifying cumulant-coefficient relations for n=4...")

np.random.seed(777)
for trial in range(5):
    roots = np.sort(np.random.randn(4) * 2)
    roots -= np.mean(roots)

    if min(np.diff(roots)) < 0.3:
        continue

    a2 = sum(roots[i]*roots[j] for i in range(4) for j in range(i+1,4))
    e3_val = sum(roots[i]*roots[j]*roots[k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4))
    a3 = -e3_val
    a4 = roots[0]*roots[1]*roots[2]*roots[3]

    # Compute kappa from roots using our function
    kappa_from_func = finite_free_cumulants(roots, 4)

    # Compute kappa from formulas
    k2_formula = -2*a2/3
    k3_formula = -2*a3
    k4_formula = (-32/3)*(a4 - a2**2/12)

    print(f"  Trial {trial}:")
    print(f"    kappa_2: from_func={kappa_from_func[2]:.8f}, formula={k2_formula:.8f}, match={abs(kappa_from_func[2]-k2_formula)<1e-8}")
    print(f"    kappa_3: from_func={kappa_from_func[3]:.8f}, formula={k3_formula:.8f}, match={abs(kappa_from_func[3]-k3_formula)<1e-8}")
    print(f"    kappa_4: from_func={kappa_from_func[4]:.8f}, formula={k4_formula:.8f}, match={abs(kappa_from_func[4]-k4_formula)<1e-8}")

# Inverse relations: express a2, a3, a4 in terms of kappa
# k2 = -2*a2/3 => a2 = -3*k2/2
# k3 = -2*a3 => a3 = -k3/2
# k4 = (-32/3)*(a4 - a2^2/12) => a4 = a2^2/12 + k4*(-3/32) = a2^2/12 - 3*k4/32
# a4 = (-3*k2/2)^2/12 - 3*k4/32 = 9*k2^2/48 - 3*k4/32 = 3*k2^2/16 - 3*k4/32

print("\n  Inverse: a2 = -3*k2/2, a3 = -k3/2, a4 = 3*k2^2/16 - 3*k4/32")

# Verify inverse relations
for trial in range(3):
    roots = np.sort(np.random.randn(4) * 2)
    roots -= np.mean(roots)
    if min(np.diff(roots)) < 0.3:
        continue

    a2 = sum(roots[i]*roots[j] for i in range(4) for j in range(i+1,4))
    e3_val = sum(roots[i]*roots[j]*roots[k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4))
    a3 = -e3_val
    a4 = roots[0]*roots[1]*roots[2]*roots[3]

    kappa = finite_free_cumulants(roots, 4)
    k2, k3, k4 = kappa[2], kappa[3], kappa[4]

    a2_rec = -3*k2/2
    a3_rec = -k3/2
    a4_rec = 3*k2**2/16 - 3*k4/32

    print(f"  Trial {trial}: a2 err={abs(a2-a2_rec):.2e}, a3 err={abs(a3-a3_rec):.2e}, a4 err={abs(a4-a4_rec):.2e}")


# ============================================================
# PART B: CHECK C_n VALUES
# ============================================================
print("\n" + "=" * 70)
print("PART B: CHECK C_n = coefficient of k2 in 1/Phi_n when k3=...=kn=0")
print("=" * 70)

# When k3 = k4 = ... = kn = 0, the polynomial is a "free semicircular" type
# (only kappa_2 nonzero). For centered polynomials with only kappa_2:
# n=2: roots are {a, -a} with kappa_2 = 2a^2, Phi_2 = 1/(2a^2) = 1/kappa_2
#       so 1/Phi_2 = kappa_2. C_2 = 1.
# n=3: roots are equally spaced {-s, 0, s} with e2 = -s^2, e3=0
#       kappa_2 = -e2 = s^2
#       Phi_3 = 9/(2s^2) (from equal-gap formula)
#       1/Phi_3 = 2s^2/9 = (2/9)*kappa_2
#       C_3 = 2/9

# For general n, let's compute numerically.
# "Only kappa_2 nonzero" means the polynomial has the same root distribution
# as an equally-weighted convex combination of Dirac masses at equally spaced points,
# suitably normalized. For n=4 with only kappa_2, we need e3=0 and e4 = a2^2/12
# (from a4 = 3*k2^2/16 - 3*k4/32 with k4=0, so a4 = 3*k2^2/16 = 3*(-3k2/2... no
# Wait: a2 = -3*k2/2, a4 = 3*k2^2/16 (when k4=0)
# So p(x) = x^4 + (-3k2/2)*x^2 + 0*x + 3*k2^2/16

print("\n  Computing C_n for n=2,3,4,5...")

for n_val in [2, 3, 4, 5]:
    # Find the roots of the "pure kappa_2" polynomial
    # and compute C_n = (1/Phi_n) / kappa_2

    C_n_values = []

    for k2_val in [0.5, 1.0, 2.0, 3.0, 5.0]:
        # Construct polynomial with only kappa_2 = k2_val
        # Need to find a2, a3=0, ..., a_n from kappa relations

        if n_val == 2:
            # p(x) = x^2 + a2 where a2 = -k2_val (since kappa_2 = -2*a2/C(2,2)*... )
            # For n=2: tilde_a_2 = a_2/C(2,2) = a_2
            # kappa_2 = -2*(tilde_a_2) = -2*a_2
            # So a_2 = -k2_val/2
            a2 = -k2_val / 2
            coeffs = [1.0, 0.0, a2]
        elif n_val == 3:
            # kappa_2 = -3*tilde_a_2 = -3*(a_2/C(3,2)) = -3*a_2/3 = -a_2
            # So a_2 = -k2_val
            a2 = -k2_val
            coeffs = [1.0, 0.0, a2, 0.0]
        elif n_val == 4:
            # kappa_2 = -2*a_2/3 => a_2 = -3*k2_val/2
            # kappa_4 = 0 => a_4 = 3*k2_val^2/16
            a2 = -3 * k2_val / 2
            a4 = 3 * k2_val ** 2 / 16
            coeffs = [1.0, 0.0, a2, 0.0, a4]
        elif n_val == 5:
            # For n=5 centered, kappa_2 only:
            # tilde_a_2 = a_2/C(5,2) = a_2/10
            # kappa_2 = -5*tilde_a_2 = -5*a_2/10 = -a_2/2
            # So a_2 = -2*k2_val
            #
            # kappa_3 = 0: need tilde_a_3 = 0 => a_3 = 0
            # kappa_4 = 0: need tilde_a_4 - 3*tilde_a_2^2 = 0
            #   tilde_a_4 = a_4/C(5,4) = a_4/5
            #   3*tilde_a_2^2 = 3*(a_2/10)^2 = 3*a_2^2/100
            #   a_4/5 = 3*a_2^2/100 => a_4 = 15*a_2^2/100 = 3*a_2^2/20
            # kappa_5 = 0: tilde_a_5 = a_5/C(5,5) = a_5
            #   Need: tilde_a_5 - 10*tilde_a_3*tilde_a_2 + ... actually the formula is complex.
            #   For centered with only kappa_2: all higher tilde_a should satisfy
            #   certain moment-cumulant constraints.
            #   kappa_5 involves tilde_a_5 and lower terms, but since tilde_a_3=0 and tilde_a_1=0:
            #   kappa_5 = n^4/24 * (tilde_a_5 - 10*tilde_a_3*tilde_a_2 + ...)
            #   With tilde_a_3 = 0: kappa_5 = (5^4/24)*(tilde_a_5 + ... terms with tilde_a_3 or tilde_a_1)
            #   Most cross terms vanish. Need: tilde_a_5 = 0 (since no odd cumulants)
            #   So a_5 = 0.
            a2_5 = -2 * k2_val
            a4_5 = 3 * a2_5 ** 2 / 20
            coeffs = [1.0, 0.0, a2_5, 0.0, a4_5, 0.0]
        else:
            continue

        roots = np.roots(coeffs)

        if np.max(np.abs(np.imag(roots))) > 1e-6:
            continue  # Not real-rooted

        roots = np.sort(np.real(roots))
        if min(np.diff(roots)) < 1e-8:
            continue  # Not simple

        phi = phi_n(roots)
        C_n = (1.0 / phi) / k2_val
        C_n_values.append(C_n)

    if C_n_values:
        C_n_avg = np.mean(C_n_values)
        C_n_std = np.std(C_n_values)

        # Try to identify as simple fraction
        best_frac = None
        best_err = 1e10
        for num in range(-20, 21):
            for den in range(1, 100):
                if num == 0 and den == 1:
                    continue
                frac = num / den
                err = abs(C_n_avg - frac)
                if err < best_err:
                    best_err = err
                    best_frac = (num, den)

        from math import gcd
        g = gcd(abs(best_frac[0]), best_frac[1])
        best_frac = (best_frac[0]//g, best_frac[1]//g)

        conjectured = 2 / (n_val * (n_val - 1))

        print(f"  n={n_val}: C_n = {C_n_avg:.10f} (std={C_n_std:.2e})")
        print(f"           Best fraction: {best_frac[0]}/{best_frac[1]} = {best_frac[0]/best_frac[1]:.10f}")
        print(f"           Conjectured 2/(n(n-1)) = {conjectured:.10f}")
        print(f"           MATCH conjectured: {abs(C_n_avg - conjectured) < 1e-6}")
    else:
        print(f"  n={n_val}: Could not compute (no valid samples)")


# ============================================================
# PART B2: VERIFY C_3 = 2/9 or 1/3?
# ============================================================
print("\n" + "=" * 70)
print("PART B2: DETAILED CHECK OF C_3")
print("=" * 70)

# For n=3, equally spaced roots {-s, 0, s}:
# Phi_3 = 9/(2s^2) (known formula)
# kappa_2 = s^2 (for n=3, kappa_2 = -e_2 = -((-s)*0 + (-s)*s + 0*s) = -(0 - s^2 + 0) = s^2)
# 1/Phi_3 = 2s^2/9 = (2/9)*kappa_2
# So C_3 = 2/9
# But 2/(3*2) = 2/6 = 1/3 != 2/9
# The conjecture C_n = 2/(n(n-1)) predicts C_3 = 1/3, but actual C_3 = 2/9.

s_test = 2.0
roots_test = np.array([-s_test, 0.0, s_test])
phi_test = phi_n(roots_test)
k2_test = s_test**2
print(f"  n=3, roots={roots_test}, Phi_3={phi_test:.10f}, 1/Phi_3={1/phi_test:.10f}")
print(f"  kappa_2={k2_test:.10f}, C_3 = (1/Phi_3)/k2 = {1/(phi_test*k2_test):.10f}")
print(f"  Expected C_3 = 2/9 = {2/9:.10f}")
print(f"  Conjectured C_3 = 2/(3*2) = 1/3 = {1/3:.10f}")
print(f"  DISCREPANCY: C_3 = 2/9 ≠ 1/3 = 2/(n(n-1))")


# ============================================================
# PART C: VERIFY THE R_4 FORMULA
# ============================================================
print("\n" + "=" * 70)
print("PART C: VERIFY CLAIMED R_4 FORMULA")
print("=" * 70)

def R4_claimed(k2, k3, k4):
    """The claimed R_4 formula from HANDOFF."""
    num = -16*k2**3*k3**2 - 4*k2**2*k4**2 + 20*k2*k3**2*k4 - 8*k3**4 - k4**3
    den = 24*(16*k2**5 - 8*k2**2*k3**2 - k2*k4**2 + 2*k3**2*k4)
    return num / den

def inv_phi_from_roots(roots):
    """Compute 1/Phi_n directly from roots."""
    return 1.0 / phi_n(roots)

# For each random centered quartic:
# 1. Compute roots and 1/Phi_4
# 2. Compute kappa_2, kappa_3, kappa_4
# 3. Check if 1/Phi_4 = C_4*k2 + R_4(k2, k3, k4)

np.random.seed(54321)
print("\n  Testing R_4 formula against direct computation...")
C4_value = None  # Will be determined

max_err = 0.0
max_rel_err = 0.0
n_tested = 0

for trial in range(200):
    roots = np.sort(np.random.randn(4) * 2)
    roots -= np.mean(roots)

    if min(np.diff(roots)) < 0.3:
        continue

    inv_phi = inv_phi_from_roots(roots)
    kappa = finite_free_cumulants(roots, 4)
    k2, k3, k4 = kappa[2], kappa[3], kappa[4]

    if k2 < 0.1:  # skip very small k2
        continue

    # Determine C_4 from pure kappa_2 case (k3=k4=0 limit)
    # We already computed C_4 above. Let's test with the exact value.
    # From Part B, C_4 should be numerically determined.

    # Test with C_4 = 1/12 (as claimed in HANDOFF)
    C4_test = 1.0/12.0

    try:
        R4 = R4_claimed(k2, k3, k4)
        predicted = C4_test * k2 + R4
        err = abs(predicted - inv_phi)
        rel_err = err / abs(inv_phi) if abs(inv_phi) > 1e-10 else err

        max_err = max(max_err, err)
        max_rel_err = max(max_rel_err, rel_err)
        n_tested += 1

        if rel_err > 0.01:
            print(f"  Trial {trial}: 1/Phi_4={inv_phi:.8f}, C4*k2+R4={predicted:.8f}, "
                  f"rel_err={rel_err:.2e}")
            print(f"    k2={k2:.6f}, k3={k3:.6f}, k4={k4:.6f}")
    except:
        pass

print(f"\n  Tested {n_tested} polynomials with C_4 = 1/12:")
print(f"  Max absolute error: {max_err:.2e}")
print(f"  Max relative error: {max_rel_err:.2e}")
if max_rel_err < 1e-6:
    print(f"  RESULT: R_4 formula with C_4=1/12 is CORRECT")
else:
    print(f"  RESULT: R_4 formula with C_4=1/12 has ERRORS")

    # Try other C_4 values
    for C4_try in [2/12, 1/6, 2/9, 1/8, 1/10]:
        max_rel_err2 = 0.0
        for trial in range(200):
            roots = np.sort(np.random.randn(4) * 2)
            roots -= np.mean(roots)
            if min(np.diff(roots)) < 0.3:
                continue
            inv_phi = inv_phi_from_roots(roots)
            kappa = finite_free_cumulants(roots, 4)
            k2, k3, k4 = kappa[2], kappa[3], kappa[4]
            if k2 < 0.1:
                continue
            try:
                R4 = R4_claimed(k2, k3, k4)
                predicted = C4_try * k2 + R4
                rel_err = abs(predicted - inv_phi) / abs(inv_phi)
                max_rel_err2 = max(max_rel_err2, rel_err)
            except:
                pass
        print(f"    C_4 = {C4_try:.6f}: max_rel_err = {max_rel_err2:.2e}")


# ============================================================
# PART C2: DERIVE 1/Phi_4 FROM SCRATCH
# ============================================================
print("\n" + "=" * 70)
print("PART C2: DERIVE 1/Phi_4 FROM SCRATCH (no assumptions)")
print("=" * 70)

# Strategy: compute 1/Phi_4 = disc/N_4 in terms of cumulants
# We have N_4 from the fit above, and disc from sympy.
# Let's substitute a_k -> kappa_k relations.

# We know:
# a_2 = -3*k2/2  (coefficient of x^2)
# a_3 = -k3/2    (coefficient of x)
# a_4 = 3*k2^2/16 - 3*k4/32  (constant term)

# Let's substitute into disc and N_4 symbolically

k2_s, k3_s, k4_s = symbols('k2 k3 k4', real=True)

# a_2, a_3, a_4 in terms of kappas
a2_expr = -Rational(3,2)*k2_s
a3_expr = -Rational(1,2)*k3_s
a4_expr = Rational(3,16)*k2_s**2 - Rational(3,32)*k4_s

# disc of x^4 + a2*x^2 + a3*x + a4
# Using the formula computed by sympy above
# Let's recompute
print("  Computing disc in terms of cumulants...")

x_s = Symbol('x')
p4_kappa = x_s**4 + a2_expr*x_s**2 + a3_expr*x_s + a4_expr
p4_poly_kappa = Poly(p4_kappa, x_s)
disc_kappa = p4_poly_kappa.discriminant()
disc_kappa_expanded = expand(disc_kappa)
print(f"  disc = {disc_kappa_expanded}")

# Factor and simplify
disc_kappa_factored = factor(disc_kappa_expanded)
print(f"  disc (factored) = {disc_kappa_factored}")

# Now compute N_4 in terms of kappas
# We need to substitute the identified N_4 polynomial
# From the numerical fit, N_4 is a polynomial in (a2, a3, a4)
# which we need to convert to kappas.

# Actually, let me redo the computation more carefully.
# Let me compute N_4 directly: N_4 = Phi_4 * disc

# Strategy: express 1/Phi_4 = disc/N_4 in cumulants.
# Both disc and N_4 are polynomials in (a_2, a_3, a_4).
# We substitute a_k = f(kappa).
# Then 1/Phi_4 = disc_kappa / N4_kappa.

# We need N_4. Let me compute it symbolically using the approach:
# Phi_4 = sum_i H(lambda_i)^2
# and use Newton's identity or Cauchy-Binet.

# Alternative: Since we fitted N_4 numerically, let's identify the exact coefficients
# and substitute.

# From the fit, we have coefficients. Let me verify them more carefully
# using high-precision arithmetic.

# Let's use sympy for exact computation.
a2_s_coeff, a3_s_coeff, a4_s_coeff = symbols('a2 a3 a4')

# For p(x) = x^4 + a2*x^2 + a3*x + a4, compute Phi_4 * disc symbolically
# using resultant/subresultant methods.

# p'(x) = 4x^3 + 2*a2*x + a3
# p''(x) = 12*x^2 + 2*a2

# H(lambda_i) = p''(lambda_i)/(2*p'(lambda_i))
# Phi_4 = sum_i (p''(lambda_i)/(2*p'(lambda_i)))^2
# = (1/4) sum_i (p''(lambda_i))^2 / p'(lambda_i)^2

# Key identity: sum_i f(lambda_i)/p'(lambda_i) = -Res(f, p)/lc(p)^{deg_f - n + 1}
# when deg(f) < deg(p) via partial fractions.
# But we need sum_i f(lambda_i)^2/p'(lambda_i)^2 which is different.

# Use the identity:
# sum_i f(lambda_i)^2/p'(lambda_i)^2 = sum_i [f(lambda_i)/p'(lambda_i)] * [f(lambda_i)/p'(lambda_i)]
# This is related to the squared resolvent.

# Alternative approach: let's just compute with explicit roots using sympy.
# For a general centered quartic with roots a, b, c, d where d = -(a+b+c).

print("\n  Computing N_4 with explicit roots using sympy...")
print("  (This may take 30-60 seconds...)")

a_sym, b_sym, c_sym = symbols('a b c')
d_sym = -(a_sym + b_sym + c_sym)

roots_list = [a_sym, b_sym, c_sym, d_sym]

# H values
H_values = []
for i in range(4):
    h = sum(sp.Integer(1) / (roots_list[i] - roots_list[j]) for j in range(4) if j != i)
    H_values.append(h)

# Phi_4 = sum H^2
Phi4_sym = sum(h**2 for h in H_values)

# Discriminant = prod_{i<j} (r_i - r_j)^2
disc_sym_roots = sp.Integer(1)
for i in range(4):
    for j in range(i+1, 4):
        disc_sym_roots *= (roots_list[i] - roots_list[j])**2

# N_4 = Phi_4 * disc
N4_sym = Phi4_sym * disc_sym_roots

print("  Expanding N_4 (this is the slow step)...")
N4_expanded = expand(N4_sym)

# Now express in terms of elementary symmetric polynomials
# e1 = a+b+c+d = 0
# e2 = sum_{i<j} ri*rj
# e3 = sum_{i<j<k} ri*rj*rk
# e4 = a*b*c*d

# Compute e_k in terms of a, b, c (with d=-(a+b+c))
e1_val = sp.Integer(0)
e2_expr = a_sym*b_sym + a_sym*c_sym + a_sym*d_sym + b_sym*c_sym + b_sym*d_sym + c_sym*d_sym
e2_expr = expand(e2_expr)
e3_expr = a_sym*b_sym*c_sym + a_sym*b_sym*d_sym + a_sym*c_sym*d_sym + b_sym*c_sym*d_sym
e3_expr = expand(e3_expr)
e4_expr = a_sym*b_sym*c_sym*d_sym
e4_expr = expand(e4_expr)

print(f"  e2 = {e2_expr}")
print(f"  e3 = {e3_expr}")
print(f"  e4 = {e4_expr}")

# The coefficients of p(x) are:
# p(x) = x^4 + 0*x^3 + e2*x^2 + (-e3)*x + e4
# Actually: p(x) = x^4 - e1*x^3 + e2*x^2 - e3*x + e4
# With e1=0: p(x) = x^4 + e2*x^2 - e3*x + e4
# So a_2 = e2, a_3 = -e3, a_4 = e4

# We need to express N4_expanded as a polynomial in e2, e3, e4.
# Since N4 is symmetric in all 4 roots, it should be expressible in terms of e_k.
# But our variable d is constrained as d = -(a+b+c), making N4 a polynomial
# in (a, b, c) that's symmetric under S_4.

# Let me check: is N4 actually symmetric?
# Phi_4 is symmetric (sum over all roots), disc is symmetric, so N_4 is symmetric. Yes.

# For a symmetric polynomial in 4 variables with the constraint r1+r2+r3+r4=0,
# we can express it in terms of e2, e3, e4 (power sums or elementary symmetric, with e1=0).

# Strategy: numerical substitution to express N4 in terms of e2, e3, e4.
# Use the monomials with weight 10 (from Part A2).

# Let me use specific numerical values and solve
print("\n  Identifying N_4 as polynomial in (e2, e3, e4)...")

# Sample specific (a,b,c) values and compute N_4, e2, e3, e4
sample_data = []
for _ in range(30):
    av = np.random.randn() * 2
    bv = np.random.randn() * 2
    cv = np.random.randn() * 2
    dv = -(av + bv + cv)

    roots_v = sorted([av, bv, cv, dv])
    if min(np.diff(roots_v)) < 0.3:
        continue

    e2v = sum(roots_v[i]*roots_v[j] for i in range(4) for j in range(i+1,4))
    e3v = sum(roots_v[i]*roots_v[j]*roots_v[k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4))
    e4v = roots_v[0]*roots_v[1]*roots_v[2]*roots_v[3]

    phi = phi_n(roots_v)
    d_val = disc_from_roots(roots_v)
    N4v = phi * d_val

    sample_data.append((e2v, e3v, e4v, N4v))

print(f"  Got {len(sample_data)} samples")

# Solve for N_4 coefficients in terms of the weight-10 monomials
A_exact = np.zeros((len(sample_data), len(monomials_10)))
b_exact = np.zeros(len(sample_data))

for idx, (e2v, e3v, e4v, N4v) in enumerate(sample_data):
    for j, (a_exp, b_exp, c_exp) in enumerate(monomials_10):
        # Here the monomials are in (a2, a3, a4) = (e2, -e3, e4)
        # Wait, we need to be careful. The monomials should be in e2, e3, e4
        # where e2, e3, e4 are the ELEMENTARY SYMMETRIC polynomials.
        # But a_2 = e2, a_3 = -e3, a_4 = e4
        # The N_4 should be in terms of the COEFFICIENTS (a_2, a_3, a_4)
        # which equal (e2, -e3, e4).
        # But for the monomial basis, let's use e2, e3, e4 directly.
        A_exact[idx, j] = e2v**a_exp * e3v**b_exp * e4v**c_exp
    b_exact[idx] = N4v

coeffs_exact, _, _, _ = np.linalg.lstsq(A_exact, b_exact, rcond=None)

# Verify
predicted_exact = A_exact @ coeffs_exact
max_rel = np.max(np.abs(predicted_exact - b_exact) / (np.abs(b_exact) + 1e-30))
print(f"  Fit quality (max rel error): {max_rel:.2e}")

# Print N_4 polynomial
print(f"\n  N_4 (in terms of e2, e3, e4) =")
N4_terms = []
for j, (a_exp, b_exp, c_exp) in enumerate(monomials_10):
    coeff = coeffs_exact[j]
    if abs(coeff) < 1e-4:
        continue

    # Identify as rational
    identified = False
    for denom_try in range(1, 500):
        val = coeff * denom_try
        if abs(val - round(val)) < 0.005:
            num = round(val)
            g = gcd(abs(num), denom_try)
            num //= g
            den = denom_try // g
            if den == 1:
                coeff_str = f"{num}"
            else:
                coeff_str = f"{num}/{den}"
            identified = True
            N4_terms.append((num, den, a_exp, b_exp, c_exp))
            break
    if not identified:
        coeff_str = f"{coeff:.8f}"
        N4_terms.append((coeff, 1, a_exp, b_exp, c_exp))

    mon = ""
    if a_exp > 0:
        mon += f"e2^{a_exp}" if a_exp > 1 else "e2"
    if b_exp > 0:
        mon += ("*" if mon else "") + (f"e3^{b_exp}" if b_exp > 1 else "e3")
    if c_exp > 0:
        mon += ("*" if mon else "") + (f"e4^{c_exp}" if c_exp > 1 else "e4")
    if not mon:
        mon = "1"

    print(f"    {coeff_str:>10s} * {mon}")

# Now substitute kappa relations to get 1/Phi_4 = disc/N_4 in kappa variables
print("\n  Now computing 1/Phi_4 = disc/N_4 in terms of kappa...")

# Use sympy for exact computation
# disc in terms of e2, e3, e4:
e2_sym, e3_sym, e4_sym = symbols('e2 e3 e4')

# Compute disc symbolically
disc_from_coeffs = Poly(x_s**4 + e2_sym*x_s**2 - e3_sym*x_s + e4_sym, x_s).discriminant()
disc_from_coeffs = expand(disc_from_coeffs)
print(f"\n  disc (in e2, e3, e4) = {disc_from_coeffs}")

# Note: p(x) = x^4 + e2*x^2 - e3*x + e4 (coefficient of x is -e3, not +e3)
# The cumulant relations use COEFFICIENTS a_2, a_3, a_4:
# a_2 = e2 (coeff of x^2)
# a_3 = -e3 (coeff of x)
# a_4 = e4 (constant)

# Substitution: e2 = -3*k2/2, -e3 = -k3/2 => e3 = k3/2, e4 = 3*k2^2/16 - 3*k4/32
e2_to_k = -Rational(3,2)*k2_s
e3_to_k = Rational(1,2)*k3_s   # e3 = k3/2 (since a_3 = -e3 = -k3/2)
e4_to_k = Rational(3,16)*k2_s**2 - Rational(3,32)*k4_s

print(f"\n  Substitution: e2 = {e2_to_k}, e3 = {e3_to_k}, e4 = {e4_to_k}")

# Substitute into disc
disc_kappa2 = disc_from_coeffs.subs([(e2_sym, e2_to_k), (e3_sym, e3_to_k), (e4_sym, e4_to_k)])
disc_kappa2 = expand(disc_kappa2)
print(f"\n  disc (in kappas) = {disc_kappa2}")

# For N_4, use the identified rational coefficients
# Build N_4 symbolically in terms of e2, e3, e4
N4_sym_poly = sp.Integer(0)
for j, (a_exp, b_exp, c_exp) in enumerate(monomials_10):
    coeff = coeffs_exact[j]
    if abs(coeff) < 1e-4:
        continue
    # Try to get exact rational
    for denom_try in range(1, 500):
        val = coeff * denom_try
        if abs(val - round(val)) < 0.005:
            num = round(val)
            g = gcd(abs(num), denom_try)
            exact_coeff = Rational(num//g, denom_try//g)
            break
    else:
        exact_coeff = Rational(round(coeff * 1000), 1000)  # fallback

    N4_sym_poly += exact_coeff * e2_sym**a_exp * e3_sym**b_exp * e4_sym**c_exp

print(f"\n  N_4 (symbolic, in e2,e3,e4) = {N4_sym_poly}")

# Substitute kappa relations into N_4
N4_kappa = N4_sym_poly.subs([(e2_sym, e2_to_k), (e3_sym, e3_to_k), (e4_sym, e4_to_k)])
N4_kappa = expand(N4_kappa)
print(f"\n  N_4 (in kappas) = {N4_kappa}")

# Now 1/Phi_4 = disc / N_4
# Compute disc/N_4 and check if it equals C_4*k2 + R_4
print("\n  Computing 1/Phi_4 = disc/N_4...")

inv_phi4 = cancel(disc_kappa2 / N4_kappa)
print(f"\n  1/Phi_4 = {inv_phi4}")

# Partial fraction decomposition w.r.t. k2
# Or just evaluate at k3=k4=0 to get C_4
inv_phi4_at_origin = inv_phi4.subs([(k3_s, 0), (k4_s, 0)])
inv_phi4_at_origin = cancel(inv_phi4_at_origin)
print(f"\n  1/Phi_4 at k3=k4=0: {inv_phi4_at_origin}")
print(f"  This should be C_4 * k2")

# Check linearity: is it C_4 * k2?
C4_extracted = cancel(inv_phi4_at_origin / k2_s)
print(f"  C_4 = {C4_extracted}")

# Compare with claimed values
print(f"\n  CLAIMED C_4 = 1/12 = {Rational(1,12)}")
print(f"  COMPUTED C_4 = {C4_extracted}")
print(f"  CONJECTURED 2/(4*3) = 1/6 = {Rational(1,6)}")


# ============================================================
# PART D: ADVERSARIAL SUPERADDITIVITY TESTING
# ============================================================
print("\n" + "=" * 70)
print("PART D: ADVERSARIAL SUPERADDITIVITY TESTING FOR R_4")
print("=" * 70)

# First determine the correct formula
# We'll use the DIRECT computation (1/Phi_4 from roots) and check superadditivity
# of 1/Phi_n, which should be equivalent to superadditivity of R_n when the
# linear part is additive.

# Direct superadditivity test of 1/Phi_4 under MSS convolution
np.random.seed(42)
n_val = 4
n_trials = 10000
violations = 0
min_gap = float('inf')
worst_case = None

print(f"\n  Running {n_trials} random trials for 1/Phi_4 superadditivity...")

for trial in range(n_trials):
    # Generate random centered quartics with distinct real roots
    roots_p = np.sort(np.random.randn(n_val) * 2)
    roots_p -= np.mean(roots_p)
    roots_q = np.sort(np.random.randn(n_val) * 2)
    roots_q -= np.mean(roots_q)

    if min(np.diff(roots_p)) < 0.1 or min(np.diff(roots_q)) < 0.1:
        continue

    try:
        roots_r = mss_convolve(roots_p, roots_q)

        if np.max(np.abs(np.imag(np.roots(np.poly(roots_r))))) > 0.01:
            continue

        phi_p = phi_n(roots_p)
        phi_q = phi_n(roots_q)
        phi_r = phi_n(roots_r)

        gap = 1/phi_r - (1/phi_p + 1/phi_q)

        if gap < -1e-8:
            violations += 1
            print(f"  VIOLATION #{violations}: gap={gap:.2e}")
            print(f"    p_roots={roots_p}, q_roots={roots_q}")

        if gap < min_gap:
            min_gap = gap
            worst_case = (roots_p.copy(), roots_q.copy(), gap)
    except:
        pass

print(f"\n  Results: {violations} violations in ~{n_trials} trials")
print(f"  Minimum gap: {min_gap:.6e}")
if worst_case:
    print(f"  Worst case: p={worst_case[0]}, q={worst_case[1]}")

# Adversarial cases
print("\n  --- Adversarial: extreme kappa ratios ---")
for ratio in [100, 1000, 0.01, 0.001]:
    violations_adv = 0
    for _ in range(500):
        scale_p = np.sqrt(ratio)
        scale_q = 1.0 / np.sqrt(ratio)
        roots_p = np.sort(np.random.randn(4)) * scale_p
        roots_p -= np.mean(roots_p)
        roots_q = np.sort(np.random.randn(4)) * scale_q
        roots_q -= np.mean(roots_q)

        if min(np.diff(roots_p)) < 0.01 or min(np.diff(roots_q)) < 0.01:
            continue

        try:
            roots_r = mss_convolve(roots_p, roots_q)
            phi_p = phi_n(roots_p)
            phi_q = phi_n(roots_q)
            phi_r = phi_n(roots_r)
            gap = 1/phi_r - (1/phi_p + 1/phi_q)
            if gap < -1e-8:
                violations_adv += 1
        except:
            pass
    print(f"  Ratio {ratio}: {violations_adv} violations in 500 trials")

print("\n  --- Adversarial: large k3/k4 relative to k2 ---")
violations_extreme = 0
for _ in range(500):
    # Create roots with small k2 but large k3/k4
    # This means roots are clustered but have large skewness/kurtosis
    roots_p = np.array([-0.1, 0.0, 0.1, 3.0])
    roots_p -= np.mean(roots_p)
    roots_p += np.random.randn(4) * 0.01
    roots_p -= np.mean(roots_p)
    roots_p = np.sort(roots_p)

    roots_q = np.sort(np.random.randn(4) * 2)
    roots_q -= np.mean(roots_q)

    if min(np.diff(roots_p)) < 0.01 or min(np.diff(roots_q)) < 0.01:
        continue

    try:
        roots_r = mss_convolve(roots_p, roots_q)
        phi_p = phi_n(roots_p)
        phi_q = phi_n(roots_q)
        phi_r = phi_n(roots_r)
        gap = 1/phi_r - (1/phi_p + 1/phi_q)
        if gap < -1e-8:
            violations_extreme += 1
    except:
        pass
print(f"  Large k3/k4: {violations_extreme} violations in 500 trials")

print("\n  --- Adversarial: sign variations ---")
for signs in [(1,1), (1,-1), (-1,1), (-1,-1)]:
    violations_sign = 0
    for _ in range(500):
        roots_p = np.sort(np.random.randn(4) * 2)
        roots_p -= np.mean(roots_p)

        # Flip sign of k3 or k4 by reflecting roots
        if signs[0] < 0:
            roots_p = -roots_p[::-1]  # negate and reverse to flip k3 sign

        roots_q = np.sort(np.random.randn(4) * 2)
        roots_q -= np.mean(roots_q)
        if signs[1] < 0:
            roots_q = -roots_q[::-1]

        if min(np.diff(roots_p)) < 0.1 or min(np.diff(roots_q)) < 0.1:
            continue

        try:
            roots_r = mss_convolve(roots_p, roots_q)
            phi_p = phi_n(roots_p)
            phi_q = phi_n(roots_q)
            phi_r = phi_n(roots_r)
            gap = 1/phi_r - (1/phi_p + 1/phi_q)
            if gap < -1e-8:
                violations_sign += 1
        except:
            pass
    print(f"  Signs {signs}: {violations_sign} violations in 500 trials")


# ============================================================
# PART E: DOMAIN ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("PART E: DOMAIN ANALYSIS FOR R_4")
print("=" * 70)

# What is the domain of valid (k2, k3, k4)?
# A centered quartic x^4 + a2*x^2 + a3*x + a4 has 4 distinct real roots iff disc > 0.
# In cumulant space: disc(k2, k3, k4) > 0.

# Check: is the denominator of R_4 always positive on this domain?
print("\n  Testing: is denominator of R_4 always positive when disc > 0?")

np.random.seed(123456)
den_negative = 0
n_valid = 0

for _ in range(10000):
    roots = np.sort(np.random.randn(4) * 2)
    roots -= np.mean(roots)
    if min(np.diff(roots)) < 0.05:
        continue

    kappa = finite_free_cumulants(roots, 4)
    k2, k3, k4 = kappa[2], kappa[3], kappa[4]

    if k2 <= 0:
        continue

    n_valid += 1
    den = 16*k2**5 - 8*k2**2*k3**2 - k2*k4**2 + 2*k3**2*k4

    if den <= 0:
        den_negative += 1
        if den_negative <= 5:
            print(f"  NEGATIVE DEN: k2={k2:.6f}, k3={k3:.6f}, k4={k4:.6f}, den={den:.6e}")
            print(f"    roots={roots}")

print(f"\n  Result: {den_negative} negative denominators out of {n_valid} valid samples")
if den_negative == 0:
    print(f"  Denominator of R_4 appears to be ALWAYS POSITIVE on the valid domain")
else:
    print(f"  FINDING: Denominator of R_4 can be NEGATIVE!")

# Check sign of R_4 itself
print("\n  Testing: sign of R_4 on the domain")
R4_positive = 0
R4_negative = 0
R4_zero = 0

for _ in range(10000):
    roots = np.sort(np.random.randn(4) * 2)
    roots -= np.mean(roots)
    if min(np.diff(roots)) < 0.1:
        continue

    kappa = finite_free_cumulants(roots, 4)
    k2, k3, k4 = kappa[2], kappa[3], kappa[4]

    if k2 <= 0:
        continue

    try:
        r4 = R4_claimed(k2, k3, k4)
        if r4 > 1e-10:
            R4_positive += 1
        elif r4 < -1e-10:
            R4_negative += 1
        else:
            R4_zero += 1
    except:
        pass

print(f"  R_4 > 0: {R4_positive}, R_4 < 0: {R4_negative}, R_4 ≈ 0: {R4_zero}")

# Boundary behavior
print("\n  Boundary: what happens as k3 -> 0 and k4 -> 0?")
for k2_val in [1.0, 2.0, 5.0]:
    for eps in [0.1, 0.01, 0.001, 0.0001]:
        try:
            r4 = R4_claimed(k2_val, eps, eps)
            print(f"  k2={k2_val}, k3=k4={eps}: R_4 = {r4:.8e}")
        except:
            print(f"  k2={k2_val}, k3=k4={eps}: UNDEFINED")


# ============================================================
# PART F: VERIFY R_4 SUPERADDITIVITY DIRECTLY
# ============================================================
print("\n" + "=" * 70)
print("PART F: DIRECT R_4 SUPERADDITIVITY TEST")
print("=" * 70)

# R_4 superadditivity: R_4(k_p + k_q) >= R_4(k_p) + R_4(k_q)
# where k_p = (k2_p, k3_p, k4_p) are cumulants of p.

# We need VALID cumulant vectors. Not all (k2, k3, k4) correspond to
# polynomials with simple real roots. We sample from actual polynomials.

np.random.seed(98765)
violations_R4 = 0
min_gap_R4 = float('inf')
n_tested_R4 = 0

for trial in range(10000):
    roots_p = np.sort(np.random.randn(4) * 2)
    roots_p -= np.mean(roots_p)
    roots_q = np.sort(np.random.randn(4) * 2)
    roots_q -= np.mean(roots_q)

    if min(np.diff(roots_p)) < 0.1 or min(np.diff(roots_q)) < 0.1:
        continue

    kp = finite_free_cumulants(roots_p, 4)
    kq = finite_free_cumulants(roots_q, 4)

    k2_p, k3_p, k4_p = kp[2], kp[3], kp[4]
    k2_q, k3_q, k4_q = kq[2], kq[3], kq[4]

    # Sum cumulants
    k2_r = k2_p + k2_q
    k3_r = k3_p + k3_q
    k4_r = k4_p + k4_q

    try:
        R4_sum = R4_claimed(k2_r, k3_r, k4_r)
        R4_p = R4_claimed(k2_p, k3_p, k4_p)
        R4_q = R4_claimed(k2_q, k3_q, k4_q)

        gap = R4_sum - (R4_p + R4_q)
        n_tested_R4 += 1

        if gap < -1e-8:
            violations_R4 += 1
            if violations_R4 <= 5:
                print(f"  VIOLATION: gap={gap:.2e}")
                print(f"    k_p=({k2_p:.4f}, {k3_p:.4f}, {k4_p:.4f})")
                print(f"    k_q=({k2_q:.4f}, {k3_q:.4f}, {k4_q:.4f})")

        if gap < min_gap_R4:
            min_gap_R4 = gap
    except:
        pass

print(f"\n  Direct R_4 superadditivity: {violations_R4} violations in {n_tested_R4} valid tests")
print(f"  Minimum gap: {min_gap_R4:.6e}")

# IMPORTANT: The cumulant sum (k_p + k_q) may not correspond to any polynomial with
# simple real roots! We should check this.
print("\n  Checking: do summed cumulants always give valid polynomials?")
n_invalid_sum = 0
for trial in range(1000):
    roots_p = np.sort(np.random.randn(4) * 2)
    roots_p -= np.mean(roots_p)
    roots_q = np.sort(np.random.randn(4) * 2)
    roots_q -= np.mean(roots_q)

    if min(np.diff(roots_p)) < 0.1 or min(np.diff(roots_q)) < 0.1:
        continue

    kp = finite_free_cumulants(roots_p, 4)
    kq = finite_free_cumulants(roots_q, 4)

    k2_r = kp[2] + kq[2]
    k3_r = kp[3] + kq[3]
    k4_r = kp[4] + kq[4]

    # Reconstruct polynomial from cumulants
    a2_r = -3*k2_r/2
    a3_r = -k3_r/2
    a4_r = 3*k2_r**2/16 - 3*k4_r/32

    coeffs_r = [1.0, 0.0, a2_r, a3_r, a4_r]
    roots_r = np.roots(coeffs_r)

    if np.max(np.abs(np.imag(roots_r))) > 1e-6:
        n_invalid_sum += 1
    elif min(np.diff(np.sort(np.real(roots_r)))) < 1e-8:
        n_invalid_sum += 1

print(f"  Invalid summed cumulant polynomials: {n_invalid_sum} out of 1000")
print(f"  NOTE: Cumulant additivity under MSS ensures the SUM is always valid")
print(f"  (because MSS convolution preserves real-rootedness by ADMITTED A)")

print("\n" + "=" * 70)
print("SUMMARY OF ALL FINDINGS")
print("=" * 70)
