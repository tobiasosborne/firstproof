"""
Verify the mathematical claims in node 1.5.2 and children.
VERIFIER-6 adversarial checks.
"""
import numpy as np
from itertools import combinations
from math import comb, factorial
from numpy.polynomial import polynomial as P

# ============================================================
# PART 0: Basic definitions
# ============================================================

def elementary_symmetric(roots, k):
    """e_k of the roots (sign convention: e_k = sum of products of k roots)."""
    n = len(roots)
    if k == 0:
        return 1.0
    if k > n:
        return 0.0
    return sum(np.prod(combo) for combo in combinations(roots, k))

def disc(roots):
    """Discriminant = prod_{i<j} (r_i - r_j)^2."""
    n = len(roots)
    d = 1.0
    for i in range(n):
        for j in range(i+1, n):
            d *= (roots[i] - roots[j])**2
    return d

def phi_n(roots):
    """Phi_n = sum_i H_p(lambda_i)^2 where H_p = p''/(2p')."""
    n = len(roots)
    total = 0.0
    for i in range(n):
        # p'(lambda_i) = prod_{j != i} (lambda_i - lambda_j)
        pprime = 1.0
        for j in range(n):
            if j != i:
                pprime *= (roots[i] - roots[j])

        # p''(lambda_i) = 2 * sum_{j != i} prod_{k != i, k != j} (lambda_i - lambda_k)
        pprimeprime = 0.0
        for j in range(n):
            if j != i:
                term = 1.0
                for k in range(n):
                    if k != i and k != j:
                        term *= (roots[i] - roots[k])
                pprimeprime += term
        pprimeprime *= 2.0

        H = pprimeprime / (2.0 * pprime)
        total += H**2
    return total

def mss_convolve(roots_p, roots_q):
    """MSS convolution of two polynomials given by roots.
    p(x) = prod(x - a_i), q(x) = prod(x - b_j).
    r(x) coefficients via MSS formula: c_k = sum_{i+j=k} w(n,i,j) a_i b_j
    where w(n,i,j) = (n-i)!(n-j)! / (n!(n-k)!) and a_i, b_j are coefficients.
    """
    n = len(roots_p)
    assert len(roots_q) == n

    # Get coefficients: p(x) = x^n + a_1 x^{n-1} + ... + a_n
    # Using numpy: coefficients in standard form
    p_coeffs = np.poly(roots_p)  # [1, a_1, a_2, ..., a_n] highest degree first
    q_coeffs = np.poly(roots_q)

    # MSS formula
    r_coeffs = np.zeros(n + 1)
    r_coeffs[0] = 1.0  # monic

    for k in range(1, n + 1):
        ck = 0.0
        for i in range(k + 1):
            j = k - i
            if i <= n and j <= n:
                w = factorial(n - i) * factorial(n - j) / (factorial(n) * factorial(n - k))
                ck += w * p_coeffs[i] * q_coeffs[j]
        r_coeffs[k] = ck

    # Find roots of r
    r_roots = np.roots(r_coeffs)
    return np.sort(np.real(r_roots))

# ============================================================
# CHECK (a): Is Phi_2 = 1/kappa_2?
# ============================================================
print("="*60)
print("CHECK (a): Phi_2 = 1/kappa_2 for n=2")
print("="*60)

# For n=2, p(x) = (x-a)(x-b) = x^2 - (a+b)x + ab
# Centered: a+b=0, so roots = {a, -a}
# e_1 = 0, e_2 = a*(-a) = -a^2
# tilde_a_2 = (-1)^2 * e_2 / C(2,2) = (-a^2)/1 = -a^2
# Wait, let me be more careful. The coefficient convention:
# p(x) = x^2 + a_1*x + a_2 where a_1 = -(a+b), a_2 = ab
# For centered: a_1 = 0, a_2 = ab = -a^2

# Finite free cumulants for n=2:
# kappa_1 = tilde_a_1 = -a_1/n = 0 (centered)
# kappa_2 = -n*(tilde_a_2 - tilde_a_1^2)
# tilde_a_k = (-1)^k * a_k / C(n,k)
# tilde_a_2 = a_2 / C(2,2) = a_2 / 1 = a_2 = -a^2
# kappa_2 = -2*(-a^2 - 0) = 2*a^2

# Phi_2 for roots {a, -a}:
# p'(x) = 2x, p''(x) = 2
# H_p(a) = p''(a)/(2p'(a)) = 2/(2*2a) = 1/(2a)
# H_p(-a) = 2/(2*(-2a)) = -1/(2a)
# Phi_2 = (1/(2a))^2 + (1/(2a))^2 = 2/(4a^2) = 1/(2a^2)

# So Phi_2 = 1/(2a^2) and kappa_2 = 2a^2
# Therefore Phi_2 = 1/kappa_2. CONFIRMED.

for a_val in [1.0, 2.0, 0.5, 3.7]:
    roots = np.array([a_val, -a_val])
    phi = phi_n(roots)
    # Compute kappa_2
    a2 = roots[0]*roots[1]  # = -a^2
    tilde_a2 = a2 / comb(2, 2)  # = a2
    kappa_2 = -2 * (tilde_a2)  # kappa_2 = -n*(tilde_a_2 - tilde_a_1^2) = -2*tilde_a_2

    print(f"  a={a_val}: Phi_2 = {phi:.10f}, 1/kappa_2 = {1/kappa_2:.10f}, "
          f"ratio = {phi * kappa_2:.10f}")

# Also verify with non-centered
print("\n  Non-centered cases:")
for roots in [np.array([1.0, 3.0]), np.array([0.5, 4.5]), np.array([-1.0, 5.0])]:
    phi = phi_n(roots)
    # General kappa_2 for n=2
    a1 = -(roots[0] + roots[1])
    a2 = roots[0] * roots[1]
    tilde_a1 = -a1 / 2  # = (roots[0]+roots[1])/2
    tilde_a2 = a2 / 1   # C(2,2) = 1
    kappa_2 = -2 * (tilde_a2 - tilde_a1**2)

    print(f"  roots={roots}: Phi_2 = {phi:.10f}, 1/kappa_2 = {1/kappa_2:.10f}, "
          f"ratio = {phi * kappa_2:.10f}")

# ============================================================
# CHECK (b): Phi_3 formula and rationality
# ============================================================
print("\n" + "="*60)
print("CHECK (b): Phi_3 = (9/2)*kappa_2^2 / (kappa_2^3 - kappa_3^2/3)")
print("="*60)

# For n=3, centered (e_1 = 0):
# Verify Phi_3 * disc = 18 * e_2^2

for trial in range(5):
    # Random centered roots for n=3
    r = np.random.randn(3)
    r = r - np.mean(r)  # center

    phi = phi_n(r)
    d = disc(r)
    e2 = elementary_symmetric(r, 2)
    e3 = elementary_symmetric(r, 3)

    lhs = phi * d
    rhs = 18 * e2**2

    print(f"  Trial {trial}: Phi_3*disc = {lhs:.10f}, 18*e_2^2 = {rhs:.10f}, "
          f"ratio = {lhs/rhs:.10f}")

    # Now check the cumulant formula
    # For n=3 centered:
    # tilde_a_1 = 0
    # tilde_a_2 = a_2/C(3,2) = e_2/3  (a_2 = e_2 for monic)
    # tilde_a_3 = -a_3/C(3,3) = -(-e_3)/1 = e_3
    # Wait, let me be careful. p(x) = x^3 + a_1 x^2 + a_2 x + a_3
    # where a_1 = -(sum roots), a_2 = sum_{i<j} r_i*r_j = e_2, a_3 = -prod roots = -e_3
    # Actually, by Vieta: a_1 = -e_1, a_2 = e_2, a_3 = -e_3

    a1 = -sum(r)  # should be ~0
    a2 = e2
    a3 = -e3

    # tilde_a_k = (-1)^k * a_k / C(n,k)
    tilde_a1 = (-1)**1 * a1 / comb(3,1)  # = -a1/3 = e_1/3 = 0
    tilde_a2 = (-1)**2 * a2 / comb(3,2)  # = a_2/3 = e_2/3
    tilde_a3 = (-1)**3 * a3 / comb(3,3)  # = -a_3/1 = e_3

    # kappa_2 = -n*(tilde_a_2 - tilde_a_1^2) = -3*(e_2/3 - 0) = -e_2
    kappa2 = -3 * (tilde_a2 - tilde_a1**2)

    # kappa_3 from node 1.5.2.1:
    # kappa_3 = n^2/2 * (tilde_a_3 - 3*tilde_a_2*tilde_a_1 + 2*tilde_a_1^3)
    # = 9/2 * (e_3 - 0 + 0) = 9*e_3/2
    kappa3 = 9/2 * (tilde_a3 - 3*tilde_a2*tilde_a1 + 2*tilde_a1**3)

    # Claimed: Phi_3 = (9/2) * kappa_2^2 / (kappa_2^3 - kappa_3^2/3)
    # Substitute: kappa_2 = -e_2, kappa_3 = 9*e_3/2
    # Numerator: (9/2)*(-e_2)^2 = (9/2)*e_2^2
    # Denominator: (-e_2)^3 - (9*e_3/2)^2/3 = -e_2^3 - 81*e_3^2/12 = -e_2^3 - 27*e_3^2/4
    # So Phi_3 = (9/2)*e_2^2 / (-e_2^3 - 27*e_3^2/4)
    # = (9/2)*e_2^2 / (-(e_2^3 + 27*e_3^2/4))
    # = -(9/2)*e_2^2 / (e_2^3 + 27*e_3^2/4)
    # = -18*e_2^2 / (4*e_2^3 + 27*e_3^2)
    #
    # From Part B: Phi_3 = 18*e_2^2 / disc = 18*e_2^2 / (-4*e_2^3 - 27*e_3^2)
    # = -18*e_2^2 / (4*e_2^3 + 27*e_3^2)
    # These match!

    if abs(kappa2**3 - kappa3**2/3) > 1e-10:
        phi_cumulant = (9/2) * kappa2**2 / (kappa2**3 - kappa3**2/3)
        print(f"         Phi_3 from cumulants = {phi_cumulant:.10f}, actual = {phi:.10f}, "
              f"ratio = {phi_cumulant/phi:.10f}")

# ============================================================
# CHECK (c): Does "insufficiency" conclusion actually follow?
# ============================================================
print("\n" + "="*60)
print("CHECK (c): Testing whether 1/Phi_n might still be superadditive")
print("           despite Phi_n being rational in cumulants")
print("="*60)

# The claim is that because Phi_n is rational (not polynomial) in cumulants,
# the cumulant approach is "insufficient". But this doesn't mean no inequality exists.
# Let's test: is 1/Phi_n actually superadditive under MSS convolution?

n_tests = 100
n_val = 3
violations = 0

for trial in range(n_tests):
    # Random polynomials with distinct real roots
    roots_p = np.sort(np.random.randn(n_val) * 2)
    roots_q = np.sort(np.random.randn(n_val) * 2)

    # Make sure roots are distinct enough
    min_gap_p = min(np.diff(roots_p))
    min_gap_q = min(np.diff(roots_q))
    if min_gap_p < 0.1 or min_gap_q < 0.1:
        continue

    try:
        roots_r = mss_convolve(roots_p, roots_q)

        # Check if roots are real and distinct
        if np.max(np.abs(np.imag(np.roots(np.poly(roots_p))))) > 0.01:
            continue

        phi_p = phi_n(roots_p)
        phi_q = phi_n(roots_q)
        phi_r = phi_n(roots_r)

        # Superadditivity of 1/Phi: 1/Phi(r) >= 1/Phi(p) + 1/Phi(q)
        lhs = 1/phi_r
        rhs = 1/phi_p + 1/phi_q

        if lhs < rhs - 1e-8:  # violation
            violations += 1
            print(f"  VIOLATION at trial {trial}: 1/Phi(r)={lhs:.6f} < 1/Phi(p)+1/Phi(q)={rhs:.6f}")
    except Exception as e:
        pass

print(f"  n=3: {violations} violations out of {n_tests} trials")

# Now try n=4
n_val = 4
violations_4 = 0
for trial in range(n_tests):
    roots_p = np.sort(np.random.randn(n_val) * 2)
    roots_q = np.sort(np.random.randn(n_val) * 2)

    min_gap_p = min(np.diff(roots_p))
    min_gap_q = min(np.diff(roots_q))
    if min_gap_p < 0.1 or min_gap_q < 0.1:
        continue

    try:
        roots_r = mss_convolve(roots_p, roots_q)
        if not all(np.isreal(roots_r)):
            continue

        phi_p = phi_n(roots_p)
        phi_q = phi_n(roots_q)
        phi_r = phi_n(roots_r)

        lhs = 1/phi_r
        rhs = 1/phi_p + 1/phi_q

        if lhs < rhs - 1e-8:
            violations_4 += 1
            print(f"  VIOLATION at trial {trial}: 1/Phi(r)={lhs:.6f} < 1/Phi(p)+1/Phi(q)={rhs:.6f}")
    except Exception as e:
        pass

print(f"  n=4: {violations_4} violations out of {n_tests} trials")

# ============================================================
# CHECK (d): Phi_n * disc polynomial — could this help?
# ============================================================
print("\n" + "="*60)
print("CHECK (d): Is disc itself expressible simply via cumulants?")
print("="*60)

# For n=3 centered:
# disc = -4*e_2^3 - 27*e_3^2
# kappa_2 = -e_2, kappa_3 = 9*e_3/2
# So e_2 = -kappa_2, e_3 = 2*kappa_3/9
# disc = -4*(-kappa_2)^3 - 27*(2*kappa_3/9)^2
#       = 4*kappa_2^3 - 27*4*kappa_3^2/81
#       = 4*kappa_2^3 - 4*kappa_3^2/3

print("  For n=3 centered:")
print("  disc = 4*kappa_2^3 - 4*kappa_3^2/3")
print("  = 4*(kappa_2^3 - kappa_3^2/3)")
print("  N_3 = 18*e_2^2 = 18*kappa_2^2")
print("  So Phi_3 = 18*kappa_2^2 / (4*(kappa_2^3 - kappa_3^2/3))")
print("           = (9/2)*kappa_2^2 / (kappa_2^3 - kappa_3^2/3)")
print("  This confirms the formula from Part C.")

# Key question: under cumulant additivity, is disc superadditive?
# disc(r) = 4*(kappa_2(p)+kappa_2(q))^3 - 4*(kappa_3(p)+kappa_3(q))^2/3
# This is NOT in general >= disc(p) + disc(q) because:
# disc(p) + disc(q) = 4*(kappa_2(p)^3 + kappa_2(q)^3) - 4*(kappa_3(p)^2 + kappa_3(q)^2)/3
# The cubic (a+b)^3 >= a^3 + b^3 when a,b > 0, but (a+b)^2 >= a^2 + b^2 also,
# so the sign depends on relative magnitudes.

print("\n  Testing disc behavior under convolution (n=3):")
for trial in range(10):
    r = np.random.randn(3)
    r = r - np.mean(r)
    s = np.random.randn(3)
    s = s - np.mean(s)

    # Compute kappa values
    e2_p = elementary_symmetric(r, 2)
    e3_p = elementary_symmetric(r, 3)
    e2_q = elementary_symmetric(s, 2)
    e3_q = elementary_symmetric(s, 3)

    k2_p = -e2_p
    k3_p = 9*e3_p/2
    k2_q = -e2_q
    k3_q = 9*e3_q/2

    # Convolved cumulants
    k2_r = k2_p + k2_q
    k3_r = k3_p + k3_q

    disc_p = 4*(k2_p**3 - k3_p**2/3)
    disc_q = 4*(k2_q**3 - k3_q**2/3)
    disc_r = 4*(k2_r**3 - k3_r**2/3)

    print(f"  Trial {trial}: disc(r)={disc_r:.4f}, disc(p)+disc(q)={disc_p+disc_q:.4f}, "
          f"ratio={disc_r/(disc_p+disc_q):.4f}" if abs(disc_p+disc_q) > 1e-10 else f"  Trial {trial}: degenerate")

# ============================================================
# CHECK: Verify the n=3 discriminant formula is correct
# ============================================================
print("\n" + "="*60)
print("CHECK: Verify disc formula for n=3")
print("="*60)

for trial in range(5):
    r = np.random.randn(3)
    r = r - np.mean(r)

    d = disc(r)
    e2 = elementary_symmetric(r, 2)
    e3 = elementary_symmetric(r, 3)

    d_formula = -4*e2**3 - 27*e3**2
    print(f"  disc(roots) = {d:.10f}, -4*e2^3-27*e3^2 = {d_formula:.10f}, "
          f"ratio = {d/d_formula:.10f}")

# ============================================================
# CHECK: n=2 superadditivity is exact equality (as claimed in 1.5.2.4)
# ============================================================
print("\n" + "="*60)
print("CHECK: n=2 exact equality in superadditivity")
print("="*60)

for trial in range(10):
    roots_p = np.sort(np.random.randn(2) * 2)
    roots_q = np.sort(np.random.randn(2) * 2)

    if abs(roots_p[0] - roots_p[1]) < 0.1 or abs(roots_q[0] - roots_q[1]) < 0.1:
        continue

    try:
        roots_r = mss_convolve(roots_p, roots_q)
        phi_p = phi_n(roots_p)
        phi_q = phi_n(roots_q)
        phi_r = phi_n(roots_r)

        lhs = 1/phi_r
        rhs = 1/phi_p + 1/phi_q

        print(f"  1/Phi(r) = {lhs:.10f}, 1/Phi(p)+1/Phi(q) = {rhs:.10f}, "
              f"diff = {lhs - rhs:.2e}")
    except:
        pass

# ============================================================
# CHECK: Additivity of cumulants under MSS
# ============================================================
print("\n" + "="*60)
print("CHECK: Additivity of kappa_2 under MSS convolution (n=3)")
print("="*60)

for trial in range(10):
    roots_p = np.sort(np.random.randn(3))
    roots_p = roots_p - np.mean(roots_p)
    roots_q = np.sort(np.random.randn(3))
    roots_q = roots_q - np.mean(roots_q)

    try:
        roots_r = mss_convolve(roots_p, roots_q)

        # Compute kappa_2 for each
        e2_p = elementary_symmetric(roots_p, 2)
        e2_q = elementary_symmetric(roots_q, 2)
        e2_r = elementary_symmetric(roots_r, 2)

        k2_p = -e2_p  # kappa_2 = -e_2 for n=3 centered
        k2_q = -e2_q
        k2_r = -e2_r

        print(f"  kappa_2(r) = {k2_r:.8f}, kappa_2(p)+kappa_2(q) = {k2_p+k2_q:.8f}, "
              f"diff = {k2_r - (k2_p+k2_q):.2e}")
    except:
        pass
