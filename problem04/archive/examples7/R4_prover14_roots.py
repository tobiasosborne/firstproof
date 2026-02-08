"""
PROVER-14: Root Geometry Approach to Fisher Superadditivity

STRATEGY: Work directly with roots, not cumulants.
For p with roots mu_1 < ... < mu_n:
  H_p(mu_i) = sum_{j!=i} 1/(mu_i - mu_j)
  Phi_n(p) = sum_i H_p(mu_i)^2

The conjecture: 1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q)
where r = p boxplus_n q.

STEP 1: Understand Phi_n as a function of root gaps.
"""
import numpy as np
from itertools import combinations
from numpy.polynomial import polynomial as P
import sympy as sp

print("="*70)
print("PROVER-14: ROOT GEOMETRY EXPLORATION")
print("="*70)

# =============================================================
# Part 1: Phi_n in terms of roots
# =============================================================
print("\n--- Part 1: Phi_n for small n ---\n")

def H_values(roots):
    """Compute H_p(lambda_i) for each root."""
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H

def Phi_n(roots):
    """Compute Phi_n(p) = sum_i H_p(lambda_i)^2."""
    H = H_values(roots)
    return np.sum(H**2)

def inv_Phi_n(roots):
    """Compute 1/Phi_n(p)."""
    phi = Phi_n(roots)
    if phi < 1e-15:
        return float('inf')
    return 1.0 / phi

# n=2: explicit formula
print("n=2:")
print("  roots = (a, b), gap d = b - a > 0")
print("  H_p(a) = 1/(a-b) = -1/d")
print("  H_p(b) = 1/(b-a) = 1/d")
print("  Phi_2 = 1/d^2 + 1/d^2 = 2/d^2")
print("  1/Phi_2 = d^2/2")
print()

# Verify
for d in [1.0, 2.0, 0.5]:
    roots = np.array([0, d])
    phi = Phi_n(roots)
    print(f"  d={d}: Phi_2 = {phi:.6f}, expected {2/d**2:.6f}")

# n=3: symbolic
print("\n\nn=3 (symbolic):")
a, b, c = sp.symbols('a b c', real=True)
# roots a < b < c
H_a = 1/(a-b) + 1/(a-c)
H_b = 1/(b-a) + 1/(b-c)
H_c = 1/(c-a) + 1/(c-b)
Phi_3_sym = sp.expand(H_a**2 + H_b**2 + H_c**2)
print(f"  Phi_3 = {Phi_3_sym}")

# Check sum H_i = 0
sum_H = sp.simplify(H_a + H_b + H_c)
print(f"  sum_i H(lambda_i) = {sum_H}")

# =============================================================
# Part 2: MSS convolution via roots
# =============================================================
print("\n\n--- Part 2: MSS convolution in root coordinates ---\n")

def mss_convolve_n(p_coeffs, q_coeffs, n):
    """
    MSS finite free additive convolution.
    p_coeffs and q_coeffs are coefficient vectors [a_0, a_1, ..., a_n]
    where p(x) = sum_k a_k x^{n-k}, a_0 = 1.
    """
    from math import factorial, comb
    r_coeffs = np.zeros(n+1)
    for k in range(n+1):
        ck = 0
        for i in range(k+1):
            j = k - i
            if i <= n and j <= n:
                coeff = factorial(n-i) * factorial(n-j) / (factorial(n) * factorial(n-k))
                ck += coeff * p_coeffs[i] * q_coeffs[j]
        r_coeffs[k] = ck
    return r_coeffs

def roots_to_monic_coeffs(roots):
    """Convert roots to monic polynomial coefficients [a_0, ..., a_n] where a_0=1."""
    n = len(roots)
    # p(x) = prod(x - r_i) = x^n - e_1 x^{n-1} + e_2 x^{n-2} - ...
    coeffs = [1.0]
    for k in range(1, n+1):
        # a_k = (-1)^k * e_k(roots)
        ek = sum(np.prod(list(combo)) for combo in combinations(roots, k))
        coeffs.append((-1)**k * ek)
    return np.array(coeffs)

def coeffs_to_poly(coeffs):
    """Convert [a_0, ..., a_n] to numpy polynomial (in standard form for np.roots)."""
    # coeffs[k] is coefficient of x^{n-k}
    # numpy wants highest degree first
    return coeffs

def mss_convolve_roots(roots_p, roots_q):
    """Compute roots of p boxplus_n q given roots of p and q."""
    n = len(roots_p)
    assert len(roots_q) == n
    p_coeffs = roots_to_monic_coeffs(roots_p)
    q_coeffs = roots_to_monic_coeffs(roots_q)
    r_coeffs = mss_convolve_n(p_coeffs, q_coeffs, n)
    # Find roots
    r_roots = np.roots(r_coeffs)
    # Should be all real
    r_roots = np.sort(np.real(r_roots))
    return r_roots

# Test: n=2
print("Test n=2:")
p_roots = np.array([-1.0, 1.0])
q_roots = np.array([-2.0, 2.0])
r_roots = mss_convolve_roots(p_roots, q_roots)
print(f"  p roots: {p_roots}")
print(f"  q roots: {q_roots}")
print(f"  r = p boxplus_2 q roots: {r_roots}")
# For n=2: r(x) = x^2 - (e1_p + e1_q)*x + ...
# Actually for n=2, boxplus_2 just shifts by mean difference

# Test n=3
print("\nTest n=3:")
p_roots = np.array([-2.0, 0.0, 2.0])
q_roots = np.array([-1.0, 0.5, 1.5])
r_roots = mss_convolve_roots(p_roots, q_roots)
print(f"  p roots: {p_roots}")
print(f"  q roots: {q_roots}")
print(f"  r = p boxplus_3 q roots: {r_roots}")

# =============================================================
# Part 3: Verify conjecture in root coordinates
# =============================================================
print("\n\n--- Part 3: Numerical verification in root coordinates ---\n")

np.random.seed(42)
n_trials = 10000
violations = 0
valid = 0
min_gap = float('inf')

for trial in range(n_trials):
    n = 3
    # Generate random roots with distinct values
    p_roots = np.sort(np.random.randn(n) * 2)
    q_roots = np.sort(np.random.randn(n) * 2)

    # Ensure roots are distinct (gap > 0.1)
    if np.min(np.diff(p_roots)) < 0.1 or np.min(np.diff(q_roots)) < 0.1:
        continue

    try:
        r_roots = mss_convolve_roots(p_roots, q_roots)
        if np.min(np.diff(r_roots)) < 1e-8:
            continue

        inv_phi_p = inv_Phi_n(p_roots)
        inv_phi_q = inv_Phi_n(q_roots)
        inv_phi_r = inv_Phi_n(r_roots)

        gap = inv_phi_r - inv_phi_p - inv_phi_q
        valid += 1
        min_gap = min(min_gap, gap)

        if gap < -1e-8:
            violations += 1
            print(f"  VIOLATION at trial {trial}: gap = {gap:.2e}")
            print(f"    p_roots = {p_roots}")
            print(f"    q_roots = {q_roots}")
            print(f"    r_roots = {r_roots}")
    except:
        continue

print(f"n=3: {violations} violations in {valid} valid trials")
print(f"  min gap = {min_gap:.6e}")

# Also n=4
np.random.seed(123)
violations_4 = 0
valid_4 = 0
min_gap_4 = float('inf')

for trial in range(n_trials):
    n = 4
    p_roots = np.sort(np.random.randn(n) * 2)
    q_roots = np.sort(np.random.randn(n) * 2)

    if np.min(np.diff(p_roots)) < 0.1 or np.min(np.diff(q_roots)) < 0.1:
        continue

    try:
        r_roots = mss_convolve_roots(p_roots, q_roots)
        if np.min(np.diff(r_roots)) < 1e-8:
            continue

        inv_phi_p = inv_Phi_n(p_roots)
        inv_phi_q = inv_Phi_n(q_roots)
        inv_phi_r = inv_Phi_n(r_roots)

        gap = inv_phi_r - inv_phi_p - inv_phi_q
        valid_4 += 1
        min_gap_4 = min(min_gap_4, gap)

        if gap < -1e-8:
            violations_4 += 1
    except:
        continue

print(f"\nn=4: {violations_4} violations in {valid_4} valid trials")
print(f"  min gap = {min_gap_4:.6e}")

# =============================================================
# Part 4: Study the H-vector and its properties
# =============================================================
print("\n\n--- Part 4: H-vector structure ---\n")

# Key fact: sum_i H_p(lambda_i) = 0
# So ||H||^2 = Phi_n and H lives on the hyperplane sum = 0
# The vector H has n-1 degrees of freedom

# For n=3, study the H-vector explicitly
print("n=3 H-vectors for various root configs:")
configs = [
    np.array([-2.0, 0.0, 2.0]),
    np.array([-3.0, 0.0, 1.0]),
    np.array([-1.0, 0.0, 3.0]),
    np.array([-1.0, -0.5, 0.5]),
]
for roots in configs:
    H = H_values(roots)
    phi = np.sum(H**2)
    print(f"  roots = {roots}")
    print(f"  H = [{', '.join(f'{h:.4f}' for h in H)}]")
    print(f"  sum(H) = {np.sum(H):.2e}, Phi = {phi:.4f}, 1/Phi = {1/phi:.6f}")
    print(f"  gaps = {np.diff(roots)}")
    print()

# =============================================================
# Part 5: Electrostatic energy interpretation
# =============================================================
print("\n--- Part 5: Electrostatic interpretation ---\n")

print("""
Electrostatic interpretation:
- Place unit positive charges at positions lambda_1, ..., lambda_n on a line
- The electric field at lambda_i (due to other charges) is E_i = sum_{j!=i} 1/(lambda_i - lambda_j)
- This is exactly H_p(lambda_i)!
- Phi_n = sum_i E_i^2 = sum_i H_i^2 = total "field energy"

The TOTAL electrostatic energy is:
  W = sum_{i<j} log(1/|lambda_i - lambda_j|) = -sum_{i<j} log|lambda_i - lambda_j|

Note: dW/d(lambda_i) = -sum_{j!=i} 1/(lambda_i - lambda_j) = -H_p(lambda_i)

So H = -grad(W) and Phi_n = ||grad W||^2.

The conjecture says: the field energy of the convolution is SMALLER
(or rather, 1/field_energy is LARGER) than the sum of individual 1/field_energies.
""")

# Verify the gradient relationship
for roots in configs[:2]:
    H = H_values(roots)
    n = len(roots)
    # Compute gradient of W
    eps = 1e-7
    grad_W = np.zeros(n)
    def W_func(r):
        val = 0
        for i in range(len(r)):
            for j in range(i+1, len(r)):
                val -= np.log(abs(r[i] - r[j]))
        return val

    for i in range(n):
        r_plus = roots.copy(); r_plus[i] += eps
        r_minus = roots.copy(); r_minus[i] -= eps
        grad_W[i] = (W_func(r_plus) - W_func(r_minus)) / (2*eps)

    print(f"  roots = {roots}")
    print(f"  H = {H}")
    print(f"  -grad(W) = {-grad_W}")
    print(f"  Match: {np.allclose(H, -grad_W)}")
    print()

# =============================================================
# Part 6: Root gap behavior under MSS convolution
# =============================================================
print("\n--- Part 6: Root gaps under MSS convolution ---\n")

np.random.seed(777)
print("Studying how root gaps of r relate to those of p and q:\n")

for trial in range(5):
    n = 4
    p_roots = np.sort(np.random.randn(n) * 2)
    q_roots = np.sort(np.random.randn(n) * 2)

    while np.min(np.diff(p_roots)) < 0.3:
        p_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(q_roots)) < 0.3:
        q_roots = np.sort(np.random.randn(n) * 2)

    try:
        r_roots = mss_convolve_roots(p_roots, q_roots)

        p_gaps = np.diff(p_roots)
        q_gaps = np.diff(q_roots)
        r_gaps = np.diff(r_roots)

        print(f"  Trial {trial}:")
        print(f"    p_gaps = {p_gaps}")
        print(f"    q_gaps = {q_gaps}")
        print(f"    r_gaps = {r_gaps}")
        print(f"    p_spread = {p_roots[-1]-p_roots[0]:.4f}, q_spread = {q_roots[-1]-q_roots[0]:.4f}, r_spread = {r_roots[-1]-r_roots[0]:.4f}")
        print(f"    Phi_p={Phi_n(p_roots):.4f}, Phi_q={Phi_n(q_roots):.4f}, Phi_r={Phi_n(r_roots):.4f}")
        print(f"    1/Phi_p={1/Phi_n(p_roots):.6f}, 1/Phi_q={1/Phi_n(q_roots):.6f}, 1/Phi_r={1/Phi_n(r_roots):.6f}")
        gap = 1/Phi_n(r_roots) - 1/Phi_n(p_roots) - 1/Phi_n(q_roots)
        print(f"    Gap = {gap:.6e}")
        print()
    except:
        continue

print("\nDone with Part 6.")
