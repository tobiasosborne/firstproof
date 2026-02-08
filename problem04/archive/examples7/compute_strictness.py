"""
Compute Phi_n for specific polynomials and their finite free convolution
to verify strict inequality for n >= 3.

Definitions:
- P_n: monic real-rooted polynomials of degree n with simple roots
- Phi_n(p) = sum_{i=1}^n H_p(lambda_i)^2
  where H_p(lambda_i) = sum_{j != i} 1/(lambda_i - lambda_j)
- Finite free convolution: r = p ⊞_n q with
  c_k = sum_{i+j=k} [(n-i)!(n-j)! / (n!(n-k)!)] a_i b_j
"""

import numpy as np
from math import factorial, comb
from itertools import product as iter_product

def poly_coeffs_from_roots(roots):
    """Get monic polynomial coefficients [a_0=1, a_1, ..., a_n] from roots.
    p(x) = prod(x - r_i) = x^n + a_1 x^{n-1} + ... + a_n
    """
    # numpy poly from roots gives [leading, ..., constant]
    coeffs = np.poly(roots)
    # Normalize: coeffs[0] should be 1 (monic)
    assert abs(coeffs[0] - 1.0) < 1e-10
    return coeffs  # [1, a_1, ..., a_n]

def discrete_hilbert(roots, i):
    """H_p(lambda_i) = sum_{j != i} 1/(lambda_i - lambda_j)"""
    n = len(roots)
    return sum(1.0 / (roots[i] - roots[j]) for j in range(n) if j != i)

def phi_n(roots):
    """Phi_n(p) = sum_{i=1}^n H_p(lambda_i)^2"""
    n = len(roots)
    return sum(discrete_hilbert(roots, i)**2 for i in range(n))

def finite_free_convolution(a_coeffs, b_coeffs, n):
    """
    Compute r = p ⊞_n q.
    a_coeffs, b_coeffs: coefficient vectors [a_0=1, a_1, ..., a_n]
    c_k = sum_{i+j=k} [(n-i)!(n-j)! / (n!(n-k)!)] a_i b_j
    Returns coefficient vector [c_0=1, c_1, ..., c_n]
    """
    c = np.zeros(n + 1)
    for k in range(n + 1):
        s = 0.0
        for i in range(k + 1):
            j = k - i
            if i <= n and j <= n:
                coeff = (factorial(n - i) * factorial(n - j)) / (factorial(n) * factorial(n - k))
                s += coeff * a_coeffs[i] * b_coeffs[j]
        c[k] = s
    return c

def roots_of_poly(coeffs):
    """Find roots of polynomial given by [c_0, c_1, ..., c_n]"""
    return np.sort(np.roots(coeffs)).real

print("=" * 70)
print("STRICTNESS ANALYSIS FOR n >= 3")
print("=" * 70)

# ===== n = 2 (baseline: expect equality or near-equality) =====
print("\n--- n = 2 ---")
p2_roots = np.array([0.0, 2.0])
q2_roots = np.array([1.0, 4.0])
p2_coeffs = poly_coeffs_from_roots(p2_roots)
q2_coeffs = poly_coeffs_from_roots(q2_roots)
r2_coeffs = finite_free_convolution(p2_coeffs, q2_coeffs, 2)
r2_roots = roots_of_poly(r2_coeffs)

phi_p2 = phi_n(p2_roots)
phi_q2 = phi_n(q2_roots)
phi_r2 = phi_n(r2_roots)

print(f"p roots: {p2_roots}, Phi_2(p) = {phi_p2:.6f}")
print(f"q roots: {q2_roots}, Phi_2(q) = {phi_q2:.6f}")
print(f"r = p ⊞_2 q coeffs: {r2_coeffs}")
print(f"r roots: {r2_roots}, Phi_2(r) = {phi_r2:.6f}")
print(f"1/Phi_2(r) = {1/phi_r2:.6f}")
print(f"1/Phi_2(p) + 1/Phi_2(q) = {1/phi_p2 + 1/phi_q2:.6f}")
print(f"Difference: {1/phi_r2 - (1/phi_p2 + 1/phi_q2):.10f}")

# ===== n = 3, Example 1 =====
print("\n--- n = 3, Example 1 ---")
p3_roots = np.array([0.0, 1.0, 3.0])
q3_roots = np.array([0.0, 2.0, 5.0])
p3_coeffs = poly_coeffs_from_roots(p3_roots)
q3_coeffs = poly_coeffs_from_roots(q3_roots)
r3_coeffs = finite_free_convolution(p3_coeffs, q3_coeffs, 3)
r3_roots = roots_of_poly(r3_coeffs)

phi_p3 = phi_n(p3_roots)
phi_q3 = phi_n(q3_roots)
phi_r3 = phi_n(r3_roots)

print(f"p roots: {p3_roots}, Phi_3(p) = {phi_p3:.6f}")
print(f"q roots: {q3_roots}, Phi_3(q) = {phi_q3:.6f}")
print(f"r = p ⊞_3 q coeffs: {r3_coeffs}")
print(f"r roots: {r3_roots}")
print(f"Phi_3(r) = {phi_r3:.6f}")

# Check simple rootedness
if len(np.unique(np.round(r3_roots, 8))) == 3:
    print(f"r has 3 distinct real roots: GOOD")
else:
    print(f"WARNING: r may not have simple roots!")

print(f"\n1/Phi_3(r) = {1/phi_r3:.6f}")
print(f"1/Phi_3(p) + 1/Phi_3(q) = {1/phi_p3 + 1/phi_q3:.6f}")
diff3 = 1/phi_r3 - (1/phi_p3 + 1/phi_q3)
print(f"Difference (should be > 0 for superadditivity): {diff3:.10f}")
print(f"Strict? {'YES' if diff3 > 1e-12 else 'NO'}")

# ===== n = 3, Example 2 =====
print("\n--- n = 3, Example 2 ---")
p3b_roots = np.array([-1.0, 0.0, 2.0])
q3b_roots = np.array([-2.0, 1.0, 3.0])
p3b_coeffs = poly_coeffs_from_roots(p3b_roots)
q3b_coeffs = poly_coeffs_from_roots(q3b_roots)
r3b_coeffs = finite_free_convolution(p3b_coeffs, q3b_coeffs, 3)
r3b_roots = roots_of_poly(r3b_coeffs)

phi_p3b = phi_n(p3b_roots)
phi_q3b = phi_n(q3b_roots)
phi_r3b = phi_n(r3b_roots)

print(f"p roots: {p3b_roots}, Phi_3(p) = {phi_p3b:.6f}")
print(f"q roots: {q3b_roots}, Phi_3(q) = {phi_q3b:.6f}")
print(f"r = p ⊞_3 q coeffs: {r3b_coeffs}")
print(f"r roots: {r3b_roots}")
print(f"Phi_3(r) = {phi_r3b:.6f}")

print(f"\n1/Phi_3(r) = {1/phi_r3b:.6f}")
print(f"1/Phi_3(p) + 1/Phi_3(q) = {1/phi_p3b + 1/phi_q3b:.6f}")
diff3b = 1/phi_r3b - (1/phi_p3b + 1/phi_q3b)
print(f"Difference: {diff3b:.10f}")
print(f"Strict? {'YES' if diff3b > 1e-12 else 'NO'}")

# ===== n = 4 =====
print("\n--- n = 4 ---")
p4_roots = np.array([0.0, 1.0, 3.0, 6.0])
q4_roots = np.array([0.0, 2.0, 5.0, 9.0])
p4_coeffs = poly_coeffs_from_roots(p4_roots)
q4_coeffs = poly_coeffs_from_roots(q4_roots)
r4_coeffs = finite_free_convolution(p4_coeffs, q4_coeffs, 4)
r4_roots = roots_of_poly(r4_coeffs)

phi_p4 = phi_n(p4_roots)
phi_q4 = phi_n(q4_roots)
phi_r4 = phi_n(r4_roots)

print(f"p roots: {p4_roots}, Phi_4(p) = {phi_p4:.6f}")
print(f"q roots: {q4_roots}, Phi_4(q) = {phi_q4:.6f}")
print(f"r = p ⊞_4 q coeffs: {r4_coeffs}")
print(f"r roots: {r4_roots}")
print(f"Phi_4(r) = {phi_r4:.6f}")

print(f"\n1/Phi_4(r) = {1/phi_r4:.6f}")
print(f"1/Phi_4(p) + 1/Phi_4(q) = {1/phi_p4 + 1/phi_q4:.6f}")
diff4 = 1/phi_r4 - (1/phi_p4 + 1/phi_q4)
print(f"Difference: {diff4:.10f}")
print(f"Strict? {'YES' if diff4 > 1e-12 else 'NO'}")

# ===== n = 5 =====
print("\n--- n = 5 ---")
p5_roots = np.array([0.0, 1.0, 3.0, 6.0, 10.0])
q5_roots = np.array([0.0, 2.0, 5.0, 9.0, 14.0])
p5_coeffs = poly_coeffs_from_roots(p5_roots)
q5_coeffs = poly_coeffs_from_roots(q5_roots)
r5_coeffs = finite_free_convolution(p5_coeffs, q5_coeffs, 5)
r5_roots = roots_of_poly(r5_coeffs)

phi_p5 = phi_n(p5_roots)
phi_q5 = phi_n(q5_roots)
phi_r5 = phi_n(r5_roots)

print(f"p roots: {p5_roots}, Phi_5(p) = {phi_p5:.6f}")
print(f"q roots: {q5_roots}, Phi_5(q) = {phi_q5:.6f}")
print(f"r roots: {r5_roots}")
print(f"Phi_5(r) = {phi_r5:.6f}")

print(f"\n1/Phi_5(r) = {1/phi_r5:.6f}")
print(f"1/Phi_5(p) + 1/Phi_5(q) = {1/phi_p5 + 1/phi_q5:.6f}")
diff5 = 1/phi_r5 - (1/phi_p5 + 1/phi_q5)
print(f"Difference: {diff5:.10f}")
print(f"Strict? {'YES' if diff5 > 1e-12 else 'NO'}")

# ===== Equality Analysis =====
print("\n" + "=" * 70)
print("EQUALITY CONDITION ANALYSIS")
print("=" * 70)

# For n = 2: check if equality holds
print("\n--- n = 2: Detailed analysis ---")
print("For n=2, p(x) = (x-a)(x-b), H_p(a) = 1/(a-b), H_p(b) = 1/(b-a)")
print("Phi_2(p) = 1/(a-b)^2 + 1/(b-a)^2 = 2/(a-b)^2")
print("So 1/Phi_2(p) = (a-b)^2/2 = gap^2/2")
print()

# For n=2 with our example
gap_p = p2_roots[1] - p2_roots[0]
gap_q = q2_roots[1] - q2_roots[0]
gap_r = r2_roots[1] - r2_roots[0]
print(f"gap(p) = {gap_p}, gap(q) = {gap_q}, gap(r) = {gap_r}")
print(f"1/Phi_2(p) = {gap_p**2/2:.6f}, 1/Phi_2(q) = {gap_q**2/2:.6f}")
print(f"1/Phi_2(r) = {gap_r**2/2:.6f}")
print(f"gap(r)^2/2 vs gap(p)^2/2 + gap(q)^2/2: {gap_r**2/2:.6f} vs {gap_p**2/2 + gap_q**2/2:.6f}")

# For n=2, the MSS convolution formula gives:
# c_0 = 1, c_1 = a_1 + b_1, c_2 = a_2 + a_1*b_1*(1/2) + b_2*(1/1) ... let me recompute
print(f"\nn=2 convolution coefficients:")
print(f"p coeffs: {p2_coeffs}")
print(f"q coeffs: {q2_coeffs}")
print(f"r coeffs: {r2_coeffs}")

# Proportionality analysis for n=3
print("\n--- n = 3: Proportionality check ---")
Hp = [discrete_hilbert(p3_roots, i) for i in range(3)]
Hq = [discrete_hilbert(q3_roots, i) for i in range(3)]
print(f"H_p vector: {Hp}")
print(f"H_q vector: {Hq}")
ratios = [Hp[i]/Hq[i] if abs(Hq[i]) > 1e-10 else float('inf') for i in range(3)]
print(f"H_p/H_q ratios: {ratios}")
print(f"Proportional? {'YES' if abs(max(ratios) - min(ratios)) < 1e-10 else 'NO'}")

# n=2 proportionality
Hp2 = [discrete_hilbert(p2_roots, i) for i in range(2)]
Hq2 = [discrete_hilbert(q2_roots, i) for i in range(2)]
print(f"\nn=2: H_p = {Hp2}, H_q = {Hq2}")
ratios2 = [Hp2[i]/Hq2[i] for i in range(2)]
print(f"H_p/H_q ratios: {ratios2}")
print(f"Proportional? {'YES' if abs(ratios2[0] - ratios2[1]) < 1e-10 else 'NO'}")

# ===== n=2 analytical verification =====
print("\n" + "=" * 70)
print("ANALYTICAL n=2 CHECK")
print("=" * 70)
print("""
For n=2: p(x) = (x-a)(x-b), q(x) = (x-c)(x-d)
  a_1 = -(a+b), a_2 = ab
  b_1 = -(c+d), b_2 = cd

MSS convolution coefficients for n=2:
  c_0 = 1
  c_1 = a_1 + b_1 = -(a+b+c+d)
  c_2 = (2!*2!)/(2!*0!) * a_0*b_2 * ... Let me use the formula directly.

c_k = sum_{i+j=k} [(n-i)!(n-j)! / (n!(n-k)!)] a_i b_j  with n=2

c_0: i=0,j=0: [2!*2!/(2!*2!)] * 1*1 = 1
c_1: i=0,j=1: [2!*1!/(2!*1!)] * 1*b_1 = b_1
    + i=1,j=0: [1!*2!/(2!*1!)] * a_1*1 = a_1
  So c_1 = a_1 + b_1

c_2: i=0,j=2: [2!*0!/(2!*0!)] * 1*b_2 = b_2
    + i=1,j=1: [1!*1!/(2!*0!)] * a_1*b_1 = a_1*b_1/2
    + i=2,j=0: [0!*2!/(2!*0!)] * a_2*1 = a_2
  So c_2 = a_2 + a_1*b_1/2 + b_2
""")

a, b = 0.0, 2.0
c, d = 1.0, 4.0
a1 = -(a+b)
a2 = a*b
b1 = -(c+d)
b2 = c*d
c1_manual = a1 + b1
c2_manual = a2 + a1*b1/2 + b2
print(f"Manual: c_1 = {c1_manual}, c_2 = {c2_manual}")
print(f"Computed: c_1 = {r2_coeffs[1]}, c_2 = {r2_coeffs[2]}")

# r(x) = x^2 + c_1*x + c_2, roots: (-c_1 +/- sqrt(c_1^2 - 4*c_2))/2
disc = c1_manual**2 - 4*c2_manual
print(f"Discriminant: {disc}")
r_root1 = (-c1_manual - np.sqrt(disc))/2
r_root2 = (-c1_manual + np.sqrt(disc))/2
gap_r_manual = r_root2 - r_root1
print(f"r roots: [{r_root1}, {r_root2}], gap = {gap_r_manual}")

# 1/Phi_2(r) = gap^2/2
# 1/Phi_2(p) + 1/Phi_2(q) = (b-a)^2/2 + (d-c)^2/2
lhs = gap_r_manual**2 / 2
rhs = (b-a)**2/2 + (d-c)**2/2
print(f"\n1/Phi_2(r) = {lhs:.6f}")
print(f"1/Phi_2(p) + 1/Phi_2(q) = {rhs:.6f}")
print(f"Difference: {lhs - rhs:.10f}")

# Analytically for n=2:
# gap(r) = sqrt(c_1^2 - 4c_2)
# c_1 = -(a+b+c+d), c_2 = ab + (a+b)(c+d)/2 + cd
# c_1^2 - 4c_2 = (a+b+c+d)^2 - 4[ab + (a+b)(c+d)/2 + cd]
# = (a+b)^2 + 2(a+b)(c+d) + (c+d)^2 - 4ab - 2(a+b)(c+d) - 4cd
# = (a+b)^2 - 4ab + (c+d)^2 - 4cd
# = (a-b)^2 + (c-d)^2
# So gap(r) = sqrt((b-a)^2 + (d-c)^2)
# And 1/Phi_2(r) = [(b-a)^2 + (d-c)^2]/2 = 1/Phi_2(p) + 1/Phi_2(q). EQUALITY!
print(f"\nANALYTICAL PROOF for n=2:")
print(f"gap(r)^2 = (b-a)^2 + (d-c)^2 = {(b-a)**2} + {(d-c)**2} = {(b-a)**2 + (d-c)**2}")
print(f"So 1/Phi_2(r) = gap(r)^2/2 = [(b-a)^2 + (d-c)^2]/2 = 1/Phi_2(p) + 1/Phi_2(q)")
print("EQUALITY ALWAYS HOLDS for n=2!")

# ===== More n=3 examples to show strictness is generic =====
print("\n" + "=" * 70)
print("SYSTEMATIC n=3 STRICTNESS CHECK")
print("=" * 70)

np.random.seed(42)
strict_count = 0
total_count = 0
min_diff = float('inf')

for trial in range(20):
    # Generate random roots (sorted, distinct)
    p_roots = np.sort(np.random.uniform(-10, 10, 3))
    while min(np.diff(p_roots)) < 0.5:
        p_roots = np.sort(np.random.uniform(-10, 10, 3))

    q_roots = np.sort(np.random.uniform(-10, 10, 3))
    while min(np.diff(q_roots)) < 0.5:
        q_roots = np.sort(np.random.uniform(-10, 10, 3))

    p_c = poly_coeffs_from_roots(p_roots)
    q_c = poly_coeffs_from_roots(q_roots)
    r_c = finite_free_convolution(p_c, q_c, 3)
    r_roots_trial = roots_of_poly(r_c)

    # Check r has 3 distinct real roots
    if np.max(np.abs(r_roots_trial.imag)) > 1e-6:
        continue
    r_roots_trial = np.sort(r_roots_trial.real)
    if min(np.diff(r_roots_trial)) < 1e-8:
        continue

    total_count += 1
    pp = phi_n(p_roots)
    pq = phi_n(q_roots)
    pr = phi_n(r_roots_trial)

    diff_val = 1/pr - (1/pp + 1/pq)
    if diff_val > 1e-12:
        strict_count += 1
    min_diff = min(min_diff, diff_val)

print(f"Tested {total_count} random pairs for n=3")
print(f"Strict inequality: {strict_count}/{total_count}")
print(f"Minimum difference observed: {min_diff:.10e}")
if strict_count == total_count:
    print("ALL examples show STRICT inequality for n=3!")
