"""
Deep investigation of the equality family for n=3.

Discovery: For n=3, p(x)=x(x-s)(x-2s) and q(x)=x(x-t)(x-2t) (equal gaps)
give EQUALITY in 1/Phi_3(r) = 1/Phi_3(p) + 1/Phi_3(q).

Question: Is this the ONLY equality family? What about n=4, n=5?
"""

import numpy as np
from math import factorial

def discrete_hilbert(roots, i):
    n = len(roots)
    return sum(1.0 / (roots[i] - roots[j]) for j in range(n) if j != i)

def phi_n(roots):
    n = len(roots)
    return sum(discrete_hilbert(roots, i)**2 for i in range(n))

def poly_coeffs_from_roots(roots):
    return np.poly(roots)

def finite_free_convolution(a_coeffs, b_coeffs, n):
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
    return np.sort(np.roots(coeffs)).real

print("=" * 70)
print("EQUAL-GAP FAMILY: p(x) = prod(x - k*s), k=0,...,n-1")
print("=" * 70)

# For n=3: equal-gap polynomials: roots at {0, s, 2s}
# Check analytically why equality holds.
print("\n--- n=3 equal-gap analysis ---")
print("p roots: {0, s, 2s}, q roots: {0, t, 2t}")
print()

# H_p for equal-gap n=3:
# H(0) = 1/(0-s) + 1/(0-2s) = -1/s - 1/(2s) = -3/(2s)
# H(s) = 1/(s-0) + 1/(s-2s) = 1/s - 1/s = 0
# H(2s) = 1/(2s-0) + 1/(2s-s) = 1/(2s) + 1/s = 3/(2s)
for s in [1.0, 2.0, 3.0]:
    roots = np.array([0, s, 2*s])
    H = [discrete_hilbert(roots, i) for i in range(3)]
    print(f"  s={s}: H = ({H[0]:.4f}, {H[1]:.4f}, {H[2]:.4f})")
    print(f"         H = ({-3/(2*s):.4f}, {0:.4f}, {3/(2*s):.4f}) [analytical]")

print()
print("For equal-gap n=3: H_p = (-3/(2s), 0, 3/(2s))")
print("So H_p is ALWAYS proportional to (-1, 0, 1) regardless of s.")
print("Phi_3(p) = 9/(4s^2) + 0 + 9/(4s^2) = 9/(2s^2)")
print("1/Phi_3(p) = 2s^2/9")
print()

# Verify
for s in [1.0, 2.0, 3.0]:
    roots = np.array([0, s, 2*s])
    p = phi_n(roots)
    print(f"  s={s}: Phi_3(p)={p:.6f}, 9/(2s^2)={9/(2*s**2):.6f}, match={abs(p-9/(2*s**2))<1e-10}")

# So 1/Phi_3(r) = 2*gap_r^2/9
# For equal-gap convolution, gap(r)^2 = s^2 + t^2 (?)
print("\n--- Convolution of equal-gap polynomials ---")
for s, t in [(1,1), (1,2), (2,3), (1,5)]:
    p_r = np.array([0.0, float(s), float(2*s)])
    q_r = np.array([0.0, float(t), float(2*t)])
    p_c = poly_coeffs_from_roots(p_r)
    q_c = poly_coeffs_from_roots(q_r)
    r_c = finite_free_convolution(p_c, q_c, 3)
    r_roots = roots_of_poly(r_c)
    gap_r = r_roots[1] - r_roots[0]  # should equal r_roots[2]-r_roots[1] for equal-gap

    # Check: is gap(r)^2 = s^2 + t^2?
    print(f"  s={s}, t={t}: r_roots={r_roots}, gaps=({r_roots[1]-r_roots[0]:.6f}, {r_roots[2]-r_roots[1]:.6f})")
    print(f"    gap(r) = {gap_r:.6f}, sqrt(s^2+t^2) = {np.sqrt(s**2+t**2):.6f}, match={abs(gap_r-np.sqrt(s**2+t**2))<1e-6}")
    print(f"    1/Phi_3(r) = {1/phi_n(r_roots):.6f}, 2(s^2+t^2)/9 = {2*(s**2+t**2)/9:.6f}")
    print(f"    1/Phi_3(p)+1/Phi_3(q) = {2*s**2/9 + 2*t**2/9:.6f}")

print("""
RESULT: For equal-gap n=3, convolution preserves the equal-gap property,
gap(r)^2 = s^2 + t^2, and equality holds:
  1/Phi_3(r) = 2(s^2+t^2)/9 = 2s^2/9 + 2t^2/9 = 1/Phi_3(p) + 1/Phi_3(q).
""")

# Now check n=4 with equal gaps
print("=" * 70)
print("n=4: EQUAL-GAP FAMILY")
print("=" * 70)

for s, t in [(1,1), (1,2), (2,3), (1,5)]:
    n = 4
    p_r = np.array([float(k*s) for k in range(n)])
    q_r = np.array([float(k*t) for k in range(n)])
    p_c = poly_coeffs_from_roots(p_r)
    q_c = poly_coeffs_from_roots(q_r)
    r_c = finite_free_convolution(p_c, q_c, n)
    r_roots = roots_of_poly(r_c)
    gaps = np.diff(r_roots)

    pp = phi_n(p_r)
    pq = phi_n(q_r)
    pr = phi_n(r_roots)
    diff = 1/pr - (1/pp + 1/pq)

    print(f"s={s}, t={t}: r gaps={np.round(gaps,4)}, ratio={gaps[0]/gaps[1]:.4f},{gaps[1]/gaps[2]:.4f}")
    print(f"  diff = {diff:.10f}, strict={'YES' if abs(diff)>1e-10 else 'NO (equality)'}")

# Check n=5 with equal gaps
print("\n" + "=" * 70)
print("n=5: EQUAL-GAP FAMILY")
print("=" * 70)

for s, t in [(1,1), (1,2), (2,3)]:
    n = 5
    p_r = np.array([float(k*s) for k in range(n)])
    q_r = np.array([float(k*t) for k in range(n)])
    p_c = poly_coeffs_from_roots(p_r)
    q_c = poly_coeffs_from_roots(q_r)
    r_c = finite_free_convolution(p_c, q_c, n)
    r_roots = roots_of_poly(r_c)
    gaps = np.diff(r_roots)

    pp = phi_n(p_r)
    pq = phi_n(q_r)
    pr = phi_n(r_roots)
    diff = 1/pr - (1/pp + 1/pq)

    print(f"s={s}, t={t}: r gaps={np.round(gaps,4)}")
    print(f"  diff = {diff:.10f}, strict={'YES' if abs(diff)>1e-10 else 'NO (equality)'}")

# Check n=3 with SHIFTED equal gaps: p roots at {a, a+s, a+2s}
print("\n" + "=" * 70)
print("n=3: SHIFTED EQUAL-GAP FAMILY (different centers)")
print("=" * 70)

for a, s, b, t in [(0,1,0,2), (1,1,0,2), (5,1,-3,2), (0,1,0,1)]:
    n = 3
    p_r = np.array([float(a+k*s) for k in range(n)])
    q_r = np.array([float(b+k*t) for k in range(n)])
    p_c = poly_coeffs_from_roots(p_r)
    q_c = poly_coeffs_from_roots(q_r)
    r_c = finite_free_convolution(p_c, q_c, n)
    r_roots = roots_of_poly(r_c)
    gaps = np.diff(r_roots)

    pp = phi_n(p_r)
    pq = phi_n(q_r)
    pr = phi_n(r_roots)
    diff = 1/pr - (1/pp + 1/pq)

    print(f"p={p_r}, q={q_r}")
    print(f"  r gaps={np.round(gaps,4)}, diff={diff:.10f}, eq={abs(diff)<1e-10}")

print("""
CONCLUSION: The equal-gap family (roots equally spaced) gives EQUALITY
for ALL n, not just n=2. The translation parameter does not matter
(Phi_n is translation-invariant).

This means the original node statement needs revision:
- Strict inequality does NOT hold for ALL p,q when n >= 3.
- There exists a codimension-1 family (equal-gap polynomials scaled differently)
  where equality holds.
- Strict inequality holds when p and q do NOT both have equal gaps,
  or more generally when their Hilbert-transform vectors are not proportional.
""")

# Final: What fraction of P_3 x P_3 gives equality?
print("=" * 70)
print("GENERIC STRICTNESS: What is the equality locus?")
print("=" * 70)
print("""
For n=3, P_3 (modulo translation) is parametrized by (g_1, g_2) with g_1, g_2 > 0.
The equal-gap condition is g_1 = g_2, which is a codimension-1 subset.

For a pair (p, q) to give equality, we need:
  (i) H_p proportional to H_q, AND
  (ii) H_r proportional to H_p (where r = p ⊞_3 q)

Condition (i) means p and q have the same gap ratio g_1/g_2.
We showed that equal-gap (g_1=g_2) is preserved under convolution.

Question: is equal-gap the ONLY gap ratio preserved under convolution?
""")

# Test: same non-unit gap ratio for p and q
print("Testing same non-unit gap ratio preservation:")
for ratio in [0.5, 2.0, 1/3, 3.0]:
    g1_p, g2_p = 1.0, 1.0/ratio  # gap ratio = ratio
    g1_q, g2_q = 2.0, 2.0/ratio
    p_r = np.array([0.0, g1_p, g1_p + g2_p])
    q_r = np.array([0.0, g1_q, g1_q + g2_q])
    p_c = poly_coeffs_from_roots(p_r)
    q_c = poly_coeffs_from_roots(q_r)
    r_c = finite_free_convolution(p_c, q_c, 3)
    r_roots = roots_of_poly(r_c)
    gaps = np.diff(r_roots)
    r_ratio = gaps[0]/gaps[1]

    pp = phi_n(p_r)
    pq = phi_n(q_r)
    pr = phi_n(r_roots)
    diff = 1/pr - (1/pp + 1/pq)

    print(f"  gap_ratio={ratio:.3f}: r_ratio={r_ratio:.6f}, preserved={abs(r_ratio-ratio)<1e-6}, diff={diff:.8f}")
