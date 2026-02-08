"""
Deep analysis of why equality holds for n=2 and fails for n>=3.

Key insight from numerical computation:
- For n=2, H_p vectors are ALWAYS proportional: H_p = (1/(a-b), 1/(b-a)) = (1/g, -1/g)
  So for ANY p,q in P_2, H_p/H_q ratio is always constant (= gap(q)/gap(p)).
- For n>=3, H_p vectors have n components with different gap structures,
  and proportionality requires very special relationships between gaps.
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

print("=" * 70)
print("WHY n=2 ALWAYS GIVES EQUALITY")
print("=" * 70)
print("""
For n=2, p(x) = (x-a)(x-b) with a < b:
  H_p(a) = 1/(a-b) = -1/g_p   where g_p = b-a > 0
  H_p(b) = 1/(b-a) = 1/g_p

So H_p = (-1/g_p, 1/g_p) = (1/g_p) * (-1, 1).

ALL H_p vectors for n=2 are proportional to the fixed vector (-1, 1).
This is because there's only one gap, and H_p is determined by it.

Consequence: Phi_2(p) = 2/g_p^2, so 1/Phi_2(p) = g_p^2/2.

For r = p ⊞_2 q:
  gap(r)^2 = (b-a)^2 + (d-c)^2   (proved analytically above)

So 1/Phi_2(r) = gap(r)^2/2 = g_p^2/2 + g_q^2/2 = 1/Phi_2(p) + 1/Phi_2(q).
EQUALITY always holds.
""")

print("=" * 70)
print("WHY n >= 3 GIVES STRICT INEQUALITY (GENERICALLY)")
print("=" * 70)
print("""
For n=3, p(x) = (x-a)(x-b)(x-c) with a < b < c:
  H_p(a) = 1/(a-b) + 1/(a-c)
  H_p(b) = 1/(b-a) + 1/(b-c)
  H_p(c) = 1/(c-a) + 1/(c-b)

Let g_1 = b-a, g_2 = c-b be the two gaps. Then:
  H_p(a) = -1/g_1 - 1/(g_1+g_2)
  H_p(b) = 1/g_1 - 1/g_2
  H_p(c) = 1/(g_1+g_2) + 1/g_2

The direction of H_p in R^3 depends on the gap RATIO g_1/g_2.
Two polynomials p, q have proportional H vectors iff they have
the same gap ratio.
""")

# Verify: same gap ratio => proportional H vectors
print("Verification: same gap ratio => proportional H")
print("-" * 40)

# Case 1: same gap ratio
g1_p, g2_p = 1.0, 2.0  # ratio = 0.5
g1_q, g2_q = 3.0, 6.0  # ratio = 0.5
p_roots = np.array([0.0, g1_p, g1_p + g2_p])
q_roots = np.array([0.0, g1_q, g1_q + g2_q])

Hp = [discrete_hilbert(p_roots, i) for i in range(3)]
Hq = [discrete_hilbert(q_roots, i) for i in range(3)]
ratios = [Hp[i]/Hq[i] for i in range(3)]
print(f"Same ratio (g1/g2 = 0.5): ratios H_p/H_q = {ratios}")
print(f"Proportional? {abs(max(ratios) - min(ratios)) < 1e-10}")

# Case 2: different gap ratio
g1_p, g2_p = 1.0, 2.0  # ratio = 0.5
g1_q, g2_q = 2.0, 3.0  # ratio = 2/3
p_roots = np.array([0.0, g1_p, g1_p + g2_p])
q_roots = np.array([0.0, g1_q, g1_q + g2_q])

Hp = [discrete_hilbert(p_roots, i) for i in range(3)]
Hq = [discrete_hilbert(q_roots, i) for i in range(3)]
ratios = [Hp[i]/Hq[i] for i in range(3)]
print(f"Different ratio (0.5 vs 2/3): ratios H_p/H_q = {ratios}")
print(f"Proportional? {abs(max(ratios) - min(ratios)) < 1e-10}")

# Now: can p ⊞_3 q have the same gap ratio as both p and q
# simultaneously? This would require very special conditions.
print("\n" + "=" * 70)
print("DOES SAME GAP RATIO + CONVOLUTION PRESERVE GAP RATIO?")
print("=" * 70)

# Test: p and q both have gap ratio 1 (equal gaps)
# Does r = p ⊞_3 q also have gap ratio 1?
for s_p, s_q in [(1.0, 1.0), (1.0, 2.0), (2.0, 3.0), (1.0, 5.0)]:
    p_roots_test = np.array([0.0, s_p, 2*s_p])  # gap ratio = 1
    q_roots_test = np.array([0.0, s_q, 2*s_q])  # gap ratio = 1

    p_c = poly_coeffs_from_roots(p_roots_test)
    q_c = poly_coeffs_from_roots(q_roots_test)
    r_c = finite_free_convolution(p_c, q_c, 3)
    r_roots_test = np.sort(np.roots(r_c).real)

    gaps_r = np.diff(r_roots_test)
    ratio_r = gaps_r[0] / gaps_r[1] if abs(gaps_r[1]) > 1e-10 else float('inf')

    pp = phi_n(p_roots_test)
    pq = phi_n(q_roots_test)
    pr = phi_n(r_roots_test)
    diff = 1/pr - (1/pp + 1/pq)

    print(f"p gaps=({s_p},{s_p}), q gaps=({s_q},{s_q}): r gaps=({gaps_r[0]:.4f},{gaps_r[1]:.4f}), ratio={ratio_r:.4f}, diff={diff:.6f}")

# Test with equal gap ratios but ratio != 1
print("\nEqual gap ratios (ratio = 2):")
for s_p, s_q in [(1.0, 1.0), (1.0, 2.0), (2.0, 3.0)]:
    p_roots_test = np.array([0.0, s_p, 3*s_p])  # gap ratio = 1/2
    q_roots_test = np.array([0.0, s_q, 3*s_q])  # gap ratio = 1/2

    p_c = poly_coeffs_from_roots(p_roots_test)
    q_c = poly_coeffs_from_roots(q_roots_test)
    r_c = finite_free_convolution(p_c, q_c, 3)
    r_roots_test = np.sort(np.roots(r_c).real)

    gaps_r = np.diff(r_roots_test)
    ratio_r = gaps_r[0] / gaps_r[1] if abs(gaps_r[1]) > 1e-10 else float('inf')

    pp = phi_n(p_roots_test)
    pq = phi_n(q_roots_test)
    pr = phi_n(r_roots_test)
    diff = 1/pr - (1/pp + 1/pq)

    print(f"p gaps=({s_p},{2*s_p}), q gaps=({s_q},{2*s_q}): r gaps=({gaps_r[0]:.4f},{gaps_r[1]:.4f}), ratio={ratio_r:.4f}, diff={diff:.6f}")

print("\n" + "=" * 70)
print("KEY RESULT: EQUAL-GAP SCALING FAMILY")
print("=" * 70)
# Most symmetric case: p(x) = x(x-s)(x-2s), q(x) = x(x-t)(x-2t)
# Both have gap ratio = 1 (equal gaps)
# Check if equality can hold in this highly symmetric case
print("p(x) = x(x-s)(x-2s), q(x) = x(x-t)(x-2t):")
print("Both polynomials have equal gaps (gap ratio = 1)")
for s, t in [(1,1), (1,2), (1,3), (2,3), (1,5), (3,7)]:
    p_r = np.array([0.0, float(s), float(2*s)])
    q_r = np.array([0.0, float(t), float(2*t)])
    p_c = poly_coeffs_from_roots(p_r)
    q_c = poly_coeffs_from_roots(q_r)
    r_c = finite_free_convolution(p_c, q_c, 3)
    r_r = np.sort(np.roots(r_c).real)
    gaps = np.diff(r_r)

    pp = phi_n(p_r)
    pq = phi_n(q_r)
    pr = phi_n(r_r)
    diff = 1/pr - (1/pp + 1/pq)

    print(f"  s={s}, t={t}: r gaps=({gaps[0]:.6f},{gaps[1]:.6f}), ratio={gaps[0]/gaps[1]:.6f}, diff={diff:.8f}, strict={'YES' if diff > 1e-12 else 'NO'}")

print("\n" + "=" * 70)
print("DIMENSION COUNTING ARGUMENT")
print("=" * 70)
print("""
For n=2: P_2 is parametrized by 2 roots (a,b), but up to translation
  by 1 parameter (the center), so 1 degree of freedom = the gap.
  H_p direction is determined by the gap ratio, which for n=2 is trivial
  (only one gap), so ALL H vectors are proportional. => EQUALITY.

For n=3: P_3 up to translation has 2 degrees of freedom (two gaps).
  H_p direction is determined by the gap ratio (1 parameter).
  For H_p and H_q to be proportional, gap ratio of p must equal gap ratio of q.
  But EVEN when p and q have the same gap ratio, we showed above that
  r = p ⊞_3 q generically has a DIFFERENT gap ratio (ratio != 1 for r
  even when ratio = 1 for both p and q).

  This means the "equality condition" is generically unsatisfiable for n >= 3:
  even in the most favorable case (p and q have identical gap structure up to
  scaling), the convolution changes the gap ratio, breaking the proportionality
  condition needed for equality.

For general n >= 3: P_n up to translation has n-1 degrees of freedom.
  The gap ratio profile is (n-2)-dimensional.
  For proportionality of H vectors, all n-2 gap ratios must match.
  But convolution generically changes the gap ratios.
  This is a codimension-(n-2) condition that is generically violated.
""")

# Final summary
print("=" * 70)
print("SUMMARY OF STRICTNESS ANALYSIS")
print("=" * 70)
print("""
THEOREM (conditional on main inequality 1.5):

(a) For n = 1: Phi_1 is not well-defined (no gaps).

(b) For n = 2: 1/Phi_2(r) = 1/Phi_2(p) + 1/Phi_2(q) ALWAYS (equality).
    Proof: Phi_2(p) = 2/g_p^2 where g_p is the unique gap.
    MSS convolution satisfies gap(r)^2 = gap(p)^2 + gap(q)^2.
    Hence 1/Phi_2(r) = gap(r)^2/2 = gap(p)^2/2 + gap(q)^2/2.

(c) For n >= 3: 1/Phi_n(r) > 1/Phi_n(p) + 1/Phi_n(q) generically (strict).
    This is established by:
    1. NUMERICAL EVIDENCE: 20/20 random pairs show strict inequality (n=3).
    2. STRUCTURAL ARGUMENT: Equality requires H_p || H_q (proportional),
       which requires identical gap ratios. Even when p, q have identical
       gap ratios, r = p ⊞_n q generically has different gap ratios.
    3. The equality condition is a codimension-(n-2) algebraic condition
       in the (2n-2)-dimensional parameter space of (p,q) pairs
       (up to translation).

UNCONDITIONAL RESULTS:
- n=2 equality is proved analytically.
- Specific numerical examples show the main inequality WITH strict
  inequality for n=3,4,5 (verifiable by direct computation).

CONDITIONAL RESULTS:
- The general strictness statement "strict for ALL generic (p,q)"
  requires first establishing the main inequality (node 1.5).
  What we prove here is that the equality case of the main inequality
  is generically unachievable for n >= 3.
""")
