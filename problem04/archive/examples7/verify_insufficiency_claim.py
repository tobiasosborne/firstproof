"""
VERIFIER-6: Deep check on the "insufficiency" claim.

The node claims that because Phi_n is rational (not polynomial) in cumulants,
the cumulant approach is "insufficient" for proving superadditivity.

This is a LOGICAL claim, not a computational one. Let me stress-test it:

1. For n=3, Phi_3 = (9/2)*k2^2 / (k2^3 - k3^2/3) where k2, k3 are additive.
   Question: Is 1/Phi_3 superadditive? i.e., is
   f(k2_p+k2_q, k3_p+k3_q) >= f(k2_p, k3_p) + f(k2_q, k3_q)
   where f(k2, k3) = 1/Phi_3 = (2/9) * (k2^3 - k3^2/3) / k2^2
                    = (2/9) * (k2 - k3^2/(3*k2^2))

2. f(k2, k3) = (2/9)*(k2 - k3^2/(3*k2^2))

   So superadditivity of 1/Phi_3 is:
   (k2_p+k2_q) - (k3_p+k3_q)^2/(3*(k2_p+k2_q)^2)
   >= k2_p - k3_p^2/(3*k2_p^2) + k2_q - k3_q^2/(3*k2_q^2)

   The k2 terms cancel (additive): k2_p+k2_q = k2_p+k2_q.
   So we need:
   -(k3_p+k3_q)^2 / (k2_p+k2_q)^2 >= -k3_p^2/k2_p^2 - k3_q^2/k2_q^2

   Multiply by -1 (flip inequality):
   (k3_p+k3_q)^2 / (k2_p+k2_q)^2 <= k3_p^2/k2_p^2 + k3_q^2/k2_q^2

   This is saying: (a+b)^2/(c+d)^2 <= a^2/c^2 + b^2/d^2
   where a=k3_p, b=k3_q, c=k2_p, d=k2_q.

   Is this true? Let's check...
"""
import numpy as np

print("="*60)
print("DEEP CHECK: Is the n=3 superadditivity provable from cumulants alone?")
print("="*60)

print("\n1/Phi_3 = (2/9)*(k2 - k3^2/(3*k2^2))")
print("Superadditivity reduces to: (k3_p+k3_q)^2/(k2_p+k2_q)^2 <= k3_p^2/k2_p^2 + k3_q^2/k2_q^2")
print()

# Test this inequality
def test_ineq(a, b, c, d):
    """Test (a+b)^2/(c+d)^2 <= a^2/c^2 + b^2/d^2"""
    lhs = (a+b)**2 / (c+d)**2
    rhs = a**2/c**2 + b**2/d**2
    return lhs, rhs, lhs <= rhs + 1e-12

# This is actually a known inequality related to Cauchy-Schwarz!
# By Cauchy-Schwarz: (a/c + b/d)^2 * (c^2 + d^2) >= (a+b)^2
# Hmm, that's not quite the same...

# Actually, let's think about it differently.
# Let u = k3_p/k2_p, v = k3_q/k2_q, s = k2_p/(k2_p+k2_q), t = k2_q/(k2_p+k2_q)
# Then s+t = 1, and the LHS is (s*u + t*v)^2 and the RHS is s*u^2/s + t*v^2/t = u^2 + v^2... no

# Let me just directly verify:
# (a+b)^2/(c+d)^2 <= a^2/c^2 + b^2/d^2  where c, d > 0

# Actually multiply out. Let x = a/c, y = b/d, alpha = c/(c+d), beta = d/(c+d).
# LHS = (alpha*x + beta*y)^2
# RHS = x^2 + y^2
# So we need (alpha*x + beta*y)^2 <= x^2 + y^2

# Since alpha, beta in (0,1) and alpha+beta=1:
# (alpha*x + beta*y)^2 <= (by Jensen, since x^2 is convex) alpha*x^2 + beta*y^2 <= x^2 + y^2
# The last step because alpha < 1 and beta < 1.

# Wait, this PROVES the inequality! Let me verify numerically.

print("Testing (a+b)^2/(c+d)^2 <= a^2/c^2 + b^2/d^2 for c,d > 0:")
violations = 0
for trial in range(10000):
    a = np.random.randn()
    b = np.random.randn()
    c = abs(np.random.randn()) + 0.01  # c > 0
    d = abs(np.random.randn()) + 0.01  # d > 0

    lhs, rhs, ok = test_ineq(a, b, c, d)
    if not ok:
        violations += 1
        print(f"  VIOLATION: a={a:.4f}, b={b:.4f}, c={c:.4f}, d={d:.4f}: "
              f"LHS={lhs:.6f}, RHS={rhs:.6f}")

print(f"  {violations} violations out of 10000 trials")

# PROOF that (a+b)^2/(c+d)^2 <= a^2/c^2 + b^2/d^2:
# Let x = a/c, y = b/d, s = c/(c+d), t = d/(c+d), so s+t=1.
# LHS = (sx + ty)^2 = s^2*x^2 + 2st*xy + t^2*y^2
# RHS = x^2 + y^2
# RHS - LHS = x^2(1-s^2) + y^2(1-t^2) - 2st*xy
#            = x^2*t*(1+s) + y^2*s*(1+t) - 2st*xy
#            (using 1-s^2 = (1-s)(1+s) = t(1+s))
#            = t(1+s)*x^2 + s(1+t)*y^2 - 2st*xy
# Since s+t=1: 1+s = 2s+t, 1+t = s+2t... Hmm, let me just use AM-GM or Cauchy-Schwarz.
#
# By Cauchy-Schwarz: (sx + ty)^2 <= (s^2/s + t^2/t)(s*x^2 + t*y^2)
# Wait, that gives (s+t)(sx^2+ty^2) = sx^2+ty^2 which is <= x^2+y^2 since s,t < 1.
# Actually Cauchy-Schwarz gives (sx+ty)^2 <= (s+t)(sx^2+ty^2) = sx^2+ty^2 <= x^2+y^2. YES!

print("\nPROOF: By Cauchy-Schwarz (or Jensen for convex x^2):")
print("  (sx+ty)^2 <= (s+t)(sx^2+ty^2) = sx^2+ty^2 <= x^2+y^2")
print("  where x=a/c, y=b/d, s=c/(c+d), t=d/(c+d), s+t=1, 0<s,t<1")
print()

# BUT WAIT: This only works when k2_p, k2_q > 0.
# For real-rooted polynomials with distinct roots, is k2 > 0?
# k2 = -e_2 (for n=3 centered). e_2 = sum_{i<j} lambda_i*lambda_j.
# For centered roots, e_1 = 0 and e_2 = -sum(lambda_i^2)/2 (since (sum lambda_i)^2 = 0 = sum lambda_i^2 + 2*e_2)
# So e_2 = -sum(lambda_i^2)/2 < 0, hence k2 = -e_2 = sum(lambda_i^2)/2 > 0. GOOD.

print("For centered n=3: k2 = sum(lambda_i^2)/2 > 0 always. CHECK.")
print()

# So for n=3, the superadditivity of 1/Phi_3 IS provable from the cumulant structure!
# The node's claim that the cumulant approach is "insufficient" is WRONG for n=3!

print("="*60)
print("CONCLUSION FOR n=3:")
print("="*60)
print("1/Phi_3 = (2/9)*(k2 - k3^2/(3*k2^2))")
print("         = (2/9)*k2 - (2/27)*k3^2/k2^2")
print()
print("Superadditivity of 1/Phi_3 with additive k2, k3 reduces to:")
print("  (k3_p+k3_q)^2/(k2_p+k2_q)^2 <= k3_p^2/k2_p^2 + k3_q^2/k2_q^2")
print()
print("This follows from Cauchy-Schwarz (Jensen's inequality for x^2).")
print("The cumulant approach IS SUFFICIENT for n=3!")
print()
print("The node's claim that 'a DIFFERENT technique is needed' is INCORRECT")
print("for n=3. The fact that Phi_n is rational does NOT prevent proving")
print("superadditivity from cumulant additivity.")

# ============================================================
# Now check: does this extend to n=4?
# ============================================================
print("\n" + "="*60)
print("CHECK: Can similar argument work for n=4?")
print("="*60)

# For n=4, centered, the structure is more complex.
# We need to know Phi_4 in terms of kappa_2, kappa_3, kappa_4.
# From Part B: Phi_4 * disc = -8*e_2^5 - 36*e_2^2*e_3^2 - 64*e_2^3*e_4 - 432*e_3^2*e_4 + 384*e_2*e_4^2

# The cumulant-to-e_k transformation for n=4:
# tilde_a_k = (-1)^k * a_k / C(4,k)
# a_1 = -e_1 = 0, a_2 = e_2, a_3 = -e_3, a_4 = e_4
# tilde_a_1 = 0, tilde_a_2 = e_2/C(4,2) = e_2/6
# tilde_a_3 = e_3/C(4,3) = e_3/4
# tilde_a_4 = e_4/C(4,4) = e_4

# The relationship between kappa and e is more involved for n=4.
# Let's just check numerically whether 1/Phi_4 is a superadditive function
# of the additive cumulant vector.

print("Testing superadditivity of 1/Phi_4 numerically (1000 trials)...")

from itertools import combinations

def elementary_symmetric(roots, k):
    n = len(roots)
    if k == 0: return 1.0
    if k > n: return 0.0
    return sum(np.prod(combo) for combo in combinations(roots, k))

def disc(roots):
    n = len(roots)
    d = 1.0
    for i in range(n):
        for j in range(i+1, n):
            d *= (roots[i] - roots[j])**2
    return d

def phi_n(roots):
    n = len(roots)
    total = 0.0
    for i in range(n):
        pprime = 1.0
        for j in range(n):
            if j != i:
                pprime *= (roots[i] - roots[j])
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

from math import factorial, comb

def mss_convolve(roots_p, roots_q):
    n = len(roots_p)
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

violations = 0
total_valid = 0
for trial in range(1000):
    roots_p = np.sort(np.random.randn(4) * 2)
    roots_q = np.sort(np.random.randn(4) * 2)

    min_gap_p = min(np.diff(roots_p))
    min_gap_q = min(np.diff(roots_q))
    if min_gap_p < 0.2 or min_gap_q < 0.2:
        continue

    try:
        roots_r = mss_convolve(roots_p, roots_q)
        # Check roots are real
        if np.max(np.abs(np.imag(np.roots(np.poly(roots_p))))) > 0.01:
            continue

        phi_p = phi_n(roots_p)
        phi_q = phi_n(roots_q)
        phi_r = phi_n(roots_r)

        if phi_p < 1e-10 or phi_q < 1e-10 or phi_r < 1e-10:
            continue

        total_valid += 1
        lhs = 1/phi_r
        rhs = 1/phi_p + 1/phi_q

        if lhs < rhs - 1e-6:
            violations += 1
            print(f"  VIOLATION: 1/Phi(r)={lhs:.6f} < 1/Phi(p)+1/Phi(q)={rhs:.6f}")
    except:
        pass

print(f"  {violations} violations out of {total_valid} valid trials")

# ============================================================
# Check: Is the n=3 proof above actually correct for non-centered case?
# ============================================================
print("\n" + "="*60)
print("CHECK: n=3 non-centered superadditivity")
print("="*60)

violations = 0
total_valid = 0
for trial in range(1000):
    roots_p = np.sort(np.random.randn(3) * 2)
    roots_q = np.sort(np.random.randn(3) * 2)
    # NOT centering

    min_gap_p = min(np.diff(roots_p))
    min_gap_q = min(np.diff(roots_q))
    if min_gap_p < 0.2 or min_gap_q < 0.2:
        continue

    try:
        roots_r = mss_convolve(roots_p, roots_q)
        phi_p = phi_n(roots_p)
        phi_q = phi_n(roots_q)
        phi_r = phi_n(roots_r)

        if phi_p < 1e-10 or phi_q < 1e-10 or phi_r < 1e-10:
            continue

        total_valid += 1
        lhs = 1/phi_r
        rhs = 1/phi_p + 1/phi_q

        if lhs < rhs - 1e-6:
            violations += 1
            print(f"  VIOLATION: 1/Phi(r)={lhs:.6f} < 1/Phi(p)+1/Phi(q)={rhs:.6f}")
    except:
        pass

print(f"  n=3 non-centered: {violations} violations out of {total_valid} valid trials")

# ============================================================
# Additional: Does Phi_n depend on kappa_1? (translation invariance)
# ============================================================
print("\n" + "="*60)
print("CHECK: Phi_n translation invariance (kappa_1 independence)")
print("="*60)

for shift in [0, 1, -2, 10, -100]:
    roots = np.array([1.0, 3.0, 5.0])
    roots_shifted = roots + shift
    print(f"  Shift={shift}: Phi_3(original)={phi_n(roots):.10f}, "
          f"Phi_3(shifted)={phi_n(roots_shifted):.10f}")
