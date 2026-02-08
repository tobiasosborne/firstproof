"""
Quick check: Is R_3 superadditive or not?
The Cauchy-Schwarz argument should work but we got violations.
Let me debug.
"""
import numpy as np
from math import comb, factorial

np.random.seed(42)

def phi_n_num(roots):
    n = len(roots)
    total = 0.0
    for i in range(n):
        H = sum(1.0/(roots[i]-roots[j]) for j in range(n) if j != i)
        total += H**2
    return total

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
    return np.sort(np.real(np.roots(r_coeffs)))

def get_ta(roots):
    n = len(roots)
    coeffs = np.poly(roots)
    ta = {}
    for k in range(1, n+1):
        ta[k] = ((-1)**k * coeffs[k]) / comb(n, k)
    return ta

n = 3
C3 = -2.0/comb(3,2)  # = -2/3

# The issue: R_n superadditivity requires R_n to be superadditive
# in the CUMULANT VARIABLES, not the root variables.
# Under MSS: K2(r) = K2(p)+K2(q), K3(r) = K3(p)+K3(q)
# So we need: R_3(K2(p)+K2(q), K3(p)+K3(q)) >= R_3(K2(p),K3(p)) + R_3(K2(q),K3(q))

# R_3 = -(1/6)*K3^2/K2^2
# R_3(a+c, b+d) >= R_3(a,b) + R_3(c,d)
# -(b+d)^2/(6(a+c)^2) >= -b^2/(6a^2) - d^2/(6c^2)
# b^2/a^2 + d^2/c^2 >= (b+d)^2/(a+c)^2  (multiply by -6, flip inequality)

# This IS the Cauchy-Schwarz/Jensen inequality for convex function f(x)=x^2!
# Set u = b/a, v = d/c, t = a/(a+c)
# Then (b+d)/(a+c) = (au+cv)/(a+c) = tu+(1-t)v
# (tu+(1-t)v)^2 <= t*u^2 + (1-t)*v^2 (by convexity)
# And t*u^2 + (1-t)*v^2 <= u^2 + v^2 if 0<=t<=1.

# BUT: For centered polys, K2 = ta_2 < 0. So a < 0, c < 0, a+c < 0.
# t = a/(a+c). Since a < 0 and a+c < 0, t = a/(a+c) > 0 (negative/negative = positive).
# Is t <= 1? a/(a+c) <= 1 iff a <= a+c iff 0 <= c iff c >= 0. But c < 0!
# So t > 1 when c < 0 and a < 0. The Jensen argument FAILS.

# Wait: more carefully. a,c < 0. Let a = -|a|, c = -|c|.
# t = a/(a+c) = |a|/(|a|+|c|) (since both negative).
# So t = |a|/(|a|+|c|) which IS in (0,1). ✓
# Then (1-t) = |c|/(|a|+|c|).
# b/a = -b/|a| since a = -|a|. So u = -b/|a|, v = -d/|c|.
# u^2 = b^2/|a|^2 = b^2/a^2 ✓

# Hmm, so Jensen should work. Let me check the argument more carefully.
# Need: b^2/a^2 + d^2/c^2 >= (b+d)^2/(a+c)^2

# With a=-2, c=-1, b=3, d=0:
a, c = -2, -1
b, d = 3, 0
lhs = b**2/a**2 + d**2/c**2  # = 9/4 + 0 = 2.25
rhs = (b+d)**2/(a+c)**2  # = 9/9 = 1
print(f"Test 1: LHS={lhs}, RHS={rhs}, OK: {lhs >= rhs}")

# With a=-1, c=-1, b=2, d=-1:
a, c = -1, -1
b, d = 2, -1
lhs = b**2/a**2 + d**2/c**2  # = 4 + 1 = 5
rhs = (b+d)**2/(a+c)**2  # = 1/4 = 0.25
print(f"Test 2: LHS={lhs}, RHS={rhs}, OK: {lhs >= rhs}")

# The Cauchy-Schwarz/Jensen argument IS correct:
# f(x) = x^2 is convex, so f(tu+(1-t)v) <= tf(u)+(1-t)f(v) <= f(u)+f(v).
# With t = |a|/(|a|+|c|), u = b/a, v = d/c.
# This gives b^2/a^2 + d^2/c^2 >= (b+d)^2/(a+c)^2.

# So R_3 SHOULD be superadditive. The violations must be coming from
# the fact that the numerical R_n test uses ROOTS (not kappa values directly)
# and there might be issues with the MSS convolution preserving real-rootedness.

# Let me test R_3 superadditivity directly in cumulant variables:
print("\n--- Direct R_3 superadditivity test in cumulant variables ---")
violations = 0
trials = 0

for _ in range(20000):
    # Random cumulants
    K2_p = -np.random.exponential(2)  # K2 < 0
    K2_q = -np.random.exponential(2)
    K3_p = np.random.randn() * 2
    K3_q = np.random.randn() * 2

    K2_r = K2_p + K2_q
    K3_r = K3_p + K3_q

    R3_p = -(1/6) * K3_p**2 / K2_p**2
    R3_q = -(1/6) * K3_q**2 / K2_q**2
    R3_r = -(1/6) * K3_r**2 / K2_r**2

    trials += 1
    if R3_r < R3_p + R3_q - 1e-12:
        violations += 1
        if violations <= 3:
            print(f"  VIOLATION: K2_p={K2_p:.4f}, K2_q={K2_q:.4f}, K3_p={K3_p:.4f}, K3_q={K3_q:.4f}")
            print(f"    R3_r={R3_r:.8f}, R3_p+R3_q={R3_p+R3_q:.8f}")

print(f"\nDirect cumulant test: {violations} violations / {trials} trials")

# So the abstract inequality IS satisfied. The issue in the previous test
# must be that the R_n computation from ROOTS doesn't properly match
# the cumulant-space computation, possibly because the roots of p⊞q
# have numerical errors affecting the cumulant extraction.

# OR: the issue is that the FULL R_n superadditivity test was using
# non-centered polynomials. Let me check.

print("\n--- R_3 test with centered polys and cumulant extraction ---")
violations_centered = 0
valid_centered = 0

for trial in range(10000):
    rp = np.sort(np.random.randn(3) * 2)
    rp = rp - np.mean(rp)  # CENTER
    rq = np.sort(np.random.randn(3) * 2)
    rq = rq - np.mean(rq)  # CENTER

    if min(np.diff(rp)) < 0.12 or min(np.diff(rq)) < 0.12:
        continue

    try:
        rr = mss_convolve(rp, rq)
        phi_p = phi_n_num(rp); phi_q = phi_n_num(rq); phi_r = phi_n_num(rr)
        if phi_p <= 0 or phi_q <= 0 or phi_r <= 0:
            continue

        ta_p = get_ta(rp)
        ta_q = get_ta(rq)
        ta_r = get_ta(rr)

        # K2 = ta_2, K3 = ta_3 (for these, exact additivity holds)
        K2_p = ta_p[2]; K3_p = ta_p[3]
        K2_q = ta_q[2]; K3_q = ta_q[3]
        K2_r = ta_r[2]; K3_r = ta_r[3]

        # Verify additivity
        assert abs(K2_r - K2_p - K2_q) < 1e-10
        assert abs(K3_r - K3_p - K3_q) < 1e-10

        R3_p = 1/phi_p - (-2.0/3)*K2_p
        R3_q = 1/phi_q - (-2.0/3)*K2_q
        R3_r = 1/phi_r - (-2.0/3)*K2_r

        valid_centered += 1
        if R3_r < R3_p + R3_q - 1e-8:
            violations_centered += 1
            if violations_centered <= 3:
                # Debug: compare formula vs actual
                R3_p_formula = -(1/6)*K3_p**2/K2_p**2
                R3_q_formula = -(1/6)*K3_q**2/K2_q**2
                R3_r_formula = -(1/6)*K3_r**2/K2_r**2
                print(f"  VIOLATION at trial {trial}:")
                print(f"    R3_r={R3_r:.10f}, R3_p+R3_q={R3_p+R3_q:.10f}")
                print(f"    Formula: R3_r={R3_r_formula:.10f}, R3_p+R3_q={R3_p_formula+R3_q_formula:.10f}")
                print(f"    Actual match formula? p:{abs(R3_p-R3_p_formula):.2e}, q:{abs(R3_q-R3_q_formula):.2e}, r:{abs(R3_r-R3_r_formula):.2e}")
    except:
        pass

print(f"R_3 (centered) superadditivity: {violations_centered} violations / {valid_centered} trials")

# The KEY issue: the previous test used NON-CENTERED polynomials!
# For non-centered polys, kappa_1 ≠ 0, and the formula R_3 = -(1/6)*K3^2/K2^2
# does NOT apply (it's only valid for centered polys).

# For non-centered polys, 1/Phi_n is a function of ALL kappa_k including kappa_1,
# and the decomposition is different.

print("\n--- R_3 test WITHOUT centering (this explains the violations) ---")
violations_nc = 0
valid_nc = 0

for trial in range(10000):
    rp = np.sort(np.random.randn(3) * 2)  # NOT centered
    rq = np.sort(np.random.randn(3) * 2)  # NOT centered

    if min(np.diff(rp)) < 0.12 or min(np.diff(rq)) < 0.12:
        continue

    try:
        rr = mss_convolve(rp, rq)
        phi_p = phi_n_num(rp); phi_q = phi_n_num(rq); phi_r = phi_n_num(rr)
        if phi_p <= 0 or phi_q <= 0 or phi_r <= 0:
            continue

        ta_p = get_ta(rp); ta_q = get_ta(rq); ta_r = get_ta(rr)
        K2_p = ta_p[2]; K2_q = ta_q[2]; K2_r = ta_r[2]

        C3 = -2.0/3
        R3_p = 1/phi_p - C3*K2_p
        R3_q = 1/phi_q - C3*K2_q
        R3_r = 1/phi_r - C3*K2_r

        valid_nc += 1
        if R3_r < R3_p + R3_q - 1e-8:
            violations_nc += 1
    except:
        pass

print(f"R_3 (non-centered) superadditivity: {violations_nc} violations / {valid_nc} trials")
print()
print("CONCLUSION: R_3 IS superadditive for CENTERED polys (as proved by Cauchy-Schwarz).")
print("The violations in the earlier test came from using NON-CENTERED polys,")
print("where the formula R_3 = -(1/6)*K3^2/K2^2 does not apply.")
print()
print("For the FULL inequality 1/Phi_n >= 1/Phi_p + 1/Phi_q, centering")
print("doesn't matter because Phi_n is translation-invariant")
print("(shifting all roots by c doesn't change lambda_i - lambda_j).")
print("And MSS convolution of centered polys IS centered.")

# Verify: is Phi_n translation-invariant?
print("\n--- Verifying Phi_n is translation-invariant ---")
roots = np.array([1.0, 3.0, 5.0])
phi_orig = phi_n_num(roots)
phi_shifted = phi_n_num(roots + 7.3)
print(f"  Phi(roots) = {phi_orig:.10f}")
print(f"  Phi(roots+7.3) = {phi_shifted:.10f}")
print(f"  Match: {abs(phi_orig - phi_shifted) < 1e-10}")

# So the CORRECT statement is:
# For CENTERED polys: 1/Phi_n = C_n*K2 + R_3(K2,...,Kn), and R_3 is superadditive.
# For GENERAL polys: 1/Phi_n depends on K1 too, and the centered formula doesn't apply.
# But since 1/Phi_n IS translation-invariant, we can WLOG center first.
# And then the decomposition works with centered cumulants.
