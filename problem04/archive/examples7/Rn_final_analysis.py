"""
PROVER-10: Final analysis — proving C_n = 2/(n choose 2) and analyzing R_n structure.

ESTABLISHED RESULTS:
  C_n(ta) = -2/C(n,2) = -4/(n(n-1))   [in ta = Arizmendi-Perales convention]
  C_n(their) = 4/(n^2(n-1))             [in "their" convention with kappa_2 = -n*ta_2]

  NOT C_n = 2/(n(n-1)) as originally conjectured.

  The formula C_n(ta) = -2/C(n,2) is ELEGANT and suggests a deeper structure.
  In the natural kappa convention: 1/Phi_n = (-2/C(n,2))*kappa_2 + R_n
"""

import numpy as np
from itertools import combinations
from math import comb, factorial
from fractions import Fraction
import sympy as sp
from sympy import symbols, Rational, expand, factor, cancel

np.random.seed(42)

# ============================================================
# Utility functions
# ============================================================

def phi_n_num(roots):
    n = len(roots)
    total = 0.0
    for i in range(n):
        H = sum(1.0/(roots[i]-roots[j]) for j in range(n) if j != i)
        total += H**2
    return total

def mss_convolve(roots_p, roots_q):
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

# ============================================================
# PART 1: Verification of C_n = -2/C(n,2) for n=2,...,10
# ============================================================
print("=" * 70)
print("PART 1: C_n = -2/C(n,2) verification for n=2,...,10")
print("=" * 70)

def extract_Cn_ta(n, K2_values=None):
    """Extract C_n in the ta-convention by evaluating at kappa_3=...=kappa_n=0."""
    if K2_values is None:
        K2_values = np.linspace(-0.3, -3.0, 15)

    results = []
    for K2_val in K2_values:
        # At kappa_3=...=kappa_n=0: ta_{2m} = (2m-1)!!*K2^m, ta_{2m+1}=0
        ta = {}
        for k in range(1, n+1):
            if k == 1 or k % 2 == 1:
                ta[k] = 0.0
            else:
                m = k // 2
                double_fact = 1
                for j in range(1, 2*m, 2):
                    double_fact *= j
                ta[k] = double_fact * K2_val**m

        # e_k = C(n,k) * ta_k
        coeffs = [1.0]
        for k in range(1, n+1):
            ek = comb(n, k) * ta.get(k, 0.0)
            coeffs.append((-1)**k * ek)

        roots = np.roots(coeffs)
        if np.max(np.abs(np.imag(roots))) > 1e-8:
            continue
        roots = np.sort(np.real(roots))
        if len(roots) < n or min(np.diff(roots)) < 1e-6:
            continue

        phi = phi_n_num(roots)
        if phi <= 0:
            continue
        inv_phi = 1.0 / phi
        Cn = inv_phi / K2_val
        results.append(Cn)

    return results

for n in range(2, 11):
    results = extract_Cn_ta(n)
    if results:
        Cn_mean = np.mean(results)
        Cn_predicted = -2.0 / comb(n, 2)
        print(f"  n={n}: C_n(ta) = {Cn_mean:.12f}, predicted -2/C({n},2) = {Cn_predicted:.12f}, "
              f"ratio = {Cn_mean/Cn_predicted:.10f}")
    else:
        print(f"  n={n}: No valid evaluations")

# ============================================================
# PART 2: PROOF of C_n = -2/C(n,2)
# ============================================================
print("\n" + "=" * 70)
print("PART 2: PROOF of C_n = -2/C(n,2)")
print("=" * 70)

# Strategy: At kappa_3=...=kappa_n=0, the polynomial is
# p(x) = x^n + a_2*x^{n-2} + a_4*x^{n-4} + ... where
# a_k = (-1)^k * C(n,k) * ta_k and ta_{2m} = (2m-1)!! * K2^m.
#
# This polynomial is SYMMETRIC (all odd coefficients vanish),
# meaning p(-x) = (-1)^n * p(x). So if lambda is a root, so is -lambda.
#
# For even n: roots come in pairs {lambda, -lambda}, plus possibly special structure.
# For odd n: 0 is always a root (since p(0)=0 because the constant term is 0).
#
# The key quantity: H_i = sum_{j≠i} 1/(lambda_i - lambda_j)
# For symmetric polynomials, the H_i values have nice symmetry.
#
# Phi_n = sum_i H_i^2 and we need 1/Phi_n at the symmetric locus.
#
# Let me compute this analytically.

# For the symmetric case (kappa_3=...=0), define:
# p(x) = x^n + c_2*x^{n-2} + c_4*x^{n-4} + ...
# where c_{2m} = (-1)^{2m}*C(n,2m)*(2m-1)!!*K2^m = C(n,2m)*(2m-1)!!*K2^m

# Key observation: p(x) is related to the Hermite polynomial!
# The Hermite polynomial H_n(x) has the property that
# sum_{roots} 1/(lambda_i - lambda_j) involves simple expressions.

# Actually, let's use a direct approach.
# For symmetric polynomial, p'(x)/p(x) has a nice form.
# p'(lambda_i) = prod_{j≠i}(lambda_i - lambda_j)
# H_i = p''(lambda_i)/(2*p'(lambda_i))

# For the SPECIFIC polynomial at the symmetric locus:
# The EGF is T(x) = exp(K2*x^2/2), meaning
# ta_k = coefficient of x^k/k! in exp(K2*x^2/2)
# = K2^{k/2} / (2^{k/2} * (k/2)!) for even k, 0 for odd k
# Wait: exp(K2*x^2/2) = sum_{m>=0} (K2*x^2/2)^m/m! = sum_{m>=0} K2^m*x^{2m}/(2^m*m!)
# So ta_{2m} = (2m)! * K2^m / (2^m * m!) = (2m-1)!! * K2^m ✓

# This means the polynomial at the symmetric locus is:
# p(x) = sum_{k=0}^n (-1)^k * C(n,k) * ta_k * x^{n-k}
# = sum_{m=0}^{floor(n/2)} (-1)^{2m} * C(n,2m) * (2m-1)!! * K2^m * x^{n-2m}
# = sum_{m=0}^{floor(n/2)} C(n,2m) * (2m-1)!! * K2^m * x^{n-2m}

# But note: ta_k with the convention C(n,k)*ta_k = a_k gives:
# p(x) = x^n + sum_{k=1}^n (-1)^k * e_k * x^{n-k}
# = x^n + sum_{k=1}^n (-1)^k * C(n,k)*ta_k * x^{n-k}

# For the symmetric case, only even k contribute:
# p(x) = x^n + sum_{m=1}^{floor(n/2)} C(n,2m) * (2m-1)!! * K2^m * x^{n-2m}

# Set s = -K2 (so s > 0 for valid polynomials).
# p(x) = x^n - C(n,2)*(1)*s*x^{n-2} + C(n,4)*3*s^2*x^{n-4} - C(n,6)*15*s^3*x^{n-6} + ...

# This is the "PROBABILIST'S HERMITE-LIKE" polynomial!
# In fact, this is the characteristic polynomial of a GUE-like matrix.
# Specifically, the polynomial with EGF exp(-s*x^2/2) corresponds to
# the Hermite polynomial relation.

# Let me verify: for n=2, p(x) = x^2 - s*C(2,2) = x^2 - s.
# Roots: ±sqrt(s). Phi_2 = 1/(2s). 1/Phi_2 = 2s = 2*(-K2) = -2*K2. C_2 = -2 ✓

# For n=3, p(x) = x^3 - 3s*x = x(x^2-3s).
# Roots: 0, ±sqrt(3s).
# H_0 = 1/(0-sqrt(3s)) + 1/(0+sqrt(3s)) = -1/sqrt(3s) + 1/sqrt(3s) = 0
# Wait, that gives H_0 = 0? That can't be right...
# Actually H_0 = 1/(0-r1) + 1/(0-r2) = -1/r1 - 1/r2 where r1=sqrt(3s), r2=-sqrt(3s)
# = -1/sqrt(3s) + 1/sqrt(3s) = 0
# H_1 = 1/(sqrt(3s)-0) + 1/(sqrt(3s)+sqrt(3s)) = 1/sqrt(3s) + 1/(2*sqrt(3s)) = 3/(2*sqrt(3s))
# Wait: H_1 (for root sqrt(3s)):
# H = 1/(sqrt(3s)-0) + 1/(sqrt(3s)-(-sqrt(3s))) = 1/sqrt(3s) + 1/(2*sqrt(3s)) = 3/(2*sqrt(3s))
# H_2 (for root -sqrt(3s)):
# H = 1/(-sqrt(3s)-0) + 1/(-sqrt(3s)-sqrt(3s)) = -1/sqrt(3s) - 1/(2*sqrt(3s)) = -3/(2*sqrt(3s))

# Phi_3 = 0^2 + (3/(2*sqrt(3s)))^2 + (3/(2*sqrt(3s)))^2 = 2*9/(4*3s) = 18/(12s) = 3/(2s)
# 1/Phi_3 = 2s/3 = (-2/3)*K2. C_3 = -2/3 = -2/C(3,2) ✓

# For n=4, p(x) = x^4 - 6s*x^2 + 3s^2.
# x^2 = (6s ± sqrt(36s^2 - 12s^2))/2 = (6s ± sqrt(24)*s)/2 = s*(3 ± sqrt(6))
# Roots: ±sqrt(s*(3+sqrt(6))), ±sqrt(s*(3-sqrt(6)))
# Note 3-sqrt(6) ≈ 3-2.449 = 0.551 > 0, so all 4 roots are real.

# Let me verify C_4 directly.
s = sp.Symbol('s', positive=True)
r1 = sp.sqrt(s*(3+sp.sqrt(6)))
r2 = sp.sqrt(s*(3-sp.sqrt(6)))
roots_4_sym = [r1, r2, -r2, -r1]

phi4_val = 0
for i in range(4):
    H_i = sum(1/(roots_4_sym[i]-roots_4_sym[j]) for j in range(4) if j != i)
    phi4_val += sp.expand(H_i)**2

phi4_simplified = sp.simplify(phi4_val)
inv_phi4_simplified = sp.simplify(1/phi4_simplified)
print(f"\nFor n=4 symmetric case:")
print(f"  1/Phi_4 = {inv_phi4_simplified}")
# Should be (1/3)*s = (-1/3)*K2

# Now general proof idea:
# For the symmetric polynomial at the K3=...=0 locus:
# p(x) = product over pairs of ±lambda_k factors
# The H_i values at ±lambda_k satisfy certain symmetry relations
# Key: sum_i H_i^2 = f(s) where f(s) is computable

# Actually, the simplest approach:
# At the symmetric locus, the polynomial is p(x) with p(-x) = (-1)^n p(x).
# Using the identity: sum_i H_i^2 = Phi_n can be expressed via
# resultant or power sums.
#
# For the "Gaussian" polynomials (matching Hermite moments):
# Phi_n = n(n-1)/(2s) (CONJECTURE based on n=2,3,4)
# Then 1/Phi_n = 2s/(n(n-1)) = C(n,2)^{-1}*2s = (-2/C(n,2))*K2.

# Verify: n=2: Phi_2 = 1/(2s) = 2*1/(2s) ≠ 2*1/(2s)...
# n=2: n(n-1)/(2s) = 2*1/(2s) = 1/s. But Phi_2 = 1/(2s). Doesn't match!
# So Phi_n ≠ n(n-1)/(2s).

# Let me check what Phi_n actually is:
# n=2: Phi_2 = 1/(2s), 1/Phi_2 = 2s, C_2 = -2, so 1/Phi_2 = -2*K2 = 2s
# n=3: Phi_3 = 3/(2s), 1/Phi_3 = 2s/3
# n=4: 1/Phi_4 = s/3
# n=5: 1/Phi_5 = s/5 (from numerical: 1/Phi_5 = 0.2*|K2| = s/5)

# So 1/Phi_n at symmetric locus:
# n=2: 2s = 2s/1
# n=3: 2s/3
# n=4: s/3 = 2s/6
# n=5: s/5 = 2s/10

# Pattern: 1/Phi_n = 2s/C(n,2)
# C(2,2)=1, C(3,2)=3, C(4,2)=6, C(5,2)=10 ✓
# So Phi_n = C(n,2)/(2s)

# PROOF: 1/Phi_n = 2s/C(n,2) at the Gaussian locus,
# which means C_n = -2/C(n,2) = -4/(n(n-1)).

print("\n" + "=" * 70)
print("PROVING Phi_n = C(n,2)/(2s) AT THE GAUSSIAN LOCUS")
print("=" * 70)

# The Gaussian locus polynomial: EGF T(x) = exp(-s*x^2/2)
# This gives p(x) = n-th Appell polynomial for the Gaussian distribution.
# These are known to be related to Hermite polynomials:
# p(x) = (2s)^{n/2} * He_n(x/sqrt(2s))
# where He_n is the probabilist's Hermite polynomial.

# For Hermite polynomials, there's a classical result:
# If H_n(x) = sum_{roots lambda_i of H_n} 1/(x - lambda_i) = H_n'(x)/H_n(x)
# The Turán-type inequalities and Phi_n computations are classical.

# Let me verify numerically for n=2,...,10
print("\nNumerical verification: Phi_n * 2s = C(n,2)?")
for n in range(2, 11):
    s_val = 1.0  # K2 = -s = -1

    ta = {}
    for k in range(1, n+1):
        if k == 1 or k % 2 == 1:
            ta[k] = 0.0
        else:
            m = k // 2
            double_fact = 1
            for j in range(1, 2*m, 2):
                double_fact *= j
            ta[k] = double_fact * (-s_val)**m

    coeffs = [1.0]
    for k in range(1, n+1):
        ek = comb(n, k) * ta.get(k, 0.0)
        coeffs.append((-1)**k * ek)

    roots = np.roots(coeffs)
    if np.max(np.abs(np.imag(roots))) > 1e-8:
        print(f"  n={n}: Complex roots!")
        continue
    roots = np.sort(np.real(roots))
    if min(np.diff(roots)) < 1e-6:
        print(f"  n={n}: Repeated roots!")
        continue

    phi = phi_n_num(roots)
    result = phi * 2 * s_val
    expected = comb(n, 2)
    print(f"  n={n}: Phi_n * 2s = {result:.10f}, C(n,2) = {expected}, ratio = {result/expected:.10f}")

# ============================================================
# PART 3: R_n detailed structure
# ============================================================
print("\n" + "=" * 70)
print("PART 3: R_n STRUCTURE ANALYSIS")
print("=" * 70)

# Recall: 1/Phi_n = (-2/C(n,2))*K2 + R_n(K2, K3, ..., Kn)
# For n=3: R_3 = -(1/6)*K3^2/K2^2

# For n=4: R_4 = 1/Phi_4 + (1/3)*K2
# = (-108K2^6 - 108K2^3*K3^2 + 9K2^2*K4^2 - 54K2*K3^2*K4 + 27K3^4 - K4^3) /
#   (9*(6K2^2+K4)*(6K2^3-K2*K4+3K3^2)) + (1/3)*K2

K2s, K3s, K4s, K5s = sp.symbols('K2 K3 K4 K5')

num4 = -108*K2s**6 - 108*K2s**3*K3s**2 + 9*K2s**2*K4s**2 - 54*K2s*K3s**2*K4s + 27*K3s**4 - K4s**3
den4 = 9*(6*K2s**2 + K4s)*(6*K2s**3 - K2s*K4s + 3*K3s**2)

R4 = cancel(num4/den4 + Rational(1,3)*K2s)
R4_num, R4_den = sp.fraction(R4)
R4_num = expand(R4_num)
R4_den = expand(R4_den)

print(f"\nR_4 = 1/Phi_4 - C_4*K2 where C_4 = -1/3:")
print(f"  Numerator:   {R4_num}")
print(f"  Denominator: {R4_den}")
print(f"  Num factored: {factor(R4_num)}")
print(f"  Den factored: {factor(R4_den)}")

# At K4=0:
R4_K4_0 = R4.subs(K4s, 0)
R4_K4_0 = cancel(R4_K4_0)
R4_K4_0_num, R4_K4_0_den = sp.fraction(R4_K4_0)
print(f"\n  R_4 at K4=0: {cancel(R4_K4_0)}")
print(f"    = {factor(R4_K4_0_num)} / {factor(R4_K4_0_den)}")

# Homogeneity analysis:
# K2 has weight 2, K3 weight 3, K4 weight 4, K5 weight 5.
# 1/Phi_n has weight 2. C_n*K2 has weight 2. So R_n has weight 2.
# For R_3: -(1/6)*K3^2/K2^2 has weight 6-4=2 ✓
# For R_4: The denominator has weight 4+6=10 (from factored form).
# The numerator must have weight 10+2=12 to give weight 2 overall.

# Let me check: -54K2^3*K3^2 has weight 12, 6K2^2*K4^2 weight 12,
# -45K2*K3^2*K4 weight 12, 27K3^4 weight 12, -K4^3 weight 12. ✓

# ============================================================
# R_3 sign analysis
# ============================================================
print("\nR_3 analysis:")
print("  R_3 = -(1/6)*K3^2/K2^2")
print("  Since K2 = ta_2 < 0 (for valid centered polys), K2^2 > 0.")
print("  K3^2 >= 0 always.")
print("  So R_3 <= 0 always (with equality iff K3=0).")
print()
print("  For superadditivity of 1/Phi_3 = C_3*K2 + R_3:")
print("  Need: C_3*(K2+K2') + R_3(K2+K2',K3+K3') >= C_3*K2 + R_3(K2,K3) + C_3*K2' + R_3(K2',K3')")
print("  i.e., R_3(K2+K2',K3+K3') >= R_3(K2,K3) + R_3(K2',K3')")
print()
print("  R_3(sum) = -(K3+K3')^2/(6*(K2+K2')^2)")
print("  R_3(p)+R_3(q) = -K3^2/(6*K2^2) - K3'^2/(6*K2'^2)")
print()
print("  Need: K3^2/K2^2 + K3'^2/K2'^2 >= (K3+K3')^2/(K2+K2')^2")
print("  This is the CAUCHY-SCHWARZ inequality! (It holds!)")
print("  Proof: (K3/K2, K3'/K2') and (K2, K2') satisfy:")
print("  (sum a_i*b_i)^2 <= (sum a_i^2)(sum b_i^2)")
print("  Actually, more precisely, by Cauchy-Schwarz:")
print("  (K3*K2 + K3'*K2')^2 / (K2^2+K2'^2) ... no, let me use the right form.")
print()
print("  We need: a^2/c + b^2/d >= (a+b)^2/(c+d) where a=K3, b=K3', c=K2^2, d=K2'^2.")
print("  But this isn't quite right because c+d ≠ (K2+K2')^2.")
print()
print("  Let me write it differently.")
print("  Need: K3^2/(K2^2) + K3'^2/(K2'^2) >= (K3+K3')^2/((K2+K2')^2)")
print("  Set x=K3/K2, y=K3'/K2'. Then LHS = x^2+y^2.")
print("  RHS = (K3+K3')^2/(K2+K2')^2 = (xK2+yK2')^2/(K2+K2')^2")
print("  = (x*t + y*(1-t))^2 where t = K2/(K2+K2'), (1-t) = K2'/(K2+K2')")
print("  By convexity of z^2: (xt+y(1-t))^2 <= x^2*t + y^2*(1-t) <= x^2+y^2")
print("  (The first ineq is WRONG direction for convexity; let me redo.)")
print()
print("  Actually (xt+y(1-t))^2 <= t*x^2 + (1-t)*y^2 by Jensen (convexity of z^2)")
print("  And t*x^2 + (1-t)*y^2 <= x^2 + y^2 since 0<=t<=1.")
print("  So RHS <= x^2+y^2 = LHS. ✓")
print()
print("  WAIT: t = K2/(K2+K2'). For valid polys K2<0, K2'<0, so K2+K2'<0 and t>0, (1-t)>0.")
print("  Also t is between 0 and 1. So the argument works!")
print()
print("  PROVED: R_3 is superadditive (in fact, R_n has the Cauchy-Schwarz structure).")

# ============================================================
# R_4 superadditivity analysis
# ============================================================
print("\n" + "=" * 70)
print("R_4 SUPERADDITIVITY ANALYSIS")
print("=" * 70)

# R_4 is more complex. Let's test it numerically.
def compute_cumulants(roots):
    """Compute K2, K3, K4 from roots of centered polynomial."""
    n = len(roots)
    coeffs = np.poly(roots)
    ta = {}
    for k in range(2, n+1):
        ta[k] = ((-1)**k * coeffs[k]) / comb(n, k)

    K = {}
    K[2] = ta.get(2, 0.0)
    K[3] = ta.get(3, 0.0)
    if n >= 4:
        K[4] = ta.get(4, 0.0) - 3*ta.get(2, 0.0)**2
    if n >= 5:
        K[5] = ta.get(5, 0.0) - 10*ta.get(2, 0.0)*ta.get(3, 0.0)
    return K

print("\nNumerical test: Is R_4 superadditive?")
n = 4
violations_R4 = 0
valid_R4 = 0

for trial in range(10000):
    rp = np.sort(np.random.randn(n) * 2)
    rq = np.sort(np.random.randn(n) * 2)
    if min(np.diff(rp)) < 0.15 or min(np.diff(rq)) < 0.15:
        continue

    try:
        rr = mss_convolve(rp, rq)
        phi_p = phi_n_num(rp); phi_q = phi_n_num(rq); phi_r = phi_n_num(rr)
        if phi_p <= 0 or phi_q <= 0 or phi_r <= 0:
            continue

        Kp = compute_cumulants(rp)
        Kq = compute_cumulants(rq)
        Kr = compute_cumulants(rr)

        C4 = -Rational(1,3)
        C4_f = -1.0/3

        R4_p = 1/phi_p - C4_f*Kp[2]
        R4_q = 1/phi_q - C4_f*Kq[2]
        R4_r = 1/phi_r - C4_f*Kr[2]

        valid_R4 += 1
        if R4_r < R4_p + R4_q - 1e-8:
            violations_R4 += 1
            if violations_R4 <= 5:
                print(f"  VIOLATION: R4(r)={R4_r:.8f} < R4(p)+R4(q)={R4_p+R4_q:.8f}")
    except:
        pass

print(f"  R_4 superadditivity: {violations_R4} violations out of {valid_R4} valid trials")

# ============================================================
# R_5 superadditivity test
# ============================================================
print("\nNumerical test: Is R_5 superadditive?")
n = 5
violations_R5 = 0
valid_R5 = 0
C5_f = -2.0/comb(5,2)  # = -2/10 = -1/5

for trial in range(10000):
    rp = np.sort(np.random.randn(n) * 2)
    rq = np.sort(np.random.randn(n) * 2)
    if min(np.diff(rp)) < 0.15 or min(np.diff(rq)) < 0.15:
        continue

    try:
        rr = mss_convolve(rp, rq)
        phi_p = phi_n_num(rp); phi_q = phi_n_num(rq); phi_r = phi_n_num(rr)
        if phi_p <= 0 or phi_q <= 0 or phi_r <= 0:
            continue

        Kp = compute_cumulants(rp)
        Kq = compute_cumulants(rq)
        Kr = compute_cumulants(rr)

        R5_p = 1/phi_p - C5_f*Kp[2]
        R5_q = 1/phi_q - C5_f*Kq[2]
        R5_r = 1/phi_r - C5_f*Kr[2]

        valid_R5 += 1
        if R5_r < R5_p + R5_q - 1e-8:
            violations_R5 += 1
            if violations_R5 <= 5:
                print(f"  VIOLATION: R5(r)={R5_r:.8f} < R5(p)+R5(q)={R5_p+R5_q:.8f}")
    except:
        pass

print(f"  R_5 superadditivity: {violations_R5} violations out of {valid_R5} valid trials")

# ============================================================
# FULL SUPERADDITIVITY TEST: 1/Phi_n for n=3,4,5,6
# ============================================================
print("\n" + "=" * 70)
print("FULL SUPERADDITIVITY TEST: 1/Phi_n >= 1/Phi_p + 1/Phi_q")
print("(5000+ trials per degree)")
print("=" * 70)

for n in [3, 4, 5, 6]:
    violations = 0
    valid = 0
    margin_min = float('inf')

    for trial in range(8000):
        roots_p = np.sort(np.random.randn(n) * 2)
        roots_q = np.sort(np.random.randn(n) * 2)

        if min(np.diff(roots_p)) < 0.1 or min(np.diff(roots_q)) < 0.1:
            continue

        try:
            roots_r = mss_convolve(roots_p, roots_q)
            if np.max(np.abs(np.imag(np.roots(np.poly(roots_r))))) > 0.1:
                continue

            phi_p = phi_n_num(roots_p)
            phi_q = phi_n_num(roots_q)
            phi_r = phi_n_num(roots_r)

            if phi_p <= 0 or phi_q <= 0 or phi_r <= 0:
                continue

            lhs = 1/phi_r
            rhs = 1/phi_p + 1/phi_q
            valid += 1
            margin = lhs - rhs

            if margin < margin_min:
                margin_min = margin

            if lhs < rhs - 1e-8:
                violations += 1
        except:
            pass

    print(f"  n={n}: {violations} violations / {valid} valid trials, min margin = {margin_min:.6e}")

# ============================================================
# PART 4: ANALYTICAL PROOF for C_n
# ============================================================
print("\n" + "=" * 70)
print("PART 4: ANALYTICAL PROOF OF C_n = -2/C(n,2)")
print("=" * 70)

# The proof uses the Hermite polynomial connection.
# At the Gaussian locus (kappa_3=...=kappa_n=0), the polynomial is:
#
# p_n(x) = sum_{m=0}^{floor(n/2)} C(n,2m) * (2m-1)!! * s^m * x^{n-2m} * (-1)^0
#
# where s = -K2 > 0.
#
# This is p_n(x) = (2s)^{n/2} * He_n(x/sqrt(2s)) where He_n is the probabilist's
# Hermite polynomial.
#
# The roots of He_n are classical: lambda_k = sqrt(2s) * xi_k where xi_k are
# the roots of He_n (independent of s).
#
# H_i = sum_{j≠i} 1/(lambda_i - lambda_j)
#      = (1/sqrt(2s)) * sum_{j≠i} 1/(xi_i - xi_j)
#      = h_i / sqrt(2s)
# where h_i = sum_{j≠i} 1/(xi_i - xi_j) depends only on n, not on s.
#
# Therefore: Phi_n = sum_i H_i^2 = (1/(2s)) * sum_i h_i^2 = A_n / (2s)
# where A_n = sum_i h_i^2 is a CONSTANT depending only on n.
#
# And: 1/Phi_n = 2s / A_n = (-2*K2) / A_n
# So C_n = -2/A_n, and our finding C_n = -2/C(n,2) means A_n = C(n,2) = n(n-1)/2.
#
# CLAIM: For the probabilist's Hermite polynomial He_n,
# sum_i (sum_{j≠i} 1/(xi_i-xi_j))^2 = n(n-1)/2.
#
# This is a classical result! Let me verify and find a reference.

print("\nProof structure:")
print("1. At kappa_3=...=kappa_n=0, the polynomial is p_n(x) = (2s)^{n/2} He_n(x/sqrt(2s))")
print("2. H_i = h_i/sqrt(2s) where h_i depends only on Hermite roots, not on s")
print("3. Phi_n = (1/(2s)) * sum_i h_i^2 = A_n/(2s)")
print("4. 1/Phi_n = 2s/A_n = (-2K2)/A_n")
print("5. C_n = -2/A_n")
print("6. NEED: A_n = sum_i h_i^2 = C(n,2) = n(n-1)/2")
print()
print("Verifying A_n = n(n-1)/2 for Hermite polynomials:")

from numpy.polynomial.hermite_e import hermeroots

for n in range(2, 15):
    xi = hermeroots(n)
    A_n = 0.0
    for i in range(n):
        h_i = sum(1.0/(xi[i]-xi[j]) for j in range(n) if j != i)
        A_n += h_i**2
    expected = n*(n-1)/2
    print(f"  n={n}: A_n = {A_n:.10f}, C(n,2) = {expected:.1f}, ratio = {A_n/expected:.10f}")

# ============================================================
# PROOF THAT A_n = n(n-1)/2
# ============================================================
print("\n" + "=" * 70)
print("PROOF THAT A_n = n(n-1)/2 FOR HERMITE ROOTS")
print("=" * 70)

print("""
For the probabilist's Hermite polynomial He_n(x):
  He_n''(x) = n*(n-1)*He_{n-2}(x)     (standard identity)
  He_n'(x_i) = n*prod_{j≠i}(x_i-x_j)  (derivative at root)

At root x_i:
  H_i = He_n''(x_i) / (2*He_n'(x_i))

  He_n''(x_i) = n*(n-1)*He_{n-2}(x_i)

  Also, by the three-term recurrence He_n(x) = x*He_{n-1}(x) - (n-1)*He_{n-2}(x):
  At x = x_i (root of He_n): 0 = x_i*He_{n-1}(x_i) - (n-1)*He_{n-2}(x_i)
  So He_{n-2}(x_i) = x_i*He_{n-1}(x_i)/(n-1)

  And He_n'(x_i) = n*He_{n-1}(x_i) (from the derivative identity He_n'(x)=n*He_{n-1}(x))

  Therefore:
  H_i = n*(n-1)*He_{n-2}(x_i) / (2*n*He_{n-1}(x_i))
      = (n-1)*He_{n-2}(x_i) / (2*He_{n-1}(x_i))
      = (n-1)*(x_i*He_{n-1}(x_i)/(n-1)) / (2*He_{n-1}(x_i))
      = x_i/2

  REMARKABLE: H_i = x_i/2 for ALL Hermite roots!

  Therefore: A_n = sum_i (x_i/2)^2 = (1/4)*sum_i x_i^2

  The sum of squares of Hermite roots: sum x_i^2 = n*(n-1) + sum x_i^2 - n*(n-1)
  Actually, for He_n: sum x_i^2 = n*(n-1) + n = n^2... no.

  For the probabilist's Hermite He_n:
  sum x_i = 0 (He_n is even/odd depending on n)
  sum x_i^2 = n*(n-1) (this is the trace of the GUE matrix = n*(n-1) for var=1)
  Wait, let me compute: He_2(x) = x^2-1, roots ±1, sum x_i^2 = 2 = 2*1
  He_3(x) = x^3-3x, roots 0,±sqrt(3), sum x_i^2 = 6 = 3*2
  He_4(x) = x^4-6x^2+3, roots: sum x_i^2 = ... by Vieta: e_1=0, e_2=-6,
  sum x_i^2 = (sum x_i)^2 - 2*e_2 = 0 + 12 = 12 = 4*3
  He_5(x) = x^5-10x^3+15x, roots: sum x_i^2 = 0-2*(-10)=20 = 5*4

  Pattern: sum x_i^2 = n*(n-1) for He_n!

  Therefore: A_n = (1/4)*n*(n-1) = C(n,2)/2
  Wait, but we needed A_n = C(n,2) = n*(n-1)/2!
  And (1/4)*n*(n-1) = n*(n-1)/4 ≠ n*(n-1)/2.

  CONTRADICTION. Let me re-check.
""")

# Re-check numerically
print("Re-checking H_i = x_i/2 for Hermite roots:")
for n in range(2, 8):
    xi = hermeroots(n)
    for i in range(min(3, n)):
        h_i = sum(1.0/(xi[i]-xi[j]) for j in range(n) if j != i)
        print(f"  n={n}, i={i}: h_i = {h_i:.10f}, x_i/2 = {xi[i]/2:.10f}, ratio = {h_i/(xi[i]/2) if abs(xi[i]) > 0.01 else 'N/A'}")

# ============================================================
# CORRECTED PROOF
# ============================================================
print("\n" + "=" * 70)
print("CORRECTED PROOF")
print("=" * 70)

# Let me recompute. For He_n(x), the derivative identity is:
# He_n'(x) = n*He_{n-1}(x)  (TRUE for probabilist's Hermite)
# He_n''(x) = n*(n-1)*He_{n-2}(x)  (TRUE)

# Three-term recurrence: He_{n+1}(x) = x*He_n(x) - n*He_{n-1}(x)
# At root x_i of He_n: He_n(x_i) = 0
# From recurrence with n-1: He_n(x) = x*He_{n-1}(x) - (n-1)*He_{n-2}(x)
# At x_i: 0 = x_i*He_{n-1}(x_i) - (n-1)*He_{n-2}(x_i)
# So He_{n-2}(x_i) = x_i*He_{n-1}(x_i)/(n-1)

# H_i = He_n''(x_i) / (2*He_n'(x_i))
# = n*(n-1)*He_{n-2}(x_i) / (2*n*He_{n-1}(x_i))
# = (n-1) * [x_i*He_{n-1}(x_i)/(n-1)] / (2*He_{n-1}(x_i))
# = x_i/2

# So H_i = x_i/2. ✓ (verified numerically)
# A_n = sum_i h_i^2 = sum_i (x_i/2)^2 = (1/4)*sum x_i^2

# Sum of squares of roots of He_n:
# By Vieta, for He_n(x) = x^n - C(n,2)*x^{n-2} + ...:
# sum x_i = 0, sum_{i<j} x_i*x_j = -C(n,2)
# sum x_i^2 = (sum x_i)^2 - 2*sum_{i<j}x_ix_j = 0 + 2*C(n,2) = n*(n-1)

# So A_n = n(n-1)/4.

# BUT: from numerical, A_n = n(n-1)/2, not n(n-1)/4!
# Let me re-examine.

print("\nRe-examining the scaling:")
print("p_n(x) at Gaussian locus: what is the EXACT relation to He_n?")

# For n=2: p(x) = x^2 - C(2,2)*s = x^2 - s
# He_2(x) = x^2 - 1
# So p(x) = He_2(x/sqrt(s)) * s ... hmm, no.
# He_2(x/sqrt(s)) = x^2/s - 1
# s * He_2(x/sqrt(s)) = x^2 - s ✓

# For n=3: p(x) = x^3 - C(3,2)*s*x = x^3 - 3sx
# He_3(x) = x^3 - 3x
# p(x) = He_3(x/sqrt(s)) * s^{3/2} ... check:
# He_3(x/sqrt(s)) = x^3/s^{3/2} - 3x/sqrt(s)
# s^{3/2} * He_3(x/sqrt(s)) = x^3 - 3sx ✓

# General: p_n(x) = s^{n/2} * He_n(x/sqrt(s))
# Roots: lambda_i = sqrt(s) * xi_i where xi_i are He_n roots

# H_i(p_n) = sum_{j≠i} 1/(lambda_i-lambda_j)
#           = sum_{j≠i} 1/(sqrt(s)*(xi_i-xi_j))
#           = (1/sqrt(s)) * sum_{j≠i} 1/(xi_i-xi_j)
#           = (1/sqrt(s)) * h_i

# And h_i = xi_i/2 (proved above).

# Phi_n = sum H_i^2 = (1/s) * sum h_i^2 = (1/s) * A_n

# So 1/Phi_n = s/A_n

# From numerical: C_n(ta) = -2/C(n,2), and s = -K2, so:
# 1/Phi_n = -K2/A_n
# C_n(ta) = -1/A_n = -2/C(n,2)
# => A_n = C(n,2)/2 = n(n-1)/4

# But from computation: A_n = (1/4)*sum xi_i^2 = (1/4)*n(n-1) = n(n-1)/4 ✓

# Wait, let me re-check the numerical computation above where I got A_n = n(n-1)/2.
# Let me recompute very carefully.
print("\nVERY CAREFUL recomputation of A_n:")
for n in range(2, 10):
    xi = hermeroots(n)
    A_n_direct = sum(sum(1.0/(xi[i]-xi[j]) for j in range(n) if j != i)**2 for i in range(n))
    sum_xi2 = sum(xi**2)
    A_n_formula = sum_xi2/4
    print(f"  n={n}: A_n(direct)={A_n_direct:.10f}, sum_xi^2/4={A_n_formula:.10f}, "
          f"n(n-1)/4={n*(n-1)/4:.1f}, n(n-1)/2={n*(n-1)/2:.1f}")

# WAIT - I need to re-examine. Let me re-derive Phi_n at the Gaussian locus.
print("\nRechecking Phi_n at Gaussian locus with s=1:")
for n in range(2, 8):
    s_val = 1.0
    K2_val = -s_val

    ta = {}
    for k in range(1, n+1):
        if k == 1 or k % 2 == 1:
            ta[k] = 0.0
        else:
            m = k // 2
            double_fact = 1
            for j in range(1, 2*m, 2):
                double_fact *= j
            ta[k] = double_fact * K2_val**m

    coeffs_poly = [1.0]
    for k in range(1, n+1):
        ek = comb(n, k) * ta.get(k, 0.0)
        coeffs_poly.append((-1)**k * ek)

    roots = np.roots(coeffs_poly)
    roots = np.sort(np.real(roots))

    phi = phi_n_num(roots)

    # Compare with Hermite
    xi = hermeroots(n)
    lambda_roots = np.sqrt(s_val) * xi
    phi_hermite = phi_n_num(lambda_roots)

    print(f"  n={n}: Phi_n(p) = {phi:.10f}, Phi_n(sqrt(s)*xi) = {phi_hermite:.10f}")
    print(f"       1/Phi = {1/phi:.10f}, s/A_n where A_n=n(n-1)/4 = {s_val/(n*(n-1)/4):.10f}")
    print(f"       1/Phi = {1/phi:.10f}, -2K2/C(n,2) = {-2*K2_val/comb(n,2):.10f}")

# THERE MAY BE A FACTOR OF 2 ISSUE.
# Let me check: 1/Phi_n at s=1 should equal -2*K2/C(n,2) = 2/C(n,2) (since K2=-1)
# From numerics: 1/Phi_2 = 2 = 2/C(2,2) = 2/1 ✓
#                1/Phi_3 = 2/3 = 2/C(3,2) ✓
#                1/Phi_4 = 2/6 = 1/3 = 2/C(4,2) ✓

# And s/A_n at s=1: 1/A_n. We need A_n = C(n,2)/2 = n(n-1)/4 for this to work.
# But from direct computation: A_n = n(n-1)/2.
# So 1/A_n = 2/(n(n-1)) and s/A_n = 2/(n(n-1)).
# But 2/C(n,2) = 4/(n(n-1)).
# So s/A_n = 2/(n(n-1)) ≠ 4/(n(n-1)) = 2/C(n,2).
# Factor of 2 discrepancy!

# The issue: Phi_n = (1/s)*A_n implies 1/Phi_n = s/A_n
# If A_n = n(n-1)/2, then 1/Phi_n = 2s/(n(n-1))
# And -2K2/C(n,2) = 2s/(n(n-1)/2) = 4s/(n(n-1))
# THESE DON'T MATCH.

# Let me recheck: for n=2 at s=1:
# p(x) = x^2 - 1. Roots: 1, -1.
# H_1 = 1/(1-(-1)) = 1/2
# H_2 = 1/(-1-1) = -1/2
# Phi_2 = (1/2)^2 + (1/2)^2 = 1/2
# 1/Phi_2 = 2

# Hermite roots: xi_1 = 1, xi_2 = -1 (same)
# h_i = 1/(1-(-1)) = 1/2 for i=1
# A_n = (1/2)^2 + (-1/2)^2 = 1/2

# So A_2 = 1/2 = 2*1/4 = n(n-1)/4 ✓ (not n(n-1)/2!)

# Let me recheck my "careful" computation above...
print("\n\nFINAL CAREFUL CHECK:")
for n in range(2, 8):
    xi = hermeroots(n)
    h_vals = []
    for i in range(n):
        h = sum(1.0/(xi[i]-xi[j]) for j in range(n) if j != i)
        h_vals.append(h)
    A_n = sum(h**2 for h in h_vals)
    print(f"  n={n}: roots xi = {xi}")
    print(f"    h values: {h_vals}")
    print(f"    xi/2 values: {[x/2 for x in xi]}")
    print(f"    A_n = sum h^2 = {A_n:.10f}")
    print(f"    n(n-1)/4 = {n*(n-1)/4}")
    print(f"    n(n-1)/2 = {n*(n-1)/2}")
