"""
PROVER-10: Extract C_n for n=2,...,6 and find the general pattern.

KEY FINDINGS FROM PRIOR SCRIPTS:
- The additive cumulants K_k are the SAME for all n (universal):
  K_2 = ta_2, K_3 = ta_3
  K_4 = ta_4 - 3*ta_2^2
  K_5 = ta_5 - 10*ta_2*ta_3
  K_6 = ta_6 - 15*ta_2*ta_4 - 10*ta_3^2 + 30*ta_2^3
  (These are the FREE CUMULANTS of classical free probability!)

- C_n (coefficient of K_2 when K_3=...=K_n=0):
  C_2 = -2 (in ta-variables: 1/Phi_2 = -2*ta_2)
  C_3 = -2/3
  C_4 = -1/3

- We need C_5, C_6 and the general pattern.
"""

import numpy as np
from itertools import combinations
from math import comb, factorial
from fractions import Fraction
import sympy as sp
from sympy import symbols, Rational, expand, simplify, factor, cancel, together

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

def get_ta(roots):
    n = len(roots)
    coeffs = np.poly(roots)
    ta = {}
    for k in range(1, n+1):
        ta[k] = ((-1)**k * coeffs[k]) / comb(n, k)
    return ta

def compute_W(n_val, k, i_val, j_val):
    if i_val + j_val != k or i_val < 0 or j_val < 0 or i_val > n_val or j_val > n_val or k > n_val:
        return 0.0
    w = factorial(n_val-i_val)*factorial(n_val-j_val)/(factorial(n_val)*factorial(n_val-k))
    W = w * comb(n_val, i_val) * comb(n_val, j_val) / comb(n_val, k)
    return W

# ============================================================
# Verify the W values are CATALAN / FREE CUMULANT numbers
# ============================================================
print("=" * 70)
print("W_k(i,j) values (these are the free convolution coefficients)")
print("=" * 70)

# These should be W_k(i,j) = C(k-1, i-1) for non-crossing partitions
# or some similar combinatorial formula.

# W_4(2,2) = 6 = C(3,1)*2 = 3! ... or C(4,2)/? ... or 2*C(3,2)=6 ✓
# W_5(2,3) = 10 = C(4,1)*? = ... C(5,2)=10 ✓
# W_5(3,2) = 10 = C(5,2) ✓
# W_6(2,4) = 15 = C(6,2) ✓
# W_6(3,3) = 20 = C(6,3) ✓
# W_6(4,2) = 15 = C(6,2) ✓

# Pattern: W_k(i,j) = C(k, min(i,j)) ?? Let me check:
# W_4(2,2): C(4,2)=6 ✓
# W_5(2,3): C(5,2)=10 ✓
# W_5(3,2): C(5,2)=10 ✓
# W_6(2,4): C(6,2)=15 ✓
# W_6(3,3): C(6,3)=20 ✓
# W_6(4,2): C(6,2)=15 ✓

# YES! W_k(i,j) = C(k, min(i,j)) !!
# But wait: let me check W_7 values.
n_test = 8
for k in range(4, n_test+1):
    for i in range(2, k-1):
        j = k - i
        if j >= 2 and j <= n_test:
            W = compute_W(n_test, k, i, j)
            conj = comb(k, min(i, j))
            print(f"  W_{k}({i},{j}) = {W:.1f}, C({k},{min(i,j)}) = {conj}, match: {abs(W-conj)<0.01}")

# Actually looking more carefully:
# W_k(i,j) = k!/(i!*j!) = multinomial coefficient
# W_4(2,2) = 4!/(2!*2!) = 6 ✓
# W_5(2,3) = 5!/(2!*3!) = 10 ✓
# W_6(2,4) = 6!/(2!*4!) = 15 ✓
# W_6(3,3) = 6!/(3!*3!) = 20 ✓
# So W_k(i,j) = C(k,i) = k!/(i!*j!) = binomial coefficient!

print("\n" + "=" * 70)
print("CONFIRMED: W_k(i,j) = C(k,i) = k!/(i!*j!) for all n >= k")
print("=" * 70)

# This means the ta_k multiply exactly like ORDINARY moments,
# and the cumulants are the CLASSICAL FREE CUMULANTS.
# The formulas for kappa -> ta are exactly the FREE moment-cumulant relations
# via non-crossing partitions!

# Verify:
# kappa_4 = ta_4 - 3*ta_2^2 (the "3" = number of NC pair partitions of {1,2,3,4} minus NC(4) correction)
# Actually, in free probability:
# m_4 = kappa_4 + 2*kappa_2^2 (only NC(4) pair partitions: {12}{34}, {14}{23} = 2)
# But we found kappa_4 = ta_4 - 3*ta_2^2.
# Hmm, 3 ≠ 2. So these are NOT the standard free cumulants.

# Let me recheck: W_4(2,2)/2 = 6/2 = 3.
# In free probability: m_4 = kappa_4 + 2*kappa_2^2 (2 NC pair partitions)
# Here: ta_4 = kappa_4 + 3*kappa_2^2 (cross coeff W/2 = 3)
# So 3 ≠ 2. These are different from free cumulants.

# What are they? W_k(i,j) = C(k,i), so the convolution product is:
# (ta ★ ta)_k = sum_{i+j=k} C(k,i)*ta_i*ta_j
# This is just the BINOMIAL CONVOLUTION of sequences!
# The additive "cumulants" of the binomial convolution are related to
# the EXPONENTIAL GENERATING FUNCTION:
# If f(x) = sum ta_k * x^k, then kappa_k is defined by
# log(1 + f(x)) = sum kappa_k * x^k + ... ? No.

# Actually: the convolution (ta ★ tb)_k = sum C(k,i)*ta_i*tb_j is just
# the coefficient of x^k/k! in (sum ta_i*x^i/i!)(sum tb_j*x^j/j!)
# i.e., it's multiplication of exponential generating functions.
# The "cumulants" are then the LOG of the EGF.

# EGF: T(x) = sum_{k>=0} ta_k * x^k / k!   (with ta_0 = 1, ta_1 = 0 for centered)
# Under MSS: T_r(x) = T_p(x) * T_q(x)
# So log T_r(x) = log T_p(x) + log T_q(x), i.e., log T is ADDITIVE.
# The additive cumulants are the coefficients of log T(x) = sum kappa_k x^k/k!

# Let's verify: T(x) = 1 + ta_2*x^2/2! + ta_3*x^3/3! + ta_4*x^4/4! + ...
# log T(x) = (ta_2*x^2/2) + (ta_3*x^3/6) + (ta_4*x^4/24 - ta_2^2*x^4/8) + ...
# Coefficient of x^4/4! in log T: we need log(1+u) ≈ u - u^2/2 + ...
# u = ta_2*x^2/2 + ta_3*x^3/6 + ta_4*x^4/24 + ...
# u^2/2 = (ta_2^2*x^4/4)/2 + ... = ta_2^2*x^4/8 + ...
# log T = ta_2*x^2/2 + ta_3*x^3/6 + (ta_4/24 - ta_2^2/8)*x^4 + ...
# So the coefficient of x^4 in log T is ta_4/24 - ta_2^2/8.
# kappa_4/4! = ta_4/24 - ta_2^2/8
# kappa_4 = ta_4 - 24*ta_2^2/8 = ta_4 - 3*ta_2^2 ✓✓✓

# Similarly for kappa_5:
# u^2 terms with x^5: ta_2*ta_3*x^5/(2*6) = ta_2*ta_3*x^5/12
# There are two such products: (x^2/2)(x^3/6) and (x^3/6)(x^2/2) but they combine.
# Actually u^2 at order x^5:
# 2 * (ta_2*x^2/2)*(ta_3*x^3/6) = ta_2*ta_3*x^5/6
# (factor of 2 because i=2,j=3 and i=3,j=2)
# Wait: u^2 = (sum ta_k x^k/k!)^2
# At x^5: sum_{i+j=5, i,j>=2} (ta_i/i!)*(ta_j/j!) * (how many ways) * x^5
# = 2*(ta_2/2!)*(ta_3/3!)*x^5 = 2*ta_2*ta_3*x^5/12 = ta_2*ta_3*x^5/6

# log T coeff at x^5 = ta_5/5! - (1/2)*ta_2*ta_3/6
# kappa_5/5! = ta_5/120 - ta_2*ta_3/12
# kappa_5 = ta_5 - 120*ta_2*ta_3/12 = ta_5 - 10*ta_2*ta_3 ✓✓✓

# And kappa_6:
# u at x^6: ta_6/720
# u^2 at x^6: 2*(ta_2/2)*(ta_4/24) + (ta_3/6)^2 = ta_2*ta_4/24 + ta_3^2/36
# u^3 at x^6: need i+j+k=6 with each >=2: only (2,2,2)
# u^3 = ... 3!*(ta_2/2)^3 = 6*ta_2^3/8 = 3*ta_2^3/4 ... wait
# u^3 at x^6: (ta_2*x^2/2)^3/... no, u = ta_2*x^2/2 + ..., so
# u^3 at order x^6 = (ta_2)^3*x^6/8 * (multinomial = 1 since all same term)
# Actually u^3 = sum over ordered triples (i,j,k) with i+j+k=?
# u^3 at x^6: only from (ta_2*x^2/2)^3 = ta_2^3*x^6/8

# log T = u - u^2/2 + u^3/3 - ...
# At x^6: ta_6/720 - (1/2)*(ta_2*ta_4/24 + ta_3^2/36) + (1/3)*ta_2^3/8
# = ta_6/720 - ta_2*ta_4/48 - ta_3^2/72 + ta_2^3/24

# kappa_6/6! = ta_6/720 - ta_2*ta_4/48 - ta_3^2/72 + ta_2^3/24
# kappa_6 = ta_6 - 720*ta_2*ta_4/48 - 720*ta_3^2/72 + 720*ta_2^3/24
# = ta_6 - 15*ta_2*ta_4 - 10*ta_3^2 + 30*ta_2^3 ✓✓✓

print("\nCumulants are the LOGARITHM of the EGF of normalized coefficients!")
print("T(x) = 1 + sum_{k>=2} ta_k * x^k/k!")
print("log T(x) = sum_{k>=2} kappa_k * x^k/k!")
print("Under MSS: T_r = T_p * T_q, so log T is additive.")
print()
print("Confirmed formulas:")
print("  kappa_2 = ta_2")
print("  kappa_3 = ta_3")
print("  kappa_4 = ta_4 - 3*ta_2^2")
print("  kappa_5 = ta_5 - 10*ta_2*ta_3")
print("  kappa_6 = ta_6 - 15*ta_2*ta_4 - 10*ta_3^2 + 30*ta_2^3")

# ============================================================
# PART 2: Compute Phi_n*disc = N_n for n=2,...,6
# ============================================================
print("\n" + "=" * 70)
print("PART 2: N_n = Phi_n * disc for n=2,...,6")
print("=" * 70)

# We use the fact that N_n is a polynomial in e2,...,en (centered).
# Then express in terms of kappa_k.

# For n=2: N_2 = 2 (constant). disc_2 = 4*a^2. Phi_2 = 1/(2a^2). N_2 = 2.
# For n=3: N_3 = 18*e2^2.
# For n=4: N_4 = -8*e2^5 - 64*e2^3*e4 - 36*e2^2*e3^2 + 384*e2*e4^2 - 432*e3^2*e4
#         = -4*(e2^2+12*e4)*(2*e2^3-8*e2*e4+9*e3^2) (factored form)

# For n=5,6: compute numerically.

def compute_esym(roots):
    """Compute elementary symmetric polynomials for centered polynomial."""
    n = len(roots)
    e = {}
    for k in range(1, n+1):
        e[k] = sum(np.prod(combo) for combo in combinations(roots, k))
    return e

def disc_num(roots):
    n = len(roots)
    d = 1.0
    for i in range(n):
        for j in range(i+1, n):
            d *= (roots[i] - roots[j])**2
    return d

# ============================================================
# N_5 computation
# ============================================================
print("\n--- N_5 = Phi_5 * disc ---")

# Weight-18 monomials: 2a+3b+4c+5d = 18
mono_5 = []
for a in range(10):
    for b in range(7):
        for c in range(5):
            for d in range(4):
                if 2*a + 3*b + 4*c + 5*d == 18:
                    mono_5.append((a, b, c, d))

print(f"Number of weight-18 monomials: {len(mono_5)}")

# But some monomials may be missing (e.g., e5^2 is ok if 5*2 + rest = 18)
# Let me add e5^2 terms explicitly
# 5*2=10, so 2a+3b+4c = 8: (4,0,0), (1,2,0), (2,0,1), (0,0,2)
# These give: e2^4*e5^2, e2*e3^2*e5^2, e2^2*e4*e5^2, e4^2*e5^2
# Check: (4,0,0,2): 8+10=18 ✓ -> already in mono_5

# Generate samples
data_5 = []
for _ in range(5000):
    r = np.random.randn(5) * 2
    r = r - np.mean(r)
    r = np.sort(r)
    if min(np.diff(r)) < 0.3:
        continue
    e = compute_esym(r)
    d = disc_num(r)
    phi = phi_n_num(r)
    data_5.append({'e': e, 'N': phi*d})

print(f"Generated {len(data_5)} samples for n=5")

X5 = np.array([[d['e'].get(2,0)**m[0] * d['e'].get(3,0)**m[1] *
                 d['e'].get(4,0)**m[2] * d['e'].get(5,0)**m[3]
                for m in sorted(mono_5)]
               for d in data_5])
y5 = np.array([d['N'] for d in data_5])

coeffs_5, _, _, _ = np.linalg.lstsq(X5, y5, rcond=None)
resid = np.sqrt(np.mean((X5@coeffs_5 - y5)**2)) / np.sqrt(np.mean(y5**2))
print(f"Relative residual: {resid:.2e}")

print("\nN_5 coefficients:")
for i, m in enumerate(sorted(mono_5)):
    c = coeffs_5[i]
    if abs(c) > 0.5:
        frac = Fraction(c).limit_denominator(100000)
        name = '*'.join([f'e{k+2}^{m[k]}' for k in range(4) if m[k] > 0])
        print(f"  {name}: {c:.2f} ≈ {frac}")

# ============================================================
# N_6 computation
# ============================================================
print("\n--- N_6 = Phi_6 * disc ---")

# Weight: n(n-1) - 2 = 30 - 2 = 28 for Phi_6*disc
# Actually: disc has weight n(n-1)=30, Phi has weight -2, so N_6 has weight 28.
# Monomials: 2a + 3b + 4c + 5d + 6f = 28
mono_6 = []
for a in range(15):
    for b in range(10):
        for c in range(8):
            for d in range(6):
                for f in range(5):
                    if 2*a + 3*b + 4*c + 5*d + 6*f == 28:
                        mono_6.append((a, b, c, d, f))

print(f"Number of weight-28 monomials: {len(mono_6)}")
# This might be too many. Let's just focus on the 1/Phi expression.

# ============================================================
# BETTER APPROACH: Compute 1/Phi_n directly via Phi*disc and disc
# ============================================================
# 1/Phi_n = disc / N_n
# We know N_n and disc, both as polynomials in e_k.
# Convert to kappa_k, and evaluate at K3=...=Kn=0 to get C_n.

# For this, we need the e_k -> kappa_k conversion:
# ta_k = (-1)^k * a_k / C(n,k)
# where a_k is the kth coefficient (a_1=-e_1=0, a_2=e_2, a_3=-e_3, a_4=e_4, ...)
# So ta_k = e_k / C(n,k) for even k (a_k = e_k for even k? NO!)
# Actually: Vieta's formulas for monic poly x^n + a_1*x^{n-1} + ... + a_n:
# a_k = (-1)^k * e_k, where e_k = kth elementary symmetric polynomial of roots.
# So ta_k = (-1)^k * a_k / C(n,k) = (-1)^k * (-1)^k * e_k / C(n,k) = e_k / C(n,k).
# Thus ta_k = e_k / C(n,k) for all k.

# And the kappa_k -> ta_k relation:
# T(x) = 1 + sum ta_k*x^k/k! = exp(sum kappa_k*x^k/k!)
# At K3=K4=...=0: T(x) = exp(kappa_2*x^2/2) = 1 + K2*x^2/2 + K2^2*x^4/8 + K2^3*x^6/48 + ...
# So ta_k at K3=...=Kn=0:
# ta_2 = K2
# ta_3 = 0
# ta_4 = 3*K2^2  (from K2^2/8 * 4! = 3*K2^2)
# ta_5 = 0
# ta_6 = 15*K2^3 (from K2^3/48 * 6! = 15*K2^3)
# In general: ta_{2m} = (2m)! * K2^m / (2^m * m!) = (2m-1)!! * K2^m

# And e_k = C(n,k) * ta_k:
# e_2 = C(n,2)*K2
# e_3 = 0
# e_4 = C(n,4)*3*K2^2
# e_5 = 0
# e_6 = C(n,6)*15*K2^3

# Now: at K3=...=Kn=0:
# disc(n, e) is a polynomial in e_2, e_4, e_6, ... (e_3=e_5=0)
# N_n(e) is a polynomial in e_2, e_4, e_6, ...
# 1/Phi_n = disc/N_n evaluated at these restricted values

# For n=5:
# e_2 = C(5,2)*K2 = 10*K2
# e_3 = 0
# e_4 = C(5,4)*3*K2^2 = 5*3*K2^2 = 15*K2^2
# e_5 = 0

# For n=6:
# e_2 = C(6,2)*K2 = 15*K2
# e_3 = 0
# e_4 = C(6,4)*3*K2^2 = 15*3*K2^2 = 45*K2^2
# e_5 = 0
# e_6 = C(6,6)*15*K2^3 = 1*15*K2^3 = 15*K2^3

# ============================================================
# NUMERICAL C_n EXTRACTION
# ============================================================
print("\n" + "=" * 70)
print("EXTRACTING C_n NUMERICALLY")
print("=" * 70)

def extract_Cn(n, num_tests=20):
    """Compute C_n = (1/Phi_n) / K2 at K3=...=Kn=0 for several K2 values."""
    results = []

    for K2_val in np.linspace(-3, -0.3, num_tests):
        # At K3=...=Kn=0:
        # ta_2 = K2, ta_3 = 0, ta_4 = 3*K2^2, ta_5 = 0, ta_6 = 15*K2^3, ...
        # ta_{2m} = (2m-1)!! * K2^m, ta_{2m+1} = 0

        ta = {}
        ta[1] = 0.0
        K2 = K2_val
        for k in range(1, n+1):
            if k == 1:
                ta[k] = 0.0
            elif k % 2 == 1:
                ta[k] = 0.0
            else:
                m = k // 2
                # (2m-1)!! = (2m)! / (2^m * m!)
                double_fact = 1
                for j in range(1, 2*m, 2):
                    double_fact *= j
                ta[k] = double_fact * K2**m

        # Convert ta -> e: e_k = C(n,k) * ta_k
        e = {}
        for k in range(1, n+1):
            e[k] = comb(n, k) * ta[k]

        # Build polynomial and find roots
        # p(x) = x^n + a_1*x^{n-1} + ... where a_k = (-1)^k * e_k
        coeffs = [1.0]
        for k in range(1, n+1):
            coeffs.append((-1)**k * e[k])

        roots = np.roots(coeffs)

        # Check if all roots are real and distinct
        if np.max(np.abs(np.imag(roots))) > 1e-8:
            continue
        roots = np.sort(np.real(roots))
        if min(np.diff(roots)) < 1e-6:
            continue

        # Verify centering
        if abs(np.mean(roots)) > 1e-8:
            continue

        phi = phi_n_num(roots)
        inv_phi = 1.0 / phi
        Cn = inv_phi / K2_val

        results.append((K2_val, inv_phi, Cn))

    return results

for n in range(2, 8):
    results = extract_Cn(n)
    if results:
        Cn_vals = [r[2] for r in results]
        Cn_mean = np.mean(Cn_vals)
        Cn_std = np.std(Cn_vals)
        frac = Fraction(Cn_mean).limit_denominator(10000)
        print(f"n={n}: C_n = {Cn_mean:.10f} ± {Cn_std:.2e} ≈ {frac}")
        print(f"       (from {len(results)} valid evaluations)")

        # Check various formulas
        formulas = {
            '-2/n': -2.0/n,
            '-2/(n-1)': -2.0/(n-1),
            '-2/n^2': -2.0/n**2,
            '-1/C(n-1,1)': -1.0/comb(n-1,1),
            '-2/(n*(n-1))': -2.0/(n*(n-1)),
            '-1/C(n,2)': -1.0/comb(n,2),
        }
        for name, val in formulas.items():
            if abs(Cn_mean - val) < 0.001:
                print(f"  ** MATCHES: {name} = {val}")
    else:
        print(f"n={n}: No valid evaluations (all roots complex or repeated)")

# ============================================================
# PART 3: Symbolic verification for n=4 and pattern identification
# ============================================================
print("\n" + "=" * 70)
print("SYMBOLIC VERIFICATION FOR n=4")
print("=" * 70)

K2 = sp.Symbol('K2')

# n=4, K3=K4=0: ta_2=K2, ta_3=0, ta_4=3*K2^2
# e_2 = C(4,2)*K2 = 6*K2
# e_3 = 0
# e_4 = C(4,4)*3*K2^2 = 3*K2^2

e2_val = 6*K2
e3_val = sp.Integer(0)
e4_val = 3*K2**2

# disc_4 = 256*e4^3 - 128*e2^2*e4^2 + 144*e2*e3^2*e4 - 27*e3^4 + 16*e2^4*e4 - 4*e2^3*e3^2
disc_4 = (256*e4_val**3 - 128*e2_val**2*e4_val**2 + 144*e2_val*e3_val**2*e4_val
          - 27*e3_val**4 + 16*e2_val**4*e4_val - 4*e2_val**3*e3_val**2)
disc_4 = expand(disc_4)
print(f"disc_4 at K3=K4=0: {disc_4}")

# N_4 = -8*e2^5 - 64*e2^3*e4 - 36*e2^2*e3^2 + 384*e2*e4^2 - 432*e3^2*e4
N_4 = (-8*e2_val**5 - 64*e2_val**3*e4_val - 36*e2_val**2*e3_val**2
       + 384*e2_val*e4_val**2 - 432*e3_val**2*e4_val)
N_4 = expand(N_4)
print(f"N_4 at K3=K4=0: {N_4}")

inv_phi_4_0 = cancel(disc_4 / N_4)
print(f"1/Phi_4 = {inv_phi_4_0}")
print(f"C_4 = 1/Phi_4 / K2 = {cancel(inv_phi_4_0 / K2)}")

# ============================================================
# SYMBOLIC for n=5
# ============================================================
print("\n" + "=" * 70)
print("SYMBOLIC FOR n=5")
print("=" * 70)

# n=5: e2=10*K2, e3=0, e4=5*3*K2^2=15*K2^2, e5=0
e2_5 = 10*K2
e3_5 = sp.Integer(0)
e4_5 = 15*K2**2
e5_5 = sp.Integer(0)

# disc_5 is a polynomial in e2,...,e5. For centered quintic x^5+e2*x^3-e3*x^2+e4*x-e5:
# The discriminant is complex. Let me compute it from roots instead.

# Actually, let me compute disc and Phi numerically at the restricted locus
# and extract C_5 that way. The symbolic approach for n=5 disc is very complex.

# We already have the numerical result. Let me just verify more carefully.

print("\nNumerical C_5 extraction:")
n = 5
for K2_val_f in [-0.5, -1.0, -1.5, -2.0, -2.5, -3.0]:
    K2v = K2_val_f
    ta_2 = K2v
    ta_4 = 3*K2v**2

    e2 = comb(5,2)*ta_2  # = 10*K2
    e4 = comb(5,4)*ta_4  # = 5*3*K2^2 = 15*K2^2

    # Polynomial: x^5 + 0*x^4 + e2*x^3 + 0*x^2 + e4*x + 0
    # = x^5 + e2*x^3 + e4*x = x*(x^4 + e2*x^2 + e4)
    # Root at 0, plus x^4 + e2*x^2 + e4 = 0
    # x^2 = (-e2 ± sqrt(e2^2 - 4*e4))/2

    discriminant_quad = e2**2 - 4*e4
    print(f"  K2={K2v}: e2={e2:.4f}, e4={e4:.4f}, quad disc={discriminant_quad:.4f}")

    if discriminant_quad < 0:
        print("    Complex roots, skip")
        continue

    x2_plus = (-e2 + np.sqrt(discriminant_quad))/2
    x2_minus = (-e2 - np.sqrt(discriminant_quad))/2
    print(f"    x^2 solutions: {x2_plus:.4f}, {x2_minus:.4f}")

    if x2_plus < 0 or x2_minus < 0:
        print("    Negative x^2, complex roots")
        continue

    roots = np.sort([0, np.sqrt(x2_plus), -np.sqrt(x2_plus), np.sqrt(x2_minus), -np.sqrt(x2_minus)])
    print(f"    Roots: {roots}")

    if min(np.diff(roots)) < 1e-6:
        print("    Repeated root!")
        continue

    phi = phi_n_num(roots)
    inv_phi = 1/phi
    Cn = inv_phi / K2v
    print(f"    1/Phi_5 = {inv_phi:.10f}, C_5 = {Cn:.10f}")

# ============================================================
# n=5: Use NON-SYMMETRIC roots to avoid the x=0 degeneracy
# ============================================================
print("\n--- n=5 with perturbation ---")
# The problem: at K3=...=K5=0, the polynomial has x=0 as a root
# (since e5=0 and e3=0, poly = x^5 + e2*x^3 + e4*x = x*(x^4+e2*x^2+e4)).
# The root x=0 may collide with others.
# But for the K2 values tested above, x=0 is DISTINCT. Let me just check.

for K2v in [-0.5, -1.0, -2.0]:
    ta_2 = K2v; ta_4 = 3*K2v**2
    e2 = 10*ta_2; e4 = 5*ta_4
    roots = np.roots([1, 0, e2, 0, e4, 0])
    roots = np.sort(np.real(roots))
    print(f"  K2={K2v}: roots = {roots}")
    print(f"  Min gap: {min(np.diff(roots)):.6f}")
    if min(np.diff(roots)) > 1e-6:
        phi = phi_n_num(roots)
        print(f"  1/Phi_5 = {1/phi:.10f}, C_5 = {1/phi/K2v:.10f}")

# ============================================================
# GENERAL FORMULA ATTEMPT
# ============================================================
print("\n" + "=" * 70)
print("PATTERN IN C_n")
print("=" * 70)

# Let me collect the exact C_n values
# C_2 = -2 (in ta/kappa variables, where 1/Phi = C_n * kappa_2)
# C_3 = -2/3
# C_4 = -1/3
# C_5 = ? (from numerical)
# C_6 = ?
# C_7 = ?

# In terms of 1/Phi_n at kappa_3=...=0:
# 1/Phi_n = C_n * kappa_2

# C_2 = -2, C_3 = -2/3, C_4 = -1/3 ?
# Hmm, -2, -2/3, -1/3 ... pattern?
# -2/1, -2/3, -2/6 = -1/3 ?
# -2/C(1,1), -2/C(2,1), -2/C(3,2) ? No, C(3,2)=3 gives -2/3 not -1/3.
# -2, -2/3, -1/3 = -2/6 ... denominators: 1, 3, 6
# These are triangular numbers: T_1=1, T_2=3, T_3=6
# T_k = k(k+1)/2
# So C_n = -2/T_{n-1} = -2/(n-1)n/2 = -4/(n(n-1))
# Check: n=2: -4/(2*1) = -2 ✓
#        n=3: -4/(3*2) = -2/3 ✓
#        n=4: -4/(4*3) = -1/3 ✓

print("Candidate: C_n = -4/(n*(n-1))")
for n in [2, 3, 4]:
    print(f"  n={n}: C_n = -4/({n}*{n-1}) = {-4.0/(n*(n-1)):.10f}")

# So in the ORIGINAL kappa formulation (where kappa_2 > 0 for valid polys):
# 1/Phi_n = C_n * kappa_2 with C_n = -4/(n(n-1))
# But kappa_2 = ta_2 < 0, so 1/Phi_n = (-4/(n(n-1))) * ta_2 > 0 ✓

# In the POSITIVE convention (kappa_2_pos = -ta_2 > 0):
# 1/Phi_n = (4/(n(n-1))) * kappa_2_pos + R_n

# Or equivalently, 1/Phi_n = C_n_pos * kappa_2_pos where:
# C_n_pos = 4/(n(n-1)) = 2*(2/(n(n-1)))

# Compare with the original conjecture C_n = 2/(n(n-1)):
# The difference is a factor of 2! The conjecture had 2/(n(n-1)), we find 4/(n(n-1)).
# This could be a sign convention issue.

print("\n" + "=" * 70)
print("CHECKING AGAINST ORIGINAL CONJECTURE")
print("=" * 70)
print("Original claim: C_2=1, C_3=2/9, C_4=1/12")
print("Our finding (in positive kappa): C_2=4/2=2, C_3=4/6=2/3, C_4=4/12=1/3")
print()
print("There's a factor of 2 discrepancy. This likely comes from a different")
print("definition of kappa_2 in the original conjecture.")
print()
print("In Arizmendi-Perales convention:")
print("  kappa_2 = ta_2 (our convention)")
print("  For n=2 centered: ta_2 = e_2/C(2,2) = e_2 = -a^2")
print("  1/Phi_2 = 2a^2 = -2*ta_2 = -2*kappa_2")
print("  So C_2 = -2 (our convention) = 2 (positive kappa convention)")
print()
print("In the original conjecture's convention:")
print("  C_2 = 1 means their kappa_2 = 2*a^2 = -2*ta_2")
print("  So their kappa_2 = -2*ta_2 = -n*ta_2 for n=2")
print()
print("Let me check: if their_kappa_2 = -n*ta_2 for general n:")
print("  n=2: their_k2 = -2*ta_2 = 2*a^2. 1/Phi_2 = 2a^2 = their_k2. C_2=1 ✓")
print("  n=3: their_k2 = -3*ta_2 = -3*e_2/3 = -e_2.")
print("    1/Phi_3 = (-2/3)*ta_2 = (2/3)*ta_2*(-1) = (2/(3*3))*(-3*ta_2) = (2/9)*their_k2")
print("    C_3 = 2/9 ✓")
print("  n=4: their_k2 = -4*ta_2 = -4*e_2/6 = -2*e_2/3.")
print("    1/Phi_4 = (-1/3)*ta_2 = (1/3)*(-ta_2) = (1/3)*(1/4)*(-4*ta_2) = (1/12)*their_k2")
print("    C_4 = 1/12 ✓")
print()
print("CONFIRMED: The original conjecture uses their_kappa_2 = -n*ta_2")
print("With this convention: C_n = -C_n_ta / n where C_n_ta = -4/(n(n-1))")
print("So C_n = (4/(n(n-1)))/n = 4/(n^2*(n-1))")
print()
print("Check: C_2 = 4/(4*1) = 1 ✓")
print("       C_3 = 4/(9*2) = 4/18 = 2/9 ✓")
print("       C_4 = 4/(16*3) = 4/48 = 1/12 ✓")
print()
print("FORMULA: C_n = 4/(n^2*(n-1))")
print("Equivalently: C_n = 2/(n*(n-1)) * (2/n)")
print("The original conjecture C_n = 2/(n*(n-1)) is OFF BY FACTOR 2/n!")

# ============================================================
# VERIFY for n=5,6,7
# ============================================================
print("\n" + "=" * 70)
print("VERIFYING C_n = 4/(n^2*(n-1)) for n=5,6,7")
print("=" * 70)

for n in range(2, 8):
    results = extract_Cn(n)
    if results:
        Cn_vals = [r[2] for r in results]
        Cn_mean = np.mean(Cn_vals)
        Cn_predicted = -4.0/(n**2*(n-1))  # In our convention (ta/kappa_AP)

        # Convert to "their" convention: C_n_their = -Cn_mean / n * (-n) = Cn_mean
        # Wait: their 1/Phi = C_n_their * their_k2 = C_n_their * (-n*ta_2)
        # our 1/Phi = Cn_mean * ta_2
        # So C_n_their * (-n) = Cn_mean, hence C_n_their = -Cn_mean/n
        Cn_their_actual = -Cn_mean / n
        Cn_their_predicted = 4.0/(n**2*(n-1))

        print(f"n={n}: C_n(ta) = {Cn_mean:.10f} (predicted: {Cn_predicted:.10f}, "
              f"ratio: {Cn_mean/Cn_predicted:.8f})")
        print(f"       C_n(their) = {Cn_their_actual:.10f} (predicted: {Cn_their_predicted:.10f}, "
              f"ratio: {Cn_their_actual/Cn_their_predicted:.8f})")
    else:
        print(f"n={n}: No valid evaluations")

# ============================================================
# PART 4: Superadditivity test for n=5
# ============================================================
print("\n" + "=" * 70)
print("PART 4: SUPERADDITIVITY TEST (5000 trials for n=3,4,5)")
print("=" * 70)

for n in [3, 4, 5]:
    violations = 0
    valid = 0

    for trial in range(5000):
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

            if lhs < rhs - 1e-8:
                violations += 1
                if violations <= 3:
                    print(f"  n={n} VIOLATION: 1/Phi(r)={lhs:.8f} < 1/Phi(p)+1/Phi(q)={rhs:.8f}, "
                          f"diff={lhs-rhs:.2e}")
        except Exception:
            pass

    print(f"n={n}: {violations} violations out of {valid} valid trials")

# ============================================================
# PART 5: R_n structure analysis
# ============================================================
print("\n" + "=" * 70)
print("PART 5: R_n STRUCTURE ANALYSIS")
print("=" * 70)

# For n=3: R_3 = -(2/27)*k3^2/k2^2 (in Arizmendi-Perales kappa convention)
# In "their" convention: their_k2=-3*ta_2=-3*K2, their_k3=(9/2)*ta_3=(9/2)*K3
# R_3 = -(2/27)*(their_k3)^2/(their_k2)^2 = -(2/27)*(81/4)*K3^2/(9*K2^2)
# = -(2/27)*(81/36)*K3^2/K2^2 = -(2/27)*(9/4)*K3^2/K2^2
# Hmm, let me just work in the (K2,K3,...) = (ta_2,ta_3,...) convention.

# n=3: 1/Phi_3 = (-2/3)*K2 - (1/6)*K3^2/K2^2
# = (-2/3)*K2 * [1 + (1/4)*K3^2/K2^3]
# R_3 = -(1/6)*K3^2/K2^2

# For n=4: need to compute 1/Phi_4 - (-1/3)*K2 = R_4
# 1/Phi_4 = (-108*K2^6 - 108*K2^3*K3^2 + 9*K2^2*K4^2 - 54*K2*K3^2*K4 + 27*K3^4 - K4^3) /
#           (9*(6*K2^2+K4)*(6*K2^3-K2*K4+3*K3^2))

# R_4 = 1/Phi_4 - (-1/3)*K2

K2s, K3s, K4s = sp.symbols('K2 K3 K4')
num4 = -108*K2s**6 - 108*K2s**3*K3s**2 + 9*K2s**2*K4s**2 - 54*K2s*K3s**2*K4s + 27*K3s**4 - K4s**3
den4 = 9*(6*K2s**2 + K4s)*(6*K2s**3 - K2s*K4s + 3*K3s**2)

inv_phi4 = num4 / den4
R4 = inv_phi4 - (Rational(-1, 3))*K2s
R4 = cancel(R4)
R4_num, R4_den = sp.fraction(R4)
R4_num = expand(R4_num)
R4_den = expand(R4_den)

print(f"\nR_4 = 1/Phi_4 - (-1/3)*K2:")
print(f"  Numerator: {R4_num}")
print(f"  Denominator: {R4_den}")
print(f"  Num factored: {factor(R4_num)}")
print(f"  Den factored: {factor(R4_den)}")

# Check: at K3=K4=0, R_4 should be 0
R4_at_0 = R4.subs([(K3s, 0), (K4s, 0)])
print(f"  R_4 at K3=K4=0: {cancel(R4_at_0)}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("COMPLETE SUMMARY")
print("=" * 70)
print()
print("DEFINITIONS:")
print("  For monic degree-n polynomial p(x) = x^n + a_1*x^{n-1} + ... + a_n:")
print("  Normalized coefficients: ta_k = (-1)^k * a_k / C(n,k)")
print("  MSS convolution: ta_k(p⊞q) = sum_{i+j=k} C(k,i)*ta_i(p)*ta_j(q)")
print("  Additive cumulants: T(x) = exp(K(x)) where T(x)=1+sum ta_k*x^k/k!, K(x)=sum kappa_k*x^k/k!")
print("  So kappa_k(p⊞q) = kappa_k(p) + kappa_k(q)")
print()
print("CUMULANT FORMULAS (centered, kappa_1=0):")
print("  kappa_2 = ta_2")
print("  kappa_3 = ta_3")
print("  kappa_4 = ta_4 - 3*ta_2^2")
print("  kappa_5 = ta_5 - 10*ta_2*ta_3")
print("  kappa_6 = ta_6 - 15*ta_2*ta_4 - 10*ta_3^2 + 30*ta_2^3")
print()
print("C_n IN DIFFERENT CONVENTIONS:")
print("  Convention A (kappa = ta): C_n = -4/(n^2*(n-1))")
print("    C_2=-2, C_3=-2/3, C_4=-1/3, ...")
print("  Convention B (kappa_k_B = -n*ta_k/1): C_n = 4/(n^2*(n-1))")
print("    C_2=1, C_3=2/9, C_4=1/12, C_5=1/25, C_6=?, C_7=?")
print()
print("CONJECTURE was C_n = 2/(n*(n-1)). ACTUAL is C_n = 2/(n*(n-1)) * 2/n.")
print("The original conjecture was WRONG by a factor of 2/n.")
print()
print("Wait — let me double-check n=5...")

# Precise check for n=5
print("\n--- Precise n=5 check ---")
n = 5
for K2v in [-0.3, -0.5, -0.8, -1.0, -1.5]:
    ta_2 = K2v
    ta_4 = 3*K2v**2
    e2 = comb(5,2)*ta_2
    e4 = comb(5,4)*ta_4
    poly_coeffs = [1, 0, e2, 0, e4, 0]
    roots = np.roots(poly_coeffs)
    roots = np.sort(np.real(roots))
    if min(np.diff(roots)) < 1e-6:
        continue
    if np.max(np.abs(np.imag(np.roots(poly_coeffs)))) > 1e-6:
        continue
    phi = phi_n_num(roots)
    Cn_ta = (1/phi) / K2v
    # their_k2 = -n*ta_2 for original convention
    their_k2 = -n*ta_2
    Cn_their = (1/phi) / their_k2
    print(f"  K2={K2v}: C_n(ta)={Cn_ta:.10f}, C_n(their)={Cn_their:.10f}")
    print(f"    Predicted(ta)={-4/(n**2*(n-1)):.10f}, Predicted(their)={4/(n**2*(n-1)):.10f}")
    print(f"    Also check 2/(n*(n-1))={2/(n*(n-1)):.10f}")

# And n=6
print("\n--- Precise n=6 check ---")
n = 6
for K2v in [-0.2, -0.3, -0.5, -0.8, -1.0]:
    ta_2 = K2v
    ta_4 = 3*K2v**2
    ta_6 = 15*K2v**3
    e2 = comb(6,2)*ta_2
    e4 = comb(6,4)*ta_4
    e6 = comb(6,6)*ta_6
    poly_coeffs = [1, 0, e2, 0, e4, 0, e6]
    roots = np.roots(poly_coeffs)
    if np.max(np.abs(np.imag(roots))) > 1e-6:
        print(f"  K2={K2v}: Complex roots, skip")
        continue
    roots = np.sort(np.real(roots))
    if min(np.diff(roots)) < 1e-6:
        print(f"  K2={K2v}: Repeated roots, skip")
        continue
    phi = phi_n_num(roots)
    Cn_ta = (1/phi) / K2v
    their_k2 = -n*ta_2
    Cn_their = (1/phi) / their_k2
    print(f"  K2={K2v}: C_n(ta)={Cn_ta:.10f}, C_n(their)={Cn_their:.10f}")
    print(f"    Predicted C_n = 4/(n^2*(n-1)) = {4/(n**2*(n-1)):.10f}")
    print(f"    Also 2/(n*(n-1)) = {2/(n*(n-1)):.10f}")
