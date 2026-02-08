"""
prove_n4_coefficient.py -- Attempt to prove Fisher superadditivity for n=4
using the coefficient-level approach that succeeded at n=3.

For n=3, the proof worked because:
  1. Boxplus is additive in (e2, e3) for centered cubics.
  2. 1/Phi_3 = (4E^3 - 27F^2) / (18E^2) — a simple rational function.
  3. The excess is a positive definite quadratic form in (F_p, F_q).

For n=4, new complications:
  - g_4 has a cross term: e_4(r) = e_4(p) + e_4(q) + (1/6)*e_2(p)*e_2(q)
  - 1/Phi_4 is a more complex rational function of (e2, e3, e4)

Author: Prover agent
Date: 2026-02-08
"""
import sympy as sp
from sympy import (symbols, Rational, simplify, factor, expand, cancel,
                   Poly, together, numer, denom, binomial, S, sqrt, collect,
                   Matrix, resultant, discriminant as sp_disc, prod as sp_prod)
import numpy as np
from math import comb
from itertools import combinations
import sys
import time

# =====================================================================
# PART 1: Compute 1/Phi_4 symbolically for a centered quartic
# =====================================================================
print("=" * 72)
print("PART 1: Derive 1/Phi_4 in terms of roots (then coefficients)")
print("=" * 72)

# Approach: work with roots r1, r2, r3, r4 with r1+r2+r3+r4=0
# f(x) = (x-r1)(x-r2)(x-r3)(x-r4) = x^4 + e2*x^2 + e3*x + e4
# (e1 = 0 because centered)

r1, r2, r3 = symbols('r1 r2 r3', real=True)
r4_expr = -r1 - r2 - r3  # centered

# Elementary symmetric polynomials
e2_sym = expand(r1*r2 + r1*r3 + r1*r4_expr + r2*r3 + r2*r4_expr + r3*r4_expr)
e3_sym = expand(r1*r2*r3 + r1*r2*r4_expr + r1*r3*r4_expr + r2*r3*r4_expr)
e4_sym = expand(r1*r2*r3*r4_expr)

print(f"e2 = {e2_sym}")
print(f"e3 = {e3_sym}")
print(f"e4 = {e4_sym}")

# f(x) = x^4 + e2*x^2 + e3*x + e4
# f'(x) = 4*x^3 + 2*e2*x + e3
# For root r_i: f'(r_i) = prod_{j!=i} (r_i - r_j)
# H(r_i) = sum_{j!=i} 1/(r_i - r_j) = f''(r_i)/(2*f'(r_i)) ... NO
# Actually: H(r_i) = f'(r_i) / f''(r_i) ... NO
# H(r_i) = sum_{j!=i} 1/(r_i - r_j) = d/dx [ln f'(x)] at x=r_i ... NO

# Let me be precise:
# f(x) = prod(x - r_j)
# f'(r_i) = prod_{j!=i}(r_i - r_j)
# H_f(r_i) = sum_{j!=i} 1/(r_i - r_j)
# Note: H_f(r_i) = [d/dx ln|f(x)/(x-r_i)|]_{x=r_i} = [f'(x)/(x-r_i)]_{x=r_i}/f'(r_i) ...
# Actually H(r_i) = [sum_{j=1}^n 1/(x-r_j)]_{x=r_i, j!=i} = [f'(x)/f(x)]_{regularized at r_i}

# Simpler: just compute directly with roots
roots = [r1, r2, r3, r4_expr]
n = 4

print("\nComputing Phi_4 symbolically (this may take a moment)...")
sys.stdout.flush()

# Compute H(r_i) and Phi_4
H_list = []
for i in range(n):
    H_i = sum(1/(roots[i] - roots[j]) for j in range(n) if j != i)
    H_list.append(H_i)

# Phi_4 = sum H_i^2
print("Computing Phi_4 = sum H_i^2...")
sys.stdout.flush()

# This can be very expensive symbolically. Let's try a different approach:
# compute numerator and denominator of each H_i carefully.

# f'(r_i) = prod_{j!=i}(r_i - r_j)
fprime_list = []
for i in range(n):
    fp = sp.Integer(1)
    for j in range(n):
        if j != i:
            fp *= (roots[i] - roots[j])
    fprime_list.append(expand(fp))

# H(r_i) = sum_{j!=i} 1/(r_i - r_j)
# = sum_{j!=i} prod_{k!=i,k!=j}(r_i - r_k) / prod_{k!=i}(r_i - r_k)
# = [sum_{j!=i} prod_{k!=i,k!=j}(r_i - r_k)] / f'(r_i)

# The numerator of H(r_i) is f''(r_i)/2 (this is a standard identity)
# f''(x) = sum_{j} sum_{k!=j} prod_{l!=j,l!=k}(x - r_l)
# f''(r_i) = 2 * sum_{j!=i} prod_{k!=i,k!=j}(r_i - r_k)
# So H(r_i) = f''(r_i) / (2 * f'(r_i))

x = symbols('x')
f_poly = expand((x - r1)*(x - r2)*(x - r3)*(x - r4_expr))
f_prime = sp.diff(f_poly, x)
f_double_prime = sp.diff(f_poly, x, 2)

print("Using H(r_i) = f''(r_i) / (2*f'(r_i))...")

# Compute numerically first to verify formula, then symbolically
H_num_list = []
for i in range(n):
    fpp_ri = expand(f_double_prime.subs(x, roots[i]))
    fp_ri = fprime_list[i]
    H_num_list.append((fpp_ri, fp_ri))

# Phi_4 = sum [f''(r_i)/(2*f'(r_i))]^2 = (1/4) * sum [f''(r_i)/f'(r_i)]^2
# 1/Phi_4 = 4 / sum [f''(r_i)/f'(r_i)]^2

# Let's compute the discriminant. disc(f) = prod_{i<j} (r_i - r_j)^2
disc_expr = sp.Integer(1)
for i in range(n):
    for j in range(i+1, n):
        disc_expr *= (roots[i] - roots[j])**2

# Note: prod f'(r_i) = (-1)^{n(n-1)/2} * disc(f) = ... for n=4: (-1)^6 = 1
# So prod f'(r_i) = disc(f) ... no, disc(f) = prod_{i<j}(r_i-r_j)^2
# while prod f'(r_i) = prod_i prod_{j!=i}(r_i-r_j) = prod_{i<j}(r_i-r_j)^2 * (-1)^...
# Actually: prod_i f'(r_i) = prod_i prod_{j!=i}(r_i-r_j)
# Each pair (i,j) appears twice: as (r_i-r_j) in f'(r_i) and (r_j-r_i) in f'(r_j)
# So prod = prod_{i<j}(r_i-r_j)(r_j-r_i) * ... hmm this is complicated.
# Let's just not go there. Let me try to compute 1/Phi_4 differently.

print("\n--- Attempting symbolic computation via resultant/discriminant ---")

# For a degree-n polynomial, there's a known formula:
# sum_i 1/f'(r_i)^2 = [result involving disc and other invariants]
# But for 1/Phi_4 = 1/sum(H_i^2), we need something different.

# Let me try numerics first to understand the structure.
print("\n--- Switching to NUMERICAL exploration first ---")

# =====================================================================
# PART 2: Numerical computation of 1/Phi_4
# =====================================================================
print("\n" + "=" * 72)
print("PART 2: Numerical computation and verification")
print("=" * 72)

def elem_sym(roots, k):
    """Compute k-th elementary symmetric polynomial of roots."""
    n = len(roots)
    if k == 0: return 1.0
    if k > n: return 0.0
    return sum(np.prod(list(combo)) for combo in combinations(roots, k))

def Phi_n_numerical(roots):
    """Compute Phi_n = sum_i H(r_i)^2."""
    n = len(roots)
    total = 0.0
    for i in range(n):
        H_i = sum(1.0/(roots[i]-roots[j]) for j in range(n) if j != i)
        total += H_i**2
    return total

def boxplus_mss(roots_p, roots_q):
    """MSS boxplus: g_k = sum_{i+j=k} C(n-j,i)/C(n,i) * e_i(p) * e_j(q)."""
    n = len(roots_p)
    assert len(roots_q) == n
    ep = [elem_sym(roots_p, k) for k in range(n+1)]
    eq = [elem_sym(roots_q, k) for k in range(n+1)]
    g = [0.0]*(n+1)
    for k in range(n+1):
        for i in range(k+1):
            j = k - i
            if i <= n and j <= n and comb(n, i) > 0:
                w = comb(n-j, i) / comb(n, i)
                g[k] += w * ep[i] * eq[j]
    # Build polynomial coefficients: x^n - g1*x^{n-1} + g2*x^{n-2} - ...
    coeffs = [1.0]
    for k in range(1, n+1):
        coeffs.append((-1)**k * g[k])
    roots_r = np.sort(np.real(np.roots(coeffs)))
    return roots_r

# Test basic computation
print("\nTest 1: 1/Phi_4 for some specific quartics")
test_cases = [
    ("equally spaced", np.array([-3.0, -1.0, 1.0, 3.0])),
    ("asymmetric", np.array([-2.0, -0.5, 0.5, 2.0])),
    ("generic", np.array([-3.0, -1.0, 0.5, 3.5])),
    ("near equal gaps", np.array([-1.5, -0.5, 0.5, 1.5])),
]

for desc, roots_arr in test_cases:
    roots_arr = roots_arr - np.mean(roots_arr)  # center
    phi = Phi_n_numerical(roots_arr)
    inv_phi = 1.0/phi
    e2 = elem_sym(roots_arr, 2)
    e3 = elem_sym(roots_arr, 3)
    e4 = elem_sym(roots_arr, 4)
    print(f"  {desc:20s}: e2={e2:.4f}, e3={e3:.4f}, e4={e4:.4f}, "
          f"1/Phi_4={inv_phi:.6f}")

# =====================================================================
# PART 3: Verify MSS boxplus for n=4 centered
# =====================================================================
print("\n" + "=" * 72)
print("PART 3: MSS boxplus coefficients for n=4 centered")
print("=" * 72)

# For n=4, centered (e_1=0):
# g_k = sum_{i+j=k} C(4-j,i)/C(4,i) * e_i(p) * e_j(q)
print("Computing MSS coefficients for n=4:")
for k in range(5):
    print(f"\n  g_{k}:")
    for i in range(k+1):
        j = k - i
        if i <= 4 and j <= 4 and comb(4, i) > 0:
            w = Rational(comb(4-j, i), comb(4, i))
            print(f"    i={i},j={j}: C({4-j},{i})/C(4,{i}) = {comb(4-j,i)}/{comb(4,i)} = {w}  * e_{i}(p)*e_{j}(q)")

print("\nFor centered (e_1 = 0):")
print("  g_0 = 1")
print("  g_1 = 0  (only term with e_1)")
print("  g_2 = e_2(p) + e_2(q)  [cross term has e_1 = 0]")
print("  g_3 = e_3(p) + e_3(q)  [cross terms have e_1 = 0]")

# g_4 in detail:
# (i,j): (0,4), (1,3), (2,2), (3,1), (4,0)
# (0,4): C(0,0)/C(4,0) = 1  -> e_4(q)
# (1,3): C(1,1)/C(4,1) = 1/4 -> e_1(p)*e_3(q) = 0
# (2,2): C(2,2)/C(4,2) = 1/6 -> e_2(p)*e_2(q)
# (3,1): C(3,3)/C(4,3) = 1/4 -> e_3(p)*e_1(q) = 0
# (4,0): C(4,4)/C(4,4) = 1  -> e_4(p)

print("  g_4 = e_4(p) + e_4(q) + (1/6)*e_2(p)*e_2(q)")
print("        [cross terms with e_1 vanish, but the e_2*e_2 term survives!]")

# Verify numerically
print("\nNumerical verification of boxplus formula:")
np.random.seed(42)
for trial in range(5):
    rp = np.sort(np.random.randn(4)*2)
    rq = np.sort(np.random.randn(4)*2)
    while np.min(np.diff(rp)) < 0.2: rp = np.sort(np.random.randn(4)*2)
    while np.min(np.diff(rq)) < 0.2: rq = np.sort(np.random.randn(4)*2)
    rp -= np.mean(rp); rq -= np.mean(rq)
    rr = boxplus_mss(rp, rq)

    ep2, ep3, ep4 = elem_sym(rp, 2), elem_sym(rp, 3), elem_sym(rp, 4)
    eq2, eq3, eq4 = elem_sym(rq, 2), elem_sym(rq, 3), elem_sym(rq, 4)
    er2, er3, er4 = elem_sym(rr, 2), elem_sym(rr, 3), elem_sym(rr, 4)

    g2_pred = ep2 + eq2
    g3_pred = ep3 + eq3
    g4_pred = ep4 + eq4 + (1.0/6)*ep2*eq2

    print(f"  Trial {trial}: e2 err={abs(er2-g2_pred):.1e}, "
          f"e3 err={abs(er3-g3_pred):.1e}, e4 err={abs(er4-g4_pred):.1e}")

# =====================================================================
# PART 4: Symbolic computation of 1/Phi_4 via coefficient parameterization
# =====================================================================
print("\n" + "=" * 72)
print("PART 4: Symbolic 1/Phi_4 in terms of (e2, e3, e4)")
print("=" * 72)

# For f(x) = x^4 + e2*x^2 + e3*x + e4 (centered quartic):
# f'(x) = 4x^3 + 2*e2*x + e3
# f''(x) = 12x^2 + 2*e2
# H(r_i) = f''(r_i)/(2*f'(r_i)) = (12*r_i^2 + 2*e2)/(2*(4*r_i^3 + 2*e2*r_i + e3))
#         = (6*r_i^2 + e2)/(4*r_i^3 + 2*e2*r_i + e3)
# But f(r_i)=0 means r_i^4 = -e2*r_i^2 - e3*r_i - e4
# So f'(r_i) = 4*r_i^3 + 2*e2*r_i + e3

# Phi_4 = sum_i H_i^2 = sum_i [(6*r_i^2 + e2) / (4*r_i^3 + 2*e2*r_i + e3)]^2

# The key identity: Phi_4 can be expressed as a rational function of (e2, e3, e4)
# via Newton's identities relating power sums to elementary symmetric polynomials.

# However, a direct symbolic computation is feasible using resultants.

# Alternative approach: compute Phi_4 numerically for many points and
# try to identify the rational function.

# Let's try: compute sum H_i^2 * (prod f'(r_i))^2 as a symmetric polynomial
# and express in terms of e2, e3, e4.

print("Attempting symbolic computation of 1/Phi_4...")
print("Strategy: compute Phi_4 as rational function of (e2, e3, e4)")
sys.stdout.flush()

e2_s, e3_s, e4_s = symbols('e2 e3 e4')

# For centered quartic: f(x) = x^4 + e2*x^2 + e3*x + e4
# f'(x) = 4x^3 + 2*e2*x + e3
# f''(x) = 12x^2 + 2*e2

# Power sums p_k = sum r_i^k can be expressed via Newton's identities:
# p_1 = 0 (centered)
# p_2 = -2*e2
# p_3 = -3*e3
# p_4 = 2*e2^2 - 4*e4
# p_5 = 5*e2*e3
# p_6 = -2*e2^3 + 6*e2*e4 + 3*e3^2

# We need: sum H_i^2 = sum [(6*r_i^2 + e2)^2 / f'(r_i)^2]
# This is NOT a symmetric function directly... unless we write it as
# sum numerator_i^2 / prod f'(r_j)^2 * prod_{j!=i} f'(r_j)^2

# Actually, let's use a cleaner approach.
# sum H_i^2 = sum [f''(r_i)]^2 / [4 f'(r_i)^2]

# Let S = sum [f''(r_i)]^2 / f'(r_i)^2
# Then Phi_4 = S/4, and 1/Phi_4 = 4/S.

# S = sum [f''(r_i)/f'(r_i)]^2

# Now f''(r_i)/f'(r_i) = d/dx[ln f'(x)]|_{x=r_i}
# But f'(x) = 4*prod(x - s_j) where s_j are roots of f'(x).
# So f''(x)/f'(x) = sum_j 1/(x - s_j).

# Hmm, this doesn't simplify easily. Let me just compute numerically
# and try to identify patterns.

print("\n--- Numerical identification of 1/Phi_4 formula ---")

def inv_Phi_from_roots(roots):
    """Compute 1/Phi_n from roots."""
    return 1.0 / Phi_n_numerical(roots)

# Let's parameterize centered quartics by (e2, e3, e4) and
# try to understand the functional form.

# First question: for e3 = 0 (symmetric quartics), what is 1/Phi_4?
print("\nCase 1: e3 = 0 (palindromic quartics, roots {-a, -b, b, a})")
for a, b in [(3, 1), (2, 1), (4, 1), (3, 2), (5, 1), (2, 0.5)]:
    roots_arr = np.array([-a, -b, b, a])
    e2 = elem_sym(roots_arr, 2)
    e4 = elem_sym(roots_arr, 4)
    inv_phi = inv_Phi_from_roots(roots_arr)
    disc = 256*e4**3 - 128*e2**2*e4**2 + 144*e2*0**2*e4 - 27*0**4 + 16*e2**4*e4 - 4*e2**3*0**2
    # disc for e3=0: 256*e4^3 - 128*e2^2*e4^2 + 16*e2^4*e4
    disc2 = 16*e4*(16*e4**2 - 8*e2**2*e4 + e2**4)
    print(f"  a={a:4.1f}, b={b:4.1f}: e2={e2:8.3f}, e4={e4:8.3f}, "
          f"1/Phi_4={inv_phi:10.6f}")

# =====================================================================
# PART 5: Compute 1/Phi_4 using resultant/subresultant technique
# =====================================================================
print("\n" + "=" * 72)
print("PART 5: Exact symbolic computation of Phi_4")
print("=" * 72)

# Strategy: use the fact that for f(x) = x^4 + e2*x^2 + e3*x + e4,
# the sum S = sum_i [f''(r_i)]^2 / [f'(r_i)]^2 can be computed via
# the formula: sum_i g(r_i)/f'(r_i) = [resultant-related computation]

# Actually, use partial fractions:
# g(x)/f(x) = sum_i g(r_i)/((x-r_i)*f'(r_i))
# The sum sum_i g(r_i)/f'(r_i) is the coefficient of 1/x in the Laurent
# expansion of g(x)/f(x) at infinity, multiplied by the leading coefficient.
# If deg(g) < deg(f) = n, then sum_i g(r_i)/f'(r_i) = [x^{n-1}] coefficient
# in the polynomial division... actually it's the negative of the residue sum.

# For our problem:
# sum_i [f''(r_i)]^2 / f'(r_i)^2 is harder because it's squared.
#
# Let's use: sum_i h(r_i)^2 = sum_i h(r_i)*h(r_i)
# where h(r_i) = f''(r_i)/f'(r_i).
#
# We can compute sum_i p(r_i)/f'(r_i) for any polynomial p using resultants.
# But h(r_i)^2 involves 1/f'(r_i)^2.

# Alternative: direct Newton identity computation.
# H_i = f''(r_i)/(2*f'(r_i))
# = (12*r_i^2 + 2*e2) / (2*(4*r_i^3 + 2*e2*r_i + e3))
# = (6*r_i^2 + e2) / (4*r_i^3 + 2*e2*r_i + e3)

# Using f(r_i) = 0: r_i^4 = -e2*r_i^2 - e3*r_i - e4

# The denominator is f'(r_i). We know that:
# 1/Phi_4 = 4 / sum_i [f''(r_i)/f'(r_i)]^2

# Let's compute sum_i [f''(r_i)/f'(r_i)]^2 using the "trace of the square"
# of the matrix [f''(r_i)/f'(r_i)].

# Actually, the cleanest approach: compute everything via power sums.

# Key formula: sum_i r_i^a / f'(r_i)^2 can be computed!
# We use: sum_i r_i^a / f'(r_i)^2 involves the discriminant.

# Actually, let me use a different approach. Let me compute via sympy's
# polynomial resultant machinery.

print("Computing via polynomial arithmetic modulo f(x)...")
sys.stdout.flush()

x_s = symbols('x')
f_s = x_s**4 + e2_s*x_s**2 + e3_s*x_s + e4_s
fp_s = sp.diff(f_s, x_s)  # 4x^3 + 2*e2*x + e3
fpp_s = sp.diff(f_s, x_s, 2)  # 12x^2 + 2*e2

# We need sum_i [fpp(r_i)]^2 / fp(r_i)^2

# Let's compute sum_i g(r_i)/fp(r_i) for polynomial g using:
# g(x)/f(x) has partial fraction sum_i g(r_i)/((x-r_i)*fp(r_i))
# So sum_i g(r_i)/fp(r_i) = [coefficient of x^{-1} in g(x)/f(x)]
# = coefficient of x^{n-1} in the Euclidean division...
# No, for deg(g) < deg(f), the partial fraction sum is:
# sum_i g(r_i)/fp(r_i) = sum of residues of g(x)/f(x) at roots of f
# By the residue theorem (or direct computation), if deg(g) < deg(f) = n:
# sum_i g(r_i)/fp(r_i) = 0 if deg(g) < n-1
# sum_i g(r_i)/fp(r_i) = leading_coeff(g)/leading_coeff(f) if deg(g) = n-1

# More precisely: if g(x)/f(x) = q(x) + r(x)/f(x) with deg(r) < n,
# then sum_i g(r_i)/fp(r_i) is the sum of residues = [x^{n-1}] coefficient of
# (g mod f) divided by leading coeff of f...

# For our problem, we need sum_i [fpp(r_i)]^2 / fp(r_i)^2.
# This is NOT sum g/fp but sum g/fp^2. Different beast.

# Key identity:
# sum_i g(r_i)/fp(r_i)^2 = -d/de4 [sum_i g(r_i)/fp(r_i)]
# because d(fp(r_i))/de4 = d/de4[4*r_i^3+2*e2*r_i+e3] = 0 (fp doesn't depend on e4)
# BUT r_i depends on e4! So dr_i/de4 = -1/fp(r_i) (implicit function theorem).
# Then d/de4[g(r_i)/fp(r_i)] = [g'(r_i)*(-1/fp(r_i))*fp(r_i) - g(r_i)*fpp(r_i)*(-1/fp(r_i))]/ fp(r_i)^2
# Hmm, this is getting complicated.

# Let me just use the DIRECT resultant approach.
# Given f(x), compute:
# R(y) = Res_x(f(x), y*fp(x) - fpp(x))
# The roots of R(y) are y_i = fpp(r_i)/fp(r_i) = 2*H(r_i).
# Then Phi_4 = (1/4) * sum y_i^2 = (1/4) * (p_1^2 - 2*p_2)
# where p_k are power sums of y_i.
# And p_1, p_2 are related to coefficients of R(y) via Newton's identities.

print("Computing resultant R(y) = Res_x(f(x), y*fp(x) - fpp(x))...")
sys.stdout.flush()

y_s = symbols('y')
g_s = y_s * fp_s - fpp_s  # y*(4x^3+2*e2*x+e3) - (12x^2+2*e2)

# R(y) = resultant of f(x) and g(x) with respect to x
# This is a polynomial in y whose roots are fpp(r_i)/fp(r_i) for i=1..4

t0 = time.time()
R_y = resultant(f_s, g_s, x_s)
t1 = time.time()
print(f"Resultant computed in {t1-t0:.1f}s")

R_y_expanded = expand(R_y)
R_poly = Poly(R_y_expanded, y_s)
print(f"R(y) degree: {R_poly.degree()}")
print(f"R(y) coefficients (from leading to constant):")
for i, c in enumerate(R_poly.all_coeffs()):
    cf = factor(c)
    print(f"  y^{R_poly.degree()-i}: {cf}")

# The roots y_i = fpp(r_i)/fp(r_i) = 2*H(r_i)
# So sum y_i^2 = 4*sum H_i^2 = 4*Phi_4
# sum y_i = -coeff(y^3)/coeff(y^4) from Vieta's
# sum y_i^2 = (sum y_i)^2 - 2*sum_{i<j} y_i*y_j
# sum_{i<j} y_i*y_j = coeff(y^2)/coeff(y^4) from Vieta's

coeffs_R = R_poly.all_coeffs()
# R(y) = a4*y^4 + a3*y^3 + a2*y^2 + a1*y + a0
deg_R = R_poly.degree()
if deg_R == 4:
    a4 = coeffs_R[0]
    a3 = coeffs_R[1]
    a2 = coeffs_R[2]
    a1 = coeffs_R[3]
    a0 = coeffs_R[4]

    # Vieta's: sum y_i = -a3/a4
    # sum_{i<j} y_i*y_j = a2/a4
    # sum y_i^2 = (-a3/a4)^2 - 2*(a2/a4) = (a3^2 - 2*a4*a2)/a4^2

    sum_y_sq = cancel((a3**2 - 2*a4*a2) / a4**2)
    # 4*Phi_4 = sum y_i^2
    # Phi_4 = sum_y_sq / 4
    # 1/Phi_4 = 4 / sum_y_sq

    Phi_4_formula = cancel(sum_y_sq / 4)
    inv_Phi_4_formula = cancel(4 / sum_y_sq)

    num_inv_Phi, den_inv_Phi = sp.fraction(inv_Phi_4_formula)

    print(f"\nsum y_i^2 = (a3^2 - 2*a4*a2) / a4^2")
    print(f"\nPhi_4 = sum_y_sq/4 = {factor(sp.fraction(Phi_4_formula)[0])} / {factor(sp.fraction(Phi_4_formula)[1])}")
    print(f"\n1/Phi_4 = 4/sum_y_sq")

    num_f = factor(num_inv_Phi)
    den_f = factor(den_inv_Phi)
    print(f"\n1/Phi_4 = {num_f}")
    print(f"          / {den_f}")
else:
    print(f"Unexpected degree {deg_R} for R(y)")
    sum_y_sq = None

# =====================================================================
# PART 6: Verify 1/Phi_4 formula numerically
# =====================================================================
print("\n" + "=" * 72)
print("PART 6: Numerical verification of 1/Phi_4 formula")
print("=" * 72)

if sum_y_sq is not None:
    inv_Phi_4_func = sp.lambdify((e2_s, e3_s, e4_s), inv_Phi_4_formula, 'numpy')

    print("\nVerifying formula against direct root computation:")
    np.random.seed(123)
    all_match = True
    for trial in range(20):
        roots_arr = np.sort(np.random.randn(4)*2)
        roots_arr -= np.mean(roots_arr)
        # Ensure distinct roots
        while np.min(np.diff(roots_arr)) < 0.3:
            roots_arr = np.sort(np.random.randn(4)*2)
            roots_arr -= np.mean(roots_arr)

        e2_v = elem_sym(roots_arr, 2)
        e3_v = elem_sym(roots_arr, 3)
        e4_v = elem_sym(roots_arr, 4)

        inv_phi_direct = 1.0 / Phi_n_numerical(roots_arr)
        inv_phi_formula_val = float(inv_Phi_4_func(e2_v, e3_v, e4_v))

        match = abs(inv_phi_direct - inv_phi_formula_val) < 1e-8 * abs(inv_phi_direct)
        if not match:
            all_match = False
            print(f"  Trial {trial}: MISMATCH! direct={inv_phi_direct:.10f}, "
                  f"formula={inv_phi_formula_val:.10f}")
        elif trial < 5:
            print(f"  Trial {trial}: direct={inv_phi_direct:.10f}, "
                  f"formula={inv_phi_formula_val:.10f} OK")

    if all_match:
        print("  ... all 20 trials match!")
    else:
        print("  FORMULA VERIFICATION FAILED!")

# =====================================================================
# PART 7: Compute the excess for n=4
# =====================================================================
print("\n" + "=" * 72)
print("PART 7: Compute the excess 1/Phi_4(r) - 1/Phi_4(p) - 1/Phi_4(q)")
print("=" * 72)

if sum_y_sq is not None:
    # Parameters for p and q
    a2, a3, a4_p = symbols('a2 a3 a4', real=True)  # e2(p), e3(p), e4(p)
    b2, b3, b4_q = symbols('b2 b3 b4', real=True)  # e2(q), e3(q), e4(q)

    # MSS boxplus for n=4 centered:
    r2 = a2 + b2           # e2(r) = e2(p) + e2(q)
    r3 = a3 + b3           # e3(r) = e3(p) + e3(q)
    r4 = a4_p + b4_q + Rational(1,6)*a2*b2  # e4(r) = e4(p) + e4(q) + (1/6)*e2(p)*e2(q)

    # 1/Phi_4 as function of (e2, e3, e4)
    inv_Phi_p = inv_Phi_4_formula.subs({e2_s: a2, e3_s: a3, e4_s: a4_p})
    inv_Phi_q = inv_Phi_4_formula.subs({e2_s: b2, e3_s: b3, e4_s: b4_q})
    inv_Phi_r = inv_Phi_4_formula.subs({e2_s: r2, e3_s: r3, e4_s: r4})

    print("Computing excess = 1/Phi_4(r) - 1/Phi_4(p) - 1/Phi_4(q)...")
    print("This is a rational function of 6 variables (a2,a3,a4,b2,b3,b4).")
    sys.stdout.flush()

    # This will be a complex rational expression. Let's first simplify
    # in special cases.

    # Special case 1: a3 = b3 = 0 (symmetric quartics)
    print("\n--- Special Case 1: a3 = b3 = 0 (symmetric quartics) ---")
    sys.stdout.flush()

    inv_Phi_p_sym = inv_Phi_4_formula.subs({e2_s: a2, e3_s: 0, e4_s: a4_p})
    inv_Phi_q_sym = inv_Phi_4_formula.subs({e2_s: b2, e3_s: 0, e4_s: b4_q})
    inv_Phi_r_sym = inv_Phi_4_formula.subs({e2_s: a2+b2, e3_s: 0,
                                             e4_s: a4_p + b4_q + Rational(1,6)*a2*b2})

    excess_sym = together(inv_Phi_r_sym - inv_Phi_p_sym - inv_Phi_q_sym)
    num_excess_sym, den_excess_sym = sp.fraction(excess_sym)
    num_excess_sym = expand(num_excess_sym)

    print(f"Excess numerator degree in (a4, b4): checking...")
    num_poly_a4 = Poly(num_excess_sym, a4_p, b4_q)
    print(f"  Total degree in (a4, b4): {num_poly_a4.total_degree()}")
    print(f"Excess denominator: {factor(den_excess_sym)}")

    # Collect as polynomial in a4, b4
    print("\nExcess as polynomial in (a4, b4) for fixed a2, b2:")
    for mon, coeff in sorted(num_poly_a4.as_dict().items()):
        coeff_f = factor(coeff)
        print(f"  a4^{mon[0]} * b4^{mon[1]}: {coeff_f}")

    # Special case 2: p = q (self-convolution)
    print("\n--- Special Case 2: p = q (self-convolution) ---")
    sys.stdout.flush()

    E2, E3, E4 = symbols('E2 E3 E4', real=True)

    inv_Phi_single = inv_Phi_4_formula.subs({e2_s: E2, e3_s: E3, e4_s: E4})
    # r = p boxplus p: e2(r) = 2*E2, e3(r) = 2*E3, e4(r) = 2*E4 + (1/6)*E2^2
    inv_Phi_double = inv_Phi_4_formula.subs({e2_s: 2*E2, e3_s: 2*E3,
                                              e4_s: 2*E4 + Rational(1,6)*E2**2})

    excess_self = together(inv_Phi_double - 2*inv_Phi_single)
    num_excess_self, den_excess_self = sp.fraction(excess_self)
    num_excess_self = expand(num_excess_self)

    print(f"Numerator of self-conv excess (factored):")
    nf = factor(num_excess_self)
    print(f"  {nf}")
    print(f"Denominator: {factor(den_excess_self)}")

# =====================================================================
# PART 8: Extensive numerical testing
# =====================================================================
print("\n" + "=" * 72)
print("PART 8: Extensive numerical testing of Fisher superadditivity at n=4")
print("=" * 72)

np.random.seed(42)

def generate_centered_quartic():
    """Generate random centered quartic with distinct real roots."""
    while True:
        roots = np.sort(np.random.randn(4) * np.random.uniform(0.5, 5))
        roots -= np.mean(roots)
        if np.min(np.diff(roots)) > 0.01:
            return roots

n_trials = 100000
violations = 0
min_excess = float('inf')
min_excess_info = None
excesses = []

print(f"Running {n_trials} random trials...")
sys.stdout.flush()
t0 = time.time()

for trial in range(n_trials):
    rp = generate_centered_quartic()
    rq = generate_centered_quartic()

    rr = boxplus_mss(rp, rq)

    inv_p = 1.0/Phi_n_numerical(rp)
    inv_q = 1.0/Phi_n_numerical(rq)
    inv_r = 1.0/Phi_n_numerical(rr)

    excess = inv_r - inv_p - inv_q
    excesses.append(excess)

    if excess < min_excess:
        min_excess = excess
        min_excess_info = (rp.copy(), rq.copy(), rr.copy(), inv_p, inv_q, inv_r)

    if excess < -1e-10:
        violations += 1

t1 = time.time()
print(f"Time: {t1-t0:.1f}s")
print(f"Trials: {n_trials}")
print(f"Violations: {violations}")
print(f"Min excess: {min_excess:.6e}")

if min_excess_info:
    rp, rq, rr, inv_p, inv_q, inv_r = min_excess_info
    print(f"\nMinimum excess case:")
    print(f"  p roots: {rp}")
    print(f"  q roots: {rq}")
    print(f"  r roots: {rr}")
    print(f"  1/Phi(p)={inv_p:.8f}, 1/Phi(q)={inv_q:.8f}, 1/Phi(r)={inv_r:.8f}")
    print(f"  excess = {inv_r - inv_p - inv_q:.6e}")
    print(f"  e2(p)={elem_sym(rp,2):.4f}, e3(p)={elem_sym(rp,3):.4f}, e4(p)={elem_sym(rp,4):.4f}")
    print(f"  e2(q)={elem_sym(rq,2):.4f}, e3(q)={elem_sym(rq,3):.4f}, e4(q)={elem_sym(rq,4):.4f}")

# Self-convolution test
print("\n--- Self-convolution test: 1/Phi(p⊞p) >= 2/Phi(p) ---")
n_self = 50000
violations_self = 0
min_excess_self = float('inf')

for trial in range(n_self):
    rp = generate_centered_quartic()
    rr = boxplus_mss(rp, rp)
    excess = 1.0/Phi_n_numerical(rr) - 2.0/Phi_n_numerical(rp)
    if excess < min_excess_self:
        min_excess_self = excess
    if excess < -1e-10:
        violations_self += 1

print(f"Self-conv trials: {n_self}")
print(f"Self-conv violations: {violations_self}")
print(f"Self-conv min excess: {min_excess_self:.6e}")

# Equally-spaced test
print("\n--- Equally-spaced roots: p = (-3d, -d, d, 3d) ---")
for d1, d2 in [(1, 1), (1, 2), (0.5, 3), (2, 0.7), (1, 0.1)]:
    rp = np.array([-3*d1, -d1, d1, 3*d1])
    rq = np.array([-3*d2, -d2, d2, 3*d2])
    rr = boxplus_mss(rp, rq)
    excess = 1.0/Phi_n_numerical(rr) - 1.0/Phi_n_numerical(rp) - 1.0/Phi_n_numerical(rq)
    print(f"  d_p={d1}, d_q={d2}: excess = {excess:.6e}")

# =====================================================================
# PART 9: Analyze what controls the minimum
# =====================================================================
print("\n" + "=" * 72)
print("PART 9: Structure of the excess — what controls the minimum?")
print("=" * 72)

# Fix e2(p) and e2(q), vary e3 and e4
# For centered quartic with roots -a, -b, b, a:
# e2 = -(a^2+b^2), e3 = 0, e4 = a^2*b^2

print("\n--- Varying e3 with fixed e2, e4 ---")
# Create quartic with specific (e2, e3, e4) by finding roots
# x^4 + e2*x^2 + e3*x + e4 = 0

# For p: fix e2(p) = -5 (so sum of squared root-gaps > 0), vary e3(p)
# For q: fix similarly

# Use equally-spaced as base, then perturb
base_p = np.array([-3.0, -1.0, 1.0, 3.0])  # e2=-10, e3=0, e4=9
base_q = np.array([-2.0, -0.5, 0.5, 2.0])  # different scale

print("Perturbing p by shifting one root (breaking symmetry):")
for shift in [0, 0.1, 0.5, 1.0, 2.0]:
    rp = base_p.copy()
    rp[1] += shift  # shift second root
    rp -= np.mean(rp)  # recenter
    rq = base_q.copy()
    rr = boxplus_mss(rp, rq)
    excess = 1.0/Phi_n_numerical(rr) - 1.0/Phi_n_numerical(rp) - 1.0/Phi_n_numerical(rq)
    e3_p = elem_sym(rp, 3)
    print(f"  shift={shift:.1f}, e3(p)={e3_p:8.4f}: excess = {excess:.6e}")

# =====================================================================
# PART 10: Try to understand 1/Phi_4 structure better
# =====================================================================
print("\n" + "=" * 72)
print("PART 10: Decompose 1/Phi_4 in terms of discriminant and other invariants")
print("=" * 72)

# For n=3: 1/Phi_3 = disc/(18*e2^2) = (-4e2^3 - 27*e3^2)/(18*e2^2)
# Can we find a similar relation for n=4?

# The discriminant of x^4 + px^2 + qx + r is:
# disc = 256*r^3 - 128*p^2*r^2 + 144*p*q^2*r - 27*q^4 + 16*p^4*r - 4*p^3*q^2

print("Computing discriminant relation...")
disc_4 = 256*e4_s**3 - 128*e2_s**2*e4_s**2 + 144*e2_s*e3_s**2*e4_s \
         - 27*e3_s**4 + 16*e2_s**4*e4_s - 4*e2_s**3*e3_s**2

if sum_y_sq is not None:
    # Check if 1/Phi_4 = disc_4 / (something)
    # Numerically test the ratio disc_4 / (1/Phi_4) for several quartics
    disc_4_func = sp.lambdify((e2_s, e3_s, e4_s), disc_4, 'numpy')

    print("\nTesting disc_4 / inv_Phi_4 ratio:")
    for trial in range(10):
        roots_arr = generate_centered_quartic()
        e2_v = elem_sym(roots_arr, 2)
        e3_v = elem_sym(roots_arr, 3)
        e4_v = elem_sym(roots_arr, 4)
        inv_phi = 1.0 / Phi_n_numerical(roots_arr)
        disc_val = disc_4_func(e2_v, e3_v, e4_v)
        if abs(inv_phi) > 1e-15:
            ratio = disc_val / inv_phi
            print(f"  Trial {trial}: disc/invPhi = {ratio:.4f}, e2={e2_v:.3f}, "
                  f"e3={e3_v:.3f}, e4={e4_v:.3f}")

# =====================================================================
# PART 11: Print the exact 1/Phi_4 formula
# =====================================================================
print("\n" + "=" * 72)
print("PART 11: Exact 1/Phi_4 formula")
print("=" * 72)

if sum_y_sq is not None:
    print(f"\n1/Phi_4 numerator (factored): {factor(num_inv_Phi)}")
    print(f"\n1/Phi_4 denominator (factored): {factor(den_inv_Phi)}")

    # Collect numerator in terms of e3
    num_collected = collect(expand(num_inv_Phi), e3_s)
    print(f"\n1/Phi_4 numerator collected by e3:")
    print(f"  {num_collected}")

    # Check: is numerator = discriminant?
    diff_disc = expand(num_inv_Phi - disc_4)
    if diff_disc == 0:
        print("\n*** NUMERATOR OF 1/Phi_4 = DISCRIMINANT! ***")
    else:
        # Check if proportional
        # Try ratio
        ratio_sym = cancel(num_inv_Phi / disc_4)
        print(f"\nRatio num(1/Phi_4)/disc_4 = {ratio_sym}")

# =====================================================================
# PART 12: Analyze the excess in the symmetric case more carefully
# =====================================================================
print("\n" + "=" * 72)
print("PART 12: Deeper analysis of the excess")
print("=" * 72)

if sum_y_sq is not None:
    # Let's try the full excess computation
    print("Computing full excess (6-variable rational function)...")
    print("This may take a while...")
    sys.stdout.flush()

    t0 = time.time()
    try:
        excess_full = together(inv_Phi_r - inv_Phi_p - inv_Phi_q)
        t1 = time.time()
        print(f"Computed in {t1-t0:.1f}s")

        num_excess, den_excess = sp.fraction(excess_full)
        num_excess = expand(num_excess)
        den_excess = expand(den_excess)

        print(f"\nDenominator factored: {factor(den_excess)}")

        # Count terms in numerator
        terms = num_excess.as_ordered_terms()
        print(f"Numerator has {len(terms)} terms")

        # Is the denominator always positive?
        den_f = factor(den_excess)
        print(f"Denominator = {den_f}")

    except Exception as ex:
        print(f"Full excess computation failed: {ex}")
        print("Falling back to numerical analysis only.")

# =====================================================================
# PART 13: The key question - can excess be negative?
# =====================================================================
print("\n" + "=" * 72)
print("PART 13: Summary of numerical findings")
print("=" * 72)

excesses_arr = np.array(excesses)
print(f"\nExcess statistics over {len(excesses_arr)} trials:")
print(f"  Min:    {excesses_arr.min():.6e}")
print(f"  Max:    {excesses_arr.max():.6e}")
print(f"  Mean:   {excesses_arr.mean():.6e}")
print(f"  Median: {np.median(excesses_arr):.6e}")
print(f"  Std:    {excesses_arr.std():.6e}")
print(f"  Count < 0: {np.sum(excesses_arr < 0)}")
print(f"  Count < -1e-10: {np.sum(excesses_arr < -1e-10)}")

# Test with very extreme parameters
print("\n--- Extreme parameter tests ---")
extreme_violations = 0
n_extreme = 50000
min_extreme = float('inf')

for trial in range(n_extreme):
    # Mix very different scales
    scale_p = 10**np.random.uniform(-2, 2)
    scale_q = 10**np.random.uniform(-2, 2)
    rp = generate_centered_quartic() * scale_p
    rq = generate_centered_quartic() * scale_q

    try:
        rr = boxplus_mss(rp, rq)
        if not np.all(np.isreal(rr)) or np.any(np.isnan(rr)):
            continue
        excess = 1.0/Phi_n_numerical(rr) - 1.0/Phi_n_numerical(rp) - 1.0/Phi_n_numerical(rq)
        if np.isfinite(excess):
            if excess < min_extreme:
                min_extreme = excess
            if excess < -1e-10:
                extreme_violations += 1
    except:
        continue

print(f"Extreme trials: {n_extreme}")
print(f"Extreme violations: {extreme_violations}")
print(f"Extreme min excess: {min_extreme:.6e}")

# Near-degenerate test (roots almost colliding)
print("\n--- Near-degenerate tests (roots close together) ---")
n_degen = 20000
violations_degen = 0
min_degen = float('inf')

for trial in range(n_degen):
    eps = 10**np.random.uniform(-6, -1)
    rp = np.array([-2, -eps, eps, 2.0])
    rp -= np.mean(rp)
    rq = generate_centered_quartic()

    try:
        rr = boxplus_mss(rp, rq)
        excess = 1.0/Phi_n_numerical(rr) - 1.0/Phi_n_numerical(rp) - 1.0/Phi_n_numerical(rq)
        if np.isfinite(excess):
            if excess < min_degen:
                min_degen = excess
            if excess < -1e-10:
                violations_degen += 1
    except:
        continue

print(f"Near-degenerate trials: {n_degen}")
print(f"Near-degenerate violations: {violations_degen}")
print(f"Near-degenerate min excess: {min_degen:.6e}")

# =====================================================================
# PART 14: Analyze when excess is near zero
# =====================================================================
print("\n" + "=" * 72)
print("PART 14: When is the excess near zero?")
print("=" * 72)

# Find cases where excess is very small
print("\nSearching for near-equality cases...")
near_zero_cases = []
np.random.seed(99)

for trial in range(200000):
    rp = generate_centered_quartic()
    rq = generate_centered_quartic()

    try:
        rr = boxplus_mss(rp, rq)
        excess = 1.0/Phi_n_numerical(rr) - 1.0/Phi_n_numerical(rp) - 1.0/Phi_n_numerical(rq)
        inv_phi_p = 1.0/Phi_n_numerical(rp)
        rel_excess = excess / (inv_phi_p + 1.0/Phi_n_numerical(rq))

        if rel_excess < 0.001:  # within 0.1% of equality
            near_zero_cases.append((rp.copy(), rq.copy(), excess, rel_excess))
    except:
        continue

print(f"Found {len(near_zero_cases)} near-equality cases (rel excess < 0.1%)")

if near_zero_cases:
    # Sort by relative excess
    near_zero_cases.sort(key=lambda x: x[3])
    print("\nTop 10 nearest-to-equality:")
    for i, (rp, rq, exc, rel_exc) in enumerate(near_zero_cases[:10]):
        e2_p, e3_p = elem_sym(rp, 2), elem_sym(rp, 3)
        e2_q, e3_q = elem_sym(rq, 2), elem_sym(rq, 3)
        e4_p, e4_q = elem_sym(rp, 4), elem_sym(rq, 4)

        # Check if roots are equally spaced
        gaps_p = np.diff(rp)
        gaps_q = np.diff(rq)
        gap_var_p = np.std(gaps_p)/np.mean(gaps_p) if np.mean(gaps_p) > 0 else 0
        gap_var_q = np.std(gaps_q)/np.mean(gaps_q) if np.mean(gaps_q) > 0 else 0

        print(f"  {i}: rel_exc={rel_exc:.2e}, gap_var_p={gap_var_p:.4f}, "
              f"gap_var_q={gap_var_q:.4f}, e3_p={e3_p:.4f}, e3_q={e3_q:.4f}")

# Check: is equality achieved when both have equally-spaced roots?
print("\n--- Equality case: equally-spaced roots ---")
for d1 in [0.5, 1.0, 2.0]:
    for d2 in [0.5, 1.0, 2.0]:
        rp = np.array([-3*d1, -d1, d1, 3*d1])
        rq = np.array([-3*d2, -d2, d2, 3*d2])
        rr = boxplus_mss(rp, rq)
        excess = 1.0/Phi_n_numerical(rr) - 1.0/Phi_n_numerical(rp) - 1.0/Phi_n_numerical(rq)
        print(f"  d_p={d1}, d_q={d2}: excess = {excess:.6e}, "
              f"e3(p)={elem_sym(rp,3):.1e}, e4(p)={elem_sym(rp,4):.4f}")

# Check with symmetric but not equally-spaced
print("\n--- Symmetric (e3=0) but not equally-spaced ---")
for a, b in [(3.0, 1.0), (2.0, 0.5), (4.0, 1.0), (3.0, 0.1)]:
    rp = np.array([-a, -b, b, a])
    for c, d_val in [(2.0, 0.5), (3.0, 1.0), (1.0, 0.3)]:
        rq = np.array([-c, -d_val, d_val, c])
        rr = boxplus_mss(rp, rq)
        excess = 1.0/Phi_n_numerical(rr) - 1.0/Phi_n_numerical(rp) - 1.0/Phi_n_numerical(rq)
        print(f"  p=({a},{b}), q=({c},{d_val}): excess = {excess:.6e}")

# =====================================================================
# FINAL SUMMARY
# =====================================================================
print("\n" + "=" * 72)
print("FINAL SUMMARY")
print("=" * 72)
print(f"""
NUMERICAL RESULTS FOR n=4 FISHER SUPERADDITIVITY:
==================================================
Total random trials: {n_trials + n_extreme + n_degen}
Total violations: {violations + extreme_violations + violations_degen}
Minimum excess found: {min(min_excess, min_extreme, min_degen):.6e}

CONCLUSION: Fisher superadditivity 1/Phi_4(p⊞q) >= 1/Phi_4(p) + 1/Phi_4(q)
appears to hold for ALL centered quartics tested.

KEY DIFFERENCES FROM n=3:
1. MSS boxplus is NOT purely additive: g_4 has cross term (1/6)*e2(p)*e2(q)
2. 1/Phi_4 is a more complex rational function of (e2, e3, e4)
3. The excess is NOT a simple quadratic form in free parameters

The coefficient-level approach that gave an elegant proof for n=3
is significantly more complex for n=4 due to the cross term in g_4
and the higher complexity of the 1/Phi_4 formula.
""")
