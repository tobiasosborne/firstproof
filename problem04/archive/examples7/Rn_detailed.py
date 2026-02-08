"""
PROVER-10: Detailed computation of 1/Phi_n in terms of finite free cumulants.

Key discovery from Part 1: The simple regression C_n = 2/(n(n-1)) fails for n>=4.
This is because the decomposition 1/Phi_n = C_n*k2 + R_n is as a RATIONAL function,
not polynomial. The "C_n" we seek is the coefficient of k2 in the Laurent expansion
of 1/Phi_n around k3=...=kn=0.

Strategy: For each n, compute 1/Phi_n as a rational function of cumulants exactly.
"""

import numpy as np
from itertools import combinations
from math import comb, factorial
from fractions import Fraction
import sympy as sp
from sympy import (symbols, Rational, expand, simplify, factor, cancel,
                   together, collect, Poly, groebner, degree, Symbol,
                   numer, denom, sqrt, Integer, oo, series, O)

np.random.seed(42)

# ============================================================
# PART A: Exact symbolic computation via SymPy for n=2,3,4
# ============================================================

def phi_n_sym(roots):
    """Symbolic Phi_n."""
    n = len(roots)
    total = 0
    for i in range(n):
        H_i = sum(1/(roots[i] - roots[j]) for j in range(n) if j != i)
        total += H_i**2
    return total

def disc_sym(roots):
    """Symbolic discriminant."""
    n = len(roots)
    d = sp.Integer(1)
    for i in range(n):
        for j in range(i+1, n):
            d *= (roots[i] - roots[j])**2
    return d

# ============================================================
# n = 2: exact
# ============================================================
print("=" * 70)
print("n = 2")
print("=" * 70)

a = sp.Symbol('a', positive=True)
roots_2 = [a, -a]
phi_2 = phi_n_sym(roots_2)
phi_2 = sp.simplify(phi_2)
print(f"Phi_2 = {phi_2}")
print(f"1/Phi_2 = {sp.simplify(1/phi_2)}")
# kappa_2 = 2*a^2 for n=2 centered
# 1/Phi_2 = 2*a^2 = kappa_2
# C_2 = 1, R_2 = 0
print("C_2 = 1, R_2 = 0")

# ============================================================
# n = 3: exact in cumulants
# ============================================================
print("\n" + "=" * 70)
print("n = 3")
print("=" * 70)

# For n=3 centered, we already know from the verify script:
# Phi_3 = (9/2)*k2^2 / (k2^3 - k3^2/3)
# 1/Phi_3 = (2/9)*(k2^3 - k3^2/3)/k2^2
#          = (2/9)*k2 - (2/27)*k3^2/k2^2

k2, k3, k4, k5, k6 = symbols('k2 k3 k4 k5 k6')

inv_phi_3 = Rational(2,9) * k2 - Rational(2,27) * k3**2 / k2**2
print(f"1/Phi_3 = {inv_phi_3}")
print(f"C_3 = 2/9 = {Rational(2,9)}")
print(f"R_3 = {-Rational(2,27) * k3**2 / k2**2}")
print(f"2/(3*2) = {Rational(2,6)} = {Rational(1,3)} != 2/9")
print(f"Actually 2/(n*(n-1)) = 2/6 = 1/3 for n=3, but C_3 = 2/9")
print(f"So the conjecture C_n = 2/(n*(n-1)) is WRONG even for n=3!")

# Wait, let me re-check. The original claim was C_3 = 2/9.
# And 2/(3*2) = 1/3, which is NOT 2/9.
# So C_3 = 2/9 ≠ 1/3 = 2/(n*(n-1)).
# Something is off with the claimed conjecture.

# Let me re-examine: maybe the conjecture is C_n = 2/(n^2) ?
# For n=2: 2/4 = 1/2 ≠ 1. No.
# For n=3: 2/9 ✓

# Or maybe C_n depends on the normalization of cumulants.
# Let me go back to first principles.

# ============================================================
# n = 3: Derive from scratch
# ============================================================
print("\n--- n=3: Derivation from scratch ---")

# Centered cubic x^3 + e2*x - e3 = 0 (Vieta: a1=0, a2=e2, a3=-e3)
# Roots r1, r2, r3 with r1+r2+r3=0, r1*r2+r1*r3+r2*r3=e2, r1*r2*r3=e3

e2s, e3s = symbols('e2 e3')

# Known: Phi_3 * disc = 18 * e2^2
# disc(n=3, centered) = -4*e2^3 - 27*e3^2
# So Phi_3 = 18*e2^2 / (-4*e2^3 - 27*e3^2)

# 1/Phi_3 = (-4*e2^3 - 27*e3^2) / (18*e2^2)
#          = -2*e2/9 - 3*e3^2/(2*e2^2)
inv_phi3_e = -Rational(2,9)*e2s - Rational(3,2)*e3s**2/e2s**2
print(f"1/Phi_3 (in e2,e3) = {inv_phi3_e}")

# Now convert to cumulants. The finite free cumulants for n=3 centered:
# Normalized coefficients: tilde_a_2 = e2/C(3,2) = e2/3
#                          tilde_a_3 = e3/C(3,3) = e3  (careful: tilde_a_3 = (-1)^3*a3/C(3,3))
# a3 = -e3, so tilde_a_3 = (-1)^3*(-e3)/1 = e3. OK.

# Cumulant formulas for n=3 centered (from Arizmendi-Perales):
# kappa_2 = -n*tilde_a_2 = -3*(e2/3) = -e2
# kappa_3 = (n^2/2)*tilde_a_3 = (9/2)*e3

# So: e2 = -k2, e3 = 2*k3/9
# 1/Phi_3 = -2/9*(-k2) - 3/2*(2*k3/9)^2/(-k2)^2
#          = 2*k2/9 - 3/2 * 4*k3^2/81 / k2^2
#          = 2*k2/9 - (6*k3^2)/(81*k2^2)
#          = 2*k2/9 - 2*k3^2/(27*k2^2)
print(f"\nConverted: 1/Phi_3 = (2/9)*k2 - (2/27)*k3^2/k2^2  ✓")

# So C_3 = 2/9, R_3 = -(2/27)*k3^2/k2^2.
# Note C_3 = 2/9 = 2/n^2 for n=3.
# For n=2: C_2 = 1 = 2/2 = 2/(n*(n-1)) or 1/1 = ...

# Hmm, 2/n^2: n=2 gives 1/2 ≠ 1. So C_n = 2/n^2 doesn't work either.
# Let me just compute C_4 exactly and look for the pattern.

# ============================================================
# n = 4: Exact symbolic computation
# ============================================================
print("\n" + "=" * 70)
print("n = 4: Exact computation")
print("=" * 70)

# From Part 2, we found: Phi_4 * disc = -8*e2^5 - 64*e2^3*e4 - 36*e2^2*e3^2 + 384*e2*e4^2 - 432*e3^2*e4
# Let's verify this symbolically and then express in cumulants.

# First verify the numerical coefficients by checking a specific case
print("Verifying Phi_4*disc formula...")

# Test with specific roots
test_roots = np.array([3.0, 1.0, -0.5, -3.5])
assert abs(sum(test_roots)) < 1e-10

e2_t = sum(test_roots[i]*test_roots[j] for i in range(4) for j in range(i+1,4))
e3_t = sum(test_roots[i]*test_roots[j]*test_roots[k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4))
e4_t = test_roots[0]*test_roots[1]*test_roots[2]*test_roots[3]

disc_t = 1.0
for i in range(4):
    for j in range(i+1, 4):
        disc_t *= (test_roots[i] - test_roots[j])**2

def phi_n_num(roots):
    n = len(roots)
    total = 0.0
    for i in range(n):
        H = sum(1.0/(roots[i]-roots[j]) for j in range(n) if j != i)
        total += H**2
    return total

phi_t = phi_n_num(test_roots)
phi_disc_actual = phi_t * disc_t
phi_disc_formula = (-8*e2_t**5 - 64*e2_t**3*e4_t - 36*e2_t**2*e3_t**2
                    + 384*e2_t*e4_t**2 - 432*e3_t**2*e4_t)
print(f"  Actual: {phi_disc_actual:.6f}")
print(f"  Formula: {phi_disc_formula:.6f}")
print(f"  Match: {abs(phi_disc_actual - phi_disc_formula) < 1e-4}")

# Good! Now express in cumulants.
# For n=4 centered:
# tilde_a_2 = a2/C(4,2) = e2/6
# tilde_a_3 = (-1)^3*a3/C(4,3) = -(-e3)/4 = e3/4
# tilde_a_4 = (-1)^4*a4/C(4,4) = e4/1 = e4

# Cumulant formulas (n=4 centered):
# kappa_2 = -4*tilde_a_2 = -4*e2/6 = -2*e2/3
# kappa_3 = (16/2)*tilde_a_3 = 8*e3/4 = 2*e3

# For kappa_4, need the centered NC(4) formula:
# tilde_a_4 = kappa_4/(n^3) + (# of NC pair partitions of {1,2,3,4})*kappa_2^2/n^2
# NC pair partitions of {1,2,3,4}: {12}{34}, {14}{23} = 2
# Wait, actually I need to think about this more carefully.
# The moment-cumulant relation for the normalized coefficients:
# For FINITE free convolution (Marcus, Spielman, Srivastava), the relation is:
#
# p(x) = sum_{k=0}^n (-1)^k C(n,k) a_k x^{n-k}
# where a_k = sum_{pi in NC(k)} prod_{V in pi} kappa_{|V|} * [n]_{k-|pi|} / [n]_k
# Wait, this doesn't look right.

# Let me look up the correct formula. The FINITE free cumulants (Arizmendi-Perales 2018)
# are defined via the relation:
# a_k = sum_{pi in NC(k)} prod_{V in pi} kappa_{|V|}
# where a_k are the BOOLEAN cumulants of the normalized measure mu_p.
# No, that's also not right.

# OK, let me just compute kappa_4 NUMERICALLY by finding the right formula.
# I'll use the approach: generate many 4-tuples, and find the formula that makes
# kappa_4 additive under MSS convolution.

print("\n--- Finding correct cumulant formulas for n=4 ---")

# The key constraint: kappa_k must be ADDITIVE under MSS convolution.
# kappa_k(p ⊞ q) = kappa_k(p) + kappa_k(q)

# For n=4 centered, express things in terms of ta_2, ta_3, ta_4:
# We know kappa_2 = -4*ta_2 (additive ✓)
# We know kappa_3 = 8*ta_3 (need to verify additivity)

# For kappa_4, must find a formula f(ta_2, ta_3, ta_4) that's additive.
# The general form: kappa_4 = alpha*ta_4 + beta*ta_2^2 + gamma*ta_3*?? + ...

# Under MSS convolution: ta_k(p⊞q) = ta_k(p) + ta_k(q)
# (This is the KEY property of normalized coefficients!)
# So kappa_k = linear function of ta_k, ta_i*ta_j, etc.?
# NO. ta_k is additive, so ANY function of ta_2,...,ta_n is just
# kappa_k = c_k * ta_k (if we want additivity and it depends only on ta_k).
# But kappa_4 might be a LINEAR function of ta_2, ta_3, ta_4:
# kappa_4 = alpha*ta_4 + beta*ta_2^2 + gamma*ta_2 + ... NO, this isn't additive if nonlinear.

# Wait: ta_k IS additive. So kappa_k = sum of linear and nonlinear terms.
# For ADDITIVITY: kappa_k(p⊞q) = kappa_k(p) + kappa_k(q)
# If kappa_4 = alpha*ta_4 + beta*ta_2^2, then
# kappa_4(p⊞q) = alpha*(ta_4(p)+ta_4(q)) + beta*(ta_2(p)+ta_2(q))^2
#               = alpha*ta_4(p) + alpha*ta_4(q) + beta*ta_2(p)^2 + beta*ta_2(q)^2 + 2*beta*ta_2(p)*ta_2(q)
#               ≠ kappa_4(p) + kappa_4(q) unless beta=0.

# So for additivity, kappa_k must be a LINEAR function of the ta_j:
# kappa_k = sum_j c_{k,j} * ta_j
# This means kappa_k = c_k * ta_k (+ possibly lower order ta_j).

# Actually, looking at the Arizmendi-Perales paper more carefully:
# The finite free cumulants are defined so that the NORMALIZED COEFFICIENTS
# themselves ARE the free cumulants (up to scaling)!
# Specifically: a_k = kappa_k for the finite free cumulants.
# No wait, that can't be right since we know kappa_2 = -n*ta_2.

# Let me just check: ARE the ta_k themselves additive?
print("Testing additivity of ta_k under MSS...")

from math import comb as comb_

def compute_ta(roots):
    n = len(roots)
    coeffs = np.poly(roots)  # [1, a1, ..., an]
    ta = [0.0] * (n+1)
    ta[0] = 1.0
    for k in range(1, n+1):
        ta[k] = ((-1)**k * coeffs[k]) / comb_(n, k)
    return ta

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

n = 4
for trial in range(5):
    rp = np.random.randn(n) * 2
    rp = rp - np.mean(rp)
    rq = np.random.randn(n) * 2
    rq = rq - np.mean(rq)

    rr = mss_convolve(rp, rq)

    ta_p = compute_ta(rp)
    ta_q = compute_ta(rq)
    ta_r = compute_ta(rr)

    for k in range(2, n+1):
        diff = ta_r[k] - (ta_p[k] + ta_q[k])
        print(f"  Trial {trial}, ta_{k}: p={ta_p[k]:.6f}, q={ta_q[k]:.6f}, "
              f"r={ta_r[k]:.6f}, p+q={ta_p[k]+ta_q[k]:.6f}, diff={diff:.2e}")

# ============================================================
# If ta_k are additive, then we can directly express 1/Phi_n in ta_k
# and the "cumulants" ARE the ta_k (up to scaling).
# ============================================================

print("\n" + "=" * 70)
print("KEY INSIGHT: If ta_k are additive, work directly with ta_k")
print("=" * 70)

# For n=4 centered: e2 = 6*ta_2, e3 = 4*ta_3, e4 = ta_4
# (from: ta_2 = e2/6, ta_3 = e3/4, ta_4 = e4)

# Phi_4*disc = -8*e2^5 - 64*e2^3*e4 - 36*e2^2*e3^2 + 384*e2*e4^2 - 432*e3^2*e4
# disc = 256*e4^3 - 128*e2^2*e4^2 + 144*e2*e3^2*e4 - 27*e3^4 + 16*e2^4*e4 - 4*e2^3*e3^2

# In terms of ta_k: e2 = 6*ta_2, e3 = 4*ta_3, e4 = ta_4
# Let t2 = ta_2, t3 = ta_3, t4 = ta_4

t2, t3, t4, t5, t6 = symbols('t2 t3 t4 t5 t6')

e2_4 = 6*t2
e3_4 = 4*t3
e4_4 = t4

N4 = -8*e2_4**5 - 64*e2_4**3*e4_4 - 36*e2_4**2*e3_4**2 + 384*e2_4*e4_4**2 - 432*e3_4**2*e4_4
N4 = expand(N4)
print(f"\nNumerator Phi_4*disc in ta_k:")
print(f"  = {N4}")

D4 = (256*e4_4**3 - 128*e2_4**2*e4_4**2 + 144*e2_4*e3_4**2*e4_4
      - 27*e3_4**4 + 16*e2_4**4*e4_4 - 4*e2_4**3*e3_4**2)
D4 = expand(D4)
print(f"\nDisc in ta_k:")
print(f"  = {D4}")

# 1/Phi_4 = disc / (Phi_4*disc) = D4 / N4
inv_phi_4 = D4 / N4
print(f"\n1/Phi_4 = disc / (Phi_4*disc)")
print("  (rational function in t2, t3, t4)")

# Simplify
inv_phi_4_simplified = cancel(inv_phi_4)
num_4, den_4 = sp.fraction(inv_phi_4_simplified)
num_4 = expand(num_4)
den_4 = expand(den_4)
print(f"\n1/Phi_4 simplified:")
print(f"  Numerator: {num_4}")
print(f"  Denominator: {den_4}")

# Factor the numerator and denominator
print(f"\n  Numerator factored: {factor(num_4)}")
print(f"  Denominator factored: {factor(den_4)}")

# Now compute the expansion around t3=t4=0:
# 1/Phi_4 = f(t2) + O(t3, t4) terms
# At t3=t4=0: 1/Phi_4 = D4(t2,0,0) / N4(t2,0,0)

D4_0 = D4.subs([(t3, 0), (t4, 0)])
N4_0 = N4.subs([(t3, 0), (t4, 0)])
print(f"\nAt t3=t4=0:")
print(f"  Disc = {D4_0}")
print(f"  Phi*disc = {N4_0}")
inv_phi_4_0 = cancel(D4_0 / N4_0)
print(f"  1/Phi_4 = {inv_phi_4_0}")
# This should give us C_4_ta * t2 (the "C_4" in terms of ta_2).

# ============================================================
# n = 3 in ta_k for comparison
# ============================================================
print("\n--- n=3 in ta_k for comparison ---")
# ta_2 = e2/3, ta_3 = e3
# e2 = 3*t2, e3 = t3
# 1/Phi_3 = (-4*e2^3 - 27*e3^2)/(18*e2^2) = (-4*(3t2)^3 - 27*t3^2)/(18*(3t2)^2)
# = (-108*t2^3 - 27*t3^2)/(162*t2^2)
# = (-108*t2^3)/(162*t2^2) - 27*t3^2/(162*t2^2)
# = -2*t2/3 - t3^2/(6*t2^2)
inv_phi3_ta = (-108*t2**3 - 27*t3**2) / (162*t2**2)
inv_phi3_ta = cancel(inv_phi3_ta)
print(f"1/Phi_3 = {inv_phi3_ta}")
print(f"        = {expand(inv_phi3_ta)}")
# Hmm, this gives -2*t2/3 - t3^2/(6*t2^2)
# But we know 1/Phi_3 should be POSITIVE for real-rooted polys.
# For centered cubic with real roots: e2 < 0 (since e2 = sum_{i<j} r_i*r_j
# and with r1+r2+r3=0, this is -(r1^2+r2^2+r3^2)/2 < 0)
# So ta_2 = e2/3 < 0 as well, hence t2 < 0.
# -2*t2/3 > 0 when t2 < 0. ✓
# -t3^2/(6*t2^2) < 0 always. But the total should be > 0.

# At t3=0: 1/Phi_3 = -2*t2/3. Since t2 < 0, this is > 0. ✓
# So "C_3_ta" = -2/3 (as a coefficient of t2, where t2 < 0 for valid polys)

print(f"\nC_3_ta = -2/3 (coefficient of t2 in 1/Phi_3)")
print(f"C_2_ta = -1 (coefficient of t2 in 1/Phi_2)")

# For n=2: ta_2 = e2/C(2,2) = e2 = a*(-a) = -a^2
# 1/Phi_2 = 2*a^2 = -2*ta_2 ... wait
# kappa_2 = -2*ta_2 = 2*a^2, 1/Phi_2 = kappa_2 = -2*ta_2
# So in ta_2: 1/Phi_2 = -2*ta_2, C_2_ta = -2

# Hmm, let me re-check:
# n=2: ta_2 = e2/C(2,2) = e2/1 = e2, and e2 = a*(-a) = -a^2
# 1/Phi_2 = 2*a^2 = -2*e2 = -2*ta_2
# Wait, kappa_2 = -n*ta_2 = -2*ta_2 = -2*(-a^2) = 2*a^2. ✓
# 1/Phi_2 = kappa_2 = -2*ta_2

# So C_2_ta = -2, meaning 1/Phi_2 = -2*ta_2.
# For n=3: 1/Phi_3 = (-2/3)*ta_2 - (1/6)*ta_3^2/ta_2^2

# Pattern so far: C_n_ta = -2/n ?
# n=2: -2/2 = -1 ≠ -2. No.
# n=2: C_2_ta = -2, n=3: C_3_ta = -2/3
# Ratio: 3*C_3_ta / C_2_ta = 3*(-2/3)/(-2) = 1. Hmm.
# Maybe C_n_ta = -2/??? Let me just compute for n=4.

# ============================================================
# n=4: Leading term in ta expansion
# ============================================================
print("\n" + "=" * 70)
print("n = 4: Computing C_4_ta (coefficient of t2 at t3=t4=0)")
print("=" * 70)

# From above:
# Phi_4*disc at t3=t4=0:
N4_t3t4_0 = N4.subs([(t3, 0), (t4, 0)])
D4_t3t4_0 = D4.subs([(t3, 0), (t4, 0)])
print(f"  N4(t3=t4=0) = {N4_t3t4_0}")
print(f"  D4(t3=t4=0) = {D4_t3t4_0}")

# If both are 0, we have 0/0 and need L'Hopital or series expansion.
# Let's check: at t3=t4=0, e2=6*t2, e3=0, e4=0
# disc = -4*e2^3*e3^2 - 27*e3^4 + 16*e2^4*e4 + ... = 0 when e3=e4=0 (all terms have e3 or e4)
# Wait, disc = 256*e4^3 - 128*e2^2*e4^2 + 144*e2*e3^2*e4 - 27*e3^4 + 16*e2^4*e4 - 4*e2^3*e3^2
# At e3=e4=0: all terms are 0. So disc = 0.
# Similarly, N4 = Phi_4*disc at e3=e4=0 is also 0.
# This makes sense: at e3=e4=0, the polynomial is x^4+e2*x^2 = x^2(x^2+e2),
# which has REPEATED ROOTS (x=0 is double). So disc=0 and Phi is undefined.

# This means we CANNOT simply set t3=t4=0. We need a different approach.

# Instead, let's parametrize around "almost repeated root" case differently.
# Or better: express 1/Phi_4 as a rational function of t2,t3,t4 and analyze its structure.

print("\nCannot set t3=t4=0 (gives 0/0). Need full rational expression.")

# Let's compute 1/Phi_4 = D4/N4 and analyze.
# Factor out common factors from N4 and D4.
from sympy import gcd

g = sp.gcd(N4, D4)
print(f"\nGCD of numerator and denominator: {factor(g)}")

N4_red = cancel(N4/g)
D4_red = cancel(D4/g)
print(f"Reduced numerator: {expand(N4_red)}")
print(f"Reduced denominator: {expand(D4_red)}")

# ============================================================
# n = 4: Try a different approach - series in t4
# ============================================================
print("\n--- Series expansion in t4 around t4=0 ---")

# 1/Phi_4 = D4/N4 as functions of t2, t3, t4
# Let's compute this as a power series in t4 with coefficients in t2, t3

# First, separate by powers of t4
N4_poly = Poly(N4, t4)
D4_poly = Poly(D4, t4)

print(f"N4 as poly in t4:")
for power, coeff in N4_poly.as_dict().items():
    print(f"  t4^{power[0]}: {coeff}")

print(f"\nD4 as poly in t4:")
for power, coeff in D4_poly.as_dict().items():
    print(f"  t4^{power[0]}: {coeff}")

# ============================================================
# n = 4: More direct approach
# ============================================================
print("\n" + "=" * 70)
print("n = 4: Direct approach via Phi as function of e-sym")
print("=" * 70)

# Since t3=t4=0 gives repeated roots, let's think about what "C_n" means.
# The claim is: 1/Phi_n = C_n * kappa_2 + R_n
# where kappa_2 is additive.
# But if 1/Phi_n is a rational function of cumulants (not polynomial),
# then C_n is the "linear part" in some appropriate sense.

# For n=3: 1/Phi_3 = (2/9)*k2 + R_3 where R_3 = -(2/27)*k3^2/k2^2
# Here C_3 = 2/9 and R_3 depends on k2 and k3.

# For n=4: we need to express 1/Phi_4 similarly.
# But first we need the CORRECT cumulant definition.
# Let me just work with the ta_k (normalized coefficients) which ARE additive.

# 1/Phi_4 = D4(t2,t3,t4) / N4(t2,t3,t4)
# We want the "linear in t2" contribution when t3,t4 are "small".

# Since the rational function blows up/degenerates at t3=t4=0,
# let's instead set t4=0 and look at 1/Phi_4 as function of t2, t3.

N4_t4_0 = expand(N4.subs(t4, 0))
D4_t4_0 = expand(D4.subs(t4, 0))
print(f"At t4=0:")
print(f"  N4 = {N4_t4_0}")
print(f"  D4 = {D4_t4_0}")

inv_phi4_t4_0 = cancel(D4_t4_0 / N4_t4_0)
print(f"  1/Phi_4 = {inv_phi4_t4_0}")

# Factor
num_t4_0, den_t4_0 = sp.fraction(inv_phi4_t4_0)
print(f"  Num: {factor(num_t4_0)}")
print(f"  Den: {factor(den_t4_0)}")

# ============================================================
# A DIFFERENT APPROACH: Express 1/Phi directly via power sums
# ============================================================
print("\n" + "=" * 70)
print("DIFFERENT APPROACH: 1/Phi_n via power sums p_k = sum lambda_i^k")
print("=" * 70)

# H_i = sum_{j≠i} 1/(lambda_i - lambda_j)
# This is related to the logarithmic derivative of the polynomial.
# Actually, if p(x) = prod_i (x - lambda_i), then p'(x) = sum_i prod_{j≠i}(x-lambda_j)
# So p'(lambda_i) = prod_{j≠i}(lambda_i - lambda_j)
# And H_i = p''(lambda_i) / (2*p'(lambda_i))

# Phi_n = sum_i [p''(lambda_i)/(2*p'(lambda_i))]^2
# = (1/4) * sum_i [p''(lambda_i)/p'(lambda_i)]^2

# This can also be written using:
# H_i = d/dx [log p'(x)] evaluated at x = lambda_i
# No, d/dx [log p(x)] = p'(x)/p(x) = sum_j 1/(x - lambda_j)
# So H_i = sum_{j≠i} 1/(lambda_i - lambda_j) = lim_{x->lambda_i} [sum_j 1/(x-lambda_j) - 1/(x-lambda_i)]
# = lim_{x->lambda_i} [p'(x)/p(x) - 1/(x-lambda_i)]
# This is (p'(lambda_i)/p(lambda_i) - ...) which involves the "regularized" value.
# Actually p'(lambda_i) = prod_{j≠i}(lambda_i-lambda_j), so:
# sum_j 1/(lambda_i - lambda_j) = p''(lambda_i)/(2*prod_{j≠i}(lambda_i-lambda_j))

# This is getting complex. Let me just compute for small n numerically
# to determine the EXACT rational functions.

# For n=4, let's use higher precision numerics with mpmath.
print("\nUsing mpmath for high-precision n=4 computation...")

import mpmath
mpmath.mp.dps = 50  # 50 decimal places

def phi_n_mp(roots):
    n = len(roots)
    total = mpmath.mpf(0)
    for i in range(n):
        H = sum(mpmath.mpf(1)/(roots[i]-roots[j]) for j in range(n) if j != i)
        total += H**2
    return total

def disc_mp(roots):
    n = len(roots)
    d = mpmath.mpf(1)
    for i in range(n):
        for j in range(i+1, n):
            d *= (roots[i]-roots[j])**2
    return d

# Generate a specific nice n=4 example
# Roots: 3, 1, -1, -3 (sum = 0)
r_nice = [mpmath.mpf(3), mpmath.mpf(1), mpmath.mpf(-1), mpmath.mpf(-3)]
e2_nice = sum(r_nice[i]*r_nice[j] for i in range(4) for j in range(i+1,4))
e3_nice = sum(r_nice[i]*r_nice[j]*r_nice[k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4))
e4_nice = r_nice[0]*r_nice[1]*r_nice[2]*r_nice[3]
print(f"Nice roots: 3, 1, -1, -3")
print(f"  e2 = {e2_nice}, e3 = {e3_nice}, e4 = {e4_nice}")

phi_nice = phi_n_mp(r_nice)
disc_nice = disc_mp(r_nice)
print(f"  Phi_4 = {phi_nice}")
print(f"  disc = {disc_nice}")
print(f"  1/Phi_4 = {1/phi_nice}")

ta2_nice = e2_nice / 6
ta3_nice = e3_nice / 4
ta4_nice = e4_nice
print(f"  ta_2 = {ta2_nice}, ta_3 = {ta3_nice}, ta_4 = {ta4_nice}")

# For this example: roots are {-3,-1,1,3}, symmetric
# e2 = (-3)(-1)+(-3)(1)+(-3)(3)+(-1)(1)+(-1)(3)+(1)(3)
#     = 3 + (-3) + (-9) + (-1) + (-3) + 3 = -10
# e3 = (-3)(-1)(1)+(-3)(-1)(3)+(-3)(1)(3)+(-1)(1)(3) = 3+9-9-3 = 0
# e4 = (-3)(-1)(1)(3) = 9
# So ta_2 = -10/6 = -5/3, ta_3 = 0, ta_4 = 9

# Another example: roots -2, -1, 1, 2
r_nice2 = [mpmath.mpf(-2), mpmath.mpf(-1), mpmath.mpf(1), mpmath.mpf(2)]
e2_n2 = sum(r_nice2[i]*r_nice2[j] for i in range(4) for j in range(i+1,4))
e3_n2 = sum(r_nice2[i]*r_nice2[j]*r_nice2[k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4))
e4_n2 = r_nice2[0]*r_nice2[1]*r_nice2[2]*r_nice2[3]
print(f"\nNice roots 2: -2, -1, 1, 2")
print(f"  e2 = {e2_n2}, e3 = {e3_n2}, e4 = {e4_n2}")
phi_n2 = phi_n_mp(r_nice2)
print(f"  1/Phi_4 = {1/phi_n2}")

# Verify formula: Phi_4*disc = -8e2^5 - 64e2^3*e4 - 36e2^2*e3^2 + 384e2*e4^2 - 432e3^2*e4
pd_n2 = -8*e2_n2**5 - 64*e2_n2**3*e4_n2 - 36*e2_n2**2*e3_n2**2 + 384*e2_n2*e4_n2**2 - 432*e3_n2**2*e4_n2
print(f"  Phi_4*disc (formula) = {pd_n2}")
print(f"  Phi_4*disc (actual) = {phi_n2 * disc_mp(r_nice2)}")

# ============================================================
# Express 1/Phi_4 symbolically using the known N_4 = Phi_4*disc
# ============================================================
print("\n" + "=" * 70)
print("SYMBOLIC 1/Phi_4")
print("=" * 70)

# 1/Phi_4 = disc / N_4 where N_4 = Phi_4 * disc

# In e-sym:
e2_s, e3_s, e4_s = symbols('e2 e3 e4')
N4_e = -8*e2_s**5 - 64*e2_s**3*e4_s - 36*e2_s**2*e3_s**2 + 384*e2_s*e4_s**2 - 432*e3_s**2*e4_s
D4_e = 256*e4_s**3 - 128*e2_s**2*e4_s**2 + 144*e2_s*e3_s**2*e4_s - 27*e3_s**4 + 16*e2_s**4*e4_s - 4*e2_s**3*e3_s**2

inv_phi4_e = cancel(D4_e / N4_e)
print(f"1/Phi_4 = {inv_phi4_e}")

num4_e, den4_e = sp.fraction(inv_phi4_e)
print(f"\nNumerator: {factor(num4_e)}")
print(f"Denominator: {factor(den4_e)}")

# Let's try partial fraction / collect by e3 powers
# First, factor out what we can
print(f"\nN4 factored: {factor(N4_e)}")
print(f"D4 factored: {factor(D4_e)}")

# ============================================================
# Verify: For n=3, N_3 = Phi_3 * disc. What is N_3?
# ============================================================
print("\n--- Verifying N_3 ---")
# For n=3 centered: N_3 = 18*e2^2
# disc_3 = -4*e2^3 - 27*e3^2
# 1/Phi_3 = disc/N_3 = (-4*e2^3 - 27*e3^2)/(18*e2^2)

# Now for n=4: N_4 = -8*e2^5 - 64*e2^3*e4 - 36*e2^2*e3^2 + 384*e2*e4^2 - 432*e3^2*e4
# Try to factor: common factor of 4?
print(f"  N_4/4 = {expand(N4_e/4)}")
# -2*e2^5 - 16*e2^3*e4 - 9*e2^2*e3^2 + 96*e2*e4^2 - 108*e3^2*e4
# Hmm, not obviously nice.

# Let me try: is N_4 = c * e2^2 * something?
# -8*e2^5: has e2^5
# -36*e2^2*e3^2: has e2^2
# -64*e2^3*e4: has e2^3
# 384*e2*e4^2: has e2
# -432*e3^2*e4: has e3^2, no e2
# So no common e2 factor.

# ============================================================
# For n=5: compute N_5 = Phi_5 * disc numerically
# ============================================================
print("\n" + "=" * 70)
print("n = 5: Computing N_5 = Phi_5 * disc")
print("=" * 70)

# Weight analysis for n=5:
# disc has weight n(n-1) = 20
# Phi has weight -2
# N_5 = Phi_5 * disc has weight 18
# Monomials: e2^a * e3^b * e4^c * e5^d with 2a+3b+4c+5d = 18
# Enumerate:
mono_5 = []
for a in range(10):
    for b in range(7):
        for c in range(5):
            for d in range(4):
                if 2*a + 3*b + 4*c + 5*d == 18:
                    mono_5.append((a, b, c, d))

print(f"Weight-18 monomials in e2,e3,e4,e5: {len(mono_5)}")
for m in sorted(mono_5):
    print(f"  e2^{m[0]} * e3^{m[1]} * e4^{m[2]} * e5^{m[3]}")

# Generate samples and fit
def compute_n5_data(roots):
    n = 5
    assert len(roots) == n and abs(sum(roots)) < 1e-8
    e2 = sum(roots[i]*roots[j] for i in range(n) for j in range(i+1,n))
    e3 = sum(roots[i]*roots[j]*roots[k] for i in range(n) for j in range(i+1,n) for k in range(j+1,n))
    e4 = sum(roots[i]*roots[j]*roots[k]*roots[l] for i in range(n) for j in range(i+1,n) for k in range(j+1,n) for l in range(k+1,n))
    e5 = np.prod(roots)

    disc_val = 1.0
    for i in range(n):
        for j in range(i+1, n):
            disc_val *= (roots[i] - roots[j])**2

    phi_val = phi_n_num(roots)
    N5 = phi_val * disc_val

    return e2, e3, e4, e5, N5

data_5 = []
for _ in range(3000):
    r = np.random.randn(5) * 2
    r = r - np.mean(r)
    r = np.sort(r)
    if min(np.diff(r)) < 0.3:
        continue
    e2, e3, e4, e5, N5 = compute_n5_data(r)
    data_5.append((e2, e3, e4, e5, N5))

print(f"Generated {len(data_5)} valid n=5 samples")

# Build design matrix
X5 = np.array([[e2**m[0] * e3**m[1] * e4**m[2] * e5**m[3]
                for m in sorted(mono_5)]
               for e2, e3, e4, e5, _ in data_5])
y5 = np.array([N5 for _, _, _, _, N5 in data_5])

coeffs_5, _, rank_5, _ = np.linalg.lstsq(X5, y5, rcond=None)
print(f"Rank of design matrix: {rank_5}")

print("\nN_5 = Phi_5 * disc coefficients:")
for i, m in enumerate(sorted(mono_5)):
    c = coeffs_5[i]
    if abs(c) > 0.01:
        frac = Fraction(c).limit_denominator(10000)
        print(f"  e2^{m[0]}*e3^{m[1]}*e4^{m[2]}*e5^{m[3]}: {c:.6f} ≈ {frac}")

# Verify residual
resid_5 = np.sqrt(np.mean((X5 @ coeffs_5 - y5)**2)) / np.sqrt(np.mean(y5**2))
print(f"\nRelative residual: {resid_5:.2e}")

print("\n" + "=" * 70)
print("SUMMARY OF Phi_n * disc = N_n FORMULAS")
print("=" * 70)
print("n=2: N_2 = 2  (Phi_2 = 1/(2*a^2), disc = 4*a^2, N_2 = 2)")
print("n=3: N_3 = 18*e2^2")
print("n=4: N_4 = -8*e2^5 - 64*e2^3*e4 - 36*e2^2*e3^2 + 384*e2*e4^2 - 432*e3^2*e4")

# ============================================================
# KEY QUESTION: What is C_n in the decomposition?
# ============================================================
print("\n" + "=" * 70)
print("KEY: Re-examining the C_n decomposition")
print("=" * 70)

# The original conjecture says: 1/Phi_n(p ⊞ q) >= 1/Phi_n(p) + 1/Phi_n(q)
# and decomposes 1/Phi_n = C_n * kappa_2 + R_n
#
# For n=2: 1/Phi_2 = kappa_2 (C_2=1, R_2=0)
# For n=3: 1/Phi_3 = (2/9)*k2 - (2/27)*k3^2/k2^2
# where k2, k3 are the Arizmendi-Perales finite free cumulants.
#
# In terms of ta_k (which ARE additive under MSS):
# n=2: 1/Phi_2 = -2*ta_2 = kappa_2 = -n*ta_2
# n=3: 1/Phi_3 = (-2/3)*ta_2 - (1/6)*ta_3^2/ta_2^2
#
# So C_3_ta = -2/3 = -n*C_3_kappa = -3*(2/9) = -2/3 ✓
# And C_2_ta = -2 = -n*C_2_kappa = -2*1 = -2 ✓
#
# Wait: kappa_2 = -n*ta_2, so C_n_kappa = C_n_ta / (-n) = (-2/n)/(-n) = 2/n^2.
# n=2: 2/4 = 1/2 ≠ 1. WRONG!
# Let me re-check: 1/Phi_2 = -2*ta_2 = (-2)*ta_2
# 1/Phi_2 = C_2_kappa * kappa_2 = C_2_kappa * (-2*ta_2)
# So C_2_kappa = (-2*ta_2)/(-2*ta_2) = 1. ✓
#
# 1/Phi_3 = (-2/3)*ta_2 - (1/6)*ta_3^2/ta_2^2
# = (-2/3)*(-k2/3) - (1/6)*(2*k3/9)^2/(-k2/3)^2
# = (2/9)*k2 - (1/6)*(4*k3^2/81)/(k2^2/9)
# = (2/9)*k2 - (1/6)*(4*9*k3^2)/(81*k2^2)
# = (2/9)*k2 - (1/6)*(36*k3^2)/(81*k2^2)
# = (2/9)*k2 - (1/6)*(4*k3^2)/(9*k2^2)
# = (2/9)*k2 - (2*k3^2)/(27*k2^2) ✓

# So the "C_n" in the kappa formulation is C_n = 2/(n^2-n)? Let me check:
# n=2: 2/(4-2) = 1 ✓
# n=3: 2/(9-3) = 2/6 = 1/3 ≠ 2/9.
# So NO, that's not right either.

# Actually for n=3: C_3 = 2/9 = 2/n^2 for n=3 ✓ but 2/n^2 for n=2 gives 1/2 ≠ 1.

# Hmm. The pattern is C_2=1, C_3=2/9. Let me compute C_4 properly.
# For n=4, 1/Phi_4 = disc/(Phi_4*disc) is a rational function.
# Let me compute it in kappa variables.

# For n=4 centered, the cumulants are:
# kappa_2 = -4*ta_2 = -4*e2/6 = -2e2/3
# kappa_3 = 8*ta_3 = 8*e3/4 = 2*e3
# For kappa_4, I need to be careful.

# Let me determine kappa_4 by checking additivity.
# Actually, the Arizmendi-Perales cumulants for degree n polynomial have:
# kappa_k = (-1)^{k-1} * n! / (n-k)! * ta_k + (nonlinear terms in ta_j for j<k)
# For the centered case (ta_1=0), the lower-order terms involve ONLY ta_j with j>=2.

# Let me determine kappa_4 by additivity testing.
print("\n--- Determining kappa_4 formula via additivity ---")

# The correct approach: finite free cumulants are defined via
# R-transform or moment-cumulant formula on the partition lattice.
# For the FINITE free case (not asymptotic), the relation is:
#
# For p(x) = sum_{k=0}^n (-1)^k C(n,k) a_k x^{n-k} with a_0=1:
# The R-transform satisfies: a_k = sum over NC(k) of c(n, pi) * prod kappa_{|V|}
# where c(n, pi) is a correction factor depending on n.
#
# From Arizmendi-Perales (2018), Theorem 3.4:
# a_k = (1/[n]_k) * sum_{pi in NC(k)} n^{|pi|} * mu_k(pi) * prod kappa_{|V|}
# where [n]_k = n*(n-1)*...*(n-k+1) is the falling factorial
# and mu_k is some weight.
#
# Actually, let me just use the SIMPLEST definition that works:
# The ta_k ARE additive. So let's work entirely with ta_k.

# FORGET about kappa_k. Just express 1/Phi_n in terms of the ADDITIVE quantities ta_k.

print("\n" + "=" * 70)
print("WORKING WITH ADDITIVE ta_k DIRECTLY")
print("=" * 70)

# For n=2 centered: ta_2 = e2. 1/Phi_2 = -2*ta_2.
# For n=3 centered: ta_2 = e2/3, ta_3 = e3.
#   1/Phi_3 = -(2/3)*ta_2 - (1/6)*ta_3^2/ta_2^2
# For n=4 centered: ta_2 = e2/6, ta_3 = e3/4, ta_4 = e4.
#   1/Phi_4 = disc_4(in ta) / N_4(in ta)

# Let's compute this properly for n=4.
# e2 = 6*t2, e3 = 4*t3, e4 = t4

# disc_4 = 256*t4^3 - 128*(6t2)^2*t4^2 + 144*(6t2)*(4t3)^2*t4 - 27*(4t3)^4 + 16*(6t2)^4*t4 - 4*(6t2)^3*(4t3)^2
D4_ta = (256*t4**3 - 128*36*t2**2*t4**2 + 144*6*16*t2*t3**2*t4
         - 27*256*t3**4 + 16*1296*t2**4*t4 - 4*216*16*t2**3*t3**2)
D4_ta = expand(D4_ta)

# N4(in ta) = -8*(6t2)^5 - 64*(6t2)^3*t4 - 36*(6t2)^2*(4t3)^2 + 384*(6t2)*t4^2 - 432*(4t3)^2*t4
N4_ta = (-8*7776*t2**5 - 64*216*t2**3*t4 - 36*36*16*t2**2*t3**2
         + 384*6*t2*t4**2 - 432*16*t3**2*t4)
N4_ta = expand(N4_ta)

print(f"D4(ta) = disc_4 in ta:")
print(f"  {D4_ta}")
print(f"\nN4(ta) = Phi_4*disc in ta:")
print(f"  {N4_ta}")

inv_phi4_ta = cancel(D4_ta / N4_ta)
num4_ta, den4_ta = sp.fraction(inv_phi4_ta)
num4_ta = expand(num4_ta)
den4_ta = expand(den4_ta)
print(f"\n1/Phi_4 in ta = num/den")
print(f"  Num: {num4_ta}")
print(f"  Den: {den4_ta}")
print(f"\n  Num factored: {factor(num4_ta)}")
print(f"  Den factored: {factor(den4_ta)}")

# ============================================================
# NUMERICAL VERIFICATION of the ta-based expressions
# ============================================================
print("\n" + "=" * 70)
print("NUMERICAL VERIFICATION")
print("=" * 70)

# Check that 1/Phi_4 matches the formula
for trial in range(5):
    roots = np.random.randn(4) * 2
    roots = roots - np.mean(roots)
    roots = np.sort(roots)
    if min(np.diff(roots)) < 0.3:
        continue

    e2_v = sum(roots[i]*roots[j] for i in range(4) for j in range(i+1,4))
    e3_v = sum(roots[i]*roots[j]*roots[k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4))
    e4_v = roots[0]*roots[1]*roots[2]*roots[3]

    t2_v = e2_v / 6
    t3_v = e3_v / 4
    t4_v = e4_v

    inv_phi_actual = 1.0 / phi_n_num(roots)

    # Evaluate the formula
    inv_phi_formula = float(inv_phi4_ta.subs([(t2, t2_v), (t3, t3_v), (t4, t4_v)]))

    print(f"  Trial: 1/Phi actual={inv_phi_actual:.10f}, formula={inv_phi_formula:.10f}, "
          f"ratio={inv_phi_actual/inv_phi_formula:.10f}")
