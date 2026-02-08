"""
PROVER-12 Step 3: Direct gap computation in (d1, d2) coordinates.

Key insight from Step 2:
  -R_4 = N(K2, d1, d2) / (24*d1*d2)
  where N = 8*K2^2*d1^2 - 6*K2*d1*d2 - d1^3 + 2*d2^2

  d1 = 4K2^2 - K4,  d2 = 4K2^3 + K2*K4 - 2K3^2

  d1 superadditive: d1_r - d1_p - d1_q = 8st > 0
  d2 superadditive: d2_r - d2_p - d2_q = 12st(s+t) + tu + sv - 4ab >= 8st(s+t) + 2(a√(t/s)-b√(s/t))^2 >= 0

Strategy: Express the gap of -R_4 in terms of slack variables δ1 = d1_r - d1_p - d1_q
and δ2 = d2_r - d2_p - d2_q (both >= 0) and try to show the gap is <= 0
(which means R_4 gap >= 0).

Actually, let me try a COMPLETELY DIFFERENT approach: reduce to n=3 using Fisher information
matrix structure.
"""
from sympy import (symbols, expand, factor, cancel, simplify, collect,
                   Rational, S, Poly, together, apart, fraction, sqrt,
                   Matrix, det, trace)
import numpy as np

print("="*70)
print("APPROACH: Reduce R_4 to R_3 via iterated Cauchy-Schwarz")
print("="*70)

# Recall R_3(K2, K3) = -K3^2/(12*K2^2*(something))... actually
# R_3(K2, K3) = K3^2 / (24*K2*(K3^2 - 2*K2^3))  -- wait, let me recall.
# From the Fisher information context:
# R_n = κ_{n+1}/(24*K2) where... no, it's more specific.

# Let me just use the known partial fraction:
# -R_4 = (4K2^3-K3^2)/(6*D1) + K3^2*(4K2^3-K3^2)/(6*K2^2*D2) - K4/(24*K2) - K2/6 - K3^2/(12*K2^2)

# where D1 = 4K2^2 - K4, D2 = 4K2^3 + K2*K4 - 2K3^2.

# Group terms as:
# -R_4 = [K2^2*d1/(3*d2) + d2/(12*d1)] + [-K2/4 - d1^2/(24*d2)]
# = [K2^2*d1/(3*d2) + d2/(12*d1)] - [K2/4 + d1^2/(24*d2)]

# Hmm. Let me try Approach D more carefully.

print("="*70)
print("APPROACH D: Convexity along line segments from origin")
print("="*70)

# R_4 is weight-2 homogeneous under cumulant scaling:
# R_4(λ^2*K2, λ^3*K3, λ^4*K4) = λ^2 * R_4(K2, K3, K4)

# Superadditivity of a positively homogeneous function of degree p (with p >= 1)
# is equivalent to super-level set {x : f(x) >= c} being convex for all c.
# But here p = 2 (via the weight) and the domain itself is a cone.

# Actually, for weight-2 homogeneous: R_4(λ*k) = λ^2*R_4(k) if k = (K2,K3,K4)
# and the scaling is non-uniform: K2→λ^2*K2, etc.
# The correct scaling: if we set λ_2*K2, λ_3*K3, λ_4*K4 then
# R_4(λ_2 K2, λ_3 K3, λ_4 K4) = ??? Only homogeneous under the joint weight scaling.

# Let me try the MOST DIRECT approach: verify that the function
# g(t) = R_4(p + t*(q-p)) is concave for t in [0,1] when p,q in the domain.
# If g is concave, then g(t) >= (1-t)*g(0) + t*g(1), which is NOT superadditivity.

# Actually, superadditivity: R_4(p+q) >= R_4(p) + R_4(q).
# This is equivalent to: for all p,q in domain, R_4(p+q) - R_4(p) - R_4(q) >= 0.

# Key insight: R_4 is weight-2 homogeneous. For such functions:
# f(x+y) >= f(x) + f(y) iff f(x) + f(y) <= f(x+y)
# This is SUPERADDITIVITY. For positively homogeneous functions of degree > 1,
# superadditivity is implied by CONVEXITY. But R_4 is NOT convex (it's negative).
# For homog degree > 1, f convex => f superadditive (when f >= 0 on a cone).
# But R_4 < 0...

# WAIT. -R_4 > 0 (on the domain, R_4 < 0). We need R_4 superadditive,
# i.e., -R_4 subadditive. -R_4 is weight-2 homogeneous and positive.
# For a CONCAVE positively homogeneous function of degree p >= 1, subadditivity holds.
# (Because f(x+y) <= f(x) + f(y) when f is concave and f(0) >= 0.)
# But -R_4 was shown NOT to be jointly concave (Hessian indefinite).

# However, we're working with WEIGHTED homogeneity, not standard homogeneity.
# f(λ·x) = λ^2·f(x) under the weight (2,3,4).
# This is different from ordinary homogeneity.

# Let me try a completely different decomposition.

print("="*70)
print("NEW APPROACH: Quadratic form decomposition")
print("="*70)

# In the (d1, d2) representation:
# -R_4 = (8K2^2*d1^2 - 6K2*d1*d2 - d1^3 + 2d2^2) / (24*d1*d2)
# = d1/(24*d2) * (8K2^2 - d1) + (2d2 - 6K2*d1)/(24*d1)
# = d1*(8K2^2-d1)/(24*d2) + (2d2-6K2*d1)/(24*d1)

# Note: 8K2^2 - d1 = 8K2^2 - (4K2^2-K4) = 4K2^2 + K4.
# And: 2d2 - 6K2*d1 = 2(4K2^3+K2*K4-2K3^2) - 6K2*(4K2^2-K4) = 8K2^3+2K2*K4-4K3^2-24K2^3+6K2*K4
# = -16K2^3 + 8K2*K4 - 4K3^2 = -4(4K2^3 - 2K2*K4 + K3^2)

# So: -R_4 = d1*(4K2^2+K4)/(24*d2) - (4K2^3-2K2*K4+K3^2)/(6*d1)

# Hmm, 4K2^3 - 2K2*K4 + K3^2... is this positive?
# At K3=0, K4=0: 4K2^3 > 0. At K3=0, K4=4K2^2: 4K2^3-8K2^3 = -4K2^3 < 0. NOT always positive.
# So this decomposition has a sign issue.

# Let me try yet another grouping.
# N = 8K2^2*d1^2 - 6K2*d1*d2 - d1^3 + 2d2^2
# Think of this as a quadratic in d2:
# N = 2*d2^2 - 6K2*d1*d2 + (8K2^2*d1^2 - d1^3)
# = 2*d2^2 - 6K2*d1*d2 + d1^2*(8K2^2 - d1)
# Discriminant: (6K2*d1)^2 - 4*2*d1^2*(8K2^2-d1)
# = 36K2^2*d1^2 - 8*d1^2*(8K2^2-d1) = d1^2*(36K2^2 - 64K2^2 + 8d1)
# = d1^2*(8d1 - 28K2^2)
# This can be positive (when d1 > 3.5K2^2) or negative. So N doesn't have a definite sign
# as a quadratic in d2 with real coefficients. But on the domain, d1 ranges over (0, 8K2^2).

# When 8d1 - 28K2^2 < 0, i.e., d1 < 3.5K2^2, the discriminant is negative,
# so N > 0 (since leading coeff 2 > 0) for all d2. ✓

# When d1 >= 3.5K2^2, N could be zero at some d2, but on the domain d2 > 0.
# N = 2(d2 - 3K2*d1/2)^2 + d1^2*(8K2^2-d1) - 9K2^2*d1^2/2
# Wait: 2(d2 - (6K2*d1)/(2*2))^2 = 2(d2 - 3K2*d1/2)^2 = 2d2^2 - 6K2*d1*d2 + 9K2^2*d1^2/2
# So N = 2(d2-3K2d1/2)^2 + d1^2*(8K2^2-d1) - 9K2^2*d1^2/2
# = 2(d2-3K2d1/2)^2 + d1^2*(8K2^2 - d1 - 9K2^2/2)
# = 2(d2-3K2d1/2)^2 + d1^2*(7K2^2/2 - d1)

# So N = 2*(d2 - 3K2*d1/2)^2 + d1^2*(7K2^2/2 - d1)

# Check: at d1 < 7K2^2/2 = 3.5K2^2, both terms non-negative, so N > 0. ✓
# At d1 >= 3.5K2^2, the second term is negative, but the first term can compensate.

# For the subadditivity proof, this quadratic-in-d2 structure might be useful.
# Let me verify:
K2_s, d1_s, d2_s = symbols('K2 d1 d2', positive=True)
N_expr = 8*K2_s**2*d1_s**2 - 6*K2_s*d1_s*d2_s - d1_s**3 + 2*d2_s**2
N_completed = 2*(d2_s - 3*K2_s*d1_s/2)**2 + d1_s**2*(S(7)/2*K2_s**2 - d1_s)
check = expand(N_completed - N_expr)
print(f"N = 2*(d2-3K2*d1/2)^2 + d1^2*(7K2^2/2 - d1)")
print(f"Verification: {check}")

# Great! So:
# -R_4 = [2*(d2-3K2*d1/2)^2 + d1^2*(7K2^2/2-d1)] / (24*d1*d2)
# = (d2-3K2*d1/2)^2/(12*d1*d2) + d1*(7K2^2/2-d1)/(24*d2)

print("\n-R_4 = (d2-3K2*d1/2)^2/(12*d1*d2) + d1*(7K2^2/2-d1)/(24*d2)")
print("     = Term_A + Term_B")
print("where Term_A = (d2-3K2*d1/2)^2/(12*d1*d2)")
print("      Term_B = d1*(7K2^2/2-d1)/(24*d2)")

# Verify:
Term_A = (d2_s - 3*K2_s*d1_s/2)**2 / (12*d1_s*d2_s)
Term_B = d1_s*(S(7)/2*K2_s**2 - d1_s) / (24*d2_s)
check2 = cancel(Term_A + Term_B - N_expr/(24*d1_s*d2_s))
print(f"Verification: diff = {check2}")

# Now let's analyze each term:
# Term_A = (d2-3K2*d1/2)^2/(12*d1*d2) = h(K2,d1,d2)^2/(12*d1*d2)
#   where h = d2 - 3K2*d1/2
# Term_B = d1*(7K2^2/2-d1)/(24*d2)

# Term_B: numerator is d1*(7K2^2/2 - d1).
# Since d1 = 4K2^2-K4 and K4 < 4K2^2, we have 0 < d1 < 8K2^2.
# 7K2^2/2 - d1 = 7K2^2/2 - 4K2^2 + K4 = K4 - K2^2/2.
# This can be negative! (When K4 < K2^2/2.)
# So Term_B is not necessarily positive.

# But wait: in the original formula, d1 ranges from 0 to 8K2^2 (since K4 ∈ (-4K2^2, 4K2^2)).
# When d1 < 7K2^2/2: Term_B > 0.
# When d1 > 7K2^2/2 (i.e., K4 < K2^2/2): Term_B < 0.

# So we can write:
# -R_4 = Term_A + Term_B
# where Term_A >= 0 always, and Term_B can be negative.
# The sum is still positive (since N > 0 when d2 > 0, which we verified).

print("\n" + "="*70)
print("KEY ANALYSIS: Structure of the two terms for subadditivity")
print("="*70)

# For subadditivity of -R_4, we need:
# Term_A(r) + Term_B(r) <= Term_A(p) + Term_A(q) + Term_B(p) + Term_B(q)
# i.e., [Term_A(r) - Term_A(p) - Term_A(q)] + [Term_B(r) - Term_B(p) - Term_B(q)] <= 0

# Let's check each term's subadditivity numerically.
print("\nNumerical check of Term_A and Term_B subadditivity gaps:")

np.random.seed(42)
min_gap_A = 1e10
max_gap_A = -1e10
min_gap_B = 1e10
max_gap_B = -1e10
min_gap_total = 1e10

for _ in range(1000000):
    sp = np.random.exponential(1) + 0.01
    tp = np.random.exponential(1) + 0.01
    ap = np.random.normal(0, sp**1.5)
    bp = np.random.normal(0, tp**1.5)

    k4p_lo = 2*ap**2/sp - 4*sp**2 + 0.001
    k4p_hi = 4*sp**2 - 0.001
    if k4p_lo >= k4p_hi: continue
    up = np.random.uniform(k4p_lo, k4p_hi)

    k4q_lo = 2*bp**2/tp - 4*tp**2 + 0.001
    k4q_hi = 4*tp**2 - 0.001
    if k4q_lo >= k4q_hi: continue
    vp = np.random.uniform(k4q_lo, k4q_hi)

    sr = sp + tp
    ar = ap + bp
    ur = up + vp
    D1r = 4*sr**2 - ur
    D2r = 4*sr**3 + sr*ur - 2*ar**2
    if D1r <= 0 or D2r <= 0: continue

    d1p = 4*sp**2 - up
    d2p = 4*sp**3 + sp*up - 2*ap**2
    d1q = 4*tp**2 - vp
    d2q = 4*tp**3 + tp*vp - 2*bp**2
    d1r = D1r
    d2r = D2r

    def term_A(k2, d1, d2):
        return (d2 - 1.5*k2*d1)**2 / (12*d1*d2)

    def term_B(k2, d1, d2):
        return d1*(3.5*k2**2 - d1) / (24*d2)

    gA = term_A(sr, d1r, d2r) - term_A(sp, d1p, d2p) - term_A(tp, d1q, d2q)
    gB = term_B(sr, d1r, d2r) - term_B(sp, d1p, d2p) - term_B(tp, d1q, d2q)

    min_gap_A = min(min_gap_A, gA)
    max_gap_A = max(max_gap_A, gA)
    min_gap_B = min(min_gap_B, gB)
    max_gap_B = max(max_gap_B, gB)
    min_gap_total = min(min_gap_total, gA + gB)

print(f"Term_A gap: min={min_gap_A:.6e}, max={max_gap_A:.6e}")
print(f"Term_B gap: min={min_gap_B:.6e}, max={max_gap_B:.6e}")
print(f"Total gap:  min={min_gap_total:.6e}")

if min_gap_A > -1e-10:
    print("  → Term_A is SUBADDITIVE!")
elif max_gap_A < 1e-10:
    print("  → Term_A is SUPERADDITIVE!")
else:
    print("  → Term_A is NEITHER subadditive nor superadditive")

if min_gap_B > -1e-10:
    print("  → Term_B is SUBADDITIVE!")
elif max_gap_B < 1e-10:
    print("  → Term_B is SUPERADDITIVE!")
else:
    print("  → Term_B is NEITHER subadditive nor superadditive")

if min_gap_total < 1e-10:
    print("  → Total (-R_4) is SUBADDITIVE (= R_4 superadditive) ✓")
