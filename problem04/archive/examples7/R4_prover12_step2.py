"""
PROVER-12 Step 2: Work with the (d1,d2) coordinate decomposition.

KEY FINDING from Step 1:
  -R_4 = (8*K2^2*d1^2 - 6*K2*d1*d2 - d1^3 + 2*d2^2) / (24*d1*d2)

where d1 = D1 = 4K2^2 - K4 > 0
      d2 = D2 = 4K2^3 + K2*K4 - 2K3^2 > 0

This means:
  -R_4 = (8*K2^2*d1)/(24*d2) - (6*K2)/(24) - d1^2/(24*d2) + (2*d2)/(24*d1)
       = K2^2*d1/(3*d2) - K2/4 - d1^2/(24*d2) + d2/(12*d1)

Wait, let me redo this more carefully.
"""
from sympy import (symbols, expand, factor, cancel, simplify, collect,
                   Rational, S, Poly, together, apart, fraction, sqrt)
import numpy as np

K2, K3, K4 = symbols('K2 K3 K4', positive=True)
d1, d2 = symbols('d1 d2', positive=True)

# -R_4 in (d1, d2):
neg_R4_d = (8*K2**2*d1**2 - 6*K2*d1*d2 - d1**3 + 2*d2**2) / (24*d1*d2)

# Separate into terms:
# = 8*K2^2*d1/(24*d2) - 6*K2/(24) - d1^2/(24*d2) + 2*d2/(24*d1)
# = K2^2*d1/(3*d2) - K2/4 - d1^2/(24*d2) + d2/(12*d1)
print("Partial fraction in d1, d2:")
print("  -R_4 = K2^2*d1/(3*d2) - K2/4 - d1^2/(24*d2) + d2/(12*d1)")

# Verify:
check = cancel(K2**2*d1/(3*d2) - K2/4 - d1**2/(24*d2) + d2/(12*d1) - neg_R4_d)
print(f"  Verification: {check}")

# Now let's substitute back: d1 = 4K2^2 - K4, d2 = 4K2^3 + K2*K4 - 2*K3^2
# Note: d2 = K2*(4K2^2 + K4) - 2K3^2 = K2*(D1 + 2*K4 + K4... wait)
# d2 = 4K2^3 + K2*K4 - 2K3^2
# d1 = 4K2^2 - K4
# So d2 = K2*d1 + 2*K2*K4 - 2*K3^2... no:
# K2*d1 = 4K2^3 - K2*K4
# d2 - K2*d1 = (4K2^3 + K2*K4 - 2K3^2) - (4K2^3 - K2*K4) = 2K2*K4 - 2K3^2
# So d2 = K2*d1 + 2*(K2*K4 - K3^2)

# Also: K4 = 4K2^2 - d1, so:
# d2 = K2*d1 + 2*(K2*(4K2^2-d1) - K3^2)
# = K2*d1 + 8K2^3 - 2K2*d1 - 2K3^2
# = -K2*d1 + 8K2^3 - 2K3^2
# = K2*(8K2^2 - d1) - 2K3^2

# Now consider the ratio d2/d1 and d1/d2 that appear.
# Let's define:
#   p = d1/(K2^2), q = d2/(K2^3)
# Then:
#   d1 = K2^2 * p, d2 = K2^3 * q
#   0 < p < 8 (since d1 = 4K2^2-K4 and K4 > -4K2^2 from D2>0 when K3=0, but could be larger)
#   Actually p = (4K2^2-K4)/K2^2 = 4 - K4/K2^2 = 4 - β, so 0 < p < 8 (since β < 4 and β > -4)
#   And q = d2/K2^3 = 4 + K4/K2^2 - 2K3^2/K2^3 = 4 + β - 2α where α = K3^2/K2^3
#   = 8 - p - 2α

# Then: -R_4 = K2^2*(K2^2*p)/(3*K2^3*q) - K2/4 - (K2^2*p)^2/(24*K2^3*q) + K2^3*q/(12*K2^2*p)
#            = K2*p/(3*q) - K2/4 - K2*p^2/(24*q) + K2*q/(12*p)
#            = K2 * [p/(3q) - 1/4 - p^2/(24q) + q/(12p)]
#            = K2 * [(8p^2 - 6pq - p^3 + 2q^2)/(24*p*q)]
# Good, consistent with the original formula scaled by K2.

# For superadditivity of R_4, we need:
# R_4(p+q) >= R_4(p) + R_4(q)
# i.e., for -R_4 (which we're analyzing as the subadditive target):
# -R_4(p+q) <= -R_4(p) + -R_4(q)

# In the (K2, d1, d2) coordinates:
# K2_r = s+t, d1_r = 4(s+t)^2 - (u+v), d2_r = 4(s+t)^3 + (s+t)(u+v) - 2(a+b)^2
# K2_p = s, d1_p = 4s^2 - u, d2_p = 4s^3 + su - 2a^2
# K2_q = t, d1_q = 4t^2 - v, d2_q = 4t^3 + tv - 2b^2

# The key question: is (8*K2^2*d1^2 - 6*K2*d1*d2 - d1^3 + 2*d2^2)/(24*d1*d2)
# subadditive as a function of (K2, d1, d2)?

# Let me think about this differently. We have 4 terms:
# A = K2^2*d1/(3*d2) = (1/3) * K2^2*d1/d2
# B = -K2/4  [linear, subadditive]
# C = -d1^2/(24*d2) = -(1/24) * d1^2/d2  [note the negative sign!]
# D = d2/(12*d1) = (1/12) * d2/d1

# For A = K2^2*d1/d2: this looks like it could be convex (since K2^2 is convex,
# d1/d2 might be too). If convex, it's superadditive, which is the WRONG direction.

# For C = -d1^2/d2: -f^2/g. If g is concave and f convex, this could be concave.

# For D = d2/d1: ratio of concave/convex = concave? Not necessarily.

# ACTUALLY, wait. Let me reconsider the whole approach.
# The core insight is:
# -R_4 = N(K2, d1, d2) / (24*d1*d2)
# where N = 8K2^2*d1^2 - 6K2*d1*d2 - d1^3 + 2d2^2

# The denominator 24*d1*d2 is NOT linear in the original variables (K2,K3,K4).
# But d1 and d2 ARE linear in K4 and K3^2 respectively (fixing K2).

# CRITICAL OBSERVATION: d1 and d2 are actually SUPERADDITIVE!
# d1_r = 4(s+t)^2 - (u+v) and d1_p + d1_q = 4s^2-u + 4t^2-v
# d1_r - (d1_p + d1_q) = 4(s+t)^2 - 4s^2 - 4t^2 = 8st > 0
# So d1_r >= d1_p + d1_q. ✓

# d2_r = 4(s+t)^3 + (s+t)(u+v) - 2(a+b)^2
# d2_p + d2_q = 4s^3+su-2a^2 + 4t^3+tv-2b^2
# d2_r - (d2_p+d2_q) = 4[(s+t)^3-s^3-t^3] + (s+t)(u+v)-su-tv - 2[(a+b)^2-a^2-b^2]
# = 4[3s^2t+3st^2] + tu+sv - 4ab
# = 12st(s+t) + tu + sv - 4ab

# Is this >= 0? By AM-GM: 4ab <= 2(a^2/s + b^2/t)*... no, need to be more careful.
# Actually by Cauchy-Schwarz: (a+b)^2 <= ??? This depends on the domain constraints.

# Let me check numerically.
print("\n" + "="*70)
print("NUMERICAL: Check superadditivity of d2")
print("="*70)

np.random.seed(42)
min_d2_gap = 1e10
for _ in range(500000):
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

    d1p = 4*sp**2 - up
    d1q = 4*tp**2 - vp
    d1r = 4*sr**2 - ur

    d2p = 4*sp**3 + sp*up - 2*ap**2
    d2q = 4*tp**3 + tp*vp - 2*bp**2
    d2r = 4*sr**3 + sr*ur - 2*ar**2

    if d1r <= 0 or d2r <= 0: continue
    if d1p <= 0 or d2p <= 0: continue
    if d1q <= 0 or d2q <= 0: continue

    d2_gap = d2r - d2p - d2q
    min_d2_gap = min(min_d2_gap, d2_gap)

print(f"min(d2_r - d2_p - d2_q) = {min_d2_gap:.6f}")
if min_d2_gap > -1e-10:
    print("d2 IS superadditive!")
else:
    print("d2 is NOT superadditive (min gap negative)")

# Let me also look at what d2_r - d2_p - d2_q equals analytically:
# = 12st(s+t) + tu + sv - 4ab
# The issue is the -4ab term. But from D2 constraints:
# 2a^2 < 4s^3 + su, so a^2 < 2s^3 + su/2
# Similarly b^2 < 2t^3 + tv/2
# By AM-GM: 4|ab| <= 4*sqrt((2s^3+su/2)(2t^3+tv/2))
# This isn't obviously smaller than 12st(s+t)+tu+sv

# But note that a and b can have OPPOSITE signs, making -4ab positive!
# The worst case is when a and b have the same sign.
# When a,b > 0: -4ab < 0. Need 12st(s+t)+tu+sv > 4ab.

# From constraints: a^2 < 2s^3 + su/2 < 2s^3 + 2s^3 = 4s^3 (using u < 4s^2)
# So |a| < 2s^{3/2}. Similarly |b| < 2t^{3/2}.
# Then 4ab < 4*2s^{3/2}*2t^{3/2} = 16*(st)^{3/2}
# And 12st(s+t) >= 12*2*st*sqrt(st) = 24*(st)^{3/2}  (by AM-GM: s+t >= 2sqrt(st))
# So 12st(s+t) > 16*(st)^{3/2} when 12st(s+t) > 16(st)^{3/2}
# i.e., 12*(s+t) > 16*sqrt(st), i.e., 3(s+t) > 4*sqrt(st)
# By AM-GM: s+t >= 2*sqrt(st), so 3(s+t) >= 6*sqrt(st) > 4*sqrt(st). ✓

# Wait... 12st(s+t) + tu + sv >= 12st(s+t) > 16(st)^{3/2} > 4ab.
# So d2 IS superadditive! And we can prove it.

print("\n" + "="*70)
print("ANALYTICAL PROOF: d2 is superadditive")
print("="*70)
print("""
d2_r - d2_p - d2_q = 12st(s+t) + tu + sv - 4ab

Claim: this is >= 0 on the domain.

Proof: From D2 > 0 and D1 > 0, we have |a| < 2*s^{3/2} and |b| < 2*t^{3/2}.
(Since a^2 < 2s^3 + su/2 < 2s^3 + 2s^3 = 4s^3.)

Case 1: ab <= 0. Then -4ab >= 0 and 12st(s+t)+tu+sv > 0. Done.

Case 2: ab > 0. Then 4ab <= 4*|a|*|b| < 4*(2s^{3/2})*(2t^{3/2}) = 16*(st)^{3/2}.
  Also: 12st(s+t) >= 12st*2*sqrt(st) = 24*(st)^{3/2} > 16*(st)^{3/2}.
  And tu + sv > 0 (since u > 2a^2/s - 4s^2 and the positive parts dominate...
  WAIT: u and v can be negative!)

Hmm, tu + sv could be negative. Let me reconsider.
""")

# u and v can be negative (K4 can be negative). Let me check the worst case.
# u > 2a^2/s - 4s^2 (from D2 > 0: 4s^3 + su - 2a^2 > 0 ⟹ u > 2a^2/s - 4s^2)
# The most negative u can be is approximately -4s^2 (when a=0).
# And v > 2b^2/t - 4t^2.
# So tu > t*(2a^2/s - 4s^2) = 2ta^2/s - 4s^2*t
# And sv > s*(2b^2/t - 4t^2) = 2sb^2/t - 4st^2
# So tu + sv > 2ta^2/s + 2sb^2/t - 4st(s+t)

# Total: 12st(s+t) + tu + sv - 4ab > 12st(s+t) + 2ta^2/s + 2sb^2/t - 4st(s+t) - 4ab
# = 8st(s+t) + 2ta^2/s + 2sb^2/t - 4ab
# = 8st(s+t) + 2(ta/sqrt(s) - b*sqrt(s))^2/t + 2b^2*s/t + 2sb^2/t - 2*(ta/sqrt(s))^2/t - 4ab
# Hmm, getting complicated. Let me try completing the square differently.

# 2ta^2/s + 2sb^2/t - 4ab = 2*(sqrt(t/s)*a - sqrt(s/t)*b)^2 + 4ab - 4ab
# Wait: 2ta^2/s + 2sb^2/t - 4ab = 2*(a*sqrt(t/s))^2 + 2*(b*sqrt(s/t))^2 - 4ab
# = 2*(a*sqrt(t/s) - b*sqrt(s/t))^2 + 4ab*sqrt(t/s)*sqrt(s/t) - 4ab + ... no.
# Let p=a*sqrt(t/s), q=b*sqrt(s/t). Then:
# 2p^2 + 2q^2 - 4ab = 2p^2 + 2q^2 - 4*p*sqrt(s/t)*q*sqrt(t/s) = 2p^2+2q^2-4pq = 2(p-q)^2
# YES! So 2ta^2/s + 2sb^2/t - 4ab = 2(a*sqrt(t/s) - b*sqrt(s/t))^2 >= 0.

print("KEY IDENTITY:")
print("  2t*a^2/s + 2s*b^2/t - 4ab = 2*(a*sqrt(t/s) - b*sqrt(s/t))^2 >= 0")
print()
print("Therefore:")
print("  d2_r - d2_p - d2_q = 12st(s+t) + tu + sv - 4ab")
print("  > 8st(s+t) + 2t*a^2/s + 2s*b^2/t - 4ab")
print("  = 8st(s+t) + 2*(a*sqrt(t/s) - b*sqrt(s/t))^2")
print("  >= 0")
print()
print("PROVED: d2 is superadditive on the domain.")

# Let me also verify this bound numerically.
print("\n" + "="*70)
print("NUMERICAL VERIFICATION of the bound")
print("="*70)

np.random.seed(42)
min_bound = 1e10
for _ in range(500000):
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

    d2_gap = 12*sp*tp*(sp+tp) + tp*up + sp*vp - 4*ap*bp
    # The claimed lower bound:
    bound = 8*sp*tp*sr + 2*tp*ap**2/sp + 2*sp*bp**2/tp - 4*ap*bp

    min_bound = min(min_bound, bound)

print(f"min(8st(s+t) + 2ta^2/s + 2sb^2/t - 4ab) = {min_bound:.6f}")
print(f"This should be >= 0: {'YES' if min_bound > -1e-10 else 'NO'}")

# Now I have: d1 superadditive (gap = 8st > 0), d2 superadditive (gap > 0).
# -R_4 = N(K2, d1, d2)/(24*d1*d2), need to show subadditivity.
# With N = 8K2^2*d1^2 - 6K2*d1*d2 - d1^3 + 2*d2^2.

# Since d1, d2 are superadditive, they are NOT easy substitution variables
# for proving subadditivity. The denominator d1*d2 is superadditive (product of
# superadditive positives), which goes the WRONG way for subadditivity.

# Let me think about this differently. Maybe I should work directly with
# the 4-term decomposition and prove subadditivity of specific combinations.

print("\n" + "="*70)
print("STEP 6: Try 3-term decomposition")
print("="*70)

# From the (d1, d2) form:
# -R_4 = K2^2*d1/(3*d2) - K2/4 - d1^2/(24*d2) + d2/(12*d1)

# Rewrite: let r = d2/(K2*d1). Then d2 = r*K2*d1, so:
# -R_4 = K2^2*d1/(3*r*K2*d1) - K2/4 - d1^2/(24*r*K2*d1) + r*K2*d1/(12*d1)
#       = K2/(3r) - K2/4 - d1/(24*r*K2) + r*K2/12
#       = K2*[1/(3r) - 1/4 + r/12] - d1/(24*r*K2)

# Hmm, this is getting circular. Let me try a completely different approach.

# KEY IDEA: Use the homogeneity!
# R_4(K2,K3,K4) = K2 * h(α, β) where α = K3^2/K2^3, β = K4/K2^2.
# h(α,β) = (-16α - 4β^2 + 20αβ - 8α^2 - β^3)/(24*(4-β)*(4+β-2α))

# Superadditivity of R_4 = superadditivity of K2*h(K3^2/K2^3, K4/K2^2).
# This is a PERSPECTIVE function f(x) = x*g(y/x) where x = K2 and
# y = (K3^2/K2^2, K4/K2) ??? Not quite...

# Actually, let's be more careful. If we set:
#   x = K2 (weight)
#   u = K3^2/K2^3 (= α)
#   v = K4/K2^2 (= β)
# Then R_4 = x * h(u, v).
# Superadditivity means: (x_p+x_q)*h(u_r,v_r) >= x_p*h(u_p,v_p) + x_q*h(u_q,v_q)
# where u_r = (a_p+a_q)^2/x_r^3, v_r = (K4_p+K4_q)/x_r^2.
# But u_r is NOT a weighted average of u_p and u_q! (Because of the square in the numerator.)

# So the perspective trick doesn't directly apply. However...

# Let me try yet another parameterization. Set:
#   w1 = K2  (weight, positive)
#   w2 = K3/K2  (ratio)
#   w3 = K4/K2  (ratio)
# Then K3 = w1*w2, K4 = w1*w3.
# R_4(w1, w1*w2, w1*w3) = w1 * h(w2^2/w1, w3/w1)  ... still messy

# OK, different approach entirely. Let me use the SUBSTITUTION approach.
# Parameterize using POSITIVE SLACK variables:
#   d1 = 4K2^2 - K4 > 0
#   e = K3 (unconstrained real)
#   d2 = 4K2^3 + K2*K4 - 2*K3^2 = K2*(8K2^2-d1) - 2*e^2 > 0

# Then K4 = 4K2^2 - d1, K3 = e, K3^2 = e^2.
# Variables: (K2, e, d1) with K2 > 0, d1 > 0, K2*(8K2^2-d1) - 2e^2 > 0
# (The last constraint is d2 > 0.)

# And: -R_4 = (8K2^2*d1^2 - 6K2*d1*d2 + 2*d2^2 - d1^3) / (24*d1*d2)
# where d2 = K2*(8K2^2-d1) - 2e^2.

# Superadditivity: K2_r = s+t, e_r = a+b, d1_r = 4(s+t)^2 - (u+v)
# d1 is a function of K2 and K4, NOT directly of our independent slack.

# I'm going in circles. Let me just try the direct computation approach:
# compute the gap numerator polynomial in 6 variables, factor it, and look for structure.

print("Computing gap numerator directly in 6 variables...")
s, t, a, b, u, v = symbols('s t a b u v')

def neg_R4(K2, K3, K4):
    num = 16*K2**3*K3**2 + 4*K2**2*K4**2 - 20*K2*K3**2*K4 + 8*K3**4 + K4**3
    D1 = 4*K2**2 - K4
    D2 = 4*K2**3 + K2*K4 - 2*K3**2
    return num, D1, D2

# Subadditivity gap of -R_4:
# We need -R_4(p+q) - [-R_4(p)] - [-R_4(q)] <= 0
# i.e., num_r/(D1_r*D2_r) - num_p/(D1_p*D2_p) - num_q/(D1_q*D2_q) <= 0
# (ignoring the factor 24)
# Clearing denominators (all positive):
# num_r * D1_p*D2_p * D1_q*D2_q - num_p * D1_r*D2_r * D1_q*D2_q - num_q * D1_r*D2_r * D1_p*D2_p
# should be <= 0 (for subadditivity of -R_4, equivalently superadditivity of R_4)

num_r, D1_r, D2_r = neg_R4(s+t, a+b, u+v)
num_p, D1_p, D2_p = neg_R4(s, a, u)
num_q, D1_q, D2_q = neg_R4(t, b, v)

print("Expanding gap numerator (slow)...")
# This is going to be a HUGE polynomial. Let me try a different route.
# Instead of raw expansion, work with specific substitutions.

# Actually, let me try the (d1, d2) approach more seriously.
# In (d1, d2) coords: N(K2, d1, d2) = 8K2^2*d1^2 - 6K2*d1*d2 - d1^3 + 2d2^2
# -R_4 = N/(24*d1*d2)

# Define:
# For point p: K2=s, d1_p=4s^2-u, d2_p=4s^3+su-2a^2
# For point q: K2=t, d1_q=4t^2-v, d2_q=4t^3+tv-2b^2
# For sum: K2=s+t, d1_r=4(s+t)^2-(u+v), d2_r=4(s+t)^3+(s+t)(u+v)-2(a+b)^2

# Subadditivity gap of -R_4 = N_r/(24*d1_r*d2_r) - N_p/(24*d1_p*d2_p) - N_q/(24*d1_q*d2_q)
# = [N_r*d1_p*d2_p*d1_q*d2_q - N_p*d1_r*d2_r*d1_q*d2_q - N_q*d1_r*d2_r*d1_p*d2_p] / (24*prod)

# For SUBadditivity we need the numerator to be <= 0.
# For SUPERadditivity of R_4 we need R_4 gap >= 0, which means the numerator of
# the gap in -R_4 must be <= 0.

# Wait, let me reconsider. R_4 = -neg_R4_num/(24*D1*D2) = -(neg_R4_num)/(24*D1*D2).
# Since D1, D2 > 0, we have R_4 < 0 when neg_R4_num > 0.
# Superadditivity: R_4(r) >= R_4(p) + R_4(q)
# -neg_R4_num_r/(24*D1_r*D2_r) >= -neg_R4_num_p/(24*D1_p*D2_p) + ...
# Multiply through by -24:
# neg_R4_num_r/(D1_r*D2_r) <= neg_R4_num_p/(D1_p*D2_p) + neg_R4_num_q/(D1_q*D2_q)
# Clear (positive) denominators:
# neg_R4_num_r * D1_p*D2_p * D1_q*D2_q <= neg_R4_num_p * D1_r*D2_r * D1_q*D2_q + neg_R4_num_q * D1_r*D2_r * D1_p*D2_p
# Gap that should be <= 0:
# G = neg_R4_num_r * D1_p*D2_p * D1_q*D2_q - neg_R4_num_p * D1_r*D2_r * D1_q*D2_q - neg_R4_num_q * D1_r*D2_r * D1_p*D2_p

# This is a monster polynomial. But using the (d1,d2) representation where
# N = 8K2^2*d1^2 - 6K2*d1*d2 - d1^3 + 2d2^2, we can write:
# G_scaled = N_r*d1_p*d2_p*d1_q*d2_q - N_p*d1_r*d2_r*d1_q*d2_q - N_q*d1_r*d2_r*d1_p*d2_p
# and we need G_scaled <= 0, where N_r = 8(s+t)^2*d1_r^2 - 6(s+t)*d1_r*d2_r - d1_r^3 + 2*d2_r^2, etc.

# Even this is huge. Let me try a NUMERICAL EXPLORATION first to understand
# where the gap is tight (close to 0).

print("\n" + "="*70)
print("NUMERICAL: Where is the superadditivity gap tightest?")
print("="*70)

np.random.seed(42)
tight_cases = []
n_trials = 2000000
valid = 0

for _ in range(n_trials):
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

    valid += 1

    def R4_eval(k2, k3, k4):
        n = -16*k2**3*k3**2 - 4*k2**2*k4**2 + 20*k2*k3**2*k4 - 8*k3**4 - k4**3
        d = 24*(4*k2**2-k4)*(4*k2**3+k2*k4-2*k3**2)
        return n/d

    rp = R4_eval(sp, ap, up)
    rq = R4_eval(tp, bp, vp)
    rr = R4_eval(sr, ar, ur)
    gap = rr - rp - rq

    if gap < 1e-4:
        # Compute scaled variables
        alpha_p = ap**2/sp**3
        beta_p = up/sp**2
        alpha_q = bp**2/tp**3
        beta_q = vp/tp**2
        tight_cases.append((gap, sp, tp, ap, bp, up, vp, alpha_p, beta_p, alpha_q, beta_q))

tight_cases.sort()
print(f"Valid trials: {valid}")
print(f"Cases with gap < 1e-4: {len(tight_cases)}")
print(f"\nTightest 10 cases:")
print(f"{'gap':>12} {'s':>8} {'t':>8} {'a':>8} {'b':>8} {'u':>8} {'v':>8} {'α_p':>8} {'β_p':>8} {'α_q':>8} {'β_q':>8}")
for c in tight_cases[:10]:
    print(f"{c[0]:12.2e} {c[1]:8.4f} {c[2]:8.4f} {c[3]:8.4f} {c[4]:8.4f} {c[5]:8.4f} {c[6]:8.4f} {c[7]:8.4f} {c[8]:8.4f} {c[9]:8.4f} {c[10]:8.4f}")

# Check patterns
print("\nAnalyzing tight cases...")
if tight_cases:
    for c in tight_cases[:5]:
        gap, sp, tp, ap, bp, up, vp = c[:7]
        sr = sp + tp
        # What are the domain parameters close to?
        d1p = 4*sp**2 - up
        d2p = 4*sp**3 + sp*up - 2*ap**2
        d1q = 4*tp**2 - vp
        d2q = 4*tp**3 + tp*vp - 2*bp**2
        print(f"  gap={gap:.2e}: s/t={sp/tp:.3f}, a*b/|a||b|={'N/A' if ap*bp==0 else (ap*bp)/(abs(ap)*abs(bp)):.3f}, d1p/d1q={d1p/d1q:.3f}, d2p/d2q={d2p/d2q:.3f}")
