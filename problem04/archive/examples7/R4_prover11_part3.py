"""
PROVER-11 Part 3: SOS decomposition for k3=0 numerator.

Num(a,b) = 128*s^5*t^2*(s+t)*b + 16*s^4*t*(s+t)*b^2
          + 128*s^2*t^5*(s+t)*a + 16*s*t^4*(s+t)*a^2
          - 16*s^2*t^2*(s+3t)*(3s+t)*ab
          + 8*s^2*t*(s-t)*a*b^2 - 8*s*t^2*(s-t)*a^2*b
          + s^2*a*b^3 + (s^2+t^2)*a^2*b^2 + t^2*a^3*b

where a = 4s^2 - u > 0, b = 4t^2 - v > 0, s, t > 0.

Additionally: u + v < 4(s+t)^2, so a + b > 4(s+t)^2 - 4s^2 - 4t^2 = 8st,
i.e., a + b > 8st.

Key insight: the negative term -16*s^2*t^2*(s+3t)*(3s+t)*ab must be absorbed
by the positive terms. Let's try to find an explicit decomposition.

Strategy: Group terms and apply AM-GM / Schur-type inequalities.
"""
from sympy import (symbols, expand, factor, collect, cancel, Poly,
                   Rational, simplify, S, Matrix, det, sqrt, solve,
                   numer, denom)
import numpy as np

s, t, a, b = symbols('s t a b', positive=True)

# Let's try a weighted AM-GM approach.
# The key negative term is: -16*s^2*t^2*(3s+t)*(s+3t)*ab
# = -16*s^2*t^2*(3s^2+10st+3t^2)*ab

# Positive linear terms: 128*s^5*t^2*(s+t)*b and 128*s^2*t^5*(s+t)*a
# These are "large" when a or b is small.

# Let me try: write Num = Q(a,b) where we group cleverly.
#
# Terms with only b (no a): 128*s^5*t^2*(s+t)*b + 16*s^4*t*(s+t)*b^2
# Terms with only a (no b): 128*s^2*t^5*(s+t)*a + 16*s*t^4*(s+t)*a^2
# Cross terms: -16*s^2*t^2*(s+3t)*(3s+t)*ab + ...higher order in a,b
#
# Idea: use the constraint a+b > 8st to help.
# Write a = 8st*alpha + ..., but this seems hard.
#
# Alternative: try writing as a quadratic form in (a, b) with
# polynomial coefficients, ignoring higher-order terms, and show the
# quadratic part dominates.
#
# Actually, let's try: view as quadratic in (sqrt(a)*b, a*sqrt(b), etc)

# Let me try a direct approach: write Num as sum of non-negative terms.
# The key is that s^2*a*b^3 + t^2*a^3*b = ab(s^2*b^2 + t^2*a^2)
# and (s^2+t^2)*a^2*b^2 = a^2*b^2*(s^2+t^2).
# These higher-order terms in a,b can absorb the negative cross term
# when a,b are large enough.

# Let me try substituting a = alpha*s, b = beta*t (dimensionless)
# Then a + b > 8st means alpha*s + beta*t > 8st

# Actually let's try a completely different approach.
# Go back to original variables and try to write the gap as
# a perfect square plus non-negative terms.

# APPROACH: Cauchy-Schwarz on the original fraction
# R4(s,0,u) = -u^2/(24*s*(4s^2-u))
# = -u^2/(24*s*A) where A = 4s^2 - u > 0
#
# So -R4(s,0,u) = u^2/(24*s*A)
#
# Superadditivity of R4 <==> subadditivity of -R4 = u^2/(24*s*A)
# i.e., (u+v)^2/(24*(s+t)*C) <= u^2/(24*s*A) + v^2/(24*s*B)
# where A = 4s^2-u, B = 4t^2-v, C = 4(s+t)^2-(u+v)
#
# Equivalently: (u+v)^2/((s+t)*C) <= u^2/(s*A) + v^2/(t*B)

print("APPROACH: Cauchy-Schwarz on u^2/(s*A)")
print("="*70)
print()
print("Need: (u+v)^2 / ((s+t)*C) <= u^2/(s*A) + v^2/(t*B)")
print("where A=4s^2-u, B=4t^2-v, C=4(s+t)^2-u-v")
print()

# Let's write x = u/sqrt(s*A), y = v/sqrt(t*B)
# Then u^2/(s*A) = x^2, v^2/(t*B) = y^2
# We need: (u+v)^2/((s+t)*C) <= x^2 + y^2
#
# By Cauchy-Schwarz: (u+v)^2 <= (s*A + t*B)*(u^2/(s*A) + v^2/(t*B))
# So it suffices to show: (s*A + t*B)/((s+t)*C) <= 1
# i.e., s*A + t*B <= (s+t)*C
#
# s*A + t*B = s*(4s^2-u) + t*(4t^2-v) = 4s^3 - su + 4t^3 - tv
# (s+t)*C = (s+t)*(4(s+t)^2 - u - v) = 4(s+t)^3 - (s+t)(u+v)
#
# Difference: (s+t)*C - s*A - t*B
# = 4(s+t)^3 - (s+t)(u+v) - 4s^3 + su - 4t^3 + tv
# = 4[(s+t)^3 - s^3 - t^3] - (s+t)(u+v) + su + tv
# = 4[3s^2*t + 3s*t^2] - (s+t)(u+v) + su + tv
# = 12st(s+t) - su - tu - sv - tv + su + tv
# = 12st(s+t) - tu - sv
# = 12st(s+t) - t*u - s*v

# Is 12st(s+t) - tu - sv >= 0?
# Use u < 4s^2, v < 4t^2:
# 12st(s+t) - tu - sv > 12st(s+t) - 4s^2*t - 4t^2*s
# = 12s^2*t + 12st^2 - 4s^2*t - 4st^2
# = 8s^2*t + 8st^2 = 8st(s+t) > 0  ✓

# But we need the exact bound, not just positivity.
# Actually we need: 12st(s+t) - tu - sv >= 0 on the domain.
# This is NOT guaranteed because u, v could be negative!
# Wait - u = k4_p can be negative. Let me check the domain.

print("Key computation: (s+t)*C - s*A - t*B = 12st(s+t) - tu - sv")
print()
print("If u, v can be negative, this could be larger than 12st(s+t).")
print("If u, v are large and positive... need u < 4s^2, v < 4t^2.")
print()

# Can u be negative? k4 is the 4th cumulant, which can be negative.
# But we need 4*k2^3 + k2*k4 - 2*k3^2 > 0. With k3=0: k4 > -4*k2^2.
# So u ∈ (-4s^2, 4s^2), v ∈ (-4t^2, 4t^2).
# Actually, 4k2^3 + k2*k4 > 0 means k4 > -4k2^2. Yes.

# For the negative case u = -4s^2 + eps:
# 12st(s+t) - t*(-4s^2) - sv = 12st(s+t) + 4s^2*t - sv
# = t*(12s^2+12st+4s^2) - sv = t*(16s^2+12st) - sv
# This is certainly > 0 for v < 4t^2: sv < 4st^2, and t*(16s^2+12st) > 4st^2.

# For the positive case u -> 4s^2, v -> 4t^2:
# 12st(s+t) - t*4s^2 - s*4t^2 = 12s^2*t+12st^2-4s^2*t-4st^2 = 8st(s+t) > 0

# So 12st(s+t) - tu - sv >= 8st(s+t) > 0 on the ENTIRE domain!
# Wait, I need to verify this more carefully.

# 12st(s+t) - tu - sv
# = 12st(s+t) - tu - sv
# Since u <= 4s^2 - eps and v <= 4t^2 - eps (strict for simple roots):
# >= 12st(s+t) - t*4s^2 - s*4t^2  (when both are at max)
# = 12s^2*t + 12st^2 - 4s^2*t - 4st^2 = 8s^2*t + 8st^2 = 8st(s+t) > 0

# But what if u is very negative? Then -tu is positive (t > 0),
# so the expression gets BIGGER. So the minimum is at u=4s^2, v=4t^2.

# CONFIRMED: 12st(s+t) - tu - sv >= 8st(s+t) > 0 on the domain.
# (Equality not achieved since u < 4s^2 strictly.)

print("RESULT: (s+t)*C - s*A - t*B = 12st(s+t) - tu - sv >= 8st(s+t) > 0")
print("Since u < 4s^2 and v < 4t^2, the minimum is at the boundary.")
print()

# Therefore s*A + t*B <= (s+t)*C, and by Cauchy-Schwarz:
# (u+v)^2 <= (s*A + t*B) * (u^2/(s*A) + v^2/(t*B))
# <= (s+t)*C * (u^2/(s*A) + v^2/(t*B))
# So (u+v)^2/((s+t)*C) <= u^2/(s*A) + v^2/(t*B)

# BUT WAIT: Cauchy-Schwarz gives (u+v)^2 <= (s*A + t*B)*(u^2/(s*A) + v^2/(t*B))
# This is just the Cauchy-Schwarz inequality (a1*b1 + a2*b2)^2 <= (a1^2+a2^2)(b1^2+b2^2)
# with a1 = sqrt(s*A), b1 = u/sqrt(s*A), a2 = sqrt(t*B), b2 = v/sqrt(t*B).
# Check: a1*b1 + a2*b2 = u + v. ✓
# (a1^2+a2^2) = s*A + t*B.
# (b1^2+b2^2) = u^2/(s*A) + v^2/(t*B).
# ✓

# So we get:
# (u+v)^2/((s+t)*C) <= (u+v)^2/(s*A+t*B) <= u^2/(s*A) + v^2/(t*B)
#
# The first inequality uses s*A+t*B <= (s+t)*C.
# The second is Cauchy-Schwarz.

print("PROOF OF k3=0 CASE:")
print("="*70)
print()
print("Claim: R4(s+t, 0, u+v) >= R4(s, 0, u) + R4(t, 0, v)")
print()
print("Proof:")
print("  R4(K2, 0, K4) = -K4^2 / (24*K2*(4*K2^2 - K4))")
print("  Let A = 4s^2 - u > 0, B = 4t^2 - v > 0, C = 4(s+t)^2 - (u+v) > 0")
print("  -R4(K2,0,K4) = K4^2 / (24*K2*(4K2^2-K4))")
print()
print("  Need: -R4(s+t,0,u+v) <= -R4(s,0,u) + -R4(t,0,v)")
print("  i.e., (u+v)^2/((s+t)*C) <= u^2/(s*A) + v^2/(t*B)")
print()
print("  Step 1: s*A + t*B <= (s+t)*C")
print("    (s+t)*C - s*A - t*B = 12st(s+t) - tu - sv")
print("    Since u < 4s^2 and v < 4t^2:")
print("    12st(s+t) - tu - sv > 12st(s+t) - 4s^2*t - 4st^2 = 8st(s+t) > 0")
print()
print("  Step 2: By Cauchy-Schwarz (with weights sqrt(s*A), sqrt(t*B)):")
print("    (u+v)^2 <= (s*A + t*B) * (u^2/(s*A) + v^2/(t*B))")
print()
print("  Combining Steps 1 and 2:")
print("    (u+v)^2/((s+t)*C) <= (u+v)^2/(s*A+t*B) <= u^2/(s*A) + v^2/(t*B)")
print("  QED.")

print()
print("="*70)
print("SYMBOLIC VERIFICATION")
print("="*70)

# Verify Step 1 symbolically
u_sym, v_sym = symbols('u_sym v_sym', real=True)
A_sym = 4*s**2 - u_sym
B_sym = 4*t**2 - v_sym
C_sym = 4*(s+t)**2 - u_sym - v_sym

diff_expr = expand((s+t)*C_sym - s*A_sym - t*B_sym)
print(f"\n(s+t)*C - s*A - t*B = {diff_expr}")
print(f"  = {factor(diff_expr)}")

# Verify it equals 12st(s+t) - tu - sv
expected = 12*s*t*(s+t) - t*u_sym - s*v_sym
print(f"\nExpected: {expand(expected)}")
print(f"Match: {simplify(diff_expr - expected) == 0}")

# Now verify numerically
print("\n" + "="*70)
print("NUMERICAL VERIFICATION (k3=0 case)")
print("="*70)
import numpy as np
np.random.seed(42)

violations = 0
n_trials = 50000
for _ in range(n_trials):
    s_val = np.random.exponential(2)
    t_val = np.random.exponential(2)
    # u in (-4s^2, 4s^2), v in (-4t^2, 4t^2)
    u_val = np.random.uniform(-4*s_val**2 + 0.01, 4*s_val**2 - 0.01)
    v_val = np.random.uniform(-4*t_val**2 + 0.01, 4*t_val**2 - 0.01)
    # Also need u+v < 4(s+t)^2
    if u_val + v_val >= 4*(s_val + t_val)**2 - 0.01:
        continue

    def R4_num(K2, K4):
        A = 4*K2**2 - K4
        D = K2 * A  # = K2*(4K2^2 - K4)
        return -K4**2 / (24 * D)

    r_sum = R4_num(s_val + t_val, u_val + v_val)
    r_p = R4_num(s_val, u_val)
    r_q = R4_num(t_val, v_val)
    gap = r_sum - r_p - r_q

    if gap < -1e-10:
        violations += 1
        print(f"VIOLATION: s={s_val:.4f}, t={t_val:.4f}, u={u_val:.4f}, v={v_val:.4f}, gap={gap:.6e}")

print(f"\n{violations} violations in {n_trials} trials")
if violations == 0:
    print("k3=0 CASE: NUMERICALLY VERIFIED")
