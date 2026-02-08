"""
PROVER-11 Part 4: Fix domain constraints and re-verify k3=0 case.

When k3=0:
- Factor 1: 4*K2^2 - K4 > 0  =>  K4 < 4*K2^2
- Factor 2: 4*K2^3 + K2*K4 > 0  =>  K4 > -4*K2^2

So u ∈ (-4s^2, 4s^2) and v ∈ (-4t^2, 4t^2).

But for the sum: K2_r = s+t, K4_r = u+v
- Factor 1: 4*(s+t)^2 - (u+v) > 0
- Factor 2: 4*(s+t)^3 + (s+t)*(u+v) > 0  =>  u+v > -4*(s+t)^2

These are all satisfied by the individual constraints... wait, factor 2
for sum: u+v > -4*(s+t)^2. Since u > -4s^2 and v > -4t^2:
u+v > -4(s^2+t^2) > -4(s+t)^2. ✓

And factor 1: u+v < 4s^2 + 4t^2 < 4(s+t)^2. ✓ (since 4s^2+4t^2 < 4(s+t)^2 = 4s^2+8st+4t^2)

OK so the domain is correct. Let me re-check why there are violations.
Actually, R4(K2, 0, K4) = -K4^2/(24*K2*(4*K2^2-K4)). But the full denominator
is 24*(4*K2^2-K4)*(4*K2^3+K2*K4). With K3=0, the second factor is
4*K2^3+K2*K4 = K2*(4*K2^2+K4).

So R4(K2, 0, K4) = (-4*K2^2*K4^2 - K4^3) / (24*(4*K2^2-K4)*K2*(4*K2^2+K4))
                 = -K4^2*(4*K2^2+K4) / (24*K2*(4*K2^2-K4)*(4*K2^2+K4))

Wait, let me recompute from scratch.
"""
from sympy import (symbols, expand, factor, collect, cancel, Poly,
                   Rational, simplify, S, Matrix, det, sqrt, solve,
                   numer, denom)
import numpy as np

s, t = symbols('s t', positive=True)
u, v = symbols('u v', real=True)

def R4_full(K2, K3, K4):
    num = -16*K2**3*K3**2 - 4*K2**2*K4**2 + 20*K2*K3**2*K4 - 8*K3**4 - K4**3
    den = 24*(4*K2**2 - K4)*(4*K2**3 + K2*K4 - 2*K3**2)
    return num, den

# k3=0 case
num_k30 = -4*s**2*u**2 - u**3
den_k30 = 24*(4*s**2 - u)*(4*s**3 + s*u)
print("R4(s, 0, u):")
print(f"  num = {expand(num_k30)} = {factor(num_k30)}")
print(f"  den = {expand(den_k30)} = {factor(den_k30)}")
print(f"  = {cancel(num_k30/den_k30)}")

# Factor: num = -u^2*(4s^2 + u), den = 24*s*(4s^2-u)*(4s^2+u)
# So R4 = -u^2*(4s^2+u) / (24*s*(4s^2-u)*(4s^2+u)) = -u^2 / (24*s*(4s^2-u))
# PROVIDED 4s^2 + u ≠ 0, which is guaranteed since u > -4s^2.
# Hmm wait, if u is very close to -4s^2, then 4s^2+u ≈ 0.
# Actually when K3=0: den = 24*K2*(4*K2^2-K4)*(4*K2^2+K4)
# The factor (4*K2^2+K4) is positive when K4 > -4*K2^2.
# When K4 = -4*K2^2 exactly, the denominator is 0. So K4 > -4*K2^2 strictly.

# After cancellation:
# R4(s,0,u) = -u^2 / (24*s*(4s^2 - u))
# This is correct.

print(f"\nSimplified: R4(s,0,u) = -u^2/(24*s*(4*s^2-u))")
print(f"  Verified: {cancel(num_k30/den_k30)}")

# Now, the superadditivity gap
# Delta = R4(s+t, 0, u+v) - R4(s, 0, u) - R4(t, 0, v)
# = -(u+v)^2/(24*(s+t)*(4*(s+t)^2-(u+v))) + u^2/(24*s*(4*s^2-u)) + v^2/(24*t*(4*t^2-v))
# Need this >= 0

# Let's verify numerically with proper domain
np.random.seed(42)
violations = 0
n_trials = 100000
min_gap = 1e10
for _ in range(n_trials):
    s_val = np.random.exponential(1) + 0.01
    t_val = np.random.exponential(1) + 0.01
    # u in (-4s^2, 4s^2), EXCLUDING boundaries
    u_val = np.random.uniform(-4*s_val**2 + 0.001, 4*s_val**2 - 0.001)
    v_val = np.random.uniform(-4*t_val**2 + 0.001, 4*t_val**2 - 0.001)

    # Check sum constraints
    sr = s_val + t_val
    ur = u_val + v_val
    if ur >= 4*sr**2 - 0.001 or ur <= -4*sr**2 + 0.001:
        continue

    def R4_num(K2, K4):
        return -K4**2 / (24 * K2 * (4*K2**2 - K4))

    r_sum = R4_num(sr, ur)
    r_p = R4_num(s_val, u_val)
    r_q = R4_num(t_val, v_val)
    gap = r_sum - r_p - r_q
    min_gap = min(min_gap, gap)

    if gap < -1e-10:
        violations += 1
        if violations <= 5:
            print(f"VIOLATION: s={s_val:.6f}, t={t_val:.6f}, u={u_val:.6f}, v={v_val:.6f}")
            print(f"  4s^2={4*s_val**2:.4f}, 4t^2={4*t_val**2:.4f}")
            print(f"  R_sum={r_sum:.8f}, R_p={r_p:.8f}, R_q={r_q:.8f}, gap={gap:.6e}")

print(f"\n{violations} violations in {n_trials} trials, min gap = {min_gap:.6e}")

# Now let me check: is my Cauchy-Schwarz argument actually correct?
# Need: (u+v)^2/((s+t)*C) <= u^2/(s*A) + v^2/(t*B)
# where A = 4s^2-u, B = 4t^2-v, C = 4(s+t)^2-(u+v)

# The Cauchy-Schwarz gives:
# (u+v)^2 <= (s*A + t*B) * (u^2/(s*A) + v^2/(t*B))
# And we showed s*A + t*B <= (s+t)*C
# So (u+v)^2 <= (s+t)*C * (u^2/(s*A) + v^2/(t*B))
# Hence (u+v)^2/((s+t)*C) <= u^2/(s*A) + v^2/(t*B)

# Wait, the direction: s*A + t*B <= (s+t)*C means
# 1/(s*A+t*B) >= 1/((s+t)*C)
# And from CS: (u+v)^2/(s*A+t*B) <= u^2/(s*A) + v^2/(t*B)
# So: (u+v)^2/((s+t)*C) <= (u+v)^2/(s*A+t*B) <= u^2/(s*A) + v^2/(t*B)  ✓

# Let me verify the key inequality symbolically one more time
print("\n\nVERIFYING s*A + t*B <= (s+t)*C:")
A_s = 4*s**2 - u
B_s = 4*t**2 - v
C_s = 4*(s+t)**2 - u - v

diff = expand((s+t)*C_s - s*A_s - t*B_s)
print(f"  (s+t)*C - s*A - t*B = {diff}")
print(f"  = {factor(diff)}")
# Should be 12st(s+t) - tu - sv
# = t*(12s^2 + 12st - u) - sv
# Since u < 4s^2: 12s^2 + 12st - u > 12s^2 + 12st - 4s^2 = 8s^2 + 12st > 0
# And v < 4t^2: sv < 4st^2
# So diff > t*(8s^2+12st) - 4st^2 = 8s^2*t + 12st^2 - 4st^2 = 8s^2*t + 8st^2 = 8st(s+t) > 0

# But I got violations! Let me look at the violation more carefully.
print("\n\nDEBUGGING VIOLATIONS:")
# Example: s=0.0140, t=1.4298, u=0.0015, v=-4.5390
s_val, t_val, u_val, v_val = 0.0140, 1.4298, 0.0015, -4.5390
# Check domain:
print(f"  u < 4s^2? {u_val} < {4*s_val**2} = {u_val < 4*s_val**2}")
print(f"  u > -4s^2? {u_val} > {-4*s_val**2} = {u_val > -4*s_val**2}")
print(f"  v < 4t^2? {v_val} < {4*t_val**2} = {v_val < 4*t_val**2}")
print(f"  v > -4t^2? {v_val} > {-4*t_val**2} = {v_val > -4*t_val**2}")

# WAIT: v = -4.5390 and 4t^2 = 4*1.4298^2 = 8.177. So -4t^2 = -8.177.
# v = -4.5390 > -8.177. ✓
# But the ACTUAL domain for polynomials also requires disc > 0.
# For k3=0: disc_4 = 16*k2^4*k4 - 128*k2^2*k4^2 + 256*k4^3 (with e2,e3,e4 substituted)
# Actually we need to check if these are actually valid cumulant vectors.

# The domain from the denominator factors is:
# (1) 4*K2^2 - K4 > 0
# (2) 4*K2^3 + K2*K4 - 2*K3^2 > 0, i.e., K4 > -4*K2^2 + 2*K3^2/K2
# With K3=0: K4 > -4*K2^2, i.e., u > -4s^2, v > -4t^2

# But the discriminant of the quartic also needs to be positive!
# disc = 16*e2^4*e4 - 4*e2^3*e3^2 - 128*e2^2*e4^2 + 144*e2*e3^2*e4 - 27*e3^4 + 256*e4^3
# With e2 = -3*K2/2, e3 = K3/2, e4 = -3*K4/32 + 3*K2^2/16
# With K3=0: disc = function of K2, K4

# But the statement is about cumulant superadditivity, not about discriminant.
# The superadditivity should hold for ALL cumulant vectors in the domain,
# not just those coming from real-rooted polynomials.

# Actually, re-reading the problem statement:
# "R_4(k2,k3,k4) = ... Domain constraints (both factors in denominator are positive):
# k2 > 0, k4 < 4*k2^2, 2*k3^2 < 4*k2^3 + k2*k4"
# So the domain is defined by the denominator being positive, not by disc > 0.

# So my domain IS correct. Let me look at what's happening with the "violations".

r_sum = -u_val**2 / (24 * (s_val+t_val) * (4*(s_val+t_val)**2 - (u_val+v_val)))

# Wait - is the formula correct? Let me be more careful.
print(f"\n  R4(s,0,u) formula: -u^2/(24*s*(4s^2-u))")
print(f"  But the ORIGINAL R4 has denominator 24*(4K2^2-K4)*(4K2^3+K2*K4)")
print(f"  = 24*K2*(4K2^2-K4)*(4K2^2+K4)")
print(f"  Numerator: -4K2^2*K4^2 - K4^3 = -K4^2*(4K2^2+K4)")
print(f"  So R4 = -K4^2*(4K2^2+K4) / (24*K2*(4K2^2-K4)*(4K2^2+K4))")
print(f"  = -K4^2 / (24*K2*(4K2^2-K4))")
print(f"  This cancellation requires 4K2^2+K4 != 0.")

# Check: when v = -4.539 and t = 1.4298:
# 4t^2 + v = 4*1.4298^2 + (-4.539) = 8.177 - 4.539 = 3.638 > 0 ✓

# So the cancellation is fine. Let me recompute more carefully.
r_p = -(u_val**2) / (24*s_val*(4*s_val**2 - u_val))
r_q = -(v_val**2) / (24*t_val*(4*t_val**2 - v_val))
sr = s_val + t_val
ur = u_val + v_val
r_sum = -(ur**2) / (24*sr*(4*sr**2 - ur))
gap = r_sum - r_p - r_q
print(f"\n  R_p = {r_p:.10f}")
print(f"  R_q = {r_q:.10f}")
print(f"  R_sum = {r_sum:.10f}")
print(f"  Gap = {gap:.10e}")
print(f"  4s^2-u = {4*s_val**2-u_val:.6f}")
print(f"  4t^2-v = {4*t_val**2-v_val:.6f}")
print(f"  4(s+t)^2-(u+v) = {4*sr**2-ur:.6f}")

# Check CS step:
A_val = 4*s_val**2 - u_val
B_val = 4*t_val**2 - v_val
C_val = 4*sr**2 - ur
sA_tB = s_val*A_val + t_val*B_val
stC = sr * C_val
print(f"\n  s*A + t*B = {sA_tB:.6f}")
print(f"  (s+t)*C  = {stC:.6f}")
print(f"  s*A+t*B <= (s+t)*C? {sA_tB <= stC}")
print(f"  12st(s+t)-tu-sv = {12*s_val*t_val*sr - t_val*u_val - s_val*v_val:.6f}")
