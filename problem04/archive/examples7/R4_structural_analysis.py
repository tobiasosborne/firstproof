"""
PROVER-9: Structural analysis of R_4 for superadditivity proof.

Key findings from Step 1:
- R_4 = P(k2,k3,k4) / Q(k2,k3,k4)
- P = -16*k2^3*k3^2 - 4*k2^2*k4^2 + 20*k2*k3^2*k4 - 8*k3^4 - k4^3
- Q = 24*(4*k2^2 - k4)*(4*k2^3 + k2*k4 - 2*k3^2)
- Domain: k2 > 0, Q > 0 (positive denominator)

Strategy:
- For the n=3 case, R_3 = -2*k3^2/(27*k2^2), the superadditivity proof used:
  k3_p^2/k2_p^2 + k3_q^2/k2_q^2 >= (k3_p+k3_q)^2/(k2_p+k2_q)^2
  This followed from Jensen/convexity of x^2.

- For n=4, we need: R_4(k_p + k_q) >= R_4(k_p) + R_4(k_q)
  where k_p = (k2_p, k3_p, k4_p), k_q = (k2_q, k3_q, k4_q)
  and k_p + k_q = (k2_p+k2_q, k3_p+k3_q, k4_p+k4_q)

Let's investigate the structure.
"""
from sympy import (symbols, Rational, simplify, factor, expand, cancel,
                   Matrix, det, S, Poly, together, numer, denom, collect,
                   sqrt, pprint, groebner, resultant, solve, diff,
                   Symbol, oo, limit, series, Expr, degree, eye,
                   zeros as sym_zeros, ones, diag, BlockMatrix, tensorproduct,
                   Sum, Product, binomial, factorial, Abs, sign, Piecewise,
                   GreaterThan, StrictGreaterThan, Le, Ge, Lt, Gt, Ne, Eq,
                   Interval, FiniteSet, Union, Intersection, EmptySet)
import sympy

k2, k3, k4 = symbols('k2 k3 k4', positive=True)
s, t = symbols('s t', positive=True)
a, b, c, d, e, f = symbols('a b c d e f', real=True)

# R_4 numerator and denominator
P = -16*k2**3*k3**2 - 4*k2**2*k4**2 + 20*k2*k3**2*k4 - 8*k3**4 - k4**3
# Note: denominator factor (4*k2^2 - k4) should be positive for valid polynomials
# because for real-rooted quartic, the discriminant must be positive
Q_factor1 = 4*k2**2 - k4
Q_factor2 = 4*k2**3 + k2*k4 - 2*k3**2
Q = 24 * Q_factor1 * Q_factor2

print("=" * 70)
print("STRUCTURAL ANALYSIS OF R_4")
print("=" * 70)

# First, let's understand the homogeneity
# k2 has weight 2, k3 has weight 3, k4 has weight 4
# P terms: k2^3*k3^2 (6+6=12), k2^2*k4^2 (4+8=12), k2*k3^2*k4 (2+6+4=12), k3^4 (12), k4^3 (12)
# So P is homogeneous of weight 12
# Q_factor1: k2^2 (4), k4 (4) -> homogeneous weight 4
# Q_factor2: k2^3 (6), k2*k4 (2+4=6), k3^2 (6) -> homogeneous weight 6
# Q: weight 4+6 = 10
# R_4 = P/Q has weight 12-10 = 2
# This matches: 1/Phi_4 has weight 2 (same as k2), and R_4 = 1/Phi_4 - k2/12

print("\nHomogeneity check: R_4 has cumulant weight 2. Confirmed.")

# Introduce scaling: k3 = k2^{3/2} * u, k4 = k2^2 * v
# Then R_4 = k2 * f(u, v) for some function f
# This is because R_4 has weight 2

u, v = symbols('u v', real=True)

P_scaled = expand(P.subs({k3: k2**(S(3)/2)*u, k4: k2**2*v}))
Q_scaled = expand(Q.subs({k3: k2**(S(3)/2)*u, k4: k2**2*v}))

# Factor out powers of k2
# P should have k2^6 (since weight 12, and k2 has weight 2)
# Q should have k2^5 (since weight 10)
P_reduced = simplify(P_scaled / k2**6)
Q_reduced = simplify(Q_scaled / k2**5)

print(f"\nAfter scaling k3 = k2^(3/2)*u, k4 = k2^2*v:")
print(f"  P = k2^6 * p(u,v) where p = {P_reduced}")
print(f"  Q = k2^5 * q(u,v) where q = {Q_reduced}")
print(f"\n  R_4 = k2 * p(u,v)/q(u,v)")

# So f(u,v) = p(u,v)/q(u,v)
p_uv = P_reduced
q_uv = Q_reduced

print(f"\n  f(u,v) = ({p_uv}) / ({q_uv})")

# Factor p and q
p_f = factor(p_uv)
q_f = factor(q_uv)
print(f"\n  p(u,v) factored = {p_f}")
print(f"  q(u,v) factored = {q_f}")

print("\n\n" + "=" * 70)
print("DOMAIN CONSTRAINTS")
print("=" * 70)

# For a real-rooted centered quartic x^4 + e2*x^2 - e3*x + e4:
# disc > 0 requires specific constraints on e2, e3, e4
# In cumulant space:
# disc = 27*(32*k2^6 - 32*k2^3*k3^2 - 6*k2^2*k4^2 + 24*k2*k3^2*k4 - 8*k3^4 - k4^3)/128
# = 27*k2^6*(32 - 32*u^2 - 6*v^2 + 24*u^2*v - 8*u^4 - v^3)/128

# disc > 0 on the valid domain (real-rooted with simple roots)
# N_4 is related to Phi_4 * disc, and Phi_4 > 0, so N_4 > 0 iff disc > 0 when Phi_4 > 0.

# Actually: Phi_4 = N_4/disc. For Phi_4 > 0, we need N_4/disc > 0.
# Since disc > 0 for simple roots, N_4 > 0.
# N_4 = 81*k2^5*(4*k2^2 - k4)*(4*k2^3 + k2*k4 - 2*k3^2)/16
# = 81*k2^5 * (4 - v) * (4 + v - 2*u^2) / 16 (in scaled variables)

print("\nN_4 in scaled variables:")
print(f"  N_4 = 81*k2^5/16 * (4-v)*(4+v-2*u^2)")
print("  For N_4 > 0 (since k2 > 0): need (4-v)*(4+v-2*u^2) > 0")
print("  This means: either both factors positive or both negative.")
print()

# The discriminant in scaled variables:
disc_scaled = 32 - 32*u**2 - 6*v**2 + 24*u**2*v - 8*u**4 - v**3
print(f"disc (scaled) = {disc_scaled}")
print("  Need disc > 0 for simple roots.")

# Important constraint: 4 - v > 0, i.e., k4 < 4*k2^2
# This is a necessary condition.
# For centered real-rooted polynomials:
# kappa_4 = excess kurtosis * kappa_2^2 (roughly)
# For real-rooted polynomials, kappa_4 can be positive or negative,
# but bounded by 4*k2^2

print("\nKey domain: v < 4 AND 4+v-2*u^2 > 0 (since both factors must be positive)")
print("  i.e., k4 < 4*k2^2 AND k4 > 2*k3^2/k2 - 4*k2^2")

print("\n\n" + "=" * 70)
print("SUPERADDITIVITY IN SCALED VARIABLES")
print("=" * 70)

# Superadditivity of R_4: for k_p, k_q with k2_p, k2_q > 0,
# R_4(k_p + k_q) >= R_4(k_p) + R_4(k_q)
#
# In scaled variables:
# R_4(k2, k3, k4) = k2 * f(k3/k2^{3/2}, k4/k2^2)
#
# With k2_p = s, k2_q = t (s,t > 0):
# R_4(s+t, k3_p+k3_q, k4_p+k4_q) = (s+t)*f(u_r, v_r)
# where u_r = (k3_p+k3_q)/(s+t)^{3/2}, v_r = (k4_p+k4_q)/(s+t)^2
#
# This is complicated because the scaled variables don't add linearly.
# Let me try a different approach.

print("\nApproach 1: Direct difference analysis")
print("Need: R_4(k_p + k_q) - R_4(k_p) - R_4(k_q) >= 0")
print()

# Let's use separate symbols for the two inputs
k2p, k3p, k4p = symbols('k2p k3p k4p', real=True)
k2q, k3q, k4q = symbols('k2q k3q k4q', real=True)

def R4_func(K2, K3, K4):
    num = -16*K2**3*K3**2 - 4*K2**2*K4**2 + 20*K2*K3**2*K4 - 8*K3**4 - K4**3
    den = 24*(4*K2**2 - K4)*(4*K2**3 + K2*K4 - 2*K3**2)
    return num / den

R4p = R4_func(k2p, k3p, k4p)
R4q = R4_func(k2q, k3q, k4q)
R4r = R4_func(k2p + k2q, k3p + k3q, k4p + k4q)

# The difference we need >= 0
Delta = R4r - R4p - R4q

# This is extremely complex. Let's try special cases first.
print("\nSpecial case 1: k3 = 0 (symmetric distributions)")
Delta_k3_0 = Delta.subs({k3p: 0, k3q: 0})
Delta_k3_0 = cancel(Delta_k3_0)
print(f"\nDelta|_{'{k3=0}'} = {Delta_k3_0}")

Delta_k3_0_num = numer(Delta_k3_0)
Delta_k3_0_den = denom(Delta_k3_0)

Delta_k3_0_num_expanded = expand(Delta_k3_0_num)
Delta_k3_0_den_expanded = expand(Delta_k3_0_den)
print(f"\nNum = {Delta_k3_0_num_expanded}")
print(f"Den = {Delta_k3_0_den_expanded}")

# Factor
print(f"\nNum factored = {factor(Delta_k3_0_num)}")
print(f"Den factored = {factor(Delta_k3_0_den)}")

print("\n\nSpecial case 2: k4 = 0")
Delta_k4_0 = Delta.subs({k4p: 0, k4q: 0})
Delta_k4_0 = cancel(Delta_k4_0)
print(f"\nDelta|_{'{k4=0}'} = {Delta_k4_0}")

Delta_k4_0_num = expand(numer(Delta_k4_0))
Delta_k4_0_den = expand(denom(Delta_k4_0))
print(f"\nNum = {Delta_k4_0_num}")
print(f"\nNum factored = {factor(numer(Delta_k4_0))}")
print(f"Den factored = {factor(denom(Delta_k4_0))}")

print("\n\nSpecial case 3: k4_p = k4_q = 0 AND k3_q = 0 (one-parameter slice)")
Delta_1d = Delta.subs({k4p: 0, k4q: 0, k3q: 0})
Delta_1d = cancel(Delta_1d)
print(f"\nDelta|_{'{k4=0,k3_q=0}'} = {Delta_1d}")
print(f"Num factored = {factor(numer(Delta_1d))}")
print(f"Den factored = {factor(denom(Delta_1d))}")
