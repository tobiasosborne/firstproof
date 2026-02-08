"""
PROVER-11 Part 7: Decompose R4 and analyze f2.

R4 = f1 + f2 where:
  f1 = -K4^2 / (24*K2*D1)           [k3=0 piece]
  f2 = -2*K3^2*G / (24*K2*D1*D2)    [k3-dependent piece]

where G = 8*K2^4 - 10*K2^2*K4 + 4*K2*K3^2 + K4^2

Superadditivity of R4 = superadditivity of f1 + superadditivity of f2
                        + cross terms.

Actually, since superadditivity is NOT preserved under addition in general,
we cannot just prove superadditivity of f1 and f2 separately.

Better approach: work with the full Delta = R4_r - R4_p - R4_q directly.

Alternative: maybe a generalized Cauchy-Schwarz on the full fraction.

Let me try the approach: write -R4 = K4^2/(24*K2*D1) + 2*K3^2*G/(24*K2*D1*D2)
and try to show -R4 is subadditive.
"""
from sympy import (symbols, expand, factor, collect, cancel, Poly,
                   Rational, simplify, S, numer, denom, apart)
import numpy as np

K2, K3, K4 = symbols('K2 K3 K4', positive=True)

D1 = 4*K2**2 - K4
D2 = 4*K2**3 + K2*K4 - 2*K3**2

# Let's examine f2 more closely
G = 8*K2**4 - 10*K2**2*K4 + 4*K2*K3**2 + K4**2

print("G = 8*K2^4 - 10*K2^2*K4 + 4*K2*K3^2 + K4^2")
print("\nIs G always positive on the domain?")

# Try completing the square in K4:
# G = K4^2 - 10*K2^2*K4 + 8*K2^4 + 4*K2*K3^2
# = (K4 - 5*K2^2)^2 - 25*K2^4 + 8*K2^4 + 4*K2*K3^2
# = (K4 - 5*K2^2)^2 - 17*K2^4 + 4*K2*K3^2
# This is NOT always positive since 17*K2^4 >> 4*K2*K3^2 typically.

# Actually, on the domain K4 < 4*K2^2:
# G = K4^2 - 10*K2^2*K4 + 8*K2^4 + 4*K2*K3^2
# At K4 = 4*K2^2: G = 16*K2^4 - 40*K2^4 + 8*K2^4 + 4*K2*K3^2 = -16*K2^4 + 4*K2*K3^2
# This can be negative.
# At K4 = 0: G = 8*K2^4 + 4*K2*K3^2 > 0
# At K4 = -4*K2^2: G = 16*K2^4 + 40*K2^4 + 8*K2^4 + 4*K2*K3^2 = 64*K2^4 + 4*K2*K3^2 > 0

# So G can change sign. This makes f2 complex.
# When G > 0: f2 < 0 (same sign as f1)
# When G < 0: f2 > 0 (opposite sign to f1)

print("G can be negative (e.g., at K4 close to 4K2^2 and small K3).")
print("So f2 can be positive, i.e., f2 partially cancels f1.")

# Let me try a completely different approach.
# Instead of decomposing R4, try the generalized Cauchy-Schwarz directly.
#
# Approach: partial fractions of -R4.
# -R4 = (16*K2^3*K3^2 + 4*K2^2*K4^2 - 20*K2*K3^2*K4 + 8*K3^4 + K4^3)
#        / (24*D1*D2)
#
# Let's do partial fractions in K4, treating D1 and D2 as linear in K4.
# D1 = -K4 + 4*K2^2, linear in K4
# D2 = K2*K4 + (4*K2^3 - 2*K3^2), linear in K4

print("\n" + "="*70)
print("APPROACH: Partial fractions in K4")
print("="*70)

# -R4_num = 16*K2^3*K3^2 + 4*K2^2*K4^2 - 20*K2*K3^2*K4 + 8*K3^4 + K4^3
neg_R4_num = 16*K2**3*K3**2 + 4*K2**2*K4**2 - 20*K2*K3**2*K4 + 8*K3**4 + K4**3

# Let's compute -R4 = neg_R4_num / (24*D1*D2)
# Partial fractions: A_coeff/(24*D1) + B_coeff/(24*D2) + C_poly/24
neg_R4 = neg_R4_num / (24*D1*D2)
pf = apart(neg_R4, K4)
print(f"\nPartial fractions of -R4 in K4:")
print(f"  {pf}")

# Hmm, this might be messy. Let me try another route.

print("\n" + "="*70)
print("APPROACH: Homogeneous form + 2-variable reduction")
print("="*70)

# Use weight-2 homogeneity: R4(K2, K3, K4) = K2 * f(u, v)
# where u = K3/K2^{3/2}, v = K4/K2^2
# Domain: v < 4, 2u^2 < v + 4

u, v = symbols('u v', real=True)

# f(u,v) = p(u,v)/q(u,v) with
# p = -8u^4 + 20u^2*v - 16u^2 - v^3 - 4v^2
# q = 24*(v-4)*(2u^2-v-4) = 24*(4-v)*(v+4-2u^2) [note signs]

p_uv = -8*u**4 + 20*u**2*v - 16*u**2 - v**3 - 4*v**2
q_uv = 24*(4 - v)*(v + 4 - 2*u**2)  # positive on domain

print(f"f(u,v) = p/q where:")
print(f"  p = {p_uv}")
print(f"  q = {q_uv}")

# Now, superadditivity of R4 means:
# (s+t)*f(u_r, v_r) >= s*f(u_p, v_p) + t*f(u_q, v_q)
# where u_r = (K3_p+K3_q)/(s+t)^{3/2}, v_r = (u_p*s^2+v_p*... complicated

# This is getting messy because the scaled variables don't add simply.
# Let me go back to the direct approach.

print("\n" + "="*70)
print("APPROACH: Direct analysis of full Delta numerator")
print("="*70)

# The key question: can we extend the Cauchy-Schwarz argument from k3=0?
# For k3=0: -R4 = K4^2/(24*K2*D1), and we used CS on K4^2/(K2*D1).
#
# For general case: -R4 = (16*K2^3*K3^2 + 4*K2^2*K4^2 - 20*K2*K3^2*K4 + 8*K3^4 + K4^3)
#                         / (24*(4*K2^2-K4)*(4*K2^3+K2*K4-2*K3^2))
#
# The numerator of -R4 is:
# N = 16*K2^3*K3^2 + 4*K2^2*K4^2 - 20*K2*K3^2*K4 + 8*K3^4 + K4^3
#
# Let me try to write N as a quadratic form in (K3^2, K4) with K2-dependent coeffs.
# Actually, N has K3^4 and K4^3 terms, so it's not a quadratic form.
# But in terms of (K3, K4): N has terms K3^4, K3^2*K4, K3^2, K4^2, K4^3.
# It IS a quadratic in K3^2 with K4-dependent coefficients:
# N = 8*(K3^2)^2 + (16*K2^3 - 20*K2*K4)*K3^2 + (4*K2^2*K4^2 + K4^3)

# Or viewing as quadratic in K3^2:
# N = 8*X^2 + (16*K2^3-20*K2*K4)*X + K4^2*(4*K2^2+K4)  where X = K3^2

X = symbols('X', positive=True)
N_in_X = 8*X**2 + (16*K2**3 - 20*K2*K4)*X + K4**2*(4*K2**2 + K4)
print(f"\nN as quadratic in X = K3^2:")
print(f"  N = {N_in_X}")

# Discriminant of this quadratic in X:
disc_X = (16*K2**3 - 20*K2*K4)**2 - 4*8*K4**2*(4*K2**2 + K4)
disc_X_exp = expand(disc_X)
print(f"\nDiscriminant = {disc_X_exp}")
print(f"  factored = {factor(disc_X_exp)}")

# If discriminant < 0: N > 0 always (since leading coeff 8 > 0)
# If discriminant >= 0: N has real roots in X = K3^2, but we need X >= 0.

# disc = 256*K2^6 - 640*K2^4*K4 + 400*K2^2*K4^2 - 128*K2^2*K4^2 - 32*K4^3
# = 256*K2^6 - 640*K2^4*K4 + 272*K2^2*K4^2 - 32*K4^3
disc_X_exp2 = expand(256*K2**6 - 640*K2**4*K4 + 272*K2**2*K4**2 - 32*K4**3)
print(f"  expanded = {disc_X_exp2}")
print(f"  match: {simplify(disc_X_exp - disc_X_exp2) == 0}")

disc_X_fac = factor(disc_X_exp)
print(f"  factored = {disc_X_fac}")
