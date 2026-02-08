"""
PROVER-11 Part 6: Full R_4 case - structural analysis.

R4(K2, K3, K4) = P(K2,K3,K4) / (24 * D1 * D2)
where:
  P = -16*K2^3*K3^2 - 4*K2^2*K4^2 + 20*K2*K3^2*K4 - 8*K3^4 - K4^3
  D1 = 4*K2^2 - K4
  D2 = 4*K2^3 + K2*K4 - 2*K3^2

For the n=3 case, R3 = -2*K3^2/(27*K2^2) = -(2/27) * (K3/K2)^2,
so -R3 = (2/27) * K3^2/K2^2 and subadditivity followed from
convexity of x^2 (Jensen).

For R4, let me try to decompose R4 into a "nice" part and show each
piece is superadditive.

Key idea: 1/Phi_4 = K2/12 + R4.
So R4 = 1/Phi_4 - K2/12.
The K2/12 part is LINEAR in cumulants, hence additive.
So superadditivity of R4 <==> superadditivity of 1/Phi_4.

The 1/Phi_4 formula:
1/Phi_4 = (32*K2^6 - 32*K2^3*K3^2 - 6*K2^2*K4^2 + 24*K2*K3^2*K4 - 8*K3^4 - K4^3)
          / (384*K2^5 - 192*K2^2*K3^2 - 24*K2*K4^2 + 48*K3^2*K4)

Let me factor the denominator of 1/Phi_4 differently.
"""
from sympy import (symbols, expand, factor, collect, cancel, Poly,
                   Rational, simplify, S, numer, denom)

K2, K3, K4 = symbols('K2 K3 K4')

# 1/Phi_4
num_inv = 32*K2**6 - 32*K2**3*K3**2 - 6*K2**2*K4**2 + 24*K2*K3**2*K4 - 8*K3**4 - K4**3
den_inv = 384*K2**5 - 192*K2**2*K3**2 - 24*K2*K4**2 + 48*K3**2*K4

print("1/Phi_4 numerator:")
print(f"  = {num_inv}")
print(f"  factored = {factor(num_inv)}")

print(f"\n1/Phi_4 denominator:")
print(f"  = {den_inv}")
print(f"  factored = {factor(den_inv)}")

# Now compute R4 = 1/Phi_4 - K2/12
R4_num = expand(num_inv - K2*den_inv/12)
R4_den = den_inv
print(f"\nR4 numerator (before simplification):")
print(f"  = {expand(R4_num)}")
print(f"  factored = {factor(R4_num)}")

# Cross-multiply: R4 = R4_num / R4_den
R4_cancel = cancel(R4_num / R4_den)
print(f"\nR4 (cancelled) = {R4_cancel}")

# Let me try a different decomposition.
# Write R4 = f1(K2, K4) + f2(K2, K3, K4) where f1 is the k3=0 part
# f1 = R4(K2, 0, K4) = -K4^2/(24*K2*(4*K2^2-K4))
# f2 = R4(K2, K3, K4) - R4(K2, 0, K4) = part that depends on K3

# R4 = P/(24*D1*D2)
P_full = -16*K2**3*K3**2 - 4*K2**2*K4**2 + 20*K2*K3**2*K4 - 8*K3**4 - K4**3
D1 = 4*K2**2 - K4
D2 = 4*K2**3 + K2*K4 - 2*K3**2

# R4(K2, 0, K4) = (-4K2^2*K4^2 - K4^3) / (24*D1*(4K2^3+K2*K4))
# = -K4^2*(4K2^2+K4) / (24*(4K2^2-K4)*K2*(4K2^2+K4))
# = -K4^2 / (24*K2*(4K2^2-K4))
# = -K4^2 / (24*K2*D1)

# f2 = R4 - f1 = P/(24*D1*D2) + K4^2/(24*K2*D1)
# = [P*K2 + K4^2*D2] / (24*K2*D1*D2)

f2_num = expand(P_full * K2 + K4**2 * D2)
print(f"\nf2 numerator = P*K2 + K4^2*D2:")
print(f"  = {f2_num}")
print(f"  factored = {factor(f2_num)}")

# Factor out K3^2
f2_num_div = cancel(f2_num / K3**2)
print(f"\nf2_num / K3^2 = {f2_num_div}")
print(f"  factored = {factor(f2_num_div)}")
