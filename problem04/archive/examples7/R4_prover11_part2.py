"""
PROVER-11 Part 2: k3=0 numerator decomposition.

We need to show Num >= 0 where:
Num = 16*s^6*v^2 + 48*s^5*t*v^2 + ... (big polynomial in s,t,u,v)

Domain: s,t > 0, u < 4s^2, v < 4t^2, u+v < 4(s+t)^2

Strategy: substitute a = 4s^2 - u > 0, b = 4t^2 - v > 0, then
express Num in terms of s, t, a, b (all positive).
"""
from sympy import (symbols, expand, factor, collect, cancel, Poly,
                   Rational, simplify, S, Matrix, det, sqrt)

s, t, a, b = symbols('s t a b', positive=True)
u, v = symbols('u v', real=True)

# The numerator
Num = (16*s**6*v**2 + 48*s**5*t*v**2 + 48*s**4*t**2*v**2
       - 8*s**4*u*v**2 - 4*s**4*v**3
       - 32*s**3*t**3*u*v - 8*s**3*t*u*v**2
       + 48*s**2*t**4*u**2 - 12*s**2*t**2*u**2*v
       - 12*s**2*t**2*u*v**2 + s**2*u**2*v**2 + s**2*u*v**3
       + 48*s*t**5*u**2 - 8*s*t**3*u**2*v
       + 16*t**6*u**2 - 4*t**4*u**3 - 8*t**4*u**2*v
       + t**2*u**3*v + t**2*u**2*v**2)

# Substitute u = 4s^2 - a, v = 4t^2 - b
Num_ab = expand(Num.subs({u: 4*s**2 - a, v: 4*t**2 - b}))
print("Num in (s, t, a, b) variables:")

# Collect by monomials in (a, b)
Num_poly = Poly(Num_ab, a, b)
print("\nAs polynomial in (a, b):")
for monom, coeff in sorted(Num_poly.as_dict().items()):
    coeff_fac = factor(coeff)
    print(f"  a^{monom[0]} * b^{monom[1]}: {coeff_fac}")

print("\n" + "="*70)
print("Analyzing signs of each coefficient...")
print("="*70)
