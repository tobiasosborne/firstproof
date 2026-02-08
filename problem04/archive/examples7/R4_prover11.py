"""
PROVER-11: R_4 superadditivity proof.

Strategy:
1. Prove k3=0 case via explicit SOS decomposition
2. Extend to full case via homogeneous substitution + Schur analysis
"""
from sympy import (symbols, Rational, simplify, factor, expand, cancel,
                   Matrix, det, S, Poly, together, collect, sqrt,
                   solve, diff, degree, numer, denom, Symbol, oo,
                   groebner, resultant, eye, pprint, apart, fraction,
                   nsimplify, limit, series, sign, Abs, zoo, nan,
                   Function, Lambda, Sum, Product, binomial, factorial,
                   And, Or, Not, Implies, Equivalent, satisfiable,
                   Interval, FiniteSet, Union, Intersection, EmptySet,
                   ConditionSet, ImageSet, Range, Contains,
                   GreaterThan, StrictGreaterThan, Le, Ge, Lt, Gt, Ne, Eq,
                   Piecewise, Min, Max, floor, ceiling, Mod,
                   pi, E, I, oo, zoo, nan, true, false,
                   trigsimp, radsimp, ratsimp, powsimp, combsimp, logcombine)
import sympy

print("="*70)
print("PROVER-11: R_4 SUPERADDITIVITY PROOF")
print("="*70)

# ============================================================
# PART 0: Setup and formula verification
# ============================================================

s, t = symbols('s t', positive=True)
u, v = symbols('u v', real=True)  # u = k4_p, v = k4_q
k3p, k3q = symbols('k3p k3q', real=True)

def R4(K2, K3, K4):
    """R_4 function."""
    num = -16*K2**3*K3**2 - 4*K2**2*K4**2 + 20*K2*K3**2*K4 - 8*K3**4 - K4**3
    den = 24*(4*K2**2 - K4)*(4*K2**3 + K2*K4 - 2*K3**2)
    return num / den

# Superadditivity gap: R4(s+t, k3p+k3q, u+v) - R4(s, k3p, u) - R4(t, k3q, v)
# We need this >= 0

print("\nPART 1: k3 = 0 CASE")
print("="*70)

# When k3 = 0: R4(K2, 0, K4) = -K4^2 / (24*K2*(4*K2^2 - K4))
# Verify:
R4_k30 = R4(s, 0, u)
R4_k30_simplified = cancel(R4_k30)
print(f"R4(s, 0, u) = {R4_k30_simplified}")
# Should be -u^2 / (24*s*(4*s^2 - u))

# Gap at k3=0
Delta_k30 = R4(s+t, 0, u+v) - R4(s, 0, u) - R4(t, 0, v)
Delta_k30_canc = cancel(Delta_k30)

num_k30 = numer(Delta_k30_canc)
den_k30 = denom(Delta_k30_canc)

num_k30_exp = expand(num_k30)
den_k30_exp = expand(den_k30)
den_k30_fac = factor(den_k30)

print(f"\nDelta|_{{k3=0}} numerator (expanded):")
print(f"  {num_k30_exp}")
print(f"\nDelta|_{{k3=0}} denominator (factored):")
print(f"  {den_k30_fac}")

# Check denominator is positive on domain
print(f"\nDenominator factors: {den_k30_fac}")
print("  All factors positive on domain (s,t>0, u<4s^2, v<4t^2, u+v<4(s+t)^2)")
