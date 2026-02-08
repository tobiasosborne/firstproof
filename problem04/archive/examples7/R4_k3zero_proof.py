"""
PROVER-9: Attempt to prove R_4|_{k3=0} superadditivity.

We established:
- Delta|_{k3=0} = Num / Den
- Den > 0 on domain
- Need: Num >= 0

Num = 16*s^6*v^2 + 48*s^5*t*v^2 + 48*s^4*t^2*v^2 - 8*s^4*u*v^2 - 4*s^4*v^3
    - 32*s^3*t^3*u*v - 8*s^3*t*u*v^2 + 48*s^2*t^4*u^2 - 12*s^2*t^2*u^2*v
    - 12*s^2*t^2*u*v^2 + s^2*u^2*v^2 + s^2*u*v^3 + 48*s*t^5*u^2
    - 8*s*t^3*u^2*v + 16*t^6*u^2 - 4*t^4*u^3 - 8*t^4*u^2*v + t^2*u^3*v + t^2*u^2*v^2

where s=k2_p, t=k2_q, u=k4_p, v=k4_q

Domain: s,t > 0, u < 4*s^2, v < 4*t^2, u+v < 4*(s+t)^2

Let's try substituting u = w_p*s^2, v = w_q*t^2 where -4 < w_p,w_q < 4.
"""
from sympy import (symbols, Rational, simplify, factor, expand, cancel,
                   Matrix, det, S, Poly, together, collect, sqrt, groebner,
                   resultant, solve, diff, degree, eye, numer, denom,
                   Symbol)
import sympy

s, t = symbols('s t', positive=True)
wp, wq = symbols('wp wq', real=True)

# Substitute u = wp*s^2, v = wq*t^2
u_expr = wp * s**2
v_expr = wq * t**2

# Build the numerator
u, v = symbols('u v', real=True)
Num = (16*s**6*v**2 + 48*s**5*t*v**2 + 48*s**4*t**2*v**2 - 8*s**4*u*v**2 - 4*s**4*v**3
    - 32*s**3*t**3*u*v - 8*s**3*t*u*v**2 + 48*s**2*t**4*u**2 - 12*s**2*t**2*u**2*v
    - 12*s**2*t**2*u*v**2 + s**2*u**2*v**2 + s**2*u*v**3 + 48*s*t**5*u**2
    - 8*s*t**3*u**2*v + 16*t**6*u**2 - 4*t**4*u**3 - 8*t**4*u**2*v + t**2*u**3*v + t**2*u**2*v**2)

Num_subst = expand(Num.subs({u: u_expr, v: v_expr}))
print("Num after substituting u=wp*s^2, v=wq*t^2:")

# Factor out s^4*t^4 (if possible)
# The minimum degree term: for smallest total degree in s,t
# 16*s^6*v^2 = 16*s^6*wq^2*t^4 → s^6*t^4
# s^2*u*v^3 = s^2*(wp*s^2)*(wq*t^2)^3 = wp*wq^3*s^4*t^6 → s^4*t^6
# All terms should be at least s^4*t^4

# Let's extract the common factor
Num_reduced = simplify(Num_subst / (s**4 * t**4))
print(f"\nNum = s^4*t^4 * F(s,t,wp,wq)")
print(f"F = {expand(Num_reduced)}")

# Now let r = s/t (ratio), so s = r*t
r = symbols('r', positive=True)
F_in_r = expand(Num_reduced.subs(s, r*t))
# Factor out powers of t
F_final = simplify(F_in_r / t**0)  # should be in r, wp, wq only
print(f"\nF(r, wp, wq) = {expand(F_final)}")

# Actually let's just substitute s=r, t=1 (WLOG by scaling)
F_r = expand(Num_reduced.subs(t, 1))
print(f"\nWith t=1, s=r: F(r, wp, wq) = {expand(F_r)}")

# Let's collect by powers of wp and wq
F_r_collected_wq = collect(expand(F_r), wq)
print(f"\nCollected by wq: {F_r_collected_wq}")

print("\n" + "="*70)
print("APPROACH: Schur-convexity / quadratic form in (wp, wq)")
print("="*70)

# F is a polynomial in wp, wq of degree at most 3
# Let's write F as a function of wp, wq for fixed r=s/t > 0

# Actually F should be at most degree 3 in each of wp, wq
# Let's identify the structure

# Group by total degree in (wp, wq):
# degree 0 in (wp,wq): none (when wp=wq=0, numerator should be 0)
# degree 1: should also be 0
# degree 2: quadratic form
# degree 3: cubic terms

# When wp=0 and wq=0: u=0, v=0, so R4_p = R4_q = 0 and R4_r = 0. gap=0. ✓

# Let's try: is F always = sum of squares times positive things?

# First, verify F as a polynomial in (wp, wq) for symbolic r
from sympy import Poly as SPoly

F_r_sym = F_r.subs(s, r)
F_poly = SPoly(F_r_sym, wp, wq)
print(f"\nF as polynomial in (wp, wq):")
print(f"  Monomials and coefficients:")
for monom, coeff in F_poly.as_dict().items():
    print(f"    wp^{monom[0]} * wq^{monom[1]}: {coeff}")

print("\n" + "="*70)
print("APPROACH: Factor the k3=0 numerator as a positive expression")
print("="*70)

# Let me try a different approach. Go back to the original variables.
# Num(s,t,u,v) where u < 4s^2, v < 4t^2
# Write u = 4s^2 - a, v = 4t^2 - b where a, b > 0

a_sym, b_sym = symbols('a b', positive=True)
u_new = 4*s**2 - a_sym
v_new = 4*t**2 - b_sym

Num_new = expand(Num.subs({u: u_new, v: v_new}))
print(f"\nNum with u=4s^2-a, v=4t^2-b:")
Num_new_collected = collect(expand(Num_new), [a_sym, b_sym])
# This might be cleaner
print(f"  = {Num_new_collected}")

# Factor
Num_new_f = factor(Num_new)
print(f"\nFactored: {Num_new_f}")
