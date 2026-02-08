"""
PROVER-9: Perspective function / Schur-convexity approach.

For k3=0: R_4 = -k4^2 / (24*k2*(4*k2^2 - k4))

Let w = k4/k2^2 (normalized cumulant). Then R_4 = k2 * g(w) where
g(w) = -w^2/(24*(4-w)).

The perspective of g: P(k2, k4) = k2 * g(k4/k2^2)

For SUPERADDITIVITY of P:
P(s+t, u+v) >= P(s,u) + P(t,v)

This is NOT the standard perspective function (which would be t*g(x/t)).
The key difference: the argument to g is k4/k2^2 (quadratic denominator).

Let me try a substitution. Let m = k2, x = k4.
P(m, x) = -x^2 / (24*m*(4*m^2 - x))

Define Q(m, x) = -P(m,x) = x^2 / (24*m*(4*m^2 - x))

We need Q SUBADDITIVE: Q(s+t, u+v) <= Q(s,u) + Q(t,v)

Idea: Write Q(m,x) = (x^2/(4m^2-x)) * (1/(24m))
     = (x/(4m^2-x)) * (x/(24m))

Let's try Q = x * h(m,x) where h = x/(24*m*(4*m^2-x))

Another idea: Q(m,x) = x^2 / (24*(4*m^3 - m*x))
Note 4*m^3 - m*x = m*(4*m^2-x), so the denominator is a degree-3 polynomial in m.

CRITICAL INSIGHT: In the n=3 case, R_3 = -2*k3^2/(27*k2^2).
The proof used that x -> x^2/y^2 is jointly convex for y > 0.
More precisely: (a+c)^2/(b+d)^2 <= a^2/b^2 + c^2/d^2 when b,d > 0.

This is the "power mean inequality" or "Cauchy-Schwarz in Engel/Titu form".

For n=4, k3=0: we need (u+v)^2/((s+t)*(4*(s+t)^2-(u+v))) <= u^2/(s*(4*s^2-u)) + v^2/(t*(4*t^2-v))

This has the form: f^2/phi <= f_p^2/phi_p + f_q^2/phi_q
where phi(m,x) = m*(4*m^2-x) = 4*m^3 - m*x.

If phi were LINEAR and positive, this would follow from Cauchy-Schwarz (Titu/Engel).
phi is NOT linear.

But: phi = m*(4*m^2-x) is "jointly concave" in (m,x) for m > 0, x < 4*m^2?
Let's check...
"""
from sympy import symbols, diff, Matrix, det, simplify, factor, expand, Rational, cancel
import numpy as np

m, x = symbols('m x', real=True)

phi = m*(4*m**2 - x)
print("phi(m, x) = m*(4*m^2 - x) = 4*m^3 - m*x")

# Hessian of phi
phi_mm = diff(phi, m, 2)
phi_mx = diff(phi, m, x)
phi_xx = diff(phi, x, 2)
print(f"phi_mm = {phi_mm}")   # = 24*m
print(f"phi_mx = {phi_mx}")   # = -1
print(f"phi_xx = {phi_xx}")   # = 0
# Hessian = [[24m, -1], [-1, 0]], det = 0 - 1 = -1 < 0
# phi is NEITHER convex nor concave.
print("phi is neither convex nor concave (indefinite Hessian).\n")

# ======================================================================
# NEW APPROACH: Direct algebraic proof for k3=0 using Schur complement
# ======================================================================

print("="*70)
print("APPROACH: Prove numerator of Delta|_{k3=0} is non-negative")
print("="*70)

# The numerator (from sympy) with s=k2p, t=k2q, u=k4p, v=k4q:
# Factor common terms from the sympy factored denominator:
# Den = 24*s*t*(s+t)*(4*s^2-u)*(4*t^2-v)*(4*(s+t)^2-u-v)
# This is positive on domain. So need Num >= 0.

# Let me write the numerator as a polynomial in u, v (for fixed s, t > 0).
# From the output:
# Num = 16*s^6*v^2 + 48*s^5*t*v^2 + 48*s^4*t^2*v^2 - 8*s^4*u*v^2 - 4*s^4*v^3
#     - 32*s^3*t^3*u*v - 8*s^3*t*u*v^2 + 48*s^2*t^4*u^2 - 12*s^2*t^2*u^2*v
#     - 12*s^2*t^2*u*v^2 + s^2*u^2*v^2 + s^2*u*v^3 + 48*s*t^5*u^2
#     - 8*s*t^3*u^2*v + 16*t^6*u^2 - 4*t^4*u^3 - 8*t^4*u^2*v + t^2*u^3*v + t^2*u^2*v^2

# Group as quadratic form in (u^2, uv, v^2, u, v) ... no, it's degree 3 in (u,v).

# The numerator has degree 3 in (u,v). Let me write it as:
# Num = A(s,t)*v^3 + B(s,t)*u^3 + ... (degree 3 terms)
#     + C(s,t)*u^2*v^2 + D(s,t)*u^2*v + E(s,t)*u*v^2 + ... (degree 2 cross terms)
#     + ... (degree 2 same variable terms)

s, t, u, v = symbols('s t u v', real=True)

Num = (16*s**6*v**2 + 48*s**5*t*v**2 + 48*s**4*t**2*v**2 - 8*s**4*u*v**2 - 4*s**4*v**3
    - 32*s**3*t**3*u*v - 8*s**3*t*u*v**2 + 48*s**2*t**4*u**2 - 12*s**2*t**2*u**2*v
    - 12*s**2*t**2*u*v**2 + s**2*u**2*v**2 + s**2*u*v**3 + 48*s*t**5*u**2
    - 8*s*t**3*u**2*v + 16*t**6*u**2 - 4*t**4*u**3 - 8*t**4*u**2*v + t**2*u**3*v + t**2*u**2*v**2)

# Check symmetry: does Num(s,t,u,v) = Num(t,s,v,u)?
Num_swapped = Num.subs({s: t, t: s, u: v, v: u})
print(f"\nNum symmetric under (s,t,u,v) <-> (t,s,v,u)? {simplify(Num - Num_swapped) == 0}")

# Good! This symmetry is expected since R_4(k_p) and R_4(k_q) are interchangeable.

# Let me try writing Num as a sum of terms, each non-negative.
# Key idea: group terms that form "squares" or "products of positives"

# On domain: u < 4*s^2, v < 4*t^2, u+v < 4*(s+t)^2
# So 4*s^2 - u > 0, 4*t^2 - v > 0

# Let a = 4*s^2 - u > 0, b = 4*t^2 - v > 0
# Then u = 4*s^2 - a, v = 4*t^2 - b

a, b = symbols('a b', positive=True)
Num_ab = expand(Num.subs({u: 4*s**2 - a, v: 4*t**2 - b}))
print(f"\nNum in (a,b) where a=4s^2-u, b=4t^2-v:")

# Collect by powers of a, b
from sympy import collect, Poly
Num_ab_poly = Poly(Num_ab, a, b)
print("\nMonomials in (a, b):")
for monom, coeff in sorted(Num_ab_poly.as_dict().items()):
    coeff_f = factor(coeff)
    print(f"  a^{monom[0]} * b^{monom[1]}: {coeff_f}")

# Check: when a=0, b=0 (i.e., u=4s^2, v=4t^2), what happens?
# This is the boundary of the domain.
Num_boundary = Num.subs({u: 4*s**2, v: 4*t**2})
print(f"\nNum at boundary (u=4s^2, v=4t^2): {expand(Num_boundary)}")
print(f"Factored: {factor(Num_boundary)}")

# What about when u=0, v=0?
Num_origin = Num.subs({u: 0, v: 0})
print(f"\nNum at origin (u=v=0): {expand(Num_origin)}")

# At u=v=0, R4_p = R4_q = R4_r = 0, so gap = 0. Num should be 0. ✓

# ======================================================================
# TRY: Write Num as sum of non-negative terms on the domain
# ======================================================================
print("\n\n" + "="*70)
print("SOS-LIKE DECOMPOSITION ATTEMPT")
print("="*70)

# Idea: factor Num = u^2 * P1(s,t,u,v) + v^2 * P2(s,t,u,v) + u*v * P3(s,t,u,v)
# Since Num vanishes at u=v=0, we can write Num = u^2*A + v^2*B + u*v*C
# Check degrees: Num has degree 3 in (u,v).
# So: A has degree 1, B has degree 1, C has degree 1 in (u,v).

# Actually, let's collect by u, v more carefully
Num_poly_uv = Poly(Num, u, v)
print("Num as polynomial in (u, v):")
for monom, coeff in sorted(Num_poly_uv.as_dict().items()):
    coeff_s = factor(coeff)
    print(f"  u^{monom[0]} * v^{monom[1]}: {coeff_s}")

print("\n\nWriting Num = v^2*(alpha + beta*u + gamma*v) + u^2*(delta + eps*v + zeta*u)")
print("  + u*v*(eta)  ... no, this doesn't separate cleanly.")

# Let me try a different grouping using the domain constraint
# On domain: u < 4s^2 so (4s^2-u) > 0, similarly (4t^2-v) > 0

# Factor attempt: Maybe Num = (4s^2-u)*A + (4t^2-v)*B + C where A,B,C >= 0?

# Or maybe: write as a quadratic form in the vector (u, v)
# Num = [u, v] * M * [u, v]^T + cubic terms
# where M has entries depending on s, t

# The quadratic (in u,v) part:
Q_part = (48*s**2*t**4 * u**2 + (-32*s**3*t**3) * u*v + 16*s**6 * v**2
         + 48*s**5*t * v**2 + 48*s**4*t**2 * v**2
         + 48*s*t**5 * u**2 + 16*t**6 * u**2)

Q_coeff_u2 = 48*s**2*t**4 + 48*s*t**5 + 16*t**6
Q_coeff_uv = -32*s**3*t**3
Q_coeff_v2 = 16*s**6 + 48*s**5*t + 48*s**4*t**2

print(f"\nQuadratic part in (u,v):")
print(f"  Coeff of u^2: {factor(Q_coeff_u2)} = {factor(Q_coeff_u2)}")
print(f"  Coeff of uv: {Q_coeff_uv}")
print(f"  Coeff of v^2: {factor(Q_coeff_v2)} = {factor(Q_coeff_v2)}")

# Q matrix
from sympy import Matrix
Q_mat = Matrix([[Q_coeff_u2, Q_coeff_uv/2], [Q_coeff_uv/2, Q_coeff_v2]])
Q_det = simplify(det(Q_mat))
print(f"\nQ matrix determinant: {factor(Q_det)}")
Q_tr = Q_mat[0,0] + Q_mat[1,1]
print(f"Q matrix trace: {factor(Q_tr)}")

# If det > 0 and trace > 0, the quadratic form is positive definite
# det = (48*s^2*t^4+48*s*t^5+16*t^6)*(16*s^6+48*s^5*t+48*s^4*t^2) - (16*s^3*t^3)^2
# = 16*t^4*(3s^2+3st+t^2)*16*s^4*(s^2+3st+3t^2) - 256*s^6*t^6
# = 256*s^4*t^4*(3s^2+3st+t^2)*(s^2+3st+3t^2) - 256*s^6*t^6
# = 256*s^4*t^4*[(3s^2+3st+t^2)(s^2+3st+3t^2) - s^2*t^2]

product = expand((3*s**2+3*s*t+t**2)*(s**2+3*s*t+3*t**2))
print(f"\n(3s^2+3st+t^2)(s^2+3st+3t^2) = {product}")
print(f"Minus s^2*t^2: {expand(product - s**2*t**2)}")
print(f"Factored: {factor(product - s**2*t**2)}")

print("\n\n" + "="*70)
print("CHECKING IF THE QUADRATIC FORM DOMINATES THE CUBIC REMAINDER")
print("="*70)

# The remaining (cubic in u,v) terms:
Cubic = Num - Q_part.subs(Q_coeff_u2, Q_coeff_u2).subs(Q_coeff_uv, Q_coeff_uv).subs(Q_coeff_v2, Q_coeff_v2)
# Actually, let me just compute it directly
Remainder = expand(Num - (Q_coeff_u2*u**2 + Q_coeff_uv*u*v + Q_coeff_v2*v**2))
print(f"\nRemainder R = Num - Q(u,v):")
Rem_poly = Poly(Remainder, u, v)
for monom, coeff in sorted(Rem_poly.as_dict().items()):
    print(f"  u^{monom[0]} * v^{monom[1]}: {factor(coeff)}")
