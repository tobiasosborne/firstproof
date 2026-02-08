"""
PROVER-11 Part 9: Analytical proof of full R_4 superadditivity.

Key decomposition:
  -R4 = (4K2^3-K3^2)/(6*D1) + K3^2*(4K2^3-K3^2)/(6*K2^2*D2) - K4/(24*K2) - K2/6 - K3^2/(12*K2^2)

Rewrite setting W = 4K2^3 - K3^2 (positive on domain, proved in Part 8):
  -R4 = W/(6*D1) + K3^2*W/(6*K2^2*D2) - K4/(24*K2) - K2/6 - K3^2/(12*K2^2)

Group: -R4 = W*(1/D1 + K3^2/(K2^2*D2))/6 - K4/(24*K2) - K2/6 - K3^2/(12*K2^2)

Hmm, this is getting complicated. Let me try a completely different approach.

APPROACH: Write 1/Phi_4 directly and prove its superadditivity.
Since R4 = 1/Phi_4 - K2/12 and K2/12 is linear (hence additive),
superadditivity of R4 <==> superadditivity of 1/Phi_4.

1/Phi_4 = N4_disc / (N4_poly)  where these are specific expressions.

But actually, let me try the SIMPLEST possible approach:
show that -R4 = h(K3, K4; K2) where for fixed K2, h is a convex
function of (K3, K4), and then use the subadditivity lemma:
if K2 -> h(K3, K4; K2) has the right structure...

Actually, the cleanest approach might be: DIRECT Cauchy-Schwarz
on the representation -R4 = X^T M X / D where X is a vector of
monomials and M is positive definite.

Let me try to write -R4 * 24 * D1 * D2 as a quadratic form.
"""
from sympy import (symbols, expand, factor, collect, cancel,
                   Rational, simplify, S, Matrix, det, Poly, numer, denom)
import numpy as np

K2, K3, K4 = symbols('K2 K3 K4')
s, t = symbols('s t', positive=True)
k3p, k3q, k4p, k4q = symbols('k3p k3q k4p k4q', real=True)

# Let me try the approach of writing R4 in a form amenable to
# the Cauchy-Schwarz / Titu's lemma generalization.
#
# Recall for n=3: -R3 = (2/27)*K3^2/K2^2.
# The subadditivity proof was: K3^2/K2^2 is convex in (K3, K2)
# (by Cauchy-Schwarz / homogeneity argument).
#
# For n=4, -R4 has more complex structure.
# Let me try: write the numerator of -R4 as a sum of terms,
# each of which gives a fraction that is subadditive.

# -R4 = [16K2^3*K3^2 + 4K2^2*K4^2 - 20K2*K3^2*K4 + 8K3^4 + K4^3]
#        / [24*(4K2^2-K4)*(4K2^3+K2*K4-2K3^2)]

# Numerator N of -R4:
N = 16*K2**3*K3**2 + 4*K2**2*K4**2 - 20*K2*K3**2*K4 + 8*K3**4 + K4**3

# Let me try to express N as sum of squares or products of positive things.
# N = 8K3^4 + K3^2*(16K2^3-20K2*K4) + K4^2*(4K2^2+K4)
# Viewing as quadratic in K3^2:
# N = 8*(K3^2)^2 + (16K2^3-20K2*K4)*K3^2 + K4^2*(4K2^2+K4)
# discriminant = (16K2^3-20K2*K4)^2 - 32*K4^2*(4K2^2+K4)
#              = 16*(K2^2-2K4)*(4K2^2-K4)^2  (from Part 7!)

# When K2^2 >= 2K4: disc >= 0, N has real roots in K3^2.
# When K2^2 < 2K4: disc < 0, N > 0 (positive definite in K3^2).
# On domain K4 < 4K2^2, so K4/K2^2 < 4. Both K2^2 >= 2K4 and K2^2 < 2K4 occur.

# Let me try: complete the square in K3^2.
# N = 8*(K3^2 + (16K2^3-20K2*K4)/16)^2 - (16K2^3-20K2*K4)^2/32 + K4^2*(4K2^2+K4)
# = 8*(K3^2 + K2*(4K2^2-5K4)/4)^2 + K4^2*(4K2^2+K4) - K2^2*(4K2^2-5K4)^2/2

# Remainder: K4^2*(4K2^2+K4) - K2^2*(4K2^2-5K4)^2/2
rem = expand(K4**2*(4*K2**2+K4) - K2**2*(4*K2**2-5*K4)**2/2)
print(f"Remainder = {rem}")
print(f"  factored = {factor(rem)}")

# Let me try a direct matrix approach for the full gap.
# The gap Delta = R4_r - R4_p - R4_q
# = P_r/(24*D1r*D2r) - P_p/(24*D1p*D2p) - P_q/(24*D1q*D2q)
# Common denominator: 24*D1p*D2p*D1q*D2q*D1r*D2r
# Numerator: P_r*D1p*D2p*D1q*D2q - P_p*D1q*D2q*D1r*D2r - P_q*D1p*D2p*D1r*D2r

# This is a massive polynomial. Let me try a smarter route.

# CLEANEST APPROACH: Epigraph / perspective-like
#
# Note that R4 has cumulant weight 2. So R4(lambda*K) = lambda^2 * R4(K)
# for all lambda > 0 (where K = (K2, K3, K4) scaled by cumulant weights).
# Wait: K2 has weight 2, K3 weight 3, K4 weight 4.
# R4(c*K2, c^{3/2}*K3, c^2*K4) = c * R4(K2, K3, K4)... let me check.
# Actually R4 has weight 2 in the cumulant weighting, meaning:
# R4(c*K2, c^{3/2}*K3, c^2*K4) = c * R4(K2, K3, K4)
# (since num has weight 12, den has weight 10, ratio has weight 2,
#  and scaling by c gives c^{12}/c^{10} = c^2, then dividing by c^1 from K2 unit... hmm)

# Actually from the formula R4 = K2 * f(u,v) with u = K3/K2^{3/2}, v = K4/K2^2:
# R4 is degree 1 in K2 (after homogeneous rescaling). ✓

# So the superadditivity gap is:
# (s+t)*f(u_r, v_r) - s*f(u_p, v_p) - t*f(u_q, v_q)
# where u_r, v_r are the scaled cumulants of the sum.

# KEY INSIGHT: This has the structure of a GENERALIZED PERSPECTIVE function!
# g(K2, K3, K4) = K2 * f(K3/K2^{3/2}, K4/K2^2)
# If f is CONCAVE (as a function of (u,v)), then g is superadditive.
# This is because g is the "perspective" of f with weight K2.

# Wait, is this true? For standard perspective: if f is concave, then
# t*f(x/t) is concave (jointly in (x,t)). This gives superadditivity:
# (s+t)*f((x_p+x_q)/(s+t)) >= s*f(x_p/s) + t*f(x_q/t) [by concavity of perspective]
# But our scaling is NOT linear: u = K3/K2^{3/2}, not K3/K2.

# Let me check: is f(u,v) concave on the domain?
print("\n" + "="*70)
print("CHECKING CONCAVITY OF f(u,v)")
print("="*70)

u, v = symbols('u v', real=True)
p_uv = -8*u**4 + 20*u**2*v - 16*u**2 - v**3 - 4*v**2
q_uv = 24*(4 - v)*(v + 4 - 2*u**2)

f_uv = p_uv / q_uv

f_uu = simplify(f_uv.diff(u, 2))
f_uv_mixed = simplify(f_uv.diff(u).diff(v))
f_vv = simplify(f_uv.diff(v, 2))

# Evaluate Hessian at a few points
import numpy as np
from sympy import N as Neval

test_points = [(0, 0), (0, 1), (0, -1), (0.5, 0), (0.5, 1), (1, 0)]
print("\nHessian eigenvalues at sample points (u, v):")
for u_val, v_val in test_points:
    # Check domain: v < 4 and 2u^2 < v+4
    if v_val >= 4 or 2*u_val**2 >= v_val + 4:
        print(f"  ({u_val}, {v_val}): OUTSIDE DOMAIN")
        continue
    h11 = float(f_uu.subs({u: u_val, v: v_val}))
    h12 = float(f_uv_mixed.subs({u: u_val, v: v_val}))
    h22 = float(f_vv.subs({u: u_val, v: v_val}))
    H = np.array([[h11, h12], [h12, h22]])
    eigs = np.linalg.eigvalsh(H)
    is_neg_semidef = all(e <= 1e-10 for e in eigs)
    print(f"  ({u_val}, {v_val}): eigs = [{eigs[0]:.6f}, {eigs[1]:.6f}] {'CONCAVE' if is_neg_semidef else 'NOT CONCAVE'}")
