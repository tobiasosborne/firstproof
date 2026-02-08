"""
PROVER-11 Part 10: Direct generalized CS approach for full R_4.

From the partial fraction decomposition (verified):
  -R4 = W/(6*D1) + K3^2*W/(6*K2^2*D2) - K4/(24*K2) - K2/6 - K3^2/(12*K2^2)

where W = 4*K2^3 - K3^2, D1 = 4*K2^2 - K4, D2 = 4*K2^3 + K2*K4 - 2*K3^2.

Regroup:
  -R4 = W/(6*D1) + K3^2*W/(6*K2^2*D2) + K2/6 + K3^2/(12*K2^2) + K4/(24*K2)

Wait, signs. Let me be careful.
  -R4 = +W/(6*D1) + K3^2*W/(6*K2^2*D2) - K4/(24*K2) - K2/6 - K3^2/(12*K2^2)

Subadditivity of -R4 means: -R4(r) <= -R4(p) + -R4(q)
i.e., each term at r <= sum of that term at p and q.

Terms:
(a) W/(6*D1): need Wr/(6*D1r) <= Wp/(6*D1p) + Wq/(6*D1q)
(b) K3^2*W/(6*K2^2*D2): need K3r^2*Wr/(6*K2r^2*D2r) <= same for p + q
(c) -K4/(24*K2): need -K4r/(24*K2r) <= -K4p/(24*K2p) + -K4q/(24*K2q)
    i.e., K4r/K2r >= K4p/K2p + K4q/K2q. NO! K4r/K2r = (K4p+K4q)/(K2p+K2q)
    which is between K4p/K2p and K4q/K2q, not their sum. So this is WRONG.
(d) -K2/6: linear, additive, gap = 0 ✓
(e) -K3^2/(12*K2^2): need K3r^2/K2r^2 <= K3p^2/K2p^2 + K3q^2/K2q^2
    This is the n=3 result. ✓

So term (c) is NOT individually subadditive. We can't decompose term by term.

Let me try a DIFFERENT approach entirely: the direct Schur-convexity method.

APPROACH E: Lagrange multiplier / parametric family
For the full gap, introduce w = k4/k2^2 (scaled k4) and mu = k3/k2^{3/2} (scaled k3).
Then R4 = k2 * f(mu, w) where f is the reduced function.

Superadditivity becomes:
  (s+t)*f(mu_r, w_r) >= s*f(mu_p, w_p) + t*f(mu_q, w_q)

where mu_r = (k3p+k3q)/(s+t)^{3/2} and w_r = (k4p+k4q)/(s+t)^2.

Now set k3p = s^{3/2}*mu_p, k3q = t^{3/2}*mu_q, k4p = s^2*w_p, k4q = t^2*w_q.
Then:
  mu_r = (s^{3/2}*mu_p + t^{3/2}*mu_q) / (s+t)^{3/2}
  w_r = (s^2*w_p + t^2*w_q) / (s+t)^2

The superadditivity condition becomes:
  (s+t)*f(mu_r, w_r) - s*f(mu_p, w_p) - t*f(mu_q, w_q) >= 0

Let me try WLOG s+t = 1 (by homogeneity). Set s = alpha, t = 1-alpha.
Then the gap = f(alpha^{3/2}*a + (1-alpha)^{3/2}*b, alpha^2*c + (1-alpha)^2*d)
              - alpha*f(a, c) - (1-alpha)*f(b, d)

where a=mu_p, b=mu_q, c=w_p, d=w_q.

This is hard to analyze in general.

APPROACH F: Reduction to POLYNOMIAL inequality
The gap numerator is a polynomial in (s, t, k3p, k3q, k4p, k4q).
The denominator is a product of positive terms.
So we need: gap numerator >= 0.

Let me compute the gap numerator symbolically for small cases.
"""
from sympy import (symbols, expand, factor, collect, cancel,
                   Rational, simplify, S, Matrix, det, Poly, numer, denom)

s, t = symbols('s t', positive=True)
k3p, k3q = symbols('k3p k3q', real=True)
k4p, k4q = symbols('k4p k4q', real=True)

def R4(K2, K3, K4):
    num = -16*K2**3*K3**2 - 4*K2**2*K4**2 + 20*K2*K3**2*K4 - 8*K3**4 - K4**3
    den = 24*(4*K2**2 - K4)*(4*K2**3 + K2*K4 - 2*K3**2)
    return num, den

# Delta = R4_r - R4_p - R4_q
# R4_i = P_i / Q_i
# Delta = P_r/Q_r - P_p/Q_p - P_q/Q_q
# = (P_r*Q_p*Q_q - P_p*Q_q*Q_r - P_q*Q_p*Q_r) / (Q_r*Q_p*Q_q)

Pr, Qr = R4(s+t, k3p+k3q, k4p+k4q)
Pp, Qp = R4(s, k3p, k4p)
Pq, Qq = R4(t, k3q, k4q)

# The numerator of Delta:
# GapNum = Pr*Qp*Qq - Pp*Qq*Qr - Pq*Qp*Qr

# This is huge. Let me try the k4=0 slice first.
print("SLICE: k4p = k4q = 0")
print("="*70)

Pr0, Qr0 = R4(s+t, k3p+k3q, 0)
Pp0, Qp0 = R4(s, k3p, 0)
Pq0, Qq0 = R4(t, k3q, 0)

# Simplify
R4r0 = cancel(Pr0/Qr0)
R4p0 = cancel(Pp0/Qp0)
R4q0 = cancel(Pq0/Qq0)

print(f"R4(s, k3p, 0) = {R4p0}")
print(f"R4(t, k3q, 0) = {R4q0}")

Delta_k40 = cancel(R4r0 - R4p0 - R4q0)
num_k40 = numer(Delta_k40)
den_k40 = denom(Delta_k40)
print(f"\nDelta|_{{k4=0}} num = {expand(num_k40)}")
print(f"  factored = {factor(num_k40)}")
print(f"\nDelta|_{{k4=0}} den factored = {factor(den_k40)}")

# At k4=0: R4(K2, K3, 0) = (-16K2^3*K3^2 - 8K3^4)/(24*4K2^2*(4K2^3-2K3^2))
# = -8K3^2*(2K2^3+K3^2) / (192*K2^2*(2K2^3-K3^2))
# = -K3^2*(2K2^3+K3^2) / (24*K2^2*(2K2^3-K3^2))

# Hmm complex. Let me try ANOTHER specific slice: k3q = k4p = k4q = 0
# Then: Delta = R4(s+t, k3p, 0) - R4(s, k3p, 0) - R4(t, 0, 0)
# R4(t, 0, 0) = 0. So Delta = R4(s+t, k3p, 0) - R4(s, k3p, 0)

print("\n\nSLICE: k3q = k4p = k4q = 0")
Delta_simple = cancel(R4r0.subs({k3q:0, k4p:0, k4q:0}) - R4p0.subs({k4p:0}))
# R4(t,0,0) = 0
print(f"Delta = R4(s+t, k3p, 0) - R4(s, k3p, 0) - 0")
print(f"  = {cancel(Delta_simple)}")
num_simple = numer(Delta_simple)
print(f"  num factored = {factor(num_simple)}")
print(f"  den factored = {factor(denom(Delta_simple))}")

# Actually, let me think about this differently.
# In the n=3 case, 1/Phi_3 = (2/9)*k2 - (2/27)*k3^2/k2^2
# The term k3^2/k2^2 is subadditive by Cauchy-Schwarz.
#
# For n=4, what if we write:
# 1/Phi_4 = c1*k2 + c2*k3^2/k2^2 + c3*k4/k2 + c4*k4^2/(k2*D1) + c5*k3^2*k4/(k2^2*D2) + ...
# and try to find a representation where each term is super/sub-additive?

# Let me try yet another approach: the SUBSTITUTION approach.
# Set k3 = k2^{3/2}*mu, k4 = k2^2*w, and work with 1/Phi_4 as k2 * F(mu, w).
# Then superadditivity of 1/Phi_4 = (s+t)*F(mu_r, w_r) >= s*F(mu_p, w_p) + t*F(mu_q, w_q).
#
# This fails if F is not concave (which we showed).
# But maybe F can be decomposed as F = F_concave + F_linear?

# Actually, 1/Phi_4 = k2/12 + R4, and k2/12 is linear.
# So 1/Phi_4 superadditive <==> R4 superadditive.
# And R4 = k2 * f(mu, w) where f = p/(24*(4-w)*(w+4-2mu^2)).
# Hessian of f is indefinite, so f is NOT concave.

# APPROACH G: Direct polynomial certificate via SYMMETRY.
# The gap is symmetric under (s,k3p,k4p) <-> (t,k3q,k4q).
# Perhaps use this symmetry + AM-GM to find a certificate.

# Let me compute the gap for a one-parameter family and try to understand
# the structure.
print("\n\n" + "="*70)
print("ONE-PARAMETER FAMILY: s=t, k3p=-k3q (antisymmetric k3)")
print("="*70)

# With s=t and k3p = -k3q = m: k3r = 0
# k4p = k4q = w
# R4(2s, 0, 2w) vs 2*R4(s, m, w)
K = symbols('K', positive=True)
m, w = symbols('m w', real=True)

R4r_sym = cancel(R4(2*K, 0, 2*w)[0] / R4(2*K, 0, 2*w)[1])
R4p_sym = cancel(R4(K, m, w)[0] / R4(K, m, w)[1])

Delta_sym = cancel(R4r_sym - 2*R4p_sym)
print(f"Delta (s=t=K, k3p=-k3q=m, k4p=k4q=w):")
print(f"  = {Delta_sym}")
num_sym = expand(numer(Delta_sym))
den_sym = factor(denom(Delta_sym))
print(f"  num = {num_sym}")
print(f"  num factored = {factor(num_sym)}")
print(f"  den factored = {den_sym}")
