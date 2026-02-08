"""
PROVER-12 Step 1: Compute the gap numerator symbolically and analyze structure.

Goal: Δ = R_4(p+q) - R_4(p) - R_4(q) >= 0

Strategy: clear denominators, get gap numerator as a polynomial in 6 variables,
then try to show it's non-negative on the domain.
"""
from sympy import (symbols, expand, factor, cancel, simplify, collect,
                   numer, denom, Rational, S, Poly, groebner, resultant,
                   sqrt, together, apart, fraction)
import sympy

s, t = symbols('s t', positive=True)       # k2_p, k2_q
a, b = symbols('a b', real=True)            # k3_p, k3_q
u, v = symbols('u v', real=True)            # k4_p, k4_q

def R4_sym(K2, K3, K4):
    """R_4 as symbolic expression"""
    num = -16*K2**3*K3**2 - 4*K2**2*K4**2 + 20*K2*K3**2*K4 - 8*K3**4 - K4**3
    D1 = 4*K2**2 - K4
    D2 = 4*K2**3 + K2*K4 - 2*K3**2
    den = 24*D1*D2
    return num, den

print("="*70)
print("STEP 1: Compute gap numerator")
print("="*70)

# R_4 at sum
num_r, den_r = R4_sym(s+t, a+b, u+v)
# R_4 at p
num_p, den_p = R4_sym(s, a, u)
# R_4 at q
num_q, den_q = R4_sym(t, b, v)

# Gap = num_r/den_r - num_p/den_p - num_q/den_q
# = (num_r * den_p * den_q - num_p * den_r * den_q - num_q * den_r * den_p) / (den_r * den_p * den_q)

# The common denominator is positive on the domain, so we need:
# gap_num = num_r * den_p * den_q - num_p * den_r * den_q - num_q * den_r * den_p >= 0

print("Computing gap numerator (this may take a moment)...")
gap_num_raw = expand(num_r * den_p * den_q - num_p * den_r * den_q - num_q * den_r * den_p)
print(f"Gap numerator has {len(gap_num_raw.as_ordered_terms())} terms")

# Let's also look at the gap using the partial fraction decomposition.
# From PROVER-11:
# -R_4 = (4*K2^3 - K3^2)/(6*D1) + K3^2*(4*K2^3 - K3^2)/(6*K2^2*D2) - K4/(24*K2) - (2*K2^3 + K3^2)/(12*K2^2)

# The last term -(2*K2^3 + K3^2)/(12*K2^2) = -K2/6 - K3^2/(12*K2^2)
# -K2/6 is linear hence additive (contributes 0 to gap)
# -K3^2/(12*K2^2) is the same form as R_3 (proved subadditive)

# So we can write:
# -R_4 = Term1 + Term2 + Term3 + Term4 + Term5
# where:
#   Term1 = (4K2^3-K3^2)/(6*D1)           [depends on K2,K3,K4]
#   Term2 = K3^2*(4K2^3-K3^2)/(6*K2^2*D2) [depends on K2,K3,K4]
#   Term3 = -K4/(24*K2)                    [depends on K2,K4]
#   Term4 = -K2/6                          [linear, additive]
#   Term5 = -K3^2/(12*K2^2)               [depends on K2,K3]

# For superadditivity of R_4, we need SUBadditivity of -R_4:
# -R_4(p+q) <= -R_4(p) + -R_4(q)
# i.e., for each term Ti, we need the gap Ti(p+q) - Ti(p) - Ti(q) <= 0.
# But PROVER-11 showed term-by-term doesn't work because T3 isn't subadditive.

# NEW IDEA: Regroup the terms differently!
# Let's combine T1+T3 and see if THAT combination is subadditive.

print("\n" + "="*70)
print("STEP 2: Alternative decomposition of -R_4")
print("="*70)

K2, K3, K4 = symbols('K2 K3 K4')
D1 = 4*K2**2 - K4
D2 = 4*K2**3 + K2*K4 - 2*K3**2

# Let me try a different partial fraction. Instead of decomposing in K4,
# let me try to express -R_4 as:
# f(K2, K4)/D1 + g(K2, K3)/D2 + h(K2, K3)
# where f, g, h are chosen to make each piece separately subadditive.

# Original: -R_4 = (16K2^3*K3^2 + 4K2^2*K4^2 - 20K2*K3^2*K4 + 8K3^4 + K4^3) / (24*D1*D2)
neg_R4_num = 16*K2**3*K3**2 + 4*K2**2*K4**2 - 20*K2*K3**2*K4 + 8*K3**4 + K4**3

# Try to write neg_R4_num = A(K2,K3,K4)*D2 + B(K2,K3,K4)*D1
# This is polynomial division in the ideal (D1, D2).

# Actually, let me try a cleaner decomposition.
# Note that D1 = 4K2^2 - K4 and D2 = 4K2^3 + K2*K4 - 2*K3^2
# So K4 = 4K2^2 - D1 and K3^2 = (4K2^3 + K2*K4 - D2)/2 = (4K2^3 + K2*(4K2^2-D1) - D2)/2
#   = (8K2^3 - K2*D1 - D2)/2

# Substitute K4 = 4K2^2 - D1 and K3^2 = (8K2^3 - K2*D1 - D2)/2 into the numerator
d1, d2 = symbols('d1 d2', positive=True)
K4_sub = 4*K2**2 - d1
K3sq_sub = (8*K2**3 - K2*d1 - d2) / 2

neg_R4_num_sub = neg_R4_num.subs(K4, K4_sub).subs(K3**2, K3sq_sub)
neg_R4_num_sub = expand(neg_R4_num_sub)
print(f"neg_R4_num in (d1, d2) coords: {neg_R4_num_sub}")
neg_R4_num_sub_collected = collect(expand(neg_R4_num_sub), [d1, d2])
print(f"Collected: {neg_R4_num_sub_collected}")

# The denominator is 24*d1*d2, so:
# -R_4 = neg_R4_num_sub / (24*d1*d2)
# Let me expand and group by powers of d1 and d2
P_sub = Poly(expand(neg_R4_num_sub), d1, d2)
print(f"\nAs polynomial in (d1, d2):")
for monom, coeff in P_sub.as_dict().items():
    i, j = monom
    print(f"  d1^{i} * d2^{j} : {factor(coeff)}")

print("\n" + "="*70)
print("STEP 3: Try homogeneous rescaling")
print("="*70)

# R_4 is degree-2 homogeneous: R_4(λ^2 K2, λ^3 K3, λ^4 K4) = λ^2 R_4(K2,K3,K4)
# (Since cumulants scale as K_n ~ λ^n under λ-rescaling)
# Actually check: numerator has degree 3+2+3+2+4 mix...
# Let me check homogeneity by assigning weights: w(K2)=2, w(K3)=3, w(K4)=4
# num = -16K2^3*K3^2 - 4K2^2*K4^2 + 20K2*K3^2*K4 - 8K3^4 - K4^3
# weights: 6+6=12, 4+8=12, 2+6+4=12, 12, 12 ✓ all weight 12
# den = 24*(4K2^2-K4)*(4K2^3+K2*K4-2K3^2)
# D1 weight: 4 (each term), D2 weight: 6 (each term)
# den weight: 4+6=10
# So R_4 has weight 12-10 = 2. ✓

# Use this: set K2 = 1 (rescale), define x = K3, y = K4.
# Then R_4(1, x, y) = f(x,y), and by homogeneity:
# R_4(K2, K3, K4) = K2 * f(K3/K2^{3/2}, K4/K2^2)

# Wait, that's exactly what PROVER-11 tried (joint concavity of f) and it FAILED.
# But let me think differently about this...

# Instead of concavity of f, maybe we can use a DIFFERENT parameterization
# that makes superadditivity transparent.

print("Exploring the scaled variables approach...")

# Let α = K3^2/(K2^3), β = K4/K2^2
# Domain: β < 4, α < (4+β)/2 (from D2 > 0: 2K3^2 < 4K2^3 + K2*K4 ⟹ 2α < 4+β)
# Then: R_4(K2, K3, K4) = K2 * r(α, β)
# where r(α,β) = (-16α - 4β^2 + 20αβ - 8α^2 - β^3) / (24*(4-β)*(4+β-2α))

# Wait, let me be careful. K3^2 = α*K2^3, K4 = β*K2^2.
# num = K2^6 * (-16α - 4β^2 + 20αβ - 8α^2 - β^3)  ... let me check
# -16*K2^3*(K3^2) = -16*K2^3*(α*K2^3) = -16α*K2^6  ✓
# -4*K2^2*K4^2 = -4*K2^2*(β*K2^2)^2 = -4β^2*K2^6   ✓
# 20*K2*K3^2*K4 = 20*K2*(α*K2^3)*(β*K2^2) = 20αβ*K2^6  ✓
# -8*K3^4 = -8*(α*K2^3)^2 = -8α^2*K2^6  ✓
# -K4^3 = -(β*K2^2)^3 = -β^3*K2^6  ✓
# den = 24*(4K2^2-K4)*(4K2^3+K2*K4-2K3^2)
#     = 24*(K2^2*(4-β))*(K2^3*(4+β-2α))
#     = 24*K2^5*(4-β)*(4+β-2α)
# So R_4 = K2^6/(24*K2^5) * (-16α - 4β^2 + 20αβ - 8α^2 - β^3)/((4-β)*(4+β-2α))
# = K2/24 * (-16α - 4β^2 + 20αβ - 8α^2 - β^3)/((4-β)*(4+β-2α))

# So R_4(K2,K3,K4) = K2 * r(α,β)
# where r(α,β) = (-16α - 4β^2 + 20αβ - 8α^2 - β^3)/(24*(4-β)*(4+β-2α))
# and α = K3^2/K2^3, β = K4/K2^2.

# For superadditivity: R_4(p+q) >= R_4(p) + R_4(q)
# (s+t)*r(α_r, β_r) >= s*r(α_p, β_p) + t*r(α_q, β_q)
# where α_r = (a+b)^2/(s+t)^3, β_r = (u+v)/(s+t)^2.

# This is a PERSPECTIVE inequality! If r were concave (it's not in general),
# this would follow. But maybe we can decompose r into pieces where the
# perspective argument works for each piece...

print("\nLet me try a different decomposition using substitution variables.")
print("Setting u_p = 4*s^2 - k4_p (>0) and w_p = 4*s^3 + s*k4_p - 2*a^2 (>0)")
print("i.e., D1_p = u_p, D2_p = w_p")

# Actually, let me try Approach B more carefully.

print("\n" + "="*70)
print("STEP 4: Approach B - Decompose into k3=0 part + correction")
print("="*70)

# From the partial fractions:
# -R_4 = T1 + T2 + T3 + T4
# At k3=0:
# T1|_{k3=0} = 4K2^3/(6*(4K2^2-K4)) = 2K2^3/(3*(4K2^2-K4))
#            = 2K2^3/(3*D1)
# T2|_{k3=0} = 0 (has K3^2 factor)
# T3|_{k3=0} = -K4/(24*K2)
# T4|_{k3=0} = -K2/6
# So -R_4|_{k3=0} = 2K2^3/(3*D1) - K4/(24*K2) - K2/6
# Let me verify: -R_4|_{k3=0} = K4^2/(24*K2*D1)
R4_k30 = -K4**2/(24*K2*(4*K2**2-K4))
neg_R4_k30 = K4**2/(24*K2*D1)
check = cancel(neg_R4_k30 - (S(2)*K2**3/(3*D1) - K4/(24*K2) - K2/6))
print(f"Check partial fractions at k3=0: {check}")
# Hmm, that should be 0. Let me compute.
pf_k30 = S(2)*K2**3/(3*D1) - K4/(24*K2) - K2/6
diff_k30 = cancel(pf_k30 - neg_R4_k30)
print(f"  pf(k3=0) - neg_R4(k3=0) = {diff_k30}")

# So the k3-dependent correction is:
# Δ(-R_4) = -R_4 - (-R_4|_{k3=0})
# = [T1 + T2 + T3 + T4] - [T1|_{k3=0} + 0 + T3 + T4]
# = (T1 - T1|_{k3=0}) + T2
# = [(4K2^3-K3^2)/(6*D1) - 4K2^3/(6*D1_{k3=0})]  Wait, D1 doesn't depend on K3!
# So T1 - T1|_{k3=0} = (4K2^3-K3^2-4K2^3)/(6*D1) = -K3^2/(6*D1)
# And T2 stays.

# So: correction = -K3^2/(6*D1) + K3^2*(4K2^3-K3^2)/(6*K2^2*D2)
# = K3^2/(6) * [-1/D1 + (4K2^3-K3^2)/(K2^2*D2)]
# = K3^2/(6) * [-K2^2*D2 + D1*(4K2^3-K3^2)] / (K2^2*D1*D2)

corr_num = expand(-K2**2*D2 + D1*(4*K2**3-K3**2))
print(f"\nCorrection numerator (inside K3^2/(6*K2^2*D1*D2)):")
print(f"  {corr_num}")
print(f"  factored: {factor(corr_num)}")

# So:
# -R_4 = K4^2/(24*K2*D1) + K3^2 * corr_num / (6*K2^2*D1*D2)
# The first part is the k3=0 piece (PROVED subadditive).
# The second part is the k3 correction.

# Let me verify this identity.
neg_R4_full = (16*K2**3*K3**2 + 4*K2**2*K4**2 - 20*K2*K3**2*K4 + 8*K3**4 + K4**3) / (24*D1*D2)
correction = K3**2 * corr_num / (6*K2**2*D1*D2)
verify = cancel(neg_R4_full - neg_R4_k30 - correction)
print(f"\nVerify: -R_4 = (k3=0 part) + correction? Diff = {verify}")

print("\n" + "="*70)
print("STEP 5: Analyze the correction term")
print("="*70)
print(f"corr_num = {expand(corr_num)}")
# Substitute D1, D2 back
corr_num_expanded = expand(corr_num.subs([(D1, 4*K2**2-K4), (D2, 4*K2**3+K2*K4-2*K3**2)]))
print(f"corr_num (expanded) = {corr_num_expanded}")
print(f"corr_num (factored) = {factor(corr_num_expanded)}")
