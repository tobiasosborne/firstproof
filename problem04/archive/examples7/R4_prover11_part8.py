"""
PROVER-11 Part 8: Partial fraction approach for the full case.

The partial fraction decomposition of -R4 in K4 is:
  -R4 = (4K2^3 - K3^2)/(6*(4K2^2 - K4))
       + K3^2*(4K2^3 - K3^2)/(6*K2^2*(4K2^3 + K2*K4 - 2K3^2))
       - K4/(24*K2)
       - (2K2^3 + K3^2)/(12*K2^2)

Let me verify this and analyze each term.

Note: 4K2^3 - K3^2 = K2*(D2 - K2*K4)/K2 + stuff... let me just check.
On the domain, D2 = 4K2^3 + K2*K4 - 2K3^2 > 0.
When K4 = 0: D2 = 4K2^3 - 2K3^2 > 0, so K3^2 < 2K2^3.
Then 4K2^3 - K3^2 > 4K2^3 - 2K2^3 = 2K2^3 > 0.
More generally: 4K2^3 - K3^2 > D2/2 + K2*K4/2... depends on sign of K4.

Actually from D2 > 0: 2K3^2 < 4K2^3 + K2*K4.
So K3^2 < 2K2^3 + K2*K4/2.
And 4K2^3 - K3^2 > 4K2^3 - 2K2^3 - K2*K4/2 = 2K2^3 - K2*K4/2 = K2*(2K2^2 - K4/2).
Since K4 < 4K2^2: K4/2 < 2K2^2, so 2K2^2 - K4/2 > 0. ✓
So 4K2^3 - K3^2 > 0 on the domain.
"""
from sympy import (symbols, expand, factor, collect, cancel, apart,
                   Rational, simplify, S, numer, denom)
import numpy as np

K2, K3, K4 = symbols('K2 K3 K4')

# Verify the partial fraction decomposition
D1 = 4*K2**2 - K4
D2 = 4*K2**3 + K2*K4 - 2*K3**2

# Terms:
T1 = (4*K2**3 - K3**2) / (6*D1)
T2 = K3**2*(4*K2**3 - K3**2) / (6*K2**2*D2)
T3 = -K4/(24*K2)
T4 = -(2*K2**3 + K3**2)/(12*K2**2)

neg_R4_pf = T1 + T2 + T3 + T4

# Original -R4
P = -16*K2**3*K3**2 - 4*K2**2*K4**2 + 20*K2*K3**2*K4 - 8*K3**4 - K4**3
neg_R4_orig = -P / (24*D1*D2)

diff = cancel(neg_R4_pf - neg_R4_orig)
print(f"Verification: partial fractions - original = {diff}")
assert diff == 0, "Partial fraction mismatch!"
print("VERIFIED: partial fraction decomposition is correct.\n")

# Now let's rewrite:
# T3 + T4 = -K4/(24*K2) - (2*K2^3 + K3^2)/(12*K2^2)
# = -K4/(24*K2) - K2/6 - K3^2/(12*K2^2)
# These are:
# -K2/6: linear in K2, hence additive (gap = 0)
# -K4/(24*K2): we need to check subadditivity
# -K3^2/(12*K2^2): this is exactly like the n=3 case!

# So: -R4 = T1 + T2 + [-K2/6] + [-K4/(24*K2)] + [-K3^2/(12*K2^2)]
# Wait, let me regroup:
# -R4 = (4K2^3-K3^2)/(6*D1) + K3^2*(4K2^3-K3^2)/(6*K2^2*D2) - K4/(24*K2) - K2/6 - K3^2/(12*K2^2)

# Linear part: -K2/6 (additive, gap = 0)
# Cauchy-Schwarz-amenable: -K3^2/(12*K2^2) (like n=3 case, known subadditive)
# Term -K4/(24*K2): is x/y subadditive? K4/(K2) is linear/linear.
#   K4_r/(K2_r) = (K4_p+K4_q)/(K2_p+K2_q) -- this is a weighted average of K4_p/K2_p and K4_q/K2_q
#   So -K4/(24*K2) evaluated at r is between -K4_p/(24*K2_p) and -K4_q/(24*K2_q).
#   NOT subadditive in general!

# Let me reconsider. Maybe combine T3 with T1 or T4.
# Actually let me try: combine T1 and T3.
# T1 + T3 = (4K2^3-K3^2)/(6*D1) - K4/(24*K2)
# = (4K2^3-K3^2)/(6*(4K2^2-K4)) - K4/(24*K2)
# = [4*K2*(4K2^3-K3^2) - K4*(4K2^2-K4)] / (24*K2*(4K2^2-K4))
comb_13_num = expand(4*K2*(4*K2**3-K3**2) - K4*(4*K2**2-K4))
print(f"T1+T3 numerator: {comb_13_num}")
print(f"  factored: {factor(comb_13_num)}")
print(f"  = {comb_13_num}")

# = 16K2^4 - 4K2*K3^2 - 4K2^2*K4 + K4^2
# Denominator: 24*K2*(4K2^2-K4) = 24*K2*D1

# So T1+T3 = (16K2^4 - 4K2*K3^2 - 4K2^2*K4 + K4^2) / (24*K2*D1)
# = (K4-4K2^2+...)... let me factor numerator:
# 16K2^4 - 4K2^2*K4 + K4^2 = (4K2^2)^2 - 2*(4K2^2)*(K4/2) + (K4/2... no
# Actually: 16K2^4 - 4K2^2*K4 + K4^2 - 4K2*K3^2
# = (4K2^2-K4)^2/... hmm
# (4K2^2-K4)^2 = 16K2^4 - 8K2^2*K4 + K4^2
# So 16K2^4 - 4K2^2*K4 + K4^2 = (4K2^2-K4)^2 + 4K2^2*K4
# And numerator = (4K2^2-K4)^2 + 4K2^2*K4 - 4K2*K3^2
# = D1^2 + 4*K2*(K2*K4 - K3^2)

# From D2 > 0: K2*K4 > 2K3^2 - 4K2^3. But K2*K4 - K3^2 could be negative.

# Let me try yet another regrouping.
# Write -R4 = A/(D1) + B/(D2) + linear_terms
# where A, B are to be determined to make each piece subadditive.

# Actually, let me try the SIMPLEST approach: direct numerical test +
# try Schur-convexity on the full gap polynomial.

print("\n" + "="*70)
print("FULL CASE: Numerical verification")
print("="*70)

np.random.seed(123)
violations = 0
n_trials = 500000
valid = 0

for _ in range(n_trials):
    sp = np.random.exponential(1) + 0.01
    tp = np.random.exponential(1) + 0.01

    # k3 values - arbitrary real
    k3p_val = np.random.normal(0, sp**1.5)
    k3q_val = np.random.normal(0, tp**1.5)

    # k4 constraints:
    # k4 < 4*k2^2 (D1 > 0)
    # k4 > (2*k3^2 - 4*k2^3)/k2 = 2*k3^2/k2 - 4*k2^2 (D2 > 0)
    k4p_lo = 2*k3p_val**2/sp - 4*sp**2 + 0.001
    k4p_hi = 4*sp**2 - 0.001
    if k4p_lo >= k4p_hi:
        continue
    k4p_val = np.random.uniform(k4p_lo, k4p_hi)

    k4q_lo = 2*k3q_val**2/tp - 4*tp**2 + 0.001
    k4q_hi = 4*tp**2 - 0.001
    if k4q_lo >= k4q_hi:
        continue
    k4q_val = np.random.uniform(k4q_lo, k4q_hi)

    # Check sum domain
    sr = sp + tp
    k3r = k3p_val + k3q_val
    k4r = k4p_val + k4q_val
    D1r = 4*sr**2 - k4r
    D2r = 4*sr**3 + sr*k4r - 2*k3r**2
    if D1r <= 0 or D2r <= 0:
        continue

    valid += 1

    def R4_eval(k2, k3, k4):
        n = -16*k2**3*k3**2 - 4*k2**2*k4**2 + 20*k2*k3**2*k4 - 8*k3**4 - k4**3
        d = 24*(4*k2**2-k4)*(4*k2**3+k2*k4-2*k3**2)
        return n/d

    rp = R4_eval(sp, k3p_val, k4p_val)
    rq = R4_eval(tp, k3q_val, k4q_val)
    rr = R4_eval(sr, k3r, k4r)
    gap = rr - rp - rq

    if gap < -1e-9:
        violations += 1
        if violations <= 3:
            print(f"VIOLATION at trial {_}: gap={gap:.6e}")
            print(f"  sp={sp:.4f}, tp={tp:.4f}")

print(f"\n{violations} violations in {valid} valid trials")
if violations == 0:
    print("FULL CASE: NUMERICALLY VERIFIED (0 violations)")
