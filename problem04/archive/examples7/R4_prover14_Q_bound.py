"""
PROVER-14 Part 6: Proving Q_r <= max(Q_p, Q_q) and consequences

Q = Sm2 * S2 = [sum_{i<j} 1/(lambda_i - lambda_j)^2] * [sum_{i<j} (lambda_i - lambda_j)^2]

This is scale-invariant (dimensionless). By Cauchy-Schwarz, Q >= [n(n-1)/2]^2.
Equality iff all pairwise distances are equal (impossible for n >= 3).

FACT: For n=2, Q = 1 always (Q_p = Q_q = Q_r = 1). The conjecture gives equality.

For n >= 3, Q depends on the "shape" of the root configuration.
More precisely, Q depends on the ratios d_i/d_j of consecutive gaps.

KEY CLAIM (numerically verified): Q_r <= max(Q_p, Q_q).
If proved, this combined with S2 additivity gives the main conjecture.

PROOF IDEA: Show that MSS convolution "regularizes" the gap structure,
making pairwise distances more uniform. Since Q is minimized at uniform
gaps, regularization decreases Q.

ALTERNATIVE: The main conjecture Q_r <= Q_p*Q_q*(alpha+beta)/(alpha*Q_q+beta*Q_p)
is equivalent to the original Fisher superadditivity.
"""
import numpy as np
from itertools import combinations
from math import factorial
import sympy as sp

print("="*70)
print("PROVER-14 Part 6: Q-BOUND AND PROOF")
print("="*70)

def H_values(roots):
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H

def Phi_n(roots):
    H = H_values(roots)
    return np.sum(H**2)

def Sm2(roots):
    n = len(roots)
    return sum(1/(roots[i]-roots[j])**2 for i in range(n) for j in range(i+1, n))

def S2(roots):
    n = len(roots)
    return sum((roots[i]-roots[j])**2 for i in range(n) for j in range(i+1, n))

def Q(roots):
    return Sm2(roots) * S2(roots)

def roots_to_monic_coeffs(roots):
    n = len(roots)
    coeffs = [1.0]
    for k in range(1, n+1):
        ek = sum(np.prod(list(combo)) for combo in combinations(roots, k))
        coeffs.append((-1)**k * ek)
    return np.array(coeffs)

def mss_convolve_n(p_coeffs, q_coeffs, n):
    r_coeffs = np.zeros(n+1)
    for k in range(n+1):
        ck = 0
        for i in range(k+1):
            j = k - i
            if i <= n and j <= n:
                coeff = factorial(n-i) * factorial(n-j) / (factorial(n) * factorial(n-k))
                ck += coeff * p_coeffs[i] * q_coeffs[j]
        r_coeffs[k] = ck
    return r_coeffs

def mss_convolve_roots(roots_p, roots_q):
    n = len(roots_p)
    p_coeffs = roots_to_monic_coeffs(roots_p)
    q_coeffs = roots_to_monic_coeffs(roots_q)
    r_coeffs = mss_convolve_n(p_coeffs, q_coeffs, n)
    r_roots = np.roots(r_coeffs)
    r_roots = np.sort(np.real(r_roots))
    return r_roots

# =============================================================
# Part 1: Q for n=2
# =============================================================
print("\n--- Part 1: Q for n=2 ---\n")

for d in [1.0, 2.0, 0.5, 3.14]:
    roots = np.array([0, d])
    print(f"  d={d}: Q = {Q(roots):.10f}")

print("\n  For n=2: Q = (1/d^2) * d^2 = 1 always. CONFIRMED.")

# =============================================================
# Part 2: Q for n=3
# =============================================================
print("\n--- Part 2: Q for n=3 ---\n")

# For n=3 with gaps s, t:
# Sm2 = 1/s^2 + 1/t^2 + 1/(s+t)^2
# S2 = s^2 + t^2 + (s+t)^2 = 2(s^2 + st + t^2)
# Q = [1/s^2 + 1/t^2 + 1/(s+t)^2] * [2(s^2 + st + t^2)]

s_sym, t_sym = sp.symbols('s t', positive=True)
Sm2_3 = 1/s_sym**2 + 1/t_sym**2 + 1/(s_sym+t_sym)**2
S2_3 = s_sym**2 + t_sym**2 + (s_sym+t_sym)**2
Q_3 = sp.simplify(sp.expand(Sm2_3 * S2_3))
print(f"  Q_3(s,t) = {Q_3}")

# Factor
Q_3_factored = sp.factor(Q_3)
print(f"  Q_3(s,t) factored = {Q_3_factored}")

# Substitute t = r*s to get Q as function of ratio only
r_sym = sp.Symbol('r', positive=True)
Q_3_ratio = Q_3.subs(t_sym, r_sym*s_sym)
Q_3_ratio_simplified = sp.simplify(Q_3_ratio)
print(f"  Q_3(s, r*s) = {Q_3_ratio_simplified}")

# This should be independent of s (Q is scale-invariant)
# Check: does s cancel?
Q_3_ratio_check = sp.simplify(Q_3_ratio_simplified * s_sym**0)  # Hmm, just check numerically

for s_val in [1.0, 2.0, 0.5]:
    for r_val in [0.5, 1.0, 2.0]:
        t_val = r_val * s_val
        roots = np.array([-(2*s_val+t_val)/3, (s_val-t_val)/3, (s_val+2*t_val)/3])
        print(f"    s={s_val}, r=t/s={r_val}: Q = {Q(roots):.6f}")

# So Q depends only on the gap ratio r = t/s.
# Let f(r) = Q(s, r*s) (should be independent of s)
print(f"\n  Q_3 as function of r=t/s:")
for r_val in np.linspace(0.1, 10, 20):
    roots = np.array([-(2+r_val)/3, (1-r_val)/3, (1+2*r_val)/3])
    print(f"    r={r_val:.2f}: Q = {Q(roots):.6f}")

# Minimum of Q_3
# dQ/dr = 0 at r = 1 (by symmetry Q(r) = Q(1/r))
print(f"\n  Q_3 at r=1 (equal gaps): {Q(np.array([-1, 0, 1])):.6f}")
# Q_3 at r=1: Sm2 = 1+1+1/4 = 9/4, S2 = 1+1+4 = 6, Q = 54/4 = 13.5
# Actually for centered roots (-1, 0, 1): gaps = 1, 1
# Sm2 = 1/1 + 1/1 + 1/4 = 2.25, S2 = 1 + 1 + 4 = 6, Q = 13.5

# What's the minimum of Q_3?
print(f"  Q_3 min = 13.5 at r=1 (equal gaps)")
# By CS: Q >= 9 (since n(n-1)/2 = 3, Q >= 9). But actual min is 13.5.

# =============================================================
# Part 3: What does Q_r <= max(Q_p, Q_q) imply?
# =============================================================
print("\n\n--- Part 3: Q_r <= max(Q_p, Q_q) implications ---\n")

print("""
CLAIM: Q_r <= max(Q_p, Q_q) where Q = Sm2 * S2.

This means: MSS convolution does NOT increase the shape factor Q.

CLAIM implies the MAIN CONJECTURE:
  We need 1/Sm2(r) >= 1/Sm2(p) + 1/Sm2(q).
  Using S2 = alpha + beta (additive):
    1/Sm2(r) = S2(r)/Q_r = (alpha+beta)/Q_r
    1/Sm2(p) = alpha/Q_p
    1/Sm2(q) = beta/Q_q

  Need: (alpha+beta)/Q_r >= alpha/Q_p + beta/Q_q

  If Q_r <= max(Q_p, Q_q), WLOG Q_p >= Q_q, so Q_r <= Q_p.
  Then: (alpha+beta)/Q_r >= (alpha+beta)/Q_p = alpha/Q_p + beta/Q_p
  And: alpha/Q_p + beta/Q_p >= alpha/Q_p + beta/Q_q  iff Q_p >= Q_q. YES!

  Wait: beta/Q_p >= beta/Q_q iff Q_q >= Q_p, which contradicts our assumption.

  So this direction doesn't work directly. Let me reconsider.
""")

print("""
CORRECTION: Q_r <= max(Q_p, Q_q) does NOT directly imply the conjecture.

Let Q_p >= Q_q (WLOG). Then Q_r <= Q_p.
  LHS = (alpha+beta)/Q_r >= (alpha+beta)/Q_p
  RHS = alpha/Q_p + beta/Q_q >= alpha/Q_p + beta/Q_p = (alpha+beta)/Q_p

So LHS >= (alpha+beta)/Q_p <= RHS. This doesn't help.

The issue is that Q_r <= Q_p makes 1/Q_r >= 1/Q_p, which helps the LHS,
but the RHS has beta/Q_q and if Q_q < Q_p, then beta/Q_q > beta/Q_p.

We need a DIFFERENT structural property of MSS convolution.
""")

# =============================================================
# Part 4: Try a weighted version
# =============================================================
print("\n--- Part 4: Weighted Q bound ---\n")

print("""
The exact condition is:
  Q_r <= Q_p * Q_q * (alpha+beta) / (alpha * Q_q + beta * Q_p)

where alpha = S2(p), beta = S2(q).

The RHS is the WEIGHTED HARMONIC MEAN of Q_p and Q_q with weights alpha, beta.

Define: W(Q_p, Q_q; alpha, beta) = Q_p * Q_q * (alpha+beta) / (alpha*Q_q + beta*Q_p)

Properties:
- W(Q, Q; alpha, beta) = Q (if Q_p = Q_q = Q, then W = Q)
- W >= min(Q_p, Q_q) (weighted harmonic mean >= min)
- W is increasing in both Q_p and Q_q

So the condition Q_r <= W is WEAKER than Q_r <= min(Q_p, Q_q).
In fact: W = Q_p*Q_q / weighted_avg >= min(Q_p, Q_q).
Actually W = (alpha+beta) / (alpha/Q_p + beta/Q_q), which is >= min(Q_p, Q_q).

The exact condition is: Q_r <= weighted_harmonic_mean(Q_p, Q_q; alpha, beta).
""")

# Verify: Q_r <= W always
np.random.seed(42)
violations_W = 0
min_ratio_W = float('inf')
for trial in range(20000):
    n = 4
    p_roots = np.sort(np.random.randn(n) * 2)
    q_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(p_roots)) < 0.1:
        p_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(q_roots)) < 0.1:
        q_roots = np.sort(np.random.randn(n) * 2)
    try:
        r_roots = mss_convolve_roots(p_roots, q_roots)
        if np.min(np.diff(r_roots)) < 1e-6:
            continue
        Qp = Q(p_roots)
        Qq = Q(q_roots)
        Qr = Q(r_roots)
        alpha = S2(p_roots)
        beta = S2(q_roots)
        W = Qp * Qq * (alpha+beta) / (alpha*Qq + beta*Qp)
        ratio = Qr / W
        min_ratio_W = min(min_ratio_W, ratio)
        if Qr > W + 1e-6:
            violations_W += 1
    except:
        continue
print(f"Violations of Q_r <= W: {violations_W}, min(Q_r/W) = {min_ratio_W:.6f}")
print("(W = weighted harmonic mean of Q_p, Q_q with weights S2_p, S2_q)")

# =============================================================
# Part 5: What determines Q under MSS?
# =============================================================
print("\n\n--- Part 5: Q under MSS convolution ---\n")

# For n=3, Q depends only on the gap ratio r = t/s.
# Under MSS convolution, r_p, r_q -> r_r.
# Study this map.

def gap_ratio_n3(roots):
    """Gap ratio t/s for centered n=3 roots."""
    gaps = np.diff(roots)
    return gaps[1] / gaps[0]

np.random.seed(42)
print("Gap ratio transformation under MSS for n=3:")
print("(r_p, r_q) -> r_r\n")

data = []
for trial in range(20):
    n = 3
    p_roots = np.sort(np.random.randn(n) * 2)
    q_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(p_roots)) < 0.3:
        p_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(q_roots)) < 0.3:
        q_roots = np.sort(np.random.randn(n) * 2)

    # Center both
    p_roots -= np.mean(p_roots)
    q_roots -= np.mean(q_roots)

    try:
        r_roots = mss_convolve_roots(p_roots, q_roots)
        r_roots -= np.mean(r_roots)

        r_p = gap_ratio_n3(p_roots)
        r_q = gap_ratio_n3(q_roots)
        r_r = gap_ratio_n3(r_roots)

        data.append((r_p, r_q, r_r, Q(p_roots), Q(q_roots), Q(r_roots)))
        print(f"  r_p={r_p:.4f}, r_q={r_q:.4f} -> r_r={r_r:.4f} (Q: {Q(p_roots):.2f},{Q(q_roots):.2f}->{Q(r_roots):.2f})")
    except:
        continue

# =============================================================
# Part 6: n=3 symbolic MSS
# =============================================================
print("\n\n--- Part 6: n=3 symbolic MSS ---\n")

# Centered n=3 polynomial: p(x) = x^3 + px + q
# where p = e2 (negative of sum of products of pairs)
# and q = -e3 (product of roots)
# With roots a, b, c = -(a+b) (centered):
# e2 = ab + ac + bc = ab + (a+b)*c ... with c = -(a+b):
# e2 = ab - a(a+b) - b(a+b) = ab - a^2 - ab - ab - b^2 = -a^2 - ab - b^2
# e3 = abc = ab(-(a+b)) = -ab(a+b)
#
# With gap parametrization: a = -(2s+t)/3, b = (s-t)/3:
# e2 = -(a^2 + ab + b^2) = -((2s+t)^2/9 + (2s+t)(s-t)/(-9) + (s-t)^2/9)
# Hmm, let me just use sympy.

s, t, u, v = sp.symbols('s t u v', positive=True)

# p roots: -(2s+t)/3, (s-t)/3, (s+2t)/3
a_p = -(2*s+t)/3
b_p = (s-t)/3
c_p = (s+2*t)/3

# q roots: -(2u+v)/3, (u-v)/3, (u+2v)/3
a_q = -(2*u+v)/3
b_q = (u-v)/3
c_q = (u+2*v)/3

# Coefficients of p: x^3 + e2_p x + e3_p
e1_p = a_p + b_p + c_p  # should be 0
e2_p = sp.expand(a_p*b_p + a_p*c_p + b_p*c_p)
e3_p = sp.expand(-a_p*b_p*c_p)
print(f"e1_p = {sp.simplify(e1_p)}")
print(f"e2_p = {e2_p}")
print(f"e3_p = {e3_p}")

e1_q = a_q + b_q + c_q
e2_q = sp.expand(a_q*b_q + a_q*c_q + b_q*c_q)
e3_q = sp.expand(-a_q*b_q*c_q)
print(f"\ne2_q = {e2_q}")
print(f"e3_q = {e3_q}")

# MSS convolution for n=3, centered:
# p(x) = x^3 + 0*x^2 + e2_p*x + e3_p -> coeffs [1, 0, e2_p, e3_p]
# q(x) = x^3 + 0*x^2 + e2_q*x + e3_q -> coeffs [1, 0, e2_q, e3_q]

# c_0 = 1, c_1 = 0 (both centered)
# c_2 = (2!*1!)/(3!*1!) * 0*0 + (1!*2!)/(3!*1!) * 0*0 + e2_p + e2_q  ... let me be careful

# c_k = sum_{i+j=k} [(n-i)!(n-j)!/(n!(n-k)!)] * a_i * b_j where n=3
# c_2: i=0,j=2: 3!*1!/(3!*1!) * 1 * e2_q = e2_q
#       i=1,j=1: 2!*2!/(3!*1!) * 0 * 0 = 0
#       i=2,j=0: 1!*3!/(3!*1!) * e2_p * 1 = e2_p
# c_2 = e2_p + e2_q

# c_3: i=0,j=3: 3!*0!/(3!*0!) * 1 * e3_q = e3_q
#       i=1,j=2: 2!*1!/(3!*0!) * 0 * e2_q = 0
#       i=2,j=1: 1!*2!/(3!*0!) * e2_p * 0 = 0
#       i=3,j=0: 0!*3!/(3!*0!) * e3_p * 1 = e3_p
# c_3 = e3_p + e3_q

print("\n\nMSS for centered n=3:")
print("  r(x) = x^3 + (e2_p + e2_q)*x + (e3_p + e3_q)")
print("  i.e., e2 and e3 are ADDITIVE!")

c2_r = e2_p + e2_q
c3_r = e3_p + e3_q

print(f"\n  e2_p = {e2_p} = -(s^2+st+t^2)/3")
print(f"  e3_p = {e3_p}")
# Simplify e3
e3_p_simplified = sp.factor(e3_p)
print(f"  e3_p factored = {e3_p_simplified}")

# e2_r = e2_p + e2_q
# Note: e2 = -(s^2+st+t^2)/3 is related to S2:
# S2 = 2(s^2+st+t^2), so e2 = -S2/6
# e2_r = e2_p + e2_q = -(S2_p + S2_q)/6 = -S2_r/6
# This confirms S2 additivity!

print("\n  Since e2 = -(s^2+st+t^2)/3 = -S2/6,")
print("  e2_r = e2_p + e2_q confirms S2 additivity.\n")

# Now: r(x) = x^3 + e2_r*x + e3_r
# The roots of r satisfy: lambda_1 + lambda_2 + lambda_3 = 0 (centered)
# lambda_1*lambda_2 + lambda_1*lambda_3 + lambda_2*lambda_3 = e2_r
# lambda_1*lambda_2*lambda_3 = -e3_r

# The gaps of r: s_r, t_r satisfy:
# S2_r = 2(s_r^2 + s_r*t_r + t_r^2) = -6*e2_r = 6*(e2_p + e2_q)... wait.
# e2_r = -(s_r^2+s_r*t_r+t_r^2)/3, so s_r^2+s_r*t_r+t_r^2 = -3*e2_r

# And e3_r is related to the skewness:
# e3_r = (2s_r+t_r)(s_r-t_r)(s_r+2t_r)/27
e3_check = sp.expand((2*s+t)*(s-t)*(s+2*t)/27)
print(f"  e3 as function of gaps = {e3_check}")
print(f"  compare with e3_p = {e3_p_simplified}")
# e3_p = -(2s+t)(s-t)(s+2t)/27 ... let me check signs
# e3 = -lambda_1*lambda_2*lambda_3 = -(-(2s+t)/3)*((s-t)/3)*((s+2t)/3)
# = (2s+t)(s-t)(s+2t)/27
print(f"  e3 = (2s+t)(s-t)(s+2t)/27: {sp.expand(e3_check - e3_p) == 0}")

# =============================================================
# Part 7: n=3 in terms of (sigma, tau) = (e2, e3)
# =============================================================
print("\n\n--- Part 7: n=3 Phi in terms of (sigma, tau) ---\n")

# For centered n=3: sigma = e2, tau = e3
# sigma = -(s^2+st+t^2)/3 < 0
# tau = (2s+t)(s-t)(s+2t)/27

# Sm2 = (s^2+st+t^2)^2 / (s^2*t^2*(s+t)^2)
# S2 = 2*(s^2+st+t^2) = -6*sigma
# Q = Sm2 * S2 = 2*(s^2+st+t^2)^3 / (s^2*t^2*(s+t)^2)

# Can express Q in terms of sigma and tau?
# sigma^3 = -(s^2+st+t^2)^3 / 27
# tau^2 = (2s+t)^2*(s-t)^2*(s+2t)^2 / 729

# The discriminant: Delta = -4*sigma^3 - 27*tau^2
# = 4*(s^2+st+t^2)^3/27 - 27*(2s+t)^2*(s-t)^2*(s+2t)^2/729
# = 4*(s^2+st+t^2)^3/27 - (2s+t)^2*(s-t)^2*(s+2t)^2/27
# = [4*(s^2+st+t^2)^3 - (2s+t)^2*(s-t)^2*(s+2t)^2] / 27

# Also: s^2*t^2*(s+t)^2 can be expressed in terms of sigma and tau?
# We have: s + t = spread_12 (but this is related to the gap structure)
# Actually, (s*t*(s+t))^2 = s^2*t^2*(s+t)^2
# And: s*t = product of consecutive gaps
# s+t = sum of consecutive gaps = total spread

# The key: the DISCRIMINANT Delta = -4*sigma^3 - 27*tau^2 of x^3 + sigma*x + tau = 0
# determines whether roots are real and distinct.

sigma_sym = sp.Symbol('sigma', negative=True)
tau_sym = sp.Symbol('tau', real=True)

# For the cubic x^3 + sigma*x + tau:
# Delta = -4*sigma^3 - 27*tau^2 > 0 for distinct real roots.

# We need Q = Sm2 * S2 in terms of sigma and tau.
# S2 = -6*sigma
# Sm2 = (s^2+st+t^2)^2 / (s^2*t^2*(s+t)^2) = 9*sigma^2 / (s^2*t^2*(s+t)^2)

# Now: s^2*t^2*(s+t)^2 in terms of sigma, tau:
# Let e2 = -(s^2+st+t^2)/3 = sigma, so s^2+st+t^2 = -3*sigma
# Let e3 = tau, so (2s+t)(s-t)(s+2t)/27 = tau

# Also: s*t*(s+t) = ?
# st = (s+t)^2 - (s^2+t^2) ... hmm
# Actually: s+t = D (total spread), st = product.
# (s+t)^2 = s^2 + 2st + t^2, so st = ((s+t)^2 - s^2 - t^2)/2
# s^2+st+t^2 = -3*sigma means s^2+t^2 = -3*sigma - st
# And (s+t)^2 = s^2 + 2st + t^2 = -3*sigma + st
# So D^2 = -3*sigma + st, giving st = D^2 + 3*sigma

# Also: e3 = (2s+t)(s-t)(s+2t)/27 = tau
# (2s+t)(s-t)(s+2t) = 27*tau
# Expand: (2s+t)(s-t) = 2s^2 - 2st + st - t^2 = 2s^2 - st - t^2
# Then *(s+2t) = 2s^3 + 4s^2t - s^2t - 2st^2 - st^2 - 2t^3
# = 2s^3 + 3s^2t - 3st^2 - 2t^3
# = 2(s^3-t^3) + 3st(s-t) = 2(s-t)(s^2+st+t^2) + 3st(s-t)
# = (s-t)[2(s^2+st+t^2) + 3st] = (s-t)(2s^2+5st+2t^2)
# = (s-t)(2s+t)(s+2t)  ... circular

# Let me instead express the discriminant:
# Delta = -4*sigma^3 - 27*tau^2
# = 4*(-3*sigma)^3/27 ... wait
# Delta of x^3 + sigma*x + tau = 0 is -4*sigma^3 - 27*tau^2
# For distinct real roots: Delta > 0, i.e., -4*sigma^3 > 27*tau^2
# Since sigma < 0: -4*sigma^3 = 4*|sigma|^3 > 0. Good.

# Now: s^2*t^2*(s+t)^2 relates to sigma and tau how?
# I can express it using the fact that the roots of the cubic are functions of sigma and tau.

# Actually, for the cubic x^3 + sigma*x + tau = 0 with roots a < b < c (a+b+c=0):
# The roots a, b, c can be written in terms of trigonometric formula:
# a = 2*sqrt(-sigma/3) * cos(theta + 2pi/3)
# b = 2*sqrt(-sigma/3) * cos(theta - 2pi/3)  -- wait, need to check order
# c = 2*sqrt(-sigma/3) * cos(theta)
# where cos(3*theta) = -tau / (2*(-sigma/3)^{3/2}) = 3*sqrt(3)*tau / (2*(-sigma)^{3/2})

# The gaps: s = b - a, t = c - b
# Using trig form... this is getting complex. Let me use a cleaner variable.

# For the centered cubic, introduce the SKEWNESS parameter:
# Let rho = 27*tau^2 / (-4*sigma^3) in [0, 1)  (rho = 0 for symmetric, rho -> 1 for degenerate)
# Then Q depends only on rho!

print("Q as function of skewness rho = 27*tau^2/(-4*sigma^3):")
for rho_val in np.linspace(0, 0.99, 20):
    # For given rho, choose sigma = -3, tau = sqrt(rho * 4 * 27 / 27) = sqrt(4*rho) * 3^{3/2}/3^{3/2}
    # Actually tau^2 = rho * (-4*sigma^3) / 27 = rho * 4 * 27 / 27 = 4*rho (for sigma=-3)
    sigma_val = -3.0
    tau_val = np.sqrt(rho_val * (-4*sigma_val**3) / 27)

    # Find roots of x^3 + sigma*x + tau = 0
    roots_cubic = np.sort(np.real(np.roots([1, 0, sigma_val, tau_val])))
    if np.min(np.diff(roots_cubic)) < 1e-8:
        continue

    q_val = Q(roots_cubic)
    print(f"  rho={rho_val:.4f}: Q = {q_val:.6f}")

# =============================================================
# Part 8: n=3 MSS in (sigma, tau) coordinates
# =============================================================
print("\n\n--- Part 8: n=3 MSS in (sigma, tau) ---\n")

print("""
For n=3 centered polynomials:
  p(x) = x^3 + sigma_p * x + tau_p
  q(x) = x^3 + sigma_q * x + tau_q
  r(x) = x^3 + (sigma_p+sigma_q)*x + (tau_p+tau_q)

So MSS convolution for centered n=3 is just COORDINATE-WISE ADDITION
of (sigma, tau).

This means: the conjecture for n=3 reduces to proving that
  1/Sm2(sigma, tau) is SUPERADDITIVE as a function of (sigma, tau).

Since Sm2 depends on the roots (which depend on sigma, tau),
we need: 1/Sm2(sigma_p+sigma_q, tau_p+tau_q) >= 1/Sm2(sigma_p,tau_p) + 1/Sm2(sigma_q,tau_q)

This is EXACTLY the statement that 1/Sm2 is a SUPERADDITIVE function
of (sigma, tau) on the domain {(sigma,tau): Delta > 0}.
""")

# Let's compute 1/Sm2 as a function of (sigma, tau) and check concavity/superadditivity
print("Checking: is 1/Sm2 CONCAVE in (sigma, tau)?\n")
print("(Concavity implies superadditivity.)\n")

def Sm2_from_sigma_tau(sigma, tau):
    """Sm2 for centered cubic with coefficients sigma, tau."""
    roots = np.sort(np.real(np.roots([1, 0, sigma, tau])))
    return Sm2(roots)

# Check concavity via Hessian
sigma_val = -3.0
tau_val = 0.5
eps = 1e-4

f_00 = 1.0/Sm2_from_sigma_tau(sigma_val, tau_val)

# Partial derivatives
f_p0 = 1.0/Sm2_from_sigma_tau(sigma_val+eps, tau_val)
f_m0 = 1.0/Sm2_from_sigma_tau(sigma_val-eps, tau_val)
f_0p = 1.0/Sm2_from_sigma_tau(sigma_val, tau_val+eps)
f_0m = 1.0/Sm2_from_sigma_tau(sigma_val, tau_val-eps)
f_pp = 1.0/Sm2_from_sigma_tau(sigma_val+eps, tau_val+eps)
f_pm = 1.0/Sm2_from_sigma_tau(sigma_val+eps, tau_val-eps)
f_mp = 1.0/Sm2_from_sigma_tau(sigma_val-eps, tau_val+eps)
f_mm = 1.0/Sm2_from_sigma_tau(sigma_val-eps, tau_val-eps)

H11 = (f_p0 - 2*f_00 + f_m0) / eps**2
H22 = (f_0p - 2*f_00 + f_0m) / eps**2
H12 = (f_pp - f_pm - f_mp + f_mm) / (4*eps**2)

print(f"At (sigma, tau) = ({sigma_val}, {tau_val}):")
print(f"  Hessian of 1/Sm2:")
print(f"    H11 = {H11:.6f}")
print(f"    H22 = {H22:.6f}")
print(f"    H12 = {H12:.6f}")
print(f"    det(H) = {H11*H22 - H12**2:.6f}")
print(f"    1/Sm2 is {'CONCAVE' if H11 < 0 and H11*H22-H12**2 > 0 else 'NOT CONCAVE'} here")

# Check at multiple points
print("\nHessian eigenvalues of 1/Sm2 at various (sigma, tau):")
for sigma_val in [-1.0, -3.0, -10.0]:
    for tau_val in [0.0, 0.5, 1.0]:
        # Check domain
        disc = -4*sigma_val**3 - 27*tau_val**2
        if disc <= 0:
            continue
        try:
            f_00 = 1.0/Sm2_from_sigma_tau(sigma_val, tau_val)
            f_p0 = 1.0/Sm2_from_sigma_tau(sigma_val+eps, tau_val)
            f_m0 = 1.0/Sm2_from_sigma_tau(sigma_val-eps, tau_val)
            f_0p = 1.0/Sm2_from_sigma_tau(sigma_val, tau_val+eps)
            f_0m = 1.0/Sm2_from_sigma_tau(sigma_val, tau_val-eps)
            f_pp = 1.0/Sm2_from_sigma_tau(sigma_val+eps, tau_val+eps)
            f_pm = 1.0/Sm2_from_sigma_tau(sigma_val+eps, tau_val-eps)
            f_mp = 1.0/Sm2_from_sigma_tau(sigma_val-eps, tau_val+eps)
            f_mm = 1.0/Sm2_from_sigma_tau(sigma_val-eps, tau_val-eps)

            H11 = (f_p0 - 2*f_00 + f_m0) / eps**2
            H22 = (f_0p - 2*f_00 + f_0m) / eps**2
            H12 = (f_pp - f_pm - f_mp + f_mm) / (4*eps**2)

            eigvals = np.linalg.eigvalsh([[H11, H12], [H12, H22]])
            concave = all(e <= 1e-6 for e in eigvals)
            print(f"  ({sigma_val:.1f}, {tau_val:.1f}): eigenvalues = [{eigvals[0]:.6f}, {eigvals[1]:.6f}], {'CONCAVE' if concave else 'NOT CONCAVE'}")
        except:
            print(f"  ({sigma_val:.1f}, {tau_val:.1f}): computation failed")

print("\n\nDone.")
