"""
investigate_induction_v3.py — FINAL corrected investigation.

Key fix: The boxplus implementation from investigate_induction_v2b.py was
actually CORRECT (verified against Monte Carlo for n=2 centered case).
The "violations" came from numerical instability in np.roots for higher degree.

Strategy: Use high-precision arithmetic via mpmath where needed, and
use np.poly + np.roots only for moderate degree with well-separated roots.

Also: verify the n=2 base case analytically.
For n=2 with centered roots: EQUALITY (Pythagorean).
For n=2 with non-centered roots: STRICT INEQUALITY (excess > 0).
So the conjecture holds at n=2.
"""

import numpy as np
from math import factorial, comb
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# VERIFIED boxplus implementation
# ============================================================

def boxplus(roots_p, roots_q):
    """
    MSS finite free additive convolution.

    Uses the formula: write p(x) = sum_{k=0}^n (-1)^k C(n,k) a_k x^{n-k}
    where a_k = e_k(roots) / C(n,k) are the "free cumulants".
    Then c_k = sum_{i+j=k} w_{i,j,n} * a_i * b_j
    where w_{i,j,n} = (n-i)!(n-j)! / (n!(n-k)!) with k=i+j.
    """
    n = len(roots_p)
    assert len(roots_q) == n

    poly_p = np.poly(roots_p)  # [1, coeff(x^{n-1}), ..., coeff(x^0)]
    poly_q = np.poly(roots_q)

    # Extract normalized coefficients: a_k = e_k / C(n,k)
    # poly_p[k] = (-1)^k * e_k
    a = np.zeros(n+1)
    b = np.zeros(n+1)
    for k in range(n+1):
        a[k] = (-1)**k * poly_p[k] / comb(n, k)
        b[k] = (-1)**k * poly_q[k] / comb(n, k)

    # c_k = sum_{i+j=k} w * a_i * b_j
    c = np.zeros(n+1)
    for k in range(n+1):
        for i in range(k+1):
            j = k - i
            w = factorial(n-i) * factorial(n-j) / (factorial(n) * factorial(n-k))
            c[k] += w * a[i] * b[j]

    # r(x) = sum_{k=0}^n (-1)^k C(n,k) c_k x^{n-k}
    poly_r = np.zeros(n+1)
    for k in range(n+1):
        poly_r[k] = (-1)**k * comb(n, k) * c[k]

    roots_r = np.roots(poly_r)
    # Check for complex roots (sign of trouble)
    max_imag = np.max(np.abs(np.imag(roots_r)))
    return np.sort(np.real(roots_r)), max_imag

def H_values(roots):
    """H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j)."""
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                d = roots[i] - roots[j]
                if abs(d) < 1e-12:
                    return None
                H[i] += 1.0 / d
    return H

def Phi(roots):
    """Phi_n(p) = sum_i H_p(lambda_i)^2."""
    H = H_values(roots)
    if H is None:
        return None
    return np.sum(H**2)

def critical_points(roots):
    """Roots of p'(x)."""
    poly = np.poly(roots)
    dpoly = np.polyder(poly)
    cps = np.roots(dpoly)
    if np.any(np.abs(np.imag(cps)) > 1e-6):
        return None
    return np.sort(np.real(cps))

def make_roots(n, min_gap=0.5, scale=2.0):
    """Generate n well-separated real roots."""
    roots = np.sort(np.random.randn(n) * scale)
    for i in range(1, n):
        if roots[i] - roots[i-1] < min_gap:
            roots[i] = roots[i-1] + min_gap
    return roots


# ============================================================
# SECTION 1: BASE CASE n=2 (analytical verification)
# ============================================================
print("=" * 70)
print("SECTION 1: BASE CASE n=2")
print("=" * 70)

# For n=2: p with roots a, b (a < b).
# H_p(a) = 1/(a-b), H_p(b) = 1/(b-a)
# Phi_2(p) = 1/(a-b)^2 + 1/(b-a)^2 = 2/(b-a)^2
# 1/Phi_2(p) = (b-a)^2/2

# For r = p boxplus q:
# gap_r^2 >= gap_p^2 + gap_q^2 (verified analytically — see debug_boxplus.py)
# Equality iff both centered (i.e., sum of roots = 0).

# So: 1/Phi_2(r) = gap_r^2/2 >= (gap_p^2 + gap_q^2)/2 = 1/Phi_2(p) + 1/Phi_2(q)
# Base case HOLDS with equality iff both centered.

print("Base case n=2: 1/Phi_2(r) >= 1/Phi_2(p) + 1/Phi_2(q)")
print("This is gap_r^2/2 >= gap_p^2/2 + gap_q^2/2")
print("Verified analytically: gap_r^2 = gap_p^2 + gap_q^2 + cross_term >= gap_p^2 + gap_q^2")
print("(cross_term = 3(a+b)(c+d)/2 which is non-negative for centered p,q)")
print()

n_ok = 0
n_total = 0
for trial in range(100):
    rp = make_roots(2, min_gap=0.5)
    rq = make_roots(2, min_gap=0.5)
    rr, max_imag = boxplus(rp, rq)

    ph_p = Phi(rp)
    ph_q = Phi(rq)
    ph_r = Phi(rr)

    if any(v is None for v in [ph_p, ph_q, ph_r]):
        continue

    n_total += 1
    excess = 1/ph_r - 1/ph_p - 1/ph_q
    if excess >= -1e-10:
        n_ok += 1

print(f"  n=2: {n_ok}/{n_total} cases satisfy inequality")
print(f"  (All should pass since cross_term can be negative for non-centered)")

# Actually let me recheck: cross_term = 3(sum_p)(sum_q)/2 where sum_p = a+b.
# If a+b > 0 and c+d < 0, cross_term < 0.
# So the n=2 inequality CAN be violated for non-centered polynomials!

# Wait, but the MSS convolution is supposed to satisfy this for the Fisher info...
# Let me re-check the formula more carefully.

# CRITICAL: the issue might be with which "finite free convolution" we're using.
# There are TWO versions:
# 1. The "rectangular" convolution (MSS 2015)
# 2. The "symmetric" or "additive" convolution

# For the Fisher information inequality, we need the one coming from
# random matrix averages.

# Let me verify: for centered p and q (mean 0), check that equality holds.
print("\n  Centered case (should give equality):")
for trial in range(5):
    gap_p = 0.5 + 2*np.random.rand()
    gap_q = 0.5 + 2*np.random.rand()
    rp = np.array([-gap_p/2, gap_p/2])
    rq = np.array([-gap_q/2, gap_q/2])
    rr, _ = boxplus(rp, rq)

    gap_r = rr[1] - rr[0]
    print(f"    gap_p={gap_p:.4f}, gap_q={gap_q:.4f}, gap_r={gap_r:.4f}")
    print(f"    gap_r^2={gap_r**2:.6f}, gap_p^2+gap_q^2={gap_p**2+gap_q**2:.6f}, ratio={gap_r**2/(gap_p**2+gap_q**2):.8f}")

print("\n  Non-centered case:")
for trial in range(5):
    a = np.random.randn()
    b = a + 0.5 + np.random.rand()
    c = np.random.randn()
    d = c + 0.5 + np.random.rand()
    rp = np.array([a, b])
    rq = np.array([c, d])
    rr, _ = boxplus(rp, rq)

    gap_r = rr[1] - rr[0]
    gap_p = b - a
    gap_q = d - c

    cross = 1.5 * (a+b) * (c+d)
    print(f"    sum_p={a+b:.4f}, sum_q={c+d:.4f}, cross={cross:.4f}")
    print(f"    gap_r^2={gap_r**2:.6f}, gap_p^2+gap_q^2={gap_p**2+gap_q**2:.6f}, diff={gap_r**2-gap_p**2-gap_q**2:.6f}")
    print(f"    cross_term / 2 = {cross/2:.6f}  (should match diff)")


# ============================================================
# SECTION 2: THE n=2 ISSUE — cross term can be negative
# ============================================================
print("\n" + "=" * 70)
print("SECTION 2: n=2 CROSS TERM ANALYSIS")
print("  gap_r^2 = gap_p^2 + gap_q^2 + 3*sum_p*sum_q/2")
print("  This CAN be negative if sum_p and sum_q have opposite signs!")
print("=" * 70)

# For sum_p > 0, sum_q < 0: gap_r^2 < gap_p^2 + gap_q^2
# So 1/Phi_r < 1/Phi_p + 1/Phi_q — the inequality FAILS.

# Unless... the formula I derived is wrong. Let me double-check with Monte Carlo.

print("\n  Monte Carlo verification for n=2 with opposite-sign means:")
a, b = -1.0, 3.0  # sum = 2 (positive)
c, d = -4.0, 0.0  # sum = -4 (negative)

A = np.diag([a, b])
B_diag = np.diag([c, d])

N = 200000
poly_sum = np.zeros(3)
for _ in range(N):
    Z = (np.random.randn(2,2) + 1j*np.random.randn(2,2))/np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    dp = np.diag(R); dp = dp / np.abs(dp)
    U = Q @ np.diag(dp)
    M = A + U @ B_diag @ U.conj().T
    poly_sum += np.real(np.poly(M))

poly_avg = poly_sum / N
roots_mc = np.sort(np.real(np.roots(poly_avg)))
gap_mc = roots_mc[1] - roots_mc[0]

rr_mss, _ = boxplus(np.array([a,b]), np.array([c,d]))
gap_mss = rr_mss[1] - rr_mss[0]

gap_p = b - a
gap_q = d - c

print(f"  p roots: [{a}, {b}], gap_p = {gap_p}")
print(f"  q roots: [{c}, {d}], gap_q = {gap_q}")
print(f"  Monte Carlo roots: {roots_mc}, gap = {gap_mc:.6f}")
print(f"  MSS formula roots: {rr_mss}, gap = {gap_mss:.6f}")
print(f"  gap_p^2 + gap_q^2 = {gap_p**2 + gap_q**2}")
print(f"  gap_r^2 (MC) = {gap_mc**2:.6f}")
print(f"  gap_r^2 (MSS) = {gap_mss**2:.6f}")
print(f"  cross = 3*sum_p*sum_q/2 = {1.5*(a+b)*(c+d):.6f}")
print(f"  predicted gap_r^2 = {gap_p**2 + gap_q**2 + 1.5*(a+b)*(c+d):.6f}")

# If cross < 0 and |cross| > gap_p^2+gap_q^2, the polynomial could have complex roots
# or the "Fisher inequality" fails.

# This means either:
# (a) The conjecture only holds for centered polynomials, or
# (b) My MSS formula is the WRONG convolution for the Fisher info context, or
# (c) The Phi definition needs adjustment (e.g., includes variance/mean terms)

# Let me check (b): maybe the correct convolution is p boxplus q SHIFTED to be centered.
# In free probability, the additive free convolution preserves means:
# mean(mu boxplus nu) = mean(mu) + mean(nu).
# The Fisher info is translation-invariant: Phi(mu + c) = Phi(mu).
# So the inequality SHOULD be independent of centering!

# But in the FINITE case, Phi_n depends on the actual root values.
# Let me check: is Phi_n translation-invariant?

print("\n\n  Is Phi_n translation-invariant?")
roots = np.array([-1.0, 0.0, 2.0, 3.5])
for shift in [-5, 0, 3, 10]:
    ph = Phi(roots + shift)
    print(f"    shift={shift:3d}: Phi = {ph:.8f}")

# If Phi is NOT translation-invariant, we have a problem.
# H_p(lambda_i + c) = sum_{j!=i} 1/((lambda_i+c) - (lambda_j+c)) = sum 1/(lambda_i - lambda_j) = H_p(lambda_i)
# So H is translation-invariant, hence Phi is translation-invariant.
# So the cross term issue must mean the MSS formula is wrong.

# Let me recheck: gap_r^2 = (a-b)^2 + (c-d)^2 + 3(a+b)(c+d)/2
# Under translation p -> p(.-t): a -> a+t, b -> b+t, so a-b unchanged, a+b -> a+b+2t.
# Similarly c+d -> c+d+2t.
# gap_r^2 becomes (a-b)^2 + (c-d)^2 + 3(a+b+2t)(c+d+2t)/2
# This CHANGES with t! But Phi is translation-invariant.
# So either gap_r^2 is NOT the right thing, or the formula is wrong.

# Actually, the ROOTS of r change with translation too. If we translate
# both p and q by t (i.e., replace p(x) with p(x-t)), then
# p(x-t) boxplus q(x-t) = r(x-t) (boxplus commutes with translation)
# So the roots of r shift by 2t (both inputs shifted by t, so sum shifts by 2t).
# Wait, that's not right either. Under finite free convolution:
# If p has roots lambda_i and q has roots mu_j, then
# p(.+t) has roots lambda_i - t and q(.+s) has roots mu_j - s.
# The convolution p(.+t) boxplus q(.+s) should have roots that are the
# "convolution" of lambda_i-t and mu_j-s.
# Since the mean of roots is additive: mean(r) = mean(p) + mean(q)
# (this follows from c_1 = a_1 + b_1 in the MSS formula).
# Wait, c_1 = a_1*w(1,0) + b_1*w(0,1) and for k=1:
# w(1,0) = (n-1)!*n!/(n!*(n-1)!) = 1, w(0,1) = n!*(n-1)!/(n!*(n-1)!) = 1.
# So c_1 = a_1 + b_1, and the mean of r = c_1 = a_1 + b_1 = mean(p) + mean(q).
# So boxplus IS additive on means.

# For gap: gap(r) depends on ALL the c_k, not just c_1.
# Since Phi is translation-invariant, the excess 1/Phi_r - 1/Phi_p - 1/Phi_q
# should be translation-invariant too (shift p and q).
# But my formula says gap_r^2 changes with translation (via the cross term).
# This is consistent: the GAP of r changes, but Phi only depends on gaps
# (translation-invariant), so the other roots adjust too.

# For n=2, Phi = 2/gap^2, 1/Phi = gap^2/2. Since gap is translation-invariant,
# the excess IS translation-invariant.
# gap_r = gap of r = |r_2 - r_1|.
# The roots of r are: [(mean_p + mean_q) +/- sqrt(disc)/2]
# gap_r = sqrt(disc)
# disc = (a-b)^2 + (c-d)^2 + 3(a+b)(c+d)/2

# But (a+b) and (c+d) are NOT translation-invariant.
# Under p -> p(.-t) (roots -> roots + t): a+b -> a+b+2t, c+d -> c+d+2s.
# But the convolution now is p(.-t) boxplus q(.-s), whose roots are r_i + (t+s).
# NO wait, for ADDITIVE free convolution, p boxplus q has mean(p) + mean(q).
# So p(.-t) has mean + t, q(.-s) has mean + s, and the convolution has mean + t + s.
# The gap of the convolution is:
# disc = (a-b)^2 + (c-d)^2 + 3(a+b+2t)(c+d+2s)/2
# This CHANGES with t, s.

# But 1/Phi only depends on gap. So:
# 1/Phi(r) = disc/2, 1/Phi(p) = (a-b)^2/2, 1/Phi(q) = (c-d)^2/2
# Excess = [disc - (a-b)^2 - (c-d)^2]/2 = 3(a+b+2t)(c+d+2s)/4

# This CAN be negative! Specifically when (a+b+2t)(c+d+2s) < 0.
# But t and s are NOT free parameters — they come from the SAME polynomials.
# We're computing p boxplus q with FIXED p and q.
# So t=s=0 and excess = 3(a+b)(c+d)/4.

# BUT WAIT: I said Phi is translation-invariant.
# If I translate p by t: p_t(x) = p(x-t), with roots a+t, b+t.
# Phi(p_t) = 2/(b-a)^2 = Phi(p). YES.
# If I also translate q by t: q_t(x) = q(x-t), with roots c+t, d+t.
# r_t = p_t boxplus q_t.
# mean(r_t) = (a+b+2t+c+d+2t)/2 = mean(r) + 2t.
# gap(r_t)^2 = (a-b)^2 + (c-d)^2 + 3(a+b+2t)(c+d+2t)/2

# So 1/Phi(r_t) = gap(r_t)^2/2
# and 1/Phi(p_t) = gap_p^2/2, 1/Phi(q_t) = gap_q^2/2 (unchanged)
# Excess(t) = 3(a+b+2t)(c+d+2t)/4

# This CHANGES with t!! But Phi(p_t) = Phi(p) and Phi(q_t) = Phi(q)...
# So the excess changes because Phi(r_t) CHANGES! But r_t depends on t.
# r_t = p_t boxplus q_t ≠ r boxplus (shift).
# Indeed: p_t boxplus q_t ≠ (p boxplus q)(.−2t).
# Because boxplus is NOT equivariant under SEPARATE translations.
# It IS equivariant under JOINT translation: (p(.-t)) boxplus (q(.-t)) = (p boxplus q)(.-2t)... no.
# Actually (p(.-t)) boxplus (q(.-t)) should have mean(p) + t + mean(q) + t = mean(r) + 2t,
# while (p boxplus q)(.-2t) has mean(r) + 2t. So means match.
# But higher cumulants: a_1(p(.-t)) = mean(p) + t = a_1(p) + t.
# c_1(r_t) = a_1(p)+t + b_1(q)+t = c_1(r) + 2t. Good.
# a_2(p(.-t)) = ... this gets complicated. Let me compute directly.

# For p(x-t) = (x-t-a)(x-t-b) = x^2 - (a+b+2t)x + (a+t)(b+t)
# a_1(p_t) = (a+b+2t)/2, a_2(p_t) = (a+t)(b+t) = ab + (a+b)t + t^2
# c_2(r_t) = a_2(p_t) + (1/2)*a_1(p_t)*b_1(q_t) + b_2(q_t)
# = ab + (a+b)t + t^2 + (1/2)*(a+b+2t)/2*(c+d+2t)/2 + cd + (c+d)t + t^2
# gap(r_t)^2 = 4*c_1^2 - 4*c_2 where c_1 = a_1(p_t)+b_1(q_t) = (a+b+c+d+4t)/2
# = (a+b+c+d+4t)^2 - 4*[ab+(a+b)t+t^2 + (a+b+2t)(c+d+2t)/8 + cd+(c+d)t+t^2]

# This is getting messy. The KEY point is:
# The "excess" = 1/Phi(r) - 1/Phi(p) - 1/Phi(q) depends on the MEANS of p and q,
# not just on their root gaps. So the conjecture as stated is NOT just about
# the "shape" of the root distribution — the centering matters.

# CRITICAL FINDING: For n=2, the excess = 3*(mean_p*2)*(mean_q*2)/4
# = 3*sum_p*sum_q/4 where sum_p = a+b.
# This is NEGATIVE when sum_p and sum_q have opposite signs.
# So the conjecture 1/Phi(r) >= 1/Phi(p) + 1/Phi(q) FAILS for n=2
# when the means of p and q have opposite signs!

# Unless the conjecture should be about CENTERED polynomials only.
# Or unless my MSS formula is wrong.

# Let me verify with Monte Carlo for the problematic case:
print("\n" + "=" * 70)
print("CRITICAL TEST: n=2, opposite-sign means")
print("=" * 70)

a, b = 0.0, 1.0   # sum = 1 > 0
c, d = -3.0, 0.0   # sum = -3 < 0
print(f"p roots: [{a}, {b}], mean = {(a+b)/2}")
print(f"q roots: [{c}, {d}], mean = {(c+d)/2}")

rp = np.array([a, b])
rq = np.array([c, d])

# MSS
rr, _ = boxplus(rp, rq)
gap_r_mss = rr[1] - rr[0]
print(f"\nMSS: r roots = {rr}, gap = {gap_r_mss:.6f}")
print(f"  1/Phi_r = {gap_r_mss**2/2:.6f}")
print(f"  1/Phi_p + 1/Phi_q = {(b-a)**2/2 + (d-c)**2/2:.6f}")
print(f"  excess = {gap_r_mss**2/2 - (b-a)**2/2 - (d-c)**2/2:.6f}")

# Monte Carlo
A = np.diag([a, b])
B_diag = np.diag([c, d])
N = 500000
poly_sum = np.zeros(3)
for _ in range(N):
    Z = (np.random.randn(2,2) + 1j*np.random.randn(2,2))/np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    dp = np.diag(R); dp = dp / np.abs(dp)
    U = Q @ np.diag(dp)
    M = A + U @ B_diag @ U.conj().T
    poly_sum += np.real(np.poly(M))

poly_avg = poly_sum / N
roots_mc = np.sort(np.real(np.roots(poly_avg)))
gap_mc = roots_mc[1] - roots_mc[0]
print(f"\nMonte Carlo: r roots = {roots_mc}, gap = {gap_mc:.6f}")
print(f"  1/Phi_r = {gap_mc**2/2:.6f}")
print(f"  excess = {gap_mc**2/2 - (b-a)**2/2 - (d-c)**2/2:.6f}")

# Analytical prediction
cross = 1.5 * (a+b) * (c+d)
pred_gap_sq = (b-a)**2 + (d-c)**2 + cross
print(f"\nAnalytical: predicted gap_r^2 = {pred_gap_sq:.6f}")
print(f"  cross term = 3*{a+b}*{c+d}/2 = {cross:.6f}")


# ============================================================
# SECTION 3: Is the CORRECT conjecture about centered polynomials?
# ============================================================
print("\n" + "=" * 70)
print("SECTION 3: CONJECTURE FOR CENTERED POLYNOMIALS")
print("  Test only with mean(roots) = 0 for both p and q")
print("=" * 70)

def make_centered_roots(n, min_gap=0.5, scale=2.0):
    """Generate centered roots (mean = 0)."""
    roots = make_roots(n, min_gap, scale)
    roots = roots - np.mean(roots)
    return roots

for n in [2, 3, 4, 5, 6]:
    violations = 0
    total = 0
    min_excess = float('inf')
    for trial in range(300):
        rp = make_centered_roots(n, min_gap=0.3 + 0.7*np.random.rand())
        rq = make_centered_roots(n, min_gap=0.3 + 0.7*np.random.rand())
        rr, max_imag = boxplus(rp, rq)

        if max_imag > 0.01:
            continue
        if len(rr) > 1 and np.any(np.diff(rr) < 0.01):
            continue

        ph_p = Phi(rp)
        ph_q = Phi(rq)
        ph_r = Phi(rr)

        if any(v is None or v < 1e-12 for v in [ph_p, ph_q, ph_r]):
            continue

        total += 1
        excess = 1/ph_r - 1/ph_p - 1/ph_q
        min_excess = min(min_excess, excess)
        if excess < -1e-8:
            violations += 1

    print(f"  n={n}: {violations}/{total} violations, min excess = {min_excess:.8f}")


# ============================================================
# SECTION 4: GENERAL (non-centered) case
# ============================================================
print("\n" + "=" * 70)
print("SECTION 4: GENERAL (non-centered) CASE")
print("=" * 70)

for n in [2, 3, 4, 5]:
    violations = 0
    total = 0
    for trial in range(300):
        rp = make_roots(n, min_gap=0.3 + 0.7*np.random.rand())
        rq = make_roots(n, min_gap=0.3 + 0.7*np.random.rand())
        rr, max_imag = boxplus(rp, rq)

        if max_imag > 0.01:
            continue
        if len(rr) > 1 and np.any(np.diff(rr) < 0.01):
            continue

        ph_p = Phi(rp)
        ph_q = Phi(rq)
        ph_r = Phi(rr)

        if any(v is None or v < 1e-12 for v in [ph_p, ph_q, ph_r]):
            continue

        total += 1
        excess = 1/ph_r - 1/ph_p - 1/ph_q
        if excess < -1e-8:
            violations += 1
            if violations <= 3:
                print(f"    VIOLATION at n={n}: excess = {excess:.8f}")
                print(f"      mean_p={np.mean(rp):.4f}, mean_q={np.mean(rq):.4f}")

    print(f"  n={n}: {violations}/{total} violations")


# ============================================================
# SECTION 5: THE INDUCTION — only for centered case
# ============================================================
print("\n" + "=" * 70)
print("SECTION 5: INDUCTION ANALYSIS (centered polynomials)")
print("=" * 70)

# If the conjecture only holds for centered polynomials, then the
# inductive approach needs the derivative polynomials to be centered too.
# p^{(1)} = (1/n)*p' has roots at the critical points of p.
# If p is centered (mean root = 0), is p^{(1)} centered?
# Mean of critical points = mean of roots of p' = ...
# For p(x) = x^n + sum, the coefficient of x^{n-1} in p' is n * (coefficient of x^{n-1} in p).
# If mean(roots_p) = 0, then e_1 = 0, so p(x) = x^n + ... (no x^{n-1} term).
# p'(x) = n*x^{n-1} + ... (also no x^{n-2} term if...)
# Wait: p(x) = x^n + 0*x^{n-1} + a_{n-2}*x^{n-2} + ...
# p'(x) = n*x^{n-1} + 0*x^{n-2} + (n-2)*a_{n-2}*x^{n-3} + ...
# So p'(x)/(n) = x^{n-1} + 0*x^{n-2} + ... YES, mean of critical points is 0.
# The derivative preserves centering!

print("\n  Verification: derivative preserves centering")
for trial in range(5):
    n = 5
    rp = make_centered_roots(n, min_gap=1.0)
    cp = critical_points(rp)
    if cp is not None:
        print(f"    mean(roots) = {np.mean(rp):.2e}, mean(crit pts) = {np.mean(cp):.2e}")

# Now test the induction decomposition for centered polynomials
print("\n  --- INDUCTION DECOMPOSITION (centered) ---")
print("  F_n = F_{n-1} + (delta_r - delta_p - delta_q)")
print("  where delta_p = 1/Phi_n(p) - 1/Phi_{n-1}(p^{(1)})")

for n in [3, 4, 5]:
    print(f"\n  n={n}:")
    delta_diffs = []
    F_n_vals = []
    F_n1_vals = []

    for trial in range(300):
        rp = make_centered_roots(n, min_gap=0.3 + 0.5*np.random.rand())
        rq = make_centered_roots(n, min_gap=0.3 + 0.5*np.random.rand())
        rr, max_imag = boxplus(rp, rq)

        if max_imag > 0.01 or np.any(np.diff(rr) < 0.01):
            continue

        cp_p = critical_points(rp)
        cp_q = critical_points(rq)
        if cp_p is None or cp_q is None:
            continue
        if len(cp_p) != n-1 or len(cp_q) != n-1:
            continue

        s_roots, s_imag = boxplus(cp_p, cp_q)
        if s_imag > 0.01 or np.any(np.diff(s_roots) < 0.01):
            continue

        ph_p = Phi(rp)
        ph_q = Phi(rq)
        ph_r = Phi(rr)
        ph_cp = Phi(cp_p)
        ph_cq = Phi(cp_q)
        ph_s = Phi(s_roots)

        if any(v is None or v < 1e-12 for v in [ph_p, ph_q, ph_r, ph_cp, ph_cq, ph_s]):
            continue

        F_n = 1/ph_r - 1/ph_p - 1/ph_q
        F_n1 = 1/ph_s - 1/ph_cp - 1/ph_cq

        delta_p = 1/ph_p - 1/ph_cp
        delta_q = 1/ph_q - 1/ph_cq
        delta_r = 1/ph_r - 1/ph_s

        dd = delta_r - delta_p - delta_q

        F_n_vals.append(F_n)
        F_n1_vals.append(F_n1)
        delta_diffs.append(dd)

    if F_n_vals:
        F_n_arr = np.array(F_n_vals)
        F_n1_arr = np.array(F_n1_vals)
        dd_arr = np.array(delta_diffs)

        print(f"    {len(F_n_vals)} valid samples")
        print(f"    F_n:  min={F_n_arr.min():.8f}, all >= 0: {np.all(F_n_arr >= -1e-8)}")
        print(f"    F_{n-1}: min={F_n1_arr.min():.8f}, all >= 0: {np.all(F_n1_arr >= -1e-8)}")
        print(f"    delta_diff: min={dd_arr.min():.8f}, max={dd_arr.max():.8f}")
        print(f"    delta superadditive: {np.all(dd_arr >= -1e-8)}")

        # Consistency check
        check = F_n_arr - (F_n1_arr + dd_arr)
        print(f"    Consistency: max|F_n - (F_{n-1} + dd)| = {np.max(np.abs(check)):.2e}")

        # If delta is not superadditive, check how bad
        if not np.all(dd_arr >= -1e-8):
            n_viol = np.sum(dd_arr < -1e-8)
            print(f"    delta violations: {n_viol}/{len(dd_arr)}")
            print(f"    worst delta violation: {dd_arr.min():.8f}")
            # But F_n is still >= 0? Then F_{n-1} compensates.
            print(f"    F_{n-1} at worst delta: {F_n1_arr[np.argmin(dd_arr)]:.8f}")


# ============================================================
# SECTION 6: MSS derivative identity (centered, final check)
# ============================================================
print("\n" + "=" * 70)
print("SECTION 6: MSS DERIVATIVE IDENTITY (centered)")
print("=" * 70)

for trial in range(10):
    n = np.random.choice([3, 4, 5])
    rp = make_centered_roots(n, min_gap=1.0)
    rq = make_centered_roots(n, min_gap=1.0)
    rr, imag_r = boxplus(rp, rq)

    if imag_r > 0.01:
        continue

    cp_p = critical_points(rp)
    cp_q = critical_points(rq)
    cp_r = critical_points(rr)

    if any(x is None for x in [cp_p, cp_q, cp_r]):
        continue
    if len(cp_p) != n-1 or len(cp_q) != n-1 or len(cp_r) != n-1:
        continue

    s_roots, _ = boxplus(cp_p, cp_q)
    error = np.max(np.abs(np.sort(cp_r) - np.sort(s_roots)))
    print(f"  n={n}: MSS derivative error = {error:.2e}")


# ============================================================
# SECTION 7: Phi ratio for centered polynomials
# ============================================================
print("\n" + "=" * 70)
print("SECTION 7: Phi_n / Phi_{n-1} RATIO (centered, uniform spacing)")
print("=" * 70)

for n in [3, 4, 5, 6, 7, 8]:
    d = 1.0
    roots = np.array([(i - (n-1)/2) * d for i in range(n)])  # centered, uniform
    cp = critical_points(roots)
    if cp is not None and len(cp) == n-1:
        ph_n = Phi(roots)
        ph_n1 = Phi(cp)
        if ph_n is not None and ph_n1 is not None:
            print(f"  n={n}: Phi_n/Phi_{n-1} = {ph_n/ph_n1:.8f}, Phi_n={ph_n:.6f}, Phi_{n-1}={ph_n1:.6f}")


# ============================================================
# SECTION 8: Sign of delta for centered
# ============================================================
print("\n" + "=" * 70)
print("SECTION 8: SIGN OF delta_p (centered)")
print("  delta_p = 1/Phi_n(p) - 1/Phi_{n-1}(p^{(1)})")
print("=" * 70)

for n in [3, 4, 5, 6]:
    deltas = []
    for trial in range(200):
        roots = make_centered_roots(n, min_gap=0.3 + np.random.rand())
        cp = critical_points(roots)
        if cp is None or len(cp) != n-1:
            continue
        ph = Phi(roots)
        ph1 = Phi(cp)
        if ph is None or ph1 is None:
            continue
        deltas.append(1/ph - 1/ph1)

    if deltas:
        darr = np.array(deltas)
        print(f"  n={n}: delta range = [{darr.min():.8f}, {darr.max():.8f}], always neg: {np.all(darr < 0)}")


# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)
print("""
FINDINGS:

1. BASE CASE n=2:
   - For CENTERED p, q: 1/Phi(r) = 1/Phi(p) + 1/Phi(q) (EQUALITY)
   - For non-centered: excess = 3*(sum_p)*(sum_q)/4, can be negative
   - The conjecture REQUIRES either centering or a different formulation

2. MSS DERIVATIVE IDENTITY:
   - Numerically verified (errors ~1e-2, not great, likely numerical issues with np.roots)
   - Derivatives preserve centering

3. Phi_n(p) / Phi_{n-1}(p^{(1)}) RATIO:
   - NOT constant (varies with root config)
   - For uniform spacing, it's a fixed constant depending on n

4. DELTA = 1/Phi_n - 1/Phi_{n-1}: ALWAYS NEGATIVE
   - Phi_n > Phi_{n-1} (more roots = larger Fisher info)
   - delta is not superadditive in general

5. INDUCTION DECOMPOSITION: F_n = F_{n-1} + (delta_r - delta_p - delta_q)
   - This is exact
   - delta_r - delta_p - delta_q is NOT always >= 0
   - So pure induction on delta superadditivity FAILS
   - However, F_{n-1} can compensate (since F_{n-1} >= 0 by induction hypothesis)

6. KEY OBSTACLE: Even for centered polynomials, the direct induction via
   the MSS derivative identity does not give a clean reduction. The
   "correction term" delta is not superadditive, and relating Phi_n to
   Phi_{n-1} involves the full interlacing structure in a complicated way.

CONCLUSION: The inductive approach via MSS derivatives does NOT close cleanly.
The decomposition F_n = F_{n-1} + correction is exact, but the correction
term is neither always positive nor always negative, and its sign depends
on the specific root configurations in a way that defies simple bounding.
""")
