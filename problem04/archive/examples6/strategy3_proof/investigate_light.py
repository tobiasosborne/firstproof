"""
Light version of structural investigation - key findings only.
"""

import numpy as np
from math import factorial, comb
from itertools import combinations

def elem_sym_poly(roots, k):
    n = len(roots)
    if k == 0: return 1.0
    if k > n: return 0.0
    return sum(np.prod([roots[i] for i in subset])
               for subset in combinations(range(n), k))

def poly_coeffs_from_roots(roots):
    n = len(roots)
    ek = [elem_sym_poly(roots, k) for k in range(n+1)]
    return [(-1)**k * ek[k] for k in range(n+1)]

def boxplus_mss(roots_p, roots_q):
    n = len(roots_p)
    a = poly_coeffs_from_roots(roots_p)
    b = poly_coeffs_from_roots(roots_q)
    c = np.zeros(n+1)
    for k in range(n+1):
        for i in range(k+1):
            j = k - i
            if i <= n and j <= n:
                coeff = factorial(n-i) * factorial(n-j) / (factorial(n) * factorial(n-k))
                c[k] += coeff * a[i] * b[j]
    return np.sort(np.real(np.roots(c))), c

def H_values(roots):
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H

def subordination(roots_p, roots_r):
    n = len(roots_p)
    roots_p_sorted = np.sort(roots_p)
    roots_r_sorted = np.sort(roots_r)
    sigma = np.arange(n)
    omega1_at_roots = roots_p_sorted[sigma]

    p_coeffs = np.poly(roots_p_sorted)
    r_coeffs = np.poly(roots_r_sorted)
    p_prime = np.polyder(p_coeffs)
    r_prime = np.polyder(r_coeffs)
    p_dprime = np.polyder(p_prime)
    r_dprime = np.polyder(r_prime)

    omega1_prime = np.zeros(n)
    omega1_dprime = np.zeros(n)

    for k in range(n):
        z0 = roots_r_sorted[k]
        w0 = omega1_at_roots[k]
        rp = np.polyval(r_prime, z0)
        rpp = np.polyval(r_dprime, z0)
        pp = np.polyval(p_prime, w0)
        ppp = np.polyval(p_dprime, w0)

        F_z = -pp * rp
        F_w = rp * pp
        omega1_prime[k] = -F_z / F_w

        F_zz = -pp * rpp
        F_zw = rpp * pp - ppp * rp
        F_ww = rp * ppp
        w1 = omega1_prime[k]
        omega1_dprime[k] = -(F_zz + 2*F_zw*w1 + F_ww*w1**2) / F_w

    return omega1_at_roots, omega1_prime, omega1_dprime, sigma


np.random.seed(42)

# ================================================================
# Test AB/h^4 as p -> q (limiting behavior)
# ================================================================
print("="*70)
print("AB/h^4 AS p -> q (equality limit)")
print("="*70)

base_roots = np.array([-2.0, 0.0, 2.0])
for epsilon in [1.0, 0.5, 0.2, 0.1, 0.05, 0.01, 0.001]:
    roots_p = base_roots + epsilon * np.array([0.1, -0.2, 0.1])
    roots_q = base_roots - epsilon * np.array([0.1, -0.2, 0.1])

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        roots_r = np.sort(np.real(roots_r))
        if np.any(np.diff(roots_r) < 0.001):
            continue

        Ph_p = sum(H_values(roots_p)**2)
        Ph_q = sum(H_values(roots_q)**2)
        Ph_r = sum(H_values(roots_r)**2)

        A = Ph_p - Ph_r
        B = Ph_q - Ph_r
        h4 = Ph_r**2

        if h4 > 0 and A > 0 and B > 0:
            print(f"  eps={epsilon:.4f}: AB/h4 = {A*B/h4:.8f}")
    except:
        pass

# ================================================================
# Key analysis: (alpha+beta)/h ratio structure
# ================================================================
print("\n" + "="*70)
print("(alpha+beta)_k / h_k RATIO STRUCTURE")
print("="*70)

for trial in range(15):
    n = np.random.choice([3, 4, 5])
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        raw = np.roots(np.poly(roots_r))
        if np.any(np.abs(np.imag(raw)) > 0.01):
            continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01):
            continue

        _, _, omega1_pp, sigma = subordination(roots_p, roots_r)
        _, _, omega2_pp, tau = subordination(roots_q, roots_r)

        h = H_values(roots_r)
        alpha = omega1_pp / 2
        beta = omega2_pp / 2

        # Check omega1' = 1
        if np.max(np.abs(omega1_pp - 2*alpha)) > 1e-10:
            continue

        if np.min(np.abs(h)) > 0.01:
            ratio_apb = (alpha + beta) / h
            ratio_a = alpha / h
            ratio_b = beta / h
            print(f"n={n}: alpha/h = {ratio_a}, beta/h = {ratio_b}")
            print(f"       (a+b)/h = {ratio_apb}")
            print(f"       <h,alpha>={np.dot(h,alpha):.6f}, <h,beta>={np.dot(h,beta):.6f}")
    except:
        pass


# ================================================================
# THE KEY INSIGHT: <h,alpha> = sum_k h_k * (u_k - h_k)
# Try to prove this >= 0 via the Cauchy-Schwarz inequality or
# by studying the structure of H_p(lambda_k) vs H_r(nu_k).
# ================================================================
print("\n\n" + "="*70)
print("KEY: COMPONENT-WISE ANALYSIS h_k * alpha_k = h_k * (u_k - h_k)")
print("="*70)

# h_k = sum_{j!=k} 1/(nu_k - nu_j)
# u_k = sum_{j!=k} 1/(lambda_k - lambda_j)
# alpha_k = u_k - h_k = sum_{j!=k} [1/(lambda_k-lambda_j) - 1/(nu_k-nu_j)]

# For the MSS convolution, omega_1(nu_k) = lambda_k (order-preserving subordination).
# omega_1 is increasing on the real line outside the poles.

# KEY STRUCTURAL PROPERTY:
# The roots of r = p boxplus_n q satisfy interlacing with respect to both p and q.
# Specifically: nu_k = lambda_k + mu_k - c_k for some corrections c_k.

# Actually, let me check: what is the relationship between the gaps?
# If gap_r_k = nu_{k+1} - nu_k, gap_p_k = lambda_{k+1} - lambda_k,
# and gap_q_k = mu_{k+1} - mu_k, does gap_r_k relate to gap_p_k + gap_q_k?

print("\nGap structure analysis:")
for trial in range(10):
    n = 4
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        roots_r = np.sort(np.real(roots_r))
        if np.any(np.diff(roots_r) < 0.01):
            continue

        gaps_p = np.diff(roots_p)
        gaps_q = np.diff(roots_q)
        gaps_r = np.diff(roots_r)

        # Check: gap_r^2 vs gap_p^2 + gap_q^2 (componentwise)
        gap_sq_sum = gaps_p**2 + gaps_q**2
        gap_r_sq = gaps_r**2

        print(f"  gap_r^2 = {gap_r_sq}")
        print(f"  gap_p^2 + gap_q^2 = {gap_sq_sum}")
        print(f"  ratio = {gap_r_sq / gap_sq_sum}")
        print()
    except:
        pass


# ================================================================
# THE REAL TEST: Can we prove <h,alpha> >= 0 from the Herglotz structure?
# ================================================================
print("="*70)
print("HERGLOTZ STRUCTURE OF OMEGA_1")
print("="*70)

# omega_1 is a Herglotz function: omega_1: C+ -> C+, omega_1(z) ~ z as z -> infinity.
# On the real line, omega_1 has a representation (Nevanlinna/Herglotz):
# omega_1(z) = az + b + integral [1/(t-z) - t/(1+t^2)] dmu(t)
# where a >= 0, b in R, mu is a positive Borel measure.
# Since omega_1(z) ~ z + O(1/z) as z -> infinity, a = 1 and integral dmu(t) = c.
# omega_1'(z) = 1 + integral 1/(t-z)^2 dmu(t) >= 1 for Im(z) = 0.

# Wait: omega_1'(z) = 1 + integral 1/(t-z)^2 dmu(t).
# At z = nu_k (real): omega_1'(nu_k) = 1 + integral 1/(t-nu_k)^2 dmu(t) = 1.
# So integral 1/(t-nu_k)^2 dmu(t) = 0.
# Since the integrand is POSITIVE, this means mu has NO support!
# That would make omega_1(z) = z + b, a linear function.
# But omega_1 is algebraic of degree n, so it can't be linear.

# CONTRADICTION? Let me reconsider.
# Actually, omega_1 might have POLES on the real line (between roots of r).
# The Herglotz representation for functions with poles is:
# omega_1(z) = az + b + sum_j c_j/(t_j - z) + integral [1/(t-z) - t/(1+t^2)] dmu_ac(t)
# where c_j > 0 are the residues at the poles t_j.

# omega_1'(z) = 1 + sum_j c_j/(t_j - z)^2 + integral 1/(t-z)^2 dmu_ac(t)
# At z = nu_k: omega_1'(nu_k) = 1 + sum_j c_j/(t_j - nu_k)^2 = 1
# So sum_j c_j/(t_j - nu_k)^2 = 0, but each term is POSITIVE.
# This again implies c_j = 0 for all j, contradicting that omega_1 has poles.

# RESOLUTION: omega_1 as defined by the subordination IS NOT globally Herglotz
# on the real line. It maps C+ to C+, but on the real line it may have poles.
# The representation is:
# omega_1(z) = z + c + sum_j c_j/(t_j - z)  for z in C+
# where c_j > 0, t_j in R (poles on the real line).

# Then omega_1'(z) = 1 + sum_j c_j/(t_j - z)^2
# At z = nu_k: 1 = 1 + sum_j c_j/(t_j - nu_k)^2
# => sum_j c_j/(t_j - nu_k)^2 = 0, impossible since c_j > 0.

# UNLESS c_j = 0 for all j. But then omega_1(z) = z + c, which means
# lambda_k = nu_k + c for all k.

# Hmm, this can't be right. Let me reconsider.

# Actually, for the FINITE polynomial case, omega_1 is a RATIONAL function
# (algebraic of degree n), not a true Herglotz function in the traditional sense.
# The Herglotz property (maps C+ to C+) holds, but the representation is different.

# For a degree-n rational function mapping C+ to C+:
# omega_1(z) = z + sum_j c_j/(t_j - z) where c_j > 0 and t_j are n-1 poles.

# omega_1'(z) = 1 + sum_j c_j/(t_j - z)^2
# This is > 1 for all z not at a pole.
# omega_1'(nu_k) should be > 1, but we computed it as = 1.

# WAIT: We computed omega_1'(nu_k) = 1 from the implicit function theorem.
# Let me verify this isn't an error.

# From F(z,w) = r'(z)p(w) - p'(w)r(z) = 0, at (nu_k, lambda_k):
# F_z = r''(nu_k)*p(lambda_k) - p'(lambda_k)*r'(nu_k) = 0 - p'(lambda_k)*r'(nu_k) = -p'(lambda_k)*r'(nu_k)
# F_w = r'(nu_k)*p'(lambda_k) - p''(lambda_k)*r(nu_k) = r'(nu_k)*p'(lambda_k) - 0 = r'(nu_k)*p'(lambda_k)
# omega_1' = -F_z/F_w = p'(lambda_k)*r'(nu_k) / (r'(nu_k)*p'(lambda_k)) = 1

# This is CORRECT: omega_1'(nu_k) = 1 at roots of r.

# But if omega_1(z) = z + sum c_j/(t_j-z), then omega_1'(z) > 1 everywhere on R \ {t_j}.
# Since omega_1'(nu_k) = 1, this representation must be wrong.

# RESOLUTION: omega_1 is NOT of the form z + sum c_j/(t_j-z) when c_j > 0.
# Actually, for the FINITE free convolution, omega_1 might NOT be a classical
# Herglotz function. It maps C+ to C+ (a finite version of this), but the
# representation omega_1(z) = z + c + integral... may have NEGATIVE c_j terms.

# Actually, in the FINITE case, the subordination function is the COMPOSITIONAL
# INVERSE of a certain function. It's algebraic but not necessarily a Stieltjes
# transform. The omega_1'(nu_k) = 1 result is specific to the finite free case.

# This is a crucial subtlety. In the INFINITE free probability case:
# omega_1'(x) > 0 on the support of mu_r (free analogue), and
# omega_1'(x) = 1 + integral ... which could be > 1 or = 1 depending on regularity.

# For the FINITE case: omega_1 is rational/algebraic, and omega_1'(nu_k) = 1 always.

# ================================================================
# ALTERNATIVE APPROACH: Direct identity for <h,alpha>
# ================================================================
print("\n" + "="*70)
print("DIRECT IDENTITY FOR <h,alpha>")
print("="*70)

# <h,alpha> = sum_k h_k * alpha_k where alpha_k = omega_1''(nu_k)/2
# = (1/2) * sum_k H_r(nu_k) * omega_1''(nu_k)

# omega_1''(nu_k) = 2*(H_p(lambda_k) - H_r(nu_k))

# So <h,alpha> = sum_k H_r(nu_k) * (H_p(lambda_k) - H_r(nu_k))
# = sum_k H_r(nu_k) * H_p(lambda_k) - Phi_r

# = <h, u> - ||h||^2 where u_k = H_p(lambda_{sigma(k)})

# The question: why is <h,u> >= ||h||^2?

# Note that both h and u are H-transform vectors evaluated at different root sets.
# h_k = H_r(nu_k) = sum_{j!=k} 1/(nu_k - nu_j)
# u_k = H_p(lambda_k) = sum_{j!=k} 1/(lambda_k - lambda_j)

# One approach: write h_k = sum_{j!=k} 1/(nu_k - nu_j) and consider the matrix
# M_{kj} = 1/(nu_k - nu_j) for j != k, M_{kk} = 0
# Then h = M * 1 (matrix times all-ones vector... no, it's the row sums)

# Actually h_k = sum_j M_{kj} where M is the Cauchy-like matrix.

# This is getting complex. Let me check numerically one more key thing:
# whether <h,alpha> = some function of the root gaps.

print("\n<h,alpha> in terms of Phi values:")
for trial in range(20):
    n = np.random.choice([3, 4, 5])
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        raw = np.roots(np.poly(roots_r))
        if np.any(np.abs(np.imag(raw)) > 0.01):
            continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01):
            continue

        _, _, omega1_pp, sigma = subordination(roots_p, roots_r)
        _, _, omega2_pp, tau = subordination(roots_q, roots_r)

        h = H_values(roots_r)
        alpha = omega1_pp / 2
        beta = omega2_pp / 2

        h_alpha = np.dot(h, alpha)
        h_beta = np.dot(h, beta)
        Ph_r = np.dot(h, h)
        Ph_p = np.dot(h + alpha, h + alpha)
        Ph_q = np.dot(h + beta, h + beta)

        A = Ph_p - Ph_r
        B = Ph_q - Ph_r

        # <h,alpha> = (A - ||alpha||^2) / 2
        # ||alpha||^2 = A - 2*<h,alpha>
        # So <h,alpha> = <h,u> - ||h||^2

        # Check: is <h,alpha> = (A - ||alpha||^2)/2?
        norm_a2 = np.dot(alpha, alpha)
        check = (A - norm_a2) / 2
        # print(f"  <h,a>={h_alpha:.6f}, (A-||a||^2)/2={check:.6f}, match={np.isclose(h_alpha, check)}")

        # Check: 2*<h,alpha> + ||alpha||^2 = A
        # <h,alpha> ranges from 0 to ... what?
        # By CS: <h,alpha>^2 <= ||h||^2 * ||alpha||^2 = Ph_r * norm_a2
        # So <h,alpha> <= sqrt(Ph_r * norm_a2)
        # And A = 2<h,alpha> + norm_a2 <= 2*sqrt(Ph_r*norm_a2) + norm_a2

        # Alternative: A = ||u||^2 - ||h||^2 = ||u-h||^2 + 2<h,u-h> - ... no
        # A = ||u||^2 - ||h||^2 = (||u|| - ||h||)(||u|| + ||h||)

        # From u = h + alpha: ||u|| = ||h + alpha|| >= ||h|| (if <h,alpha> >= 0)
        # And ||u|| <= ||h|| + ||alpha||
        # So ||h|| <= ||u|| <= ||h|| + ||alpha||
        # A = ||u||^2 - ||h||^2 >= 0 and A <= (||h||+||alpha||)^2 - ||h||^2 = 2||h||*||alpha|| + ||alpha||^2

    except:
        pass

# ================================================================
# THE MONOTONE COUPLING APPROACH
# ================================================================
print("\n" + "="*70)
print("MONOTONE COUPLING: lambda_k > nu_k or lambda_k < nu_k?")
print("="*70)

# Since omega_1(nu_k) = lambda_k and omega_1(z) ~ z + c as z -> infty,
# lambda_k = nu_k + c + [correction terms].
# If c > 0, then lambda_k > nu_k for most k.
# The value of c is related to the mean shift:
# sum lambda_k = sum nu_k + n*c + [sum of corrections]
# mean(lambda) - mean(nu) = ...

for trial in range(10):
    n = 4
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        roots_r = np.sort(np.real(roots_r))
        if np.any(np.diff(roots_r) < 0.01):
            continue

        shifts = roots_p - roots_r
        mean_p = np.mean(roots_p)
        mean_q = np.mean(roots_q)
        mean_r = np.mean(roots_r)

        print(f"  mean_p={mean_p:.4f}, mean_q={mean_q:.4f}, mean_r={mean_r:.4f}")
        print(f"  mean_p - mean_r = {mean_p - mean_r:.4f}, mean_q - mean_r = {mean_q - mean_r:.4f}")
        print(f"  lambda_k - nu_k = {shifts}")
        print(f"  Note: mean(lambda) - mean(nu) = {mean_p - mean_r:.6f} = -mean(q)_shift")
        print()
    except:
        pass
