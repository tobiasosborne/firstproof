"""
Deep structural analysis of the vectors alpha, beta, h.

Now that we know <h,alpha> >= 0 and AB >= h^4 numerically,
we investigate the STRUCTURE that makes this true.

Key questions:
1. What is the relationship between alpha and h? Are they parallel?
2. Is there a monotonicity argument for <h,alpha> >= 0?
3. What is the tightest bound: is AB = h^4 + something positive?
4. Can we express AB - h^4 in terms of known quantities?
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

def Phi_n(roots):
    return np.sum(H_values(roots)**2)

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
    p_tprime = np.polyder(p_dprime)
    r_tprime = np.polyder(r_dprime)

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

print("="*70)
print("DEEP STRUCTURAL ANALYSIS")
print("="*70)

# Collect detailed data
data = []
for trial in range(500):
    n = np.random.choice([3, 4, 5, 6])

    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5

    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5

    try:
        roots_r, c_vec = boxplus_mss(roots_p, roots_q)
        raw_roots = np.roots(c_vec)
        if np.any(np.abs(np.imag(raw_roots)) > 0.01):
            continue
        roots_r = np.sort(np.real(raw_roots))
        gaps = np.diff(roots_r)
        if np.any(gaps < 0.01):
            continue

        _, omega1_p, omega1_pp, sigma = subordination(roots_p, roots_r)
        _, omega2_p, omega2_pp, tau = subordination(roots_q, roots_r)

        h = H_values(roots_r)
        alpha = omega1_pp / 2
        beta = omega2_pp / 2
        u = H_values(roots_p)[sigma]
        v = H_values(roots_q)[tau]

        cr_err = max(np.max(np.abs(u - (h + alpha))), np.max(np.abs(v - (h + beta))))
        if cr_err > 0.01:
            continue

        h_alpha = np.dot(h, alpha)
        h_beta = np.dot(h, beta)
        norm_h2 = np.dot(h, h)
        norm_alpha2 = np.dot(alpha, alpha)
        norm_beta2 = np.dot(beta, beta)
        alpha_beta = np.dot(alpha, beta)

        A = 2*h_alpha + norm_alpha2
        B = 2*h_beta + norm_beta2

        # Cosine of angle between h and alpha
        if norm_alpha2 > 1e-15:
            cos_ha = h_alpha / np.sqrt(norm_h2 * norm_alpha2)
        else:
            cos_ha = 1.0
        if norm_beta2 > 1e-15:
            cos_hb = h_beta / np.sqrt(norm_h2 * norm_beta2)
        else:
            cos_hb = 1.0

        # Ratio <h,alpha>/||h||^2
        ratio_ha = h_alpha / norm_h2 if norm_h2 > 1e-15 else 0

        # AB decomposition
        AB = A * B
        h4 = norm_h2**2
        AB_h4 = AB - h4

        # <u,v> = <h+alpha, h+beta> = ||h||^2 + <h,alpha> + <h,beta> + <alpha,beta>
        uv = norm_h2 + h_alpha + h_beta + alpha_beta

        data.append({
            'n': n, 'h_alpha': h_alpha, 'h_beta': h_beta,
            'norm_h2': norm_h2, 'norm_alpha2': norm_alpha2, 'norm_beta2': norm_beta2,
            'alpha_beta': alpha_beta,
            'A': A, 'B': B, 'AB': AB, 'h4': h4, 'AB_h4': AB_h4,
            'cos_ha': cos_ha, 'cos_hb': cos_hb,
            'ratio_ha': ratio_ha,
            'uv': uv,
            'h': h, 'alpha': alpha, 'beta': beta,
        })
    except:
        pass

print(f"Valid trials: {len(data)}")

# Analysis 1: Is h parallel to alpha?
print("\n--- Angle between h and alpha ---")
cos_ha_vals = [d['cos_ha'] for d in data]
print(f"  min cos(h,alpha) = {min(cos_ha_vals):.6f}")
print(f"  max cos(h,alpha) = {max(cos_ha_vals):.6f}")
print(f"  mean cos(h,alpha) = {np.mean(cos_ha_vals):.6f}")
print(f"  Parallel (cos > 0.99): {sum(1 for c in cos_ha_vals if c > 0.99)}/{len(cos_ha_vals)}")
print(f"  Nearly parallel (cos > 0.9): {sum(1 for c in cos_ha_vals if c > 0.9)}/{len(cos_ha_vals)}")

# Analysis 2: Ratio <h,alpha>/||h||^2
print("\n--- Ratio <h,alpha>/||h||^2 ---")
ratios = [d['ratio_ha'] for d in data]
print(f"  min = {min(ratios):.6f}")
print(f"  max = {max(ratios):.6f}")
print(f"  mean = {np.mean(ratios):.6f}")

# Analysis 3: <h,alpha> component-wise analysis
print("\n--- Component-wise h_k * alpha_k signs ---")
all_positive = 0
mixed = 0
for d in data:
    h, alpha = d['h'], d['alpha']
    products = h * alpha
    if np.all(products >= -1e-10):
        all_positive += 1
    else:
        mixed += 1
print(f"  All h_k*alpha_k >= 0: {all_positive}/{len(data)}")
print(f"  Mixed signs: {mixed}/{len(data)}")

# Analysis 4: What makes AB >= h^4 work?
# AB = (2<h,a> + ||a||^2)(2<h,b> + ||b||^2)
# = 4<h,a><h,b> + 2<h,a>||b||^2 + 2||a||^2<h,b> + ||a||^2||b||^2
# We need this >= ||h||^4
# If <h,a> >= 0 and <h,b> >= 0, then each of the 4 terms is >= 0.

# By Cauchy-Schwarz: <h,a> <= ||h||*||a||, so
# 4<h,a><h,b> <= 4||h||^2||a||*||b||
# 2<h,a>||b||^2 <= 2||h||*||a||*||b||^2
# etc. These don't obviously sum to >= ||h||^4.

# Try: Can we show 4<h,a><h,b> >= ||h||^4?
# This requires <h,a><h,b> >= ||h||^4/4.
print("\n--- Can 4<h,a><h,b> >= ||h||^4? ---")
term1_sufficient = sum(1 for d in data if 4*d['h_alpha']*d['h_beta'] >= d['h4'] - 1e-10)
print(f"  4<h,a><h,b> >= h^4: {term1_sufficient}/{len(data)}")

# Hmm, probably not. Let's check.
ratios_4hahb = [4*d['h_alpha']*d['h_beta']/d['h4'] for d in data if d['h4'] > 1e-10]
print(f"  min(4<h,a><h,b>/h^4) = {min(ratios_4hahb):.6f}")
print(f"  mean = {np.mean(ratios_4hahb):.6f}")

# Analysis 5: Try the 2x2 matrix approach
# [[||u||^2, <u,v>], [<u,v>, ||v||^2]] is PSD
# => ||u||^2*||v||^2 >= <u,v>^2
# This gives Phi_p*Phi_q >= <u,v>^2
# We need: (Phi_p-Phi_r)(Phi_q-Phi_r) >= Phi_r^2
# i.e., Phi_p*Phi_q - (Phi_p+Phi_q)*Phi_r + Phi_r^2 >= Phi_r^2
# i.e., Phi_p*Phi_q >= (Phi_p+Phi_q)*Phi_r
# i.e., 1/Phi_r >= (Phi_p+Phi_q)/(Phi_p*Phi_q) = 1/Phi_q + 1/Phi_p. YES!
# WAIT: This IS the target inequality!

# So the question reduces to: can we show Phi_p*Phi_q >= (Phi_p+Phi_q)*Phi_r?
# Equivalently: Phi_r <= Phi_p*Phi_q/(Phi_p+Phi_q) = harmonic mean of Phi_p, Phi_q.

# Let's explore the <u,v> angle approach.
print("\n--- <u,v> analysis ---")
uv_bound_tests = []
for d in data:
    uv = d['uv']
    norm_u2 = d['norm_h2'] + d['A']
    norm_v2 = d['norm_h2'] + d['B']
    Phi_r = d['norm_h2']

    # The target: norm_u2 * norm_v2 >= (norm_u2 + norm_v2) * Phi_r
    lhs = norm_u2 * norm_v2
    rhs = (norm_u2 + norm_v2) * Phi_r
    uv_bound_tests.append(lhs >= rhs - 1e-10)

print(f"  ||u||^2*||v||^2 >= (||u||^2+||v||^2)*||h||^2: {sum(uv_bound_tests)}/{len(uv_bound_tests)}")

# This should be equivalent to the original inequality. Let's verify.
# ||u||^2 * ||v||^2 >= (||u||^2 + ||v||^2)*||h||^2
# Let a = ||u||^2/||h||^2, b = ||v||^2/||h||^2 (both > 1)
# Then: a*b >= a + b, i.e., ab - a - b >= 0, i.e., (a-1)(b-1) >= 1.
# With a = Phi_p/Phi_r, b = Phi_q/Phi_r:
# (Phi_p/Phi_r - 1)(Phi_q/Phi_r - 1) >= 1
# (Phi_p - Phi_r)(Phi_q - Phi_r) >= Phi_r^2
# AB >= h^4. YES, this is exactly the target.

print("\nSo the inequality is equivalent to (a-1)(b-1) >= 1 where a = Phi_p/Phi_r, b = Phi_q/Phi_r.")
print("This is the HM-GM-AM inequality: Phi_r <= HM(Phi_p, Phi_q)")

# Analysis 6: Explore direct proof via <h,u>
# <h,u> = ||h||^2 + <h,alpha> >= ||h||^2 (since <h,alpha> >= 0)
# <h,v> = ||h||^2 + <h,beta> >= ||h||^2
# By Cauchy-Schwarz: <h,u>^2 <= ||h||^2*||u||^2 = Phi_r * Phi_p
# So (||h||^2 + <h,alpha>)^2 <= Phi_r * Phi_p
# Similarly (||h||^2 + <h,beta>)^2 <= Phi_r * Phi_q
# Multiply: (||h||^2 + <h,alpha>)^2 * (||h||^2 + <h,beta>)^2 <= Phi_r^2 * Phi_p * Phi_q
# Since each factor >= ||h||^4:
# ||h||^4 * ||h||^4 <= Phi_r^2 * Phi_p * Phi_q
# Phi_r^2 <= Phi_p * Phi_q (which we already knew)

# We need BETTER: use the fact that <h,alpha> and <h,beta> are NOT too small.

# Analysis 7: The key identity
# <h,alpha> = <h, u-h> = <h,u> - ||h||^2
# <h,beta> = <h, v-h> = <h,v> - ||h||^2
# AB = (||u||^2 - ||h||^2)(||v||^2 - ||h||^2)
# = (||h||^2 + 2<h,alpha> + ||alpha||^2 - ||h||^2)(similar for beta)
# = (2<h,alpha> + ||alpha||^2)(2<h,beta> + ||beta||^2)

# Alternative: write u = h + alpha, ||u||^2 = ||h||^2 + 2<h,alpha> + ||alpha||^2
# Then ||u||^2 * ||v||^2 - (||u||^2 + ||v||^2)*||h||^2 + ||h||^4
# = ||u||^2*||v||^2 - ||u||^2*||h||^2 - ||v||^2*||h||^2 + ||h||^4
# = (||u||^2 - ||h||^2)(||v||^2 - ||h||^2)
# We need this >= 0, which it is since both factors are positive (A,B > 0).
# Wait, we need it >= ||h||^4, not >= 0!

# The target is (a-1)(b-1) >= 1. Can we use the Cauchy-Schwarz and <h,alpha> >= 0
# in a more clever way?

# Approach: use the MATRIX
# M = [u | v] in R^{n x 2}
# M^T M = [[||u||^2, <u,v>], [<u,v>, ||v||^2]]
# det(M^T M) = ||u||^2*||v||^2 - <u,v>^2 >= 0 (always)
# But we need to relate <u,v> to ||h||^2.
# <u,v> = <h+alpha, h+beta> = ||h||^2 + <h,alpha> + <h,beta> + <alpha,beta>

# Key insight: if we could show <u,v> >= ||h||^2 + ||h||^2 = 2||h||^2... no, that's not right.

# Let's look at the ratio <u,v>/(||u||*||v||) = cos(angle between u and v)
print("\n--- cos(u,v) analysis ---")
cos_uv_vals = []
for d in data:
    norm_u = np.sqrt(d['norm_h2'] + d['A'])
    norm_v = np.sqrt(d['norm_h2'] + d['B'])
    cos_uv = d['uv'] / (norm_u * norm_v)
    cos_uv_vals.append(cos_uv)

print(f"  min cos(u,v) = {min(cos_uv_vals):.6f}")
print(f"  max cos(u,v) = {max(cos_uv_vals):.6f}")
print(f"  mean cos(u,v) = {np.mean(cos_uv_vals):.6f}")

# Analysis 8: Check if there's a SIMPLE algebraic identity
# For n=2 we showed AB = h^4 (equality). What is the excess for n=3?
print("\n--- AB/h^4 ratio ---")
for n_val in [3, 4, 5, 6]:
    ratios = [d['AB']/d['h4'] for d in data if d['n'] == n_val and d['h4'] > 1e-10]
    if ratios:
        print(f"  n={n_val}: min={min(ratios):.6f}, mean={np.mean(ratios):.4f}, max={max(ratios):.4f}")

# Analysis 9: Cross terms
# AB = (2<h,alpha> + ||alpha||^2)(2<h,beta> + ||beta||^2)
# Let x = <h,alpha>/||h||^2, y = <h,beta>/||h||^2
# a2 = ||alpha||^2/||h||^2, b2 = ||beta||^2/||h||^2
# Then AB/h^4 = (2x + a2)(2y + b2)
# Need (2x + a2)(2y + b2) >= 1
# By CS: x^2 <= a2, y^2 <= b2 (since <h,alpha>^2 <= ||h||^2 * ||alpha||^2)

# So x = cos(theta_a) * sqrt(a2), y = cos(theta_b) * sqrt(b2)
# (2 cos_a sqrt(a2) + a2)(2 cos_b sqrt(b2) + b2) >= 1

print("\n--- Normalized quantities (x, y, a2, b2) ---")
for d in data[:5]:
    x = d['h_alpha'] / d['norm_h2']
    y = d['h_beta'] / d['norm_h2']
    a2 = d['norm_alpha2'] / d['norm_h2']
    b2 = d['norm_beta2'] / d['norm_h2']
    ab_ratio = d['AB'] / d['h4']
    print(f"  n={d['n']}: x={x:.4f}, y={y:.4f}, a2={a2:.4f}, b2={b2:.4f}, "
          f"(2x+a2)={2*x+a2:.4f}, (2y+b2)={2*y+b2:.4f}, product={ab_ratio:.6f}")

# Analysis 10: The ORTHOGONAL DECOMPOSITION approach
# Write alpha = c_a * h/||h|| + alpha_perp where alpha_perp perp h
# Then <h,alpha> = c_a * ||h||, so c_a = <h,alpha>/||h||
# ||alpha||^2 = c_a^2 + ||alpha_perp||^2
# A = 2*c_a*||h|| + c_a^2 + ||alpha_perp||^2

# Let s_a = c_a/||h|| = <h,alpha>/||h||^2, t_a^2 = ||alpha_perp||^2/||h||^2
# Then A/||h||^2 = 2*s_a + s_a^2 + t_a^2 = (1+s_a)^2 - 1 + t_a^2
# Similarly B/||h||^2 = (1+s_b)^2 - 1 + t_b^2
# AB/||h||^4 = [(1+s_a)^2 - 1 + t_a^2] * [(1+s_b)^2 - 1 + t_b^2]

# For n=2: alpha is parallel to h (cos = 1), so t_a = t_b = 0.
# AB/h4 = [(1+s_a)^2 - 1] * [(1+s_b)^2 - 1] = s_a(2+s_a) * s_b(2+s_b)
# We showed this = 1 for n=2.

print("\n--- Orthogonal decomposition s_a, t_a, s_b, t_b ---")
for d in data[:10]:
    s_a = d['h_alpha'] / d['norm_h2']
    s_b = d['h_beta'] / d['norm_h2']
    t_a2 = d['norm_alpha2'] / d['norm_h2'] - (d['h_alpha'])**2 / d['norm_h2']**2
    t_b2 = d['norm_beta2'] / d['norm_h2'] - (d['h_beta'])**2 / d['norm_h2']**2
    fac_a = (1+s_a)**2 - 1 + t_a2
    fac_b = (1+s_b)**2 - 1 + t_b2
    print(f"  n={d['n']}: s_a={s_a:.4f}, t_a^2={t_a2:.6f}, s_b={s_b:.4f}, t_b^2={t_b2:.6f}, "
          f"fac_a={fac_a:.4f}, fac_b={fac_b:.4f}, product={fac_a*fac_b:.6f}")

# KEY OBSERVATION: With the orthogonal decomposition,
# fac_a = s_a(2+s_a) + t_a^2 >= s_a(2+s_a) (since t_a^2 >= 0)
# fac_b = s_b(2+s_b) + t_b^2 >= s_b(2+s_b)

# So AB/h^4 >= s_a(2+s_a) * s_b(2+s_b)  [using only the h-parallel component]

# For n=2, s_a(2+s_a)*s_b(2+s_b) = 1 exactly.
# For n >= 3, we need s_a(2+s_a)*s_b(2+s_b) >= 1 - epsilon,
# where the t^2 terms make up the difference.

# Actually, with t_a^2, t_b^2, the product is:
# = s_a(2+s_a)*s_b(2+s_b) + s_a(2+s_a)*t_b^2 + t_a^2*s_b(2+s_b) + t_a^2*t_b^2
# Need total >= 1.

# Check: is s_a(2+s_a)*s_b(2+s_b) >= 1?
print("\n--- Is the parallel part alone sufficient? ---")
for d in data[:5]:
    s_a = d['h_alpha'] / d['norm_h2']
    s_b = d['h_beta'] / d['norm_h2']
    parallel = s_a*(2+s_a) * s_b*(2+s_b)
    print(f"  n={d['n']}: parallel_product = {parallel:.6f} {'>=1 YES' if parallel >= 1-1e-8 else '< 1 NO'}")

parallel_sufficient = sum(1 for d in data
    if (d['h_alpha']/d['norm_h2'] * (2 + d['h_alpha']/d['norm_h2'])) *
       (d['h_beta']/d['norm_h2'] * (2 + d['h_beta']/d['norm_h2'])) >= 1 - 1e-8)
print(f"  Parallel part >= 1: {parallel_sufficient}/{len(data)}")


# Analysis 11: Is there a nice expression for s_a * s_b?
# s_a = <h,alpha>/||h||^2, s_b = <h,beta>/||h||^2
# s_a * s_b = <h,alpha>*<h,beta>/||h||^4
# Recall: alpha - beta = u - v. So <h,alpha-beta> = <h,u-v> = <h,u> - <h,v>
# = (||h||^2 + <h,alpha>) - (||h||^2 + <h,beta>) = <h,alpha> - <h,beta>.
# This is tautological.

# But: alpha + beta = (u - h) + (v - h) = u + v - 2h
# <h, alpha + beta> = <h,u> + <h,v> - 2||h||^2 = <h,alpha> + <h,beta>
# Also tautological.

# New approach: <u,v> = <h+alpha, h+beta> = ||h||^2 + <h,alpha> + <h,beta> + <alpha,beta>
# = ||h||^2(1 + s_a + s_b) + <alpha,beta>

print("\n--- <alpha,beta>/||h||^2 analysis ---")
ab_ratios = [d['alpha_beta']/d['norm_h2'] for d in data]
print(f"  min = {min(ab_ratios):.6f}")
print(f"  max = {max(ab_ratios):.6f}")
print(f"  mean = {np.mean(ab_ratios):.6f}")
print(f"  <alpha,beta> >= 0: {sum(1 for x in ab_ratios if x >= -1e-8)}/{len(ab_ratios)}")

# Analysis 12: Check if the KEY bound works via AM-GM
# We want (2s_a + a2)(2s_b + b2) >= 1
# where a2 = s_a^2 + t_a^2, b2 = s_b^2 + t_b^2
# = (s_a + (s_a + t_a^2/...)) * (s_b + ...)
# This is getting complex. Let me try a completely different approach.

# APPROACH via ||u-v||
# ||u-v||^2 = ||alpha-beta||^2
# ||u||^2 + ||v||^2 - 2<u,v> = ||alpha-beta||^2
# <u,v> = (||u||^2 + ||v||^2 - ||alpha-beta||^2) / 2
# ||u||^2*||v||^2 >= <u,v>^2 = (||u||^2 + ||v||^2 - ||alpha-beta||^2)^2 / 4

# The target: ||u||^2*||v||^2 - (||u||^2+||v||^2)*||h||^2 + ||h||^4 >= 0
# Let P = Phi_p, Q = Phi_q, R = Phi_r = ||h||^2
# Target: PQ - (P+Q)R + R^2 >= 0, i.e., (P-R)(Q-R) >= R^2... wait that's >= 0 not >= R^2.
# NO: The target is (P-R)(Q-R) >= R^2.

# Hmm, let me think about this differently.
# 1/R >= 1/P + 1/Q
# R <= PQ/(P+Q)  [R is at most the harmonic mean of P and Q]
# PQ/(P+Q) - R >= 0
# PQ - R(P+Q) >= 0
# This is weaker than (P-R)(Q-R) >= R^2.
# Actually: (P-R)(Q-R) = PQ - R(P+Q) + R^2. So (P-R)(Q-R) >= R^2 iff PQ >= R(P+Q).
# And PQ >= R(P+Q) iff PQ/(P+Q) >= R iff R <= HM(P,Q).
# YES, they are equivalent.

# So the question is: can we show R <= PQ/(P+Q)?
# This is the "subharmonic" property: Phi_r <= HM(Phi_p, Phi_q).

# Given <h,alpha> >= 0, <h,beta> >= 0, we know P >= R and Q >= R.
# We need the STRONGER: PQ >= R(P+Q), i.e., 1/R >= 1/P + 1/Q.

# A possible strategy exploiting <h,alpha> >= 0:
# Since u = h + alpha, we have ||u||^2 = ||h||^2 + 2<h,alpha> + ||alpha||^2 >= ||h||^2
# P = ||u||^2 >= R + 2<h,alpha> >= R (sharp bound?)
# Similarly Q >= R + 2<h,beta>

# The lower bound: P - R >= 2<h,alpha> and Q - R >= 2<h,beta>
# So (P-R)(Q-R) >= 4<h,alpha><h,beta>
# Need: 4<h,alpha><h,beta> >= R^2.

# Check this numerically:
print("\n--- Does 4<h,alpha><h,beta> >= R^2 = ||h||^4? ---")
sufficient_count = sum(1 for d in data if 4*d['h_alpha']*d['h_beta'] >= d['h4'] - 1e-8)
print(f"  YES: {sufficient_count}/{len(data)}")

# If NOT always, then we need the ||alpha||^2 and ||beta||^2 terms.
# Let's check what the gap looks like:
gaps = [(4*d['h_alpha']*d['h_beta'] - d['h4'])/d['h4'] for d in data if d['h4'] > 1e-10]
print(f"  min(4<h,a><h,b>/h^4 - 1) = {min(gaps):.6f}")
print(f"  mean = {np.mean(gaps):.6f}")

# CRITICAL: What if we use the FULL expansion?
# (P-R)(Q-R) = (2<h,a> + ||a||^2)(2<h,b> + ||b||^2)
#            >= (2<h,a>)(2<h,b>) + (2<h,a>)(||b||^2) + (||a||^2)(2<h,b>)
#            = 4<h,a><h,b> + 2(<h,a>||b||^2 + ||a||^2<h,b>)
# By AM-GM on the last two terms:
# <h,a>||b||^2 + ||a||^2<h,b> >= 2*sqrt(<h,a><h,b>*||a||^2*||b||^2)
# = 2||a||*||b||*sqrt(<h,a><h,b>)

# This is getting complicated. Let's try a DIFFERENT approach entirely.

# ==========================================================
# THE CAUCHY-SCHWARZ WITH WEIGHTS APPROACH (NODE 1.8)
# ==========================================================

print("\n\n" + "="*70)
print("NODE 1.8: CAUCHY-SCHWARZ WITH WEIGHTS APPROACH")
print("="*70)

# Key idea: instead of plain CS on <h,u>, use WEIGHTED CS.
# Consider the matrix M = [u | v] and h.
# We have: <h,u> = ||h||^2 + <h,alpha> >= ||h||^2
# And:     <h,v> = ||h||^2 + <h,beta> >= ||h||^2

# SCHUR COMPLEMENT: The 3x3 Gram matrix
# G = [[||u||^2, <u,v>, <u,h>],
#      [<u,v>,  ||v||^2, <v,h>],
#      [<u,h>,  <v,h>,  ||h||^2]]
# is PSD (Gram matrix of u, v, h).

# The Schur complement condition says:
# ||h||^2 - [<u,h>, <v,h>] * [[||u||^2, <u,v>], [<u,v>, ||v||^2]]^{-1} * [<u,h>; <v,h>] >= 0

# Let's compute this and see what it gives.
print("\n--- 3x3 Gram matrix Schur complement ---")
for d in data[:5]:
    h, alpha, beta = d['h'], d['alpha'], d['beta']
    u = h + alpha
    v = h + beta
    R = d['norm_h2']
    P = R + d['A']
    Q = R + d['B']
    hu = R + d['h_alpha']
    hv = R + d['h_beta']
    uv = d['uv']

    G = np.array([[P, uv, hu], [uv, Q, hv], [hu, hv, R]])
    eigvals = np.linalg.eigvalsh(G)
    det_G = np.linalg.det(G)

    # Schur complement of R in G:
    # det(G) / R = R*(P*Q - uv^2) - (hu^2*Q + hv^2*P - 2*hu*hv*uv)
    # ... actually det(G) = R*(PQ - uv^2) - hu*(hu*Q - hv*uv) + hv*(hu*uv - hv*P)
    # Hmm let me just compute.

    print(f"  n={d['n']}: det(G)={det_G:.6f}, eigenvalues={np.sort(eigvals)}")
    print(f"    PQ-uv^2 = {P*Q - uv**2:.6f}")
    print(f"    (P-R)(Q-R)-R^2 = {d['AB_h4']:.6f}")

# The Gram matrix approach gives: det(G) >= 0, which is:
# R*(PQ - uv^2) - hu*(hu*Q - hv*uv) + hv*(hu*uv - hv*P) >= 0 ... messy.
# Better: det(G) = R*PQ + 2*hu*hv*uv - R*uv^2 - P*hv^2 - Q*hu^2.
# This is the mixed discriminant condition. It's hard to use directly.

# Let me try yet another approach.

# APPROACH: Use u = h + alpha, v = h + beta, and the PARALLELOGRAM identity.
# ||u||^2 + ||v||^2 = 2||h||^2 + 2<h,alpha+beta> + ||alpha||^2 + ||beta||^2
# ||u-v||^2 = ||alpha-beta||^2
# ||u+v||^2 = ||2h + alpha + beta||^2

# Hmm. Let me think about what we KNOW and what we NEED.

# KNOW:
# (K1) u = h + alpha, v = h + beta (vector decomposition)
# (K2) <h,alpha> >= 0, <h,beta> >= 0 (numerically verified, needs proof)
# (K3) ||u||^2 = Phi_p, ||v||^2 = Phi_q, ||h||^2 = Phi_r (definitions)
# (K4) alpha_k = omega_1''(nu_k)/2, beta_k = omega_2''(nu_k)/2 (subordination)
# (K5) omega_1 + omega_2 = z + 1/(nG_r(z))  (subordination identity)

# NEED: (Phi_p - Phi_r)(Phi_q - Phi_r) >= Phi_r^2

# From K2: Phi_p >= Phi_r, Phi_q >= Phi_r.

# From the subordination identity K5:
# omega_1(z) + omega_2(z) = z + F_r^{-1}(nG_r(z))... hmm, this is the free prob version.
# In the finite case: omega_1(z) + omega_2(z) - z = (some function of z).

# At z = nu_k: omega_1(nu_k) + omega_2(nu_k) = lambda_{sigma(k)} + mu_{tau(k)}
# omega_1'(nu_k) + omega_2'(nu_k) = 1 + 1 = 2
# omega_1''(nu_k) + omega_2''(nu_k) = alpha_k*2 + beta_k*2 => alpha_k + beta_k = ...

# KEY: Differentiating omega_1 + omega_2 = z + something:
# If omega_1 + omega_2 = z + phi(z), then
# omega_1' + omega_2' = 1 + phi'(z)
# At z = nu_k: 2 = 1 + phi'(nu_k), so phi'(nu_k) = 1.
# omega_1'' + omega_2'' = phi''(z)
# So alpha + beta = phi''/2 where phi = omega_1 + omega_2 - z.

# phi(z) = omega_1(z) + omega_2(z) - z.
# phi maps reals to reals (since omega_1, omega_2 do).
# phi'(nu_k) = 1 for all k.
# What is phi more explicitly?

# In the infinite free probability case:
# F_r(z) = z - 1/(nG_r(z)), and omega_i are the subordination functions satisfying
# F_r(z) = F_p(omega_1(z)) = F_q(omega_2(z)) and omega_1(z) + omega_2(z) = z + F_r(z)... hmm.

# Actually in the free prob setup (Biane):
# The R-transform satisfies: R_{mu boxplus nu}(z) = R_mu(z) + R_nu(z)
# And the subordination functions satisfy: omega_1(z) + omega_2(z) = z + R_{mu boxplus nu}(G_r(z))

# For the finite case, the relation might be different. Let's just check numerically.

print("\n--- omega_1(nu_k) + omega_2(nu_k) structure ---")
for d in data[:5]:
    h, alpha, beta = d['h'], d['alpha'], d['beta']
    sum_alpha_beta = alpha + beta
    print(f"  n={d['n']}: alpha+beta = {sum_alpha_beta}, ||alpha+beta||^2 = {np.dot(sum_alpha_beta, sum_alpha_beta):.6f}")
    print(f"    <h, alpha+beta> = {np.dot(h, sum_alpha_beta):.6f} = <h,a>+<h,b> = {d['h_alpha']+d['h_beta']:.6f}")

# FINAL KEY OBSERVATION: Check if there's a clean sufficient condition
print("\n\n--- TESTING SUFFICIENT CONDITIONS ---")

# Test: does (2<h,alpha> + ||alpha||^2) * (2<h,beta> + ||beta||^2) >= ||h||^4 follow from
# some algebraic identity involving alpha, beta, h?

# Given alpha_k = omega_1''(nu_k)/2 and the subordination structure,
# what constraints exist between alpha and h?

# From the chain rule at each point:
# H_r(nu_k) = H_p(omega_1(nu_k)) - omega_1''(nu_k)/2
# h_k = u_k - alpha_k
# So alpha_k = u_k - h_k and beta_k = v_k - h_k.

# The vectors u and v live in a specific CONE determined by the interlacing.
# u_k = H_p(lambda_k) = sum_{j!=k} 1/(lambda_k - lambda_j)
# v_k = H_q(mu_k) = sum_{j!=k} 1/(mu_k - mu_j)
# h_k = H_r(nu_k) = sum_{j!=k} 1/(nu_k - nu_j)

# The key structural fact is that ALL three sets of H-values have the same SIGN PATTERN:
# h_1 < 0, h_2 > 0 alternating (no, that's not quite right for general polynomials)

# Let's check the sign patterns
print("\n--- Sign patterns of h, u, v ---")
for d in data[:5]:
    h, alpha, beta = d['h'], d['alpha'], d['beta']
    u = h + alpha
    v = h + beta
    print(f"  n={d['n']}:")
    print(f"    h signs: {['+'if x>0 else '-' for x in h]}")
    print(f"    u signs: {['+'if x>0 else '-' for x in u]}")
    print(f"    v signs: {['+'if x>0 else '-' for x in v]}")
    print(f"    alpha signs: {['+'if x>0 else '-' for x in alpha]}")
