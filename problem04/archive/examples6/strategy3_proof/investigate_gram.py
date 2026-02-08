"""
Investigate the near-zero Gram determinant at n=3 and the parallel structure.

Key observation from deep analysis:
- For n=3, det(G) of the Gram matrix [u,v,h] is near zero!
- This means u, v, h are nearly coplanar in R^3 (they're IN R^3 for n=3).
- For n=3, u,v,h are vectors in R^3, so they span at most a 3D space.
  The Gram matrix is 3x3, and det = 0 means they're linearly dependent
  (which generically they are NOT, since generically 3 vectors in R^3 are a basis).
  But det(G) ~ 0 numerically suggests NEAR linear dependence.

Actually wait: for n=3, h, u, v are in R^3. The Gram matrix of these 3 vectors in R^3
has determinant = |det([h|u|v])|^2 (up to sign). So det(G) = 0 iff h, u, v are
linearly dependent in R^3.

Since u = h + alpha, v = h + beta, we have {h, u, v} linearly dependent iff
{h, alpha, beta} linearly dependent. Since alpha - beta = u - v, we have
beta = alpha - (u-v), so alpha and beta always lie in the plane spanned by alpha and u-v.
{h, alpha, beta} dependent iff h lies in span(alpha, beta) iff h lies in span(alpha, u-v).

For n=3 (3-dimensional vectors), this is generically true only for a 2-dimensional
subspace. But the subordination structure might force this!

Let's check: for n=3, is h always in span(alpha, beta)?
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

# ================================================================
# Test: for n=2, h, alpha, beta are in R^2.
# u = h + alpha, v = h + beta, both in R^2.
# The Gram matrix of [h, u, v] in R^2 is always singular (3 vectors in R^2).
# For n=2, h is PARALLEL to alpha and beta (cos = 1 always).
# Generically, for n=3, h, alpha, beta span R^3.
# ================================================================

print("="*70)
print("GRAM DETERMINANT AND DIMENSIONALITY")
print("="*70)

for n_val in [2, 3, 4, 5]:
    print(f"\nn={n_val}:")
    det_vals = []
    coplanar_count = 0

    for trial in range(200):
        roots_p = np.sort(np.random.randn(n_val) * 2)
        for i in range(1, n_val):
            if roots_p[i] - roots_p[i-1] < 0.5:
                roots_p[i] = roots_p[i-1] + 0.5
        roots_q = np.sort(np.random.randn(n_val) * 2)
        for i in range(1, n_val):
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

            _, _, omega1_pp, sigma = subordination(roots_p, roots_r)
            _, _, omega2_pp, tau = subordination(roots_q, roots_r)

            h = H_values(roots_r)
            alpha = omega1_pp / 2
            beta = omega2_pp / 2

            # Check linear dependence of {h, alpha, beta}
            M = np.column_stack([h, alpha, beta])
            if n_val >= 3:
                # SVD to check rank
                s = np.linalg.svd(M, compute_uv=False)
                min_sv = s[-1]
                if min_sv < 1e-8:
                    coplanar_count += 1

            # Gram matrix
            u = h + alpha
            v = h + beta
            G = np.array([
                [np.dot(h,h), np.dot(h,u), np.dot(h,v)],
                [np.dot(u,h), np.dot(u,u), np.dot(u,v)],
                [np.dot(v,h), np.dot(v,u), np.dot(v,v)],
            ])
            det_G = np.linalg.det(G)
            det_vals.append(det_G)

        except:
            pass

    if det_vals:
        print(f"  det(G) stats: min={min(det_vals):.6e}, max={max(det_vals):.6e}, mean={np.mean(det_vals):.6e}")
        print(f"  |det(G)| < 1e-8: {sum(1 for d in det_vals if abs(d) < 1e-8)}/{len(det_vals)}")
        if n_val >= 3:
            print(f"  h,alpha,beta coplanar (min_sv < 1e-8): {coplanar_count}/{len(det_vals)}")


# ================================================================
# KEY IDENTITY TEST: alpha_k + beta_k = f(nu_k, roots_r)
# ================================================================
print("\n\n" + "="*70)
print("TESTING: alpha + beta STRUCTURE")
print("="*70)

# For each k: alpha_k + beta_k = omega_1''(nu_k)/2 + omega_2''(nu_k)/2
# And h_k = H_r(nu_k) = sum_{j!=k} 1/(nu_k - nu_j)
# Also: u_k = H_p(lambda_k), v_k = H_q(mu_k)
# alpha_k = u_k - h_k, beta_k = v_k - h_k
# alpha_k + beta_k = u_k + v_k - 2h_k

# Is there a simple formula for alpha_k + beta_k in terms of roots of r only?
# Check if alpha_k + beta_k = c * h_k for some constant c... probably not.
# Check if alpha_k + beta_k = f(nu_k) for some function f.

print("\nChecking if (alpha+beta)/h is constant componentwise:")
for trial in range(10):
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
        if np.any(np.abs(np.imag(np.roots(np.poly(roots_r)))) > 0.01):
            continue
        roots_r = np.sort(np.real(roots_r))
        if np.any(np.diff(roots_r) < 0.01):
            continue

        _, _, omega1_pp, sigma = subordination(roots_p, roots_r)
        _, _, omega2_pp, tau = subordination(roots_q, roots_r)

        h = H_values(roots_r)
        alpha = omega1_pp / 2
        beta = omega2_pp / 2

        if np.min(np.abs(h)) > 0.01:
            ratio = (alpha + beta) / h
            print(f"  n={n}: (alpha+beta)/h = {ratio}")
    except:
        pass


# ================================================================
# CRITICAL TEST: <h, alpha> >= 0 PROOF ATTEMPT
# ================================================================
print("\n\n" + "="*70)
print("PROVING <h,alpha> >= 0")
print("="*70)

# alpha_k = omega_1''(nu_k)/2
# h_k = H_r(nu_k) = sum_{j!=k} 1/(nu_k - nu_j)

# From the implicit function F(z,w) = r'(z)p(w) - p'(w)r(z) = 0:
# At (nu_k, lambda_k): r(nu_k) = p(lambda_k) = 0
# omega_1''(nu_k) = -[F_zz + 2F_zw + F_ww] / F_w  (with omega_1' = 1)
# = -[-p'(lambda_k)*r''(nu_k) + 2(r''(nu_k)*p'(lambda_k) - p''(lambda_k)*r'(nu_k)) + r'(nu_k)*p''(lambda_k)] / [r'(nu_k)*p'(lambda_k)]
# = [p'(lambda_k)*r''(nu_k) - 2r''(nu_k)*p'(lambda_k) + 2p''(lambda_k)*r'(nu_k) - r'(nu_k)*p''(lambda_k)] / [r'(nu_k)*p'(lambda_k)]
# = [-r''(nu_k)*p'(lambda_k) + p''(lambda_k)*r'(nu_k)] / [r'(nu_k)*p'(lambda_k)]
# = p''(lambda_k)/p'(lambda_k) - r''(nu_k)/r'(nu_k)
# = 2*H_p(lambda_k) - 2*H_r(nu_k)
# = 2*(u_k - h_k) = 2*alpha_k

# WAIT! This is CIRCULAR: omega_1''(nu_k)/2 = u_k - h_k = alpha_k.
# That's just the definition. Not useful for proving <h,alpha> >= 0.

# Let me think about this differently.
# alpha_k = H_p(lambda_k) - H_r(nu_k) = sum_{j!=k} [1/(lambda_k - lambda_j) - 1/(nu_k - nu_j)]

# <h, alpha> = sum_k H_r(nu_k) * [H_p(lambda_k) - H_r(nu_k)]
#            = sum_k h_k * u_k - sum_k h_k^2
#            = <h, u> - ||h||^2

# So <h,alpha> >= 0 iff <h,u> >= ||h||^2.

# Now u = H_p reindexed. Can we show <h,u> >= ||h||^2?

# <h,u> = sum_k H_r(nu_k) * H_p(lambda_k)
# ||h||^2 = sum_k H_r(nu_k)^2

# So the question is: sum_k H_r(nu_k) * H_p(lambda_k) >= sum_k H_r(nu_k)^2
# i.e., sum_k H_r(nu_k) * [H_p(lambda_k) - H_r(nu_k)] >= 0

# This is NOT obvious. Let's look at it term by term.

# H_r(nu_k) = sum_{j!=k} 1/(nu_k - nu_j)
# H_p(lambda_k) = sum_{j!=k} 1/(lambda_k - lambda_j)

# The roots of r interlace with the roots of p in some way (via the convolution).
# The key structural property must come from the specific relationship between
# the roots of r and the roots of p.

# For the MSS convolution: lambda_k and nu_k are related by:
# nu_k is the k-th root of p boxplus_n q
# lambda_k is the k-th root of p
# omega_1(nu_k) = lambda_k (subordination)

# The subordination function omega_1 is real, increasing on each interval between
# the poles. omega_1'(nu_k) = 1 > 0, so omega_1 is locally increasing.

# KEY PROPERTY: omega_1 is a Herglotz function restricted to the real line.
# Between consecutive roots of r, omega_1 is a real-analytic function.
# omega_1(nu_k) = lambda_k, omega_1'(nu_k) = 1.
# On the REAL LINE, omega_1 maps the interval [nu_k, nu_{k+1}] to [lambda_k, lambda_{k+1}].
# Since omega_1' = 1 at the endpoints and omega_1 maps interval to interval,
# the intermediate value theorem says omega_1 maps [nu_k, nu_{k+1}] onto [lambda_k, lambda_{k+1}]
# (assuming the Herglotz property makes omega_1 monotone).

# This means: lambda_{k+1} - lambda_k >= nu_{k+1} - nu_k for all k?
# NOT necessarily. omega_1 is monotone but not necessarily with slope >= 1 everywhere.
# It has slope 1 at nu_k, but could have slope < 1 in between.
# Actually, omega_1(nu_k) = lambda_k, omega_1(nu_{k+1}) = lambda_{k+1}.
# Mean value theorem: there exists xi in (nu_k, nu_{k+1}) with
# omega_1'(xi) = (lambda_{k+1} - lambda_k)/(nu_{k+1} - nu_k).
# This could be > 1 or < 1.

# The gap ratio:
print("\nGap ratios lambda_{k+1}-lambda_k vs nu_{k+1}-nu_k:")
for trial in range(5):
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
        gaps_p = np.diff(roots_p)
        gaps_r = np.diff(roots_r)
        ratios = gaps_p / gaps_r
        print(f"  gaps_p/gaps_r = {ratios}")
    except:
        pass


# ================================================================
# ATTEMPT: Use convexity/concavity of omega_1
# ================================================================
print("\n\n" + "="*70)
print("OMEGA_1 SECOND DERIVATIVE AND CONVEXITY")
print("="*70)

# omega_1''(nu_k)/2 = alpha_k = u_k - h_k
# <h,alpha> = sum_k h_k * alpha_k
# = sum_k H_r(nu_k) * omega_1''(nu_k)/2

# omega_1 is a Herglotz function on C+. On the real line between roots of r,
# omega_1 is well-defined and monotone increasing.

# For a Herglotz function, there's a representation:
# omega_1(z) = z + c + integral [1/(t-z) - t/(t^2+1)] d_mu(t)
# But our omega_1 satisfies omega_1'(nu_k) = 1, which suggests the correction
# term is "second order" at each root of r.

# KEY INSIGHT: Think about sum_k h_k * omega_1''(nu_k)
# = sum_k H_r(nu_k) * omega_1''(nu_k)
# = sum_k [sum_{j!=k} 1/(nu_k - nu_j)] * omega_1''(nu_k)

# This looks like a discrete version of int H_r(x) * omega_1''(x) dx.
# We might be able to use integration by parts (discrete or continuous).

# Another approach: <h, alpha> = sum_k h_k * (u_k - h_k)
# = sum_k h_k * u_k - sum_k h_k^2

# Can we use the CAUCHY-INTERLACING theorem or MAJORIZATION?
# The roots of p and r may satisfy some majorization relation.

# For the free additive convolution in the INFINITE case (mu boxplus nu):
# The support of mu boxplus nu is [a+c, b+d] where [a,b] = supp(mu), [c,d] = supp(nu).
# So the roots of r spread MORE than the roots of p alone.
# This suggests Phi_r < Phi_p (since wider spacing -> smaller H values -> smaller Phi).

# For the finite case: the roots of r = p boxplus_n q are MORE SPREAD than those of p or q.
# This is consistent with A = Phi_p - Phi_r > 0 (roots of r more spread => smaller Fisher info).

# But we need QUANTITATIVE control on HOW MUCH more spread.

# ================================================================
# FINAL APPROACH: Try to prove the inequality using matrix methods
# ================================================================

# Consider the n x n diagonal matrices:
# D_r = diag(nu_1,...,nu_n)
# D_p = diag(lambda_1,...,lambda_n)
# D_q = diag(mu_1,...,mu_n)

# H_r(nu_k) = tr((D_r - nu_k I)^{-1}) with the k-th diagonal removed
# = sum_{j!=k} 1/(nu_j - nu_k) ... wait, sign.

# H_r(nu_k) = sum_{j!=k} 1/(nu_k - nu_j) = d/dx [log |prod_{j!=k}(x - nu_j)|] at x = nu_k

# <h,u> = sum_k H_r(nu_k) * H_p(lambda_k)

# This is a sum involving TWO sets of roots. Hard to bound directly.

# Let's try NUMERICAL optimization to find the minimum of AB/h^4.
# From the data: min(AB/h^4) = 1.007 at n=3. Can we find the exact minimum?

print("\n\n" + "="*70)
print("FINDING MINIMUM OF AB/h^4 FOR n=3")
print("="*70)

from scipy.optimize import minimize

def neg_AB_over_h4(params):
    """Minimize -AB/h^4 for n=3."""
    n = 3
    # params = [a1, a2, a3, b1, b2, b3] = roots of p and q (sorted)
    a = np.sort(params[:n])
    b = np.sort(params[n:])

    # Ensure distinct
    for i in range(1, n):
        if a[i] - a[i-1] < 0.1:
            return 1e10
        if b[i] - b[i-1] < 0.1:
            return 1e10

    try:
        roots_r, _ = boxplus_mss(a, b)
        raw = np.roots(np.poly(roots_r))
        if np.any(np.abs(np.imag(raw)) > 0.01):
            return 1e10
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01):
            return 1e10

        Ph_p = sum(H_values(a)**2)
        Ph_q = sum(H_values(b)**2)
        Ph_r = sum(H_values(roots_r)**2)

        A = Ph_p - Ph_r
        B = Ph_q - Ph_r
        h4 = Ph_r**2

        if A <= 0 or B <= 0 or h4 <= 0:
            return 1e10

        return -A*B / h4  # minimize -AB/h4 to find minimum of AB/h4
    except:
        return 1e10

# Try many random starting points
best_val = 1e10
best_params = None
for _ in range(200):
    x0 = np.sort(np.random.randn(3) * 2).tolist() + np.sort(np.random.randn(3) * 2).tolist()
    for i in [1, 2, 4, 5]:
        prev = i - 1 if i != 3 else 2
        if i < 3 and x0[i] - x0[i-1] < 0.3:
            x0[i] = x0[i-1] + 0.3
        if i >= 3 and i > 3 and x0[i] - x0[i-1] < 0.3:
            x0[i] = x0[i-1] + 0.3

    try:
        res = minimize(neg_AB_over_h4, x0, method='Nelder-Mead',
                      options={'maxiter': 5000, 'xatol': 1e-10, 'fatol': 1e-12})
        if res.fun < best_val:
            best_val = res.fun
            best_params = res.x
    except:
        pass

print(f"Minimum AB/h^4 found: {-best_val:.10f}")
if best_params is not None:
    p = np.sort(best_params[:3])
    q = np.sort(best_params[3:])
    print(f"At p roots: {p}")
    print(f"At q roots: {q}")

    roots_r, _ = boxplus_mss(p, q)
    roots_r = np.sort(np.real(roots_r))
    print(f"r roots: {roots_r}")
    print(f"Phi_p={sum(H_values(p)**2):.8f}, Phi_q={sum(H_values(q)**2):.8f}, Phi_r={sum(H_values(roots_r)**2):.8f}")

    # What does the minimum look like? Is it approached when p and q are "similar"?
    print(f"Gap ratio p/q: {np.diff(p)/np.diff(q)}")


# Check: does the minimum approach 1 as the polynomials approach each other?
print("\n--- AB/h^4 as p -> q ---")
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
            print(f"  eps={epsilon:.3f}: AB/h4 = {A*B/h4:.8f}, Phi_p={Ph_p:.6f}, Phi_q={Ph_q:.6f}, Phi_r={Ph_r:.6f}")
    except:
        pass
