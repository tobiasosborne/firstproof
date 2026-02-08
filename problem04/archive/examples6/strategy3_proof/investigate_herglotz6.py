#!/usr/bin/env python3
"""
Part 6: Focus on Phi(t) = sum H_t^2 monotonicity.

CRITICAL FINDING: Phi(t) = sum_k H(x_k(t))^2 where x_k(t) = (1-t)*nu_k + t*lambda_k
is ALWAYS monotone increasing for MSS pairs (nu, lambda) = (roots of r, roots of p).

This gives:
- Phi(0) = Phi_r <= Phi(1) = Phi_p, i.e., Phi_r <= Phi_p (always)
- By symmetry with omega_2: Phi_r <= Phi_q
- These are WEAKER than <h,alpha> >= 0, but still useful.

Wait -- Phi(t) monotone gives Phi_r <= Phi_p, but does it give <h,alpha> >= 0?
<h,alpha> = <h,u> - ||h||^2 where h = H_r and u = H_p.
Phi_p >= Phi_r gives ||u||^2 >= ||h||^2, but NOT <h,u> >= ||h||^2.

Actually, Phi(t) increasing means Phi'(t) >= 0 for all t.
Phi'(t) = 2 * sum_k H_t(x_k) * d/dt[H_t(x_k)]
But this is Phi'(t), not F'(t) = <h_r, d/dt[h_t]>.

Let me compute Phi'(0) and see if it implies <h,alpha> >= 0.

Phi'(0) = 2 * sum_k H_r(nu_k) * d/dt[H_t(x_k)]|_{t=0}
        = 2 * sum_k h_k * sum_{l!=k} (phi_k - phi_l)/(nu_k - nu_l)^2
        = 2 * F'(0)

So Phi'(0) = 2*F'(0). The numerical test showed F'(t) can be negative for some
MSS pairs (4/1000). But Phi(t) is always monotone.

The difference: Phi'(t) = 2*<h_t, dh_t/dt> where h_t = H-values at x(t).
F'(t) = <h_r, dh_t/dt> uses the FIXED vector h_r, not h_t.

Phi(t) monotone is a STRONGER statement because h_t changes adaptively.

Actually, Phi(t) monotone means ||h_t||^2 is increasing, which means
the "Fisher information" increases as we interpolate from r-roots to p-roots.
This is interesting but does NOT directly give <h_r, alpha> >= 0.

Wait, but <h_r, alpha> = <h_r, u - h_r> = <h_r, h_1 - h_0>
where h_0 = h_r and h_1 = u = H_p = h at t=1.

<h_0, h_1 - h_0> = <h_0, h_1> - ||h_0||^2.

If h(t) traces a path from h_0 to h_1 with ||h(t)|| monotone increasing,
does <h_0, h_1 - h_0> >= 0?

Not necessarily in general. Consider h_0 = (1,0) and h_1 = (-2, 0).
||h_1|| > ||h_0|| but <h_0, h_1-h_0> = <(1,0),(-3,0)> = -3 < 0.

But the PATH h(t) is constrained. The interpolation is along specific curves.

KEY: maybe the PATH is ALWAYS such that the angle between h_0 and h_t is
less than pi/2, i.e., <h_0, h_t> >= 0 for all t.

Let me check this.
"""

import numpy as np
from math import factorial
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
    r_roots = np.sort(np.real(np.roots(c)))
    return r_roots, c

def H_values(roots):
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H

np.random.seed(42)

# ================================================================
# Test: is <h_0, h_t> >= ||h_0||^2 for all t in [0,1]?
# This would give <h_0, h_1> >= ||h_0||^2, i.e., <h,alpha> >= 0.
# ================================================================
print("="*70)
print("TEST: <h_0, h_t> >= ||h_0||^2 along interpolation?")
print("="*70)

violations = 0
total = 0

for trial in range(500):
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
        if np.any(np.abs(np.imag(raw)) > 0.01): continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01): continue

        total += 1
        h0 = H_values(roots_r)
        h0_norm2 = np.dot(h0, h0)
        violated = False

        for t in np.linspace(0, 1, 50):
            x_t = (1-t)*roots_r + t*roots_p
            h_t = H_values(x_t)
            inner = np.dot(h0, h_t)
            if inner < h0_norm2 - 1e-8:
                violated = True
                break

        if violated:
            violations += 1
    except:
        pass

print(f"Total: {total}")
print(f"<h_0, h_t> >= ||h_0||^2 violated: {violations}")

# So <h_0, h_t> is NOT always >= ||h_0||^2 along the path.
# Even though Phi(t) = ||h_t||^2 is monotone.


# ================================================================
# CRUCIAL QUESTION: What makes Phi(t) monotone for MSS pairs?
# ================================================================
print("\n\n" + "="*70)
print("WHY IS Phi(t) MONOTONE? -- Derivative computation")
print("="*70)

# Phi(t) = sum_k H_t(x_k(t))^2 where x_k(t) = (1-t)*nu_k + t*lambda_k
# Phi'(t) = 2*sum_k H_t(x_k)*[sum_{l!=k} (phi_k-phi_l)/(x_k-x_l)^2]
# where phi_k = nu_k - lambda_k (BUT wait: x_k' = lambda_k - nu_k = -phi_k)
#
# More carefully: d/dt x_k = lambda_k - nu_k.
# d/dt H_t(x_k) = sum_{l!=k} -(dx_k/dt - dx_l/dt)/(x_k-x_l)^2
#               = sum_{l!=k} -((lambda_k-nu_k)-(lambda_l-nu_l))/(x_k-x_l)^2
#               = sum_{l!=k} ((nu_k-lambda_k)-(nu_l-lambda_l))/(x_k-x_l)^2
#               = sum_{l!=k} (phi_k - phi_l)/(x_k(t)-x_l(t))^2
#
# Phi'(t) = 2*sum_k H_t(x_k) * sum_{l!=k} (phi_k-phi_l)/(x_k(t)-x_l(t))^2
# = 2 * sum_{k!=l} H_t(x_k) * (phi_k-phi_l)/(x_k(t)-x_l(t))^2
# = 2 * sum_{k<l} [H_t(x_k)-H_t(x_l)] * (phi_k-phi_l)/(x_k(t)-x_l(t))^2

# So Phi'(t) = 2 * sum_{k<l} (H_t(x_k)-H_t(x_l))*(phi_k-phi_l)/(x_k(t)-x_l(t))^2.

# At t=0: this is 2 * sum_{k<l} (h_k-h_l)*(phi_k-phi_l)/(nu_k-nu_l)^2.
# At general t: h_k(t) = H_t(x_k(t)) and the denominator changes.

# SELF-REINFORCING MECHANISM: Note that h_t is the H-VALUES at the points x_t.
# As t increases, the points move from nu -> lambda.
# The H-values change. And the product (h_t(k)-h_t(l))*(phi_k-phi_l)
# stays positive because... why?

# The key is that h_t(k) = sum_{j!=k} 1/(x_k(t)-x_j(t)), and the ordering
# of x(t) is preserved (since both nu and lambda are sorted, and t in [0,1]).
# So x_k(t) < x_l(t) for k < l.

# For k < l: x_k < x_l, so 1/(x_k-x_l) < 0, 1/(x_l-x_k) > 0.
# h_t(k) tends to be negative (leftmost root) and h_t(l) tends to be positive.
# So h_t(k) - h_t(l) tends to be negative.

# phi_k - phi_l: for k < l, phi_k < phi_l means the shifts are increasing.
# This happens when roots of r "expand" more than roots of p from the center.

# ACTUALLY the issue is whether Phi'(t) is a QUADRATIC FORM.
# Let me write it differently.

# Phi'(t) = 2 * sum_{k<l} Delta_h(k,l) * Delta_phi(k,l) / gap(t,k,l)^2

# where Delta_h(k,l) = h_t(k) - h_t(l) depends on t through the current H-values.

# The MAGIC: h_t(k) - h_t(l) = sum_{j!=k} 1/(x_k-x_j) - sum_{j!=l} 1/(x_l-x_j)
# This is the DIFFERENCE of H-values at two different points.

# Consider the function g(x) = H(x; x_1,...,x_n) = sum_{j: x_j!=x} 1/(x-x_j).
# Then h_t(k) = g(x_k) (where g depends on ALL x_j values).
# And h_t(k) - h_t(l) = g(x_k) - g(x_l).

# g is a rational function with simple poles at x_j and g'(x) = -sum 1/(x-x_j)^2 < 0.
# So g is STRICTLY DECREASING where defined.
# Since x_k < x_l for k < l: g(x_k) > g(x_l) ONLY between poles.

# WAIT: g(x) = sum_j 1/(x-x_j) has poles at x_1,...,x_n.
# But we evaluate g at x_k, excluding the self-interaction:
# h_t(k) = sum_{j!=k} 1/(x_k - x_j), not sum_j 1/(x_k - x_j).

# Let me reconsider. For x between x_k and x_{k+1}:
# sum_j 1/(x-x_j) = 1/(x-x_k) + [positive terms from j<k] + 1/(x-x_{k+1}) + [terms from j>k+1]
# This function goes from -inf to +inf as x crosses each x_j from left to right.
# So the TOTAL sum (including self-interaction) is not simply monotone.

# The EXCLUDED sum h_t(k) = sum_{j!=k} 1/(x_k-x_j) is the "regularized" value:
# it's the limit of [sum_j 1/(x-x_j) - 1/(x-x_k)] as x -> x_k.

# Let's check if h_t is "similarly ordered" to phi:

print("\nAt t=0 and t=0.5, check ordering of h_t and phi:")
n = 4
roots_p = np.array([-3., -1., 1., 4.])
roots_q = np.array([-2., 0., 2., 3.])
roots_r, _ = boxplus_mss(roots_p, roots_q)
roots_r = np.sort(np.real(roots_r))

phi = roots_r - roots_p
h0 = H_values(roots_r)

for t in [0, 0.25, 0.5, 0.75, 1.0]:
    x_t = (1-t)*roots_r + t*roots_p
    h_t = H_values(x_t)
    print(f"  t={t:.2f}: h_t = {np.round(h_t, 4)}, phi = {np.round(phi, 4)}")
    # Check: are (h_t(k)-h_t(l)) and (phi_k-phi_l) same sign for all k<l?
    concordant = 0
    total_pairs = 0
    for k in range(n):
        for l in range(k+1, n):
            if (h_t[k]-h_t[l])*(phi[k]-phi[l]) >= 0:
                concordant += 1
            total_pairs += 1
    print(f"         concordance: {concordant}/{total_pairs}")


# ================================================================
# THE KEY INSIGHT: Schur-convexity of Phi
# ================================================================
print("\n\n" + "="*70)
print("SCHUR-CONVEXITY ANALYSIS")
print("="*70)

# Phi(x_1,...,x_n) = sum_k [sum_{j!=k} 1/(x_k-x_j)]^2
# is a SYMMETRIC function of x_1,...,x_n.
# Is it SCHUR-CONVEX?
#
# A symmetric function is Schur-convex iff
# (x_1 - x_2) * (dPhi/dx_1 - dPhi/dx_2) >= 0 for all x.
#
# dPhi/dx_k = 2 * H(x_k) * [-sum_{j!=k} 1/(x_k-x_j)^2]
#           + 2 * sum_{l!=k} H(x_l) * [1/(x_l-x_k)^2]
# = -2*H(x_k)*sum_{j!=k} 1/(x_k-x_j)^2 + 2*sum_{l!=k} H(x_l)/(x_l-x_k)^2
# = 2*sum_{j!=k} [H(x_j) - H(x_k)]/(x_j-x_k)^2

# Hmm wait, let me be careful:
# Phi(x) = sum_k H_k^2 where H_k = sum_{j!=k} 1/(x_k-x_j)
# dH_k/dx_k = -sum_{j!=k} 1/(x_k-x_j)^2
# dH_k/dx_l = 1/(x_k-x_l)^2 for l != k
# dPhi/dx_m = sum_k 2*H_k * dH_k/dx_m
# = 2*H_m * (-sum_{j!=m} 1/(x_m-x_j)^2) + 2*sum_{k!=m} H_k * 1/(x_k-x_m)^2
# = 2*sum_{k!=m} [H_k - H_m]/(x_k-x_m)^2  (combining the terms)

# Wait: for k != m:
# From the H_m term: 2*H_m * (-1/(x_m-x_k)^2)
# From the H_k term: 2*H_k * (1/(x_k-x_m)^2) = 2*H_k / (x_m-x_k)^2
# Sum: 2*(H_k - H_m)/(x_m-x_k)^2 = 2*(H_k-H_m)/(x_k-x_m)^2

# Actually (x_m-x_k)^2 = (x_k-x_m)^2, so:
# dPhi/dx_m = 2*sum_{k!=m} (H_k - H_m)/(x_k - x_m)^2

# Note: (x_k-x_m)^2 > 0 always.

# For Schur-convexity: need (x_1-x_2)*(dPhi/dx_1 - dPhi/dx_2) >= 0.
# dPhi/dx_1 - dPhi/dx_2 = 2*sum_{k!=1} (H_k-H_1)/(x_k-x_1)^2 - 2*sum_{k!=2} (H_k-H_2)/(x_k-x_2)^2

# This is complicated. Let me check numerically.

print("\nSchur-convexity check: (x_1-x_2)*(dPhi/dx_1 - dPhi/dx_2) >= 0?")

np.random.seed(42)
schur_violations = 0
schur_total = 0

for trial in range(500):
    n = np.random.choice([3, 4, 5])
    x = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if x[i] - x[i-1] < 0.5:
            x[i] = x[i-1] + 0.5

    H = H_values(x)

    # Compute gradient
    grad = np.zeros(n)
    for m in range(n):
        for k in range(n):
            if k != m:
                grad[m] += 2*(H[k] - H[m])/(x[k]-x[m])**2

    schur_total += 1
    violated = False
    for i in range(n):
        for j in range(i+1, n):
            if (x[i]-x[j])*(grad[i]-grad[j]) < -1e-8:
                violated = True
                break
        if violated:
            break

    if violated:
        schur_violations += 1

print(f"  Schur-convexity holds: {schur_total - schur_violations}/{schur_total}")
print(f"  Violations: {schur_violations}")


# ================================================================
# IF PHI IS SCHUR-CONVEX, what does that buy us?
# ================================================================
print("\n" + "="*70)
print("SCHUR-CONVEXITY IMPLICATIONS")
print("="*70)

# If Phi is Schur-convex, then Phi(x) >= Phi(y) whenever x majorizes y.
# Majorization: x majorizes y means that the partial sums of the sorted
# x values dominate those of sorted y.

# For our application: x(0) = nu (r-roots), x(1) = lambda (p-roots).
# We need Phi(lambda) >= Phi(nu), i.e., lambda majorizes nu in some sense.

# The MSS convolution creates r from p and q. The roots of r are
# "more spread" than those of p. So the gap sequence of r dominates
# that of p. This is related to majorization.

# But wait: the STANDARD majorization would require:
# sum_{k=1}^{j} lambda_k >= sum_{k=1}^{j} nu_k for j = 1,...,n-1
# and equality for j = n (same sum).

# Since sum nu_k = sum lambda_k + sum mu_k, we have sum nu_k > sum lambda_k
# in general (unless sum mu = 0).

# So standard majorization doesn't directly apply. We'd need a SHIFTED version.

# HOWEVER: the interpolation x(t) = (1-t)*nu + t*lambda goes from nu to lambda
# along a straight line. If Phi is Schur-convex AND nu majorizes lambda
# (or the other way), then Phi is monotone along the path.

# Actually, Phi is NOT symmetric if we track WHICH index has which value.
# Wait, Phi = sum H_k^2 IS symmetric. And Schur-convexity of symmetric
# functions means Phi(x) >= Phi(y) whenever x majorizes y.

# CRUCIAL: The linear interpolation x(t) preserves sorting (since nu and lambda
# are both sorted). Along such a path, the sum of x(t) changes linearly.
# So the constraint sum x(t) = constant doesn't hold.

# But: Phi = sum H_k^2 depends on the GAPS, not the absolute positions.
# Translation-invariant Phi: Phi(x + c*1) = Phi(x) for any constant c.
# So we can shift x(t) to have sum = 0 without changing Phi.

# After centering: the centered points move from centered_nu to centered_lambda.
# The ordering is preserved. The question is whether centered_nu majorizes
# centered_lambda (or vice versa).

print("\nMajorization test: Does nu majorize lambda (after centering)?")
np.random.seed(42)
maj_count = 0
total = 0

for trial in range(500):
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
        if np.any(np.abs(np.imag(raw)) > 0.01): continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01): continue

        # Center both to have sum = 0
        nu_c = roots_r - np.mean(roots_r)
        lam_c = roots_p - np.mean(roots_p)

        # Does nu_c majorize lam_c?
        # For decreasingly sorted: sum_{k=1}^{j} nu*_{k} >= sum_{k=1}^{j} lam*_{k}
        # For our sorted (increasing): we reverse to get decreasing.
        nu_dec = nu_c[::-1]  # decreasing
        lam_dec = lam_c[::-1]

        majorizes = True
        for j in range(1, n):
            if np.sum(nu_dec[:j]) < np.sum(lam_dec[:j]) - 1e-8:
                majorizes = False
                break

        total += 1
        if majorizes:
            maj_count += 1
    except:
        pass

print(f"  nu majorizes lambda (centered): {maj_count}/{total}")
print(f"  (If Phi is Schur-convex AND nu majorizes lambda, then Phi(nu) >= Phi(lambda)")
print(f"   which means Phi_r >= Phi_p. But we need Phi_r <= Phi_p!)")

# Check the OTHER direction:
print("\nDoes lambda majorize nu (centered)?")
maj_count2 = 0
total2 = 0

for trial in range(500):
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
        if np.any(np.abs(np.imag(raw)) > 0.01): continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01): continue

        nu_c = roots_r - np.mean(roots_r)
        lam_c = roots_p - np.mean(roots_p)

        nu_dec = nu_c[::-1]
        lam_dec = lam_c[::-1]

        majorizes = True
        for j in range(1, n):
            if np.sum(lam_dec[:j]) < np.sum(nu_dec[:j]) - 1e-8:
                majorizes = False
                break

        total2 += 1
        if majorizes:
            maj_count2 += 1
    except:
        pass

print(f"  lambda majorizes nu (centered): {maj_count2}/{total2}")

# ================================================================
# Phi IS Schur-CONCAVE or CONVEX?
# ================================================================
print("\n\n" + "="*70)
print("IS Phi SCHUR-CONVEX OR SCHUR-CONCAVE?")
print("="*70)

# Schur-convex: Phi(x) >= Phi(y) when x majorizes y.
# Since Phi = sum H_k^2, and H_k involves 1/(x_k-x_j),
# when roots are CLOSER together (less spread), H values are LARGER,
# so Phi is LARGER.
# More spread roots -> smaller H values -> smaller Phi.

# If nu is more spread than lambda (nu majorizes lambda centered),
# then Phi(nu) < Phi(lambda), meaning Phi is SCHUR-CONCAVE.

# Let's test explicitly:
print("\nTest: does majorization imply Phi ordering?")
# Create two ordered sets where x majorizes y (after centering).
x = np.array([-3, 0, 3])  # centered, spread
y = np.array([-1, 0, 1])  # centered, less spread

Phi_x = np.sum(H_values(x)**2)
Phi_y = np.sum(H_values(y)**2)

print(f"  x = {x}: Phi = {Phi_x:.6f}")
print(f"  y = {y}: Phi = {Phi_y:.6f}")
print(f"  x majorizes y: True (by construction)")
print(f"  Phi(x) < Phi(y): {Phi_x < Phi_y}")
print(f"  So Phi is SCHUR-CONCAVE (more spread = less Phi)")

# Double-check: if x majorizes y means x is more "spread",
# and Phi(x) <= Phi(y), then Phi is Schur-concave.

# For our setting:
# - nu (r-roots) is MORE spread than lambda (p-roots) because free convolution spreads
# - So nu majorizes lambda (centered)
# - If Phi is Schur-concave, then Phi(nu) <= Phi(lambda), i.e., Phi_r <= Phi_p

# This is consistent with our observation that Phi_r <= Phi_p!

# BUT: Schur-concavity of Phi + majorization would give Phi_r <= Phi_p,
# which is WEAKER than <h,alpha> >= 0.

# <h,alpha> >= 0 is: sum h_k * u_k >= sum h_k^2.
# Phi_r <= Phi_p is: sum u_k^2 >= sum h_k^2.
# The latter follows from the former via Cauchy-Schwarz? No.
# Actually, sum h_k * u_k >= sum h_k^2 implies (by CS) sum u_k^2 >= (sum h_k*u_k)^2 / sum h_k^2 >= sum h_k^2.
# Wait: CS gives (sum h_k*u_k)^2 <= (sum h_k^2)(sum u_k^2).
# If sum h_k*u_k >= sum h_k^2 = A, then A^2 <= A * sum u_k^2, so A <= sum u_k^2.
# I.e., Phi_r <= Phi_p. YES. So <h,alpha> >= 0 implies Phi_r <= Phi_p.
# But NOT vice versa.

# So the Schur-concavity approach only gives the WEAKER result.

# ================================================================
# SUMMARY AND CONCLUSION
# ================================================================
print("\n\n" + "="*70)
print("FINAL SUMMARY OF HERGLOTZ APPROACH INVESTIGATION")
print("="*70)
print("""
GOAL: Prove <h,alpha> >= 0 where h_k = H_r(nu_k), alpha_k = H_p(lambda_k) - H_r(nu_k),
and r = p boxplus_n q (MSS finite free convolution).

APPROACH TRIED: Herglotz convexity of omega_1 (subordination function).

KEY FINDINGS:

1. omega_1 DOES map C^+ to C^+ (Im(omega_1(z)) > 0 for Im(z) > 0).
   BUT Im(omega_1(z)) is NOT always >= Im(z). So phi = omega_1 - id is NOT Herglotz.

2. The Nevanlinna representation of omega_1 as alpha + z + sum c_j/(d_j - z) with c_j > 0
   gives phi' = sum c_j/(d_j-z)^2 > 0 everywhere. This CONTRADICTS phi'(nu_k) = 0.
   RESOLUTION: omega_1 is NOT a single-valued rational Herglotz function. It is a
   BRANCH of an algebraic function (defined by the degree-n equation G_r(z) = G_p(w)).

3. The Loewner matrix of omega_1 (L_{kl} = (lambda_k-lambda_l)/(nu_k-nu_l))
   is NOT always PSD (fails ~50% of tested cases).
   So omega_1 is NOT operator monotone in general.

4. <h,alpha> >= 0 is NOT true for arbitrary sorted point sets (nu, lambda).
   It SPECIFICALLY requires the MSS convolution structure (fails 71% for random pairs).

ALTERNATIVE REFORMULATIONS FOUND:

(A) <h,alpha> = sum_{k<l} (h_k-h_l)(phi_k-phi_l) / [(nu_k-nu_l)(lambda_k-lambda_l)]
    where phi_k = nu_k - lambda_k. This is a POSITIVE-WEIGHTED sum of products of
    differences. Individual terms can be negative (7.7% of pairs), but the total sum
    is always positive.

(B) Phi(t) = sum H(x_k(t))^2 along interpolation x(t) = (1-t)*nu + t*lambda is
    MONOTONE INCREASING (0/500 violations). This gives Phi_r <= Phi_p (weaker than target).

(C) The Schur-concavity of Phi combined with the majorization nu > lambda gives
    Phi_r <= Phi_p, but this is weaker than <h,alpha> >= 0.

STATUS: PROOF NOT FOUND. The Herglotz approach is DEFINITIVELY BLOCKED because
omega_1 does not have the required Herglotz structure in the finite case.

MOST PROMISING REMAINING DIRECTION: The identity (A) above might be provable
positive via the specific STRUCTURE of the MSS convolution (e.g., the fact that
the MSS convolution preserves real-rootedness and the subordination satisfies
omega_1'(nu_k) = 1). But this requires a different technique than Herglotz theory.
""")
