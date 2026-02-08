#!/usr/bin/env python3
"""
Part 5: Test the F(t) monotonicity approach more thoroughly.

KEY FINDING: F(t) = <h_r, H_t(x(t))> where x_k(t) = (1-t)*nu_k + t*lambda_k
appears to be monotone increasing in t.

F(0) = ||h||^2 = Phi_r
F(1) = <h, u> = Phi_r + <h, alpha>

If we can prove F'(t) >= 0 for all t in [0,1], then <h,alpha> >= 0.

F'(t) = sum_k h_k * d/dt[H_t(x_k(t))]

Let's derive F'(t) carefully.

x_k(t) = (1-t)*nu_k + t*lambda_k
x_k'(t) = lambda_k - nu_k = -phi_k (where phi = nu - lambda)

H_t(x_k) = sum_{l != k} 1/(x_k(t) - x_l(t))

d/dt H_t(x_k) = sum_{l != k} -(x_k' - x_l') / (x_k - x_l)^2
              = sum_{l != k} -(-phi_k + phi_l) / (x_k(t) - x_l(t))^2
              = sum_{l != k} (phi_k - phi_l) / (x_k(t) - x_l(t))^2

F'(t) = sum_k h_k * sum_{l != k} (phi_k - phi_l) / (x_k(t) - x_l(t))^2
       = sum_{k != l} h_k * (phi_k - phi_l) / (x_k(t) - x_l(t))^2
       = sum_{k < l} (h_k - h_l)(phi_k - phi_l) / (x_k(t) - x_l(t))^2

Wait, that's NOT quite right. Let me redo the antisymmetrization:
sum_{k != l} h_k * (phi_k - phi_l) / (x_k - x_l)^2
= sum_{k < l} [h_k*(phi_k-phi_l) + h_l*(phi_l-phi_k)] / (x_k-x_l)^2
= sum_{k < l} (h_k - h_l)(phi_k - phi_l) / (x_k(t) - x_l(t))^2

NOTE: Here h_k is FIXED at H_r(nu_k), but the denominator (x_k(t)-x_l(t))^2
varies with t!

At t = 0: (x_k-x_l)^2 = (nu_k-nu_l)^2
At t = 1: (x_k-x_l)^2 = (lambda_k-lambda_l)^2

And x_k(t) - x_l(t) = (1-t)(nu_k-nu_l) + t(lambda_k-lambda_l).

So F'(t) = sum_{k<l} (h_k-h_l)(phi_k-phi_l) / [(1-t)(nu_k-nu_l) + t(lambda_k-lambda_l)]^2

The denominator is always positive (since both nu and lambda are sorted).

IMPORTANT: The NUMERATOR (h_k-h_l)(phi_k-phi_l) does NOT depend on t,
but the DENOMINATOR does. So F'(t) is a WEIGHTED sum with t-dependent weights.

For F'(t) >= 0: we need the positive terms to dominate, weighted by 1/gap(t)^2.
As t varies, the weights change, but the signs of individual terms don't.

Let me verify F'(t) >= 0 for all t numerically with many trials.
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


def F_prime(t, h, phi, nu, lam):
    """Compute F'(t) = sum_{k<l} (h_k-h_l)(phi_k-phi_l) / gap(t)^2."""
    n = len(h)
    result = 0
    for k in range(n):
        for l in range(k+1, n):
            gap_t = (1-t)*(nu[k]-nu[l]) + t*(lam[k]-lam[l])
            if abs(gap_t) < 1e-15:
                return np.nan
            result += (h[k]-h[l])*(phi[k]-phi[l]) / gap_t**2
    return result


np.random.seed(42)

print("="*70)
print("TESTING F'(t) >= 0 FOR ALL t IN [0,1]")
print("="*70)

violations = 0
total = 0
min_Fprime = np.inf

for trial in range(1000):
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
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        raw = np.roots(np.poly(roots_r))
        if np.any(np.abs(np.imag(raw)) > 0.01): continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01): continue

        h = H_values(roots_r)
        phi = roots_r - roots_p

        total += 1
        local_min = np.inf

        for t in np.linspace(0, 1, 50):
            fp = F_prime(t, h, phi, roots_r, roots_p)
            if np.isnan(fp):
                continue
            if fp < local_min:
                local_min = fp
            if fp < min_Fprime:
                min_Fprime = fp

        if local_min < -1e-8:
            violations += 1

    except:
        pass

print(f"\nTotal valid trials: {total}")
print(f"Violations (F'(t) < 0 somewhere): {violations}")
print(f"Global min F'(t): {min_Fprime:.8e}")


# ================================================================
# DIFFERENT APPROACH: F(t) where h(t) also varies
# ================================================================
print("\n\n" + "="*70)
print("TESTING: Is Phi(t) = sum H_t(x_k(t))^2 monotone increasing?")
print("="*70)

violations_phi = 0
total_phi = 0

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

        total_phi += 1

        t_vals = np.linspace(0, 1, 50)
        Phi_vals = []
        for t in t_vals:
            x_t = (1-t)*roots_r + t*roots_p
            h_t = H_values(x_t)
            Phi_vals.append(np.sum(h_t**2))

        dPhi = np.diff(Phi_vals)
        if np.any(dPhi < -1e-8):
            violations_phi += 1

    except:
        pass

print(f"Total: {total_phi}")
print(f"Phi(t) NOT monotone: {violations_phi}")
print(f"Phi(t) monotone increasing: {total_phi - violations_phi}/{total_phi}")


# ================================================================
# KEY TEST: Does F'(t) >= 0 hold EVEN for non-MSS pairs (nu, lambda)?
# ================================================================
print("\n\n" + "="*70)
print("CONTROL: Does F'(t) >= 0 for RANDOM (non-MSS) pairs?")
print("="*70)

# If F'(t) >= 0 holds for ANY pair of sorted point sets nu, lambda with
# h = H(nu) and phi = nu - lambda, then the MSS structure is irrelevant!
# Let's test with random pairs NOT arising from convolution.

np.random.seed(42)
ctrl_violations = 0
ctrl_total = 0

for trial in range(1000):
    n = np.random.choice([3, 4, 5])
    nu = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if nu[i] - nu[i-1] < 0.5:
            nu[i] = nu[i-1] + 0.5

    lam = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if lam[i] - lam[i-1] < 0.5:
            lam[i] = lam[i-1] + 0.5

    h = H_values(nu)
    phi = nu - lam

    ctrl_total += 1
    violated = False
    for t in np.linspace(0, 1, 50):
        fp = F_prime(t, h, phi, nu, lam)
        if np.isnan(fp): continue
        if fp < -1e-8:
            violated = True
            break

    if violated:
        ctrl_violations += 1

print(f"Control total: {ctrl_total}")
print(f"Control violations: {ctrl_violations}")
print(f"Percentage violated: {100*ctrl_violations/ctrl_total:.1f}%")


# Also test <h,alpha> >= 0 for random pairs:
np.random.seed(42)
ctrl_ha_violations = 0
ctrl_ha_total = 0

for trial in range(1000):
    n = np.random.choice([3, 4, 5])
    nu = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if nu[i] - nu[i-1] < 0.5:
            nu[i] = nu[i-1] + 0.5

    lam = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if lam[i] - lam[i-1] < 0.5:
            lam[i] = lam[i-1] + 0.5

    h = H_values(nu)
    u = H_values(lam)
    alpha = u - h

    ctrl_ha_total += 1
    if np.dot(h, alpha) < -1e-8:
        ctrl_ha_violations += 1

print(f"\nControl <h,alpha> >= 0: violations = {ctrl_ha_violations}/{ctrl_ha_total}")
print(f"Percentage violated: {100*ctrl_ha_violations/ctrl_ha_total:.1f}%")

if ctrl_ha_violations > 0:
    print("\n*** <h,alpha> >= 0 is NOT universally true for arbitrary sorted point sets! ***")
    print("This means the MSS convolution structure IS needed for the proof.")
else:
    print("\n*** <h,alpha> >= 0 appears universal! MSS structure may NOT be needed. ***")
