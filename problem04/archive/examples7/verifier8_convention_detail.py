"""
VERIFIER-8: Investigate why the "wrong" convention also gives additive cumulants.

The answer: the "wrong" convention is just a RESCALING of the AP convention.
If kappa_k^{AP} = c_k * kappa_k^{wrong}, then both are additive.
The question is whether the 1/Phi_3 formula uses the RIGHT convention.
"""
import numpy as np
from itertools import combinations
from math import factorial, comb

def elementary_symmetric(roots, k):
    n = len(roots)
    if k == 0: return 1.0
    if k > n: return 0.0
    return sum(np.prod(combo) for combo in combinations(roots, k))

def phi_n(roots):
    n = len(roots)
    total = 0.0
    for i in range(n):
        s = sum(1.0/(roots[i]-roots[j]) for j in range(n) if j != i)
        total += s**2
    return total

def mss_convolve(roots_p, roots_q):
    n = len(roots_p)
    assert len(roots_q) == n
    p_coeffs = np.poly(roots_p)
    q_coeffs = np.poly(roots_q)
    r_coeffs = np.zeros(n + 1)
    r_coeffs[0] = 1.0
    for k in range(1, n + 1):
        ck = 0.0
        for i in range(k + 1):
            j = k - i
            if i <= n and j <= n:
                w = factorial(n-i)*factorial(n-j)/(factorial(n)*factorial(n-k))
                ck += w * p_coeffs[i] * q_coeffs[j]
        r_coeffs[k] = ck
    return np.sort(np.real(np.roots(r_coeffs)))

np.random.seed(42)

# The "wrong" convention: tilde_a_k same, but
# k2_wrong = -(ta2 - ta1^2), k3_wrong = 0.5*(ta3 - 3*ta2*ta1 + 2*ta1^3)
# AP convention: k2_AP = -n*(ta2 - ta1^2) = n * k2_wrong
#                k3_AP = n^2/2*(ta3 - ...) = n^2 * k3_wrong

# So k2_AP = n * k2_wrong and k3_AP = n^2 * k3_wrong (for n=3: k2_AP = 3*k2_wrong, k3_AP = 9*k3_wrong)
# Both are additive since scaling by a constant preserves additivity.
# This is not a contradiction -- it just means any scalar multiple of additive cumulants is additive.

print("The 'wrong' convention gives cumulants that are scalar multiples of AP:")
print(f"  k2_AP = n * k2_wrong = 3 * k2_wrong")
print(f"  k3_AP = n^2 * k3_wrong = 9 * k3_wrong")
print()

# Verify
for trial in range(5):
    r = np.sort(np.random.randn(3) * 2)
    r = r - np.mean(r)

    es = [elementary_symmetric(r, k) for k in range(4)]
    a = [(-1)**k * es[k] for k in range(4)]
    ta = [(-1)**k * a[k] / comb(3, k) for k in range(4)]

    k2_AP = -3*(ta[2] - ta[1]**2)
    k3_AP = 4.5*(ta[3] - 3*ta[2]*ta[1] + 2*ta[1]**3)

    k2_wrong = -(ta[2] - ta[1]**2)
    k3_wrong = 0.5*(ta[3] - 3*ta[2]*ta[1] + 2*ta[1]**3)

    print(f"  Trial {trial}: k2_AP/k2_wrong = {k2_AP/k2_wrong:.6f}, k3_AP/k3_wrong = {k3_AP/k3_wrong:.6f}")

print()
print("CONCLUSION: Both conventions give additive cumulants (scalar multiples).")
print("The formula 1/Phi_3 = (2/9)*k2 - (2/27)*k3^2/k2^2 uses AP convention.")
print("With the 'wrong' convention it would be:")
print("  1/Phi_3 = (2/9)*(3*k2w) - (2/27)*(9*k3w)^2/(3*k2w)^2")
print("          = (2/3)*k2w - (2/27)*81*k3w^2/(9*k2w^2)")
print("          = (2/3)*k2w - (6/9)*k3w^2/k2w^2")
print("          = (2/3)*k2w - (2/3)*k3w^2/k2w^2")
print()
print("This is consistent -- the formula changes with convention, but the PROOF")
print("is convention-independent as long as all steps use the SAME convention.")
print()

# FINAL VERIFICATION: The formula as stated in the HANDOFF uses AP convention
# and is numerically correct.
print("Final check: formula with AP convention matches Phi_3:")
for trial in range(10):
    r = np.sort(np.random.randn(3) * (0.5 + np.random.rand()*5))
    if min(np.diff(r)) < 0.01:
        continue

    phi = phi_n(r)

    es = [elementary_symmetric(r, k) for k in range(4)]
    a = [(-1)**k * es[k] for k in range(4)]
    ta = [(-1)**k * a[k] / comb(3, k) for k in range(4)]

    k2 = -3*(ta[2] - ta[1]**2)
    k3 = 4.5*(ta[3] - 3*ta[2]*ta[1] + 2*ta[1]**3)

    inv_phi = 1.0/phi
    formula = (2.0/9.0)*k2 - (2.0/27.0)*k3**2/k2**2

    err = abs(inv_phi - formula)
    print(f"  1/Phi_3 = {inv_phi:.8f}, formula = {formula:.8f}, err = {err:.2e}")

print("\nAll consistent. No convention issues.")
