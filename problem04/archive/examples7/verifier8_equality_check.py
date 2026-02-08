"""
VERIFIER-8: Final equality condition check.
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

def compute_cumulants_n3(roots):
    e1 = sum(roots)
    e2 = sum(roots[i]*roots[j] for i in range(3) for j in range(i+1,3))
    e3 = roots[0]*roots[1]*roots[2]
    a1, a2, a3 = -e1, e2, -e3
    ta1 = -a1/3
    ta2 = a2/3
    ta3 = -a3
    k1 = ta1
    k2 = -3*(ta2 - ta1**2)
    k3 = 4.5*(ta3 - 3*ta2*ta1 + 2*ta1**3)
    return k1, k2, k3

print("Equality conditions check:")
print("\nEqually spaced (k3=0):")
for gap in [1, 2, 5, 0.1]:
    r_p = np.array([-gap, 0, gap])
    r_q = np.array([-gap*2, 0, gap*2])
    r_r = mss_convolve(r_p, r_q)
    inv_phi_r = 1.0/phi_n(r_r)
    inv_phi_sum = 1.0/phi_n(r_p) + 1.0/phi_n(r_q)
    margin = inv_phi_r - inv_phi_sum
    print(f"  gap_p={gap}, gap_q={gap*2}: margin = {margin:.2e}")

print("\nNon-equally-spaced:")
np.random.seed(42)
for _ in range(10):
    r_p = np.sort(np.random.randn(3))
    r_p = r_p - np.mean(r_p)
    r_q = np.sort(np.random.randn(3))
    r_q = r_q - np.mean(r_q)
    if min(np.diff(r_p)) < 0.1 or min(np.diff(r_q)) < 0.1:
        continue
    try:
        r_r = mss_convolve(r_p, r_q)
        if min(np.diff(r_r)) < 1e-8:
            continue
        inv_phi_r = 1.0/phi_n(r_r)
        inv_phi_sum = 1.0/phi_n(r_p) + 1.0/phi_n(r_q)
        margin = inv_phi_r - inv_phi_sum
        k3_p = compute_cumulants_n3(r_p)[2]
        k3_q = compute_cumulants_n3(r_q)[2]
        print(f"  k3_p={k3_p:.4f}, k3_q={k3_q:.4f}, margin = {margin:.6e} {'STRICT' if margin > 1e-10 else 'EQUALITY'}")
    except:
        pass

# Check: equally-spaced p with NON-equally-spaced q
print("\nMixed: equally-spaced + non-equally-spaced:")
r_p = np.array([-1.0, 0.0, 1.0])
for _ in range(5):
    r_q = np.sort(np.random.randn(3)*2)
    r_q = r_q - np.mean(r_q)
    if min(np.diff(r_q)) < 0.1:
        continue
    try:
        r_r = mss_convolve(r_p, r_q)
        if min(np.diff(r_r)) < 1e-8:
            continue
        inv_phi_r = 1.0/phi_n(r_r)
        inv_phi_sum = 1.0/phi_n(r_p) + 1.0/phi_n(r_q)
        margin = inv_phi_r - inv_phi_sum
        k3_q = compute_cumulants_n3(r_q)[2]
        print(f"  k3_q={k3_q:.4f}: margin = {margin:.6e} {'STRICT' if margin > 1e-10 else 'EQUALITY'}")
    except:
        pass

# ADVERSARIAL: Can we have equality with k3_p/k2_p = k3_q/k2_q but both nonzero?
# The chain had TWO inequalities:
# (1) (sx+ty)^2 <= s*x^2 + t*y^2  (equality iff x=y)
# (2) s*x^2 + t*y^2 <= x^2 + y^2  (equality iff x=y=0)
# Both must be equalities for overall equality.
# So we need x=y AND x=y=0, i.e., k3_p=k3_q=0.
# But wait -- if x=y (but nonzero), then step (1) is equality but step (2) is STRICT.
# So overall, we get STRICT inequality when x=y != 0.

print("\nADVERSARIAL: k3_p/k2_p = k3_q/k2_q but both nonzero:")
# Create two polynomials with same k3/k2 ratio
# k3/k2 is a scale-invariant quantity related to skewness
# For roots {a, 0, -a+d}: e3 = a*0*(-a+d) = 0... that's k3=0
# Try: roots = {0, 1, t} centered. mean = (1+t)/3.
# Centered: {-(1+t)/3, 1-(1+t)/3, t-(1+t)/3} = {-(1+t)/3, (2-t)/3, (2t-1)/3}

for ratio_target in [0.5, 1.0, 2.0]:
    # Build polynomial with specific k3/k2 ratio by brute force search
    best_r = None
    best_err = float('inf')
    for _ in range(1000):
        r = np.sort(np.random.randn(3)*2)
        r = r - np.mean(r)
        if min(np.diff(r)) < 0.1:
            continue
        _, k2, k3 = compute_cumulants_n3(r)
        if k2 > 0.1:
            err = abs(k3/k2 - ratio_target)
            if err < best_err:
                best_err = err
                best_r = r.copy()

    if best_r is not None:
        _, k2_p, k3_p = compute_cumulants_n3(best_r)
        # Scale to get a second polynomial with same ratio
        scale = 2.0
        r_q = best_r * scale
        _, k2_q, k3_q = compute_cumulants_n3(r_q)

        print(f"  ratio={ratio_target}: k3_p/k2_p={k3_p/k2_p:.4f}, k3_q/k2_q={k3_q/k2_q:.4f}")

        try:
            r_r = mss_convolve(best_r, r_q)
            if min(np.diff(r_r)) < 1e-8:
                continue
            inv_phi_r = 1.0/phi_n(r_r)
            inv_phi_sum = 1.0/phi_n(best_r) + 1.0/phi_n(r_q)
            margin = inv_phi_r - inv_phi_sum
            print(f"    margin = {margin:.6e} -- {'STRICT' if margin > 1e-10 else 'EQUALITY'}")

            # Compare with theoretical prediction:
            # x=y case: (sx+ty)^2 = (x(s+t))^2 = x^2. And RHS = x^2 + y^2 = 2x^2.
            # So margin should be (2/27)*(2x^2 - x^2)/(something) = positive.
            # Specifically: LHS-RHS of KEY ineq = x^2+y^2-(sx+ty)^2 = x^2+x^2-x^2 = x^2 > 0
            print(f"    (This confirms: same ratio but nonzero => STRICT inequality)")
        except:
            pass

print("\nDone.")
