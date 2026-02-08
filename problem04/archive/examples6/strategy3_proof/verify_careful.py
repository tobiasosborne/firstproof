#!/usr/bin/env python3
"""
Careful investigation of the apparent violation.
The key question: for nearly degenerate polynomials, is the
direct superadditivity 1/Phi(r) >= 1/Phi(p) + 1/Phi(q) actually violated,
or is this a numerical artifact from inaccurate root-finding of r?
"""

import numpy as np
from math import factorial
import mpmath

mpmath.mp.dps = 100  # 100 decimal digits

np.random.seed(42)

def boxplus_n_mp(p_coeffs, q_coeffs, n):
    """Compute boxplus using mpmath for arbitrary precision."""
    a = [mpmath.mpf(str(x)) for x in p_coeffs]
    b = [mpmath.mpf(str(x)) for x in q_coeffs]
    r = [mpmath.mpf(0)] * (n + 1)
    r[0] = mpmath.mpf(1)
    for k in range(1, n + 1):
        c_k = mpmath.mpf(0)
        for i in range(0, k + 1):
            j = k - i
            if i <= n and j <= n:
                coeff = (mpmath.factorial(n - i) * mpmath.factorial(n - j)) / (mpmath.factorial(n) * mpmath.factorial(n - k))
                c_k += coeff * a[i] * b[j]
        r[k] = c_k
    return r

def poly_from_roots_mp(roots):
    """Compute polynomial coefficients from roots using mpmath."""
    coeffs = [mpmath.mpf(1)]
    for r in roots:
        new = [mpmath.mpf(0)] * (len(coeffs) + 1)
        for i, c in enumerate(coeffs):
            new[i] += c
            new[i+1] -= c * r
        coeffs = new
    return coeffs

def H_values_mp(roots):
    """Compute H_p(lambda_i) using mpmath."""
    n = len(roots)
    H = [mpmath.mpf(0)] * n
    for i in range(n):
        for j in range(n):
            if i != j:
                H[i] += mpmath.mpf(1) / (roots[i] - roots[j])
    return H

def Phi_n_mp(roots):
    """Compute Phi_n using mpmath."""
    H = H_values_mp(roots)
    return sum(h**2 for h in H)

# ============================================================
# Test the specific violation case
# ============================================================
print("=" * 60)
print("HIGH PRECISION INVESTIGATION OF VIOLATION")
print("=" * 60)

n = 4
p_roots_f = [4.56501398, 4.56613507, 4.56725615, 4.56837724]
q_roots_f = [2.37508359, 2.3763709, 2.37771781, 2.37974953]

p_roots = [mpmath.mpf(str(x)) for x in p_roots_f]
q_roots = [mpmath.mpf(str(x)) for x in q_roots_f]

print(f"p_roots = {[float(x) for x in p_roots]}")
print(f"q_roots = {[float(x) for x in q_roots]}")

# Compute coefficients
p_coeffs = poly_from_roots_mp(p_roots)
q_coeffs = poly_from_roots_mp(q_roots)

print(f"\np_coeffs = {[float(x) for x in p_coeffs]}")
print(f"q_coeffs = {[float(x) for x in q_coeffs]}")

# Compute boxplus
r_coeffs = boxplus_n_mp(p_coeffs, q_coeffs, n)
print(f"r_coeffs = {[float(x) for x in r_coeffs]}")

# Find roots of r using mpmath with high precision
try:
    r_roots = mpmath.polyroots(r_coeffs, maxsteps=500, cleanup=True, error=True)
    r_roots_vals = sorted([r[0] for r in zip(r_roots[0], r_roots[0])], key=lambda x: mpmath.re(x))
    print(f"\nr_roots = {[str(mpmath.nstr(x, 15)) for x in r_roots[0]]}")
    print(f"r_root errors = {[str(mpmath.nstr(x, 5)) for x in r_roots[1]]}")

    r_roots_sorted = sorted(r_roots[0], key=lambda x: float(mpmath.re(x)))
    r_roots_real = [mpmath.re(x) for x in r_roots_sorted]
    r_imag = [mpmath.im(x) for x in r_roots_sorted]
    print(f"Imaginary parts: {[float(x) for x in r_imag]}")

    if any(abs(float(x)) > 1e-10 for x in r_imag):
        print("WARNING: r has complex roots! The polynomial may not be real-rooted.")
    else:
        # Compute Phi values
        phi_p = Phi_n_mp(p_roots)
        phi_q = Phi_n_mp(q_roots)
        phi_r = Phi_n_mp(r_roots_real)

        print(f"\nPhi(p) = {float(phi_p):.15e}")
        print(f"Phi(q) = {float(phi_q):.15e}")
        print(f"Phi(r) = {float(phi_r):.15e}")

        A = phi_p - phi_r
        B = phi_q - phi_r
        print(f"A = Phi(p) - Phi(r) = {float(A):.15e}")
        print(f"B = Phi(q) - Phi(r) = {float(B):.15e}")

        if phi_r > 0:
            ratio = A * B / phi_r**2
            print(f"ratio = AB/Phi(r)^2 = {float(ratio):.15e}")

        superadd = 1/phi_r - 1/phi_p - 1/phi_q
        print(f"1/Phi(r) - 1/Phi(p) - 1/Phi(q) = {float(superadd):.15e}")
        print(f"Superadditivity: {'HOLDS' if superadd >= 0 else 'VIOLATED'}")

except Exception as e:
    print(f"Root finding failed: {e}")
    print("Trying alternative: compute r_roots from the formula directly...")

    # Alternative: use numpy roots but refine with mpmath
    r_coeffs_np = [float(x) for x in r_coeffs]
    r_roots_np = np.sort(np.real(np.roots(r_coeffs_np)))
    print(f"numpy r_roots = {r_roots_np}")

    # Refine each root using Newton's method in mpmath
    from mpmath import polyval as mp_polyval

    def eval_poly_mp(coeffs, x):
        """Evaluate polynomial with mpmath coefficients."""
        result = mpmath.mpf(0)
        for c in coeffs:
            result = result * x + c
        return result

    def eval_dpoly_mp(coeffs, x):
        """Evaluate derivative of polynomial."""
        n = len(coeffs) - 1
        result = mpmath.mpf(0)
        for i, c in enumerate(coeffs[:-1]):
            result = result * x + c * (n - i)
        return result

    r_roots_refined = []
    for root_init in r_roots_np:
        x = mpmath.mpf(str(root_init))
        for _ in range(200):
            fx = eval_poly_mp(r_coeffs, x)
            fpx = eval_dpoly_mp(r_coeffs, x)
            if fpx == 0:
                break
            x_new = x - fx / fpx
            if abs(x_new - x) < mpmath.mpf(10)**(-90):
                break
            x = x_new
        r_roots_refined.append(x)
        residual = eval_poly_mp(r_coeffs, x)
        print(f"  Root {float(root_init):.10f} -> {mpmath.nstr(x, 20)}, residual = {mpmath.nstr(residual, 5)}")

    r_roots_refined = sorted(r_roots_refined, key=lambda x: float(x))

    phi_p = Phi_n_mp(p_roots)
    phi_q = Phi_n_mp(q_roots)
    phi_r = Phi_n_mp(r_roots_refined)

    print(f"\nHIGH PRECISION RESULTS:")
    print(f"Phi(p) = {mpmath.nstr(phi_p, 30)}")
    print(f"Phi(q) = {mpmath.nstr(phi_q, 30)}")
    print(f"Phi(r) = {mpmath.nstr(phi_r, 30)}")

    A = phi_p - phi_r
    B = phi_q - phi_r
    print(f"A = {mpmath.nstr(A, 30)}")
    print(f"B = {mpmath.nstr(B, 30)}")

    if phi_r > 0:
        ratio = A * B / phi_r**2
        print(f"ratio = {mpmath.nstr(ratio, 30)}")

    superadd = 1/phi_r - 1/phi_p - 1/phi_q
    print(f"superadditivity excess = {mpmath.nstr(superadd, 30)}")
    print(f"Superadditivity: {'HOLDS' if superadd >= 0 else 'VIOLATED'}")

# ============================================================
# Now test with WELL-SEPARATED roots using high precision
# ============================================================
print("\n" + "=" * 60)
print("WELL-SEPARATED ROOTS, HIGH PRECISION CHECK")
print("=" * 60)

for n in [3, 4, 5]:
    np.random.seed(42)
    min_ratio = mpmath.inf
    violations = 0
    valid = 0

    for trial in range(200):
        p_roots_f = sorted(np.random.uniform(-5, 5, n))
        q_roots_f = sorted(np.random.uniform(-5, 5, n))

        # Ensure gaps >= 0.5
        for i in range(1, n):
            if p_roots_f[i] - p_roots_f[i-1] < 0.5:
                p_roots_f[i] = p_roots_f[i-1] + 0.5
            if q_roots_f[i] - q_roots_f[i-1] < 0.5:
                q_roots_f[i] = q_roots_f[i-1] + 0.5

        p_roots = [mpmath.mpf(str(x)) for x in p_roots_f]
        q_roots = [mpmath.mpf(str(x)) for x in q_roots_f]

        p_coeffs = poly_from_roots_mp(p_roots)
        q_coeffs = poly_from_roots_mp(q_roots)
        r_coeffs = boxplus_n_mp(p_coeffs, q_coeffs, n)

        # Find roots using numpy then refine with Newton
        r_coeffs_np = [float(x) for x in r_coeffs]
        r_roots_np = np.roots(r_coeffs_np)
        if np.max(np.abs(np.imag(r_roots_np))) > 1e-6:
            continue

        r_roots_np = np.sort(np.real(r_roots_np))

        def eval_poly_mp(coeffs, x):
            result = mpmath.mpf(0)
            for c in coeffs:
                result = result * x + c
            return result

        def eval_dpoly_mp(coeffs, x):
            nn = len(coeffs) - 1
            result = mpmath.mpf(0)
            for i, c in enumerate(coeffs[:-1]):
                result = result * x + c * (nn - i)
            return result

        r_roots_hp = []
        for root_init in r_roots_np:
            x = mpmath.mpf(str(root_init))
            for _ in range(100):
                fx = eval_poly_mp(r_coeffs, x)
                fpx = eval_dpoly_mp(r_coeffs, x)
                if fpx == 0: break
                x_new = x - fx / fpx
                if abs(x_new - x) < mpmath.mpf(10)**(-90): break
                x = x_new
            r_roots_hp.append(x)

        r_roots_hp = sorted(r_roots_hp, key=lambda x: float(x))

        if len(r_roots_hp) < 2:
            continue
        min_gap = min(float(r_roots_hp[i+1] - r_roots_hp[i]) for i in range(len(r_roots_hp)-1))
        if min_gap < 1e-10:
            continue

        phi_p = Phi_n_mp(p_roots)
        phi_q = Phi_n_mp(q_roots)
        phi_r = Phi_n_mp(r_roots_hp)

        A = phi_p - phi_r
        B = phi_q - phi_r

        valid += 1

        if phi_r > 0:
            ratio = A * B / phi_r**2
            if ratio < min_ratio:
                min_ratio = ratio

            if ratio < mpmath.mpf(1) - mpmath.mpf(10)**(-10):
                violations += 1
                print(f"  VIOLATION n={n} trial {trial}: ratio={float(ratio):.15e}")

        superadd = 1/phi_r - 1/phi_p - 1/phi_q
        if superadd < -mpmath.mpf(10)**(-10):
            print(f"  SUPERADD VIOLATION n={n} trial {trial}: excess={float(superadd):.15e}")

    print(f"  n={n}: min_ratio={float(min_ratio):.15e}, violations={violations}/{valid}")
