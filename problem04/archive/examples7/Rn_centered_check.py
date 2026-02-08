"""
Check R_4 and R_5 superadditivity with CENTERED polynomials.
"""
import numpy as np
from math import comb, factorial

np.random.seed(42)

def phi_n_num(roots):
    n = len(roots)
    total = 0.0
    for i in range(n):
        H = sum(1.0/(roots[i]-roots[j]) for j in range(n) if j != i)
        total += H**2
    return total

def mss_convolve(roots_p, roots_q):
    n = len(roots_p)
    p_coeffs = np.poly(roots_p)
    q_coeffs = np.poly(roots_q)
    r_coeffs = np.zeros(n + 1)
    r_coeffs[0] = 1.0
    for k in range(1, n + 1):
        ck = 0.0
        for i in range(k + 1):
            j = k - i
            if i <= n and j <= n:
                w = factorial(n - i) * factorial(n - j) / (factorial(n) * factorial(n - k))
                ck += w * p_coeffs[i] * q_coeffs[j]
        r_coeffs[k] = ck
    return np.sort(np.real(np.roots(r_coeffs)))

def get_ta(roots):
    n = len(roots)
    coeffs = np.poly(roots)
    ta = {}
    for k in range(1, n+1):
        ta[k] = ((-1)**k * coeffs[k]) / comb(n, k)
    return ta

# Test R_4 with centered polys
for n in [4, 5]:
    Cn = -2.0/comb(n,2)
    violations = 0
    valid = 0

    for trial in range(15000):
        rp = np.sort(np.random.randn(n) * 2)
        rp = rp - np.mean(rp)  # CENTER
        rq = np.sort(np.random.randn(n) * 2)
        rq = rq - np.mean(rq)  # CENTER

        if min(np.diff(rp)) < 0.12 or min(np.diff(rq)) < 0.12:
            continue

        try:
            rr = mss_convolve(rp, rq)
            phi_p = phi_n_num(rp); phi_q = phi_n_num(rq); phi_r = phi_n_num(rr)
            if phi_p <= 0 or phi_q <= 0 or phi_r <= 0:
                continue

            ta_p = get_ta(rp); ta_q = get_ta(rq); ta_r = get_ta(rr)

            R_p = 1/phi_p - Cn*ta_p[2]
            R_q = 1/phi_q - Cn*ta_q[2]
            R_r = 1/phi_r - Cn*ta_r[2]

            valid += 1
            if R_r < R_p + R_q - 1e-8:
                violations += 1
                if violations <= 3:
                    print(f"  n={n} VIOLATION: R_r={R_r:.8f}, R_p+R_q={R_p+R_q:.8f}, diff={R_r-R_p-R_q:.2e}")
        except:
            pass

    print(f"n={n} (centered): R_n superadditive? {violations} violations / {valid} trials "
          f"({'YES' if violations == 0 else 'NO'})")
