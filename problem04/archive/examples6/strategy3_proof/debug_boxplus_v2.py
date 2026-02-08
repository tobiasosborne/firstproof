"""
debug_boxplus_v2.py — Identify which convolution formula is correct.

The Monte Carlo simulation (expected characteristic polynomial of A + UBU*)
gives DIFFERENT results from my MSS coefficient formula for non-centered n=2.

HYPOTHESIS: The MSS formula I'm using might be wrong, or there's a sign/
normalization error.

Let me test multiple formulations and compare with Monte Carlo.
"""

import numpy as np
from math import factorial, comb
from itertools import combinations

np.random.seed(42)

def mc_boxplus(roots_p, roots_q, N=500000):
    """Monte Carlo: E[char poly of A + UBU*] with Haar unitary U."""
    n = len(roots_p)
    A = np.diag(roots_p)
    B = np.diag(roots_q)

    poly_sum = np.zeros(n+1)
    for _ in range(N):
        Z = (np.random.randn(n, n) + 1j * np.random.randn(n, n)) / np.sqrt(2)
        Q, R = np.linalg.qr(Z)
        dp = np.diag(R)
        dp = dp / np.abs(dp)
        U = Q @ np.diag(dp)
        M = A + U @ B @ U.conj().T
        poly_sum += np.real(np.poly(M))

    poly_avg = poly_sum / N
    roots_r = np.sort(np.real(np.roots(poly_avg)))
    return roots_r, poly_avg

def boxplus_v1(roots_p, roots_q):
    """MSS formula: c_k = sum w_{i,j} a_i b_j with w = (n-i)!(n-j)!/(n!(n-k)!)."""
    n = len(roots_p)
    poly_p = np.poly(roots_p)
    poly_q = np.poly(roots_q)

    a = np.zeros(n+1)
    b = np.zeros(n+1)
    for k in range(n+1):
        a[k] = (-1)**k * poly_p[k] / comb(n, k)
        b[k] = (-1)**k * poly_q[k] / comb(n, k)

    c = np.zeros(n+1)
    for k in range(n+1):
        for i in range(k+1):
            j = k - i
            w = factorial(n-i) * factorial(n-j) / (factorial(n) * factorial(n-k))
            c[k] += w * a[i] * b[j]

    poly_r = np.zeros(n+1)
    for k in range(n+1):
        poly_r[k] = (-1)**k * comb(n, k) * c[k]

    return np.sort(np.real(np.roots(poly_r))), poly_r

def boxplus_v2(roots_p, roots_q):
    """
    Simple additive free cumulants:
    r_k = e_k / C(n,k), and r_k(p boxplus q) = r_k(p) + r_k(q).
    """
    n = len(roots_p)
    poly_p = np.poly(roots_p)
    poly_q = np.poly(roots_q)

    poly_r = np.zeros(n+1)
    for k in range(n+1):
        # e_k(p) = (-1)^k * poly_p[k]
        # r_k(p) = e_k(p) / C(n,k)
        # r_k(r) = r_k(p) + r_k(q)
        # e_k(r) = C(n,k) * (r_k(p) + r_k(q))
        e_k_p = (-1)**k * poly_p[k]
        e_k_q = (-1)**k * poly_q[k]
        r_k_p = e_k_p / comb(n, k)
        r_k_q = e_k_q / comb(n, k)
        e_k_r = comb(n, k) * (r_k_p + r_k_q)
        poly_r[k] = (-1)**k * e_k_r

    return np.sort(np.real(np.roots(poly_r))), poly_r

def boxplus_v3(roots_p, roots_q):
    """
    MSS convolution using the CORRECT formula from Marcus-Spielman-Srivastava
    'Interlacing Families II' Theorem 4.4.

    The finite free convolution is defined via the characteristic polynomial:
    E[det(xI - A - UBU*)] over Haar-random U.

    For degree n monic polynomials p, q:
    (p boxplus q)(x) = sum_{k=0}^n sum_{j=0}^{n-k} (-1)^{n-k}
                       * C(n-k, j) / C(n, j) * p^{(j)}(x) * q_j(0)
    ... this is getting complicated. Let me use the explicit coefficient formula.

    From Marcus 2021 "Polynomial convolutions":
    Write p(x) = sum_{k=0}^n (-1)^k e_k x^{n-k} (elementary symmetric coefficients)
    Write q(x) = sum_{k=0}^n (-1)^k f_k x^{n-k}

    Then (p boxplus_n q)(x) = sum_{k=0}^n (-1)^k g_k x^{n-k}
    where g_k = sum_{i+j=k} C(n-j, i) / C(n, i) * e_i * f_j
            (equivalently: C(n-i, j) / C(n, j) * e_i * f_j by symmetry)
    """
    n = len(roots_p)
    poly_p = np.poly(roots_p)
    poly_q = np.poly(roots_q)

    e = np.zeros(n+1)
    f = np.zeros(n+1)
    for k in range(n+1):
        e[k] = (-1)**k * poly_p[k]
        f[k] = (-1)**k * poly_q[k]

    g = np.zeros(n+1)
    for k in range(n+1):
        for i in range(k+1):
            j = k - i
            if comb(n, i) > 0:
                w = comb(n-j, i) / comb(n, i)
                g[k] += w * e[i] * f[j]

    poly_r = np.zeros(n+1)
    for k in range(n+1):
        poly_r[k] = (-1)**k * g[k]

    return np.sort(np.real(np.roots(poly_r))), poly_r


# ============================================================
# TEST: Compare all formulations with Monte Carlo
# ============================================================

test_cases = [
    ("centered n=2", np.array([-1.0, 1.0]), np.array([-1.0, 1.0])),
    ("non-centered n=2", np.array([0.0, 1.0]), np.array([-3.0, 0.0])),
    ("non-centered n=2 v2", np.array([-1.0, 3.0]), np.array([-4.0, 0.0])),
    ("centered n=3", np.array([-1.0, 0.0, 1.0]), np.array([-2.0, 0.0, 2.0])),
    ("non-centered n=3", np.array([0.0, 1.0, 3.0]), np.array([-1.0, 0.0, 2.0])),
]

for name, rp, rq in test_cases:
    print(f"\n{'='*60}")
    print(f"Test: {name}")
    print(f"  p roots: {rp}")
    print(f"  q roots: {rq}")
    print(f"{'='*60}")

    mc_roots, mc_poly = mc_boxplus(rp, rq, N=300000)
    v1_roots, v1_poly = boxplus_v1(rp, rq)
    v2_roots, v2_poly = boxplus_v2(rp, rq)
    v3_roots, v3_poly = boxplus_v3(rp, rq)

    print(f"  MC roots:    {mc_roots}")
    print(f"  v1 (MSS) :   {v1_roots}")
    print(f"  v2 (simple): {v2_roots}")
    print(f"  v3 (Marcus): {v3_roots}")

    print(f"\n  MC poly:     {np.array2string(mc_poly, precision=4)}")
    print(f"  v1 (MSS) :   {np.array2string(v1_poly, precision=4)}")
    print(f"  v2 (simple): {np.array2string(v2_poly, precision=4)}")
    print(f"  v3 (Marcus): {np.array2string(v3_poly, precision=4)}")

    # Which matches MC best?
    err_v1 = np.max(np.abs(v1_poly - mc_poly))
    err_v2 = np.max(np.abs(v2_poly - mc_poly))
    err_v3 = np.max(np.abs(v3_poly - mc_poly))
    print(f"\n  Max poly error vs MC:")
    print(f"    v1 (MSS):    {err_v1:.6f}")
    print(f"    v2 (simple): {err_v2:.6f}")
    print(f"    v3 (Marcus): {err_v3:.6f}")

    best = ["v1", "v2", "v3"][[err_v1, err_v2, err_v3].index(min(err_v1, err_v2, err_v3))]
    print(f"  Best match: {best}")
