"""
DEFINITIVE numerical investigation using the CORRECT MSS formula.

From the proof document (fisher_subordination_proof.md, line 18):
c_k = sum_{i+j=k} [(n-i)!(n-j)! / (n!(n-k)!)] * a_i * b_j

where p(x) = sum_{k=0}^n a_k x^{n-k} (monic, so a_0 = 1).
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
    """Compute coefficient vector a_0,...,a_n where p(x) = sum a_k x^{n-k}."""
    n = len(roots)
    # a_k = (-1)^k * e_k
    ek = [elem_sym_poly(roots, k) for k in range(n+1)]
    return [(-1)**k * ek[k] for k in range(n+1)]


def boxplus_mss(roots_p, roots_q):
    """Compute p boxplus_n q using the EXACT MSS formula.

    c_k = sum_{i+j=k} [(n-i)!(n-j)! / (n!(n-k)!)] * a_i * b_j
    """
    n = len(roots_p)
    assert len(roots_q) == n

    a = poly_coeffs_from_roots(roots_p)
    b = poly_coeffs_from_roots(roots_q)

    c = np.zeros(n+1)
    for k in range(n+1):
        for i in range(k+1):
            j = k - i
            if i <= n and j <= n:
                coeff = factorial(n-i) * factorial(n-j) / (factorial(n) * factorial(n-k))
                c[k] += coeff * a[i] * b[j]

    # Build polynomial: p(x) = sum c_k x^{n-k} = c_0 x^n + c_1 x^{n-1} + ...
    poly_coeffs = np.array(c)
    roots_r = np.sort(np.real(np.roots(poly_coeffs)))
    return roots_r, c


def H_values(roots):
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H


def Phi_n(roots):
    H = H_values(roots)
    return np.sum(H**2)


def subordination(roots_p, roots_r):
    """Compute omega_1 properties at roots of r."""
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
        rppp = np.polyval(r_tprime, z0) if len(r_tprime) > 0 else 0

        pp = np.polyval(p_prime, w0)
        ppp = np.polyval(p_dprime, w0)
        pppp = np.polyval(p_tprime, w0) if len(p_tprime) > 0 else 0

        # F(z,w) = r'(z)*p(w) - p'(w)*r(z)
        # At root point: p(w0)=0, r(z0)=0
        F_z = -pp * rp  # + rpp*p(w0) = 0
        F_w = rp * pp   # - ppp*r(z0) = 0

        omega1_prime[k] = -F_z / F_w

        F_zz = -pp * rpp  # + rppp*p(w0) = 0
        F_zw = rpp * pp - ppp * rp
        F_ww = rp * ppp  # - pppp*r(z0) = 0

        w1 = omega1_prime[k]
        omega1_dprime[k] = -(F_zz + 2*F_zw*w1 + F_ww*w1**2) / F_w

    return omega1_at_roots, omega1_prime, omega1_dprime, sigma


# =============================================================
# Verify the MSS formula
# =============================================================

print("="*70)
print("VERIFYING MSS FORMULA")
print("="*70)

# n=2 check: should give exact equality
a, b, c, d = -1.0, 1.0, -2.0, 2.0
roots_r, c_vec = boxplus_mss(np.array([a, b]), np.array([c, d]))
print(f"n=2 test: p=({a},{b}), q=({c},{d})")
print(f"  MSS roots:   {roots_r}")
print(f"  Exact roots: [{-(a-b)**2/4 + (c-d)**2/4:.6f}... wait]")
# exact n=2 formula: mean = (a+b+c+d)/2, spread = sqrt(s^2+t^2)/2
s, t = b-a, d-c
mean_r = (a+b+c+d)/2
spread_r = np.sqrt(s**2 + t**2)/2
exact_r = np.array([mean_r - spread_r, mean_r + spread_r])
print(f"  Exact roots: {exact_r}")
print(f"  Match: {np.allclose(roots_r, exact_r)}")

# Check the coefficient formula for n=2
# a = [1, -(a+b), ab], b = [1, -(c+d), cd]
# c_0 = (2!*2!)/(2!*2!) * 1 * 1 = 1
# c_1 = (1!*2!)/(2!*1!) * (-(a+b)) * 1 + (2!*1!)/(2!*1!) * 1 * (-(c+d))
#      = (1/2) * (-(a+b)) + 1 * (-(c+d))  ... wait this can't be right
# Let me check: c_1 = sum_{i+j=1} [(n-i)!(n-j)!/(n!(n-1)!)] * a_i * b_j
#   i=0,j=1: [(2)!(1)!/(2!*1!)] * 1 * (-(c+d)) = 1 * (-(c+d))
#   i=1,j=0: [(1)!(2)!/(2!*1!)] * (-(a+b)) * 1 = 1 * (-(a+b))
# c_1 = -(a+b+c+d)  OK that's the mean, correct.

# c_2 = sum_{i+j=2} [(n-i)!(n-j)!/(n!(n-2)!)] * a_i * b_j (n=2, n-2=0)
#   i=0,j=2: [(2)!(0)!/(2!*0!)] * 1 * cd = 1 * cd
#   i=1,j=1: [(1)!(1)!/(2!*0!)] * (-(a+b)) * (-(c+d)) = (1/2) * (a+b)(c+d)
#   i=2,j=0: [(0)!(2)!/(2!*0!)] * ab * 1 = 1 * ab
# c_2 = cd + (a+b)(c+d)/2 + ab

# For p=(x-(-1))(x-1) = x^2-1: a=[1, 0, -1]
# q=(x-(-2))(x-2) = x^2-4: b=[1, 0, -4]
# c_2 = -4 + 0 + (-1) = -5
# r(x) = x^2 - 0*x + (-5) = x^2 - 5. Roots: +/-sqrt(5) ✓

# Another test
a_coeff = poly_coeffs_from_roots(np.array([0.0, 1.0]))
b_coeff = poly_coeffs_from_roots(np.array([0.0, 3.0]))
print(f"\np=x(x-1), a={a_coeff}")
print(f"q=x(x-3), b={b_coeff}")
roots_r2, c2 = boxplus_mss(np.array([0.0, 1.0]), np.array([0.0, 3.0]))
exact_r2 = np.array([0.5+1.5-np.sqrt(1+9)/2, 0.5+1.5+np.sqrt(1+9)/2])
print(f"MSS roots: {roots_r2}")
print(f"Exact:     {exact_r2}")
print(f"Match: {np.allclose(roots_r2, exact_r2)}")


# n=3 test against MC
print("\nn=3 test against Monte Carlo:")
roots_p = np.array([-2.0, 0.0, 2.0])
roots_q = np.array([-3.0, 0.0, 3.0])
roots_mss, _ = boxplus_mss(roots_p, roots_q)

# MC verification
def boxplus_mc(roots_p, roots_q, samples=500000):
    n = len(roots_p)
    A = np.diag(np.array(roots_p, dtype=float))
    B = np.diag(np.array(roots_q, dtype=float))
    sum_ek = np.zeros(n+1)
    for _ in range(samples):
        Z = (np.random.randn(n, n) + 1j * np.random.randn(n, n)) / np.sqrt(2)
        Q, R_mat = np.linalg.qr(Z)
        dd = np.diagonal(R_mat)
        ph = dd / np.abs(dd)
        U = Q * ph[np.newaxis, :]
        M = A + U @ B @ U.conj().T
        eigs = np.sort(np.linalg.eigvalsh(M))
        for k in range(n+1):
            sum_ek[k] += elem_sym_poly(eigs, k)
    avg_ek = sum_ek / samples
    coeffs = [(-1)**k * avg_ek[k] for k in range(n+1)]
    return np.sort(np.real(np.roots(coeffs))), avg_ek

np.random.seed(42)
roots_mc, ek_mc = boxplus_mc(roots_p, roots_q, 500000)
print(f"MSS formula: {roots_mss}")
print(f"MC (500k):   {roots_mc}")
print(f"Error:       {np.max(np.abs(roots_mss - roots_mc)):.6f}")


# =============================================================
# COMPREHENSIVE INEQUALITY TEST
# =============================================================

print("\n\n" + "="*70)
print("COMPREHENSIVE INEQUALITY TEST WITH CORRECT MSS FORMULA")
print("="*70)

np.random.seed(42)

results = []
failures = 0
total = 0

for trial in range(500):
    n = np.random.choice([2, 3, 4, 5])

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

        # Check real and distinct
        raw_roots = np.roots(c_vec)
        if np.any(np.abs(np.imag(raw_roots)) > 0.01):
            continue
        roots_r = np.sort(np.real(raw_roots))
        gaps = np.diff(roots_r)
        if np.any(gaps < 0.01):
            continue

        Phi_p_val = Phi_n(roots_p)
        Phi_q_val = Phi_n(roots_q)
        Phi_r_val = Phi_n(roots_r)

        ineq = 1/Phi_r_val - 1/Phi_p_val - 1/Phi_q_val
        total += 1

        if ineq < -1e-8:
            failures += 1
            if failures <= 5:
                print(f"  FAIL trial={trial}, n={n}: ineq={ineq:.10f}")
                print(f"    Phi_p={Phi_p_val:.6f}, Phi_q={Phi_q_val:.6f}, Phi_r={Phi_r_val:.6f}")
                print(f"    p={roots_p}, q={roots_q}, r={roots_r}")

        results.append({
            'n': n, 'ineq': ineq,
            'roots_p': roots_p, 'roots_q': roots_q, 'roots_r': roots_r,
            'Phi_p': Phi_p_val, 'Phi_q': Phi_q_val, 'Phi_r': Phi_r_val,
        })
    except Exception as e:
        pass

print(f"\nTotal valid: {total}, Failures: {failures}")
for n_val in [2, 3, 4, 5]:
    n_res = [r for r in results if r['n'] == n_val]
    if n_res:
        ineqs = [r['ineq'] for r in n_res]
        print(f"  n={n_val}: {len(n_res)} trials, min(ineq)={min(ineqs):.10f}, failures={sum(1 for x in ineqs if x < -1e-8)}")


# =============================================================
# STRUCTURAL ANALYSIS (NODE 1.7): alpha, beta, <h,alpha>
# =============================================================

print("\n\n" + "="*70)
print("STRUCTURAL ANALYSIS (NODE 1.7)")
print("="*70)

np.random.seed(42)

h_alpha_signs = []
h_beta_signs = []
AB_h4_signs = []
struct_results = []

for trial in range(300):
    n = np.random.choice([2, 3, 4, 5])

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

        # Subordination
        _, omega1_p, omega1_pp, sigma = subordination(roots_p, roots_r)
        _, omega2_p, omega2_pp, tau = subordination(roots_q, roots_r)

        h = H_values(roots_r)
        alpha = omega1_pp / 2
        beta = omega2_pp / 2
        u = H_values(roots_p)[sigma]
        v = H_values(roots_q)[tau]

        # Check chain rule
        cr_err = max(np.max(np.abs(u - (h + alpha))), np.max(np.abs(v - (h + beta))))
        if cr_err > 0.01:
            continue

        h_alpha = np.dot(h, alpha)
        h_beta = np.dot(h, beta)
        norm_h2 = np.dot(h, h)
        A = np.dot(u,u) - norm_h2
        B = np.dot(v,v) - norm_h2
        AB = A * B
        h4 = norm_h2**2

        h_alpha_signs.append(h_alpha >= -1e-8)
        h_beta_signs.append(h_beta >= -1e-8)
        AB_h4_signs.append(AB >= h4 - 1e-8)

        struct_results.append({
            'n': n, 'h_alpha': h_alpha, 'h_beta': h_beta,
            'A': A, 'B': B, 'AB': AB, 'h4': h4,
            'gap': AB - h4,
        })

    except:
        pass

print(f"Total valid: {len(h_alpha_signs)}")
print(f"<h,alpha> >= 0: {sum(h_alpha_signs)}/{len(h_alpha_signs)} ({100*sum(h_alpha_signs)/len(h_alpha_signs):.1f}%)")
print(f"<h,beta>  >= 0: {sum(h_beta_signs)}/{len(h_beta_signs)} ({100*sum(h_beta_signs)/len(h_beta_signs):.1f}%)")
print(f"AB >= h^4:      {sum(AB_h4_signs)}/{len(AB_h4_signs)} ({100*sum(AB_h4_signs)/len(AB_h4_signs):.1f}%)")

# By n
for n_val in [2, 3, 4, 5]:
    sr = [r for r in struct_results if r['n'] == n_val]
    if sr:
        print(f"\n  n={n_val}: {len(sr)} trials")
        print(f"    min(<h,alpha>) = {min(r['h_alpha'] for r in sr):.8f}")
        print(f"    min(<h,beta>)  = {min(r['h_beta'] for r in sr):.8f}")
        print(f"    min(AB-h4)     = {min(r['gap'] for r in sr):.10f}")
        print(f"    <h,alpha> < 0: {sum(1 for r in sr if r['h_alpha'] < -1e-8)}/{len(sr)}")
        print(f"    AB < h4:       {sum(1 for r in sr if r['gap'] < -1e-8)}/{len(sr)}")

# Show counterexamples to <h,alpha> >= 0
neg_examples = [r for r in struct_results if r['h_alpha'] < -1e-8]
if neg_examples:
    print(f"\n<h,alpha> < 0 examples (first 3):")
    for r in neg_examples[:3]:
        print(f"  n={r['n']}: <h,a>={r['h_alpha']:.6f}, <h,b>={r['h_beta']:.6f}, A={r['A']:.6f}, B={r['B']:.6f}, gap={r['gap']:.10f}")

# =============================================================
# KEY QUESTION: Does A > 0 always? (i.e., Phi_p > Phi_r)
# =============================================================
print("\n\n" + "="*70)
print("KEY QUESTION: Is A > 0 always? (Phi_p > Phi_r)")
print("="*70)
A_positive = sum(1 for r in struct_results if r['A'] > -1e-8)
B_positive = sum(1 for r in struct_results if r['B'] > -1e-8)
print(f"A > 0: {A_positive}/{len(struct_results)}")
print(f"B > 0: {B_positive}/{len(struct_results)}")
A_neg_examples = [r for r in struct_results if r['A'] < -1e-8]
if A_neg_examples:
    print(f"A < 0 examples (first 3):")
    for r in A_neg_examples[:3]:
        print(f"  n={r['n']}: A={r['A']:.8f}, B={r['B']:.8f}")
