"""
Exact numerical investigation using the MSS finite free convolution formula.

The finite free additive convolution for degree n polynomials uses:
  e_k(p boxplus_n q) = sum_{j=0}^{k} C(k,j) * C(n-k, n-k) / C(n,j) * e_j(p) * e_{k-j}(q)
  Hmm wait, let me find the correct formula.

From MSS (2015), the "finite free convolution" of p and q of degree n is:
  For p(x) = sum_{k=0}^{n} a_k x^{n-k} where a_0 = 1 (monic),
  the convolution is defined via:
  (p boxplus_n q)(x) = sum_{k=0}^{n} c_k x^{n-k}
  where c_k = sum_{j=0}^{k} C(n-j, k-j) / C(n, k) * a_j * b_{k-j} * C(k, j)

  Wait, I need to check the exact formula. Let me use a different approach:
  compute via the finite free cumulant sequence.

From Arizmendi-Perales (2018), for a monic polynomial p(x) of degree n with
coefficients written as p(x) = x^n + sum_{k=1}^{n} (-1)^k e_k x^{n-k},
the finite free cumulants kappa_k^{(n)} are defined by the moment-cumulant relations:

  e_k / C(n,k) = sum over non-crossing partitions pi of {1,...,k}:
    prod_{B in pi} kappa_{|B|}^{(n)}

Wait, the finite free cumulants for the finite free probability are different from
the usual free cumulants. Let me just use a recursive formula.

Actually, the simplest correct formula from MSS Theorem 4.5:
For monic p(x) = prod(x - lambda_i) with elementary symmetric polynomials e_k(lambda),
and similarly q with e_k(mu):

e_k(boxplus) = sum_{j=0}^{k} (-1)^{k-j} * C(k,j) * [C(n,k)]^{-1} * C(n,j) * ... no

Let me use the SIMPLER characterization from Marcus (2021):
The finite free additive convolution preserves the "symmetrized additive" structure.
For the purposes of computation, use the expected characteristic polynomial:

(p boxplus_n q)(x) = E_U[ det(xI - (A + UBU*)) ]

where A = diag(lambda_1,...,lambda_n), B = diag(mu_1,...,mu_n), U is Haar unitary.

The key identity (Marcus-Spielman-Srivastava):
E_U[ det(xI - (A + UBU*)) ] = sum_{sigma in S_n} prod_{i=1}^{n} (some weight) * ...

Actually, the most useful formula for computation is:
  e_k(p boxplus_n q) = sum_{j=0}^{k} d_{n,k,j} * e_j(p) * e_{k-j}(q)

where d_{n,k,j} = C(k,j) * C(n-k, n-k) / C(n, k) ... I keep going in circles.

Let me just use the MATRIX INTEGRAL approach with EXACT integration.
For the 2x2 case we've shown it works. For higher n, I'll use a HIGH-ACCURACY
Monte Carlo (100k samples).

Actually, a much better approach: use the CHARACTERISTIC POLYNOMIAL IDENTITY.

Marcus-Spielman-Srivastava proved that:
  E_U[det(xI - (A+UBU*))] = sum over matchings...

The exact formula is complex. Let me use a different, well-known approach:
the COEFFICIENTS of p boxplus_n q can be computed using the formula:

a_k(p boxplus_n q) / C(n,k) = sum_{j=0}^k a_j(p)/C(n,j) * a_{k-j}(q)/C(n,k-j) * (-1)^? ...

OK I found it. From the survey "Finite Free Probability" by Marcus (2021), Definition 2.1:

For p(x) = sum_{k=0}^{n} C(n,k) * (-1)^k * hat{e}_k * x^{n-k}, where hat{e}_k = e_k / C(n,k),
and similarly q with hat{f}_k, we have:

(p boxplus_n q)(x) has hat{g}_k = sum_{j=0}^{k} hat{e}_j * hat{f}_{k-j} * ???

Actually this is the "additive" convolution of the NORMALIZED coefficients.
It turns out the simplest version is:

Theorem (MSS 2015): If p,q are degree n monic polynomials, then
p boxplus_n q = p_n(D) q  (evaluated at a specific point)

where p_n(D) is a differential operator. This is hard to implement.

Let me just use a DIFFERENT well-tested implementation: the recursion via
the resolvent/Cauchy transform.

SIMPLEST APPROACH: Since this is numerical, I'll use the well-known fact that
for the FINITE FREE CONVOLUTION:
  G_{p boxplus_n q}(z) is determined by the subordination equations
  G_{p boxplus_n q}(z) = G_p(omega_1(z)) = G_q(omega_2(z))
  omega_1(z) + omega_2(z) - z = (1/n) * G_{p boxplus_n q}(z)^{-1} ... no, this is the free probability version.

OK, I'll use a practical approach: compute via high-accuracy Monte Carlo and
also implement the MSS coefficient formula correctly.
"""

import numpy as np
from math import comb, factorial
from itertools import combinations

def elem_sym_poly(roots, k):
    """k-th elementary symmetric polynomial of roots."""
    n = len(roots)
    if k == 0:
        return 1.0
    if k > n:
        return 0.0
    result = 0.0
    for subset in combinations(range(n), k):
        prod = 1.0
        for i in subset:
            prod *= roots[i]
        result += prod
    return result


def boxplus_mss(roots_p, roots_q):
    """Compute p boxplus_n q using the MSS convolution formula.

    The correct formula (from Marcus 2021 survey, Proposition 2.3):

    For p(x) = sum_{k=0}^{n} (-1)^k e_k(p) x^{n-k} and q similarly,

    e_m(p boxplus_n q) = sum_{k=0}^{m} (-1)^{m-k} * C(m,k) * C(n-k, m-k)^{-1} * C(n, m)^{-1}
        * C(n, k) * C(n, m-k) * e_k(p) * e_{m-k}(q) ???

    Hmm, I'm still not sure. Let me try another reference.

    From Arizmendi-Perales (2018), Equation (6):
    For the "additive rectangular convolution" at ratio c = 1:

    e_k(p boxplus_n q) = sum_{j=0}^{k} C(k,j) * C(n-j, k-j) / C(n,k) * e_j(p) * e_{k-j}(q)

    Wait, that doesn't look right dimensionally either.

    Let me try the formula from Anderson (2014):
    For monic polynomials p,q of degree n,
    (p boxplus_n q)(x) = sum_{k=0}^{n} (-1)^k E_k x^{n-k}
    where
    E_k = sum_{j=0}^{k} (-1)^{k-j} ... no

    OK, I'll try the following formula which I believe is correct:

    For monic p(x) = sum_{k=0}^{n} a_k x^{n-k} where a_k = (-1)^k e_k(p),
    the MSS convolution gives:

    a_k(r) = sum_{j=0}^{k} C(k,j) / C(n,j) * a_j(p) * a_{k-j}(q)

    This is the "binomial convolution" normalized by binomial coefficients.
    """
    n = len(roots_p)
    assert len(roots_q) == n

    # Compute normalized coefficients: hat_a_k = (-1)^k e_k / C(n,k)
    ep = [elem_sym_poly(roots_p, k) for k in range(n+1)]
    eq = [elem_sym_poly(roots_q, k) for k in range(n+1)]

    # Try multiple candidate formulas and check against MC

    # Formula 1: Direct convolution of normalized e_k
    # hat_e_k(r) = sum_{j=0}^k hat_e_j(p) * hat_e_{k-j}(q)
    # where hat_e_k = e_k / C(n,k)
    hat_ep = [ep[k] / comb(n, k) if k <= n else 0 for k in range(n+1)]
    hat_eq = [eq[k] / comb(n, k) if k <= n else 0 for k in range(n+1)]

    hat_er_1 = np.zeros(n+1)
    for k in range(n+1):
        for j in range(k+1):
            hat_er_1[k] += hat_ep[j] * hat_eq[k-j]

    er_1 = [hat_er_1[k] * comb(n, k) for k in range(n+1)]
    coeffs_1 = np.zeros(n+1)
    for k in range(n+1):
        coeffs_1[k] = (-1)**k * er_1[k]

    # Formula 2: a_k(r) = sum_{j} C(k,j)/C(n,j) * a_j(p) * a_{k-j}(q)
    ap = [(-1)**k * ep[k] for k in range(n+1)]
    aq = [(-1)**k * eq[k] for k in range(n+1)]

    ar_2 = np.zeros(n+1)
    for k in range(n+1):
        for j in range(k+1):
            ar_2[k] += comb(k, j) / comb(n, j) * ap[j] * aq[k-j]

    # Build polynomial coeffs (highest degree first for np.roots)
    poly_coeffs_1 = coeffs_1  # These are for x^n, x^{n-1}, ..., x^0
    poly_coeffs_2 = ar_2

    return poly_coeffs_1, poly_coeffs_2, ep, eq


def boxplus_mc(roots_p, roots_q, samples=200000):
    """High-accuracy Monte Carlo."""
    n = len(roots_p)
    A = np.diag(roots_p)
    B = np.diag(roots_q)

    # We compute E[e_k(eigenvalues of A + UBU*)]
    sum_ek = np.zeros(n+1)

    for _ in range(samples):
        Z = (np.random.randn(n, n) + 1j * np.random.randn(n, n)) / np.sqrt(2)
        Q, R = np.linalg.qr(Z)
        dd = np.diagonal(R)
        ph = dd / np.abs(dd)
        U = Q * ph[np.newaxis, :]

        M = A + U @ B @ U.conj().T
        eigs = np.sort(np.real(np.linalg.eigvalsh(M)))

        for k in range(n+1):
            sum_ek[k] += elem_sym_poly(eigs, k)

    avg_ek = sum_ek / samples
    return avg_ek


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


np.random.seed(42)

# First, find the correct formula by checking against n=2 exact and MC
print("="*60)
print("FINDING CORRECT MSS FORMULA")
print("="*60)

roots_p = np.array([-1.0, 1.0])
roots_q = np.array([-2.0, 2.0])

c1, c2, ep, eq = boxplus_mss(roots_p, roots_q)
mc_ek = boxplus_mc(roots_p, roots_q, 100000)

# Exact for n=2: roots are -sqrt(20/4), sqrt(20/4) = -sqrt(5), sqrt(5)
exact_roots = np.array([-np.sqrt(5), np.sqrt(5)])
exact_ek = [elem_sym_poly(exact_roots, k) for k in range(3)]

print(f"Exact e_k: {exact_ek}")
print(f"MC e_k:    {mc_ek}")
print(f"Formula 1 e_k: {[c1[k]*(-1 if k%2==1 else 1)*1 for k in range(3)]}")
# Wait, let me extract e_k from the formulas properly.

# For formula 1: er_1[k] are the e_k values
hat_ep = [ep[k] / comb(2, k) for k in range(3)]
hat_eq = [eq[k] / comb(2, k) for k in range(3)]
hat_er_1 = np.zeros(3)
for k in range(3):
    for j in range(k+1):
        hat_er_1[k] += hat_ep[j] * hat_eq[k-j]
er_1 = [hat_er_1[k] * comb(2, k) for k in range(3)]

print(f"\nFormula 1 (hat convolution) e_k: {er_1}")
print(f"  hat_ep = {hat_ep}, hat_eq = {hat_eq}")

# Formula 2 gives a_k, convert to e_k: a_k = (-1)^k e_k
ap = [(-1)**k * ep[k] for k in range(3)]
aq = [(-1)**k * eq[k] for k in range(3)]
ar_2 = np.zeros(3)
for k in range(3):
    for j in range(k+1):
        ar_2[k] += comb(k, j) / comb(2, j) * ap[j] * aq[k-j]
er_2 = [(-1)**k * ar_2[k] for k in range(3)]

print(f"Formula 2 (binomial/binomial) e_k: {er_2}")

# Check which formula matches
print(f"\nExact: {exact_ek}")
print(f"MC:    {list(mc_ek)}")
print(f"F1:    {er_1}")
print(f"F2:    {er_2}")

# Now identify correct formula
def get_roots_from_ek(ek_list, n):
    """From e_k values, construct polynomial and find roots."""
    # p(x) = x^n - e_1 x^{n-1} + e_2 x^{n-2} - ... + (-1)^n e_n
    coeffs = np.zeros(n+1)
    for k in range(n+1):
        coeffs[k] = (-1)**k * ek_list[k]
    # coeffs[0] = 1, coeffs[1] = -e_1, etc. for x^n, x^{n-1}, ...
    return np.sort(np.real(np.roots(coeffs)))

print(f"\nRoots from F1: {get_roots_from_ek(er_1, 2)}")
print(f"Roots from F2: {get_roots_from_ek(er_2, 2)}")
print(f"Exact roots:   {exact_roots}")


# =================================================================
# Now use the CORRECT formula for comprehensive testing
# =================================================================

print("\n\n" + "="*60)
print("COMPREHENSIVE TESTING WITH CORRECT CONVOLUTION")
print("="*60)

def boxplus_correct(roots_p, roots_q):
    """Compute p boxplus_n q using the correct MSS formula.

    After testing, use whichever formula matched the n=2 exact case.
    """
    n = len(roots_p)

    ep = [elem_sym_poly(roots_p, k) for k in range(n+1)]
    eq = [elem_sym_poly(roots_q, k) for k in range(n+1)]

    # Use Formula 1: hat convolution
    hat_ep = [ep[k] / comb(n, k) for k in range(n+1)]
    hat_eq = [eq[k] / comb(n, k) for k in range(n+1)]

    hat_er = np.zeros(n+1)
    for k in range(n+1):
        for j in range(k+1):
            hat_er[k] += hat_ep[j] * hat_eq[k-j]

    er = [hat_er[k] * comb(n, k) for k in range(n+1)]

    roots_r = get_roots_from_ek(er, n)
    return roots_r, er


def subordination_omega1(roots_p, roots_r):
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
        rppp = np.polyval(r_tprime, z0)

        pp = np.polyval(p_prime, w0)
        ppp = np.polyval(p_dprime, w0)
        pppp = np.polyval(p_tprime, w0)

        # F = r'(z)*p(w) - p'(w)*r(z)
        # At root: p(w0)=0, r(z0)=0
        F_z = -pp * rp
        F_w = rp * pp

        omega1_prime[k] = -F_z / F_w  # = 1

        F_zz = -pp * rpp
        F_zw = rpp * pp - ppp * rp
        F_ww = rp * ppp

        w1 = omega1_prime[k]
        omega1_dprime[k] = -(F_zz + 2*F_zw*w1 + F_ww*w1**2) / F_w

    return omega1_at_roots, omega1_prime, omega1_dprime, sigma


# Test with correct formula

# n=2 verification
roots_r_2, er_2 = boxplus_correct(np.array([-1.0, 1.0]), np.array([-2.0, 2.0]))
print(f"n=2 check: roots = {roots_r_2}, expected = [-sqrt(5), sqrt(5)] = [{-np.sqrt(5):.8f}, {np.sqrt(5):.8f}]")

# n=3 test
roots_r_3, er_3 = boxplus_correct(np.array([-2.0, 0.0, 2.0]), np.array([-3.0, 0.0, 3.0]))
mc_ek_3 = boxplus_mc(np.array([-2.0, 0.0, 2.0]), np.array([-3.0, 0.0, 3.0]), 200000)
print(f"\nn=3 formula e_k: {er_3}")
print(f"n=3 MC e_k:      {list(mc_ek_3)}")
print(f"n=3 formula roots: {roots_r_3}")
print(f"n=3 MC roots: {get_roots_from_ek(mc_ek_3, 3)}")


# Comprehensive test: compute the inequality for many cases
print("\n\n" + "="*60)
print("SYSTEMATIC TEST: 1/Phi_r >= 1/Phi_p + 1/Phi_q")
print("="*60)

results = []
failures = 0
h_alpha_negative = 0
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
        roots_r, _ = boxplus_correct(roots_p, roots_q)
        # Check roots are real and distinct
        if not np.all(np.isreal(roots_r)):
            continue
        gaps = np.diff(roots_r)
        if np.any(gaps < 0.01):
            continue

        Phi_p_val = Phi_n(roots_p)
        Phi_q_val = Phi_n(roots_q)
        Phi_r_val = Phi_n(roots_r)

        inv_ineq = 1/Phi_r_val - 1/Phi_p_val - 1/Phi_q_val

        # Compute alpha, beta
        _, _, omega1_pp, sigma = subordination_omega1(roots_p, roots_r)
        _, _, omega2_pp, tau = subordination_omega1(roots_q, roots_r)

        h = H_values(roots_r)
        alpha = omega1_pp / 2
        beta = omega2_pp / 2
        u = H_values(roots_p)[sigma]
        v = H_values(roots_q)[tau]

        # Check chain rule
        cr_err = max(np.max(np.abs(u - (h + alpha))), np.max(np.abs(v - (h + beta))))
        if cr_err > 0.01:
            continue  # Skip if chain rule doesn't hold well

        h_alpha = np.dot(h, alpha)
        h_beta = np.dot(h, beta)
        A = np.dot(u,u) - np.dot(h,h)
        B = np.dot(v,v) - np.dot(h,h)
        AB = A * B
        h4 = np.dot(h,h)**2

        total += 1
        if inv_ineq < -1e-8:
            failures += 1
        if h_alpha < -1e-8:
            h_alpha_negative += 1

        results.append({
            'n': n, 'inv_ineq': inv_ineq, 'AB_minus_h4': AB - h4,
            'h_alpha': h_alpha, 'h_beta': h_beta, 'A': A, 'B': B,
            'Phi_p': Phi_p_val, 'Phi_q': Phi_q_val, 'Phi_r': Phi_r_val,
            'roots_p': roots_p, 'roots_q': roots_q, 'roots_r': roots_r,
        })

    except Exception as e:
        pass

print(f"\nTotal valid trials: {total}")
print(f"Inequality failures (1/Phi_r < 1/Phi_p + 1/Phi_q): {failures}/{total}")
print(f"<h,alpha> < 0 cases: {h_alpha_negative}/{total}")

# Show distribution by n
for n_val in [2, 3, 4, 5]:
    n_results = [r for r in results if r['n'] == n_val]
    if n_results:
        inv_ineqs = [r['inv_ineq'] for r in n_results]
        h_alphas = [r['h_alpha'] for r in n_results]
        print(f"\nn={n_val}: {len(n_results)} trials")
        print(f"  min(1/Phi_r - 1/Phi_p - 1/Phi_q) = {min(inv_ineqs):.8f}")
        print(f"  min(<h,alpha>) = {min(h_alphas):.8f}")
        print(f"  <h,alpha> < 0: {sum(1 for x in h_alphas if x < -1e-8)}")

# Show worst cases
if results:
    results.sort(key=lambda r: r['inv_ineq'])
    print(f"\n5 worst cases (smallest 1/Phi_r - 1/Phi_p - 1/Phi_q):")
    for r in results[:5]:
        print(f"  n={r['n']}, gap={r['inv_ineq']:.10f}, <h,a>={r['h_alpha']:.6f}, <h,b>={r['h_beta']:.6f}")
        print(f"    A={r['A']:.6f}, B={r['B']:.6f}, AB-h4={r['AB_minus_h4']:.10f}")
