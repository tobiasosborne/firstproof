"""
Numerical investigation of nodes 1.7 and 1.8.

We compute:
- Free additive convolution r = p boxplus_n q via the MSS coefficient formula
- Subordination functions omega_1, omega_2 from the implicit equation F(z,w) = 0
- The vectors h, u, v, alpha, beta
- Inner products <h,alpha>, <h,beta>
- A, B, ||h||^4, AB - ||h||^4
- Test whether <h,alpha> >= 0, <h,beta> >= 0
- Test whether AB >= ||h||^4
"""

import numpy as np
from numpy.polynomial import polynomial as P
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

def poly_from_roots(roots):
    """Create monic polynomial from roots: prod(x - r_i)."""
    # numpy poly convention: highest degree first
    coeffs = np.poly(roots)
    return coeffs

def poly_eval(coeffs, z):
    """Evaluate polynomial with numpy convention (highest first)."""
    return np.polyval(coeffs, z)

def poly_deriv(coeffs):
    """Derivative of polynomial (numpy convention)."""
    return np.polyder(coeffs)

def H_values(roots):
    """Compute H_p(lambda_i) = p''(lambda_i)/(2p'(lambda_i)) for each root.

    Equivalently, H_p(lambda_i) = sum_{j != i} 1/(lambda_i - lambda_j).
    Wait, that's not right. Let me be careful.

    p(x) = prod(x - lambda_j)
    p'(lambda_i) = prod_{j != i} (lambda_i - lambda_j)
    p''(lambda_i) = 2 * sum_{j != i} prod_{k != i, k != j} (lambda_i - lambda_k)

    So H_p(lambda_i) = p''(lambda_i)/(2p'(lambda_i)) = sum_{j != i} 1/(lambda_i - lambda_j)
    """
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H

def Phi_n(roots):
    """Phi_n(p) = sum_i H_p(lambda_i)^2."""
    H = H_values(roots)
    return np.sum(H**2)

def free_additive_convolution_roots(roots_p, roots_q, tol=1e-10):
    """Compute roots of r = p boxplus_n q using the subordination approach.

    We use the fact that if p(x) = prod(x - lambda_i) and q(x) = prod(x - mu_j),
    then r = p boxplus_n q has roots nu_k where nu_k = lambda_{sigma(k)} + mu_{tau(k)} - c
    for some constant c related to the means.

    Actually, let's use the MSS finite free convolution directly.
    p boxplus_n q is defined via the coefficient identity involving
    the finite free cumulants. For degree n:

    The MSS convolution uses the relation between coefficients and
    expected characteristic polynomials of random matrices, but for
    numerical testing, we can compute it via the R-transform approach.

    Simpler: use the Marcus-Spielman-Srivastava characteristic polynomial approach.
    For monic polynomials of degree n, p boxplus_n q can be computed as:
    the expected characteristic polynomial E[det(xI - (A + UBU*))]
    where A has eigenvalues = roots of p, B has eigenvalues = roots of q,
    U is Haar-random unitary. For numerical testing, we can average.

    Even simpler for testing: use the free convolution via finite free cumulants.
    """
    n = len(roots_p)
    assert len(roots_q) == n

    # Use the MSS formula: (p boxplus_n q)(x) is computed via mixed coefficients
    # e_k(p boxplus_n q) = sum_{j=0}^{k} C(n-j, k-j)/C(n, k-j) * (-1)^j * ...
    # This is complex. Let's use a Monte Carlo approach for numerical testing.

    # Monte Carlo: average char poly of A + U B U* over random unitaries
    A = np.diag(roots_p)
    B = np.diag(roots_q)

    num_samples = 5000
    sum_coeffs = np.zeros(n + 1)

    for _ in range(num_samples):
        # Random unitary from Haar measure
        Z = (np.random.randn(n, n) + 1j * np.random.randn(n, n)) / np.sqrt(2)
        Q, R = np.linalg.qr(Z)
        # Fix phases
        d = np.diagonal(R)
        ph = d / np.abs(d)
        U = Q * ph[np.newaxis, :]

        M = A + U @ B @ U.conj().T
        eigenvalues = np.linalg.eigvalsh(M.real + M.real.T) / 2  # Symmetrize
        # Actually, A + UBU* is Hermitian, so eigenvalues are real
        eigenvalues = np.sort(np.real(np.linalg.eigvals(M)))
        char_coeffs = np.poly(eigenvalues)
        sum_coeffs += np.real(char_coeffs)

    avg_coeffs = sum_coeffs / num_samples
    roots_r = np.sort(np.real(np.roots(avg_coeffs)))
    return roots_r


def free_additive_convolution_exact(roots_p, roots_q):
    """Exact finite free additive convolution using the MSS formula.

    For p(x) = sum_{k=0}^{n} (-1)^k e_k(p) x^{n-k} and similarly for q,
    the coefficients of r = p boxplus_n q are:

    e_k(r) = sum_{j=0}^{k} C(k,j) * C(n-j, k-j)^{-1} * C(n, k)^{-1} * ...

    Actually the correct formula from MSS is:
    e_k(p boxplus_n q) = sum_{j=0}^{k} (-1)^{k-j} * C(k,j) / C(n,j) * e_j(p) * e_{k-j}(q) * C(n,k) / ???

    Let me use the simpler characterization:
    The finite free cumulant kappa_k(p) is defined recursively from the coefficients.
    Then kappa_k(p boxplus_n q) = kappa_k(p) + kappa_k(q).

    For a monic degree-n polynomial with coefficients:
    p(x) = x^n - c_1 x^{n-1} + c_2 x^{n-2} - ... + (-1)^n c_n

    The finite free cumulants (Marcus-Spielman-Srivastava) are:
    kappa_1 = c_1 / n  (= mean of roots)

    Actually, I'll use a cleaner formulation. The MSS paper defines:
    For p(x) = x^n + a_{n-1} x^{n-1} + ... + a_0,
    the finite free convolution coefficients via:

    a_k(p boxplus_n q) / C(n,n-k) = sum_{j} a_j(p)/C(n,j) * a_{n-k-j}(q)/C(n,n-k-j)

    No, this isn't right either. Let me just use the explicit formula for small n.
    """
    n = len(roots_p)

    # Use elementary symmetric polynomials
    # For polynomial with roots r_1,...,r_n:
    # e_k = sum_{|S|=k} prod_{i in S} r_i  (note: these are the elementary symmetric polynomials of the roots)
    # p(x) = x^n - e_1 x^{n-1} + e_2 x^{n-2} - ... + (-1)^n e_n

    from math import comb

    def elem_sym(roots, k):
        """k-th elementary symmetric polynomial of roots."""
        if k == 0:
            return 1.0
        if k > len(roots):
            return 0.0
        result = 0.0
        for subset in combinations(range(len(roots)), k):
            result += np.prod([roots[i] for i in subset])
        return result

    ep = [elem_sym(roots_p, k) for k in range(n+1)]
    eq = [elem_sym(roots_q, k) for k in range(n+1)]

    # MSS finite free convolution formula:
    # e_k(r) = sum_{j=0}^{k} C(k,j) * C(n-k, n-k) / C(n, j) * ...
    #
    # The correct formula from Marcus-Spielman-Srivastava (2015) Theorem 4.5:
    # If p(x) = sum_{k=0}^{n} (-1)^k e_k x^{n-k}, then
    # e_k(p boxplus_n q) = sum_{j=0}^{k} C(n-j, k-j) / C(n, k) * C(k, j) * e_j(p) * e_{k-j}(q) ???
    #
    # Hmm, let me use a different reference. From the survey by Arizmendi-Perales:
    # The finite free additive convolution is characterized by:
    # kappa_k^{(n)}(p boxplus_n q) = kappa_k^{(n)}(p) + kappa_k^{(n)}(q)
    # where kappa_k^{(n)} are the finite free cumulants.
    #
    # Finite free cumulants are defined by:
    # e_k / C(n,k) = sum over partitions pi of [k]: prod_{B in pi} kappa_{|B|}^{(n)} / ...
    #
    # This is getting complicated. Let me just use the Monte Carlo approach.
    pass
    return None


def subordination_omega1(roots_p, roots_r):
    """Compute omega_1 at the roots of r using the subordination relation.

    G_r(z) = G_p(omega_1(z)). At z = nu_k (root of r), G_r has a pole.
    omega_1(nu_k) must be a root of p, say lambda_{sigma(k)}.

    We find sigma by matching: omega_1 is the unique analytic function with
    Im(omega_1(z)) > 0 when Im(z) > 0, and omega_1(z) ~ z as z -> infinity.

    For the finite polynomial case, omega_1 is determined by:
    F(z,w) = r'(z) * p(w) - p'(w) * r(z) = 0

    At z = nu_k (root of r), F(nu_k, w) = r'(nu_k) * p(w) = 0,
    so w must be a root of p. The matching sigma is determined by the
    analytic continuation / interlacing.

    For our purposes, we need to find which root of p each root of r maps to,
    and compute omega_1' and omega_1'' at each root of r.
    """
    n = len(roots_p)
    roots_p_sorted = np.sort(roots_p)
    roots_r_sorted = np.sort(roots_r)

    # The subordination function omega_1 maps roots of r to roots of p
    # in an order-preserving way (since omega_1 is Herglotz/monotone on reals).
    # So sigma is the identity permutation (when both are sorted).
    sigma = np.arange(n)

    # omega_1(nu_k) = lambda_k (identity permutation)
    omega1_at_roots = roots_p_sorted[sigma]

    # Now compute omega_1'(nu_k) and omega_1''(nu_k) from the implicit function theorem
    # applied to F(z,w) = r'(z)p(w) - p'(w)r(z) = 0.

    p_coeffs = np.poly(roots_p_sorted)
    r_coeffs = np.poly(roots_r_sorted)
    p_prime_coeffs = np.polyder(p_coeffs)
    r_prime_coeffs = np.polyder(r_coeffs)
    p_double_prime_coeffs = np.polyder(p_prime_coeffs)
    r_double_prime_coeffs = np.polyder(r_prime_coeffs)

    omega1_prime = np.zeros(n)
    omega1_double_prime = np.zeros(n)

    for k in range(n):
        z0 = roots_r_sorted[k]
        w0 = omega1_at_roots[k]

        # F(z,w) = r'(z)*p(w) - p'(w)*r(z)
        # F_z = r''(z)*p(w) - p'(w)*r'(z)
        # F_w = r'(z)*p'(w) - p''(w)*r(z)

        rp_z0 = np.polyval(r_prime_coeffs, z0)
        rpp_z0 = np.polyval(r_double_prime_coeffs, z0)
        r_z0 = np.polyval(r_coeffs, z0)  # should be ~0

        pp_w0 = np.polyval(p_prime_coeffs, w0)
        ppp_w0 = np.polyval(p_double_prime_coeffs, w0)
        p_w0 = np.polyval(p_coeffs, w0)  # should be ~0

        # At (nu_k, lambda_{sigma(k)}): r(nu_k) = 0, p(lambda_{sigma(k)}) = 0
        # F_z = r''(z0)*p(w0) - p'(w0)*r'(z0) = 0 - p'(w0)*r'(z0) = -p'(w0)*r'(z0)
        # F_w = r'(z0)*p'(w0) - p''(w0)*r(z0) = r'(z0)*p'(w0) - 0 = r'(z0)*p'(w0)

        F_z = -pp_w0 * rp_z0  # + rpp_z0 * p_w0 (but p_w0 = 0)
        F_w = rp_z0 * pp_w0   # - ppp_w0 * r_z0 (but r_z0 = 0)

        # omega_1'(z) = -F_z / F_w
        omega1_prime[k] = -F_z / F_w  # = pp_w0 * rp_z0 / (rp_z0 * pp_w0) = 1

        # For omega_1''(z), we need second derivatives of F
        # Using: omega_1''(z) = -[F_zz + 2*F_zw*omega_1' + F_ww*(omega_1')^2] / F_w

        rppp_z0 = np.polyval(np.polyder(r_double_prime_coeffs), z0) if len(r_double_prime_coeffs) > 1 else 0
        pppp_w0 = np.polyval(np.polyder(p_double_prime_coeffs), w0) if len(p_double_prime_coeffs) > 1 else 0

        # F_zz = r'''(z)*p(w) + r''(z)*0 - p'(w)*r''(z) - 0  [at the root point]
        # Wait, let me be more careful with F = r'(z)*p(w) - p'(w)*r(z)
        # F_z = r''(z)*p(w) - p'(w)*r'(z)
        # F_zz = r'''(z)*p(w) - p'(w)*r''(z)
        # At root point: p(w0)=0, r(z0)=0
        F_zz = rppp_z0 * p_w0 - pp_w0 * rpp_z0  # = 0 - pp_w0 * rpp_z0
        F_zz = -pp_w0 * rpp_z0

        # F_zw = r''(z)*p'(w) - p''(w)*r'(z)
        # At root point:
        F_zw = rpp_z0 * pp_w0 - ppp_w0 * rp_z0

        # F_ww = r'(z)*p''(w) - p'''(w)*r(z)
        # At root point: r(z0)=0
        F_ww = rp_z0 * ppp_w0 - pppp_w0 * r_z0  # = rp_z0 * ppp_w0
        F_ww = rp_z0 * ppp_w0

        w1 = omega1_prime[k]
        omega1_double_prime[k] = -(F_zz + 2*F_zw*w1 + F_ww*w1**2) / F_w

    return omega1_at_roots, omega1_prime, omega1_double_prime, sigma


def run_test(roots_p, roots_q, label=""):
    """Run a complete test for given p and q."""
    n = len(roots_p)

    # Compute r = p boxplus_n q using Monte Carlo
    roots_r = free_additive_convolution_roots(roots_p, roots_q)

    # Compute H values
    H_p = H_values(roots_p)
    H_q = H_values(roots_q)
    H_r = H_values(roots_r)

    Phi_p = np.sum(H_p**2)
    Phi_q = np.sum(H_q**2)
    Phi_r = np.sum(H_r**2)

    # Compute subordination
    omega1_vals, omega1_p, omega1_pp, sigma = subordination_omega1(roots_p, roots_r)
    omega2_vals, omega2_p, omega2_pp, tau = subordination_omega1(roots_q, roots_r)

    # Vectors
    h = H_r  # h_k = H_r(nu_k)
    u = H_p[sigma]  # u_k = H_p(lambda_{sigma(k)})
    v = H_q[tau]    # v_k = H_q(mu_{tau(k)})

    alpha = omega1_pp / 2  # alpha_k = omega_1''(nu_k) / 2
    beta = omega2_pp / 2   # beta_k = omega_2''(nu_k) / 2

    # Check chain rule: u = h + alpha?
    chain_rule_error_1 = np.max(np.abs(u - (h + alpha)))
    chain_rule_error_2 = np.max(np.abs(v - (h + beta)))

    # Compute key quantities
    h_alpha = np.dot(h, alpha)
    h_beta = np.dot(h, beta)

    norm_h_sq = np.dot(h, h)
    norm_u_sq = np.dot(u, u)
    norm_v_sq = np.dot(v, v)
    norm_alpha_sq = np.dot(alpha, alpha)
    norm_beta_sq = np.dot(beta, beta)

    A = norm_u_sq - norm_h_sq  # = 2*h_alpha + norm_alpha_sq
    B = norm_v_sq - norm_h_sq  # = 2*h_beta + norm_beta_sq

    AB = A * B
    target = norm_h_sq**2

    # Also check: h_alpha = (A - norm_alpha_sq) / 2
    h_alpha_check = (A - norm_alpha_sq) / 2
    h_beta_check = (B - norm_beta_sq) / 2

    # Main inequality
    main_ineq = 1.0/Phi_r - (1.0/Phi_p + 1.0/Phi_q)

    print(f"\n{'='*60}")
    print(f"TEST: {label}")
    print(f"n = {n}")
    print(f"roots_p = {np.sort(roots_p)}")
    print(f"roots_q = {np.sort(roots_q)}")
    print(f"roots_r = {roots_r}")
    print(f"")
    print(f"Chain rule errors: {chain_rule_error_1:.6e}, {chain_rule_error_2:.6e}")
    print(f"omega1' at roots: {omega1_p}")
    print(f"omega2' at roots: {omega2_p}")
    print(f"")
    print(f"h = {h}")
    print(f"u = {u}")
    print(f"v = {v}")
    print(f"alpha = {alpha}")
    print(f"beta  = {beta}")
    print(f"")
    print(f"Phi_n(p) = {Phi_p:.8f}, Phi_n(q) = {Phi_q:.8f}, Phi_n(r) = {Phi_r:.8f}")
    print(f"||h||^2 = {norm_h_sq:.8f}  (should = Phi_r)")
    print(f"||u||^2 = {norm_u_sq:.8f}  (should = Phi_p)")
    print(f"||v||^2 = {norm_v_sq:.8f}  (should = Phi_q)")
    print(f"")
    print(f"<h,alpha> = {h_alpha:.8f}  (check: {h_alpha_check:.8f})")
    print(f"<h,beta>  = {h_beta:.8f}  (check: {h_beta_check:.8f})")
    print(f"<h,alpha> >= 0? {h_alpha >= -1e-8}")
    print(f"<h,beta>  >= 0? {h_beta >= -1e-8}")
    print(f"")
    print(f"A = Phi_p - Phi_r = {A:.8f}")
    print(f"B = Phi_q - Phi_r = {B:.8f}")
    print(f"AB = {AB:.8f}")
    print(f"||h||^4 = Phi_r^2 = {target:.8f}")
    print(f"AB - ||h||^4 = {AB - target:.8f}")
    print(f"AB >= ||h||^4? {AB >= target - 1e-8}")
    print(f"")
    print(f"Main inequality 1/Phi_r - 1/Phi_p - 1/Phi_q = {main_ineq:.8f} >= 0? {main_ineq >= -1e-8}")
    print(f"Equivalent: AB >= ||h||^4? {AB >= target - 1e-8}")

    # Extra diagnostics
    print(f"\n--- Extra diagnostics ---")
    print(f"||alpha||^2 = {norm_alpha_sq:.8f}")
    print(f"||beta||^2  = {norm_beta_sq:.8f}")
    print(f"||alpha||*||beta|| = {np.sqrt(norm_alpha_sq*norm_beta_sq):.8f}")
    print(f"<alpha,beta> = {np.dot(alpha,beta):.8f}")
    print(f"<h,alpha>*<h,beta> = {h_alpha*h_beta:.8f}")
    print(f"||h||^4/4 = {target/4:.8f}")
    print(f"<h,alpha>*<h,beta> >= ||h||^4/4? {h_alpha*h_beta >= target/4 - 1e-8}")

    # CS-type bound: <h,alpha>^2 <= ||h||^2 * ||alpha||^2
    cs_alpha = h_alpha**2 / (norm_h_sq * norm_alpha_sq) if norm_alpha_sq > 1e-15 else float('inf')
    cs_beta = h_beta**2 / (norm_h_sq * norm_beta_sq) if norm_beta_sq > 1e-15 else float('inf')
    print(f"<h,alpha>^2 / (||h||^2*||alpha||^2) = {cs_alpha:.8f}  (CS says <= 1)")
    print(f"<h,beta>^2 / (||h||^2*||beta||^2) = {cs_beta:.8f}  (CS says <= 1)")

    # Angle analysis
    if norm_alpha_sq > 1e-15 and norm_h_sq > 1e-15:
        cos_h_alpha = h_alpha / np.sqrt(norm_h_sq * norm_alpha_sq)
        print(f"cos(h,alpha) = {cos_h_alpha:.8f}")
    if norm_beta_sq > 1e-15 and norm_h_sq > 1e-15:
        cos_h_beta = h_beta / np.sqrt(norm_h_sq * norm_beta_sq)
        print(f"cos(h,beta) = {cos_h_beta:.8f}")

    return {
        'h_alpha': h_alpha, 'h_beta': h_beta,
        'A': A, 'B': B, 'AB': AB, 'h4': target,
        'gap': AB - target, 'main': main_ineq,
        'Phi_p': Phi_p, 'Phi_q': Phi_q, 'Phi_r': Phi_r,
        'alpha': alpha, 'beta': beta, 'h': h, 'u': u, 'v': v,
        'chain_err_1': chain_rule_error_1, 'chain_err_2': chain_rule_error_2,
    }


# ===== TESTS =====

print("=" * 60)
print("NUMERICAL INVESTIGATION FOR NODES 1.7 AND 1.8")
print("=" * 60)

np.random.seed(42)

# Test 1: n=2, simple case
results_n2 = run_test(np.array([-1.0, 1.0]), np.array([-2.0, 2.0]), "n=2, symmetric")

# Test 2: n=2, asymmetric
results_n2b = run_test(np.array([0.0, 1.0]), np.array([0.0, 3.0]), "n=2, asymmetric")

# Test 3: n=3
results_n3 = run_test(np.array([-2.0, 0.0, 2.0]), np.array([-3.0, 0.0, 3.0]), "n=3, symmetric")

# Test 4: n=3, asymmetric
results_n3b = run_test(np.array([-1.0, 0.0, 3.0]), np.array([-2.0, 1.0, 4.0]), "n=3, asymmetric")

# Test 5: n=4
results_n4 = run_test(np.array([-3.0, -1.0, 1.0, 3.0]), np.array([-2.0, -0.5, 0.5, 2.0]), "n=4, symmetric")

# Random tests
print("\n\n" + "=" * 60)
print("RANDOM TESTS: Sign of <h,alpha> and <h,beta>")
print("=" * 60)

h_alpha_signs = []
h_beta_signs = []
gap_signs = []
h_alpha_neg_examples = []

for trial in range(200):
    n = np.random.choice([2, 3, 4, 5])
    roots_p = np.sort(np.random.randn(n) * 2)
    # Ensure distinct
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.3:
            roots_p[i] = roots_p[i-1] + 0.3

    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.3:
            roots_q[i] = roots_q[i-1] + 0.3

    try:
        res = run_test.__wrapped__(roots_p, roots_q) if hasattr(run_test, '__wrapped__') else None

        roots_r = free_additive_convolution_roots(roots_p, roots_q)
        H_r_vals = H_values(roots_r)

        _, _, omega1_pp, sigma = subordination_omega1(roots_p, roots_r)
        _, _, omega2_pp, tau = subordination_omega1(roots_q, roots_r)

        h = H_r_vals
        alpha = omega1_pp / 2
        beta = omega2_pp / 2
        u = H_values(roots_p)[sigma]
        v = H_values(roots_q)[tau]

        h_alpha = np.dot(h, alpha)
        h_beta = np.dot(h, beta)

        A = np.dot(u, u) - np.dot(h, h)
        B = np.dot(v, v) - np.dot(h, h)
        AB = A * B
        h4 = np.dot(h, h)**2

        h_alpha_signs.append(h_alpha >= -1e-6)
        h_beta_signs.append(h_beta >= -1e-6)
        gap_signs.append(AB >= h4 - 1e-6)

        if h_alpha < -1e-6:
            h_alpha_neg_examples.append({
                'trial': trial, 'n': n,
                'roots_p': roots_p.copy(), 'roots_q': roots_q.copy(),
                'h_alpha': h_alpha, 'h_beta': h_beta,
                'A': A, 'B': B, 'AB': AB, 'h4': h4,
            })
    except Exception as e:
        pass

print(f"\nTotal trials: {len(h_alpha_signs)}")
print(f"<h,alpha> >= 0: {sum(h_alpha_signs)} / {len(h_alpha_signs)} = {sum(h_alpha_signs)/len(h_alpha_signs)*100:.1f}%")
print(f"<h,beta> >= 0: {sum(h_beta_signs)} / {len(h_beta_signs)} = {sum(h_beta_signs)/len(h_beta_signs)*100:.1f}%")
print(f"AB >= ||h||^4: {sum(gap_signs)} / {len(gap_signs)} = {sum(gap_signs)/len(gap_signs)*100:.1f}%")

if h_alpha_neg_examples:
    print(f"\n<h,alpha> < 0 counterexamples:")
    for ex in h_alpha_neg_examples[:5]:
        print(f"  Trial {ex['trial']}: n={ex['n']}, <h,alpha>={ex['h_alpha']:.6f}, A={ex['A']:.6f}, B={ex['B']:.6f}, AB={ex['AB']:.6f}, h4={ex['h4']:.6f}")
        print(f"    roots_p={ex['roots_p']}, roots_q={ex['roots_q']}")
