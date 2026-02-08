"""
PROVER-13: Heat Flow / de Bruijn Identity Approach for Fisher Superadditivity

Investigates whether the classical proof strategy for Stam's inequality
(via de Bruijn identity + heat flow) has a finite free analog.

Key questions:
1. Compute Phi_n(G_s) for the Gaussian polynomial G_s
2. Test monotonicity of Phi_n under finite free convolution with G_t
3. Search for a de Bruijn-type identity: Phi_n(p) = -dS_n(p_t)/dt|_{t=0}
4. Test d/dt (1/Phi_n(p_t)) for definite sign
"""

import numpy as np
from numpy.polynomial import polynomial as P
from scipy.special import comb as scipy_comb
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# PART 0: Core infrastructure
# ============================================================

def finite_free_convolution(p_coeffs, q_coeffs):
    """
    Compute p ⊞_n q using MSS formula.
    p_coeffs, q_coeffs are coefficient vectors [a_0, a_1, ..., a_n]
    where p(x) = sum_k a_k x^{n-k}, with a_0 = 1 (monic).
    """
    n = len(p_coeffs) - 1
    assert len(q_coeffs) - 1 == n
    assert abs(p_coeffs[0] - 1.0) < 1e-12
    assert abs(q_coeffs[0] - 1.0) < 1e-12

    from math import factorial
    c = np.zeros(n + 1)
    for k in range(n + 1):
        s = 0.0
        for i in range(k + 1):
            j = k - i
            if i <= n and j <= n:
                coeff = (factorial(n - i) * factorial(n - j)) / (factorial(n) * factorial(n - k))
                s += coeff * p_coeffs[i] * q_coeffs[j]
        c[k] = s
    return c


def coeffs_from_roots(roots):
    """
    Given roots lambda_1, ..., lambda_n, return coefficient vector
    [a_0, a_1, ..., a_n] where p(x) = sum_k a_k x^{n-k}.
    a_0 = 1 (monic).
    """
    # numpy poly gives coefficients in descending order of power
    coeffs = np.polynomial.polynomial.polyfromroots(roots)
    # polyfromroots gives ascending order: c_0 + c_1*x + ... + c_n*x^n
    # We want descending: a_0*x^n + a_1*x^{n-1} + ... + a_n
    coeffs = coeffs[::-1]
    # Normalize to monic
    coeffs = coeffs / coeffs[0]
    return coeffs


def roots_from_coeffs(coeffs):
    """
    Given coefficient vector [a_0, ..., a_n] (descending powers), return roots.
    """
    # np.roots expects descending order coefficients
    return np.sort(np.roots(coeffs)).real


def eval_poly_from_coeffs(coeffs, x):
    """Evaluate polynomial with descending-power coefficients at x."""
    n = len(coeffs) - 1
    val = 0.0
    for k, a in enumerate(coeffs):
        val += a * x**(n - k)
    return val


def H_values(roots):
    """Compute H_p(lambda_i) = sum_{j != i} 1/(lambda_i - lambda_j) for all i."""
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                diff = roots[i] - roots[j]
                if abs(diff) < 1e-15:
                    return None  # repeated root
                H[i] += 1.0 / diff
    return H


def Phi_n(roots):
    """Compute Phi_n(p) = sum_i H_p(lambda_i)^2."""
    H = H_values(roots)
    if H is None:
        return float('inf')
    return np.sum(H**2)


def inv_Phi_n(roots):
    """Compute 1/Phi_n(p)."""
    phi = Phi_n(roots)
    if phi == 0 or phi == float('inf'):
        return 0.0
    return 1.0 / phi


# ============================================================
# PART 1: Gaussian polynomial G_s
# ============================================================

def hermite_prob_coeffs(n):
    """
    Return coefficient vector for probabilist's Hermite polynomial He_n(x).
    He_n(x) = n! * sum_{m=0}^{floor(n/2)} (-1)^m / (m! * 2^m * (n-2m)!) * x^{n-2m}
    Returned as [a_0, a_1, ..., a_n] with descending powers.
    """
    from math import factorial
    coeffs = np.zeros(n + 1)
    for m in range(n // 2 + 1):
        k = 2 * m  # this contributes to x^{n-2m}, which is position 2m in descending
        coeff = factorial(n) * ((-1)**m) / (factorial(m) * (2**m) * factorial(n - 2*m))
        coeffs[k] = coeff
    return coeffs


def gaussian_poly_coeffs(n, s):
    """
    Gaussian polynomial G_s(x) = s^{n/2} * He_n(x / sqrt(s)).
    Returns coefficient vector [a_0, ..., a_n] in descending powers.

    G_s(x) = sum_k a_k(He) * s^{n/2} * (x/sqrt(s))^{n-k}
           = sum_k a_k(He) * s^{n/2} * x^{n-k} * s^{-(n-k)/2}
           = sum_k a_k(He) * s^{k/2} * x^{n-k}
    """
    he = hermite_prob_coeffs(n)
    gs = np.zeros(n + 1)
    for k in range(n + 1):
        gs[k] = he[k] * s**(k / 2.0)
    return gs


def gaussian_poly_roots(n, s):
    """Roots of G_s = sqrt(s) * (Hermite zeros)."""
    he_coeffs = hermite_prob_coeffs(n)
    he_roots = np.sort(np.roots(he_coeffs)).real
    return np.sort(he_roots * np.sqrt(s))


print("=" * 70)
print("PROVER-13: HEAT FLOW / DE BRUIJN IDENTITY APPROACH")
print("=" * 70)

# ============================================================
# TEST 1: Verify Phi_n(G_s) = n(n-1)/(4s)
# ============================================================

print("\n" + "=" * 70)
print("TEST 1: Phi_n(G_s) for Gaussian polynomial")
print("=" * 70)

for n in [3, 4, 5, 6, 8, 10]:
    for s in [0.5, 1.0, 2.0, 5.0]:
        roots = gaussian_poly_roots(n, s)
        phi = Phi_n(roots)
        predicted = n * (n - 1) / (4 * s)
        ratio = phi / predicted if predicted > 0 else float('nan')
        if abs(ratio - 1.0) > 0.01:
            print(f"  n={n}, s={s}: Phi_n = {phi:.6f}, predicted = {predicted:.6f}, RATIO = {ratio:.6f} *** MISMATCH")
        else:
            print(f"  n={n}, s={s}: Phi_n = {phi:.6f}, predicted = {predicted:.6f}, ratio = {ratio:.6f} OK")

# ============================================================
# TEST 2: Monotonicity of Phi_n under convolution with G_t
# ============================================================

print("\n" + "=" * 70)
print("TEST 2: Does Phi_n decrease under convolution with G_t?")
print("=" * 70)
print("Testing: Phi_n(p ⊞_n G_t) <= Phi_n(p) for all t > 0?")
print("(This is the finite free analog of 'adding noise decreases Fisher info')")

np.random.seed(42)

violations_monotone = 0
total_monotone = 0

for n in [3, 4, 5]:
    for trial in range(200):
        # Random polynomial with well-separated roots
        roots_p = np.sort(np.random.randn(n) * 2)
        # Ensure simple roots
        min_gap = np.min(np.diff(roots_p))
        if min_gap < 0.1:
            continue

        coeffs_p = coeffs_from_roots(roots_p)
        phi_p = Phi_n(roots_p)

        for t in [0.01, 0.1, 0.5, 1.0, 2.0]:
            coeffs_gt = gaussian_poly_coeffs(n, t)
            coeffs_pt = finite_free_convolution(coeffs_p, coeffs_gt)
            roots_pt = roots_from_coeffs(coeffs_pt)

            if np.any(np.iscomplex(roots_pt)):
                continue
            if np.min(np.diff(np.sort(roots_pt.real))) < 1e-10:
                continue

            phi_pt = Phi_n(roots_pt.real)
            total_monotone += 1

            if phi_pt > phi_p + 1e-8:
                violations_monotone += 1
                if violations_monotone <= 5:
                    print(f"  VIOLATION: n={n}, t={t:.2f}, Phi(p)={phi_p:.4f}, Phi(p⊞G_t)={phi_pt:.4f}")

print(f"\nMonotonicity test: {violations_monotone} violations in {total_monotone} trials")
if violations_monotone == 0:
    print("*** Phi_n DECREASES under convolution with Gaussian -- CONFIRMED ***")

# ============================================================
# TEST 3: Compute the "entropy" S_n(p_t) = -integral of Phi
# ============================================================

print("\n" + "=" * 70)
print("TEST 3: Heat flow trajectory and entropy candidate")
print("=" * 70)
print("For p_t = p ⊞_n G_t, compute Phi_n(p_t) as function of t")
print("Check if there exists S_n with S_n'(t) = -Phi_n(p_t)")

# For several test polynomials, trace the heat flow
for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    # Test polynomial: roots at 1, 2, ..., n
    roots_p = np.arange(1, n + 1, dtype=float)
    coeffs_p = coeffs_from_roots(roots_p)

    t_values = np.linspace(0.001, 5.0, 200)
    phi_values = []
    inv_phi_values = []

    for t in t_values:
        coeffs_gt = gaussian_poly_coeffs(n, t)
        coeffs_pt = finite_free_convolution(coeffs_p, coeffs_gt)
        roots_pt = roots_from_coeffs(coeffs_pt)
        if np.any(np.abs(np.imag(roots_pt)) > 1e-8):
            phi_values.append(float('nan'))
            inv_phi_values.append(float('nan'))
            continue
        roots_pt = np.sort(roots_pt.real)
        phi = Phi_n(roots_pt)
        phi_values.append(phi)
        inv_phi_values.append(1.0 / phi if phi > 0 else float('nan'))

    phi_values = np.array(phi_values)
    inv_phi_values = np.array(inv_phi_values)

    # Check monotonicity of Phi
    phi_diffs = np.diff(phi_values)
    phi_decreasing = np.all(phi_diffs[~np.isnan(phi_diffs)] <= 1e-8)
    print(f"  Phi_n(p_t) decreasing in t: {phi_decreasing}")

    # Check monotonicity of 1/Phi
    inv_phi_diffs = np.diff(inv_phi_values)
    inv_phi_increasing = np.all(inv_phi_diffs[~np.isnan(inv_phi_diffs)] >= -1e-8)
    print(f"  1/Phi_n(p_t) increasing in t: {inv_phi_increasing}")

    # Check if 1/Phi_n(p_t) is concave in t
    inv_phi_dd = np.diff(inv_phi_values, 2)
    inv_phi_concave = np.all(inv_phi_dd[~np.isnan(inv_phi_dd)] <= 1e-6)
    print(f"  1/Phi_n(p_t) concave in t: {inv_phi_concave}")

    # Check large-t asymptotics
    # As t -> inf, p ⊞_n G_t should look like G_t (CLT), so Phi_n ~ n(n-1)/(4t)
    # Hence 1/Phi_n ~ 4t/(n(n-1))
    for idx in [-5, -3, -1]:
        t_val = t_values[idx]
        inv_phi_val = inv_phi_values[idx]
        predicted = 4 * t_val / (n * (n - 1))
        print(f"  t={t_val:.2f}: 1/Phi = {inv_phi_val:.6f}, 4t/(n(n-1)) = {predicted:.6f}, ratio = {inv_phi_val/predicted:.4f}")

    # Numerical derivative: d/dt (1/Phi_n)
    dt = t_values[1] - t_values[0]
    d_inv_phi = np.gradient(inv_phi_values, dt)
    # Check if d/dt(1/Phi) is decreasing (i.e., 1/Phi is concave)
    d_inv_phi_diffs = np.diff(d_inv_phi)
    mask = ~np.isnan(d_inv_phi_diffs)
    if np.sum(mask) > 0:
        d_inv_phi_decreasing = np.all(d_inv_phi_diffs[mask] <= 1e-4)
        print(f"  d/dt(1/Phi_n) decreasing (=> concavity): {d_inv_phi_decreasing}")


# ============================================================
# TEST 4: Check de Bruijn-type identity
# ============================================================

print("\n" + "=" * 70)
print("TEST 4: De Bruijn-type identity")
print("=" * 70)
print("Classical: I(X) = -d/dt H(X + sqrt(t)Z)|_{t=0}")
print("Finite free analog: Phi_n(p) = ?? * (-d/dt S_n(p ⊞_n G_t))|_{t=0}")
print()

# If S_n(p) is some entropy, then by the chain rule:
# d/dt S_n(p ⊞_n G_t)|_{t=0} should involve Phi_n(p).
# Let's try S_n(p) = log(discriminant) or similar.

def discriminant(roots):
    """Compute discriminant of polynomial from roots."""
    n = len(roots)
    disc = 1.0
    for i in range(n):
        for j in range(i + 1, n):
            disc *= (roots[j] - roots[i])**2
    return disc


def log_discriminant(roots):
    """log |disc(p)| = 2 * sum_{i<j} log|lambda_j - lambda_i|."""
    n = len(roots)
    val = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            diff = abs(roots[j] - roots[i])
            if diff < 1e-15:
                return float('-inf')
            val += 2 * np.log(diff)
    return val


def sum_log_vandermonde(roots):
    """sum_{i<j} log|lambda_j - lambda_i| = log_disc/2."""
    return log_discriminant(roots) / 2


# Candidate entropies to test
def S_candidate_1(roots):
    """S_1 = (1/n) * sum_{i<j} log|lambda_i - lambda_j|  (Voiculescu-type)"""
    n = len(roots)
    return sum_log_vandermonde(roots) / n


def S_candidate_2(roots):
    """S_2 = sum_{i<j} log|lambda_i - lambda_j|  (= log_disc/2)"""
    return sum_log_vandermonde(roots)


def S_candidate_3(roots):
    """S_3 = (2/n^2) * sum_{i<j} log|lambda_i - lambda_j|  (normalized)"""
    n = len(roots)
    return 2 * sum_log_vandermonde(roots) / (n * n)


print("Testing: -dS/dt|_{t=0} vs c * Phi_n(p) for various S candidates")
print()

for n in [3, 4, 5, 6]:
    print(f"--- n = {n} ---")

    for trial_name, roots_p in [("equispaced", np.arange(1, n + 1, dtype=float)),
                                 ("random1", np.sort(np.array([0.1*i**2 for i in range(1, n+1)]))),
                                 ("random2", np.sort(np.random.RandomState(100+n).randn(n) * 3))]:
        min_gap = np.min(np.diff(roots_p))
        if min_gap < 0.05:
            continue

        coeffs_p = coeffs_from_roots(roots_p)
        phi_p = Phi_n(roots_p)

        # Compute numerical derivative of S(p ⊞ G_t) at t=0+ using small t
        dt = 1e-5
        t_vals = [dt, 2*dt, 3*dt]
        for S_name, S_func in [("S1=logV/n", S_candidate_1),
                                ("S2=logV", S_candidate_2),
                                ("S3=2logV/n^2", S_candidate_3)]:
            S_vals = []
            for t in t_vals:
                coeffs_gt = gaussian_poly_coeffs(n, t)
                coeffs_pt = finite_free_convolution(coeffs_p, coeffs_gt)
                roots_pt = roots_from_coeffs(coeffs_pt)
                if np.any(np.abs(np.imag(roots_pt)) > 1e-8):
                    S_vals.append(float('nan'))
                    continue
                roots_pt = np.sort(roots_pt.real)
                S_vals.append(S_func(roots_pt))

            # Forward difference for dS/dt at t ~ dt
            if len(S_vals) >= 2 and not any(np.isnan(S_vals)):
                dSdt = (S_vals[1] - S_vals[0]) / dt
                # Also need S at t=0
                S0 = S_func(roots_p)
                dSdt_at0 = (S_vals[0] - S0) / dt
                ratio = -dSdt_at0 / phi_p if phi_p > 1e-10 else float('nan')
                if trial_name == "equispaced":
                    print(f"  {trial_name}: {S_name}: -dS/dt|_0 = {-dSdt_at0:.6f}, Phi = {phi_p:.6f}, ratio = {ratio:.6f}")

    print()


# ============================================================
# TEST 5: Superadditivity directly from heat flow
# ============================================================

print("\n" + "=" * 70)
print("TEST 5: Superadditivity from heat flow / concavity")
print("=" * 70)
print("If 1/Phi_n(p_t) is concave in t, then superadditivity may follow")
print("from the semigroup property: (p⊞q)⊞G_t = (p⊞G_t)⊞(q⊞G_t) ??? ")
print("Actually: (p⊞q)⊞G_t = (p⊞G_{t/2})⊞(q⊞G_{t/2}) is NOT true in general.")
print("But: G_s ⊞ G_t = G_{s+t} IS true (check below)")
print()

# Verify: G_s ⊞_n G_t = G_{s+t}
print("Verifying G_s ⊞_n G_t = G_{s+t}:")
for n in [3, 4, 5, 6]:
    for s, t in [(1.0, 1.0), (0.5, 2.0), (1.0, 3.0)]:
        gs = gaussian_poly_coeffs(n, s)
        gt = gaussian_poly_coeffs(n, t)
        gst_conv = finite_free_convolution(gs, gt)
        gst_direct = gaussian_poly_coeffs(n, s + t)
        err = np.max(np.abs(gst_conv - gst_direct))
        status = "OK" if err < 1e-10 else f"ERR={err:.2e}"
        if n <= 4:
            print(f"  n={n}, s={s}, t={t}: max_err = {err:.2e} {status}")

print()

# ============================================================
# TEST 6: The key semigroup identity for the proof
# ============================================================

print("\n" + "=" * 70)
print("TEST 6: Key identity for Stam-type proof")
print("=" * 70)
print()
print("Classical Stam proof uses: X+Y = (X + sqrt(t)Z_1) + (Y + sqrt(t)Z_2) - sqrt(2t)Z_3")
print("where Z_i are independent Gaussians.")
print()
print("In finite free setting, the key is:")
print("  p ⊞ q ⊞ G_t = (p ⊞ G_s) ⊞ (q ⊞ G_{t-s})  ??? NO, G distributes differently")
print()
print("Actually the correct identity is simply:")
print("  (p ⊞ q) ⊞ G_t = p ⊞ q ⊞ G_t  (associativity)")
print("and G_s ⊞ G_t = G_{s+t}")
print()
print("The data processing approach:")
print("  Phi_n(p ⊞ G_t) <= Phi_n(p)  (monotonicity/data processing)")
print("  1/Phi_n(p ⊞ G_t) >= 1/Phi_n(p)  (inverse monotonicity)")
print()
print("For Stam, we need: 1/Phi(p⊞q) >= 1/Phi(p) + 1/Phi(q)")
print("Attempted via: at t -> inf, p ⊞ G_t ~ G_t, so 1/Phi ~ 4t/n(n-1)")
print()

# ============================================================
# TEST 7: Derivative formula for 1/Phi along heat flow
# ============================================================

print("=" * 70)
print("TEST 7: d/dt(1/Phi_n(p ⊞ G_t)) at t=0")
print("=" * 70)
print("If there's a universal formula for d/dt(1/Phi_n(p_t))|_{t=0} in terms")
print("of the roots of p, this could yield the de Bruijn analog.")
print()

np.random.seed(123)

for n in [3, 4, 5]:
    print(f"--- n = {n} ---")
    results = []
    for trial in range(20):
        roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n) * 2)
        min_gap = np.min(np.diff(roots_p))
        if min_gap < 0.2:
            continue

        coeffs_p = coeffs_from_roots(roots_p)
        phi_p = Phi_n(roots_p)

        # Numerical derivative of 1/Phi at t=0+
        dt = 1e-6
        inv_phis = []
        for t in [0, dt, 2*dt, 3*dt]:
            if t == 0:
                inv_phis.append(1.0 / phi_p)
            else:
                coeffs_gt = gaussian_poly_coeffs(n, t)
                coeffs_pt = finite_free_convolution(coeffs_p, coeffs_gt)
                roots_pt = roots_from_coeffs(coeffs_pt)
                if np.any(np.abs(np.imag(roots_pt)) > 1e-8):
                    inv_phis.append(float('nan'))
                    continue
                roots_pt = np.sort(roots_pt.real)
                inv_phis.append(1.0 / Phi_n(roots_pt))

        if any(np.isnan(inv_phis)):
            continue

        # d/dt (1/Phi) at t=0 using forward difference
        d_inv_phi = (-3*inv_phis[0] + 4*inv_phis[1] - inv_phis[2]) / (2*dt)

        # Candidate quantities it might equal:
        H = H_values(roots_p)
        # sum H_i^4 / (sum H_i^2)^2
        ratio1 = np.sum(H**4) / (np.sum(H**2)**2)
        # 1/(sum H_i^2)  (= 1/Phi)
        ratio2 = 1.0 / np.sum(H**2)
        # Various normalizations
        ratio3 = d_inv_phi * phi_p  # = d/dt(1/Phi) * Phi = dimensionless

        results.append({
            'phi': phi_p,
            'd_inv_phi': d_inv_phi,
            'dimensionless': ratio3,
            'sum_H4': np.sum(H**4),
            'sum_H2': np.sum(H**2),
        })

    if results:
        dim_vals = [r['dimensionless'] for r in results]
        phi_vals = [r['phi'] for r in results]
        d_vals = [r['d_inv_phi'] for r in results]
        h4_h2_vals = [r['sum_H4']/r['sum_H2']**2 for r in results]

        print(f"  d/dt(1/Phi)|_0:  min={min(d_vals):.6f}, max={max(d_vals):.6f}, all positive: {all(d > -1e-8 for d in d_vals)}")
        print(f"  Phi * d/dt(1/Phi): min={min(dim_vals):.6f}, max={max(dim_vals):.6f}")
        print(f"  sum(H^4)/sum(H^2)^2: min={min(h4_h2_vals):.6f}, max={max(h4_h2_vals):.6f}")
        print()

# ============================================================
# TEST 8: Free cumulant approach
# ============================================================

print("\n" + "=" * 70)
print("TEST 8: Cumulant decomposition along heat flow")
print("=" * 70)
print("Free cumulants linearize: kappa_k(p ⊞ G_t) = kappa_k(p) + kappa_k(G_t)")
print("For Gaussian: kappa_2(G_t) = t, kappa_k(G_t) = 0 for k >= 3")
print("So: kappa_2(p_t) = kappa_2(p) + t, kappa_k(p_t) = kappa_k(p) for k >= 3")
print()

# Compute finite free cumulants from roots
def power_sums(roots, max_k=None):
    """Compute power sums m_k = (1/n) * sum_i lambda_i^k."""
    n = len(roots)
    if max_k is None:
        max_k = n
    return [np.sum(roots**k) / n for k in range(max_k + 1)]


def coeffs_to_elementary_sym(coeffs):
    """
    From descending-power coefficients [1, a_1, ..., a_n],
    extract elementary symmetric polynomials e_k.
    p(x) = x^n - e_1 x^{n-1} + e_2 x^{n-2} - ...
    So a_k = (-1)^k * e_k.
    """
    n = len(coeffs) - 1
    e = np.zeros(n + 1)
    e[0] = 1
    for k in range(1, n + 1):
        e[k] = (-1)**k * coeffs[k]
    return e


def finite_free_cumulants(coeffs, max_k=None):
    """
    Compute finite free cumulants kappa_k from polynomial coefficients.
    Using the relation: a_k = (-1)^k * binom(n,k) * kappa_k / n^{k-1} approximately.

    Actually, for finite free cumulants (Arizmendi-Perales):
    The k-th finite free cumulant of degree-n polynomial p is:
    kappa_k = n^{k-1} * (-1)^k * a_k / binom(n,k)
    where p(x) = x^n + a_1 x^{n-1} + ... + a_n.
    """
    from math import comb
    n = len(coeffs) - 1
    if max_k is None:
        max_k = n
    kappas = [0.0]  # kappa_0 placeholder
    for k in range(1, min(max_k, n) + 1):
        kappas.append((-1)**k * coeffs[k] * n**(k-1) / comb(n, k))
    return kappas


# Verify cumulant additivity for ⊞
print("Verifying finite free cumulant additivity:")
np.random.seed(77)
for n in [3, 4, 5, 6]:
    roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n))
    roots_q = np.sort(np.random.randn(n) * 1.5 + np.arange(n) * 0.5)
    coeffs_p = coeffs_from_roots(roots_p)
    coeffs_q = coeffs_from_roots(roots_q)
    coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)

    kp = finite_free_cumulants(coeffs_p)
    kq = finite_free_cumulants(coeffs_q)
    kpq = finite_free_cumulants(coeffs_pq)

    print(f"  n={n}:")
    for k in range(1, n + 1):
        expected = kp[k] + kq[k]
        actual = kpq[k]
        err = abs(actual - expected)
        status = "OK" if err < 1e-8 else f"ERR={err:.2e}"
        if k <= 4:
            print(f"    kappa_{k}: p={kp[k]:.4f}, q={kq[k]:.4f}, p+q={expected:.4f}, p⊞q={actual:.4f} {status}")


# Verify Gaussian cumulants: kappa_2 = s, kappa_k = 0 for k >= 3
print("\nGaussian polynomial cumulants:")
for n in [4, 5, 6, 8]:
    for s in [1.0, 2.0]:
        gs = gaussian_poly_coeffs(n, s)
        kg = finite_free_cumulants(gs)
        print(f"  n={n}, s={s}: kappa_1={kg[1]:.6f}, kappa_2={kg[2]:.6f}", end="")
        if len(kg) > 3:
            print(f", kappa_3={kg[3]:.6f}, kappa_4={kg[4]:.6f}", end="")
        print()

# ============================================================
# TEST 9: Phi_n in terms of cumulants along heat flow
# ============================================================

print("\n" + "=" * 70)
print("TEST 9: Phi_n as function of kappa_2 (= kappa_2(p) + t)")
print("=" * 70)
print("Along heat flow: kappa_2 increases by t, all others fixed.")
print("So Phi_n(p_t) = Phi_n(kappa_2 + t, kappa_3, kappa_4, ...)")
print("1/Phi_n(p_t) = f(kappa_2 + t)")
print("If f is concave in kappa_2, then superadditivity might follow!")
print()

# For small n, the cumulant decomposition gives explicit formulas
# n=3: Phi_3 = 12*k2, so 1/Phi_3 = 1/(12*k2) -- concave? NO, convex.
# Actually 1/Phi_3 = 1/(12*k2) is convex in k2, which is bad.
# But wait -- the actual conjecture involves TWO polynomials.

# Let's compute Phi_n(p) in terms of finite free cumulants numerically
print("Phi_n in terms of cumulants (n=3):")
print("For n=3, kappa_1 doesn't affect Phi (translation invariant).")
print("Testing: Phi_3 as function of kappa_2 only (since kappa_3 doesn't exist for n=3)")

# n=3: only kappa_1, kappa_2, kappa_3
# But Phi_3 = sum H_i^2 = sum (sum 1/(l_i - l_j))^2
# For n=3 with cumulants kappa_1=0, kappa_2=k2, kappa_3=k3:
# The polynomial is x^3 - 3*k2*x - 2*k3 (centered, kappa_1=0)

for k2 in [0.5, 1.0, 2.0, 4.0]:
    for k3 in [0.0, 0.1, 0.5]:
        # p(x) = x^3 + a1*x^2 + a2*x + a3
        # kappa_1 = -a1 = 0 => a1 = 0
        # kappa_2 = -a2 * n / binom(n,2) = -a2 * 3 / 3 = -a2 => a2 = -k2
        # kappa_3 = a3 * n^2 / binom(n,3) = a3 * 9 / 1 = 9*a3 => a3 = k3/9
        # Wait, let me recompute. For finite free cumulants:
        # kappa_k = (-1)^k * a_k * n^{k-1} / binom(n,k)
        from math import comb
        n = 3
        a1 = 0  # kappa_1 = 0
        a2 = (-1)**2 * k2 * comb(n, 2) / n**(2-1)  # kappa_2 = (-1)^2 * a2 * n / binom(n,2)
        # Actually: kappa_k = (-1)^k * a_k * n^{k-1} / binom(n,k)
        # So: a_k = (-1)^k * kappa_k * binom(n,k) / n^{k-1}
        a2_correct = (-1)**2 * k2 * comb(3, 2) / 3**1  # = k2 * 3/3 = k2
        a3_correct = (-1)**3 * k3 * comb(3, 3) / 3**2  # = -k3/9

        coeffs = np.array([1.0, a1, a2_correct, a3_correct])
        roots = roots_from_coeffs(coeffs)
        if np.any(np.abs(np.imag(roots)) > 1e-8):
            continue
        roots = np.sort(roots.real)
        if np.min(np.diff(roots)) < 1e-10:
            continue
        phi = Phi_n(roots)
        print(f"  k2={k2:.1f}, k3={k3:.1f}: a2={a2_correct:.4f}, a3={a3_correct:.4f}, Phi_3 = {phi:.6f}")


# ============================================================
# TEST 10: The critical test -- concavity of 1/Phi in kappa_2
# ============================================================

print("\n" + "=" * 70)
print("TEST 10: Concavity of 1/Phi_n(p_t) in t (= kappa_2 direction)")
print("=" * 70)
print("This is the KEY test for the heat flow approach.")
print("If 1/Phi_n is concave in kappa_2, then:")
print("  1/Phi(kappa_2^p + kappa_2^q, ...) >= 1/Phi(kappa_2^p, ...) + 1/Phi(kappa_2^q, ...)")
print("  NO -- that's superlinearity, not what we need.")
print()
print("Actually for Stam, the proof uses:")
print("  1/Phi(p_t) = 1/Phi(p) + integral_0^t (d/ds)(1/Phi(p_s)) ds")
print("  and 1/Phi(q_t) = 1/Phi(q) + integral_0^t (d/ds)(1/Phi(q_s)) ds")
print("  At t = s*, p_s* ⊞ q_s* has the same Phi as (p⊞q)_{2s*}")
print("  The argument is more subtle.")
print()

# Instead, let's test the EPI-type approach directly.
# Define F(t) = 1/Phi_n((p ⊞ q) ⊞ G_t)
#        G(t) = 1/Phi_n(p ⊞ G_t) + 1/Phi_n(q ⊞ G_t)
# We want F(0) >= G(0).
# At t -> inf, F(t) ~ G(t) ~ 4t/n(n-1) (both converge to same thing).
# If F'(0) <= G'(0), then the gap starts increasing -- bad.
# If F'(0) >= G'(0), the gap starts closing from above -- promising.

print("Testing F(t) vs G(t) along heat flow:")
np.random.seed(456)

for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    for trial in range(5):
        roots_p = np.sort(np.random.randn(n) * 1.5 + np.arange(n))
        roots_q = np.sort(np.random.randn(n) * 1.5 + np.arange(n) * 0.7)
        min_gap_p = np.min(np.diff(roots_p))
        min_gap_q = np.min(np.diff(roots_q))
        if min_gap_p < 0.2 or min_gap_q < 0.2:
            continue

        coeffs_p = coeffs_from_roots(roots_p)
        coeffs_q = coeffs_from_roots(roots_q)
        coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)

        t_vals = np.linspace(0.001, 3.0, 100)
        F_vals = []
        G_vals = []

        for t in t_vals:
            gt = gaussian_poly_coeffs(n, t)

            # F(t) = 1/Phi((p⊞q)⊞G_t)
            c_pq_t = finite_free_convolution(coeffs_pq, gt)
            r_pq_t = roots_from_coeffs(c_pq_t)
            if np.any(np.abs(np.imag(r_pq_t)) > 1e-8):
                F_vals.append(float('nan'))
                G_vals.append(float('nan'))
                continue
            F_vals.append(inv_Phi_n(np.sort(r_pq_t.real)))

            # G(t) = 1/Phi(p⊞G_t) + 1/Phi(q⊞G_t)
            c_p_t = finite_free_convolution(coeffs_p, gt)
            c_q_t = finite_free_convolution(coeffs_q, gt)
            r_p_t = roots_from_coeffs(c_p_t)
            r_q_t = roots_from_coeffs(c_q_t)
            if np.any(np.abs(np.imag(r_p_t)) > 1e-8) or np.any(np.abs(np.imag(r_q_t)) > 1e-8):
                G_vals.append(float('nan'))
                continue
            G_vals.append(inv_Phi_n(np.sort(r_p_t.real)) + inv_Phi_n(np.sort(r_q_t.real)))

        F_vals = np.array(F_vals)
        G_vals = np.array(G_vals)
        gap = F_vals - G_vals

        valid = ~np.isnan(gap)
        if np.sum(valid) < 10:
            continue

        gap_valid = gap[valid]
        t_valid = t_vals[valid]

        print(f"  Trial {trial}: gap min={np.min(gap_valid):.6f}, gap at t=0.001: {gap_valid[0]:.6f}")
        print(f"    F(0.001)={F_vals[valid][0]:.6f}, G(0.001)={G_vals[valid][0]:.6f}")
        print(f"    gap decreasing? {np.all(np.diff(gap_valid) <= 1e-6)}")

        # Check if gap is monotonically decreasing
        gap_diffs = np.diff(gap_valid)
        n_decrease = np.sum(gap_diffs < -1e-10)
        n_increase = np.sum(gap_diffs > 1e-10)
        print(f"    gap changes: {n_decrease} decrease, {n_increase} increase")


# ============================================================
# SUMMARY
# ============================================================

print("\n" + "=" * 70)
print("SUMMARY OF HEAT FLOW INVESTIGATION")
print("=" * 70)
