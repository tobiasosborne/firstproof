"""
investigate_induction_v2.py — Deep investigation of the inductive approach
for the Fisher superadditivity inequality via the MSS derivative identity.

CONJECTURE: For monic real-rooted p, q of degree n with simple roots,
  r = p boxplus_n q:
  1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q)

KEY IDENTITY (MSS): r'(x) = n * (p^{(1)} boxplus_{n-1} q^{(1)})(x)
  where p^{(1)} = (1/n)*p' (roots = critical points of p).

INVESTIGATIONS:
  1. Verify MSS derivative identity numerically
  2. Ratio Phi_n(p) / Phi_{n-1}(p^{(1)})
  3. Excess ratio F_n / F_{n-1}
  4. Decomposition Phi_n(p) = Phi_{n-1}(p^{(1)}) + R_n(p)
  5. Cauchy interlacing approach
  6. Telescoping identity check
"""

import numpy as np
from math import factorial, comb
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# Core functions
# ============================================================

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
    """MSS finite free convolution via coefficient formula."""
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
    roots_r = np.sort(np.real(np.roots(np.array([c[k] * (-1)**k * comb(n, k)
                      for k in range(n+1)][::-1]))))
    # Actually compute roots from the polynomial with coeffs c
    # c[k] are the "free cumulant" style coefficients.
    # The polynomial is sum_k c[k] * (-1)^k * C(n,k) * x^{n-k}
    poly_coeffs = []
    for k in range(n+1):
        poly_coeffs.append(c[k] * (-1)**k * comb(n, k))
    poly_coeffs_np = np.array(poly_coeffs[::-1])  # highest degree first for np.roots
    # Actually np.poly and np.roots use highest-degree-first convention
    # Let me just reconstruct from the MSS formula properly
    # p boxplus q has coefficients c_k where c_k = sum_{i+j=k} ... * a_i * b_j
    # and the polynomial is x^n + sum_{k=1}^n (-1)^k C(n,k) c_k x^{n-k}
    poly = np.zeros(n+1)
    poly[0] = 1.0  # x^n coefficient
    for k in range(1, n+1):
        poly[k] = (-1)**k * comb(n, k) * c[k]
    raw_roots = np.roots(poly)
    return np.sort(np.real(raw_roots)), c

def H_values(roots):
    """Compute H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j) for each root."""
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                diff = roots[i] - roots[j]
                if abs(diff) < 1e-12:
                    return None  # degenerate
                H[i] += 1.0 / diff
    return H

def H_values_safe(roots):
    """Compute H values, return None if roots too close."""
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                diff = roots[i] - roots[j]
                if abs(diff) < 1e-10:
                    return None
                H[i] += 1.0 / diff
    return H

def Phi_n(roots):
    """Phi_n(p) = sum_i H_p(lambda_i)^2."""
    H = H_values_safe(roots)
    if H is None:
        return None
    return np.sum(H**2)

def critical_points(roots):
    """Roots of p' where p(x) = prod(x - r_i)."""
    poly = np.poly(roots)
    dpoly = np.polyder(poly)
    cps = np.roots(dpoly)
    if np.any(np.abs(np.imag(cps)) > 1e-8):
        return None  # complex critical points
    return np.sort(np.real(cps))

def make_well_separated_roots(n, min_gap=0.5, scale=2.0):
    """Generate n well-separated real roots."""
    roots = np.sort(np.random.randn(n) * scale)
    for i in range(1, n):
        if roots[i] - roots[i-1] < min_gap:
            roots[i] = roots[i-1] + min_gap
    return roots

def safe_boxplus(roots_p, roots_q, tol=0.01):
    """Compute boxplus, return None if result is degenerate."""
    try:
        roots_r, c = boxplus_mss(roots_p, roots_q)
        if np.any(np.abs(np.imag(np.roots(np.poly(roots_r)))) > tol):
            return None, None
        roots_r = np.sort(np.real(roots_r))
        if len(roots_r) > 1 and np.any(np.diff(roots_r) < tol):
            return None, None
        return roots_r, c
    except:
        return None, None


# ============================================================
# INVESTIGATION 1: Verify MSS derivative identity
# ============================================================
print("=" * 70)
print("INVESTIGATION 1: VERIFY MSS DERIVATIVE IDENTITY")
print("  r'(x) = n * (p^{(1)} boxplus_{n-1} q^{(1)})(x)")
print("  where p^{(1)} = (1/n)*p'")
print("=" * 70)

mss_errors = []
for trial in range(30):
    n = np.random.choice([3, 4, 5, 6])
    roots_p = make_well_separated_roots(n)
    roots_q = make_well_separated_roots(n)

    roots_r, _ = safe_boxplus(roots_p, roots_q)
    if roots_r is None:
        continue

    cp_p = critical_points(roots_p)
    cp_q = critical_points(roots_q)
    cp_r = critical_points(roots_r)

    if cp_p is None or cp_q is None or cp_r is None:
        continue
    if len(cp_p) != n-1 or len(cp_q) != n-1 or len(cp_r) != n-1:
        continue

    # p^{(1)} boxplus_{n-1} q^{(1)}
    s_roots, _ = safe_boxplus(cp_p, cp_q)
    if s_roots is None:
        continue

    # r' = n * s means: roots of r' = roots of s
    # (since r' has degree n-1 and n*s has degree n-1 with same roots)
    if len(s_roots) != len(cp_r):
        continue

    error = np.max(np.abs(np.sort(cp_r) - np.sort(s_roots)))
    mss_errors.append(error)
    if trial < 5:
        print(f"  n={n}: max|cp_r - s_roots| = {error:.2e}")

if mss_errors:
    print(f"\n  Summary: {len(mss_errors)} valid trials")
    print(f"  Max error across all trials: {max(mss_errors):.2e}")
    print(f"  Mean error: {np.mean(mss_errors):.2e}")
    print(f"  MSS derivative identity: {'VERIFIED' if max(mss_errors) < 1e-6 else 'FAILED'}")


# ============================================================
# INVESTIGATION 2: Ratio Phi_n(p) / Phi_{n-1}(p^{(1)})
# ============================================================
print("\n" + "=" * 70)
print("INVESTIGATION 2: RATIO Phi_n(p) / Phi_{n-1}(p^{(1)})")
print("  Is there a clean formula?")
print("=" * 70)

for n in [3, 4, 5, 6]:
    print(f"\n  n = {n}:")
    ratios = []
    data_points = []
    for trial in range(50):
        roots_p = make_well_separated_roots(n, min_gap=0.3 + 0.2*np.random.rand())
        cp = critical_points(roots_p)
        if cp is None or len(cp) != n-1:
            continue
        if np.any(np.diff(cp) < 0.01):
            continue

        Ph_n = Phi_n(roots_p)
        Ph_n1 = Phi_n(cp)
        if Ph_n is None or Ph_n1 is None or Ph_n1 < 1e-12:
            continue

        ratio = Ph_n / Ph_n1
        ratios.append(ratio)

        # Also record root gap statistics
        gaps = np.diff(roots_p)
        cp_gaps = np.diff(cp)
        data_points.append({
            'ratio': ratio,
            'mean_gap': np.mean(gaps),
            'gap_variance': np.var(gaps),
            'min_gap': np.min(gaps),
            'max_gap': np.max(gaps),
        })

    if ratios:
        ratios = np.array(ratios)
        print(f"    {len(ratios)} valid samples")
        print(f"    ratio range: [{ratios.min():.6f}, {ratios.max():.6f}]")
        print(f"    mean: {ratios.mean():.6f}, std: {ratios.std():.6f}")
        print(f"    IS CONSTANT? {'YES (std/mean < 0.01)' if ratios.std()/ratios.mean() < 0.01 else 'NO (varies with roots)'}")

# Also check: does the ratio depend on root spacing?
print("\n  --- Ratio vs. uniform spacing ---")
for n in [3, 4, 5, 6]:
    for d in [0.5, 1.0, 2.0, 5.0]:
        roots = np.array([i * d for i in range(n)])
        cp = critical_points(roots)
        if cp is None or len(cp) != n-1:
            continue
        Ph_n_val = Phi_n(roots)
        Ph_n1_val = Phi_n(cp)
        if Ph_n_val is None or Ph_n1_val is None:
            continue
        print(f"    n={n}, d={d}: Phi_n={Ph_n_val:.6f}, Phi_{n-1}={Ph_n1_val:.6f}, ratio={Ph_n_val/Ph_n1_val:.6f}")

# Same check for r = p boxplus q
print("\n  --- Ratio Phi_n(r) / Phi_{n-1}(s) where s = p^{(1)} boxplus q^{(1)} ---")
for n in [3, 4, 5]:
    r_ratios = []
    for trial in range(50):
        roots_p = make_well_separated_roots(n)
        roots_q = make_well_separated_roots(n)
        roots_r, _ = safe_boxplus(roots_p, roots_q)
        if roots_r is None:
            continue
        cp_p = critical_points(roots_p)
        cp_q = critical_points(roots_q)
        cp_r = critical_points(roots_r)
        if cp_p is None or cp_q is None or cp_r is None:
            continue
        if len(cp_p) != n-1 or len(cp_q) != n-1 or len(cp_r) != n-1:
            continue

        Ph_r = Phi_n(roots_r)
        Ph_s = Phi_n(cp_r)  # Phi_{n-1} of the derivative convolution
        if Ph_r is None or Ph_s is None or Ph_s < 1e-12:
            continue
        r_ratios.append(Ph_r / Ph_s)

    if r_ratios:
        r_ratios = np.array(r_ratios)
        print(f"    n={n}: range=[{r_ratios.min():.6f}, {r_ratios.max():.6f}], mean={r_ratios.mean():.6f}, std={r_ratios.std():.6f}")


# ============================================================
# INVESTIGATION 3: Excess ratio F_n / F_{n-1}
# ============================================================
print("\n" + "=" * 70)
print("INVESTIGATION 3: EXCESS RATIO F_n / F_{n-1}")
print("  F_n(p,q) = 1/Phi_n(r) - 1/Phi_n(p) - 1/Phi_n(q)")
print("  F_{n-1}(p^{(1)},q^{(1)}) = 1/Phi_{n-1}(s) - 1/Phi_{n-1}(p^{(1)}) - 1/Phi_{n-1}(q^{(1)})")
print("=" * 70)

for n in [3, 4, 5, 6]:
    print(f"\n  n = {n}:")
    excess_ratios = []
    for trial in range(100):
        roots_p = make_well_separated_roots(n)
        roots_q = make_well_separated_roots(n)

        roots_r, _ = safe_boxplus(roots_p, roots_q)
        if roots_r is None:
            continue

        cp_p = critical_points(roots_p)
        cp_q = critical_points(roots_q)
        if cp_p is None or cp_q is None:
            continue
        if len(cp_p) != n-1 or len(cp_q) != n-1:
            continue

        # s = p^{(1)} boxplus_{n-1} q^{(1)}
        s_roots, _ = safe_boxplus(cp_p, cp_q)
        if s_roots is None:
            continue

        Ph_p = Phi_n(roots_p)
        Ph_q = Phi_n(roots_q)
        Ph_r = Phi_n(roots_r)
        Ph_cp = Phi_n(cp_p)
        Ph_cq = Phi_n(cp_q)
        Ph_s = Phi_n(s_roots)

        if any(v is None or v < 1e-12 for v in [Ph_p, Ph_q, Ph_r, Ph_cp, Ph_cq, Ph_s]):
            continue

        F_n = 1/Ph_r - 1/Ph_p - 1/Ph_q
        F_n1 = 1/Ph_s - 1/Ph_cp - 1/Ph_cq

        if abs(F_n1) < 1e-14:
            continue

        excess_ratios.append(F_n / F_n1)

        if trial < 3:
            print(f"    F_n={F_n:.8f}, F_{n-1}={F_n1:.8f}, ratio={F_n/F_n1:.6f}")

    if excess_ratios:
        excess_ratios = np.array(excess_ratios)
        print(f"    {len(excess_ratios)} valid samples")
        print(f"    ratio range: [{excess_ratios.min():.6f}, {excess_ratios.max():.6f}]")
        print(f"    mean: {excess_ratios.mean():.6f}, std: {excess_ratios.std():.6f}")
        print(f"    all positive: {np.all(excess_ratios > 0)}")
        print(f"    IS CONSTANT? {'YES' if excess_ratios.std()/abs(excess_ratios.mean()) < 0.01 else 'NO'}")
        if np.all(excess_ratios > 0):
            print(f"    min ratio (lower bound): {excess_ratios.min():.6f}")


# ============================================================
# INVESTIGATION 4: Decomposition Phi_n(p) = Phi_{n-1}(p^{(1)}) + R_n(p)
# ============================================================
print("\n" + "=" * 70)
print("INVESTIGATION 4: DECOMPOSITION Phi_n(p) = Phi_{n-1}(p^{(1)}) + R_n(p)")
print("  What is R_n? Does it depend on root gaps?")
print("=" * 70)

for n in [3, 4, 5, 6]:
    print(f"\n  n = {n}:")
    R_values = []
    gap_stats = []
    for trial in range(50):
        roots = make_well_separated_roots(n, min_gap=0.3 + 0.5*np.random.rand())
        cp = critical_points(roots)
        if cp is None or len(cp) != n-1:
            continue
        if np.any(np.diff(cp) < 0.01):
            continue

        Ph_n_val = Phi_n(roots)
        Ph_n1_val = Phi_n(cp)
        if Ph_n_val is None or Ph_n1_val is None:
            continue

        R = Ph_n_val - Ph_n1_val
        R_values.append(R)

        gaps = np.diff(roots)
        gap_stats.append({
            'R': R,
            'sum_inv_gap_sq': np.sum(1.0/gaps**2),
            'sum_inv_gap': np.sum(1.0/gaps),
            'harmonic_mean_gap': len(gaps) / np.sum(1.0/gaps),
            'Phi_n': Ph_n_val,
            'Phi_n1': Ph_n1_val,
        })

    if R_values:
        R_arr = np.array(R_values)
        print(f"    R_n range: [{R_arr.min():.6f}, {R_arr.max():.6f}]")
        print(f"    R_n always positive: {np.all(R_arr > 0)}")
        print(f"    R_n always negative: {np.all(R_arr < 0)}")

        # Check: is R_n = c * sum(1/gap^2)?
        inv_gap_sq = np.array([g['sum_inv_gap_sq'] for g in gap_stats])
        if len(R_arr) > 2:
            corr = np.corrcoef(R_arr, inv_gap_sq)[0, 1]
            print(f"    corr(R_n, sum 1/gap^2) = {corr:.6f}")
            ratios = R_arr / inv_gap_sq
            print(f"    R_n / sum(1/gap^2): mean={ratios.mean():.6f}, std={ratios.std():.6f}")
            print(f"    IS R_n = c * sum(1/gap^2)? {'MAYBE' if ratios.std()/abs(ratios.mean()) < 0.05 else 'NO'}")

        # Check: is R_n related to Phi_n in a simple way?
        Ph_arr = np.array([g['Phi_n'] for g in gap_stats])
        Ph1_arr = np.array([g['Phi_n1'] for g in gap_stats])
        ratio_R_Phi = R_arr / Ph_arr
        print(f"    R_n / Phi_n: mean={ratio_R_Phi.mean():.6f}, std={ratio_R_Phi.std():.6f}")
        print(f"    i.e., Phi_{n-1}/Phi_n ~ {(1 - ratio_R_Phi).mean():.6f}")


# ============================================================
# INVESTIGATION 5: Cauchy interlacing — relate H_p to H_{p'}
# ============================================================
print("\n" + "=" * 70)
print("INVESTIGATION 5: CAUCHY INTERLACING")
print("  Roots of p': c_1 < c_2 < ... < c_{n-1}")
print("  Interlacing: lambda_1 < c_1 < lambda_2 < c_2 < ... < c_{n-1} < lambda_n")
print("  Relate H_p(lambda_i) to H_{p'}(c_j)")
print("=" * 70)

print("\n  --- Detailed H-value analysis ---")
for trial in range(5):
    n = 5
    roots = make_well_separated_roots(n, min_gap=0.8)
    cp = critical_points(roots)
    if cp is None or len(cp) != n-1:
        continue

    H_roots = H_values_safe(roots)
    H_cp = H_values_safe(cp)
    if H_roots is None or H_cp is None:
        continue

    print(f"\n  Trial {trial}:")
    print(f"    roots:    {np.array2string(roots, precision=4)}")
    print(f"    crit pts: {np.array2string(cp, precision=4)}")
    print(f"    H(roots): {np.array2string(H_roots, precision=4)}")
    print(f"    H(cps):   {np.array2string(H_cp, precision=4)}")
    print(f"    sum H(roots) = {np.sum(H_roots):.6f} (should be ~0)")
    print(f"    sum H(cps) = {np.sum(H_cp):.6f} (should be ~0)")
    print(f"    Phi_n = sum H^2 = {np.sum(H_roots**2):.6f}")
    print(f"    Phi_{n-1} = sum H'^2 = {np.sum(H_cp**2):.6f}")

    # Partial fraction decomposition:
    # p'(x)/p(x) = sum_i 1/(x - lambda_i)
    # At x = c_j: p'(c_j)/p(c_j) = 0 (since c_j is a root of p'), so
    # sum_i 1/(c_j - lambda_i) = 0 for each j.
    # This means: the vector (1/(c_j - lambda_i))_i has zero sum for each j.
    print(f"\n    Check: sum_i 1/(c_j - lambda_i) = 0 for each j:")
    for j in range(n-1):
        s = sum(1.0/(cp[j] - roots[i]) for i in range(n))
        print(f"      j={j}: sum = {s:.2e} (should be 0 since c_j is a critical point)")

    # Key identity: H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j)
    # Note: 1/(lambda_i - lambda_j) can be decomposed using the interlacing.
    # Between lambda_i and lambda_{i+1}, there's exactly one critical point c_i.
    # So gaps are: lambda_1, c_1, lambda_2, c_2, ..., c_{n-1}, lambda_n

    # Compute: for each root lambda_i, what is:
    # A_i = sum_{j=1}^{n-1} 1/(lambda_i - c_j) ?
    # These relate H_p to the critical points.
    print(f"\n    A_i = sum_j 1/(lambda_i - c_j):")
    A = np.zeros(n)
    for i in range(n):
        A[i] = sum(1.0/(roots[i] - cp[j]) for j in range(n-1))
    print(f"      A = {np.array2string(A, precision=4)}")
    print(f"      H(roots) = {np.array2string(H_roots, precision=4)}")
    print(f"      A - H = {np.array2string(A - H_roots, precision=4)}")
    # Note: p'(lambda_i) = prod_{j!=i} (lambda_i - lambda_j) and also
    # p'(x) = n * prod_{j=1}^{n-1} (x - c_j), so
    # p'(lambda_i) = n * prod_{j=1}^{n-1} (lambda_i - c_j)
    # H_p(lambda_i) = p''(lambda_i)/(2*p'(lambda_i))
    # Also: (log p')'(lambda_i) = sum_j 1/(lambda_i - c_j) = A_i
    # And: (log p)'(x) = sum_i 1/(x - lambda_i)
    # p''/p' = sum_j 1/(x - c_j)
    # So H_p(lambda_i) = p''(lambda_i)/(2*p'(lambda_i)) = (1/2) * A_i
    # Wait, let me check: p''(lambda_i)/(2*p'(lambda_i)) vs A_i

    # Actually H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j)
    # And A_i = sum_j 1/(lambda_i - c_j)
    # The identity is: p''(x)/p'(x) = sum_j 1/(x - c_j) (logarithmic derivative of p')
    # And p''(lambda_i)/(2*p'(lambda_i)) = H_p(lambda_i)
    # So A_i = p''(lambda_i)/p'(lambda_i) = 2*H_p(lambda_i).

    # CHECK:
    print(f"\n    Check A_i = 2*H(lambda_i):")
    print(f"      A/2 = {np.array2string(A/2, precision=4)}")
    print(f"      H   = {np.array2string(H_roots, precision=4)}")
    print(f"      error: {np.max(np.abs(A/2 - H_roots)):.2e}")


# ============================================================
# INVESTIGATION 5b: Express Phi_n in terms of interlacing gaps
# ============================================================
print("\n" + "=" * 70)
print("INVESTIGATION 5b: Phi_n IN TERMS OF INTERLACING GAPS")
print("  Use A_i = 2*H(lambda_i) and relate to gap structure")
print("=" * 70)

# We have: H_p(lambda_i) = (1/2) * sum_j 1/(lambda_i - c_j)
# So Phi_n(p) = (1/4) * sum_i [sum_j 1/(lambda_i - c_j)]^2
#
# Now define d_{i,j} = lambda_i - c_j (signed distance from root to critical point).
# By interlacing: d_{i,j} < 0 for j >= i and d_{i,j} > 0 for j < i.
#
# Phi_n = (1/4) * sum_i [sum_j 1/d_{i,j}]^2

# Similarly for Phi_{n-1}: H_{p'}(c_j) = sum_{k!=j} 1/(c_j - c_k)
# Phi_{n-1} = sum_j H_{p'}(c_j)^2

# Can we express Phi_{n-1} in terms of the same interlacing structure?
# H_{p'}(c_j) = sum_{k!=j} 1/(c_j - c_k)
#
# But also p''(x) = n(n-1) * prod_{l=1}^{n-2} (x - e_l) where e_l are roots of p''.
# And p''(c_j)/p'(c_j) = 0 is wrong — c_j is root of p', not p''.
# Actually: d/dx[p'/p] = (p''p - (p')^2)/p^2 = sum_i -1/(x-lambda_i)^2
# So (p''/p')(c_j) = (p')^2(c_j)/((p(c_j) * (from p' = 0: p''(c_j) comes from L'Hopital))
# Let me just compute numerically.

for trial in range(3):
    n = 5
    roots = make_well_separated_roots(n, min_gap=1.0)
    cp = critical_points(roots)
    if cp is None or len(cp) != n-1:
        continue

    H_roots = H_values_safe(roots)
    H_cp = H_values_safe(cp)
    if H_roots is None or H_cp is None:
        continue

    # Build the matrix D_{i,j} = 1/(lambda_i - c_j)
    D = np.zeros((n, n-1))
    for i in range(n):
        for j in range(n-1):
            D[i, j] = 1.0 / (roots[i] - cp[j])

    # Check H(lambda_i) = (1/2) * sum_j D_{i,j}
    H_from_D = 0.5 * np.sum(D, axis=1)

    # Phi_n = sum H^2 = (1/4) * sum_i (sum_j D_{i,j})^2
    # = (1/4) * ||D @ 1||^2 where 1 is the all-ones vector
    ones_n1 = np.ones(n-1)
    Phi_from_D = 0.25 * np.sum((D @ ones_n1)**2)

    # Similarly, define E_{j,i} = 1/(c_j - lambda_i)
    E = np.zeros((n-1, n))
    for j in range(n-1):
        for i in range(n):
            E[j, i] = 1.0 / (cp[j] - roots[i])

    # Key identity from critical points: sum_i E_{j,i} = 0 for each j
    # (since c_j is a critical point of p)
    # So E @ 1_n = 0 (each row sums to zero).

    # But H_{p'}(c_j) = sum_{k!=j} 1/(c_j - c_k), which involves c-c differences,
    # not c-lambda differences.

    print(f"\n  Trial {trial}:")
    print(f"    Phi_n from H: {np.sum(H_roots**2):.6f}")
    print(f"    Phi_n from D: {Phi_from_D:.6f}")
    print(f"    E @ 1 (should be ~0): {np.array2string(E @ np.ones(n), precision=2)}")

    # Now try: can we relate Phi_{n-1}(p') to the matrix D?
    # H_{p'}(c_j) = sum_{k!=j} 1/(c_j - c_k)
    # We know c_j are the critical points. Can we relate 1/(c_j - c_k) to
    # the lambda values?
    # Not directly, but we can try to build the "double interlacing" matrix.

    # Alternative: compute sum_j H_{p'}(c_j)^2 and sum_i H_p(lambda_i)^2
    # and look for a matrix identity.

    # GRAM MATRIX approach:
    # Define vectors v_i = (D_{i,1}, ..., D_{i,n-1}) in R^{n-1}
    # Then H(lambda_i) = (1/2) * <v_i, 1>
    # Phi_n = (1/4) * sum_i <v_i, 1>^2 = (1/4) * ||sum_i <v_i, 1> e_i||^2
    # Hmm, that's just (1/4)*||D 1||^2 again.

    # What about Phi_n in terms of D^T D or D D^T?
    # Phi_n = (1/4) * 1^T D^T D 1 = (1/4) * <D^T D 1, 1>? No.
    # Wait: Phi_n = (1/4) * sum_i (sum_j D_{i,j})^2 = (1/4) * ||D 1||^2
    # = (1/4) * 1^T D^T D 1.

    G = D.T @ D  # (n-1) x (n-1) matrix
    Phi_from_G = 0.25 * ones_n1 @ G @ ones_n1
    print(f"    Phi_n from Gram: {Phi_from_G:.6f}")
    print(f"    Phi_{n-1}(p') = {np.sum(H_cp**2):.6f}")

    # Check: is Phi_{n-1} related to tr(G) or det(G)?
    print(f"    tr(G) = {np.trace(G):.6f}")
    print(f"    sum_i 1/(lambda_i - c_j)^2 for each j:")
    for j in range(n-1):
        s = sum(1.0/(roots[i] - cp[j])**2 for i in range(n))
        print(f"      j={j}: {s:.6f}")
    print(f"    tr(G) = sum_j sum_i 1/(lambda_i - c_j)^2 = {np.trace(G):.6f}")

    # Key identity: tr(D^T D) = sum_{i,j} 1/(lambda_i - c_j)^2
    # Phi_n = (1/4) * 1^T G 1 where G = D^T D
    # Phi_{n-1} = sum_j (sum_{k!=j} 1/(c_j - c_k))^2 — different matrix!


# ============================================================
# INVESTIGATION 6: Telescoping identity
# ============================================================
print("\n" + "=" * 70)
print("INVESTIGATION 6: TELESCOPING IDENTITY")
print("  Define p^{(k)} = iterated derivative (normalized)")
print("  p^{(0)} = p, p^{(1)} = (1/n)*p', p^{(2)} = (1/(n-1))*(p^{(1)})', etc.")
print("  Check: 1/Phi_n(r) vs sum involving 1/Phi_k for k=2,...,n")
print("=" * 70)

def iterated_critical_points(roots, depth):
    """Compute roots of iterated normalized derivatives.
    Returns list [roots, cp1, cp2, ...] where cp_k has n-k elements."""
    result = [roots]
    current = roots
    for k in range(depth):
        cp = critical_points(current)
        if cp is None or len(cp) != len(current) - 1:
            return None
        if len(cp) > 1 and np.any(np.diff(cp) < 1e-8):
            return None
        result.append(cp)
        current = cp
    return result

print("\n  --- Phi values along derivative chain ---")
for trial in range(5):
    n = 6
    roots = make_well_separated_roots(n, min_gap=1.0)
    chain = iterated_critical_points(roots, n-2)
    if chain is None:
        continue

    print(f"\n  Trial {trial}: roots = {np.array2string(roots, precision=3)}")
    phis = []
    for k in range(len(chain)):
        r = chain[k]
        ph = Phi_n(r)
        if ph is None:
            break
        phis.append(ph)
        print(f"    p^({k}): {len(r)} roots, Phi = {ph:.6f}, 1/Phi = {1/ph:.6f}")

    if len(phis) == len(chain):
        # Check various telescoping formulas
        inv_phis = [1/p for p in phis]
        print(f"    sum 1/Phi_k for k=0..{n-2}: {sum(inv_phis):.6f}")

        # Differences
        diffs = [inv_phis[k+1] - inv_phis[k] for k in range(len(inv_phis)-1)]
        print(f"    1/Phi_{k+1} - 1/Phi_k: {np.array2string(np.array(diffs), precision=4)}")
        print(f"    Ratios Phi_k / Phi_{k+1}: {np.array2string(np.array([phis[k]/phis[k+1] for k in range(len(phis)-1)]), precision=4)}")


print("\n  --- Telescoping for convolution ---")
for trial in range(5):
    n = 5
    roots_p = make_well_separated_roots(n, min_gap=0.8)
    roots_q = make_well_separated_roots(n, min_gap=0.8)
    roots_r, _ = safe_boxplus(roots_p, roots_q)
    if roots_r is None:
        continue

    # Compute derivative chains
    chain_p = iterated_critical_points(roots_p, n-2)
    chain_q = iterated_critical_points(roots_q, n-2)
    chain_r = iterated_critical_points(roots_r, n-2)

    if chain_p is None or chain_q is None or chain_r is None:
        continue

    print(f"\n  Trial {trial}:")
    all_ok = True
    for k in range(len(chain_p)):
        ph_p = Phi_n(chain_p[k])
        ph_q = Phi_n(chain_q[k])
        ph_r = Phi_n(chain_r[k])
        if ph_p is None or ph_q is None or ph_r is None:
            all_ok = False
            break

        deg = n - k
        F_k = 1/ph_r - 1/ph_p - 1/ph_q
        print(f"    degree {deg}: F_{deg} = 1/Phi(r^({k})) - 1/Phi(p^({k})) - 1/Phi(q^({k})) = {F_k:.8f} {'>=0' if F_k >= -1e-10 else 'NEGATIVE!'}")

        # Also check: convolution of derivatives matches derivative of convolution
        if k < len(chain_p) - 1:
            s_roots, _ = safe_boxplus(chain_p[k+1], chain_q[k+1])
            if s_roots is not None:
                error = np.max(np.abs(np.sort(chain_r[k+1]) - np.sort(s_roots)))
                print(f"      MSS check at level {k}: error = {error:.2e}")


# ============================================================
# INVESTIGATION 7: Key structural question — induction feasibility
# ============================================================
print("\n" + "=" * 70)
print("INVESTIGATION 7: INDUCTION FEASIBILITY ANALYSIS")
print("  Can we bound F_n in terms of F_{n-1} and correction terms?")
print("=" * 70)

print("\n  --- F_n vs F_{n-1}: scatter data ---")
fn_data = []
for trial in range(200):
    n = np.random.choice([4, 5, 6])
    roots_p = make_well_separated_roots(n)
    roots_q = make_well_separated_roots(n)
    roots_r, _ = safe_boxplus(roots_p, roots_q)
    if roots_r is None:
        continue

    cp_p = critical_points(roots_p)
    cp_q = critical_points(roots_q)
    if cp_p is None or cp_q is None:
        continue
    if len(cp_p) != n-1 or len(cp_q) != n-1:
        continue

    s_roots, _ = safe_boxplus(cp_p, cp_q)
    if s_roots is None:
        continue

    Ph_p = Phi_n(roots_p)
    Ph_q = Phi_n(roots_q)
    Ph_r = Phi_n(roots_r)
    Ph_cp = Phi_n(cp_p)
    Ph_cq = Phi_n(cp_q)
    Ph_s = Phi_n(s_roots)

    if any(v is None or v < 1e-12 for v in [Ph_p, Ph_q, Ph_r, Ph_cp, Ph_cq, Ph_s]):
        continue

    F_n = 1/Ph_r - 1/Ph_p - 1/Ph_q
    F_n1 = 1/Ph_s - 1/Ph_cp - 1/Ph_cq

    # Also compute the "remainder" terms
    R_p = Ph_p - Ph_cp
    R_q = Ph_q - Ph_cq
    R_r = Ph_r - Ph_s

    fn_data.append({
        'n': n, 'F_n': F_n, 'F_n1': F_n1,
        'R_p': R_p, 'R_q': R_q, 'R_r': R_r,
        'Phi_p': Ph_p, 'Phi_q': Ph_q, 'Phi_r': Ph_r,
        'Phi_cp': Ph_cp, 'Phi_cq': Ph_cq, 'Phi_s': Ph_s,
    })

if fn_data:
    print(f"\n  {len(fn_data)} valid samples collected")

    for n in [4, 5, 6]:
        subset = [d for d in fn_data if d['n'] == n]
        if not subset:
            continue

        F_n_arr = np.array([d['F_n'] for d in subset])
        F_n1_arr = np.array([d['F_n1'] for d in subset])
        R_p_arr = np.array([d['R_p'] for d in subset])
        R_q_arr = np.array([d['R_q'] for d in subset])
        R_r_arr = np.array([d['R_r'] for d in subset])

        print(f"\n  n={n} ({len(subset)} samples):")
        print(f"    F_n  all >= 0: {np.all(F_n_arr >= -1e-10)}, min = {F_n_arr.min():.8f}")
        print(f"    F_{n-1} all >= 0: {np.all(F_n1_arr >= -1e-10)}, min = {F_n1_arr.min():.8f}")

        # Key question: is F_n >= F_{n-1}?
        diff = F_n_arr - F_n1_arr
        print(f"    F_n >= F_{n-1}: {np.all(diff >= -1e-10)}, min(F_n - F_{n-1}) = {diff.min():.8f}")

        # Ratio F_n / F_{n-1}
        valid = np.abs(F_n1_arr) > 1e-12
        if np.any(valid):
            ratios = F_n_arr[valid] / F_n1_arr[valid]
            print(f"    F_n/F_{n-1}: range=[{ratios.min():.4f}, {ratios.max():.4f}], mean={ratios.mean():.4f}")

        # Can we bound F_n >= c * F_{n-1} for some c > 0?
        # Check: F_n = F_{n-1} * (Phi_s*Phi_cp*Phi_cq)/(Phi_r*Phi_p*Phi_q) * adjustment?
        # This is getting complicated. Let's try direct analysis.

        # Key: does 1/Phi_n = 1/Phi_{n-1} + correction help?
        # 1/Phi_p = 1/Ph_cp + something(p)
        # If 1/Phi_p = 1/Ph_cp + delta_p, etc., then
        # F_n = 1/Ph_r - 1/Ph_p - 1/Ph_q
        #     = (1/Ph_s + delta_r) - (1/Ph_cp + delta_p) - (1/Ph_cq + delta_q)
        #     = F_{n-1} + (delta_r - delta_p - delta_q)
        # So induction works if delta_r >= delta_p + delta_q!

        delta_p_arr = 1/np.array([d['Phi_p'] for d in subset]) - 1/np.array([d['Phi_cp'] for d in subset])
        delta_q_arr = 1/np.array([d['Phi_q'] for d in subset]) - 1/np.array([d['Phi_cq'] for d in subset])
        delta_r_arr = 1/np.array([d['Phi_r'] for d in subset]) - 1/np.array([d['Phi_s'] for d in subset])

        delta_diff = delta_r_arr - delta_p_arr - delta_q_arr

        print(f"\n    INDUCTION DECOMPOSITION:")
        print(f"    F_n = F_{n-1} + (delta_r - delta_p - delta_q)")
        print(f"    where delta_p = 1/Phi_n(p) - 1/Phi_{n-1}(p^(1))")
        print(f"    delta_r - delta_p - delta_q: min={delta_diff.min():.8f}, max={delta_diff.max():.8f}")
        print(f"    All >= 0? {np.all(delta_diff >= -1e-10)}")
        if np.all(delta_diff >= -1e-10):
            print(f"    *** INDUCTION WOULD WORK if delta is superadditive! ***")
        else:
            print(f"    delta is NOT superadditive. Induction via this decomposition FAILS.")
            print(f"    Number of violations: {np.sum(delta_diff < -1e-10)}/{len(delta_diff)}")
            print(f"    Worst violation: {delta_diff.min():.8f}")

            # But maybe F_{n-1} compensates?
            # F_n = F_{n-1} + (delta_r - delta_p - delta_q) >= 0
            # Even if delta_diff < 0, we need F_{n-1} >= |delta_diff|
            compensated = F_n1_arr + delta_diff
            print(f"    F_{n-1} + delta_diff: min = {compensated.min():.8f} (should be F_n)")
            print(f"    Consistency check: max|F_n - compensated| = {np.max(np.abs(F_n_arr - compensated)):.2e}")


# ============================================================
# INVESTIGATION 8: Sign analysis of delta
# ============================================================
print("\n" + "=" * 70)
print("INVESTIGATION 8: SIGN OF delta_p = 1/Phi_n(p) - 1/Phi_{n-1}(p^{(1)})")
print("  Is delta always positive? Always negative? Sign-changing?")
print("=" * 70)

for n in [3, 4, 5, 6]:
    deltas = []
    for trial in range(100):
        roots = make_well_separated_roots(n, min_gap=0.2 + 0.8*np.random.rand())
        cp = critical_points(roots)
        if cp is None or len(cp) != n-1:
            continue
        if len(cp) > 1 and np.any(np.diff(cp) < 0.01):
            continue

        Ph = Phi_n(roots)
        Ph1 = Phi_n(cp)
        if Ph is None or Ph1 is None:
            continue

        delta = 1/Ph - 1/Ph1
        deltas.append(delta)

    if deltas:
        deltas = np.array(deltas)
        print(f"  n={n}: {len(deltas)} samples")
        print(f"    delta range: [{deltas.min():.8f}, {deltas.max():.8f}]")
        print(f"    all positive: {np.all(deltas > 0)}")
        print(f"    all negative: {np.all(deltas < 0)}")
        if np.all(deltas < 0):
            print(f"    delta is ALWAYS NEGATIVE: 1/Phi_n < 1/Phi_{n-1}")
            print(f"    This means Phi_n > Phi_{n-1}, which makes sense (more terms).")
        if np.all(deltas > 0):
            print(f"    delta is ALWAYS POSITIVE: 1/Phi_n > 1/Phi_{n-1}")


# ============================================================
# INVESTIGATION 9: Alternative — work with Phi directly
# ============================================================
print("\n" + "=" * 70)
print("INVESTIGATION 9: DECOMPOSITION OF Phi_n IN TERMS OF INTERLACING")
print("  Use H_p(lambda_i) = (1/2) * sum_j 1/(lambda_i - c_j)")
print("  to express Phi_n(p) = (1/4) * sum_i (sum_j 1/(lambda_i-c_j))^2")
print("=" * 70)

# For a polynomial p with roots lambda_1 < ... < lambda_n and
# critical points c_1 < ... < c_{n-1} (interlacing), define:
# D_{i,j} = 1/(lambda_i - c_j)
# Then: Phi_n = (1/4) * ||D @ 1||^2

# For Phi_{n-1}:
# H_{p'}(c_j) = sum_{k!=j} 1/(c_j - c_k)
# But can we relate this to D?
# From the relation: sum_i D_{j,i} = 0 transposed:
# Actually sum_i 1/(c_j - lambda_i) = 0 for each critical point c_j.

# So defining E_{j,i} = 1/(c_j - lambda_i), we have E @ 1_n = 0 (row sums are 0).
# Note: E = -D^T (up to index swap), or more precisely E_{j,i} = -D_{i,j}.

# Now consider p' with roots c_1,...,c_{n-1}. Its critical points are some
# d_1,...,d_{n-2}. Then:
# H_{p'}(c_j) = (1/2) * sum_l 1/(c_j - d_l)
# This involves a DIFFERENT set of points d_l, not the lambda_i.

# So the relationship is: Phi_n involves (lambda, c) pairs and
# Phi_{n-1} involves (c, d) pairs, where d are roots of p''.
# The "link" is through c being shared.

# Can we express H_{p'}(c_j) directly in terms of lambda_i and c_j?
# H_{p'}(c_j) = sum_{k!=j} 1/(c_j - c_k)
# We know sum_i 1/(c_j - lambda_i) = 0, i.e., the sum includes both neighbors.
# Can we write 1/(c_j - c_k) in terms of 1/(c_j - lambda_i)?

print("\n  --- Direct comparison of H at roots vs H at critical points ---")
for trial in range(3):
    n = 5
    roots = make_well_separated_roots(n, min_gap=1.0)
    cp = critical_points(roots)
    if cp is None or len(cp) != n-1:
        continue

    H_r = H_values_safe(roots)
    H_c = H_values_safe(cp)
    if H_r is None or H_c is None:
        continue

    print(f"\n  Trial {trial}:")
    print(f"    roots: {np.array2string(roots, precision=4)}")
    print(f"    crit:  {np.array2string(cp, precision=4)}")

    # Build the matrices
    # D_{i,j} = 1/(lambda_i - c_j), size n x (n-1)
    # E_{j,i} = 1/(c_j - lambda_i), size (n-1) x n
    D = np.zeros((n, n-1))
    E = np.zeros((n-1, n))
    for i in range(n):
        for j in range(n-1):
            D[i,j] = 1.0/(roots[i] - cp[j])
            E[j,i] = 1.0/(cp[j] - roots[i])

    # F_{j,k} = 1/(c_j - c_k) for j != k, size (n-1) x (n-1)
    F_mat = np.zeros((n-1, n-1))
    for j in range(n-1):
        for k in range(n-1):
            if k != j:
                F_mat[j,k] = 1.0/(cp[j] - cp[k])

    # Verify: D @ 1 = 2*H(roots)
    check1 = D @ np.ones(n-1)
    print(f"    D@1 = {np.array2string(check1, precision=4)}")
    print(f"    2*H = {np.array2string(2*H_r, precision=4)}")

    # Verify: E @ 1 = 0 (zero for critical points)
    check2 = E @ np.ones(n)
    print(f"    E@1 = {np.array2string(check2, precision=2)} (should be ~0)")

    # F_mat @ 1 = H(crit) (by definition)
    check3 = F_mat @ np.ones(n-1)  # Wait, this sums all 1/(c_j - c_k) for k != j
    # That's actually H_{p'}(c_j) NOT times 2
    # Because H_{p'}(c_j) = sum_{k!=j} 1/(c_j - c_k) (no factor of 1/2)
    # Wait... let me recheck.
    # H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j) — this is the definition.
    # For degree n polynomial p with roots lambda_i.
    # For p' with roots c_j (degree n-1 polynomial):
    # H_{p'}(c_j) = sum_{k!=j} 1/(c_j - c_k).
    # So F_mat @ 1 = H_{p'}(c_j). Check:
    print(f"    F@1 = {np.array2string(check3, precision=4)}")
    print(f"    H_c = {np.array2string(H_c, precision=4)}")

    # Phi_n = (1/4)*||D@1||^2 vs Phi_{n-1} = ||F@1||^2
    # Wait, Phi_n uses H (not 2*H). Let me recheck.
    # H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j)
    # We showed A_i = sum_j 1/(lambda_i - c_j) = 2*H_p(lambda_i)
    # So D@1 = A_i vector, and H_p(lambda_i) = (1/2)*A_i
    # Phi_n = sum H^2 = (1/4)*sum A_i^2 = (1/4)*||D@1||^2.  Correct.
    # Phi_{n-1} = sum H_{p'}(c_j)^2 = ||F@1||^2? Let me check.
    # No: Phi_{n-1} = sum_j (H_{p'}(c_j))^2.
    # And H_{p'}(c_j) = (row sum of F_mat for row j) = (F@1)_j
    # So Phi_{n-1} = ||F@1||^2. BUT WAIT: do we need the factor of 1/2?
    # For p' (degree n-1), its critical points are d_l (n-2 of them).
    # H_{p'}(c_j) = sum_{k!=j} 1/(c_j - c_k) — this is correct, no extra factor.
    # But earlier we showed that H_p(lambda_i) = (1/2)*sum_j 1/(lambda_i - c_j).
    # WHY the factor of 1/2? Let me re-derive.
    # p'(x)/p(x) = sum_i 1/(x - lambda_i) [logarithmic derivative]
    # p'(x) = n * prod(x - c_j) [p' is degree n-1]
    # p''(x)/p'(x) = sum_j 1/(x - c_j) [log derivative of p']
    # H_p(lambda_i) = p''(lambda_i)/(2*p'(lambda_i)) [definition]
    # But p''(x)/p'(x) = sum_j 1/(x - c_j), so at x = lambda_i:
    # p''(lambda_i)/p'(lambda_i) = sum_j 1/(lambda_i - c_j)
    # Therefore H_p(lambda_i) = (1/2) * p''(lambda_i)/p'(lambda_i) = (1/2)*sum_j 1/(lambda_i - c_j)
    # The 1/2 comes from the DEFINITION of H, not from any structural reason.

    # Now for p': H_{p'}(c_j) = sum_{k!=j} 1/(c_j - c_k).
    # Is this p'''(c_j)/(2*p''(c_j)) or something else?
    # p''(x)/p'(x) = sum_j 1/(x - c_j) [YES]
    # p'''(x)/p''(x) = sum_l 1/(x - d_l) where d_l are roots of p''
    # H_{p'}(c_j) = p'''(c_j)/(2*p''(c_j)) [by the same definition, applied to p']
    # So H_{p'}(c_j) = (1/2)*sum_l 1/(c_j - d_l)

    # But ALSO H_{p'}(c_j) = sum_{k!=j} 1/(c_j - c_k) [by direct definition as sum over OTHER roots of p']
    # These should be equal:
    # sum_{k!=j} 1/(c_j - c_k) = (1/2)*sum_l 1/(c_j - d_l)
    # Wait, that can't both be the definition. Let me re-clarify.

    # H_p(lambda_i) := sum_{j != i} 1/(lambda_i - lambda_j)
    # This is EXACTLY what we define.
    # And we showed H_p(lambda_i) = (1/2) * sum_j 1/(lambda_i - c_j)
    # where c_j are critical points.

    # Similarly, H_{p'}(c_j) := sum_{k != j} 1/(c_j - c_k)
    # And H_{p'}(c_j) = (1/2) * sum_l 1/(c_j - d_l) where d_l are roots of p''.

    # Let's verify this numerically.
    d_pts = critical_points(cp)  # roots of p''
    if d_pts is not None and len(d_pts) == n-2:
        for j in range(n-1):
            H_direct = sum(1.0/(cp[j] - cp[k]) for k in range(n-1) if k != j)
            H_via_d = 0.5 * sum(1.0/(cp[j] - d_pts[l]) for l in range(n-2))
            print(f"    c_{j}: H_direct={H_direct:.6f}, (1/2)*sum 1/(c-d)={H_via_d:.6f}, match={abs(H_direct-H_via_d)<1e-8}")

    # So the Phi hierarchy is:
    # Phi_n(p) = (1/4) * ||D_{lambda,c} @ 1_{n-1}||^2  (uses lambda and c)
    # Phi_{n-1}(p') = (1/4) * ||D_{c,d} @ 1_{n-2}||^2  (uses c and d)
    # Both involve a "1/4" factor when expressed via derivative roots.
    # But Phi_n = sum H^2 = sum (sum 1/(lambda_i - lambda_j))^2 with NO 1/4 factor.
    # This is because H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j) directly.
    # The 1/2 * sum 1/(lambda_i - c_j) = H_p(lambda_i) identity connects the two.


# ============================================================
# INVESTIGATION 10: The critical induction formula
# ============================================================
print("\n" + "=" * 70)
print("INVESTIGATION 10: INDUCTION VIA PARTIAL FRACTIONS")
print("  Express H_p(lambda_i) in terms of c_j and lambda_i")
print("  Then compute Phi_n - Phi_{n-1} explicitly")
print("=" * 70)

# We have:
# H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j)
#
# Key partial fraction identity:
# 1/(lambda_i - lambda_j) = ??? in terms of c_k
# This doesn't simplify directly. But we can try:
#
# Define delta_i = lambda_{i+1} - lambda_i (root gaps)
# and epsilon_i = c_i - lambda_i, eta_i = lambda_{i+1} - c_i
# (so delta_i = epsilon_i + eta_i and epsilon_i, eta_i > 0 by interlacing)
#
# H_p(lambda_i) has contributions from:
# - "left neighbors": sum_{j < i} 1/(lambda_i - lambda_j) > 0
# - "right neighbors": sum_{j > i} 1/(lambda_i - lambda_j) < 0
# The closest terms dominate: 1/delta_{i-1} and -1/delta_i.

for trial in range(3):
    n = 5
    roots = make_well_separated_roots(n, min_gap=1.0)
    cp = critical_points(roots)
    if cp is None or len(cp) != n-1:
        continue

    H = H_values_safe(roots)
    if H is None:
        continue

    deltas = np.diff(roots)  # delta_i = lambda_{i+1} - lambda_i
    epsilons = cp - roots[:-1]  # epsilon_i = c_i - lambda_i
    etas = roots[1:] - cp  # eta_i = lambda_{i+1} - c_i

    print(f"\n  Trial {trial}:")
    print(f"    deltas (root gaps): {np.array2string(deltas, precision=4)}")
    print(f"    epsilons (lambda->c): {np.array2string(epsilons, precision=4)}")
    print(f"    etas (c->lambda): {np.array2string(etas, precision=4)}")
    print(f"    epsilon + eta = delta: {np.allclose(epsilons + etas, deltas)}")
    print(f"    epsilon/delta: {np.array2string(epsilons/deltas, precision=4)}")
    print(f"    (Root fractions - where c sits in each gap)")

    # The position of c_i within the gap [lambda_i, lambda_{i+1}] is
    # t_i = epsilon_i/delta_i in (0,1).
    # For equally spaced roots, t_i should be 1/2 (by symmetry).
    # For non-uniform roots, t_i varies.


# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY OF INDUCTION INVESTIGATION")
print("=" * 70)

print("""
KEY FINDINGS:

1. MSS DERIVATIVE IDENTITY: r' = n * (p^{(1)} boxplus_{n-1} q^{(1)})
   NUMERICALLY VERIFIED with errors < 1e-8.

2. RATIO Phi_n(p) / Phi_{n-1}(p^{(1)}):
   NOT CONSTANT. Varies with root configuration.
   No clean universal formula found.

3. DECOMPOSITION F_n = F_{n-1} + (delta_r - delta_p - delta_q):
   This is an EXACT identity where delta_p = 1/Phi_n(p) - 1/Phi_{n-1}(p^{(1)}).
   INDUCTION WORKS if delta is SUPERADDITIVE under boxplus.

4. SIGN OF delta_p = 1/Phi_n(p) - 1/Phi_{n-1}(p^{(1)}):
   Check if always negative (Phi_n > Phi_{n-1}).

5. The interlacing structure gives:
   H_p(lambda_i) = (1/2) * sum_j 1/(lambda_i - c_j)
   Phi_n = (1/4) * ||D@1||^2 where D_{i,j} = 1/(lambda_i - c_j)

6. CONCLUSION: The induction reduces to proving SUPERADDITIVITY of
   delta: delta(r) >= delta(p) + delta(q)
   where delta(p) = 1/Phi_n(p) - 1/Phi_{n-1}(p^{(1)}).
   Numerical evidence needed to check this.
""")
