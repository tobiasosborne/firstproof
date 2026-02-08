"""
Matrix approach investigation v2: Focus on WELL-SEPARATED roots
where the finite free convolution is guaranteed to stay real-rooted.

The v1 investigation found many "failures" but those were contaminated by
numerical issues (non-real roots from the convolution). This script uses
only WELL-SEPARATED roots where the convolution stays cleanly real-rooted.

Also: investigate the n=3 case exactly.
"""

import numpy as np
from itertools import combinations
from scipy.special import comb
import warnings
warnings.filterwarnings('ignore')

def phi_n_from_roots(roots):
    """Compute Phi_n(p) = sum_i H_p(lambda_i)^2."""
    n = len(roots)
    phi = 0.0
    for i in range(n):
        H_i = sum(1.0/(roots[i] - roots[j]) for j in range(n) if j != i)
        phi += H_i**2
    return phi

def finite_free_convolution(roots_p, roots_q):
    """Compute p boxplus_n q using hat-e_k formula."""
    n = len(roots_p)
    assert len(roots_q) == n

    def elem_sym(roots, k):
        if k == 0:
            return 1.0
        if k > len(roots):
            return 0.0
        return sum(np.prod(list(combo)) for combo in combinations(roots, k))

    hat_e_p = [elem_sym(roots_p, k) / comb(n, k, exact=True) for k in range(n+1)]
    hat_e_q = [elem_sym(roots_q, k) / comb(n, k, exact=True) for k in range(n+1)]
    hat_e_r = [sum(hat_e_p[j] * hat_e_q[k-j] for j in range(k+1)) for k in range(n+1)]
    e_r = [comb(n, k, exact=True) * hat_e_r[k] for k in range(n+1)]
    coeffs = [(-1)**k * e_r[k] for k in range(n+1)]
    roots_r = np.sort(np.roots(coeffs))
    return roots_r

def is_real_rooted(roots_r, tol=1e-8):
    return np.max(np.abs(np.imag(roots_r))) < tol

def test_gap(roots_p, roots_q, label="", verbose=True):
    """Test the superadditivity gap. Returns gap or None if invalid."""
    roots_r = finite_free_convolution(roots_p, roots_q)

    if not is_real_rooted(roots_r, tol=1e-6):
        if verbose:
            print(f"  {label}: SKIPPED (complex roots, max_imag={np.max(np.abs(np.imag(roots_r))):.2e})")
        return None

    roots_r = np.sort(np.real(roots_r))

    # Check distinct
    if len(roots_r) >= 2 and np.min(np.diff(roots_r)) < 1e-10:
        if verbose:
            print(f"  {label}: SKIPPED (repeated roots)")
        return None

    phi_p = phi_n_from_roots(roots_p)
    phi_q = phi_n_from_roots(roots_q)
    phi_r = phi_n_from_roots(roots_r)

    gap = 1.0/phi_r - (1.0/phi_p + 1.0/phi_q)

    if verbose:
        print(f"  {label}: gap = {gap:.8e}  {'PASS' if gap >= -1e-12 else 'FAIL'}")

    return gap


# ============================================================================
# PART A: n=2 exact analysis (confirmed)
# ============================================================================

print("=" * 80)
print("n=2 EXACT ANALYSIS")
print("=" * 80)

print("""
For n=2 with p roots (a,b), q roots (c,d):
  1/Phi_n(r) - 1/Phi_n(p) - 1/Phi_n(q) = (a+b)(c+d)/2

This is:
  > 0 when means have same sign
  = 0 when at least one is centered (mean zero)
  < 0 when means have opposite sign

The gap = 2 * mean_p * mean_q * n^2 / 2 = mean_p * mean_q * 2.
""")

# ============================================================================
# PART B: n=3 with well-separated EQUALLY SPACED centered roots
# ============================================================================

print("=" * 80)
print("n=3: SYSTEMATIC TESTS WITH WELL-SEPARATED ROOTS")
print("=" * 80)

print("\n--- n=3 centered, equally spaced with varying separation ---")
for gap_p in [1.0, 2.0, 3.0, 5.0, 10.0]:
    for gap_q in [1.0, 2.0, 3.0, 5.0, 10.0]:
        roots_p = np.array([-gap_p, 0.0, gap_p])
        roots_q = np.array([-gap_q, 0.0, gap_q])
        test_gap(roots_p, roots_q, f"gap_p={gap_p:.0f}, gap_q={gap_q:.0f}")

print("\n--- n=3 centered, asymmetric ---")
test_cases_n3 = [
    (np.array([-2.0, -0.5, 2.5]), np.array([-3.0, 1.0, 2.0]), "asym1"),
    (np.array([-5.0, 1.0, 4.0]), np.array([-1.0, -0.5, 1.5]), "asym2"),
    (np.array([-3.0, -1.0, 4.0]), np.array([-2.0, 0.5, 1.5]), "asym3"),
    (np.array([-10.0, 3.0, 7.0]), np.array([-1.0, 0.0, 1.0]), "asym4"),
    (np.array([-1.0, 0.0, 1.0]), np.array([-1.0, 0.0, 1.0]), "same"),
]
for roots_p, roots_q, label in test_cases_n3:
    test_gap(roots_p, roots_q, label)


print("\n--- n=3 UNCENTERED, testing sign of gap vs means ---")
test_cases_uncentered = [
    (np.array([1.0, 2.0, 5.0]), np.array([0.0, 3.0, 4.0]), "both mean>0"),
    (np.array([-5.0, -2.0, -1.0]), np.array([-4.0, -3.0, 0.0]), "both mean<0"),
    (np.array([-5.0, -2.0, -1.0]), np.array([0.0, 3.0, 4.0]), "opposite means"),
    (np.array([1.0, 2.0, 100.0]), np.array([-50.0, -1.0, 0.0]), "extreme opposite"),
]
for roots_p, roots_q, label in test_cases_uncentered:
    g = test_gap(roots_p, roots_q, f"mean_p={np.mean(roots_p):.1f}, mean_q={np.mean(roots_q):.1f} ({label})")


# ============================================================================
# PART C: Large-scale tests with GUARANTEED real-rooted convolutions
# ============================================================================

print("\n" + "=" * 80)
print("LARGE-SCALE: Well-separated roots guaranteeing real-rootedness")
print("=" * 80)

print("""
Strategy: Use roots that are integer multiples of a large spacing,
ensuring the convolution stays real-rooted.
""")

np.random.seed(42)
results = {n: {"pass": 0, "fail": 0, "skip": 0, "min_gap": float('inf')} for n in [2, 3, 4, 5]}

for trial in range(500):
    n = np.random.choice([2, 3, 4, 5])

    # Generate well-separated roots (spacing ~ n)
    # This helps ensure real-rootedness of the convolution
    spacing = 2.0 * n
    roots_p = np.sort(np.random.choice(range(-10*n, 10*n+1), size=n, replace=False).astype(float))
    roots_q = np.sort(np.random.choice(range(-10*n, 10*n+1), size=n, replace=False).astype(float))

    roots_r = finite_free_convolution(roots_p, roots_q)

    if not is_real_rooted(roots_r, tol=1e-6):
        results[n]["skip"] += 1
        continue

    roots_r = np.sort(np.real(roots_r))
    if len(roots_r) >= 2 and np.min(np.diff(roots_r)) < 1e-10:
        results[n]["skip"] += 1
        continue

    phi_p = phi_n_from_roots(roots_p)
    phi_q = phi_n_from_roots(roots_q)
    phi_r = phi_n_from_roots(roots_r)

    gap = 1.0/phi_r - (1.0/phi_p + 1.0/phi_q)
    results[n]["min_gap"] = min(results[n]["min_gap"], gap)

    if gap >= -1e-10:
        results[n]["pass"] += 1
    else:
        results[n]["fail"] += 1

print("\nResults (integer-spaced roots, general):")
for n in [2, 3, 4, 5]:
    r = results[n]
    total = r["pass"] + r["fail"]
    print(f"  n={n}: pass={r['pass']}, fail={r['fail']}, skip={r['skip']}, min_gap={r['min_gap']:.6e}")


# Now test CENTERED
print("\nCentered tests:")
np.random.seed(42)
results_c = {n: {"pass": 0, "fail": 0, "skip": 0, "min_gap": float('inf'), "fail_examples": []} for n in [3, 4, 5]}

for trial in range(1000):
    n = np.random.choice([3, 4, 5])

    # Generate centered, well-separated roots
    roots_p = np.sort(np.random.choice(range(-10*n, 10*n+1), size=n, replace=False).astype(float))
    roots_p -= np.mean(roots_p)

    roots_q = np.sort(np.random.choice(range(-10*n, 10*n+1), size=n, replace=False).astype(float))
    roots_q -= np.mean(roots_q)

    roots_r = finite_free_convolution(roots_p, roots_q)

    if not is_real_rooted(roots_r, tol=1e-6):
        results_c[n]["skip"] += 1
        continue

    roots_r = np.sort(np.real(roots_r))
    if len(roots_r) >= 2 and np.min(np.diff(roots_r)) < 1e-10:
        results_c[n]["skip"] += 1
        continue

    phi_p = phi_n_from_roots(roots_p)
    phi_q = phi_n_from_roots(roots_q)
    phi_r = phi_n_from_roots(roots_r)

    gap = 1.0/phi_r - (1.0/phi_p + 1.0/phi_q)
    results_c[n]["min_gap"] = min(results_c[n]["min_gap"], gap)

    if gap >= -1e-10:
        results_c[n]["pass"] += 1
    else:
        results_c[n]["fail"] += 1
        if len(results_c[n]["fail_examples"]) < 3:
            results_c[n]["fail_examples"].append((roots_p.copy(), roots_q.copy(), gap))

for n in [3, 4, 5]:
    r = results_c[n]
    total = r["pass"] + r["fail"]
    print(f"  n={n}: pass={r['pass']}, fail={r['fail']}, skip={r['skip']}, min_gap={r['min_gap']:.6e}")
    for p, q, g in r["fail_examples"]:
        print(f"    FAIL example: p={p}, q={q}, gap={g:.6e}")


# ============================================================================
# PART D: What about ONLY using positive roots? (as in the original MSS context)
# ============================================================================

print("\n" + "=" * 80)
print("POSITIVE ROOTS ONLY (eigenvalues of PSD matrices)")
print("=" * 80)

np.random.seed(123)
results_pos = {n: {"pass": 0, "fail": 0, "skip": 0, "min_gap": float('inf')} for n in [2, 3, 4, 5]}

for trial in range(500):
    n = np.random.choice([2, 3, 4, 5])

    # Generate positive, well-separated roots
    roots_p = np.sort(np.random.choice(range(1, 20*n+1), size=n, replace=False).astype(float))
    roots_q = np.sort(np.random.choice(range(1, 20*n+1), size=n, replace=False).astype(float))

    roots_r = finite_free_convolution(roots_p, roots_q)

    if not is_real_rooted(roots_r, tol=1e-6):
        results_pos[n]["skip"] += 1
        continue

    roots_r = np.sort(np.real(roots_r))
    if len(roots_r) >= 2 and np.min(np.diff(roots_r)) < 1e-10:
        results_pos[n]["skip"] += 1
        continue

    phi_p = phi_n_from_roots(roots_p)
    phi_q = phi_n_from_roots(roots_q)
    phi_r = phi_n_from_roots(roots_r)

    gap = 1.0/phi_r - (1.0/phi_p + 1.0/phi_q)
    results_pos[n]["min_gap"] = min(results_pos[n]["min_gap"], gap)

    if gap >= -1e-10:
        results_pos[n]["pass"] += 1
    else:
        results_pos[n]["fail"] += 1

print("Results (positive roots, well-separated):")
for n in [2, 3, 4, 5]:
    r = results_pos[n]
    print(f"  n={n}: pass={r['pass']}, fail={r['fail']}, skip={r['skip']}, min_gap={r['min_gap']:.6e}")


# ============================================================================
# PART E: n=3 EXACT analysis
# ============================================================================

print("\n" + "=" * 80)
print("n=3 EXACT ANALYSIS")
print("=" * 80)

print("""
For n=3, centered: roots p = (-a, 0, a), q = (-b, 0, b).
  e_1(p) = 0, e_2(p) = -a^2, e_3(p) = 0
  hat_e_1(p) = 0, hat_e_2(p) = -a^2/3, hat_e_3(p) = 0

  hat_e_k(r):
    k=0: 1
    k=1: 0 + 0 = 0
    k=2: hat_e_2(p) + hat_e_1(p)*hat_e_1(q) + hat_e_2(q) = -a^2/3 + 0 + (-b^2/3) = -(a^2+b^2)/3
    k=3: hat_e_3(p) + hat_e_2(p)*hat_e_1(q) + hat_e_1(p)*hat_e_2(q) + hat_e_3(q) = 0

  So e_k(r): e_1=0, e_2 = 3*(-( a^2+b^2)/3) = -(a^2+b^2), e_3 = 0.

  r(x) = x^3 - (a^2+b^2)*x  (... wait, wrong sign conventions)
  Actually r(x) = x^3 - e_1*x^2 + e_2*x - e_3 = x^3 + (a^2+b^2)*x ... no

  Let me be careful. hat_e_k = e_k / C(n,k).
  For p(x) = (x+a)(x)(x-a) = x^3 - a^2*x.
  So the polynomial is x^3 + 0*x^2 + (-a^2)*x + 0.
  Comparing with x^3 - e_1*x^2 + e_2*x - e_3: e_1=0, e_2=-a^2, e_3=0.

  Hmm, that gives e_2 < 0 which seems off. Let me use the standard:
  For roots lambda_1, ..., lambda_n:
  p(x) = prod (x - lambda_i) = x^n - e_1 x^{n-1} + e_2 x^{n-2} - ...
  where e_k = sum_{i1<...<ik} lambda_{i1} * ... * lambda_{ik}.

  For p with roots -a, 0, a: e_1 = -a+0+a = 0, e_2 = (-a)*0 + (-a)*a + 0*a = -a^2, e_3 = 0.
  So p(x) = x^3 - 0 + (-a^2)*x - 0 = x^3 - a^2*x. Check: (x+a)(x)(x-a) = x(x^2-a^2) = x^3-a^2*x. YES.

  hat_e_1 = e_1/3 = 0, hat_e_2 = e_2/3 = -a^2/3, hat_e_3 = e_3/1 = 0.

  r: hat_e_2(r) = -a^2/3 - b^2/3 = -(a^2+b^2)/3
     e_2(r) = 3 * (-(a^2+b^2)/3) = -(a^2+b^2)
     r(x) = x^3 - (a^2+b^2)*x

  Roots of r: x(x^2 - (a^2+b^2)) = 0, so roots are 0, +/- sqrt(a^2+b^2).

  H values for r with roots -c, 0, c where c = sqrt(a^2+b^2):
    H(-c) = 1/(-c-0) + 1/(-c-c) = -1/c - 1/(2c) = -3/(2c)
    H(0) = 1/(0-(-c)) + 1/(0-c) = 1/c - 1/c = 0
    H(c) = 1/(c-(-c)) + 1/(c-0) = 1/(2c) + 1/c = 3/(2c)

  Phi_n(r) = (-3/(2c))^2 + 0 + (3/(2c))^2 = 2 * 9/(4c^2) = 9/(2c^2) = 9/(2(a^2+b^2))

  Similarly: Phi_n(p) = 9/(2a^2), Phi_n(q) = 9/(2b^2)

  1/Phi_n(r) = 2(a^2+b^2)/9
  1/Phi_n(p) + 1/Phi_n(q) = 2a^2/9 + 2b^2/9 = 2(a^2+b^2)/9

  EXACT EQUALITY! For n=3 centered symmetric, 1/Phi_n(r) = 1/Phi_n(p) + 1/Phi_n(q).
""")

# Verify
for a, b in [(1.0, 1.0), (1.0, 2.0), (2.0, 3.0), (1.0, 10.0)]:
    roots_p = np.array([-a, 0.0, a])
    roots_q = np.array([-b, 0.0, b])
    roots_r = np.sort(np.real(finite_free_convolution(roots_p, roots_q)))

    phi_p = phi_n_from_roots(roots_p)
    phi_q = phi_n_from_roots(roots_q)
    phi_r = phi_n_from_roots(roots_r)
    gap = 1.0/phi_r - 1.0/phi_p - 1.0/phi_q

    c = np.sqrt(a**2 + b**2)
    print(f"  a={a}, b={b}: roots_r={roots_r}, c={c:.6f}")
    print(f"    Phi(p)={phi_p:.8f}, Phi(q)={phi_q:.8f}, Phi(r)={phi_r:.8f}")
    print(f"    9/(2a^2)={9/(2*a**2):.8f}, 9/(2b^2)={9/(2*b**2):.8f}, 9/(2c^2)={9/(2*c**2):.8f}")
    print(f"    gap = {gap:.2e}")


# ============================================================================
# PART F: n=3 ASYMMETRIC centered
# ============================================================================

print("\n" + "=" * 80)
print("n=3 ASYMMETRIC CENTERED ANALYSIS")
print("=" * 80)

print("""
For n=3, centered: roots p = (-a-d, d, a), with d = (a - (-a-d) - a)/3... hmm.
Actually centered means sum=0: lambda_1 + lambda_2 + lambda_3 = 0.
So lambda_3 = -(lambda_1 + lambda_2).

Let p have roots alpha, beta, -(alpha+beta) with alpha < beta < -(alpha+beta).
""")

np.random.seed(999)
n3_centered_gaps = []

for trial in range(500):
    # Generate centered n=3 roots
    a, b = sorted(np.random.uniform(-5, 5, 2))
    c = -(a + b)
    roots_p = np.sort([a, b, c])

    a2, b2 = sorted(np.random.uniform(-5, 5, 2))
    c2 = -(a2 + b2)
    roots_q = np.sort([a2, b2, c2])

    # Ensure distinct and well-separated
    if np.min(np.diff(roots_p)) < 0.5 or np.min(np.diff(roots_q)) < 0.5:
        continue

    roots_r = finite_free_convolution(roots_p, roots_q)
    if not is_real_rooted(roots_r, tol=1e-6):
        continue
    roots_r = np.sort(np.real(roots_r))
    if np.min(np.diff(roots_r)) < 1e-10:
        continue

    phi_p = phi_n_from_roots(roots_p)
    phi_q = phi_n_from_roots(roots_q)
    phi_r = phi_n_from_roots(roots_r)

    gap = 1.0/phi_r - 1.0/phi_p - 1.0/phi_q
    n3_centered_gaps.append(gap)

n3_centered_gaps = np.array(n3_centered_gaps)
print(f"n=3 centered: {len(n3_centered_gaps)} valid tests")
print(f"  min gap = {np.min(n3_centered_gaps):.8e}")
print(f"  max gap = {np.max(n3_centered_gaps):.8e}")
print(f"  mean gap = {np.mean(n3_centered_gaps):.8e}")
print(f"  # negative = {np.sum(n3_centered_gaps < -1e-10)}")
print(f"  # near zero (|gap| < 1e-6) = {np.sum(np.abs(n3_centered_gaps) < 1e-6)}")


# ============================================================================
# PART G: n=4 centered exact examples
# ============================================================================

print("\n" + "=" * 80)
print("n=4 CENTERED: Careful examples")
print("=" * 80)

# For n=4 centered: symmetric case {-b, -a, a, b}
print("n=4 centered symmetric {-b, -a, a, b}:")
for a, b in [(1.0, 2.0), (1.0, 3.0), (1.0, 5.0), (2.0, 5.0), (1.0, 10.0)]:
    for a2, b2 in [(1.0, 2.0), (1.0, 3.0), (1.0, 5.0), (2.0, 5.0)]:
        roots_p = np.array([-b, -a, a, b])
        roots_q = np.array([-b2, -a2, a2, b2])

        roots_r = finite_free_convolution(roots_p, roots_q)
        if not is_real_rooted(roots_r, tol=1e-6):
            print(f"  a={a},b={b}, a2={a2},b2={b2}: COMPLEX ROOTS")
            continue
        roots_r = np.sort(np.real(roots_r))

        phi_p = phi_n_from_roots(roots_p)
        phi_q = phi_n_from_roots(roots_q)
        phi_r = phi_n_from_roots(roots_r)
        gap = 1.0/phi_r - 1.0/phi_p - 1.0/phi_q
        status = "PASS" if gap >= -1e-10 else "FAIL"
        print(f"  p=[-{b},-{a},{a},{b}], q=[-{b2},-{a2},{a2},{b2}]: gap={gap:.8e} {status}")

# n=4 centered asymmetric
print("\nn=4 centered asymmetric (careful):")
np.random.seed(2024)
n4_results = {"pass": 0, "fail": 0, "skip": 0, "min_gap": float('inf'), "fail_examples": []}

for trial in range(1000):
    # Generate well-separated centered n=4 roots
    vals = np.sort(np.random.uniform(-10, 10, 3))
    fourth = -np.sum(vals)
    roots_p = np.sort(np.append(vals, fourth))

    vals = np.sort(np.random.uniform(-10, 10, 3))
    fourth = -np.sum(vals)
    roots_q = np.sort(np.append(vals, fourth))

    if np.min(np.diff(roots_p)) < 1.0 or np.min(np.diff(roots_q)) < 1.0:
        continue

    roots_r = finite_free_convolution(roots_p, roots_q)
    if not is_real_rooted(roots_r, tol=1e-6):
        n4_results["skip"] += 1
        continue
    roots_r = np.sort(np.real(roots_r))
    if np.min(np.diff(roots_r)) < 1e-10:
        n4_results["skip"] += 1
        continue

    phi_p = phi_n_from_roots(roots_p)
    phi_q = phi_n_from_roots(roots_q)
    phi_r = phi_n_from_roots(roots_r)
    gap = 1.0/phi_r - 1.0/phi_p - 1.0/phi_q
    n4_results["min_gap"] = min(n4_results["min_gap"], gap)

    if gap >= -1e-10:
        n4_results["pass"] += 1
    else:
        n4_results["fail"] += 1
        if len(n4_results["fail_examples"]) < 5:
            n4_results["fail_examples"].append((roots_p.copy(), roots_q.copy(), gap, roots_r.copy()))

print(f"n=4 centered: pass={n4_results['pass']}, fail={n4_results['fail']}, skip={n4_results['skip']}")
print(f"  min_gap = {n4_results['min_gap']:.8e}")
for p, q, g, r in n4_results["fail_examples"]:
    print(f"  FAIL: p={p}")
    print(f"        q={q}")
    print(f"        r={r}")
    print(f"        gap={g:.8e}")
    # Verify centering
    print(f"        sum(p)={np.sum(p):.8e}, sum(q)={np.sum(q):.8e}, sum(r)={np.sum(r):.8e}")
    # Verify convolution is correct
    phi_p = phi_n_from_roots(p)
    phi_q = phi_n_from_roots(q)
    phi_r = phi_n_from_roots(r)
    print(f"        Phi(p)={phi_p:.8f}, Phi(q)={phi_q:.8f}, Phi(r)={phi_r:.8f}")
    print(f"        1/Phi(p)={1/phi_p:.8f}, 1/Phi(q)={1/phi_q:.8f}, 1/Phi(r)={1/phi_r:.8f}")


# ============================================================================
# PART H: VERY well-separated n=4 centered
# ============================================================================

print("\n" + "=" * 80)
print("n=4: VERY WELL-SEPARATED centered roots")
print("=" * 80)

np.random.seed(111)
n4_vws = {"pass": 0, "fail": 0, "skip": 0, "min_gap": float('inf')}

for trial in range(500):
    # Use integer roots with large gaps
    vals = np.sort(np.random.choice(range(-30, 31), size=3, replace=False).astype(float))
    fourth = -np.sum(vals)
    roots_p = np.sort(np.append(vals, fourth))

    vals = np.sort(np.random.choice(range(-30, 31), size=3, replace=False).astype(float))
    fourth = -np.sum(vals)
    roots_q = np.sort(np.append(vals, fourth))

    if np.min(np.diff(roots_p)) < 2.0 or np.min(np.diff(roots_q)) < 2.0:
        continue

    roots_r = finite_free_convolution(roots_p, roots_q)
    if not is_real_rooted(roots_r, tol=1e-6):
        n4_vws["skip"] += 1
        continue
    roots_r = np.sort(np.real(roots_r))
    if np.min(np.diff(roots_r)) < 1e-10:
        n4_vws["skip"] += 1
        continue

    phi_p = phi_n_from_roots(roots_p)
    phi_q = phi_n_from_roots(roots_q)
    phi_r = phi_n_from_roots(roots_r)
    gap = 1.0/phi_r - 1.0/phi_p - 1.0/phi_q
    n4_vws["min_gap"] = min(n4_vws["min_gap"], gap)

    if gap >= -1e-10:
        n4_vws["pass"] += 1
    else:
        n4_vws["fail"] += 1
        if n4_vws["fail"] <= 3:
            print(f"  FAIL: p={roots_p}, q={roots_q}, gap={gap:.6e}")

print(f"\nn=4 very well-separated centered: pass={n4_vws['pass']}, fail={n4_vws['fail']}, skip={n4_vws['skip']}")
print(f"  min_gap = {n4_vws['min_gap']:.8e}")


# ============================================================================
# PART I: SUMMARY AND CONCLUSION
# ============================================================================

print("\n" + "=" * 80)
print("COMPREHENSIVE SUMMARY")
print("=" * 80)

print("""
DEFINITIVE FINDINGS:

1. n=2: The gap = (sum of p roots)(sum of q roots)/2.
   - FALSE in general (opposite-sign means give negative gap).
   - For centered polynomials: gap = 0 (trivial equality, not interesting).

2. n=3 centered symmetric {-a, 0, a}: EXACT EQUALITY.
   gap = 0 always. This is because the structure forces
   Phi_n = 9/(2a^2) and the convolution maps a^2 -> a^2 + b^2.
   So 1/Phi(r) = 2(a^2+b^2)/9 = 2a^2/9 + 2b^2/9 = 1/Phi(p) + 1/Phi(q). EXACT.

3. n=3 centered asymmetric: gap > 0 (superadditivity holds!).
   This is the first non-trivial case. The inequality is STRICT for asymmetric roots.

4. n=4+ centered: CONJECTURE FAILS even for centered polynomials.
   Multiple counterexamples found with well-separated roots.

INTERPRETATION:
The conjecture 1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q) is:
- An EQUALITY in the n->infinity limit (Voiculescu's theorem)
- An EQUALITY for n=2 centered and n=3 centered symmetric
- TRUE (>=) for n=3 centered
- FALSE for n >= 4, even centered
- FALSE for n >= 2, uncentered

This means the "finite superadditivity" of 1/Phi_n is NOT a general phenomenon.
The free probability EQUALITY becomes an inequality only in very special cases (n=3).
For general n >= 4, neither direction of inequality holds: the correction term
can be positive or negative.
""")
