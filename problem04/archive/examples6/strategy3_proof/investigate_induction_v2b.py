"""
investigate_induction_v2b.py — CORRECTED investigation.

Issue from v2: The boxplus_mss function had a bug in root extraction.
The MSS coefficient formula gives c_k, and the polynomial is:
  r(x) = sum_{k=0}^n (-1)^k C(n,k) c_k x^{n-k}
with c_0 = 1 (monic).

Also: the CONJECTURE F_n >= 0 was being violated. Need to check if this
is a bug or a genuine failure.

Let me also verify against the correct MSS formula more carefully.
"""

import numpy as np
from math import factorial, comb
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# CORE FUNCTIONS — carefully verified
# ============================================================

def poly_from_roots(roots):
    """Monic polynomial from roots, highest degree first."""
    return np.poly(roots)

def roots_from_poly(coeffs):
    """Roots of polynomial (highest degree first)."""
    return np.roots(coeffs)

def boxplus_coeffs(roots_p, roots_q):
    """
    Compute the MSS finite free additive convolution p boxplus_n q.

    If p(x) = sum_{k=0}^n (-1)^k C(n,k) a_k x^{n-k} (a_0=1)
    and q(x) = sum_{k=0}^n (-1)^k C(n,k) b_k x^{n-k} (b_0=1)
    then r = p boxplus q has:
    r(x) = sum_{k=0}^n (-1)^k C(n,k) c_k x^{n-k}
    where c_k = sum_{i+j=k} [C(n-i,j)^{-1} * C(n,j)^{-1} * C(n,k)] * a_i * b_j

    Actually the MSS formula is:
    c_k = sum_{i+j=k} [(n-i)!(n-j)! / (n!(n-k)!)] * a_i * b_j

    where a_k = e_k(roots_p) (elementary symmetric polynomials of roots of p).
    """
    n = len(roots_p)
    assert len(roots_q) == n

    # Extract a_k = e_k(roots_p) and b_k = e_k(roots_q)
    # np.poly gives [1, -e_1, e_2, -e_3, ...] for monic polynomial
    poly_p = np.poly(roots_p)  # [1, -e_1, e_2, ...]
    poly_q = np.poly(roots_q)

    # a_k = e_k(roots_p), with sign convention: poly_p[k] = (-1)^k * C(n,k)^{-1} ... no
    # Actually np.poly gives the standard polynomial coefficients:
    # p(x) = poly_p[0]*x^n + poly_p[1]*x^{n-1} + ... + poly_p[n]
    # For monic poly with roots r_1,...,r_n:
    # p(x) = x^n - e_1 x^{n-1} + e_2 x^{n-2} - ... + (-1)^n e_n
    # So poly_p[k] = (-1)^k * e_k(roots_p)

    # The MSS "normalized" coefficients:
    # Write p(x) = sum_{k=0}^n (-1)^k C(n,k) a_k x^{n-k}
    # Then (-1)^k C(n,k) a_k = poly_p[k] = (-1)^k e_k
    # So a_k = e_k / C(n,k)

    a = np.zeros(n+1)
    b = np.zeros(n+1)
    for k in range(n+1):
        a[k] = (-1)**k * poly_p[k] / comb(n, k)  # a_k = e_k / C(n,k)
        b[k] = (-1)**k * poly_q[k] / comb(n, k)

    # MSS convolution: c_k = sum_{i+j=k} [(n-i)!(n-j)! / (n!(n-k)!)] * a_i * b_j
    c = np.zeros(n+1)
    for k in range(n+1):
        for i in range(k+1):
            j = k - i
            if i <= n and j <= n:
                coeff = factorial(n-i) * factorial(n-j) / (factorial(n) * factorial(n-k))
                c[k] += coeff * a[i] * b[j]

    # Reconstruct polynomial: r(x) = sum_{k=0}^n (-1)^k C(n,k) c_k x^{n-k}
    poly_r = np.zeros(n+1)
    for k in range(n+1):
        poly_r[k] = (-1)**k * comb(n, k) * c[k]

    return poly_r, c

def boxplus(roots_p, roots_q):
    """Compute p boxplus q and return sorted real roots."""
    poly_r, c = boxplus_coeffs(roots_p, roots_q)
    raw_roots = np.roots(poly_r)
    return np.sort(np.real(raw_roots)), c

def H_values(roots):
    """H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j)."""
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                d = roots[i] - roots[j]
                if abs(d) < 1e-12:
                    return None
                H[i] += 1.0 / d
    return H

def Phi(roots):
    """Phi_n(p) = sum_i H_p(lambda_i)^2."""
    H = H_values(roots)
    if H is None:
        return None
    return np.sum(H**2)

def critical_points(roots):
    """Roots of p'(x) where p(x) = prod(x - r_i)."""
    poly = np.poly(roots)
    dpoly = np.polyder(poly)
    cps = np.roots(dpoly)
    if np.any(np.abs(np.imag(cps)) > 1e-8):
        return None
    return np.sort(np.real(cps))

def make_roots(n, min_gap=0.5, scale=2.0):
    """Generate n well-separated real roots."""
    roots = np.sort(np.random.randn(n) * scale)
    for i in range(1, n):
        if roots[i] - roots[i-1] < min_gap:
            roots[i] = roots[i-1] + min_gap
    return roots


# ============================================================
# TEST 0: Verify boxplus implementation
# ============================================================
print("=" * 70)
print("TEST 0: VERIFY BOXPLUS IMPLEMENTATION")
print("=" * 70)

# For n=2: p(x) = x^2 - (a+b)x + ab with roots a,b
# q(x) = x^2 - (c+d)x + cd with roots c,d
# p boxplus q should have roots that are a+c, b+d... no, that's not right.
# For n=2 the convolution is just shift: p boxplus q(x) = product over pairs...
# Actually for degree 2, a_1 = (lambda_1+lambda_2)/2, a_2 = lambda_1*lambda_2/1
# Let me just check: boxplus of p(x)=x^2 with q(x)=x^2 (both roots at 0)
# should give r(x) = x^2 (roots at 0).
print("\n  Test: p=q=x^2 (roots at 0,0)")
roots_r, _ = boxplus([0.0, 1.0], [0.0, 1.0])
print(f"    boxplus roots: {roots_r}")
# For p = (x)(x-1) = x^2 - x, q = same
# a_0=1, a_1 = e_1/C(2,1) = 1/2, a_2 = e_2/C(2,2) = 0/1 = 0
# c_k = sum_{i+j=k} [(2-i)!(2-j)!/(2!(2-k)!)] a_i b_j
# c_0 = a_0*b_0 * [2!*2!/(2!*2!)] = 1
# c_1 = [1!*2!/(2!*1!) * a_1*b_0 + 2!*1!/(2!*1!) * a_0*b_1] = [1/2 * 1/2 + 1/2] = 3/4... hmm
# Let me just check empirically.
print(f"    Expected: roots of boxplus should sum to sum(p_roots)+sum(q_roots) = 2")
print(f"    Got sum: {sum(roots_r):.6f}")

# Check against known result: boxplus shifts means
# sum of roots of r = sum of roots of p + sum of roots of q
for trial in range(5):
    n = 4
    rp = make_roots(n)
    rq = make_roots(n)
    rr, _ = boxplus(rp, rq)
    print(f"  n={n}: sum_p={sum(rp):.4f}, sum_q={sum(rq):.4f}, sum_r={sum(rr):.4f}, sum_p+sum_q={sum(rp)+sum(rq):.4f}")


# ============================================================
# TEST 1: Verify MSS derivative identity MORE CAREFULLY
# ============================================================
print("\n" + "=" * 70)
print("TEST 1: MSS DERIVATIVE IDENTITY (careful)")
print("  r'(x) = n * (p^{(1)} boxplus_{n-1} q^{(1)})(x)")
print("  Equivalently: critical points of r = boxplus of critical points")
print("=" * 70)

for trial in range(15):
    n = np.random.choice([3, 4, 5])
    rp = make_roots(n, min_gap=1.0)
    rq = make_roots(n, min_gap=1.0)

    rr, _ = boxplus(rp, rq)

    # Check r has real, simple roots
    if np.any(np.diff(rr) < 0.01):
        print(f"  Trial {trial}: r has near-degenerate roots, skip")
        continue

    cp_p = critical_points(rp)
    cp_q = critical_points(rq)
    cp_r = critical_points(rr)

    if cp_p is None or cp_q is None or cp_r is None:
        continue
    if len(cp_p) != n-1 or len(cp_q) != n-1:
        continue

    # s = p^{(1)} boxplus_{n-1} q^{(1)}
    s_roots, _ = boxplus(cp_p, cp_q)

    error = np.max(np.abs(np.sort(cp_r) - np.sort(s_roots)))
    status = "OK" if error < 1e-6 else "FAIL"
    print(f"  n={n}: error = {error:.2e} [{status}]")
    if error > 0.01:
        print(f"    rp = {np.array2string(rp, precision=4)}")
        print(f"    rq = {np.array2string(rq, precision=4)}")
        print(f"    rr = {np.array2string(rr, precision=4)}")
        print(f"    cp_r = {np.array2string(cp_r, precision=4)}")
        print(f"    s    = {np.array2string(s_roots, precision=4)}")


# ============================================================
# TEST 2: Verify the CONJECTURE 1/Phi(r) >= 1/Phi(p) + 1/Phi(q)
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: CONJECTURE VERIFICATION")
print("  1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q) where r = p boxplus q")
print("=" * 70)

for n in [2, 3, 4, 5, 6]:
    violations = 0
    total = 0
    worst = 0
    for trial in range(200):
        rp = make_roots(n, min_gap=0.5)
        rq = make_roots(n, min_gap=0.5)
        rr, _ = boxplus(rp, rq)

        if np.any(np.diff(rr) < 0.01):
            continue

        ph_p = Phi(rp)
        ph_q = Phi(rq)
        ph_r = Phi(rr)

        if any(v is None or v < 1e-12 for v in [ph_p, ph_q, ph_r]):
            continue

        total += 1
        excess = 1/ph_r - 1/ph_p - 1/ph_q
        if excess < -1e-10:
            violations += 1
            if excess < worst:
                worst = excess
                worst_data = (rp.copy(), rq.copy(), rr.copy(), ph_p, ph_q, ph_r)

    print(f"  n={n}: {violations}/{total} violations, worst excess = {worst:.8f}")
    if violations > 0 and 'worst_data' in dir():
        rp, rq, rr, pp, pq, pr = worst_data
        print(f"    Worst case:")
        print(f"      p roots: {np.array2string(rp, precision=4)}")
        print(f"      q roots: {np.array2string(rq, precision=4)}")
        print(f"      r roots: {np.array2string(rr, precision=4)}")
        print(f"      1/Phi_p={1/pp:.6f}, 1/Phi_q={1/pq:.6f}, 1/Phi_r={1/pr:.6f}")


# ============================================================
# TEST 3: Check n=2 base case (should be exact equality)
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: BASE CASE n=2")
print("=" * 70)

for trial in range(10):
    rp = make_roots(2, min_gap=0.5)
    rq = make_roots(2, min_gap=0.5)
    rr, _ = boxplus(rp, rq)

    ph_p = Phi(rp)
    ph_q = Phi(rq)
    ph_r = Phi(rr)

    excess = 1/ph_r - 1/ph_p - 1/ph_q
    print(f"  p={np.array2string(rp,precision=3)}, q={np.array2string(rq,precision=3)}")
    print(f"    Phi_p={ph_p:.6f}, Phi_q={ph_q:.6f}, Phi_r={ph_r:.6f}")
    print(f"    excess = {excess:.10f}")


# ============================================================
# TEST 4: Check computation of Phi for uniform spacing
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: Phi FOR UNIFORM SPACING (sanity check)")
print("=" * 70)

# For roots 0, d, 2d, ..., (n-1)d:
# H_p(lambda_i) = sum_{j!=i} 1/((i-j)*d) = (1/d) * sum_{j!=i} 1/(i-j)
# These are known: H_i = (1/d) * [H_n - 1/(n-i) - ... + 1/i - ...] (harmonic numbers)

for n in [2, 3, 4, 5]:
    d = 1.0
    roots = np.array([i * d for i in range(n)])
    H = H_values(roots)
    ph = Phi(roots)
    print(f"  n={n}, d={d}: roots={roots}")
    print(f"    H = {np.array2string(H, precision=6)}")
    print(f"    sum(H) = {sum(H):.10f} (should be 0)")
    print(f"    Phi = {ph:.6f}")
    print(f"    1/Phi = {1/ph:.6f}")
    print(f"    (1/Phi)/d^2 = {1/(ph*d**2):.6f} (universal constant for this n)")


# ============================================================
# Now the MAIN investigation with corrected boxplus
# ============================================================
print("\n" + "=" * 70)
print("MAIN: INDUCTION ANALYSIS WITH CORRECTED BOXPLUS")
print("=" * 70)

# Investigation 2 (corrected): Ratio Phi_n / Phi_{n-1}
print("\n--- Phi_n(p) / Phi_{n-1}(p^{(1)}) for uniform spacing ---")
for n in [3, 4, 5, 6, 7, 8]:
    d = 1.0
    roots = np.array([i * d for i in range(n)])
    cp = critical_points(roots)
    if cp is not None and len(cp) == n-1:
        ph_n = Phi(roots)
        ph_n1 = Phi(cp)
        if ph_n is not None and ph_n1 is not None and ph_n1 > 0:
            print(f"  n={n}: ratio = {ph_n/ph_n1:.8f}")

# Investigation 7 (corrected): delta superadditivity
print("\n--- DELTA SUPERADDITIVITY CHECK (corrected boxplus) ---")
for n in [3, 4, 5]:
    print(f"\n  n={n}:")
    delta_data = []
    for trial in range(200):
        rp = make_roots(n, min_gap=0.5 + 0.5*np.random.rand())
        rq = make_roots(n, min_gap=0.5 + 0.5*np.random.rand())
        rr, _ = boxplus(rp, rq)

        if np.any(np.diff(rr) < 0.01):
            continue

        cp_p = critical_points(rp)
        cp_q = critical_points(rq)
        if cp_p is None or cp_q is None:
            continue
        if len(cp_p) != n-1 or len(cp_q) != n-1:
            continue

        s_roots, _ = boxplus(cp_p, cp_q)
        if np.any(np.diff(s_roots) < 0.01):
            continue

        ph_p = Phi(rp)
        ph_q = Phi(rq)
        ph_r = Phi(rr)
        ph_cp = Phi(cp_p)
        ph_cq = Phi(cp_q)
        ph_s = Phi(s_roots)

        if any(v is None or v < 1e-12 for v in [ph_p, ph_q, ph_r, ph_cp, ph_cq, ph_s]):
            continue

        F_n = 1/ph_r - 1/ph_p - 1/ph_q
        F_n1 = 1/ph_s - 1/ph_cp - 1/ph_cq

        delta_p = 1/ph_p - 1/ph_cp
        delta_q = 1/ph_q - 1/ph_cq
        delta_r = 1/ph_r - 1/ph_s

        delta_data.append({
            'F_n': F_n, 'F_n1': F_n1,
            'delta_p': delta_p, 'delta_q': delta_q, 'delta_r': delta_r,
            'delta_diff': delta_r - delta_p - delta_q,
            'ph_p': ph_p, 'ph_q': ph_q, 'ph_r': ph_r,
        })

    if delta_data:
        F_n_arr = np.array([d['F_n'] for d in delta_data])
        F_n1_arr = np.array([d['F_n1'] for d in delta_data])
        dd_arr = np.array([d['delta_diff'] for d in delta_data])

        print(f"    {len(delta_data)} valid samples")
        print(f"    F_n:  min={F_n_arr.min():.8f}, max={F_n_arr.max():.8f}")
        print(f"    F_n all >= 0: {np.all(F_n_arr >= -1e-10)}")
        n_violations = np.sum(F_n_arr < -1e-10)
        if n_violations > 0:
            print(f"    *** {n_violations} VIOLATIONS of the main conjecture! ***")
        print(f"    F_{n-1}: min={F_n1_arr.min():.8f}, max={F_n1_arr.max():.8f}")
        print(f"    delta_diff = delta_r - delta_p - delta_q:")
        print(f"      min={dd_arr.min():.8f}, max={dd_arr.max():.8f}")
        print(f"      all >= 0 (superadditive): {np.all(dd_arr >= -1e-10)}")

        # Verify: F_n = F_{n-1} + delta_diff
        check = F_n_arr - (F_n1_arr + dd_arr)
        print(f"    Consistency check: max|F_n - (F_{n-1}+dd)| = {np.max(np.abs(check)):.2e}")


# ============================================================
# Investigation: SIGN pattern of delta
# ============================================================
print("\n" + "=" * 70)
print("SIGN ANALYSIS OF delta_p = 1/Phi_n(p) - 1/Phi_{n-1}(p^{(1)})")
print("=" * 70)

for n in [3, 4, 5, 6]:
    deltas = []
    for trial in range(200):
        roots = make_roots(n, min_gap=0.3 + np.random.rand())
        cp = critical_points(roots)
        if cp is None or len(cp) != n-1:
            continue
        if len(cp) > 1 and np.any(np.diff(cp) < 0.01):
            continue
        ph = Phi(roots)
        ph1 = Phi(cp)
        if ph is None or ph1 is None:
            continue
        deltas.append(1/ph - 1/ph1)

    if deltas:
        deltas = np.array(deltas)
        print(f"  n={n}: delta range = [{deltas.min():.8f}, {deltas.max():.8f}]")
        print(f"    all negative: {np.all(deltas < 0)}")
        print(f"    Interpretation: Phi_n > Phi_{n-1} always (more terms => bigger)")


# ============================================================
# KEY QUESTION: Is the conjecture even TRUE?
# ============================================================
print("\n" + "=" * 70)
print("CRITICAL: IS THE CONJECTURE TRUE?")
print("Focus on small n where we can verify carefully")
print("=" * 70)

# n=2: EXACT equality (Cauchy-Schwarz style)
print("\n--- n=2: expect EQUALITY ---")
for trial in range(5):
    a, b = sorted(np.random.randn(2) * 2)
    if b - a < 0.5:
        b = a + 0.5
    c, d = sorted(np.random.randn(2) * 2)
    if d - c < 0.5:
        d = c + 0.5
    rp = np.array([a, b])
    rq = np.array([c, d])
    rr, _ = boxplus(rp, rq)

    # For n=2: H_p(a) = 1/(a-b), H_p(b) = 1/(b-a) = -H_p(a)
    # Phi_2(p) = 2/(b-a)^2
    # 1/Phi_2(p) = (b-a)^2/2
    gap_p = b - a
    gap_q = d - c
    gap_r = rr[1] - rr[0]

    ph_p = 2/gap_p**2
    ph_q = 2/gap_q**2
    ph_r = 2/gap_r**2

    excess = gap_r**2/2 - gap_p**2/2 - gap_q**2/2
    print(f"  gaps: p={gap_p:.4f}, q={gap_q:.4f}, r={gap_r:.4f}")
    print(f"  gap_r^2 = {gap_r**2:.6f}, gap_p^2 + gap_q^2 = {gap_p**2+gap_q**2:.6f}")
    print(f"  excess = {excess:.8f}")

# For n=2 boxplus: gap_r^2 = gap_p^2 + gap_q^2 (PYTHAGOREAN!)
# This gives EQUALITY in the Fisher inequality.

print("\n--- n=3: careful check ---")
for trial in range(20):
    rp = make_roots(3, min_gap=0.5 + np.random.rand())
    rq = make_roots(3, min_gap=0.5 + np.random.rand())
    rr, _ = boxplus(rp, rq)

    if np.any(np.diff(rr) < 0.01):
        continue

    ph_p = Phi(rp)
    ph_q = Phi(rq)
    ph_r = Phi(rr)

    if any(v is None for v in [ph_p, ph_q, ph_r]):
        continue

    excess = 1/ph_r - 1/ph_p - 1/ph_q
    print(f"  excess = {excess:.10f}  {'OK' if excess >= -1e-10 else 'VIOLATION!'}")
    if excess < -1e-10:
        print(f"    rp = {rp}")
        print(f"    rq = {rq}")
        print(f"    rr = {rr}")
        print(f"    Phi: p={ph_p:.6f}, q={ph_q:.6f}, r={ph_r:.6f}")
        print(f"    1/Phi: p={1/ph_p:.6f}, q={1/ph_q:.6f}, r={1/ph_r:.6f}")


# ============================================================
# Double-check: compute boxplus a different way
# ============================================================
print("\n" + "=" * 70)
print("CROSS-CHECK: Alternative boxplus computation")
print("=" * 70)

def boxplus_via_poly(roots_p, roots_q):
    """Alternative: directly from the polynomial coefficient formula.

    For monic polynomials p(x) = x^n + p_{n-1} x^{n-1} + ... + p_0
    and q(x) = x^n + q_{n-1} x^{n-1} + ... + q_0, define
    a_k = (-1)^k p_{n-k} / C(n,k) and similarly for b_k.

    Then c_k = sum_{i+j=k} [(n-i)!(n-j)!/(n!(n-k)!)] a_i b_j
    and r(x) = x^n + sum_{k=1}^n (-1)^k C(n,k) c_k x^{n-k}
    """
    n = len(roots_p)
    poly_p = np.poly(roots_p)  # [1, coeff of x^{n-1}, ..., coeff of x^0]
    poly_q = np.poly(roots_q)

    # a_k: the "normalized" elementary symmetric functions
    a = np.zeros(n+1)
    b = np.zeros(n+1)
    for k in range(n+1):
        # poly_p[k] is the coefficient of x^{n-k}
        # poly_p[k] = (-1)^k * e_k where e_k = C(n,k)*a_k
        # So a_k = (-1)^k * poly_p[k] / C(n,k)
        a[k] = (-1)**k * poly_p[k] / comb(n, k)
        b[k] = (-1)**k * poly_q[k] / comb(n, k)

    # c_k = sum_{i+j=k} [(n-i)!(n-j)!/(n!(n-k)!)] a_i b_j
    c = np.zeros(n+1)
    for k in range(n+1):
        for i in range(k+1):
            j = k - i
            if 0 <= j <= n:
                w = factorial(n-i) * factorial(n-j) / (factorial(n) * factorial(n-k))
                c[k] += w * a[i] * b[j]

    # r(x) = sum_{k=0}^n (-1)^k C(n,k) c_k x^{n-k}
    # = x^n + (-1)^1 C(n,1) c_1 x^{n-1} + (-1)^2 C(n,2) c_2 x^{n-2} + ...
    poly_r = np.zeros(n+1)
    for k in range(n+1):
        poly_r[k] = (-1)**k * comb(n, k) * c[k]

    raw_roots = np.roots(poly_r)
    return np.sort(np.real(raw_roots))

# Compare the two implementations
for trial in range(5):
    n = 4
    rp = make_roots(n, min_gap=1.0)
    rq = make_roots(n, min_gap=1.0)

    rr1, _ = boxplus(rp, rq)
    rr2 = boxplus_via_poly(rp, rq)

    error = np.max(np.abs(rr1 - rr2))
    print(f"  Trial {trial}: error between implementations = {error:.2e}")


# ============================================================
# ADDITIONAL: Check the Cauchy-Schwarz / variance interpretation
# ============================================================
print("\n" + "=" * 70)
print("VARIANCE INTERPRETATION")
print("  For n=2: 1/Phi_2(p) = gap^2/2 = Var(roots)")
print("  For general n: is 1/Phi_n related to a variance?")
print("=" * 70)

for n in [2, 3, 4, 5, 6]:
    for trial in range(3):
        roots = make_roots(n, min_gap=0.5 + np.random.rand())
        ph = Phi(roots)
        if ph is None:
            continue
        var = np.var(roots)
        mean_gap = np.mean(np.diff(roots))
        print(f"  n={n}: 1/Phi={1/ph:.6f}, Var={var:.6f}, mean_gap^2={mean_gap**2:.6f}, ratio(1/Phi)/Var={1/(ph*var):.6f}")


# ============================================================
# FINAL: Check if the n=2 Pythagorean property extends
# ============================================================
print("\n" + "=" * 70)
print("PYTHAGOREAN CHECK: 1/Phi(r) vs 1/Phi(p) + 1/Phi(q)")
print("  For n=2: EQUALITY (gap^2 is additive)")
print("  For n>2: ?")
print("=" * 70)

for n in [2, 3, 4, 5]:
    print(f"\n  n={n}:")
    for trial in range(10):
        rp = make_roots(n, min_gap=1.0)
        rq = make_roots(n, min_gap=1.0)
        rr, _ = boxplus(rp, rq)

        if np.any(np.diff(rr) < 0.01):
            continue

        ph_p = Phi(rp)
        ph_q = Phi(rq)
        ph_r = Phi(rr)

        if any(v is None for v in [ph_p, ph_q, ph_r]):
            continue

        lhs = 1/ph_r
        rhs = 1/ph_p + 1/ph_q
        excess = lhs - rhs
        ratio = lhs / rhs if rhs != 0 else float('inf')

        print(f"    1/Phi_r={lhs:.6f}, 1/Phi_p+1/Phi_q={rhs:.6f}, excess={excess:.8f}, ratio={ratio:.6f}")
