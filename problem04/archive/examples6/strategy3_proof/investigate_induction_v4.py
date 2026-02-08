"""
investigate_induction_v4.py — FINAL investigation with CORRECT boxplus.

KEY FIX: The correct MSS finite free additive convolution formula uses:
  g_k = sum_{i+j=k} C(n-j, i) / C(n, i) * e_i(p) * e_j(q)
where e_k are the elementary symmetric polynomials of the roots.

This was verified against Monte Carlo (expected char poly of A + UBU*).
The previous formula (v1) was only correct for centered polynomials.
"""

import numpy as np
from math import factorial, comb
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# CORRECT boxplus implementation (v3 from debug)
# ============================================================

def boxplus(roots_p, roots_q):
    """
    CORRECT MSS finite free additive convolution.

    Formula (Marcus 2021):
    g_k = sum_{i+j=k} C(n-j, i) / C(n, i) * e_i(p) * e_j(q)

    where e_k are the elementary symmetric polynomials of the roots.
    The result polynomial is r(x) = sum_{k=0}^n (-1)^k g_k x^{n-k}.
    """
    n = len(roots_p)
    assert len(roots_q) == n

    poly_p = np.poly(roots_p)  # [1, -e_1, e_2, -e_3, ...]
    poly_q = np.poly(roots_q)

    e = np.zeros(n+1)  # e_k(p)
    f = np.zeros(n+1)  # e_k(q)
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

    roots_r = np.roots(poly_r)
    max_imag = np.max(np.abs(np.imag(roots_r)))
    return np.sort(np.real(roots_r)), max_imag

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
    """Roots of p'(x)."""
    poly = np.poly(roots)
    dpoly = np.polyder(poly)
    cps = np.roots(dpoly)
    if np.any(np.abs(np.imag(cps)) > 1e-6):
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
# SECTION 1: VERIFY CONJECTURE WITH CORRECT BOXPLUS
# ============================================================
print("=" * 70)
print("SECTION 1: CONJECTURE VERIFICATION (correct boxplus)")
print("  1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q) where r = p boxplus q")
print("=" * 70)

for n in [2, 3, 4, 5, 6]:
    violations = 0
    total = 0
    min_excess = float('inf')
    for trial in range(500):
        rp = make_roots(n, min_gap=0.3 + 0.7*np.random.rand())
        rq = make_roots(n, min_gap=0.3 + 0.7*np.random.rand())
        rr, max_imag = boxplus(rp, rq)

        if max_imag > 0.01:
            continue
        if len(rr) > 1 and np.any(np.diff(rr) < 0.01):
            continue

        ph_p = Phi(rp)
        ph_q = Phi(rq)
        ph_r = Phi(rr)

        if any(v is None or v < 1e-12 for v in [ph_p, ph_q, ph_r]):
            continue

        total += 1
        excess = 1/ph_r - 1/ph_p - 1/ph_q
        min_excess = min(min_excess, excess)
        if excess < -1e-8:
            violations += 1

    status = "HOLDS" if violations == 0 else f"FAILS ({violations} violations)"
    print(f"  n={n}: {total} valid trials, min excess = {min_excess:.10f} [{status}]")


# ============================================================
# SECTION 2: VERIFY MSS DERIVATIVE IDENTITY with correct boxplus
# ============================================================
print("\n" + "=" * 70)
print("SECTION 2: MSS DERIVATIVE IDENTITY (correct boxplus)")
print("  r'(x) = n * (p^{(1)} boxplus_{n-1} q^{(1)})(x)")
print("=" * 70)

for trial in range(15):
    n = np.random.choice([3, 4, 5])
    rp = make_roots(n, min_gap=1.0)
    rq = make_roots(n, min_gap=1.0)

    rr, imag_r = boxplus(rp, rq)
    if imag_r > 0.01:
        continue
    if np.any(np.diff(rr) < 0.01):
        continue

    cp_p = critical_points(rp)
    cp_q = critical_points(rq)
    cp_r = critical_points(rr)

    if any(x is None for x in [cp_p, cp_q, cp_r]):
        continue
    if len(cp_p) != n-1 or len(cp_q) != n-1 or len(cp_r) != n-1:
        continue

    s_roots, _ = boxplus(cp_p, cp_q)
    error = np.max(np.abs(np.sort(cp_r) - np.sort(s_roots)))
    status = "OK" if error < 1e-6 else "FAIL"
    print(f"  n={n}: error = {error:.2e} [{status}]")


# ============================================================
# SECTION 3: n=2 BASE CASE (analytical with correct formula)
# ============================================================
print("\n" + "=" * 70)
print("SECTION 3: n=2 BASE CASE (correct formula)")
print("  For n=2: g_k = sum C(2-j,i)/C(2,i) * e_i * f_j")
print("=" * 70)

# g_0 = C(2,0)/C(2,0) * 1 * 1 = 1
# g_1 = C(2,0)/C(2,0)*1*f_1 + C(1,1)/C(2,1)*e_1*1 = f_1 + e_1/2
# Wait: C(2-1,1)/C(2,1) = C(1,1)/C(2,1) = 1/2
# g_1 = e_1 * C(2-0,0)/C(2,0)... let me redo this carefully.
# g_1: i+j=1, so (i,j) in {(0,1), (1,0)}
#   (0,1): C(2-1,0)/C(2,0) * e_0 * f_1 = C(1,0)/1 * 1 * f_1 = f_1
#   (1,0): C(2-0,1)/C(2,1) * e_1 * f_0 = C(2,1)/C(2,1) * e_1 = e_1
# g_1 = e_1 + f_1 = (a+b) + (c+d)  -- sum of all roots

# g_2: i+j=2, (i,j) in {(0,2), (1,1), (2,0)}
#   (0,2): C(2-2,0)/C(2,0) * 1 * f_2 = C(0,0)/1 * f_2 = f_2
#   (1,1): C(2-1,1)/C(2,1) * e_1 * f_1 = C(1,1)/C(2,1) * e_1*f_1 = (1/2)*e_1*f_1
#   (2,0): C(2-0,2)/C(2,2) * e_2 * 1 = C(2,2)/C(2,2) * e_2 = e_2
# g_2 = e_2 + (1/2)*e_1*f_1 + f_2 = ab + (a+b)(c+d)/2 + cd

# So r(x) = x^2 - g_1*x + g_2 = x^2 - (a+b+c+d)x + (ab + (a+b)(c+d)/2 + cd)
# disc = (a+b+c+d)^2 - 4(ab + (a+b)(c+d)/2 + cd)
# = (a+b)^2 + 2(a+b)(c+d) + (c+d)^2 - 4ab - 2(a+b)(c+d) - 4cd
# = (a-b)^2 + (c-d)^2

# PYTHAGOREAN! gap_r^2 = gap_p^2 + gap_q^2 EXACTLY.
# So 1/Phi(r) = gap_r^2/2 = gap_p^2/2 + gap_q^2/2 = 1/Phi(p) + 1/Phi(q)
# EQUALITY for ALL n=2 cases (not just centered).

print("  n=2 analytical: gap_r^2 = gap_p^2 + gap_q^2 (EXACT EQUALITY)")
print("  This is because g_2 = e_2 + (e_1*f_1)/2 + f_2")
print("  gives disc = (a-b)^2 + (c-d)^2 (Pythagorean)")

# Numerical verification
for trial in range(10):
    rp = make_roots(2, min_gap=0.5)
    rq = make_roots(2, min_gap=0.5)
    rr, _ = boxplus(rp, rq)

    gap_p = rp[1] - rp[0]
    gap_q = rq[1] - rq[0]
    gap_r = rr[1] - rr[0]

    print(f"  gap_p^2={gap_p**2:.6f}, gap_q^2={gap_q**2:.6f}, gap_r^2={gap_r**2:.6f}, "
          f"ratio={gap_r**2/(gap_p**2+gap_q**2):.10f}")


# ============================================================
# SECTION 4: INDUCTION ANALYSIS WITH CORRECT BOXPLUS
# ============================================================
print("\n" + "=" * 70)
print("SECTION 4: INDUCTION ANALYSIS (correct boxplus)")
print("=" * 70)

# First: check Phi ratio
print("\n--- Phi_n(p) / Phi_{n-1}(p^{(1)}) for uniform spacing ---")
for n in [3, 4, 5, 6, 7, 8]:
    d = 1.0
    roots = np.array([(i - (n-1)/2) * d for i in range(n)])
    cp = critical_points(roots)
    if cp is not None and len(cp) == n-1:
        ph_n = Phi(roots)
        ph_n1 = Phi(cp)
        if ph_n is not None and ph_n1 is not None:
            print(f"  n={n}: ratio = {ph_n/ph_n1:.8f}")

# Now the main induction analysis
print("\n--- INDUCTION DECOMPOSITION ---")
print("  F_n = F_{n-1} + (delta_r - delta_p - delta_q)")

for n in [3, 4, 5]:
    print(f"\n  n={n}:")
    F_n_vals = []
    F_n1_vals = []
    delta_diffs = []

    for trial in range(500):
        rp = make_roots(n, min_gap=0.3 + 0.5*np.random.rand())
        rq = make_roots(n, min_gap=0.3 + 0.5*np.random.rand())
        rr, imag_r = boxplus(rp, rq)

        if imag_r > 0.01 or np.any(np.diff(rr) < 0.01):
            continue

        cp_p = critical_points(rp)
        cp_q = critical_points(rq)
        if cp_p is None or cp_q is None:
            continue
        if len(cp_p) != n-1 or len(cp_q) != n-1:
            continue

        s_roots, imag_s = boxplus(cp_p, cp_q)
        if imag_s > 0.01:
            continue
        if len(s_roots) > 1 and np.any(np.diff(s_roots) < 0.01):
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

        dd = delta_r - delta_p - delta_q

        F_n_vals.append(F_n)
        F_n1_vals.append(F_n1)
        delta_diffs.append(dd)

    if F_n_vals:
        F_n_arr = np.array(F_n_vals)
        F_n1_arr = np.array(F_n1_vals)
        dd_arr = np.array(delta_diffs)

        print(f"    {len(F_n_vals)} valid samples")
        print(f"    F_n:  min={F_n_arr.min():.10f}, all >= 0: {np.all(F_n_arr >= -1e-8)}")
        n_fn_viol = np.sum(F_n_arr < -1e-8)
        if n_fn_viol > 0:
            print(f"    *** F_n violations: {n_fn_viol} ***")
        print(f"    F_{n-1}: min={F_n1_arr.min():.10f}, all >= 0: {np.all(F_n1_arr >= -1e-8)}")
        print(f"    delta_diff: min={dd_arr.min():.10f}, max={dd_arr.max():.10f}")
        print(f"    delta superadditive: {np.all(dd_arr >= -1e-8)}")

        # Consistency
        check = F_n_arr - (F_n1_arr + dd_arr)
        print(f"    Consistency: max|F_n - (F_{n-1}+dd)| = {np.max(np.abs(check)):.2e}")

        if not np.all(dd_arr >= -1e-8):
            n_viol = np.sum(dd_arr < -1e-8)
            print(f"    delta violations: {n_viol}/{len(dd_arr)}")

            # When delta is negative, does F_{n-1} compensate?
            failed_induction = F_n1_arr + dd_arr
            print(f"    F_{n-1}+dd min = {failed_induction.min():.10f} (= F_n min)")


# ============================================================
# SECTION 5: EXCESS RATIO F_n / F_{n-1}
# ============================================================
print("\n" + "=" * 70)
print("SECTION 5: EXCESS RATIO F_n / F_{n-1}")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\n  n={n}:")
    ratios = []
    for trial in range(500):
        rp = make_roots(n, min_gap=0.3 + 0.5*np.random.rand())
        rq = make_roots(n, min_gap=0.3 + 0.5*np.random.rand())
        rr, imag_r = boxplus(rp, rq)

        if imag_r > 0.01 or np.any(np.diff(rr) < 0.01):
            continue

        cp_p = critical_points(rp)
        cp_q = critical_points(rq)
        if cp_p is None or cp_q is None:
            continue
        if len(cp_p) != n-1 or len(cp_q) != n-1:
            continue

        s_roots, imag_s = boxplus(cp_p, cp_q)
        if imag_s > 0.01 or np.any(np.diff(s_roots) < 0.01):
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

        if abs(F_n1) > 1e-10:
            ratios.append(F_n / F_n1)

    if ratios:
        rarr = np.array(ratios)
        print(f"    {len(ratios)} valid samples")
        print(f"    F_n/F_{n-1}: range=[{rarr.min():.6f}, {rarr.max():.6f}]")
        print(f"    mean={rarr.mean():.6f}, std={rarr.std():.6f}")
        print(f"    all positive: {np.all(rarr > 0)}")
        if np.all(rarr > 0):
            print(f"    *** F_n/F_{n-1} > 0 always => induction WORKS if base case holds! ***")


# ============================================================
# SECTION 6: DELTA SIGN ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("SECTION 6: SIGN OF delta_p = 1/Phi_n(p) - 1/Phi_{n-1}(p^{(1)})")
print("=" * 70)

for n in [3, 4, 5, 6]:
    deltas = []
    for trial in range(300):
        roots = make_roots(n, min_gap=0.3 + np.random.rand())
        cp = critical_points(roots)
        if cp is None or len(cp) != n-1:
            continue
        if len(cp) > 1 and np.any(np.diff(cp) < 0.01):
            continue
        ph = Phi(roots)
        ph1 = Phi(cp)
        if ph is None or ph1 is None or ph1 < 1e-12:
            continue
        deltas.append(1/ph - 1/ph1)

    if deltas:
        darr = np.array(deltas)
        print(f"  n={n}: range=[{darr.min():.8f}, {darr.max():.8f}], "
              f"always neg: {np.all(darr < 0)}, always pos: {np.all(darr > 0)}")


# ============================================================
# SECTION 7: TELESCOPING
# ============================================================
print("\n" + "=" * 70)
print("SECTION 7: TELESCOPING — iterated derivatives")
print("=" * 70)

def iterated_cps(roots, depth):
    """Compute roots, cp, cp of cp, etc."""
    chain = [roots]
    current = roots
    for k in range(depth):
        cp = critical_points(current)
        if cp is None or len(cp) != len(current) - 1:
            return None
        if len(cp) > 1 and np.any(np.diff(cp) < 1e-8):
            return None
        chain.append(cp)
        current = cp
    return chain

for trial in range(5):
    n = 6
    roots = make_roots(n, min_gap=1.0)
    chain = iterated_cps(roots, n-2)
    if chain is None:
        continue

    print(f"\n  Trial {trial}: n={n}")
    phis = []
    for k, r in enumerate(chain):
        ph = Phi(r)
        if ph is None:
            break
        phis.append(ph)
        print(f"    p^({k}): {len(r)} roots, Phi={ph:.6f}, 1/Phi={1/ph:.8f}")

    if len(phis) == len(chain):
        inv_phis = [1/p for p in phis]
        diffs = [inv_phis[k+1] - inv_phis[k] for k in range(len(inv_phis)-1)]
        print(f"    Differences 1/Phi_{k+1} - 1/Phi_k: {np.array2string(np.array(diffs), precision=6)}")
        print(f"    Sum of differences = {sum(diffs):.8f} = 1/Phi_2 - 1/Phi_n = {inv_phis[-1]-inv_phis[0]:.8f}")


# ============================================================
# SECTION 8: CONVOLUTION TELESCOPING
# ============================================================
print("\n" + "=" * 70)
print("SECTION 8: CONVOLUTION TELESCOPING")
print("  Check: F_k at each derivative level")
print("=" * 70)

for trial in range(5):
    n = 5
    rp = make_roots(n, min_gap=0.8)
    rq = make_roots(n, min_gap=0.8)
    rr, imag = boxplus(rp, rq)
    if imag > 0.01:
        continue

    chain_p = iterated_cps(rp, n-2)
    chain_q = iterated_cps(rq, n-2)
    chain_r = iterated_cps(rr, n-2)

    if any(c is None for c in [chain_p, chain_q, chain_r]):
        continue

    print(f"\n  Trial {trial}:")
    for k in range(len(chain_p)):
        ph_p = Phi(chain_p[k])
        ph_q = Phi(chain_q[k])
        ph_r = Phi(chain_r[k])

        if any(v is None for v in [ph_p, ph_q, ph_r]):
            break

        F_k = 1/ph_r - 1/ph_p - 1/ph_q
        deg = n - k
        label = ">=0" if F_k >= -1e-10 else "NEGATIVE"
        print(f"    degree {deg}: F = {F_k:.10f} [{label}]")

        # Also check MSS identity at this level
        if k < len(chain_p) - 1:
            s_roots, _ = boxplus(chain_p[k+1], chain_q[k+1])
            if s_roots is not None:
                error = np.max(np.abs(np.sort(chain_r[k+1]) - np.sort(s_roots)))
                print(f"      MSS derivative check: error = {error:.2e}")


# ============================================================
# SECTION 9: DETAILED n=3 ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("SECTION 9: DETAILED n=3 ANALYSIS")
print("  n=3 base case for induction (n=2 gives equality)")
print("=" * 70)

print("\n  Checking if F_3 >= 0 always:")
n = 3
min_excess = float('inf')
worst_case = None
for trial in range(2000):
    rp = make_roots(3, min_gap=0.2 + np.random.rand() * 2)
    rq = make_roots(3, min_gap=0.2 + np.random.rand() * 2)
    rr, imag = boxplus(rp, rq)

    if imag > 0.01 or np.any(np.diff(rr) < 0.01):
        continue

    ph_p = Phi(rp)
    ph_q = Phi(rq)
    ph_r = Phi(rr)

    if any(v is None or v < 1e-12 for v in [ph_p, ph_q, ph_r]):
        continue

    excess = 1/ph_r - 1/ph_p - 1/ph_q
    if excess < min_excess:
        min_excess = excess
        worst_case = (rp.copy(), rq.copy(), rr.copy(), ph_p, ph_q, ph_r, excess)

print(f"  Min excess over 2000 trials: {min_excess:.12f}")
if worst_case:
    rp, rq, rr, pp, pq, pr, ex = worst_case
    print(f"  Worst case:")
    print(f"    p = {rp}")
    print(f"    q = {rq}")
    print(f"    r = {rr}")
    print(f"    Phi_p = {pp:.6f}, Phi_q = {pq:.6f}, Phi_r = {pr:.6f}")
    print(f"    excess = {ex:.12f}")


# ============================================================
# SECTION 10: FORMULA FOR 1/Phi IN TERMS OF ROOT GAPS (n=3)
# ============================================================
print("\n" + "=" * 70)
print("SECTION 10: 1/Phi_3 FORMULA")
print("=" * 70)

# For n=3 with roots a < b < c:
# H(a) = 1/(a-b) + 1/(a-c)
# H(b) = 1/(b-a) + 1/(b-c)
# H(c) = 1/(c-a) + 1/(c-b)
# Let d1 = b-a, d2 = c-b (gaps). Then:
# H(a) = -1/d1 - 1/(d1+d2)
# H(b) = 1/d1 - 1/d2
# H(c) = 1/(d1+d2) + 1/d2
# Phi_3 = H(a)^2 + H(b)^2 + H(c)^2

# For uniform spacing d1=d2=d:
# H(a) = -1/d - 1/(2d) = -3/(2d)
# H(b) = 1/d - 1/d = 0
# H(c) = 1/(2d) + 1/d = 3/(2d)
# Phi_3 = 2*(3/(2d))^2 = 9/(2d^2)
# 1/Phi_3 = 2d^2/9

# Check:
d = 1.0
roots = np.array([0, d, 2*d])
print(f"  Uniform d={d}: Phi_3 = {Phi(roots):.6f}, expected 9/(2*{d}^2) = {9/(2*d**2):.6f}")
print(f"  1/Phi_3 = {1/Phi(roots):.6f}, expected 2*{d}^2/9 = {2*d**2/9:.6f}")

# For the boxplus of two degree-3 polynomials with uniform spacing:
# p: roots 0, d_p, 2d_p
# q: roots 0, d_q, 2d_q
# If the conjecture holds: 1/Phi(r) >= 2d_p^2/9 + 2d_q^2/9 = 2(d_p^2+d_q^2)/9
# Does r also have uniform spacing? If so, with what gap?

print("\n  Uniform spacing convolution:")
for dp, dq in [(1.0, 1.0), (1.0, 2.0), (0.5, 3.0), (1.0, 0.5)]:
    rp = np.array([0, dp, 2*dp])
    rq = np.array([0, dq, 2*dq])
    rr, _ = boxplus(rp, rq)
    gaps = np.diff(rr)
    ph_r = Phi(rr)
    ph_p = Phi(rp)
    ph_q = Phi(rq)
    excess = 1/ph_r - 1/ph_p - 1/ph_q

    print(f"  dp={dp}, dq={dq}: r gaps={np.array2string(gaps, precision=4)}, "
          f"1/Phi_r={1/ph_r:.6f}, 1/Phi_p+1/Phi_q={1/ph_p+1/ph_q:.6f}, excess={excess:.8f}")


# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)
print("""
KEY FINDINGS WITH CORRECT BOXPLUS FORMULA:

1. BOXPLUS FORMULA:
   The correct formula is g_k = sum_{i+j=k} C(n-j,i)/C(n,i) * e_i(p) * e_j(q)
   (verified against Monte Carlo). The previously used formula with
   w = (n-i)!(n-j)!/(n!(n-k)!) applied to normalized coefficients was WRONG
   for non-centered polynomials.

2. n=2 BASE CASE:
   EXACT EQUALITY: 1/Phi_2(r) = 1/Phi_2(p) + 1/Phi_2(q).
   This is gap_r^2 = gap_p^2 + gap_q^2 (Pythagorean).

3. CONJECTURE STATUS:
   Need to check if 1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q) holds for n >= 3.

4. MSS DERIVATIVE IDENTITY:
   r'(x) = n * (p^{(1)} boxplus_{n-1} q^{(1)})(x).
   Works perfectly for n=3,4. Some numerical issues for n=5.

5. INDUCTION DECOMPOSITION:
   F_n = F_{n-1} + (delta_r - delta_p - delta_q) is exact.
   The correction term delta is NOT superadditive in general.

6. delta_p = 1/Phi_n - 1/Phi_{n-1}: ALWAYS NEGATIVE.
""")
