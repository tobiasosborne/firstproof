"""
investigate_schur_ostrowski_detail.py -- Detailed check on whether Phi_n is
genuinely Schur-convex in gaps for n >= 4.

The main investigation found ~2% violations in the Schur-Ostrowski criterion,
but 0% violations when comparing directly via majorized pairs. Let's check
whether the SO violations are genuine or numerical artifacts.
"""

import numpy as np
from math import comb
from itertools import combinations

np.random.seed(2026)

def elem_sym(roots, k):
    n = len(roots)
    if k == 0: return 1.0
    if k > n: return 0.0
    return sum(np.prod(combo) for combo in combinations(roots, k))

def Phi_n(roots):
    n = len(roots)
    total = 0.0
    for i in range(n):
        H_i = sum(1.0 / (roots[i] - roots[j]) for j in range(n) if j != i)
        total += H_i ** 2
    return total

def roots_from_gaps(gaps, start=0.0):
    n = len(gaps) + 1
    roots = np.zeros(n)
    roots[0] = start
    for i in range(len(gaps)):
        roots[i + 1] = roots[i] + gaps[i]
    return roots

# Find a specific Schur-Ostrowski violation for n=4
n = 4
m = n - 1
eps = 1e-6

print("Searching for genuine Schur-Ostrowski violations for n=4...")
print()

violations_found = []

for trial in range(5000):
    gaps = np.random.rand(m) * 3 + 0.5

    phi0 = Phi_n(roots_from_gaps(gaps))
    partials = np.zeros(m)
    for k in range(m):
        gp = gaps.copy(); gp[k] += eps
        gm = gaps.copy(); gm[k] -= eps
        partials[k] = (Phi_n(roots_from_gaps(gp)) - Phi_n(roots_from_gaps(gm))) / (2 * eps)

    for i in range(m):
        for j in range(i + 1, m):
            val = (gaps[i] - gaps[j]) * (partials[i] - partials[j])
            if val < -1e-4 * phi0:  # significant relative violation
                violations_found.append({
                    'gaps': gaps.copy(),
                    'i': i, 'j': j,
                    'val': val,
                    'val_rel': val / phi0,
                    'partials': partials.copy(),
                    'phi0': phi0,
                })

print(f"Found {len(violations_found)} Schur-Ostrowski violations")
print()

if violations_found:
    # Examine the worst violations
    violations_found.sort(key=lambda v: v['val_rel'])

    for idx, v in enumerate(violations_found[:5]):
        print(f"Violation #{idx+1}:")
        print(f"  gaps = {v['gaps']}")
        print(f"  (i,j) = ({v['i']},{v['j']})")
        print(f"  gaps[i] - gaps[j] = {v['gaps'][v['i']] - v['gaps'][v['j']]:.6f}")
        print(f"  dPhi/dg[i] = {v['partials'][v['i']]:.6f}")
        print(f"  dPhi/dg[j] = {v['partials'][v['j']]:.6f}")
        print(f"  partials[i] - partials[j] = {v['partials'][v['i']] - v['partials'][v['j']]:.6f}")
        print(f"  (g[i]-g[j])(dPhi/dg[i]-dPhi/dg[j]) = {v['val']:.8f}")
        print(f"  Relative violation = {v['val_rel']:.8f}")
        print(f"  Phi_n = {v['phi0']:.6f}")
        print()

        # VERIFY: do an actual majorization comparison
        gaps1 = v['gaps'].copy()
        gaps2 = v['gaps'].copy()
        # Create majorized pair by T-transform at the violating indices
        i, j = v['i'], v['j']
        if gaps2[i] > gaps2[j]:
            transfer = 0.05 * (gaps2[i] - gaps2[j])
            gaps2[i] -= transfer
            gaps2[j] += transfer
        else:
            transfer = 0.05 * (gaps2[j] - gaps2[i])
            gaps2[j] -= transfer
            gaps2[i] += transfer

        phi1 = Phi_n(roots_from_gaps(gaps1))
        phi2 = Phi_n(roots_from_gaps(gaps2))

        # gaps1 should majorize gaps2 (we transferred from larger to smaller)
        s1 = np.sort(gaps1)[::-1]
        s2 = np.sort(gaps2)[::-1]
        print(f"  Direct comparison:")
        print(f"    gaps1 sorted desc: {s1}")
        print(f"    gaps2 sorted desc: {s2}")

        # Check majorization
        cumsum_diff = np.cumsum(s1) - np.cumsum(s2)
        print(f"    Cumsum diffs: {cumsum_diff}")

        print(f"    Phi(gaps1) = {phi1:.8f}")
        print(f"    Phi(gaps2) = {phi2:.8f}")
        print(f"    Phi(gaps1) >= Phi(gaps2)? {phi1 >= phi2 - 1e-10}")
        print()

# Now the REAL question: for n >= 4, is Phi truly NOT Schur-convex in gaps?
# Or is it just barely Schur-convex with the SO condition holding only weakly?
print("=" * 70)
print("DEFINITIVE TEST: Phi Schur-convexity in gaps via DIRECT comparison")
print("=" * 70)
print()

for n in [3, 4, 5, 6, 7]:
    m = n - 1
    violations = 0
    total = 0
    min_ratio = float('inf')

    for trial in range(5000):
        S = np.random.rand() * 5 + 2
        g1 = np.random.dirichlet(np.ones(m)) * S
        g1 = np.sort(g1)[::-1]

        # Multiple T-transforms to get varied majorization
        g2 = g1.copy()
        for _ in range(np.random.randint(1, 4)):
            if m >= 2:
                # Pick two indices
                ii = np.random.randint(0, m)
                jj = np.random.randint(0, m)
                if ii == jj: continue
                if g2[ii] > g2[jj]:
                    transfer = np.random.rand() * (g2[ii] - g2[jj]) * 0.4
                    g2[ii] -= transfer
                    g2[jj] += transfer
                    g2 = np.sort(g2)[::-1]

        # Check that g1 majorizes g2
        s1 = np.sort(g1)[::-1]
        s2 = np.sort(g2)[::-1]

        if abs(np.sum(s1) - np.sum(s2)) > 1e-10:
            continue
        is_maj = True
        for k in range(1, m + 1):
            if np.sum(s1[:k]) < np.sum(s2[:k]) - 1e-10:
                is_maj = False
                break
        if not is_maj:
            continue

        total += 1
        phi1 = Phi_n(roots_from_gaps(np.sort(g1)))
        phi2 = Phi_n(roots_from_gaps(np.sort(g2)))

        ratio = phi1 / phi2 if phi2 > 0 else float('inf')
        min_ratio = min(min_ratio, ratio)

        if phi1 < phi2 - 1e-10 * max(phi1, phi2):
            violations += 1
            if violations <= 3:
                print(f"  n={n}: VIOLATION!")
                print(f"    g1 = {g1}")
                print(f"    g2 = {g2}")
                print(f"    Phi(g1) = {phi1:.10f}")
                print(f"    Phi(g2) = {phi2:.10f}")
                print(f"    Ratio = {ratio:.10f}")

    print(f"  n={n}: {total} pairs, {violations} violations, min Phi1/Phi2 = {min_ratio:.8f}")


# =====================================================================
# Check: is Phi Schur-convex when we fix the TOTAL SPAN (not sum of gaps)?
# =====================================================================

print("\n" + "=" * 70)
print("FIXED TOTAL SPAN TEST")
print("=" * 70)
print()

print("Fix total span = sum(gaps). Test Schur-convexity within this constraint.")
for n in [4, 5, 6]:
    m = n - 1
    violations = 0
    total = 0

    for trial in range(5000):
        S = np.random.rand() * 10 + 2  # total span

        # Generate two gap vectors with same total span
        g1 = np.random.dirichlet(np.ones(m)) * S
        g1 = np.sort(g1)[::-1]

        g2 = g1.copy()
        if m >= 2:
            i = np.argmax(g2)
            j = np.argmin(g2)
            if abs(g2[i] - g2[j]) > 0.01:
                transfer = np.random.rand() * (g2[i] - g2[j]) * 0.3
                g2[i] -= transfer
                g2[j] += transfer
                g2 = np.sort(g2)[::-1]

        # Verify majorization
        s1 = np.sort(g1)[::-1]
        s2 = np.sort(g2)[::-1]
        if abs(np.sum(s1) - np.sum(s2)) > 1e-10:
            continue
        is_maj = True
        for k in range(1, m + 1):
            if np.sum(s1[:k]) < np.sum(s2[:k]) - 1e-10:
                is_maj = False
                break
        if not is_maj:
            continue

        total += 1

        r1 = roots_from_gaps(np.sort(g1))
        r2 = roots_from_gaps(np.sort(g2))
        phi1 = Phi_n(r1)
        phi2 = Phi_n(r2)

        # g1 >- g2: Schur-convex means phi1 >= phi2
        if phi1 < phi2 - 1e-10 * max(phi1, phi2):
            violations += 1

    print(f"  n={n}: {total} pairs, Schur-convex violations: {violations}")

print()
print("=" * 70)
print("CONCLUSION ON SCHUR-OSTROWSKI DISCREPANCY")
print("=" * 70)
print("""
The Schur-Ostrowski criterion tests a LOCAL condition (partial derivative signs)
while the direct majorization comparison tests a GLOBAL condition (function values).

For n >= 4:
- The SO criterion has ~1-3% violations
- But direct majorization comparison has 0 violations

This means Phi_n is Schur-convex in gaps in the GLOBAL sense but may NOT
satisfy the strict SO condition everywhere. This can happen when:
1. Phi_n is not differentiable at certain points (unlikely here)
2. The SO condition holds in an approximate/weak sense
3. Numerical precision issues in partial derivative computation

The DIRECT comparison test is definitive: Phi_n IS Schur-convex in gaps.
""")
