"""
investigate_AplusB.py -- Does A + B >= 2*Phi_r hold?

If true, then by AM-GM: AB >= ((A+B)/2)^2 >= Phi_r^2 = ||h||^4, closing the conjecture.

Setup:
  p, q: monic real-rooted degree n, simple roots.
  r = p boxplus_n q (MSS finite free additive convolution).
  h_k = H_r(nu_k), u_k = H_p(lambda_k), v_k = H_q(mu_k)
  Phi_f = ||H_f||^2
  A = Phi_p - Phi_r = ||u||^2 - ||h||^2
  B = Phi_q - Phi_r = ||v||^2 - ||h||^2
  Target: A + B >= 2*Phi_r, equivalently Phi_p + Phi_q >= 4*Phi_r.

We use TWO boxplus formulas and verify against Monte Carlo.
"""

import numpy as np
from math import factorial, comb
from itertools import combinations
from scipy.optimize import minimize, differential_evolution
import warnings
warnings.filterwarnings('ignore')

# ================================================================
# BOXPLUS FORMULAS
# ================================================================

def elem_sym(roots, k):
    """Elementary symmetric polynomial e_k(roots)."""
    n = len(roots)
    if k == 0: return 1.0
    if k > n: return 0.0
    return sum(np.prod([roots[i] for i in subset])
               for subset in combinations(range(n), k))

def boxplus_formula1(roots_p, roots_q):
    """Original MSS coefficient formula:
    c_k = sum_{i+j=k} [(n-i)!(n-j)! / (n!(n-k)!)] * a_i * b_j
    where p(x) = sum a_k x^{n-k} (standard poly coefficients).
    """
    n = len(roots_p)
    a = np.real(np.poly(roots_p))  # [1, -e1, e2, ..., (-1)^n e_n]
    b = np.real(np.poly(roots_q))

    c = np.zeros(n+1)
    for k in range(n+1):
        for i in range(k+1):
            j = k - i
            if i <= n and j <= n:
                w = factorial(n-i) * factorial(n-j) / (factorial(n) * factorial(n-k))
                c[k] += w * a[i] * b[j]

    roots_r = np.sort(np.real(np.roots(c)))
    return roots_r, c

def boxplus_formula2(roots_p, roots_q):
    """Marcus 2021 formula using elementary symmetric polynomials:
    g_k = sum_{i+j=k} [C(n-j, i) / C(n, i)] * e_i(p) * e_j(q)
    where e_k = (-1)^k * poly_coeff[k].
    """
    n = len(roots_p)
    poly_p = np.real(np.poly(roots_p))
    poly_q = np.real(np.poly(roots_q))

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

    roots_r = np.sort(np.real(np.roots(poly_r)))
    return roots_r, poly_r

def mc_boxplus(roots_p, roots_q, N=50000):
    """Monte Carlo: E[char poly of A + UBU*] with Haar unitary U."""
    n = len(roots_p)
    A = np.diag(roots_p.astype(complex))
    B = np.diag(roots_q.astype(complex))

    poly_sum = np.zeros(n+1)
    for _ in range(N):
        Z = (np.random.randn(n, n) + 1j * np.random.randn(n, n)) / np.sqrt(2)
        Q, R = np.linalg.qr(Z)
        dp = np.diag(R)
        dp = dp / np.abs(dp)
        U = Q @ np.diag(dp)
        M = A + U @ B @ U.conj().T
        poly_sum += np.real(np.poly(M))

    return poly_sum / N

# ================================================================
# CORE ANALYSIS FUNCTIONS
# ================================================================

def H_values(roots):
    """H_f(root_k) = sum_{j!=k} 1/(root_k - root_j)."""
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H

def compute_quantities(roots_p, roots_q, use_formula=2):
    """Compute A, B, Phi_r and check A + B >= 2*Phi_r.

    use_formula: 1 = original MSS, 2 = Marcus 2021 (default).
    Returns dict or None if invalid.
    """
    n = len(roots_p)

    if use_formula == 1:
        roots_r, poly_r = boxplus_formula1(roots_p, roots_q)
    else:
        roots_r, poly_r = boxplus_formula2(roots_p, roots_q)

    # Check real-rootedness
    raw_roots = np.roots(poly_r)
    if np.any(np.abs(np.imag(raw_roots)) > 1e-6):
        return None
    roots_r = np.sort(np.real(raw_roots))

    # Check distinct roots
    if np.any(np.diff(roots_r) < 1e-6):
        return None

    h = H_values(roots_r)
    u = H_values(roots_p)
    v = H_values(roots_q)

    Phi_r = np.dot(h, h)
    Phi_p = np.dot(u, u)
    Phi_q = np.dot(v, v)

    A = Phi_p - Phi_r
    B = Phi_q - Phi_r

    if A <= -1e-10 or B <= -1e-10 or Phi_r <= 0:
        return None

    return {
        'n': n,
        'roots_p': roots_p, 'roots_q': roots_q, 'roots_r': roots_r,
        'Phi_p': Phi_p, 'Phi_q': Phi_q, 'Phi_r': Phi_r,
        'A': A, 'B': B,
        'ApB': A + B,                     # A + B
        'two_Phi_r': 2 * Phi_r,           # target lower bound
        'ratio_ApB': (A + B) / (2 * Phi_r),  # should be >= 1
        'AB': A * B,
        'h4': Phi_r**2,
        'AB_ratio': A * B / Phi_r**2,     # AB / h^4, should be >= 1
    }

def random_roots(n, scale=2.0, min_gap=0.3):
    """Generate random sorted roots with minimum gap."""
    roots = np.sort(np.random.randn(n) * scale)
    for i in range(1, n):
        if roots[i] - roots[i-1] < min_gap:
            roots[i] = roots[i-1] + min_gap
    return roots


# ================================================================
# PART 0: VERIFY FORMULAS AGAINST MONTE CARLO
# ================================================================

print("=" * 70)
print("PART 0: VERIFY BOXPLUS FORMULAS AGAINST MONTE CARLO")
print("=" * 70)

np.random.seed(42)

test_cases = [
    ("centered n=2", np.array([-1.0, 1.0]), np.array([-1.0, 1.0])),
    ("non-centered n=2", np.array([0.0, 1.0]), np.array([-3.0, 0.0])),
    ("centered n=3", np.array([-1.0, 0.0, 1.0]), np.array([-2.0, 0.0, 2.0])),
    ("non-centered n=3", np.array([0.0, 1.0, 3.0]), np.array([-1.0, 0.0, 2.0])),
    ("non-centered n=4", np.array([0.0, 1.0, 3.0, 5.0]), np.array([-2.0, -1.0, 0.0, 2.0])),
]

for name, rp, rq in test_cases:
    mc_poly = mc_boxplus(rp, rq, N=50000)
    _, f1_poly = boxplus_formula1(rp, rq)
    _, f2_poly = boxplus_formula2(rp, rq)

    err1 = np.max(np.abs(f1_poly - mc_poly))
    err2 = np.max(np.abs(f2_poly - mc_poly))

    print(f"  {name:25s} | F1 err: {err1:.6f} | F2 err: {err2:.6f} | "
          f"{'F1' if err1 < err2 else 'F2'} better")

print("\n  F1 = original MSS coefficient formula")
print("  F2 = Marcus 2021 (elementary symmetric)")
print("  (Errors > 0.01 suggest wrong formula for non-centered case)")


# ================================================================
# PART 1: NUMERICAL TEST OF A + B >= 2*Phi_r
# ================================================================

print("\n" + "=" * 70)
print("PART 1: A + B >= 2*Phi_r NUMERICAL TEST")
print("=" * 70)

np.random.seed(123)

for formula_id in [1, 2]:
    fname = "F1 (MSS coeff)" if formula_id == 1 else "F2 (Marcus 2021)"
    print(f"\n--- Using {fname} ---")

    for n in [2, 3, 4, 5, 6]:
        n_trials = 600
        violations = 0
        min_ratio = float('inf')
        max_ratio = 0
        valid = 0
        ratios_list = []

        for trial in range(n_trials):
            rp = random_roots(n)
            rq = random_roots(n)

            result = compute_quantities(rp, rq, use_formula=formula_id)
            if result is None:
                continue

            valid += 1
            r = result['ratio_ApB']
            ratios_list.append(r)
            min_ratio = min(min_ratio, r)
            max_ratio = max(max_ratio, r)

            if r < 1.0 - 1e-8:
                violations += 1

        if valid > 0:
            print(f"  n={n}: valid={valid:4d}, violations={violations}, "
                  f"min(A+B)/(2*Phi_r)={min_ratio:.8f}, "
                  f"max={max_ratio:.4f}, "
                  f"mean={np.mean(ratios_list):.4f}")
        else:
            print(f"  n={n}: no valid trials")


# ================================================================
# PART 2: ANALYTIC PROOF FOR n=2
# ================================================================

print("\n" + "=" * 70)
print("PART 2: ANALYTIC PROOF FOR n=2")
print("=" * 70)

print("""
For n=2 with p having roots {a, b} and q having roots {c, d}:
  gap_p = |a - b|, gap_q = |c - d|
  Phi_2(p) = 2 / gap_p^2  (since H = [1/(a-b), 1/(b-a)] and ||H||^2 = 2/(a-b)^2)

For r = p boxplus_2 q:
  r has roots {nu_1, nu_2} with gap_r satisfying gap_r^2 = gap_p^2 + gap_q^2
  (This is the Pythagorean property of n=2 boxplus.)
  Phi_2(r) = 2 / gap_r^2 = 2 / (gap_p^2 + gap_q^2)

So Phi_p + Phi_q >= 4*Phi_r becomes:
  2/gap_p^2 + 2/gap_q^2 >= 8/(gap_p^2 + gap_q^2)

Let x = gap_p^2 > 0, y = gap_q^2 > 0. We need:
  (1/x + 1/y) >= 4/(x + y)
  (x + y)/(xy) >= 4/(x + y)
  (x + y)^2 >= 4*xy
  x^2 + 2xy + y^2 >= 4xy
  x^2 - 2xy + y^2 >= 0
  (x - y)^2 >= 0   ALWAYS TRUE.

Equality iff x = y, i.e., gap_p = gap_q.

PROVED: Phi_p + Phi_q >= 4*Phi_r for n=2, with equality iff gap_p = gap_q.

For n=2, A + B = Phi_p + Phi_q - 2*Phi_r
  = 2/x + 2/y - 4/(x+y)
  = 2(x+y)/(xy) - 4/(x+y)
  = [2(x+y)^2 - 4xy] / [xy(x+y)]
  = 2[(x+y)^2 - 2xy] / [xy(x+y)]
  = 2[x^2 + y^2] / [xy(x+y)]
  >= 0 always.

Also: 2*Phi_r = 4/(x+y).

A + B >= 2*Phi_r <=> 2(x^2+y^2)/[xy(x+y)] >= 4/(x+y)
<=> (x^2+y^2)/(xy) >= 2 <=> (x-y)^2 >= 0. QED.
""")

# Verify numerically
print("Numerical verification for n=2:")
for trial in range(5):
    gp = np.random.uniform(0.5, 5.0)
    gq = np.random.uniform(0.5, 5.0)
    rp = np.array([-gp/2, gp/2])
    rq = np.array([-gq/2, gq/2])

    result = compute_quantities(rp, rq, use_formula=2)
    if result:
        x, y = gp**2, gq**2
        analytic_ratio = (x**2 + y**2) / (2*x*y)
        print(f"  gap_p={gp:.3f}, gap_q={gq:.3f}: "
              f"(A+B)/(2Phi_r) = {result['ratio_ApB']:.8f}, "
              f"analytic = {analytic_ratio:.8f}, "
              f"match = {abs(result['ratio_ApB'] - analytic_ratio) < 1e-6}")


# ================================================================
# PART 2b: VERIFY PYTHAGOREAN PROPERTY FOR n=2
# ================================================================

print("\n--- Verifying Pythagorean property gap_r^2 = gap_p^2 + gap_q^2 for n=2 ---")
for trial in range(5):
    rp = random_roots(2, scale=3.0)
    rq = random_roots(2, scale=3.0)

    for fid in [1, 2]:
        roots_r, _ = (boxplus_formula1 if fid == 1 else boxplus_formula2)(rp, rq)
        gp2 = (rp[1] - rp[0])**2
        gq2 = (rq[1] - rq[0])**2
        gr2 = (roots_r[1] - roots_r[0])**2
        print(f"  F{fid}: gap_p^2={gp2:.4f}, gap_q^2={gq2:.4f}, "
              f"gap_r^2={gr2:.4f}, gap_p^2+gap_q^2={gp2+gq2:.4f}, "
              f"match={abs(gr2 - gp2 - gq2) < 0.01}")


# ================================================================
# PART 3: TIGHTNESS FOR n >= 3
# ================================================================

print("\n" + "=" * 70)
print("PART 3: TIGHTNESS ANALYSIS - min (A+B)/(2*Phi_r) FOR n >= 3")
print("=" * 70)

np.random.seed(456)

for n in [3, 4, 5, 6]:
    min_ratio = float('inf')
    min_case = None
    ratios = []

    for trial in range(2000):
        rp = random_roots(n, scale=2.0, min_gap=0.2)
        rq = random_roots(n, scale=2.0, min_gap=0.2)

        result = compute_quantities(rp, rq, use_formula=2)
        if result is None:
            continue

        r = result['ratio_ApB']
        ratios.append(r)
        if r < min_ratio:
            min_ratio = r
            min_case = result

    print(f"\n  n={n}: {len(ratios)} valid trials")
    print(f"    min (A+B)/(2*Phi_r) = {min_ratio:.8f}")
    print(f"    5th percentile       = {np.percentile(ratios, 5):.6f}")
    print(f"    mean                 = {np.mean(ratios):.4f}")
    print(f"    max                  = {np.max(ratios):.4f}")

    if min_case:
        print(f"    At minimum: Phi_p={min_case['Phi_p']:.4f}, "
              f"Phi_q={min_case['Phi_q']:.4f}, Phi_r={min_case['Phi_r']:.4f}")
        # Check if roots are approximately equally spaced
        gaps_p = np.diff(min_case['roots_p'])
        gaps_q = np.diff(min_case['roots_q'])
        print(f"    p gaps: {gaps_p}")
        print(f"    q gaps: {gaps_q}")
        gap_cv_p = np.std(gaps_p) / np.mean(gaps_p) if np.mean(gaps_p) > 0 else 0
        gap_cv_q = np.std(gaps_q) / np.mean(gaps_q) if np.mean(gaps_q) > 0 else 0
        print(f"    gap CV (p) = {gap_cv_p:.4f}, gap CV (q) = {gap_cv_q:.4f}")


# ================================================================
# PART 3b: TARGETED OPTIMIZATION to find minimum
# ================================================================

print("\n--- Optimization search for minimum (A+B)/(2*Phi_r) ---")

def neg_ratio_ApB(params, n):
    """Negative of (A+B)/(2*Phi_r) for optimization (minimize this)."""
    half = n
    gaps_p = np.exp(params[:n-1])  # ensure positive gaps
    gaps_q = np.exp(params[n-1:2*(n-1)])

    rp = np.cumsum(np.concatenate([[0], gaps_p]))
    rq = np.cumsum(np.concatenate([[0], gaps_q]))

    result = compute_quantities(rp, rq, use_formula=2)
    if result is None:
        return 0  # penalty: return 0 (we're minimizing negative, so 0 is bad)

    return -result['ratio_ApB']

for n in [3, 4, 5]:
    best_ratio = float('inf')

    for attempt in range(20):
        x0 = np.random.randn(2*(n-1)) * 0.5

        try:
            res = minimize(neg_ratio_ApB, x0, args=(n,), method='Nelder-Mead',
                          options={'maxiter': 5000, 'xatol': 1e-10, 'fatol': 1e-12})
            ratio = -res.fun
            if 0 < ratio < best_ratio:
                best_ratio = ratio
        except:
            pass

    print(f"  n={n}: min (A+B)/(2*Phi_r) approx = {best_ratio:.10f}")


# ================================================================
# PART 4: EQUALLY SPACED ROOTS -- EXACT FORMULA
# ================================================================

print("\n" + "=" * 70)
print("PART 4: EQUALLY SPACED ROOTS")
print("=" * 70)

print("\nFor equally spaced p = {0, s, 2s, ...} and q = {0, t, 2t, ...}:")

for n in [2, 3, 4, 5, 6, 7, 8]:
    for s, t in [(1.0, 1.0), (1.0, 2.0), (1.0, 0.5), (2.0, 3.0)]:
        rp = np.arange(n, dtype=float) * s
        rq = np.arange(n, dtype=float) * t

        result = compute_quantities(rp, rq, use_formula=2)
        if result:
            print(f"  n={n}, s={s:.1f}, t={t:.1f}: "
                  f"(A+B)/(2Phi_r) = {result['ratio_ApB']:.8f}, "
                  f"AB/h^4 = {result['AB_ratio']:.8f}")


# ================================================================
# PART 5: THE AM-GM PROOF CHAIN
# ================================================================

print("\n" + "=" * 70)
print("PART 5: FULL PROOF CHAIN VERIFICATION")
print("=" * 70)

print("""
If A + B >= 2*Phi_r (which we are testing), then:

Step 1: A, B > 0 (already established: Phi_p > Phi_r, Phi_q > Phi_r)
Step 2: By AM-GM for positive reals: AB >= ((A+B)/2)^2
Step 3: ((A+B)/2)^2 >= (2*Phi_r/2)^2 = Phi_r^2 = ||h||^4

Therefore AB >= ||h||^4. QED (assuming A + B >= 2*Phi_r).
""")

np.random.seed(789)
print("Verifying the chain for 1000 random cases at each n:")

for n in [2, 3, 4, 5, 6]:
    chain_ok = 0
    total = 0
    amgm_tight = float('inf')  # how tight is AM-GM step
    sum_tight = float('inf')   # how tight is sum step

    for trial in range(1000):
        rp = random_roots(n, scale=2.0, min_gap=0.3)
        rq = random_roots(n, scale=2.0, min_gap=0.3)

        result = compute_quantities(rp, rq, use_formula=2)
        if result is None:
            continue

        total += 1
        A, B = result['A'], result['B']
        Phi_r = result['Phi_r']
        h4 = Phi_r**2

        # Chain check
        step1 = (A > 0 and B > 0)
        step2_val = A * B  # actual
        step2_bound = ((A + B) / 2)**2  # AM-GM bound
        step3_bound = Phi_r**2  # sum bound

        if step1 and step2_val >= step2_bound - 1e-10 and step2_bound >= step3_bound - 1e-10:
            chain_ok += 1

        if step2_bound > 0:
            amgm_tight = min(amgm_tight, step2_val / step2_bound)
        if step3_bound > 0 and A + B > 0:
            sum_tight = min(sum_tight, (A + B) / (2 * Phi_r))

    print(f"  n={n}: {chain_ok}/{total} pass, "
          f"AM-GM tightness min(AB/((A+B)/2)^2)={amgm_tight:.6f}, "
          f"Sum tightness min((A+B)/(2Phi_r))={sum_tight:.6f}")


# ================================================================
# PART 6: WHEN IS EQUALITY IN A + B = 2*Phi_r ACHIEVED?
# ================================================================

print("\n" + "=" * 70)
print("PART 6: EQUALITY CONDITIONS IN A + B = 2*Phi_r")
print("=" * 70)

print("\nFor n=2: equality iff gap_p = gap_q (proved analytically).")
print("\nFor n=3: is the infimum exactly 1, or strictly greater than 1?")

# Try equally-spaced roots with same spacing
for n in [3, 4, 5]:
    for s, t in [(1.0, 1.0), (0.5, 0.5), (2.0, 2.0)]:
        rp = np.arange(n, dtype=float) * s
        rq = np.arange(n, dtype=float) * t
        result = compute_quantities(rp, rq, use_formula=2)
        if result:
            print(f"  n={n}, equal gaps s=t={s}: (A+B)/(2Phi_r) = {result['ratio_ApB']:.10f}")

# For n=3, try p and q with IDENTICAL gap structure but different scales
print("\nn=3, same shape, different scale (gap_p/gap_q = c for all gaps):")
for c in [0.5, 0.8, 0.9, 0.95, 0.99, 1.0, 1.01, 1.05, 1.1, 1.5, 2.0]:
    rp = np.array([0.0, 1.0, 3.0])
    rq = rp * c
    result = compute_quantities(rp, rq, use_formula=2)
    if result:
        print(f"  c={c:.2f}: (A+B)/(2Phi_r) = {result['ratio_ApB']:.10f}")


# ================================================================
# PART 7: ATTEMPT ANALYTIC APPROACH FOR n >= 3
# ================================================================

print("\n" + "=" * 70)
print("PART 7: STRUCTURAL ANALYSIS")
print("=" * 70)

print("""
Key identity: A + B = Phi_p + Phi_q - 2*Phi_r.

The question is: Phi_p + Phi_q >= 4*Phi_r.

Rewrite: Phi_p + Phi_q >= 4*Phi_r
i.e., (Phi_p - Phi_r) + (Phi_q - Phi_r) >= 2*Phi_r
i.e., sum_k (u_k^2 - h_k^2) + sum_k (v_k^2 - h_k^2) >= 2 * sum_k h_k^2
i.e., sum_k (u_k^2 + v_k^2) >= 4 * sum_k h_k^2
i.e., ||u||^2 + ||v||^2 >= 4 * ||h||^2

Since h = (u + v)/2 ... NO, h is NOT the average of u and v in general.

Let's check: is h_k = (u_k + v_k)/2?
""")

np.random.seed(111)
for n in [2, 3, 4]:
    rp = random_roots(n, scale=2.0, min_gap=0.5)
    rq = random_roots(n, scale=2.0, min_gap=0.5)
    result = compute_quantities(rp, rq, use_formula=2)
    if result:
        h = H_values(result['roots_r'])
        u = H_values(result['roots_p'])
        v = H_values(result['roots_q'])
        avg = (u + v) / 2
        print(f"  n={n}: h = {h}")
        print(f"       (u+v)/2 = {avg}")
        print(f"       match? {np.allclose(h, avg, atol=1e-6)}")
        print()

print("""
h != (u+v)/2 in general. But the inequality ||u||^2 + ||v||^2 >= 4||h||^2
is equivalent to: ||u - v||^2 + ||u + v||^2 >= 4||h||^2 (parallelogram),
which is: ||u - v||^2 >= 4||h||^2 - ||u + v||^2.
This doesn't obviously simplify.

Alternative: ||u||^2 + ||v||^2 >= 4||h||^2 is just
  ||u||^2 + ||v||^2 - 2||h||^2 >= 2||h||^2
  A + B >= 2||h||^2 = 2*Phi_r

From u = h + alpha, v = h + beta:
  ||u||^2 + ||v||^2 = 2||h||^2 + 2<h,alpha+beta> + ||alpha||^2 + ||beta||^2

So A + B = 2<h,alpha+beta> + ||alpha||^2 + ||beta||^2.
Needing this >= 2*Phi_r = 2||h||^2.

So: ||alpha||^2 + ||beta||^2 + 2<h, alpha+beta> >= 2||h||^2.

Let w = alpha + beta = u + v - 2h.
Then: ||alpha||^2 + ||beta||^2 + 2<h, w> >= 2||h||^2.

Note ||alpha||^2 + ||beta||^2 >= ||w||^2/2 by Cauchy-Schwarz (|a|^2+|b|^2 >= |a+b|^2/2).
So: ||w||^2/2 + 2<h,w> >= 2||h||^2?
i.e., ||w||^2 + 4<h,w> >= 4||h||^2?
i.e., ||w||^2 + 4<h,w> + 4||h||^2 >= 8||h||^2?
i.e., ||w + 2h||^2 >= 8||h||^2?
i.e., ||u + v||^2 >= 8||h||^2?
This would require ||u+v|| >= 2*sqrt(2)*||h||, which is too strong in general.

So the triangle inequality approach is too lossy. The bound must use specific
properties of the MSS convolution structure.
""")


# ================================================================
# PART 8: COMPARISON OF FORMULAS - DO THEY GIVE SAME A+B?
# ================================================================

print("=" * 70)
print("PART 8: DO FORMULAS F1 AND F2 AGREE ON A+B/(2*Phi_r)?")
print("=" * 70)

np.random.seed(222)
max_diff = 0
for n in [2, 3, 4, 5]:
    for trial in range(100):
        rp = random_roots(n, scale=2.0, min_gap=0.4)
        rq = random_roots(n, scale=2.0, min_gap=0.4)

        r1 = compute_quantities(rp, rq, use_formula=1)
        r2 = compute_quantities(rp, rq, use_formula=2)

        if r1 is not None and r2 is not None:
            d = abs(r1['ratio_ApB'] - r2['ratio_ApB'])
            max_diff = max(max_diff, d)

print(f"  Max |F1_ratio - F2_ratio| over all tests: {max_diff:.2e}")
if max_diff < 0.01:
    print("  CONCLUSION: Both formulas agree to high precision on (A+B)/(2*Phi_r).")
else:
    print("  WARNING: Formulas disagree! Investigate further.")


# ================================================================
# PART 9: SPECIFIC STRUCTURED FAMILIES
# ================================================================

print("\n" + "=" * 70)
print("PART 9: TESTING SPECIFIC FAMILIES")
print("=" * 70)

# Family 1: p = {-s, 0, s}, q = {-t, 0, t} for n=3
print("\nFamily: p = {-s, 0, s}, q = {-t, 0, t} for various s, t:")
for s in [0.5, 1.0, 2.0, 5.0]:
    for t in [0.5, 1.0, 2.0, 5.0]:
        rp = np.array([-s, 0, s])
        rq = np.array([-t, 0, t])
        r = compute_quantities(rp, rq, use_formula=2)
        if r:
            print(f"  s={s:.1f}, t={t:.1f}: (A+B)/(2Phi_r) = {r['ratio_ApB']:.8f}, "
                  f"AB/h^4 = {r['AB_ratio']:.8f}")

# Family 2: p with one large gap, q equally spaced, for n=3
print("\nFamily: p = {0, 1, 1+L} (one large gap L), q = {0, 1, 2} for n=3:")
for L in [0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0]:
    rp = np.array([0.0, 1.0, 1.0 + L])
    rq = np.array([0.0, 1.0, 2.0])
    r = compute_quantities(rp, rq, use_formula=2)
    if r:
        print(f"  L={L:6.1f}: (A+B)/(2Phi_r) = {r['ratio_ApB']:.8f}, "
              f"AB/h^4 = {r['AB_ratio']:.8f}")

# Family 3: extreme gap ratio
print("\nFamily: p = {0, eps, 1} (tiny gap), q = {0, 1, 2} for n=3:")
for eps in [0.5, 0.1, 0.01, 0.001]:
    rp = np.array([0.0, eps, 1.0])
    rq = np.array([0.0, 1.0, 2.0])
    r = compute_quantities(rp, rq, use_formula=2)
    if r:
        print(f"  eps={eps:.4f}: (A+B)/(2Phi_r) = {r['ratio_ApB']:.8f}, "
              f"AB/h^4 = {r['AB_ratio']:.8f}")


# ================================================================
# PART 10: DECOMPOSITION OF A + B INTO CAUCHY-TYPE SUMS
# ================================================================

print("\n" + "=" * 70)
print("PART 10: CAUCHY KERNEL DECOMPOSITION")
print("=" * 70)

print("""
Phi_p = sum_k (sum_{j!=k} 1/(lambda_k - lambda_j))^2

For n=2: Phi = 2/gap^2 (one gap).
For n=3: Phi = (1/g1 + 1/g2)^2 + (-1/g1 + 1/(g1+g2))^2 + (-1/g2 - 1/(g1+g2))^2
         where g1 = lambda_2 - lambda_1, g2 = lambda_3 - lambda_2.

The n=3 case with gaps (g1, g2) for p and (h1, h2) for q:
  Phi_p = sum of squared Cauchy terms, involves 1/g1, 1/g2, 1/(g1+g2)
  Phi_q = sum of squared Cauchy terms, involves 1/h1, 1/h2, 1/(h1+h2)
  Phi_r = sum of squared Cauchy terms for r = p boxplus q.

Let's compute the n=3 boxplus gap structure more carefully.
""")

# For n=3, compute gap structure of r from p and q
print("n=3 gap analysis:")
for trial in range(5):
    rp = random_roots(3, scale=2.0, min_gap=0.5)
    rq = random_roots(3, scale=2.0, min_gap=0.5)
    result = compute_quantities(rp, rq, use_formula=2)
    if result is None:
        continue

    rr = result['roots_r']
    gp = np.diff(result['roots_p'])
    gq = np.diff(result['roots_q'])
    gr = np.diff(rr)

    # Is there a simple relation between gr and (gp, gq)?
    print(f"  gp = {gp}, gq = {gq}")
    print(f"  gr = {gr}")
    print(f"  gp^2+gq^2 = {gp**2 + gq**2}")
    print(f"  gr^2 = {gr**2}")
    print()


# ================================================================
# FINAL SUMMARY
# ================================================================

print("=" * 70)
print("FINAL SUMMARY")
print("=" * 70)
print("""
Results to be documented in findings_AplusB.md after running this script.
""")
