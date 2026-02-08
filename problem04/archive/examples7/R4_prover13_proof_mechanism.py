"""
PROVER-13 Part 3: The Proof Mechanism

ESTABLISHED FACTS:
1. dr_i/dt = H_i(r) exactly (root dynamics under heat flow p_t = p ⊞ G_t)
2. dS/dt = Phi_n(p) where S = sum_{i<j} log|r_i - r_j| (de Bruijn identity)
3. 1/Phi_n(p_t) is concave in t (numerically confirmed)
4. Phi_n(G_s) = n(n-1)/(4s)
5. G_s ⊞ G_t = G_{s+t}

ALGEBRAIC IDENTITY (established in Part 2):
  sum_{i<j} (H_i - H_j)/(r_i - r_j) = sum_i H_i^2 = Phi_n

This is remarkable: the FIRST derivative of S along heat flow equals Phi,
and the root dynamics are EXACTLY Dyson Brownian motion.

THIS SCRIPT: Investigate whether the Blachman-Stam argument works.
The classical proof of Stam's inequality is:

Method 1 (via EPI):
  1/I(X+Y) >= 1/I(X) + 1/I(Y)
  This follows from the entropy power inequality (EPI).

Method 2 (via heat flow + MMSE):
  More direct. Uses that 1/I is concave along heat flow.

Method 3 (Blachman-Stam via score functions):
  Score: rho_i(x) = d/dx log f_i(x)
  I(X+Y) = E[rho_{X+Y}(X+Y)^2]
  Key: rho_{X+Y}(x) = E[rho_X(X) | X+Y=x] (score decomposition)
  Then Cauchy-Schwarz: E[rho^2] >= (E[rho])^2 + Var(rho)

For finite free setting, the "score function" is H_p(lambda_i).
Let's investigate if there's an analog of the score decomposition.
"""

import numpy as np
from math import factorial, comb
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Core functions
# ============================================================

def finite_free_convolution(p_coeffs, q_coeffs):
    n = len(p_coeffs) - 1
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
    coeffs = np.polynomial.polynomial.polyfromroots(roots)[::-1]
    return coeffs / coeffs[0]

def roots_from_coeffs(coeffs):
    return np.sort(np.roots(coeffs)).real

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

def hermite_prob_coeffs(n):
    coeffs = np.zeros(n + 1)
    for m in range(n // 2 + 1):
        k = 2 * m
        coeff = factorial(n) * ((-1)**m) / (factorial(m) * (2**m) * factorial(n - 2*m))
        coeffs[k] = coeff
    return coeffs

def gaussian_poly_coeffs(n, s):
    he = hermite_prob_coeffs(n)
    gs = np.zeros(n + 1)
    for k in range(n + 1):
        gs[k] = he[k] * s**(k / 2.0)
    return gs


print("=" * 70)
print("PROVER-13 PART 3: PROOF MECHANISM")
print("=" * 70)

# ============================================================
# A. Compute d²S/dt² explicitly
# ============================================================

print("\n" + "=" * 70)
print("A. EXPLICIT d²S/dt² (Psi_n)")
print("=" * 70)
print()
print("Given dr_i/dt = H_i, we have:")
print("  S(t) = sum_{i<j} log|r_i(t) - r_j(t)|")
print("  dS/dt = sum_{i<j} (H_i - H_j)/(r_i - r_j) = Phi_n")
print()
print("For d²S/dt²:")
print("  d²S/dt² = sum_{i<j} d/dt [(H_i - H_j)/(r_i - r_j)]")
print("         = sum_{i<j} [(dH_i/dt - dH_j/dt)(r_i - r_j) - (H_i - H_j)(dr_i/dt - dr_j/dt)] / (r_i - r_j)^2")
print()
print("Since dr_i/dt = H_i:")
print("  d²S/dt² = sum_{i<j} [(dH_i/dt - dH_j/dt)/(r_i - r_j) - (H_i - H_j)^2/(r_i - r_j)^2]")
print()
print("And dH_i/dt = d/dt sum_{j!=i} 1/(r_i - r_j)")
print("            = sum_{j!=i} -(dr_i/dt - dr_j/dt)/(r_i - r_j)^2")
print("            = -sum_{j!=i} (H_i - H_j)/(r_i - r_j)^2")
print()

# Compute everything explicitly
for n in [3, 4, 5, 6]:
    for trial_name, roots_p in [("equi", np.arange(1, n+1, dtype=float)),
                                 ("quad", np.arange(1, n+1, dtype=float)**1.5)]:
        if np.min(np.diff(roots_p)) < 0.1:
            continue

        H = H_values(roots_p)
        phi = Phi_n(roots_p)

        # Compute dH_i/dt = -sum_{j!=i} (H_i - H_j)/(r_i - r_j)^2
        dH = np.zeros(n)
        for i in range(n):
            for j in range(n):
                if j != i:
                    dH[i] -= (H[i] - H[j]) / (roots_p[i] - roots_p[j])**2

        # Compute d²S/dt²
        d2S = 0.0
        for i in range(n):
            for j in range(i+1, n):
                term1 = (dH[i] - dH[j]) / (roots_p[i] - roots_p[j])
                term2 = (H[i] - H[j])**2 / (roots_p[i] - roots_p[j])**2
                d2S += term1 - term2

        # Also compute numerically for verification
        coeffs_p = coeffs_from_roots(roots_p)
        dt = 1e-5
        phis = []
        for k in range(5):
            if k == 0:
                phis.append(phi)
            else:
                gt = gaussian_poly_coeffs(n, k*dt)
                ct = finite_free_convolution(coeffs_p, gt)
                rt = np.sort(roots_from_coeffs(ct).real)
                phis.append(Phi_n(rt))

        d2S_num = (phis[2] - 2*phis[1] + phis[0]) / dt**2  # d²S/dt² = dPhi/dt

        Psi = -d2S  # should be positive
        Psi_num = -(phis[1] - phis[0]) / dt  # negative of dPhi/dt should be Psi
        # Actually d2S = dPhi/dt, so -d2S = -dPhi/dt = Psi

        if n <= 5:
            print(f"  n={n} ({trial_name}): d²S/dt² = {d2S:.6f} (analytic), dPhi/dt = {(phis[1]-phis[0])/dt:.6f} (num)")
            print(f"    Psi = -d²S = {Psi:.6f}, should be >= 0: {Psi >= -1e-6}")

            # Decompose Psi into its two parts
            part_A = 0.0  # sum of (dH_i - dH_j)/(r_i - r_j) terms
            part_B = 0.0  # sum of (H_i - H_j)^2/(r_i - r_j)^2 terms
            for i in range(n):
                for j in range(i+1, n):
                    part_A += (dH[i] - dH[j]) / (roots_p[i] - roots_p[j])
                    part_B += (H[i] - H[j])**2 / (roots_p[i] - roots_p[j])**2

            print(f"    Part A (sum dH terms): {part_A:.6f}")
            print(f"    Part B (sum (H-H)^2 terms): {part_B:.6f}")
            print(f"    d²S = A - B: {part_A - part_B:.6f}")
            print()


# ============================================================
# B. The crucial question: Can we derive d/dt(1/Phi) >= c?
# ============================================================

print("\n" + "=" * 70)
print("B. d/dt(1/Phi) = Psi/Phi^2 -- WHAT IS Psi?")
print("=" * 70)
print()
print("Psi = -d²S/dt² = part_B - part_A")
print("   = sum_{i<j} (H_i-H_j)^2/(r_i-r_j)^2 - sum_{i<j} (dH_i-dH_j)/(r_i-r_j)")
print()
print("If we can show Psi >= Phi^2 / (some function of cumulants),")
print("that gives d/dt(1/Phi) >= 1/(some function).")
print()

# Compute Psi / Phi^2 for many examples to find patterns
np.random.seed(42)

for n in [3, 4, 5, 6]:
    print(f"--- n = {n} ---")
    ratios = []
    for trial in range(100):
        roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n) * 2)
        if np.min(np.diff(roots_p)) < 0.2:
            continue

        H = H_values(roots_p)
        phi = Phi_n(roots_p)

        dH = np.zeros(n)
        for i in range(n):
            for j in range(n):
                if j != i:
                    dH[i] -= (H[i] - H[j]) / (roots_p[i] - roots_p[j])**2

        d2S = 0.0
        for i in range(n):
            for j in range(i+1, n):
                d2S += (dH[i] - dH[j]) / (roots_p[i] - roots_p[j])
                d2S -= (H[i] - H[j])**2 / (roots_p[i] - roots_p[j])**2

        Psi = -d2S
        ratio = Psi / phi**2 if phi > 0 else float('nan')
        ratios.append(ratio)

    if ratios:
        print(f"  Psi/Phi^2: min={min(ratios):.6f}, max={max(ratios):.6f}, mean={np.mean(ratios):.6f}")
        print(f"  This is d/dt(1/Phi)|_0.")
        # Check if min is bounded below by 4/(n(n-1))
        expected_asymptotic = 4.0 / (n * (n - 1))
        print(f"  4/(n(n-1)) = {expected_asymptotic:.6f}")
        print(f"  min ratio / (4/(n(n-1))) = {min(ratios)/expected_asymptotic:.6f}")
    print()


# ============================================================
# C. Direct approach: Cauchy-Schwarz on the root level
# ============================================================

print("\n" + "=" * 70)
print("C. CAUCHY-SCHWARZ APPROACH")
print("=" * 70)
print()
print("Phi_n(r) = sum_i H_i^2 where H_i = sum_{j!=i} 1/(r_i - r_j)")
print()
print("For p ⊞ q with roots mu_k, we want:")
print("  1/Phi(mu) >= 1/Phi(lambda) + 1/Phi(rho)")
print("where lambda = roots(p), rho = roots(q)")
print()
print("This is equivalent to:")
print("  Phi(lambda) * Phi(rho) >= Phi(mu) * [Phi(lambda) + Phi(rho)]")
print()
print("Or: Phi(mu) <= Phi(lambda)*Phi(rho)/[Phi(lambda)+Phi(rho)] = harmonic_mean/2")
print()
print("So the conjecture says: Phi(p ⊞ q) <= HARMONIC MEAN of Phi(p), Phi(q)")
print()

# Verify this formulation
np.random.seed(42)
violations = 0
total = 0

for n in [3, 4, 5, 6]:
    n_viol = 0
    n_tot = 0
    for trial in range(200):
        roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n))
        roots_q = np.sort(np.random.randn(n) * 1.5 + np.arange(n) * 0.5)
        if np.min(np.diff(roots_p)) < 0.1 or np.min(np.diff(roots_q)) < 0.1:
            continue

        coeffs_p = coeffs_from_roots(roots_p)
        coeffs_q = coeffs_from_roots(roots_q)
        coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)
        roots_pq = roots_from_coeffs(coeffs_pq)
        if np.any(np.abs(np.imag(roots_pq)) > 1e-8):
            continue

        phi_p = Phi_n(roots_p)
        phi_q = Phi_n(roots_q)
        phi_pq = Phi_n(np.sort(roots_pq.real))

        harmonic_mean = 2 * phi_p * phi_q / (phi_p + phi_q) if (phi_p + phi_q) > 0 else 0

        n_tot += 1
        total += 1
        if phi_pq > harmonic_mean + 1e-10:
            n_viol += 1
            violations += 1

    print(f"  n={n}: {n_viol} violations / {n_tot} (harmonic mean bound)")

print(f"\nTotal: {violations} / {total}")
if violations == 0:
    print("*** Phi(p⊞q) <= HarmonicMean(Phi(p), Phi(q)) CONFIRMED ***")


# ============================================================
# D. The EPI approach more carefully
# ============================================================

print("\n" + "=" * 70)
print("D. EPI APPROACH: S(p ⊞ G_t) BEHAVIOR")
print("=" * 70)
print()
print("Define N(p) = exp(2*S(p)/m) where m = binom(n,2) = n(n-1)/2")
print("(the number of log-gap terms)")
print()
print("Analog of entropy power: N(p) = exp(2*S/m)")
print("EPI analog: N(p ⊞ q) >= N(p) + N(q)?")
print()

def N_entropy_power(roots, n):
    """Finite free entropy power."""
    m = n * (n - 1) // 2  # number of pairs
    S = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            S += np.log(abs(roots[j] - roots[i]))
    return np.exp(2 * S / m)


# Test EPI
np.random.seed(42)
for n_val in [3, 4, 5, 6]:
    n = n_val
    epi_violations = 0
    epi_total = 0
    for trial in range(200):
        roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n))
        roots_q = np.sort(np.random.randn(n) * 1.5 + np.arange(n) * 0.5)
        if np.min(np.diff(roots_p)) < 0.1 or np.min(np.diff(roots_q)) < 0.1:
            continue

        coeffs_p = coeffs_from_roots(roots_p)
        coeffs_q = coeffs_from_roots(roots_q)
        coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)
        roots_pq = roots_from_coeffs(coeffs_pq)
        if np.any(np.abs(np.imag(roots_pq)) > 1e-8):
            continue

        N_p = N_entropy_power(roots_p, n)
        N_q = N_entropy_power(roots_q, n)
        N_pq = N_entropy_power(np.sort(roots_pq.real), n)

        epi_total += 1
        if N_pq < N_p + N_q - 1e-10:
            epi_violations += 1

    print(f"  n={n}: N(p⊞q) >= N(p) + N(q)? {epi_violations} violations / {epi_total}")


# Also try with different normalizations
print("\nTrying N = exp(S/m) (without the factor 2):")
for n_val in [3, 4, 5, 6]:
    n = n_val
    m = n * (n - 1) // 2
    epi_violations = 0
    epi_total = 0
    for trial in range(200):
        roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n))
        roots_q = np.sort(np.random.randn(n) * 1.5 + np.arange(n) * 0.5)
        if np.min(np.diff(roots_p)) < 0.1 or np.min(np.diff(roots_q)) < 0.1:
            continue

        coeffs_p = coeffs_from_roots(roots_p)
        coeffs_q = coeffs_from_roots(roots_q)
        coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)
        roots_pq = roots_from_coeffs(coeffs_pq)
        if np.any(np.abs(np.imag(roots_pq)) > 1e-8):
            continue

        S_p = sum([np.log(abs(roots_p[j] - roots_p[i]))
                   for i in range(n) for j in range(i+1, n)])
        S_q = sum([np.log(abs(roots_q[j] - roots_q[i]))
                   for i in range(n) for j in range(i+1, n)])
        roots_pq_sorted = np.sort(roots_pq.real)
        S_pq = sum([np.log(abs(roots_pq_sorted[j] - roots_pq_sorted[i]))
                    for i in range(n) for j in range(i+1, n)])

        N_p = np.exp(S_p / m)
        N_q = np.exp(S_q / m)
        N_pq = np.exp(S_pq / m)

        epi_total += 1
        if N_pq < N_p + N_q - 1e-10:
            epi_violations += 1

    print(f"  n={n}: N(p⊞q) >= N(p) + N(q)? {epi_violations} violations / {epi_total}")


# ============================================================
# E. The Blachman-Stam argument with finite free Fisher
# ============================================================

print("\n" + "=" * 70)
print("E. BLACHMAN-STAM ARGUMENT")
print("=" * 70)
print()
print("Classical Blachman-Stam:")
print("  X = X + 0, Y = 0 + Y")
print("  X + Y has score rho_{X+Y}(x) = E[rho_X(X) | X+Y=x]")
print("  By Cauchy-Schwarz on conditional expectations:")
print("  E[rho_{X+Y}^2] * E[rho_X^2] >= (E[rho_X^2])^2")
print("  ... leading to 1/I(X+Y) >= 1/I(X) + 1/I(Y)")
print()
print("Finite free analog:")
print("  For p ⊞ q, the roots mu_k are deterministic.")
print("  But via the Hermitian random matrix model:")
print("  mu_k = eigenvalues of A + UBU*, averaged over Haar(U)")
print("  So the Blachman argument needs a matrix-level formulation.")
print()
print("Alternative: use the heat flow directly.")
print()

# ============================================================
# F. Heat flow proof attempt
# ============================================================

print("\n" + "=" * 70)
print("F. HEAT FLOW PROOF ATTEMPT")
print("=" * 70)
print()
print("Define f_p(t) = 1/Phi_n(p ⊞ G_t).")
print("We have established:")
print("  (i)   f_p(t) is concave in t")
print("  (ii)  f_p(0) = 1/Phi_n(p)")
print("  (iii) f_p(t) ~ 4t/(n(n-1)) as t -> inf")
print()
print("Key identity: p ⊞ q ⊞ G_t = (p ⊞ G_{t/2}) ⊞ (q ⊞ G_{t/2})")
print("Wait -- is this true? Let's check!")
print()

# Check: p ⊞ q ⊞ G_t = (p ⊞ G_s) ⊞ (q ⊞ G_{t-s}) for any s?
# NO! This would require G_t to "split" between the two terms.
# The correct identity is just associativity: (p ⊞ q) ⊞ G_t = p ⊞ (q ⊞ G_t)

np.random.seed(42)
for n in [3, 4]:
    roots_p = np.arange(1, n+1, dtype=float)
    roots_q = np.arange(1, n+1, dtype=float) * 0.5
    coeffs_p = coeffs_from_roots(roots_p)
    coeffs_q = coeffs_from_roots(roots_q)

    t = 1.0
    s = t / 2

    # LHS: (p ⊞ q) ⊞ G_t
    pq = finite_free_convolution(coeffs_p, coeffs_q)
    gt = gaussian_poly_coeffs(n, t)
    lhs = finite_free_convolution(pq, gt)

    # RHS attempt: (p ⊞ G_s) ⊞ (q ⊞ G_s)
    gs = gaussian_poly_coeffs(n, s)
    p_gs = finite_free_convolution(coeffs_p, gs)
    q_gs = finite_free_convolution(coeffs_q, gs)
    rhs = finite_free_convolution(p_gs, q_gs)

    err = np.max(np.abs(lhs - rhs))
    print(f"  n={n}: |(p⊞q)⊞G_t - (p⊞G_{s})⊞(q⊞G_{s})| = {err:.6e}")
    # Compare with (p ⊞ G_t) ⊞ q = p ⊞ (q ⊞ G_t)
    p_gt = finite_free_convolution(coeffs_p, gt)
    rhs2 = finite_free_convolution(p_gt, coeffs_q)
    err2 = np.max(np.abs(lhs - rhs2))
    print(f"  n={n}: |(p⊞q)⊞G_t - (p⊞G_t)⊞q| = {err2:.6e}  (associativity)")

print()
print("CRITICAL: (p⊞q)⊞G_t != (p⊞G_{t/2})⊞(q⊞G_{t/2})")
print("This is because ⊞ is NOT the same as sum of independent RVs!")
print("It's an average over Haar measure, not a product of independent variables.")
print()
print("BUT: in the random matrix model,")
print("  A + UBU* + V(tI)V* = A + t*I + UBU*  (if adding to diagonal)")
print("So adding Gaussian noise commutes with addition: this IS associative.")
print("The identity (p⊞q)⊞G_t = (p⊞G_t)⊞q holds by associativity.")
print()

# ============================================================
# G. The CORRECT approach: interpolation via heat flow
# ============================================================

print("=" * 70)
print("G. INTERPOLATION APPROACH")
print("=" * 70)
print()
print("Define for 0 <= alpha <= 1:")
print("  p_alpha = p ⊞ G_{alpha*s}")
print("  q_alpha = q ⊞ G_{alpha*t}")
print()
print("Then r_alpha = p_alpha ⊞ q_alpha has:")
print("  Phi(r_0) = Phi(p ⊞ q)")
print("  As alpha -> inf, Phi(r_alpha) -> Phi(G_{alpha*s} ⊞ G_{alpha*t})")
print("                               = Phi(G_{alpha(s+t)})")
print("                               = n(n-1)/(4*alpha*(s+t))")
print()
print("Key insight: at alpha=0, r_0 = p ⊞ q (what we want)")
print("At alpha -> inf: 1/Phi(r_alpha) ~ 4*alpha*(s+t)/(n(n-1))")
print("  and 1/Phi(p_alpha) ~ 4*alpha*s/(n(n-1))")
print("  and 1/Phi(q_alpha) ~ 4*alpha*t/(n(n-1))")
print("  So 1/Phi(p_alpha) + 1/Phi(q_alpha) ~ 4*alpha*(s+t)/(n(n-1)) = 1/Phi(r_alpha)")
print()
print("So the Stam gap 1/Phi(r) - 1/Phi(p) - 1/Phi(q) starts at some value at alpha=0")
print("and approaches 0 as alpha -> inf.")
print()
print("If we can show the gap is DECREASING in alpha (the gap closes monotonically),")
print("then the gap at alpha=0 must be >= 0!")
print()

# Test this!
np.random.seed(42)
for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    for trial in range(5):
        roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n))
        roots_q = np.sort(np.random.randn(n) * 1.5 + np.arange(n) * 0.5)
        if np.min(np.diff(roots_p)) < 0.15 or np.min(np.diff(roots_q)) < 0.15:
            continue

        coeffs_p = coeffs_from_roots(roots_p)
        coeffs_q = coeffs_from_roots(roots_q)

        # Choose s, t as the "variance" parameters
        # s = 1/Phi(p), t = 1/Phi(q)
        phi_p = Phi_n(roots_p)
        phi_q = Phi_n(roots_q)
        s_param = 1.0 / phi_p
        t_param = 1.0 / phi_q

        alphas = np.linspace(0.001, 5.0, 100)
        gap_vals = []

        for alpha in alphas:
            gs_a = gaussian_poly_coeffs(n, alpha * s_param)
            gt_a = gaussian_poly_coeffs(n, alpha * t_param)

            p_a = finite_free_convolution(coeffs_p, gs_a)
            q_a = finite_free_convolution(coeffs_q, gt_a)
            r_a = finite_free_convolution(p_a, q_a)

            rts_p = roots_from_coeffs(p_a)
            rts_q = roots_from_coeffs(q_a)
            rts_r = roots_from_coeffs(r_a)

            if (np.any(np.abs(np.imag(rts_p)) > 1e-8) or
                np.any(np.abs(np.imag(rts_q)) > 1e-8) or
                np.any(np.abs(np.imag(rts_r)) > 1e-8)):
                gap_vals.append(float('nan'))
                continue

            inv_phi_p = 1.0 / Phi_n(np.sort(rts_p.real))
            inv_phi_q = 1.0 / Phi_n(np.sort(rts_q.real))
            inv_phi_r = 1.0 / Phi_n(np.sort(rts_r.real))

            gap = inv_phi_r - inv_phi_p - inv_phi_q
            gap_vals.append(gap)

        gap_vals = np.array(gap_vals)
        valid = ~np.isnan(gap_vals)
        if np.sum(valid) < 10:
            continue

        g = gap_vals[valid]
        a = alphas[valid]

        print(f"  Trial {trial}:")
        print(f"    gap at alpha~0: {g[0]:.6f}")
        print(f"    gap at alpha=5: {g[-1]:.6f}")
        print(f"    gap min: {np.min(g):.6f}")
        gap_diffs = np.diff(g)
        n_dec = np.sum(gap_diffs < -1e-10)
        n_inc = np.sum(gap_diffs > 1e-10)
        print(f"    gap changes: {n_dec} decrease, {n_inc} increase")
        print(f"    gap monotonically decreasing: {n_inc == 0}")


# ============================================================
# H. Alternative: FIXED noise, varying mixture
# ============================================================

print("\n" + "=" * 70)
print("H. FIXED NOISE APPROACH")
print("=" * 70)
print()
print("For FIXED t > 0, define:")
print("  gap(t) = 1/Phi(p⊞q⊞G_t) - 1/Phi(p⊞G_t) - 1/Phi(q⊞G_t)")
print()
print("We showed gap(t) -> 0 as t -> inf.")
print("If gap(t) is monotonically DECREASING, then gap(0) >= gap(t) -> 0.")
print()

np.random.seed(42)
for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    for trial in range(3):
        roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n))
        roots_q = np.sort(np.random.randn(n) * 1.5 + np.arange(n) * 0.5)
        if np.min(np.diff(roots_p)) < 0.15 or np.min(np.diff(roots_q)) < 0.15:
            continue

        coeffs_p = coeffs_from_roots(roots_p)
        coeffs_q = coeffs_from_roots(roots_q)
        coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)

        ts = np.linspace(0.001, 5.0, 200)
        gaps = []

        for t in ts:
            gt = gaussian_poly_coeffs(n, t)

            c_pq_t = finite_free_convolution(coeffs_pq, gt)
            c_p_t = finite_free_convolution(coeffs_p, gt)
            c_q_t = finite_free_convolution(coeffs_q, gt)

            r_pq = roots_from_coeffs(c_pq_t)
            r_p = roots_from_coeffs(c_p_t)
            r_q = roots_from_coeffs(c_q_t)

            if (np.any(np.abs(np.imag(r_pq)) > 1e-8) or
                np.any(np.abs(np.imag(r_p)) > 1e-8) or
                np.any(np.abs(np.imag(r_q)) > 1e-8)):
                gaps.append(float('nan'))
                continue

            inv_pq = 1.0 / Phi_n(np.sort(r_pq.real))
            inv_p = 1.0 / Phi_n(np.sort(r_p.real))
            inv_q = 1.0 / Phi_n(np.sort(r_q.real))

            gaps.append(inv_pq - inv_p - inv_q)

        gaps = np.array(gaps)
        valid = ~np.isnan(gaps)
        if np.sum(valid) < 10:
            continue

        g = gaps[valid]
        print(f"  Trial {trial}:")
        print(f"    gap at t~0: {g[0]:.6f}")
        print(f"    gap at t=5: {g[-1]:.6f}")
        print(f"    gap min: {np.min(g):.6f}")
        print(f"    gap monotonically decreasing: {np.all(np.diff(g) <= 1e-8)}")


# ============================================================
# I. What IS the derivative d/dt gap(t)?
# ============================================================

print("\n" + "=" * 70)
print("I. DERIVATIVE OF THE GAP: d/dt[1/Phi(r_t) - 1/Phi(p_t) - 1/Phi(q_t)]")
print("=" * 70)
print()
print("d/dt(1/Phi(p_t)) = Psi_n(p_t) / Phi_n(p_t)^2")
print("So d/dt(gap) = Psi(r_t)/Phi(r_t)^2 - Psi(p_t)/Phi(p_t)^2 - Psi(q_t)/Phi(q_t)^2")
print()
print("For the gap to be decreasing, we need:")
print("  Psi(r_t)/Phi(r_t)^2 <= Psi(p_t)/Phi(p_t)^2 + Psi(q_t)/Phi(q_t)^2")
print()
print("Since Psi/Phi^2 = d/dt(1/Phi), this is asking for SUBADDITIVITY of d/dt(1/Phi).")
print("That would follow if d/dt(1/Phi) is a CONCAVE function of the polynomial.")
print()

# Compute d/dt(1/Phi) for p⊞q, p, q
np.random.seed(42)
for n in [3, 4, 5]:
    print(f"\n--- n = {n} ---")
    for trial in range(5):
        roots_p = np.sort(np.random.randn(n) * 2 + np.arange(n))
        roots_q = np.sort(np.random.randn(n) * 1.5 + np.arange(n) * 0.5)
        if np.min(np.diff(roots_p)) < 0.15 or np.min(np.diff(roots_q)) < 0.15:
            continue

        coeffs_p = coeffs_from_roots(roots_p)
        coeffs_q = coeffs_from_roots(roots_q)
        coeffs_pq = finite_free_convolution(coeffs_p, coeffs_q)
        roots_pq = roots_from_coeffs(coeffs_pq)
        if np.any(np.abs(np.imag(roots_pq)) > 1e-8):
            continue
        roots_pq = np.sort(roots_pq.real)

        # Compute d/dt(1/Phi) at t=0 for each
        def compute_d_inv_phi(coeffs_base, roots_base):
            phi = Phi_n(roots_base)
            dt = 1e-5
            gt = gaussian_poly_coeffs(n, dt)
            ct = finite_free_convolution(coeffs_base, gt)
            rt = np.sort(roots_from_coeffs(ct).real)
            phi_dt = Phi_n(rt)
            return (1/phi_dt - 1/phi) / dt

        d_inv_p = compute_d_inv_phi(coeffs_p, roots_p)
        d_inv_q = compute_d_inv_phi(coeffs_q, roots_q)
        d_inv_pq = compute_d_inv_phi(coeffs_pq, roots_pq)

        print(f"  Trial {trial}:")
        print(f"    d/dt(1/Phi) at p: {d_inv_p:.6f}")
        print(f"    d/dt(1/Phi) at q: {d_inv_q:.6f}")
        print(f"    d/dt(1/Phi) at p⊞q: {d_inv_pq:.6f}")
        print(f"    sum at p,q: {d_inv_p + d_inv_q:.6f}")
        print(f"    d(p⊞q) <= d(p)+d(q)? {d_inv_pq <= d_inv_p + d_inv_q + 1e-8}")


print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)
