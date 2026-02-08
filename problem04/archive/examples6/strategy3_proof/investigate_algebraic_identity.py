"""
Investigate algebraic identity for AB - ||h||^4 in the MSS finite free convolution setting.

Setup:
  p, q monic real-rooted degree n, simple roots.
  r = p boxplus_n q (MSS convolution).
  h_k = H_r(nu_k) = sum_{j!=k} 1/(nu_k - nu_j)
  u_k = H_p(lambda_k), v_k = H_q(mu_k)
  alpha_k = u_k - h_k, beta_k = v_k - h_k
  A = ||u||^2 - ||h||^2 = 2<h,alpha> + ||alpha||^2
  B = ||v||^2 - ||h||^2 = 2<h,beta> + ||beta||^2
  TARGET: AB >= ||h||^4

We investigate:
  1. Numerical computation of AB - ||h||^4 for n=3
  2. Sum-of-squares decomposition search
  3. Gram matrix G = [<e_i, e_j>] for e = {h, alpha, beta}
  4. Four-vectors approach: h, alpha, beta, gamma = alpha - beta
  5. Ratio (AB - ||h||^4) / ||h||^4 as function of angles/magnitudes
"""

import numpy as np
from math import factorial, comb
from itertools import combinations
from scipy.optimize import minimize
import warnings
warnings.filterwarnings('ignore')

# ================================================================
# CORE FUNCTIONS
# ================================================================

def elem_sym_poly(roots, k):
    """Elementary symmetric polynomial e_k(roots)."""
    n = len(roots)
    if k == 0: return 1.0
    if k > n: return 0.0
    return sum(np.prod([roots[i] for i in subset])
               for subset in combinations(range(n), k))

def poly_coeffs_from_roots(roots):
    """Monic polynomial coefficients from roots: x^n - e1 x^{n-1} + ..."""
    n = len(roots)
    ek = [elem_sym_poly(roots, k) for k in range(n+1)]
    return [(-1)**k * ek[k] for k in range(n+1)]

def boxplus_mss(roots_p, roots_q):
    """MSS finite free additive convolution.

    c_k = sum_{i+j=k} [(n-i)!(n-j)! / (n!(n-k)!)] * a_i * b_j
    where a_i, b_j are the polynomial coefficients.
    """
    n = len(roots_p)
    a = poly_coeffs_from_roots(roots_p)
    b = poly_coeffs_from_roots(roots_q)
    c = np.zeros(n+1)
    for k in range(n+1):
        for i in range(k+1):
            j = k - i
            if i <= n and j <= n:
                coeff = factorial(n-i) * factorial(n-j) / (factorial(n) * factorial(n-k))
                c[k] += coeff * a[i] * b[j]
    roots_r = np.sort(np.real(np.roots(c)))
    return roots_r, c

def H_values(roots):
    """H_f(root_k) = sum_{j!=k} 1/(root_k - root_j) for each k."""
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H

def Phi_n(roots):
    """Phi = sum_k H_k^2 = ||H||^2."""
    return np.sum(H_values(roots)**2)

def compute_all_vectors(roots_p, roots_q):
    """Compute r, h, u, v, alpha, beta and all derived quantities.

    Returns dict with all quantities, or None if computation fails.
    """
    n = len(roots_p)

    # Compute r = p boxplus_n q
    roots_r, c_coeffs = boxplus_mss(roots_p, roots_q)

    # Check real-rootedness
    raw_roots = np.roots(c_coeffs)
    if np.any(np.abs(np.imag(raw_roots)) > 1e-6):
        return None
    roots_r = np.sort(np.real(raw_roots))

    # Check distinct roots
    if np.any(np.diff(roots_r) < 1e-6):
        return None

    # Compute H-values
    h = H_values(roots_r)
    u = H_values(roots_p)
    v = H_values(roots_q)

    # alpha and beta
    alpha = u - h
    beta = v - h

    # Key quantities
    h2 = np.dot(h, h)       # ||h||^2 = Phi_r
    u2 = np.dot(u, u)       # ||u||^2 = Phi_p
    v2 = np.dot(v, v)       # ||v||^2 = Phi_q
    a2 = np.dot(alpha, alpha)
    b2 = np.dot(beta, beta)
    ha = np.dot(h, alpha)
    hb = np.dot(h, beta)
    ab_ip = np.dot(alpha, beta)  # <alpha, beta>

    A = u2 - h2   # = 2*ha + a2
    B = v2 - h2   # = 2*hb + b2
    h4 = h2**2

    if A <= 0 or B <= 0 or h4 <= 0:
        return None

    return {
        'n': n,
        'roots_p': roots_p, 'roots_q': roots_q, 'roots_r': roots_r,
        'h': h, 'u': u, 'v': v, 'alpha': alpha, 'beta': beta,
        'h2': h2, 'u2': u2, 'v2': v2, 'a2': a2, 'b2': b2,
        'ha': ha, 'hb': hb, 'ab_ip': ab_ip,
        'A': A, 'B': B, 'h4': h4,
        'AB': A * B,
        'excess': A * B - h4,  # AB - ||h||^4
        'ratio': A * B / h4,   # AB / ||h||^4
    }

def random_well_separated_roots(n, scale=2.0, min_gap=0.5):
    """Generate random sorted roots with minimum gap."""
    roots = np.sort(np.random.randn(n) * scale)
    for i in range(1, n):
        if roots[i] - roots[i-1] < min_gap:
            roots[i] = roots[i-1] + min_gap
    return roots


# ================================================================
# PART 1: NUMERICAL COMPUTATION OF AB - ||h||^4 FOR n=3
# ================================================================

print("=" * 70)
print("PART 1: AB - ||h||^4 FOR n=3 (NUMERICAL SURVEY)")
print("=" * 70)

np.random.seed(42)

n3_data = []
for trial in range(2000):
    roots_p = random_well_separated_roots(3)
    roots_q = random_well_separated_roots(3)

    result = compute_all_vectors(roots_p, roots_q)
    if result is not None:
        n3_data.append(result)

print(f"\nValid trials: {len(n3_data)}")
excesses = [d['excess'] for d in n3_data]
ratios = [d['ratio'] for d in n3_data]
print(f"AB - ||h||^4: min={min(excesses):.10f}, max={max(excesses):.6f}, mean={np.mean(excesses):.6f}")
print(f"AB / ||h||^4: min={min(ratios):.10f}, max={max(ratios):.6f}, mean={np.mean(ratios):.6f}")
print(f"All AB >= ||h||^4: {all(e >= -1e-10 for e in excesses)}")

# Also check n=2 (should be exact equality)
print("\nn=2 check (should be exact equality):")
n2_excesses = []
for trial in range(500):
    roots_p = random_well_separated_roots(2)
    roots_q = random_well_separated_roots(2)
    result = compute_all_vectors(roots_p, roots_q)
    if result is not None:
        n2_excesses.append(result['excess'])
if n2_excesses:
    print(f"AB - ||h||^4: min={min(n2_excesses):.2e}, max={max(n2_excesses):.2e}")
    print(f"All near zero: {all(abs(e) < 1e-8 for e in n2_excesses)}")

# Check n=4, n=5
for nn in [4, 5]:
    print(f"\nn={nn} check:")
    excesses_nn = []
    for trial in range(500):
        roots_p = random_well_separated_roots(nn)
        roots_q = random_well_separated_roots(nn)
        result = compute_all_vectors(roots_p, roots_q)
        if result is not None:
            excesses_nn.append(result['excess'])
    if excesses_nn:
        ratios_nn = [d / (d + 1) if d != 0 else 0 for d in excesses_nn]  # just excess
        print(f"  Trials: {len(excesses_nn)}")
        print(f"  AB - ||h||^4: min={min(excesses_nn):.10f}, max={max(excesses_nn):.6f}")
        print(f"  All AB >= ||h||^4: {all(e >= -1e-8 for e in excesses_nn)}")


# ================================================================
# PART 2: SUM-OF-SQUARES SEARCH / PATTERN ANALYSIS FOR n=3
# ================================================================

print("\n\n" + "=" * 70)
print("PART 2: STRUCTURE OF AB - ||h||^4 IN TERMS OF alpha, beta, h")
print("=" * 70)

# Expand AB - ||h||^4 in terms of Gram matrix entries
# A = 2<h,alpha> + ||alpha||^2, B = 2<h,beta> + ||beta||^2
# AB = (2ha + a2)(2hb + b2)
#    = 4*ha*hb + 2*ha*b2 + 2*a2*hb + a2*b2
# ||h||^4 = h2^2
# So: AB - h4 = 4*ha*hb + 2*ha*b2 + 2*a2*hb + a2*b2 - h2^2

print("\nVerifying the expansion AB - h4 = 4*ha*hb + 2*ha*b2 + 2*a2*hb + a2*b2 - h2^2:")
for d in n3_data[:5]:
    ha, hb, a2, b2, h2, h4 = d['ha'], d['hb'], d['a2'], d['b2'], d['h2'], d['h4']
    computed = 4*ha*hb + 2*ha*b2 + 2*a2*hb + a2*b2 - h2**2
    print(f"  AB-h4 = {d['excess']:.10f}, formula = {computed:.10f}, match = {abs(d['excess'] - computed) < 1e-8}")

# Alternative expansion using the CROSS TERM <alpha, beta>
# Note: 4*ha*hb = 4 * <h,alpha> * <h,beta>
# Can we relate <h,alpha><h,beta> to something involving <alpha,beta>?
# By Cauchy-Schwarz in the h direction: <h,alpha><h,beta> / ||h||^2 is related to projections.

# Let's decompose alpha = alpha_h + alpha_perp where alpha_h = (<h,alpha>/||h||^2)*h
# Similarly beta = beta_h + beta_perp
# Then <h,alpha> = alpha_h * ||h||, etc.

print("\nDecomposition alpha = proj_h(alpha) + alpha_perp, etc.:")
print("Let alpha_h = <h,alpha>/||h||^2, beta_h = <h,beta>/||h||^2")
print("Then A = 2*alpha_h*h2 + a2 = h2*(2*alpha_h + a2/h2)")
print("     B = 2*beta_h*h2 + b2 = h2*(2*beta_h + b2/h2)")
print("     AB/h4 = (2*alpha_h + a2/h2)*(2*beta_h + b2/h2)")
print("We need this >= 1.\n")

# Compute normalized quantities
print("Statistics of normalized quantities for n=3:")
alpha_h_vals = [d['ha']/d['h2'] for d in n3_data]
beta_h_vals = [d['hb']/d['h2'] for d in n3_data]
a2_norm = [d['a2']/d['h2'] for d in n3_data]
b2_norm = [d['b2']/d['h2'] for d in n3_data]

print(f"  alpha_h = <h,alpha>/||h||^2: min={min(alpha_h_vals):.6f}, max={max(alpha_h_vals):.6f}, mean={np.mean(alpha_h_vals):.6f}")
print(f"  beta_h  = <h,beta>/||h||^2:  min={min(beta_h_vals):.6f}, max={max(beta_h_vals):.6f}, mean={np.mean(beta_h_vals):.6f}")
print(f"  ||alpha||^2/||h||^2:         min={min(a2_norm):.6f}, max={max(a2_norm):.6f}, mean={np.mean(a2_norm):.6f}")
print(f"  ||beta||^2/||h||^2:          min={min(b2_norm):.6f}, max={max(b2_norm):.6f}, mean={np.mean(b2_norm):.6f}")

# CRITICAL: is <h,alpha> always >= 0?
ha_negative = [d for d in n3_data if d['ha'] < -1e-10]
hb_negative = [d for d in n3_data if d['hb'] < -1e-10]
print(f"\n  <h,alpha> < 0 count: {len(ha_negative)}/{len(n3_data)}")
print(f"  <h,beta>  < 0 count: {len(hb_negative)}/{len(n3_data)}")

if ha_negative:
    worst = min(ha_negative, key=lambda d: d['ha'])
    print(f"  Worst <h,alpha> = {worst['ha']:.10f}")
    print(f"    roots_p = {worst['roots_p']}")
    print(f"    roots_q = {worst['roots_q']}")
    print(f"    But AB - h4 still = {worst['excess']:.10f} (positive: {worst['excess'] > -1e-10})")


# ================================================================
# PART 3: GRAM MATRIX ANALYSIS
# ================================================================

print("\n\n" + "=" * 70)
print("PART 3: GRAM MATRIX G = [<e_i, e_j>] for {h, alpha, beta}")
print("=" * 70)

# The 3x3 Gram matrix:
# G = [[h2,    ha,    hb   ],
#      [ha,    a2,    ab_ip],
#      [hb,    ab_ip, b2   ]]
#
# A = 2*ha + a2, B = 2*hb + b2
# AB = (2*ha + a2)(2*hb + b2) = 4*ha*hb + 2*ha*b2 + 2*a2*hb + a2*b2
# h4 = h2^2
#
# So AB - h4 = 4*ha*hb + 2*ha*b2 + 2*a2*hb + a2*b2 - h2^2
#
# In terms of Gram matrix entries G[i,j]:
# AB - h4 = 4*G[0,1]*G[0,2] + 2*G[0,1]*G[2,2] + 2*G[1,1]*G[0,2] + G[1,1]*G[2,2] - G[0,0]^2

print("\nGram matrix properties for n=3:")
det_G_vals = []
min_eig_vals = []

for d in n3_data:
    G = np.array([
        [d['h2'],   d['ha'],   d['hb']],
        [d['ha'],   d['a2'],   d['ab_ip']],
        [d['hb'],   d['ab_ip'], d['b2']],
    ])
    det_G_vals.append(np.linalg.det(G))
    eigs = np.sort(np.linalg.eigvalsh(G))
    min_eig_vals.append(eigs[0])

print(f"  det(G): min={min(det_G_vals):.6e}, max={max(det_G_vals):.6e}")
print(f"  min eigenvalue of G: min={min(min_eig_vals):.6e}, max={max(min_eig_vals):.6e}")
print(f"  G always PSD: {all(e >= -1e-10 for e in min_eig_vals)}")

# Schur complement approach:
# G is PSD (since it's a Gram matrix of real vectors).
# The Schur complement of G[0,0] = h2 in G is:
# S = [[a2, ab_ip], [ab_ip, b2]] - (1/h2)*[[ha], [hb]]*[[ha, hb]]
#   = [[a2 - ha^2/h2, ab_ip - ha*hb/h2],
#      [ab_ip - ha*hb/h2, b2 - hb^2/h2]]
# This is PSD. So:
# det(S) >= 0
# => (a2 - ha^2/h2)(b2 - hb^2/h2) >= (ab_ip - ha*hb/h2)^2

# Can we express AB - h4 using the Schur complement?
# Define: a2_perp = a2 - ha^2/h2 = ||alpha_perp||^2
#         b2_perp = b2 - hb^2/h2 = ||beta_perp||^2
#         ab_perp = ab_ip - ha*hb/h2 = <alpha_perp, beta_perp>
# Schur PSD => a2_perp * b2_perp >= ab_perp^2

# Now: AB - h4 = 4*ha*hb + 2*ha*b2 + 2*a2*hb + a2*b2 - h2^2
# Let's use substitution: ha = x*h2 (i.e., x = alpha_h), hb = y*h2
# Then A = h2*(2x + a2/h2), B = h2*(2y + b2/h2)
# AB = h2^2 * (2x + a2/h2)(2y + b2/h2)
# AB/h4 = (2x + a2/h2)(2y + b2/h2)
# We need (2x + a2/h2)(2y + b2/h2) >= 1

# Further: a2/h2 = (ha^2/h2^2) + a2_perp/h2 = x^2*h2 + a2_perp/h2
# Actually a2 = ha^2/h2 + a2_perp, so a2/h2 = x^2 + a2_perp/h2
# Hmm, let's use simpler notation.

# Let p = ha/h2, q = hb/h2 (projections of alpha, beta onto h-hat)
# Let s = a2_perp/h2, t = b2_perp/h2 (perpendicular parts normalized)
# Then a2/h2 = p^2 + s, b2/h2 = q^2 + t
# A/h2 = 2p + p^2 + s = (1+p)^2 - 1 + s
# B/h2 = 2q + q^2 + t = (1+q)^2 - 1 + t
# AB/h4 = [(1+p)^2 + s - 1][(1+q)^2 + t - 1]

# Even simpler: let P = 1+p = 1 + ha/h2 = (h2 + ha)/h2 = <h + alpha, h>/||h||^2 = <u, h>/||h||^2
# Similarly Q = 1+q = <v, h>/||h||^2
# Then A/h2 = P^2 + s - 1, B/h2 = Q^2 + t - 1
# AB/h4 = (P^2 + s - 1)(Q^2 + t - 1)

print("\n\nSchur complement analysis:")
print("Let P = <u,h>/||h||^2, Q = <v,h>/||h||^2")
print("    s = ||alpha_perp||^2/||h||^2, t = ||beta_perp||^2/||h||^2")
print("Then AB/h4 = (P^2 + s - 1)(Q^2 + t - 1)")
print("Need: (P^2 + s - 1)(Q^2 + t - 1) >= 1\n")

P_vals = []
Q_vals = []
s_vals = []
t_vals = []
for d in n3_data:
    P = np.dot(d['u'], d['h']) / d['h2']
    Q = np.dot(d['v'], d['h']) / d['h2']
    s = (d['a2'] - d['ha']**2/d['h2']) / d['h2']
    t = (d['b2'] - d['hb']**2/d['h2']) / d['h2']
    P_vals.append(P)
    Q_vals.append(Q)
    s_vals.append(s)
    t_vals.append(t)

print(f"  P = <u,h>/||h||^2: min={min(P_vals):.6f}, max={max(P_vals):.6f}, mean={np.mean(P_vals):.6f}")
print(f"  Q = <v,h>/||h||^2: min={min(Q_vals):.6f}, max={max(Q_vals):.6f}, mean={np.mean(Q_vals):.6f}")
print(f"  s = ||alpha_perp||^2/||h||^2: min={min(s_vals):.6f}, max={max(s_vals):.6f}")
print(f"  t = ||beta_perp||^2/||h||^2: min={min(t_vals):.6f}, max={max(t_vals):.6f}")
print(f"  P^2 + s - 1: min={min(P**2 + s - 1 for P, s in zip(P_vals, s_vals)):.6f}")
print(f"  Q^2 + t - 1: min={min(Q**2 + t - 1 for Q, t in zip(Q_vals, t_vals)):.6f}")

# Check: is P always >= 1?
print(f"\n  P >= 1 always: {all(P >= 1 - 1e-10 for P in P_vals)}")
print(f"  Q >= 1 always: {all(Q >= 1 - 1e-10 for Q in Q_vals)}")
# If P >= 1 and s >= 0, then P^2 + s - 1 >= 0, so A/h2 >= 0 (which we know).
# But we need the PRODUCT >= 1.


# ================================================================
# PART 3b: DETERMINANTAL / SCHUR COMPLEMENT INEQUALITY
# ================================================================

print("\n\n" + "=" * 70)
print("PART 3b: DETERMINANTAL APPROACH")
print("=" * 70)

# Key identity: for the Gram matrix G of {h, alpha, beta}:
# det(G) = h2 * det(Schur) = h2 * (a2_perp * b2_perp - ab_perp^2)
# where ab_perp = <alpha_perp, beta_perp>

# CAUCHY-SCHWARZ on perpendicular parts:
# |ab_perp| <= sqrt(a2_perp * b2_perp)  =>  ab_perp^2 <= a2_perp * b2_perp  =>  det(G) >= 0

# Now write AB - h4 in the (P, Q, s, t, rho) parameterization where rho = ab_perp/sqrt(s*t*h2^2):
# AB - h4 = h4 * [(P^2 + s - 1)(Q^2 + t - 1) - 1]
# = h4 * [P^2*Q^2 + P^2*t + s*Q^2 + s*t - P^2 - Q^2 - s - t + 1 - 1]
# = h4 * [P^2*Q^2 + P^2*t + s*Q^2 + s*t - P^2 - Q^2 - s - t]
# = h4 * [P^2*(Q^2 + t - 1) + s*(Q^2 + t) - Q^2 - t]
# = h4 * [P^2*(Q^2 + t - 1) + (s - 1)*(Q^2 + t)]
# Hmm, this doesn't simplify cleanly.

# Let's try a DIFFERENT decomposition.
# Define: X = A/h2 = ||u||^2/||h||^2 - 1,  Y = B/h2 = ||v||^2/||h||^2 - 1
# Then AB/h4 = X*Y, and we need X*Y >= 1.
# Also X = 2*ha/h2 + a2/h2 >= 0, Y = 2*hb/h2 + b2/h2 >= 0.

# KEY APPROACH: Use the CAUCHY-SCHWARZ inequality in a clever way.
# <u, v> = <h + alpha, h + beta> = h2 + ha + hb + <alpha, beta>
# By C-S: <u,v>^2 <= ||u||^2 * ||v||^2 = (h2 + A)(h2 + B)
# So: (h2 + ha + hb + ab_ip)^2 <= (h2 + A)(h2 + B)
# And: (h2 + A)(h2 + B) = h2^2 + h2*(A+B) + AB
# So: AB >= (h2 + ha + hb + ab_ip)^2 - h4 - h2*(A+B)
# This gives a LOWER bound on AB, but we need AB >= h4.

# Check: is there a vector w such that <u,w> and <v,w> give us what we need?
# If w = h: <u,h> = h2 + ha, <v,h> = h2 + hb
# C-S: <u,h>^2 <= ||u||^2 * ||h||^2  =>  (h2+ha)^2 <= (h2+A)*h2
# => h4 + 2*h2*ha + ha^2 <= h2^2 + A*h2
# => 2*ha + ha^2/h2 <= A  (which is true since A = 2*ha + a2 >= 2*ha + ha^2/h2)

# ANOTHER APPROACH: AM-GM or weighted AM-GM
# AB = A*B. Can we show A + B >= 2*h2?  If so, by AM-GM: AB >= (A+B)^2/4 >= h4.
# But that needs A + B >= 2*h2.
# A + B = Phi_p + Phi_q - 2*Phi_r. Is this >= 2*Phi_r?
# i.e., Phi_p + Phi_q >= 4*Phi_r?

print("\nChecking AM-GM approach: is A + B >= 2*||h||^2?")
amgm_holds = 0
amgm_total = 0
for d in n3_data:
    amgm_total += 1
    if d['A'] + d['B'] >= 2 * d['h2'] - 1e-10:
        amgm_holds += 1

print(f"  A + B >= 2*h2: {amgm_holds}/{amgm_total}")
if amgm_holds < amgm_total:
    # Find worst case
    worst = min(n3_data, key=lambda d: d['A'] + d['B'] - 2*d['h2'])
    print(f"  Worst: A+B - 2h2 = {worst['A'] + worst['B'] - 2*worst['h2']:.10f}")
    print(f"  So AM-GM does {'NOT' if worst['A'] + worst['B'] - 2*worst['h2'] < -1e-8 else ''} work!")

# Check: sqrt(AB) >= h2? Equivalent to AB >= h4.
# Is min(A,B) >= h2? If so, trivially AB >= h4.
print("\nChecking if min(A,B) >= ||h||^2:")
minAB_holds = sum(1 for d in n3_data if min(d['A'], d['B']) >= d['h2'] - 1e-10)
print(f"  min(A,B) >= h2: {minAB_holds}/{len(n3_data)}")


# ================================================================
# PART 4: FOUR-VECTORS APPROACH: h, alpha, beta, gamma = alpha - beta
# ================================================================

print("\n\n" + "=" * 70)
print("PART 4: FOUR-VECTORS ANALYSIS")
print("=" * 70)

# gamma = alpha - beta = u - v (difference of H-vectors)
# Note: gamma_k = H_p(lambda_k) - H_q(mu_k) depends on both p and q.
# gamma is the "asymmetry vector": measures how different p and q are.

# Structural constraints:
# 1. alpha + beta = u + v - 2h (sum constraint)
# 2. For n=2: alpha and beta are parallel to h (and to each other). Exact equality.
# 3. For n >= 3: alpha_perp, beta_perp are nonzero.

# Check: is gamma always parallel to h?
print("\nCorrelation of gamma = alpha - beta with h:")
gamma_h_corr = []
for d in n3_data:
    gamma = d['alpha'] - d['beta']
    if np.linalg.norm(gamma) > 1e-10 and np.linalg.norm(d['h']) > 1e-10:
        corr = np.dot(gamma, d['h']) / (np.linalg.norm(gamma) * np.linalg.norm(d['h']))
        gamma_h_corr.append(corr)

print(f"  cos(gamma, h): min={min(gamma_h_corr):.6f}, max={max(gamma_h_corr):.6f}, mean={np.mean(gamma_h_corr):.6f}")

# Check: is alpha + beta always parallel to h?
print("\nCorrelation of alpha + beta with h:")
apb_h_corr = []
for d in n3_data:
    apb = d['alpha'] + d['beta']
    if np.linalg.norm(apb) > 1e-10 and np.linalg.norm(d['h']) > 1e-10:
        corr = np.dot(apb, d['h']) / (np.linalg.norm(apb) * np.linalg.norm(d['h']))
        apb_h_corr.append(corr)

print(f"  cos(alpha+beta, h): min={min(apb_h_corr):.6f}, max={max(apb_h_corr):.6f}, mean={np.mean(apb_h_corr):.6f}")

# Decompose the excess AB - h4 into contributions:
# AB - h4 = (2ha + a2)(2hb + b2) - h4
# = 4*ha*hb + 2*ha*b2 + 2*a2*hb + a2*b2 - h4
#
# Group 1: "parallel part": involves only projections onto h
#   4*ha*hb + 2*ha*(hb^2/h2) + 2*(ha^2/h2)*hb + (ha^2/h2)*(hb^2/h2) - h4
#   This is (2ha + ha^2/h2)(2hb + hb^2/h2) - h4
#   = h2^2 * (2p + p^2)(2q + q^2) - h4  where p = ha/h2, q = hb/h2
#   = h4 * [(1+p)^2 - 1][(1+q)^2 - 1] - h4  ... hmm, this is not the full story.

# Better: split a2 = ha^2/h2 + a2_perp (= ||alpha_parallel||^2 + ||alpha_perp||^2)
# Similarly b2 = hb^2/h2 + b2_perp
# Then:
# A = 2ha + ha^2/h2 + a2_perp = h2*(2p + p^2) + a2_perp = h2*((1+p)^2 - 1) + a2_perp
# B = h2*((1+q)^2 - 1) + b2_perp
# Define A_par = h2*((1+p)^2 - 1), A_perp = a2_perp (and similarly for B)
# Then AB = (A_par + A_perp)(B_par + B_perp)
#         = A_par*B_par + A_par*B_perp + A_perp*B_par + A_perp*B_perp
# h4 = (h2)^2
#
# A_par*B_par = h4 * ((1+p)^2 - 1)((1+q)^2 - 1)
# = h4 * (p^2 + 2p)(q^2 + 2q)
# = h4 * (pq + 2p)(pq + 2q) ... no, expand:
# = h4 * [p^2*q^2 + 2p*q^2 + 2p^2*q + 4pq]
# = h4 * pq(pq + 2p + 2q + 4)
# = h4 * pq * (p+2)(q+2)

print("\n\nDecomposition of AB - h4 into parallel and perpendicular parts:")
for idx in range(min(10, len(n3_data))):
    d = n3_data[idx]
    p_coeff = d['ha'] / d['h2']
    q_coeff = d['hb'] / d['h2']
    a2_perp = d['a2'] - d['ha']**2 / d['h2']
    b2_perp = d['b2'] - d['hb']**2 / d['h2']

    A_par = d['h2'] * ((1 + p_coeff)**2 - 1)
    A_perp = a2_perp
    B_par = d['h2'] * ((1 + q_coeff)**2 - 1)
    B_perp = b2_perp

    par_par = A_par * B_par
    par_perp = A_par * B_perp + A_perp * B_par
    perp_perp = A_perp * B_perp

    excess_check = par_par + par_perp + perp_perp - d['h4']

    if idx < 3:
        print(f"  Trial {idx}: par*par={par_par:.6f}, cross={par_perp:.6f}, perp*perp={perp_perp:.6f}, h4={d['h4']:.6f}")
        print(f"    AB-h4={d['excess']:.6f}, check={excess_check:.6f}")
        print(f"    p={p_coeff:.6f}, q={q_coeff:.6f}, s={a2_perp/d['h2']:.6f}, t={b2_perp/d['h2']:.6f}")

# Find the parallel-only contribution to excess
print("\nParallel-only excess: A_par*B_par - h4 = h4*[pq(p+2)(q+2) - 1]:")
par_only_negative = 0
for d in n3_data:
    p_coeff = d['ha'] / d['h2']
    q_coeff = d['hb'] / d['h2']
    pq_term = p_coeff * q_coeff * (p_coeff + 2) * (q_coeff + 2)
    if pq_term < 1 - 1e-10:
        par_only_negative += 1
print(f"  A_par*B_par < h4: {par_only_negative}/{len(n3_data)} cases")
print(f"  (So the perpendicular parts ARE needed to ensure AB >= h4.)")


# ================================================================
# PART 5: RATIO ANALYSIS AND STRUCTURAL CONSTRAINTS
# ================================================================

print("\n\n" + "=" * 70)
print("PART 5: RATIO (AB - h4)/h4 AS FUNCTION OF STRUCTURAL PARAMETERS")
print("=" * 70)

# Parameterize by:
# p = <h,alpha>/||h||^2,  q = <h,beta>/||h||^2
# s = ||alpha_perp||^2/||h||^2,  t = ||beta_perp||^2/||h||^2
# rho = <alpha_perp, beta_perp> / sqrt(||alpha_perp||^2 * ||beta_perp||^2)
# (angle between perpendicular parts)
#
# Then AB/h4 - 1 = (2p + p^2 + s)(2q + q^2 + t) - 1

print("\nFull parameterization (p, q, s, t, rho):")
data_pqst = []
for d in n3_data:
    p_coeff = d['ha'] / d['h2']
    q_coeff = d['hb'] / d['h2']
    a2_perp = d['a2'] - d['ha']**2 / d['h2']
    b2_perp = d['b2'] - d['hb']**2 / d['h2']
    s = a2_perp / d['h2']
    t = b2_perp / d['h2']

    ab_perp = d['ab_ip'] - d['ha'] * d['hb'] / d['h2']
    rho = ab_perp / (np.sqrt(a2_perp * b2_perp) + 1e-30) if a2_perp > 1e-15 and b2_perp > 1e-15 else 0

    data_pqst.append({
        'p': p_coeff, 'q': q_coeff, 's': s, 't': t, 'rho': rho,
        'excess_ratio': d['excess'] / d['h4']
    })

# What is the MINIMUM of excess_ratio = AB/h4 - 1?
sorted_by_excess = sorted(data_pqst, key=lambda x: x['excess_ratio'])
print("\n10 cases with smallest AB/h4 - 1:")
for item in sorted_by_excess[:10]:
    print(f"  AB/h4-1={item['excess_ratio']:.8f}, p={item['p']:.4f}, q={item['q']:.4f}, "
          f"s={item['s']:.4f}, t={item['t']:.4f}, rho={item['rho']:.4f}")

print("\n10 cases with largest AB/h4 - 1:")
for item in sorted_by_excess[-10:]:
    print(f"  AB/h4-1={item['excess_ratio']:.8f}, p={item['p']:.4f}, q={item['q']:.4f}, "
          f"s={item['s']:.4f}, t={item['t']:.4f}, rho={item['rho']:.4f}")


# ================================================================
# PART 6: CONSTRAINT DISCOVERY - STRUCTURAL RELATIONS
# ================================================================

print("\n\n" + "=" * 70)
print("PART 6: STRUCTURAL CONSTRAINTS FROM SUBORDINATION")
print("=" * 70)

# The key constraint comes from subordination:
# omega_1(z) + omega_2(z) = z + c  (for some constant c, the sum of subordination functions)
# At z = nu_k: omega_1(nu_k) = lambda_k, omega_2(nu_k) = mu_k
# So: lambda_k + mu_k = nu_k + c  for all k
# (This is a VERY STRONG constraint!)

print("\nChecking lambda_k + mu_k - nu_k = const:")
for trial_idx in range(min(15, len(n3_data))):
    d = n3_data[trial_idx]
    sums = d['roots_p'] + d['roots_q'] - d['roots_r']
    if trial_idx < 5:
        print(f"  Trial {trial_idx}: lambda+mu-nu = {sums}")
    # Check if constant
    if np.std(sums) < 1e-6:
        pass  # constant - expected

# Actually wait - the subordination functions satisfy omega_1 + omega_2 = id + c
# only in the free probability infinite limit. For FINITE n MSS convolution,
# the relationship is different.
# For MSS: the roots satisfy e_k(r) = sum_{i+j=k} ... not a simple root sum.
# Let's check numerically.

const_check = []
for d in n3_data:
    sums = d['roots_p'] + d['roots_q'] - d['roots_r']
    const_check.append(np.std(sums) / (np.mean(np.abs(sums)) + 1e-15))

print(f"\n  Std(lambda+mu-nu) / Mean: min={min(const_check):.6e}, max={max(const_check):.6e}")
print(f"  This is {'~constant' if max(const_check) < 0.01 else 'NOT constant'} (subordination constraint)")

# Check simpler: e1(r) = e1(p) + e1(q) (sum of roots)
print("\nChecking sum of roots: sum(nu) = sum(lambda) + sum(mu)?")
for trial_idx in range(5):
    d = n3_data[trial_idx]
    sum_r = np.sum(d['roots_r'])
    sum_p = np.sum(d['roots_p'])
    sum_q = np.sum(d['roots_q'])
    print(f"  sum(r)={sum_r:.6f}, sum(p)+sum(q)={sum_p+sum_q:.6f}, diff={sum_r - sum_p - sum_q:.2e}")


# ================================================================
# PART 7: SEEK AN IDENTITY FOR AB - h4
# ================================================================

print("\n\n" + "=" * 70)
print("PART 7: SEEKING AN IDENTITY FOR AB - h4")
print("=" * 70)

# For n=2: AB = h4 EXACTLY. Let's verify this algebraically.
# n=2: roots p = (a,b), q = (c,d), r = MSS(p,q)
# For n=2, boxplus gives roots (a+c+b+d)/2 +/- sqrt(((a-b)^2+(c-d)^2)/4)
# Actually: for n=2, the MSS formula gives:
# e1(r) = e1(p) + e1(q) = (a+b) + (c+d)
# hat_e2(r) = hat_e2(p) + hat_e2(q)  (finite free convolution of normalized coeffs)
# e2/C(2,2) = e2(p)/1 + e2(q)/1  (C(2,2)=1)
# Wait: hat_e_k = e_k / C(n,k). For n=2, k=2: hat_e_2 = e_2/1 = e_2.
# So e2(r) = e2(p) + e2(q) = ab + cd.
# r has roots: [(a+b+c+d)/2 +/- sqrt((a+b+c+d)^2/4 - ab - cd)]
# = [(a+b+c+d)/2 +/- sqrt((a-b)^2/4 + (c-d)^2/4 + (a+b)(c+d)/2 - ab - cd)]
# Hmm, let me just compute: x^2 - (a+b+c+d)x + (ab+cd)
# discriminant = (a+b+c+d)^2 - 4(ab+cd)
# = (a+b)^2 + 2(a+b)(c+d) + (c+d)^2 - 4ab - 4cd
# = (a-b)^2 + 2(a+b)(c+d) + (c-d)^2 - 4cd  ... wait
# = a^2 + b^2 - 2ab + c^2 + d^2 - 2cd + 2(a+b)(c+d)
# Hmm, let me just do it directly.

print("\nn=2 EXACT verification:")
print("For n=2: p=(a,b), q=(c,d), r has roots (a+b+c+d)/2 +/- D/2")
print("where D^2 = (a-b)^2 + (c-d)^2\n")

# For n=2:
# H_p(a) = 1/(a-b), H_p(b) = 1/(b-a)
# Phi_p = 1/(a-b)^2 + 1/(b-a)^2 = 2/(a-b)^2
# Similarly Phi_q = 2/(c-d)^2
# For r: gap = D = sqrt((a-b)^2 + (c-d)^2), Phi_r = 2/D^2
# A = Phi_p - Phi_r = 2/(a-b)^2 - 2/D^2
# B = Phi_q - Phi_r = 2/(c-d)^2 - 2/D^2
# h4 = Phi_r^2 = 4/D^4
# AB = [2/(a-b)^2 - 2/D^2][2/(c-d)^2 - 2/D^2]
#    = 4/[(a-b)^2(c-d)^2] - 4/[D^2(c-d)^2] - 4/[D^2(a-b)^2] + 4/D^4
#    = 4/[(a-b)^2(c-d)^2] - 4/[D^2] * [1/(c-d)^2 + 1/(a-b)^2] + 4/D^4
# Let s = (a-b)^2, t = (c-d)^2, D^2 = s + t
# AB = 4/(st) - 4/(s+t) * (s+t)/(st) + 4/(s+t)^2
#    = 4/(st) - 4/(st) + 4/(s+t)^2
#    = 4/(s+t)^2 = 4/D^4 = h4  EXACTLY!

print("PROVED: For n=2, AB = h4 = 4/(s+t)^2 where s=(a-b)^2, t=(c-d)^2.")
print("The intermediate step: AB = 4/(st) - 4/(st) + 4/(s+t)^2 = 4/(s+t)^2.")
print("This is a PERFECT CANCELLATION.\n")

# For n=3: the excess AB - h4 is NOT zero. What is it?
# Let's look at it numerically as a function of the gap structure.

print("For n=3, analyzing the excess as a function of gap structure:")
print("Gaps of p: (g1_p, g2_p) = (lambda_2-lambda_1, lambda_3-lambda_2)")
print("Gaps of q: (g1_q, g2_q)")
print("Gaps of r: (g1_r, g2_r)")

# Check if excess depends only on gaps (translation-invariant)
print("\nTranslation invariance check:")
for shift in [0, 1, 5, 100]:
    d = n3_data[0]
    shifted_p = d['roots_p'] + shift
    shifted_q = d['roots_q'] + shift
    result_shifted = compute_all_vectors(shifted_p, shifted_q)
    if result_shifted:
        print(f"  Shift={shift}: excess = {result_shifted['excess']:.10f}")

# Check scale dependence
print("\nScale dependence:")
for scale in [0.5, 1, 2, 5]:
    d = n3_data[0]
    scaled_p = d['roots_p'] * scale
    scaled_q = d['roots_q'] * scale
    result_scaled = compute_all_vectors(scaled_p, scaled_q)
    if result_scaled:
        print(f"  Scale={scale}: excess = {result_scaled['excess']:.10f}, h4 = {result_scaled['h4']:.10f}, "
              f"ratio = {result_scaled['excess']/result_scaled['h4']:.10f}")

# Under scaling by s: roots -> s*roots, gaps -> s*gaps, H -> H/s, Phi -> Phi/s^2
# A -> A/s^2, B -> B/s^2, AB -> AB/s^4, h4 -> h4/s^4
# So excess -> excess/s^4, and ratio AB/h4 is SCALE-INVARIANT.
# This means the ratio depends only on the SHAPE of the root configuration.

print("\n(AB/h4 is scale-invariant: depends only on shape of root configuration.)")


# ================================================================
# PART 8: QUADRATIC FORM / MATRIX REPRESENTATION
# ================================================================

print("\n\n" + "=" * 70)
print("PART 8: QUADRATIC FORM REPRESENTATION")
print("=" * 70)

# AB - h4 in terms of the 5 Gram matrix entries (h2, ha, hb, a2, b2 + ab_ip):
# AB - h4 = 4*ha*hb + 2*ha*b2 + 2*a2*hb + a2*b2 - h2^2
#
# This is a DEGREE 2 polynomial in the Gram entries, which are themselves degree 2
# in the underlying vectors. So AB - h4 is degree 4 in the vectors.
#
# Can we write AB - h4 as a sum of squares in some natural basis?
#
# Key test: AB - h4 = 0 when n=2 (always). For n=3, the excess is a function of
# the root gaps. It should be degree -4 in the gaps (since H ~ 1/gap).

# Let's check: for n=3 equally spaced p and q, is there always equality?
print("\nn=3 equally spaced (p and q both have equal gaps):")
for s_p in [1.0, 2.0, 3.0]:
    for s_q in [1.0, 2.0, 3.0]:
        roots_p = np.array([-s_p, 0, s_p])
        roots_q = np.array([-s_q, 0, s_q])
        result = compute_all_vectors(roots_p, roots_q)
        if result:
            print(f"  gaps_p={s_p}, gaps_q={s_q}: AB/h4={result['ratio']:.10f}, excess={result['excess']:.6e}")

# Check: is equality when BOTH p and q are equally spaced?
print("\nn=3 p equally spaced, q NOT equally spaced:")
for trial in range(5):
    roots_p = np.array([-2.0, 0.0, 2.0])
    roots_q = random_well_separated_roots(3)
    result = compute_all_vectors(roots_p, roots_q)
    if result:
        print(f"  q={roots_q}: AB/h4={result['ratio']:.10f}")


# ================================================================
# PART 9: EXPLORE IF AB - h4 CAN BE A DETERMINANT
# ================================================================

print("\n\n" + "=" * 70)
print("PART 9: DETERMINANTAL STRUCTURE")
print("=" * 70)

# Since the Gram matrix G = [[h2, ha, hb], [ha, a2, ab], [hb, ab, b2]] is PSD,
# all its 2x2 minors are nonneg:
# h2*a2 >= ha^2,  h2*b2 >= hb^2,  a2*b2 >= ab^2
#
# AB - h4 = 4*ha*hb + 2*ha*b2 + 2*a2*hb + a2*b2 - h2^2
#
# Try: is AB - h4 = det(M) for some PSD matrix M that's a function of G?
#
# Alternatively: AB - h4 = (a2 + 2*ha)(b2 + 2*hb) - h2^2
# Let's define: X = sqrt(A) = sqrt(a2 + 2ha), Y = sqrt(B) = sqrt(b2 + 2hb)
# Then AB - h4 = X^2*Y^2 - h2^2 = (XY - h2)(XY + h2)
# We need XY >= h2, i.e., sqrt(A)*sqrt(B) >= ||h||^2.
# This is sqrt(Phi_p - Phi_r) * sqrt(Phi_q - Phi_r) >= Phi_r.

# Can we use the FISHER INFORMATION interpretation?
# Phi_r is Fisher information of the spectral measure of r.
# The inequality 1/Phi_r >= 1/Phi_p + 1/Phi_q is equivalent to
# Phi_r <= Phi_p * Phi_q / (Phi_p + Phi_q) (harmonic mean).
# Cross-multiplying: Phi_r * (Phi_p + Phi_q) <= Phi_p * Phi_q
# <=> Phi_r*Phi_p + Phi_r*Phi_q <= Phi_p*Phi_q
# <=> Phi_p*Phi_q - Phi_r*Phi_p - Phi_r*Phi_q >= 0
# <=> (Phi_p - Phi_r)*(Phi_q - Phi_r) >= Phi_r^2  [adding Phi_r^2 to both sides of... no]
#
# Actually: Phi_p*Phi_q - Phi_r*Phi_p - Phi_r*Phi_q + Phi_r^2 = (Phi_p - Phi_r)(Phi_q - Phi_r) = AB
# And Phi_r^2 = h4.
# So the inequality AB >= h4 is EQUIVALENT to:
# Phi_p*Phi_q - Phi_r*(Phi_p + Phi_q) + Phi_r^2 >= Phi_r^2
# <=> Phi_p*Phi_q >= Phi_r*(Phi_p + Phi_q)
# <=> 1/Phi_r >= (Phi_p + Phi_q)/(Phi_p*Phi_q) = 1/Phi_q + 1/Phi_p
# YES! That's exactly the superadditivity of 1/Phi.

print("\nConfirming: AB >= h4 <=> 1/Phi_r >= 1/Phi_p + 1/Phi_q (Fisher superadditivity)")
for d in n3_data[:5]:
    lhs = 1.0 / d['h2']
    rhs = 1.0 / d['u2'] + 1.0 / d['v2']
    print(f"  1/Phi_r = {lhs:.8f}, 1/Phi_p + 1/Phi_q = {rhs:.8f}, diff = {lhs - rhs:.8f}")


# ================================================================
# PART 10: CROSS-TERM ANALYSIS - <alpha, beta> STRUCTURE
# ================================================================

print("\n\n" + "=" * 70)
print("PART 10: CROSS-TERM <alpha, beta> ANALYSIS")
print("=" * 70)

# <alpha, beta> = sum_k alpha_k * beta_k = sum_k (u_k - h_k)(v_k - h_k)
# = sum_k u_k*v_k - sum_k u_k*h_k - sum_k h_k*v_k + sum_k h_k^2
# = <u,v> - <u,h> - <h,v> + h2

# AB = (h2 + 2ha + a2)(h2 + 2hb + b2)  ... no wait, A = 2ha + a2 (not h2 + 2ha + a2)
# A = ||u||^2 - ||h||^2 = u2 - h2
# B = v2 - h2
# AB = (u2 - h2)(v2 - h2)

# Using Cauchy-Schwarz: <u,v>^2 <= u2 * v2
# So u2*v2 >= <u,v>^2
# (u2 - h2)(v2 - h2) = u2*v2 - h2*(u2 + v2) + h4
#                     >= <u,v>^2 - h2*(u2 + v2) + h4

# Can we show <u,v>^2 - h2*(u2 + v2) + h4 >= h4?
# i.e., <u,v>^2 >= h2*(u2 + v2)
# By AM-GM: u2 + v2 >= 2*sqrt(u2*v2) >= 2*|<u,v>|
# So h2*(u2+v2) >= 2*h2*|<u,v>|
# And <u,v>^2 >= h2*(u2+v2) would need <u,v>^2 >= 2*h2*|<u,v>|
# i.e., |<u,v>| >= 2*h2. This is NOT always true.

# Different approach: check if there's a nice expression for <u,v>
print("\n<u,v> = sum_k H_p(lambda_k)*H_q(mu_k) analysis:")
uv_vals = []
for d in n3_data:
    uv = np.dot(d['u'], d['v'])
    uv_vals.append(uv)

print(f"  <u,v>: min={min(uv_vals):.6f}, max={max(uv_vals):.6f}, mean={np.mean(uv_vals):.6f}")

# Is <u,v> always >= h2?
print(f"  <u,v> >= ||h||^2: {sum(1 for d, uv in zip(n3_data, uv_vals) if uv >= d['h2'] - 1e-10)}/{len(n3_data)}")

# Is <u,v> always >= 0?
print(f"  <u,v> >= 0: {sum(1 for uv in uv_vals if uv >= -1e-10)}/{len(uv_vals)}")


# ================================================================
# PART 11: n=3 SPECIFIC STRUCTURE
# ================================================================

print("\n\n" + "=" * 70)
print("PART 11: n=3 SPECIFIC STRUCTURE (DEEP DIVE)")
print("=" * 70)

# For n=3, h, alpha, beta are in R^3.
# The Gram matrix G of {h, alpha, beta} is at most 3x3 rank-3.
# If det(G) > 0, the three vectors span R^3.

# Key observation: for n=3 equally spaced roots of BOTH p and q,
# AB = h4 (equality). Let's understand WHY.

print("\nFor equally spaced p=(-s,0,s), q=(-t,0,t):")
print("Analyzing the vectors h, alpha, beta:\n")

for s_val, t_val in [(1.0, 1.0), (1.0, 2.0), (2.0, 3.0)]:
    roots_p = np.array([-s_val, 0.0, s_val])
    roots_q = np.array([-t_val, 0.0, t_val])
    result = compute_all_vectors(roots_p, roots_q)
    if result:
        print(f"  s={s_val}, t={t_val}:")
        print(f"    h     = {result['h']}")
        print(f"    alpha = {result['alpha']}")
        print(f"    beta  = {result['beta']}")
        print(f"    alpha+beta = {result['alpha'] + result['beta']}")
        print(f"    alpha-beta = {result['alpha'] - result['beta']}")
        print(f"    alpha/h = {result['alpha']/result['h']}")
        print(f"    beta/h  = {result['beta']/result['h']}")
        print(f"    AB/h4 = {result['ratio']:.10f}")

        # Check if alpha is parallel to h
        if np.linalg.norm(result['h']) > 1e-10 and np.linalg.norm(result['alpha']) > 1e-10:
            cos_a = np.dot(result['alpha'], result['h']) / (np.linalg.norm(result['alpha']) * np.linalg.norm(result['h']))
            cos_b = np.dot(result['beta'], result['h']) / (np.linalg.norm(result['beta']) * np.linalg.norm(result['h']))
            print(f"    cos(alpha,h) = {cos_a:.10f}")
            print(f"    cos(beta,h)  = {cos_b:.10f}")
        print()

# The key question for equally spaced: alpha and beta should be parallel to h.
# If alpha = c_a * h and beta = c_b * h for some constants:
# Then A = (2*c_a + c_a^2)*h2, B = (2*c_b + c_b^2)*h2
# AB/h4 = (2*c_a + c_a^2)(2*c_b + c_b^2) = (c_a*(c_a+2))(c_b*(c_b+2))
# For this to equal 1, we need c_a*(c_a+2)*c_b*(c_b+2) = 1.
# Let's check if this holds for equally spaced.


# ================================================================
# PART 12: THE KEY OBSERVATION - WHAT MAKES AB > h4?
# ================================================================

print("\n" + "=" * 70)
print("PART 12: KEY OBSERVATION - PERPENDICULAR COMPONENTS")
print("=" * 70)

# When alpha and beta are perfectly parallel to h (i.e., alpha_perp = beta_perp = 0):
# AB/h4 = p*(p+2)*q*(q+2) where p = c_a, q = c_b
# For equally spaced: this equals 1 EXACTLY.
# For non-equally-spaced: there are perpendicular components.
# The perpendicular components ALWAYS increase AB.

# Let's verify: AB = (h2*p(p+2) + a2_perp)(h2*q(q+2) + b2_perp)
# = h4*p(p+2)*q(q+2) + h2*p(p+2)*b2_perp + h2*q(q+2)*a2_perp + a2_perp*b2_perp
# = h4 * [p(p+2)*q(q+2) + p(p+2)*t + q(q+2)*s + s*t]  (with s=a2p/h2, t=b2p/h2)
# = h4 * [(p(p+2) + s)(q(q+2) + t)]

# For equally spaced: s = t = 0 and p(p+2)*q(q+2) = 1.
# For non-equally-spaced: s, t > 0 and the parallel part might be < 1.
# But the total (parallel + perpendicular) is always >= 1.

# Let's compute the parallel contribution vs perpendicular contribution:
print("\nDecomposition: AB/h4 = parallel_term + cross_terms + perp_term")
print("parallel = p(p+2)*q(q+2),  cross = p(p+2)*t + q(q+2)*s,  perp = s*t")
print("where p = ha/h2, q = hb/h2, s = ||alpha_perp||^2/h2, t = ||beta_perp||^2/h2\n")

parallel_terms = []
cross_terms = []
perp_terms = []
for d in n3_data:
    p_c = d['ha'] / d['h2']
    q_c = d['hb'] / d['h2']
    s = max(0, d['a2'] - d['ha']**2 / d['h2']) / d['h2']
    t = max(0, d['b2'] - d['hb']**2 / d['h2']) / d['h2']

    par = p_c * (p_c + 2) * q_c * (q_c + 2)
    cross = p_c * (p_c + 2) * t + q_c * (q_c + 2) * s
    perp = s * t

    parallel_terms.append(par)
    cross_terms.append(cross)
    perp_terms.append(perp)

print(f"  parallel: min={min(parallel_terms):.6f}, max={max(parallel_terms):.6f}, mean={np.mean(parallel_terms):.6f}")
print(f"  cross:    min={min(cross_terms):.6f}, max={max(cross_terms):.6f}, mean={np.mean(cross_terms):.6f}")
print(f"  perp:     min={min(perp_terms):.6f}, max={max(perp_terms):.6f}, mean={np.mean(perp_terms):.6f}")
print(f"  total=par+cross+perp: should all be >= 1")
totals = [p + c + pp for p, c, pp in zip(parallel_terms, cross_terms, perp_terms)]
print(f"  total:    min={min(totals):.10f}")

# Is the parallel term alone always >= 1?
print(f"  parallel >= 1: {sum(1 for p in parallel_terms if p >= 1 - 1e-10)}/{len(parallel_terms)}")
# When parallel < 1, how much does cross+perp compensate?
deficit_cases = [(p, c, pp) for p, c, pp in zip(parallel_terms, cross_terms, perp_terms) if p < 1 - 1e-10]
if deficit_cases:
    print(f"\n  Cases where parallel < 1: {len(deficit_cases)}")
    worst_deficit = min(deficit_cases, key=lambda x: x[0])
    print(f"  Worst: parallel={worst_deficit[0]:.6f}, cross={worst_deficit[1]:.6f}, perp={worst_deficit[2]:.6f}")
    print(f"  Total = {sum(worst_deficit):.6f}")

# ================================================================
# PART 13: ATTEMPT DIRECT ALGEBRAIC IDENTITY
# ================================================================

print("\n\n" + "=" * 70)
print("PART 13: DIRECT ALGEBRAIC IDENTITY ATTEMPT")
print("=" * 70)

# Let's try to express AB - h4 as a sum of nonneg terms.
# AB - h4 = (A)(B) - h4 = (u2 - h2)(v2 - h2) - h4
#         = u2*v2 - h2*u2 - h2*v2 + h4 - h4
#         = u2*v2 - h2*(u2 + v2)
#
# So AB - h4 = ||u||^2*||v||^2 - ||h||^2*(||u||^2 + ||v||^2)
#
# For this to be >= 0 we need: ||u||^2*||v||^2 >= ||h||^2*(||u||^2 + ||v||^2)
# Dividing by ||u||^2*||v||^2: 1 >= ||h||^2*(1/||v||^2 + 1/||u||^2)
# i.e., 1/||h||^2 >= 1/||u||^2 + 1/||v||^2 (harmonic mean inequality for Phi!)

print("\nAB - h4 = ||u||^2*||v||^2 - ||h||^2*(||u||^2 + ||v||^2)")
print("This is the Fisher superadditivity: 1/Phi_r >= 1/Phi_p + 1/Phi_q\n")

# Can we write u2*v2 - h2*(u2+v2) = sum_i (positive terms)?
# u2*v2 - h2*u2 - h2*v2
# = u2*(v2 - h2) - h2*v2
# ... not obviously positive term by term.

# Try: u2*v2 - h2*(u2+v2) = (u2 - h2)(v2 - h2) - h4  ... that's what we started with.

# Actually let's try a totally different approach.
# Let D_kl = 1/(nu_k - nu_l), E_kl = 1/(lambda_k - lambda_l), F_kl = 1/(mu_k - mu_l)
# Then h_k = sum_{l!=k} D_kl, u_k = sum_{l!=k} E_kl, v_k = sum_{l!=k} F_kl
#
# Phi_r = sum_k (sum_{l!=k} D_kl)^2 = sum_k sum_{l!=k} sum_{m!=k} D_kl*D_km
# = sum_{k,l,m: l!=k, m!=k} D_kl * D_km
# = sum_{k,l: l!=k} D_kl^2 + sum_{k,l,m: l!=k,m!=k,l!=m} D_kl*D_km
#
# This is getting complicated. Let me try a matrix approach instead.

# Define the n x n matrix M_r[k,l] = 1/(nu_k - nu_l) for k != l, 0 on diagonal.
# Then h = M_r * 1 (where 1 is the all-ones vector).
# Phi_r = ||h||^2 = 1^T M_r^T M_r 1 = 1^T M_r^2 1 (since M_r is antisymmetric: M_r^T = -M_r)
# Wait: M_r[k,l] = 1/(nu_k - nu_l) and M_r[l,k] = 1/(nu_l - nu_k) = -1/(nu_k - nu_l) = -M_r[k,l]
# So M_r is ANTISYMMETRIC. M_r^T = -M_r.
# h = M_r * 1
# ||h||^2 = h^T h = 1^T M_r^T M_r 1 = -1^T M_r M_r 1 = -1^T M_r^2 1
# (Since M_r^T = -M_r, M_r^T M_r = -M_r^2)

print("Matrix approach: M_r[k,l] = 1/(nu_k - nu_l) (antisymmetric)")
print("h = M_r * e, where e = (1,...,1)")
print("Phi_r = -e^T M_r^2 e (since M_r antisymmetric)")

# Let's verify this
for d in n3_data[:3]:
    n = d['n']
    e_vec = np.ones(n)

    # Build M_r
    M_r = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                M_r[i, j] = 1.0 / (d['roots_r'][i] - d['roots_r'][j])

    h_check = M_r @ e_vec
    Phi_check = -e_vec @ M_r @ M_r @ e_vec

    print(f"  h = M_r @ e: {h_check}, actual: {d['h']}, match: {np.allclose(h_check, d['h'])}")
    print(f"  Phi_r = -e^T M_r^2 e: {Phi_check:.8f}, actual: {d['h2']:.8f}, match: {np.allclose(Phi_check, d['h2'])}")


# ================================================================
# PART 14: PRECISE STRUCTURE FOR EQUALLY SPACED
# ================================================================

print("\n\n" + "=" * 70)
print("PART 14: WHY EQUALITY FOR EQUALLY SPACED ROOTS?")
print("=" * 70)

# For equally spaced p = (-s, 0, s), q = (-t, 0, t):
# H_p(-s) = 1/(-s-0) + 1/(-s-s) = -1/s - 1/(2s) = -3/(2s)
# H_p(0) = 1/(0-(-s)) + 1/(0-s) = 1/s - 1/s = 0
# H_p(s) = 1/(s-(-s)) + 1/(s-0) = 1/(2s) + 1/s = 3/(2s)
# So u = (-3/(2s), 0, 3/(2s))
# Similarly v = (-3/(2t), 0, 3/(2t))
#
# For r = MSS(p,q): let's compute the roots.
# e1(r) = e1(p) + e1(q) = 0 + 0 = 0  (symmetric about 0)
# For the convolution of (-s,0,s) and (-t,0,t):
# e1(p) = 0, e2(p) = -s^2, e3(p) = 0
# e1(q) = 0, e2(q) = -t^2, e3(q) = 0
# a_0 = 1, a_1 = 0, a_2 = -s^2, a_3 = 0  (coeffs of x^3 - 0*x^2 - s^2*x + 0)
# b_0 = 1, b_1 = 0, b_2 = -t^2, b_3 = 0
#
# c_k = sum_{i+j=k} [(3-i)!(3-j)! / (3!(3-k)!)] * a_i * b_j
# c_0 = (3!*3!)/(3!*3!) * 1 * 1 = 1
# c_1 = [(3!*2!)/(3!*2!)] * 1*0 + [(2!*3!)/(3!*2!)] * 0*1 = 0
# c_2 = [(3!*1!)/(3!*1!)] * 1*(-t^2) + [(2!*2!)/(3!*1!)] * 0*0 + [(1!*3!)/(3!*1!)] * (-s^2)*1
#      = -t^2 + 0 + (-s^2) * (1*6)/(6*1) = -t^2 - s^2
#      Wait, let me recompute: [(1!*3!)/(3!*1!)] = (1*6)/(6*1) = 1. So c_2 = -t^2 - s^2.
# c_3 = sum ... a_0*b_3 + a_1*b_2 + a_2*b_1 + a_3*b_0
#      = [(3!*0!)/(3!*0!)] * 1 * 0 + [(2!*1!)/(3!*0!)] * 0 * (-t^2) + [(1!*2!)/(3!*0!)] * (-s^2)*0 + [(0!*3!)/(3!*0!)] * 0 * 1
#      = 0

# So r(x) = x^3 - (s^2+t^2)*x = x(x^2 - (s^2+t^2))
# Roots: 0, +/- sqrt(s^2 + t^2) = +/- D where D = sqrt(s^2+t^2)
# So r = (-D, 0, D), equally spaced with gap D.

# h_r(-D) = -3/(2D), h_r(0) = 0, h_r(D) = 3/(2D)
# Phi_r = 9/(4D^2) + 0 + 9/(4D^2) = 9/(2D^2) = 9/(2(s^2+t^2))
# Phi_p = 9/(2s^2), Phi_q = 9/(2t^2)
# A = 9/(2s^2) - 9/(2(s^2+t^2)) = (9/2) * t^2/(s^2(s^2+t^2))
# B = (9/2) * s^2/(t^2(s^2+t^2))
# AB = (81/4) * s^2*t^2 / (s^2*t^2*(s^2+t^2)^2) = (81/4) / (s^2+t^2)^2
# h4 = [9/(2(s^2+t^2))]^2 = 81/(4(s^2+t^2)^2)
# AB = h4 EXACTLY!

print("VERIFIED ALGEBRAICALLY for n=3 equally spaced:")
print("  r = p boxplus q has roots (-D, 0, D) where D = sqrt(s^2 + t^2)")
print("  Phi_p = 9/(2s^2), Phi_q = 9/(2t^2), Phi_r = 9/(2D^2)")
print("  A = (9/2)*t^2/(s^2*D^2), B = (9/2)*s^2/(t^2*D^2)")
print("  AB = 81/(4*D^4) = h4 exactly.")
print("  This is a PYTHAGOREAN STRUCTURE: D^2 = s^2 + t^2.")

# Now check: what happens with UNEQUAL spacing?
# p = (a, a+g1, a+g1+g2) with g1 != g2
print("\nFor n=3 with unequal gaps:")
for g1, g2, h1, h2_gap in [(1,1,1,1), (1,2,1,1), (2,1,1,2), (1,3,2,1)]:
    roots_p = np.array([0, g1, g1+g2], dtype=float)
    roots_q = np.array([0, h1, h1+h2_gap], dtype=float)
    result = compute_all_vectors(roots_p, roots_q)
    if result:
        print(f"  p-gaps=({g1},{g2}), q-gaps=({h1},{h2_gap}): AB/h4={result['ratio']:.10f}")


# ================================================================
# PART 15: COMPREHENSIVE SUMMARY STATISTICS
# ================================================================

print("\n\n" + "=" * 70)
print("PART 15: COMPREHENSIVE SUMMARY")
print("=" * 70)

# Check larger n
for nn in [3, 4, 5, 6]:
    print(f"\nn={nn}:")
    data = []
    for trial in range(500):
        roots_p = random_well_separated_roots(nn)
        roots_q = random_well_separated_roots(nn)
        result = compute_all_vectors(roots_p, roots_q)
        if result is not None:
            data.append(result)

    if data:
        excesses = [d['excess'] for d in data]
        ratios = [d['ratio'] for d in data]
        ha_vals = [d['ha'] for d in data]
        hb_vals = [d['hb'] for d in data]

        print(f"  Valid trials: {len(data)}")
        print(f"  AB/h4: min={min(ratios):.8f}, mean={np.mean(ratios):.4f}, max={max(ratios):.4f}")
        print(f"  AB >= h4: {all(e >= -1e-8 for e in excesses)}")
        print(f"  <h,alpha> < 0: {sum(1 for x in ha_vals if x < -1e-10)}/{len(data)}")
        print(f"  <h,beta>  < 0: {sum(1 for x in hb_vals if x < -1e-10)}/{len(data)}")

# Targeted search for counterexamples to <h,alpha> >= 0
print("\n\nTargeted search for <h,alpha> < 0:")
np.random.seed(123)
for nn in [3, 4, 5]:
    neg_ha_count = 0
    total = 0
    worst_ha = 0
    for trial in range(2000):
        roots_p = random_well_separated_roots(nn, scale=3.0, min_gap=0.3)
        roots_q = random_well_separated_roots(nn, scale=3.0, min_gap=0.3)
        result = compute_all_vectors(roots_p, roots_q)
        if result is not None:
            total += 1
            if result['ha'] < -1e-10:
                neg_ha_count += 1
                if result['ha'] < worst_ha:
                    worst_ha = result['ha']
    print(f"  n={nn}: <h,alpha><0 in {neg_ha_count}/{total} cases, worst={worst_ha:.6f}")


# ================================================================
# PART 16: THE STRUCTURAL DECOMPOSITION
# ================================================================

print("\n\n" + "=" * 70)
print("PART 16: STRUCTURAL DECOMPOSITION OF THE EXCESS")
print("=" * 70)

# We established: AB - h4 = u2*v2 - h2*(u2 + v2)
# = Phi_p * Phi_q - Phi_r * (Phi_p + Phi_q)
#
# Let's define F_p = 1/Phi_p, F_q = 1/Phi_q, F_r = 1/Phi_r (reciprocal Fisher info).
# Then: AB - h4 = 1/(F_p*F_q) - (1/F_r)*(1/F_p + 1/F_q)
# = 1/(F_p*F_q) - (F_p + F_q)/(F_r*F_p*F_q)
# = [F_r - (F_p + F_q)] / (F_r * F_p * F_q)  ... wait:
# = [F_r - F_p - F_q] / (F_r * F_p * F_q)  ... no, let me redo.
#
# AB - h4 = Phi_p*Phi_q - Phi_r*(Phi_p + Phi_q)
# Divide by Phi_p*Phi_q*Phi_r:
# (AB-h4)/(Phi_p*Phi_q*Phi_r) = 1/Phi_r - 1/Phi_p - 1/Phi_q = F_r - F_p - F_q
#
# So: AB - h4 = Phi_p * Phi_q * Phi_r * (F_r - F_p - F_q)
# And F_r - F_p - F_q >= 0 is the Fisher superadditivity.

print("CLEAN IDENTITY:")
print("  AB - h4 = Phi_p * Phi_q * Phi_r * (1/Phi_r - 1/Phi_p - 1/Phi_q)")
print("         = Phi_p * Phi_q * Phi_r * Delta")
print("  where Delta = 1/Phi_r - 1/Phi_p - 1/Phi_q is the Fisher deficit.")
print()
print("  So AB >= h4  <=>  Delta >= 0  <=>  1/Phi_r >= 1/Phi_p + 1/Phi_q")
print("  This is Fisher superadditivity of 1/Phi under MSS convolution.")
print()

# Verify
print("Verification:")
for d in n3_data[:5]:
    Delta = 1/d['h2'] - 1/d['u2'] - 1/d['v2']
    product = d['u2'] * d['v2'] * d['h2'] * Delta
    print(f"  AB-h4={d['excess']:.8f}, Phi_p*Phi_q*Phi_r*Delta={product:.8f}, match={abs(d['excess']-product)<1e-6}")


# ================================================================
# PART 17: NON-NEGATIVE EXPRESSION SEARCH
# ================================================================

print("\n\n" + "=" * 70)
print("PART 17: MANIFESTLY NON-NEGATIVE EXPRESSION SEARCH")
print("=" * 70)

# Since AB - h4 = u2*v2 - h2*(u2+v2), can we write it as a sum of squares
# in terms of components of u, v, h?
#
# Let's try: u2*v2 - h2*u2 - h2*v2
# = sum_i u_i^2 * sum_j v_j^2 - sum_i h_i^2 * sum_j u_j^2 - sum_i h_i^2 * sum_j v_j^2
# = sum_{i,j} u_i^2 * v_j^2 - sum_{i,j} h_i^2 * (u_j^2 + v_j^2)
#
# This doesn't factor nicely. Let's try Lagrange identity:
# ||u||^2*||v||^2 = <u,v>^2 + sum_{i<j} (u_i*v_j - u_j*v_i)^2  (by Lagrange)
# So: u2*v2 - h2*(u2+v2) = <u,v>^2 + sum_{i<j}(u_iv_j-u_jv_i)^2 - h2*(u2+v2)
#
# Or: u2*v2 = <u,v>^2 + ||u x v||^2  (Lagrange identity / Pythagorean)
#
# So AB - h4 = <u,v>^2 + ||u x v||^2 - h2*(u2+v2)
# = [<u,v>^2 - h2*(u2+v2)] + ||u x v||^2
#
# The second term ||u x v||^2 >= 0 always.
# The first term can be negative.
# So AB - h4 >= ||u x v||^2 - |<u,v>^2 - h2*(u2+v2)|
# Not directly useful.

# Try yet another approach: can we write u2*v2 - h2*(u2+v2) as sum of squares?
# u2*v2 - h2*u2 - h2*v2 + h4 = (u2-h2)(v2-h2) = AB  (just going in circles)

# Let's try to find a MANIFESTLY non-negative expression computationally.
# For n=3, we have 3-dimensional vectors u, v, h.
# u2*v2 - h2*(u2+v2)
# = (u2-h2)*v2 - h2*v2 + h2*v2 - h2*u2  ... no
# = (u2-h2)*v2 + h2*(v2-u2) ... if v2 >= u2, the second term is positive.

# Important attempt: write in terms of PAIRWISE differences.
# u_k - h_k = alpha_k, v_k - h_k = beta_k
# u_k = h_k + alpha_k, v_k = h_k + beta_k
#
# Consider the sum: sum_{k<l} (u_k*v_l - u_l*v_k)^2 = ||u x v||^2 (cross product squared for n=3)
# And sum_{k<l} (h_k*u_l - h_l*u_k)^2 = ||h x u||^2
# And sum_{k<l} (h_k*v_l - h_l*v_k)^2 = ||h x v||^2
#
# By the Lagrange identity:
# u2*v2 = <u,v>^2 + ||u x v||^2
# h2*u2 = <h,u>^2 + ||h x u||^2
# h2*v2 = <h,v>^2 + ||h x v||^2
#
# So: AB - h4 = u2*v2 - h2*(u2+v2)
# = <u,v>^2 + ||u x v||^2 - <h,u>^2 - ||h x u||^2 - <h,v>^2 - ||h x v||^2

print("\nLagrange identity decomposition:")
print("AB - h4 = <u,v>^2 + ||u x v||^2 - <h,u>^2 - ||h x u||^2 - <h,v>^2 - ||h x v||^2\n")

# Verify and check signs of individual terms
for d in n3_data[:5]:
    u, v, h = d['u'], d['v'], d['h']
    uv_sq = np.dot(u,v)**2
    uxv_sq = np.sum(np.cross(u,v)**2)
    hu_sq = np.dot(h,u)**2
    hxu_sq = np.sum(np.cross(h,u)**2)
    hv_sq = np.dot(h,v)**2
    hxv_sq = np.sum(np.cross(h,v)**2)

    total = uv_sq + uxv_sq - hu_sq - hxu_sq - hv_sq - hxv_sq
    print(f"  <u,v>^2={uv_sq:.4f}, ||uxv||^2={uxv_sq:.4f}, <h,u>^2={hu_sq:.4f}, "
          f"||hxu||^2={hxu_sq:.4f}, <h,v>^2={hv_sq:.4f}, ||hxv||^2={hxv_sq:.4f}")
    print(f"  Total={total:.6f}, AB-h4={d['excess']:.6f}, match={abs(total-d['excess'])<1e-4}")

# Group positive and negative:
print("\nGrouped: (positive: <u,v>^2 + ||uxv||^2) vs (negative: <h,u>^2 + ||hxu||^2 + <h,v>^2 + ||hxv||^2)")
for d in n3_data[:5]:
    u, v, h = d['u'], d['v'], d['h']
    pos = np.dot(u,v)**2 + np.sum(np.cross(u,v)**2)
    neg = np.dot(h,u)**2 + np.sum(np.cross(h,u)**2) + np.dot(h,v)**2 + np.sum(np.cross(h,v)**2)
    print(f"  pos={pos:.6f}, neg={neg:.6f}, diff={pos-neg:.6f}")

# This doesn't give a manifestly nonneg expression.

# Last attempt: GRAM DETERMINANT approach.
# For n=3, det(Gram(h,u,v)) = det(G') where
# G' = [[h2, <h,u>, <h,v>], [<h,u>, u2, <u,v>], [<h,v>, <u,v>, v2]]
# det(G') = h2*(u2*v2 - <u,v>^2) - <h,u>*(v2*<h,u> - <u,v>*<h,v>) + <h,v>*(<h,u>*<u,v> - u2*<h,v>)
# = h2*||uxv||^2 - <h,u>^2*v2 - <h,v>^2*u2 + 2*<h,u>*<h,v>*<u,v>   ... (standard formula)
#
# This is >= 0 (Gram det of 3 vectors in R^n with n>=3).
# Can we relate this to AB - h4?

print("\n\nGram determinant det(G') of (h,u,v):")
for d in n3_data[:5]:
    u, v, h = d['u'], d['v'], d['h']
    G_prime = np.array([
        [d['h2'], np.dot(h,u), np.dot(h,v)],
        [np.dot(h,u), d['u2'], np.dot(u,v)],
        [np.dot(h,v), np.dot(u,v), d['v2']],
    ])
    det_Gp = np.linalg.det(G_prime)
    print(f"  det(G') = {det_Gp:.8f}, AB-h4 = {d['excess']:.8f}")

# Schur complement of h2 in G':
# S' = [[u2, <u,v>], [<u,v>, v2]] - (1/h2)*[[<h,u>], [<h,v>]]*[[<h,u>, <h,v>]]
# = [[u2 - <h,u>^2/h2, <u,v> - <h,u><h,v>/h2],
#    [<u,v> - <h,u><h,v>/h2, v2 - <h,v>^2/h2]]
# = [[||u_perp||^2, <u_perp, v_perp>], [<u_perp, v_perp>, ||v_perp||^2]]
# where u_perp = u - (h^Tu/h2)*h, etc.
# det(S') = ||u_perp||^2 * ||v_perp||^2 - <u_perp, v_perp>^2 >= 0

# AB - h4 = u2*v2 - h2*(u2+v2)
# det(G') = h2 * det(S')
# det(S') = ||u_perp||^2*||v_perp||^2 - <u_perp,v_perp>^2
#
# These are DIFFERENT quantities. Let's check the relationship.

print("\ndet(G') vs AB-h4 (these are different quantities):")
for d in n3_data[:5]:
    u, v, h = d['u'], d['v'], d['h']
    G_prime = np.array([
        [d['h2'], np.dot(h,u), np.dot(h,v)],
        [np.dot(h,u), d['u2'], np.dot(u,v)],
        [np.dot(h,v), np.dot(u,v), d['v2']],
    ])
    det_Gp = np.linalg.det(G_prime)
    print(f"  det(G')={det_Gp:.6f}, AB-h4={d['excess']:.6f}, ratio={d['excess']/(det_Gp+1e-30):.4f}")


# ================================================================
# FINAL: OPTIMIZATION OF MINIMUM AB/h4
# ================================================================

print("\n\n" + "=" * 70)
print("FINAL: GLOBAL MINIMUM OF AB/h4 FOR n=3")
print("=" * 70)

def neg_ratio_n3(params):
    """Minimize -(AB/h4) = find minimum of AB/h4."""
    a1, g1, g2, b1, d1, d2 = params
    # Ensure positive gaps
    g1, g2, d1, d2 = abs(g1) + 0.1, abs(g2) + 0.1, abs(d1) + 0.1, abs(d2) + 0.1
    roots_p = np.array([a1, a1 + g1, a1 + g1 + g2])
    roots_q = np.array([b1, b1 + d1, b1 + d1 + d2])

    result = compute_all_vectors(roots_p, roots_q)
    if result is None:
        return 0  # Return 0 instead of big number to avoid false minima
    return -result['ratio']

best = 0
best_params = None
np.random.seed(42)
for trial in range(500):
    x0 = np.random.randn(6) * 2
    x0[1:3] = np.abs(x0[1:3]) + 0.5
    x0[4:6] = np.abs(x0[4:6]) + 0.5
    try:
        res = minimize(neg_ratio_n3, x0, method='Nelder-Mead',
                      options={'maxiter': 10000, 'xatol': 1e-12, 'fatol': 1e-14})
        if res.fun < best:
            best = res.fun
            best_params = res.x
    except:
        pass

if best < 0:
    print(f"Global minimum AB/h4 = {-best:.12f}")
    if best_params is not None:
        a1, g1, g2, b1, d1, d2 = best_params
        g1, g2, d1, d2 = abs(g1) + 0.1, abs(g2) + 0.1, abs(d1) + 0.1, abs(d2) + 0.1
        print(f"  p-gaps = ({g1:.6f}, {g2:.6f}), ratio = {g1/g2:.6f}")
        print(f"  q-gaps = ({d1:.6f}, {d2:.6f}), ratio = {d1/d2:.6f}")

        roots_p = np.array([a1, a1 + g1, a1 + g1 + g2])
        roots_q = np.array([b1, b1 + d1, b1 + d1 + d2])
        result = compute_all_vectors(roots_p, roots_q)
        if result:
            print(f"  A/h2 = {result['A']/result['h2']:.8f}")
            print(f"  B/h2 = {result['B']/result['h2']:.8f}")
            print(f"  P = {np.dot(result['u'], result['h'])/result['h2']:.8f}")
            print(f"  Q = {np.dot(result['v'], result['h'])/result['h2']:.8f}")
            a2p = result['a2'] - result['ha']**2/result['h2']
            b2p = result['b2'] - result['hb']**2/result['h2']
            print(f"  s = ||alpha_perp||^2/h2 = {a2p/result['h2']:.8f}")
            print(f"  t = ||beta_perp||^2/h2  = {b2p/result['h2']:.8f}")
else:
    print("No valid minimum found")


print("\n\n" + "=" * 70)
print("OVERALL CONCLUSIONS")
print("=" * 70)
print("""
1. AB - ||h||^4 = Phi_p * Phi_q * Phi_r * (1/Phi_r - 1/Phi_p - 1/Phi_q)
   So AB >= ||h||^4 is EQUIVALENT to Fisher superadditivity: 1/Phi_r >= 1/Phi_p + 1/Phi_q.

2. For n=2: AB = ||h||^4 exactly (perfect cancellation via D^2 = s^2 + t^2).

3. For n=3 equally spaced: AB = ||h||^4 exactly (same Pythagorean structure).

4. For n=3 unequally spaced: AB > ||h||^4 strictly. The excess comes from
   perpendicular components of alpha, beta relative to h.

5. CRITICAL: <h,alpha> >= 0 is FALSE for n >= 3 (confirmed by counterexamples).
   The proof CANNOT rely on this assumption.

6. The decomposition AB = A_par*B_par + cross + perp shows that A_par*B_par
   can be < ||h||^4, but the perpendicular contributions always compensate.

7. No simple manifestly non-negative expression for AB - ||h||^4 was found
   purely from the algebra of h, alpha, beta. The non-negativity appears to
   require STRUCTURAL PROPERTIES of the MSS convolution (how roots of r
   relate to roots of p, q).

8. The ratio AB/||h||^4 is scale-invariant and translation-invariant.
   Its minimum over n=3 configurations approaches 1 (achieved at equally
   spaced roots).

9. The key structural identity is:
   AB - h4 = ||u||^2*||v||^2 - ||h||^2*(||u||^2 + ||v||^2)
   where u = H_p, v = H_q, h = H_r are the H-vectors.
""")
