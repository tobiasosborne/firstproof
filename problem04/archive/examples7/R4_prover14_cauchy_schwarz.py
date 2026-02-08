"""
PROVER-14 Part 3: Cauchy-Schwarz / Information-theoretic approach

KEY REALIZATION: The inequality 1/Phi_r >= 1/Phi_p + 1/Phi_q is formally
identical to the CRAMER-RAO bound in information geometry:

1/I(theta_1 + theta_2) >= 1/I(theta_1) + 1/I(theta_2)

where I(theta) is the Fisher information. This holds when the "estimators"
are independent. The Fisher information of the sum is bounded by the
harmonic mean of the individual Fisher informations.

PLAN:
1. Express Phi_n as a Fisher information
2. Show that MSS convolution corresponds to addition of independent parameters
3. Apply the information-theoretic inequality

APPROACH A: Phi_n as Fisher information of the empirical measure
  mu_p = (1/n) sum delta_{lambda_i}
  Phi_n / n^2 approximates Phi*(mu) = int (H mu)^2 d mu (Voiculescu)

APPROACH B: Direct Cauchy-Schwarz on the derivative of the MSS formula
  Since (p boxplus_n q)' = n (p^{(1)} boxplus_{n-1} q^{(1)}),
  the Cauchy transform satisfies a recursive relation.
  Can we exploit this?

APPROACH C: Variational characterization of 1/Phi_n
  1/Phi_n = min over unit vectors u: (u^T H)^2 / ||H||^2 ... no
  Actually 1/Phi_n = 1/||H||^2 which is the inverse of a squared norm.

  Better: 1/||v||^2 = min_{||w||=1} 1/|<v,w>|^2  (not useful)

  1/||v||^2 is convex in v? NO, it's not convex.
  But <v, e_i>^2 / ||v||^2 = cos^2(angle) is bounded.

APPROACH D: Resolvent identity
  Key identity: If R_A(z) = (z-A)^{-1} and R_B(z) = (z-B)^{-1}, then
  R_{A+B}(z) - R_A(z) - R_B(z) involves cross terms.

  For A = diag(mu), the resolvent trace is G_p(z).
"""
import numpy as np
from itertools import combinations
from math import factorial
import sympy as sp

print("="*70)
print("PROVER-14 Part 3: CAUCHY-SCHWARZ / INFO-GEOMETRIC APPROACH")
print("="*70)

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

def roots_to_monic_coeffs(roots):
    n = len(roots)
    coeffs = [1.0]
    for k in range(1, n+1):
        ek = sum(np.prod(list(combo)) for combo in combinations(roots, k))
        coeffs.append((-1)**k * ek)
    return np.array(coeffs)

def mss_convolve_n(p_coeffs, q_coeffs, n):
    r_coeffs = np.zeros(n+1)
    for k in range(n+1):
        ck = 0
        for i in range(k+1):
            j = k - i
            if i <= n and j <= n:
                coeff = factorial(n-i) * factorial(n-j) / (factorial(n) * factorial(n-k))
                ck += coeff * p_coeffs[i] * q_coeffs[j]
        r_coeffs[k] = ck
    return r_coeffs

def mss_convolve_roots(roots_p, roots_q):
    n = len(roots_p)
    p_coeffs = roots_to_monic_coeffs(roots_p)
    q_coeffs = roots_to_monic_coeffs(roots_q)
    r_coeffs = mss_convolve_n(p_coeffs, q_coeffs, n)
    r_roots = np.roots(r_coeffs)
    r_roots = np.sort(np.real(r_roots))
    return r_roots

# =============================================================
# Part 1: Variational characterization
# =============================================================
print("\n--- Part 1: Variational characterization of 1/Phi ---\n")

print("""
Observation: 1/||v||^2 is not convex in v. But consider:

  1/Phi_n(p) = 1 / sum_i H_i^2

The inequality says this is superadditive under MSS convolution.

KEY INSIGHT: Consider the "dual" formulation.
By Cauchy-Schwarz: (sum_i a_i b_i)^2 <= (sum_i a_i^2)(sum_i b_i^2)

With a_i = H_p(mu_i) * w_i, b_i = 1/w_i for positive weights w_i:
  (sum_i H_p(mu_i))^2 <= (sum_i H_p(mu_i)^2 w_i^2)(sum_i 1/w_i^2)

But sum_i H_p(mu_i) = 0! So this gives 0 <= product, which is trivial.

Instead, try: For a fixed test function f,
  (sum_i H_p(mu_i) f(mu_i))^2 <= Phi_n(p) * sum_i f(mu_i)^2

This gives a LOWER BOUND on Phi_n(p):
  Phi_n(p) >= (sum_i H_p(mu_i) f(mu_i))^2 / sum_i f(mu_i)^2

And thus an UPPER BOUND on 1/Phi_n(p):
  1/Phi_n(p) <= sum_i f(mu_i)^2 / (sum_i H_p(mu_i) f(mu_i))^2
""")

# =============================================================
# Part 2: The ADMITTED C derivative structure
# =============================================================
print("\n--- Part 2: Derivative structure of MSS ---\n")

print("""
ADMITTED C: (p boxplus_n q)' = n * (p^{(1)} boxplus_{n-1} q^{(1)})

Let r = p boxplus_n q with roots lambda_1 < ... < lambda_n.
Let mu_1' < ... < mu_{n-1}' be the critical points of p (roots of p').
Let nu_1' < ... < nu_{n-1}' be the critical points of q (roots of q').
Let rho_1 < ... < rho_{n-1} be the critical points of r (roots of r').

Then: The polynomial with roots rho_1, ..., rho_{n-1} is
      (p^{(1)} boxplus_{n-1} q^{(1)})
(up to scaling, since r' = n * p^{(1)} boxplus_{n-1} q^{(1)}).

This means: the CRITICAL POINTS of r are obtained by MSS-convolving
the monic polynomials of critical points of p and q!

So MSS convolution is RECURSIVE on the critical-point tower:
  r = p boxplus_n q
  r' corresponds to p^{(1)} boxplus_{n-1} q^{(1)}
  r'' corresponds to p^{(2)} boxplus_{n-2} q^{(2)}
  ...

This recursive structure is CRUCIAL for understanding Phi_n.

Recall: H_p(lambda_i) = p''(lambda_i) / (2 p'(lambda_i))
      = n * sum of 1/(lambda_i - rho_j) ... wait, that's not right.

Actually, p'(lambda_i) / p''(lambda_i) ... let me think.
""")

# Let's compute H using critical points
print("Connection between H and critical points:")
roots = np.array([-2.0, 0.0, 2.0])
n = len(roots)
H_direct = H_values(roots)

# p'(x) = n * prod(x - mu_j') where mu_j' are critical points
# Find critical points
p_coeffs = roots_to_monic_coeffs(roots)
# p(x) = x^3 + a1 x^2 + a2 x + a3
# p'(x) = 3x^2 + 2a1 x + a2
p_prime_coeffs = [3*p_coeffs[0], 2*p_coeffs[1], p_coeffs[2]]
crit_points = np.sort(np.real(np.roots(p_prime_coeffs)))
print(f"  roots = {roots}")
print(f"  critical points = {crit_points}")

# H_p(lambda_i) = p''(lambda_i) / (2 p'(lambda_i))
# But also:
# p'(x) = n * prod(x - mu_j')
# So p'(lambda_i) = n * prod(lambda_i - mu_j')
# And log|p'(x)| = log(n) + sum_j log|x - mu_j'|
# [log|p'(x)|]' = sum_j 1/(x - mu_j')

# Also p''(x) = n * sum_j prod_{k!=j} (x - mu_k')
# So p''(lambda_i) / p'(lambda_i) = sum_j 1/(lambda_i - mu_j')

# Therefore: H_p(lambda_i) = (1/2) * p''(lambda_i)/p'(lambda_i)
#           = (1/2) * sum_j 1/(lambda_i - mu_j')

for i in range(n):
    H_via_crit = 0.5 * sum(1/(roots[i] - cp) for cp in crit_points)
    print(f"  H[{i}] = {H_direct[i]:.6f}, via crit pts = {H_via_crit:.6f}")

# Wait, this doesn't match. Let me recheck.
# p(x) = prod(x - lambda_i), p'(x) = sum_i prod_{j!=i} (x - lambda_j)
# At x = lambda_k: p'(lambda_k) = prod_{j!=k} (lambda_k - lambda_j)
# p''(x) = sum_i sum_{j!=i} prod_{m!=i,m!=j} (x - lambda_m)
# At x = lambda_k: p''(lambda_k) = 2 sum_{j!=k} prod_{m!=k,m!=j} (lambda_k - lambda_m)
# So p''(lambda_k) / p'(lambda_k) = 2 sum_{j!=k} prod_{m!=k,m!=j}(lambda_k-lambda_m) / prod_{j!=k}(lambda_k-lambda_j)
#                                  = 2 sum_{j!=k} 1/(lambda_k - lambda_j)
#                                  = 2 H_p(lambda_k)

print("\nSo H_p(lambda_i) = p''(lambda_i) / (2 p'(lambda_i)). Confirmed above.")

# Now use the critical point representation
# p'(x) = n * prod(x - mu_j'), so at x = lambda_i:
# p'(lambda_i) = n * prod(lambda_i - mu_j')
# p''(x) = n * d/dx prod(x - mu_j') = n * sum_j prod_{k!=j}(x - mu_k')
# p''(lambda_i) = n * sum_j prod_{k!=j}(lambda_i - mu_k')
# p''(lambda_i)/p'(lambda_i) = sum_j 1/(lambda_i - mu_j')
# H_p(lambda_i) = (1/2) sum_j 1/(lambda_i - mu_j')
print("\nSo H_p(lambda_i) = (1/2) sum_j 1/(lambda_i - mu_j')")
print("where mu_j' are the critical points of p.\n")

for i in range(n):
    H_via_crit = 0.5 * sum(1/(roots[i] - cp) for cp in crit_points)
    print(f"  H_p(lambda_{i}) = {H_direct[i]:.6f}, (1/2) sum 1/(lambda_i - mu_j') = {H_via_crit:.6f}")

print("\nHmm, these don't match. Let me recheck the identity.")
# Actually for n=3, p(x) = x^3 - 4x. p'(x) = 3x^2 - 4.
# At x = -2: p'(-2) = 12-4 = 8
# p''(x) = 6x, p''(-2) = -12
# H = p''/2p' = -12/16 = -0.75. Correct!
# Critical pts of p: 3x^2 - 4 = 0, x = +/- 2/sqrt(3) = +/- 1.1547
print(f"\n  p(x) = x^3 - 4x, p'(x) = 3x^2 - 4")
print(f"  crit pts = {crit_points}")
# sum 1/(-2 - mu_j') = 1/(-2 - 1.1547) + 1/(-2 + 1.1547) = 1/(-3.1547) + 1/(-0.8453)
# = -0.3170 + (-1.1830) = -1.5
# (1/2) * (-1.5) = -0.75. YES it matches!
print(f"  (1/2) * [1/(-2 - {crit_points[0]:.4f}) + 1/(-2 - {crit_points[1]:.4f})]")
print(f"  = (1/2) * [{1/(-2-crit_points[0]):.4f} + {1/(-2-crit_points[1]):.4f}]")
print(f"  = (1/2) * {1/(-2-crit_points[0]) + 1/(-2-crit_points[1]):.4f}")
print(f"  = {0.5 * (1/(-2-crit_points[0]) + 1/(-2-crit_points[1])):.4f}")
print(f"  H_p(-2) = {H_direct[0]:.4f}")
# They match!

# =============================================================
# Part 3: Phi_n via critical points
# =============================================================
print("\n\n--- Part 3: Phi_n via critical points ---\n")

print("""
IMPORTANT IDENTITY:
  H_p(lambda_i) = (1/2) sum_{j=1}^{n-1} 1/(lambda_i - mu_j')

where mu_j' are the critical points of p.

This means:
  Phi_n = sum_i H_p(lambda_i)^2
        = (1/4) sum_i [sum_j 1/(lambda_i - mu_j')]^2

This is a BILINEAR form in the "cauchy matrix" C_{ij} = 1/(lambda_i - mu_j').
Specifically, letting C be the n x (n-1) matrix with C_{ij} = 1/(lambda_i - mu_j'):

  H = (1/2) C * 1_{n-1}  (where 1_{n-1} is the all-ones vector)

  Phi_n = (1/4) 1_{n-1}^T C^T C 1_{n-1} = (1/4) ||C^T 1_n ... no.

Wait: H_i = (1/2) sum_j C_{ij}. So H = (1/2) C * 1.
  Phi_n = ||H||^2 = (1/4) ||C * 1||^2 = (1/4) 1^T C^T C 1.

The matrix C^T C is (n-1) x (n-1).

CAUCHY MATRIX PROPERTY:
The matrix C with C_{ij} = 1/(lambda_i - mu_j') is a CAUCHY MATRIX.
These have beautiful determinant formulas and are well-studied.

For our Cauchy matrix, by interlacing (lambda_1 < mu_1' < lambda_2 < ... < mu_{n-1}' < lambda_n),
all entries of C have definite signs.
""")

# Verify: roots interlace with critical points
print("Interlacing verification:")
for config in [np.array([-2., 0., 2.]), np.array([-3., -1., 1., 3.]),
               np.array([-4., -1., 0., 2., 5.])]:
    n = len(config)
    coeffs = roots_to_monic_coeffs(config)
    # p'(x) coefficients
    p_prime = [(n-k)*coeffs[k] for k in range(n)]
    crit = np.sort(np.real(np.roots(p_prime)))
    print(f"  roots: {config}")
    print(f"  crits: {crit}")
    interlaces = all(config[i] < crit[i] < config[i+1] for i in range(n-1))
    print(f"  interlaces: {interlaces}")
    print()

# =============================================================
# Part 4: Study the Cauchy matrix C and its norms
# =============================================================
print("\n--- Part 4: Cauchy matrix analysis ---\n")

np.random.seed(42)
for trial in range(3):
    n = 4
    roots = np.sort(np.random.randn(n) * 3)
    while np.min(np.diff(roots)) < 0.3:
        roots = np.sort(np.random.randn(n) * 3)

    coeffs = roots_to_monic_coeffs(roots)
    p_prime = [(n-k)*coeffs[k] for k in range(n)]
    crits = np.sort(np.real(np.roots(p_prime)))

    # Cauchy matrix
    C = np.zeros((n, n-1))
    for i in range(n):
        for j in range(n-1):
            C[i,j] = 1.0 / (roots[i] - crits[j])

    H_vec = 0.5 * C @ np.ones(n-1)
    H_direct_vec = H_values(roots)
    phi = np.sum(H_vec**2)

    print(f"  Trial {trial}: roots = {roots}")
    print(f"  crits = {crits}")
    print(f"  H (via C) = {H_vec}")
    print(f"  H (direct) = {H_direct_vec}")
    print(f"  Match: {np.allclose(H_vec, H_direct_vec)}")
    print(f"  Phi = {phi:.6f}")
    print(f"  singular values of C: {np.linalg.svd(C, compute_uv=False)}")
    print(f"  ||C||_F^2 = {np.sum(C**2):.6f}")
    print()

# =============================================================
# Part 5: Key idea - Sum of Cauchy-matrix contributions
# =============================================================
print("\n--- Part 5: Decomposition into Cauchy-matrix terms ---\n")

print("""
For r = p boxplus_n q, by ADMITTED C:
  r' = n * (p^{(1)} boxplus_{n-1} q^{(1)})

The critical points rho_j of r are the roots of p^{(1)} boxplus_{n-1} q^{(1)}.
So rho_j = MSS-convolve(mu'_j-pattern, nu'_j-pattern).

The Cauchy matrix for r is:
  C_r[i,j] = 1 / (lambda_i - rho_j)

And H_r(lambda_i) = (1/2) sum_j C_r[i,j] = (1/2) sum_j 1/(lambda_i - rho_j).

Similarly:
  H_p(mu_i) = (1/2) sum_j 1/(mu_i - mu_j')
  H_q(nu_i) = (1/2) sum_j 1/(nu_i - nu_j')

The roots lambda of r depend on both mu (roots of p) and nu (roots of q).
The critical points rho of r depend on mu' (crits of p) and nu' (crits of q).

This creates a MULTI-LEVEL structure:
  Level 0: roots lambda, mu, nu
  Level 1: critical points rho, mu', nu'
  Level 2: critical points of critical points, etc.
""")

# =============================================================
# Part 6: n=2 exact proof
# =============================================================
print("\n--- Part 6: n=2 exact verification (equality case) ---\n")

# For n=2: p(x) = (x-a)(x-b), q(x) = (x-c)(x-d)
# p boxplus_2 q has coefficients computed by MSS formula
# a_0 = b_0 = 1, a_1 = -(a+b), a_2 = ab, b_1 = -(c+d), b_2 = cd
# c_0 = 1
# c_1 = (2!*1!)/(2!*1!) * a_1 + (1!*2!)/(2!*1!) * b_1 = a_1 + b_1
# c_2 = (0!*2!)/(2!*0!) * a_2*b_0 + (1!*1!)/(2!*0!) * a_1*b_1 + (2!*0!)/(2!*0!) * a_0*b_2
#      = a_2 + (1/2)*a_1*b_1 + b_2

# Wait, let me recompute. c_k = sum_{i+j=k} (n-i)!(n-j)! / (n! (n-k)!) * a_i * b_j
# For n=2:
# c_0: i=0,j=0: 2!*2!/(2!*2!) = 1. c_0 = 1
# c_1: i=0,j=1: 2!*1!/(2!*1!) = 1. i=1,j=0: 1!*2!/(2!*1!) = 1. c_1 = a_1 + b_1 = -(a+b+c+d)
# c_2: i=0,j=2: 2!*0!/(2!*0!) = 1. i=1,j=1: 1!*1!/(2!*0!) = 1/2. i=2,j=0: 0!*2!/(2!*0!) = 1.
# c_2 = ab + (1/2)*(a+b)*(c+d) + cd

print("n=2 MSS convolution:")
print("  p(x) = (x-a)(x-b), q(x) = (x-c)(x-d)")
print("  r(x) = x^2 - (a+b+c+d)x + [ab + cd + (a+b)(c+d)/2]")

# For p: Phi_2 = 2/(b-a)^2, 1/Phi_2 = (b-a)^2/2
# For q: Phi_2 = 2/(d-c)^2, 1/Phi_2 = (d-c)^2/2
# For r: roots are (a+b+c+d +/- sqrt((a+b+c+d)^2 - 4*c_2)) / 2
# gap^2 = (a+b+c+d)^2 - 4*c_2 = (a+b+c+d)^2 - 4ab - 4cd - 2(a+b)(c+d)
#        = (a+b)^2 + 2(a+b)(c+d) + (c+d)^2 - 4ab - 4cd - 2(a+b)(c+d)
#        = (a+b)^2 - 4ab + (c+d)^2 - 4cd
#        = (a-b)^2 + (c-d)^2

print("  gap_r^2 = (a-b)^2 + (c-d)^2")
print("  1/Phi_r = gap_r^2/2 = [(a-b)^2 + (c-d)^2]/2")
print("  1/Phi_p + 1/Phi_q = (a-b)^2/2 + (c-d)^2/2 = [(a-b)^2 + (c-d)^2]/2")
print("  EQUALITY! 1/Phi_r = 1/Phi_p + 1/Phi_q for n=2.")
print("  This is consistent with the claim that equality holds iff n <= 2.")

# Verify
for _ in range(5):
    a, b = np.sort(np.random.randn(2))
    c, d = np.sort(np.random.randn(2))
    p_roots = np.array([a, b])
    q_roots = np.array([c, d])
    r_roots = mss_convolve_roots(p_roots, q_roots)
    lhs = 1/Phi_n(r_roots)
    rhs = 1/Phi_n(p_roots) + 1/Phi_n(q_roots)
    print(f"  1/Phi_r = {lhs:.8f}, 1/Phi_p + 1/Phi_q = {rhs:.8f}, diff = {lhs-rhs:.2e}")

# =============================================================
# Part 7: n=2 identity gives the BLUEPRINT
# =============================================================
print("\n\n--- Part 7: Blueprint from n=2 ---\n")

print("""
For n=2, we proved: gap_r^2 = gap_p^2 + gap_q^2 (PYTHAGOREAN THEOREM!)

This means: 1/Phi_r = gap_r^2/2 = gap_p^2/2 + gap_q^2/2 = 1/Phi_p + 1/Phi_q.

The Pythagorean property gap_r^2 = gap_p^2 + gap_q^2 is because:
  r has 2 roots, gap^2 = (root1 - root2)^2
  = discriminant = (sum of roots)^2 - 4*(product of roots)
  = (mean_p + mean_q)^2 - 4*(var_p + covariance + var_q)  ... no
  Actually = gap_p^2 + gap_q^2 by the explicit formula.

For n >= 3, the conjecture says gap_r^2/2 >= gap_p^2/2 + gap_q^2/2 IN THE
APPROPRIATE GENERALIZED SENSE where "gap^2/2" is replaced by 1/Phi_n.

QUESTION: Is there a "generalized Pythagorean inequality" for root gaps
under MSS convolution?

For n >= 3, Phi_n depends on ALL pairwise gaps 1/(lambda_i - lambda_j),
not just nearest-neighbor gaps. So the relevant quantity is more complex.

KEY FORMULA:
  Phi_n = sum_i (sum_{j!=i} 1/(lambda_i - lambda_j))^2
        = sum_i sum_{j!=i} sum_{k!=i} 1/((lambda_i-lambda_j)(lambda_i-lambda_k))
        = 2 * sum_i sum_{j!=i} 1/(lambda_i-lambda_j)^2
          + 2 * sum_{i} sum_{j!=i} sum_{k!=i,k!=j} 1/((lambda_i-lambda_j)(lambda_i-lambda_k))

Wait, let me expand more carefully:
  Phi_n = sum_i [sum_{j!=i} 1/(lambda_i - lambda_j)]^2
        = sum_i [sum_{j!=i} 1/d_{ij}]^2   where d_{ij} = lambda_i - lambda_j

  = sum_i [sum_{j!=i} 1/d_{ij}^2  + sum_{j!=i} sum_{k!=i,k!=j} 1/(d_{ij}*d_{ik})]

  = sum_i sum_{j!=i} 1/d_{ij}^2 + sum_i sum_{j!=i,k!=i,j!=k} 1/(d_{ij}*d_{ik})

  First sum: S_1 = sum_{i!=j} 1/(lambda_i - lambda_j)^2 = 2 sum_{i<j} 1/(lambda_i-lambda_j)^2
  Second sum: S_2 = cross terms

Actually this is getting complex. Let me just think about the right quantity.

For the ELECTROSTATIC interpretation:
  Phi_n = sum_i E_i^2 where E_i = field at charge i due to all others.

The key question: HOW does MSS convolution act on the field configuration?
""")

# =============================================================
# Part 8: Pairwise distance sums
# =============================================================
print("\n--- Part 8: Pairwise distance sums under MSS ---\n")

# Study various symmetric functions of pairwise distances
np.random.seed(42)
print("Various symmetric functions of pairwise distances d_{ij} = |lambda_i - lambda_j|:\n")

for trial in range(8):
    n = 3
    p_roots = np.sort(np.random.randn(n) * 2)
    q_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(p_roots)) < 0.2:
        p_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(q_roots)) < 0.2:
        q_roots = np.sort(np.random.randn(n) * 2)

    try:
        r_roots = mss_convolve_roots(p_roots, q_roots)
        if np.min(np.diff(r_roots)) < 1e-6:
            continue

        # S2 = sum_{i<j} d_{ij}^2
        def S2(roots):
            return sum((roots[i]-roots[j])**2 for i in range(len(roots)) for j in range(i+1, len(roots)))

        # Sm2 = sum_{i<j} 1/d_{ij}^2
        def Sm2(roots):
            return sum(1/(roots[i]-roots[j])**2 for i in range(len(roots)) for j in range(i+1, len(roots)))

        # Note: Phi_n contains Sm2 as part of its expansion
        # Actually Phi_n = sum_i (sum_{j!=i} 1/d_{ij})^2
        # while Sm2 = sum_{i<j} 1/d_{ij}^2 = (1/2) sum_{i!=j} 1/d_{ij}^2

        print(f"  Trial {trial}:")
        print(f"    S2: p={S2(p_roots):.4f}, q={S2(q_roots):.4f}, r={S2(r_roots):.4f}")
        print(f"    S2_r vs S2_p + S2_q: {S2(r_roots):.4f} vs {S2(p_roots)+S2(q_roots):.4f}")
        print(f"    Sm2: p={Sm2(p_roots):.4f}, q={Sm2(q_roots):.4f}, r={Sm2(r_roots):.4f}")
        print(f"    Phi: p={Phi_n(p_roots):.4f}, q={Phi_n(q_roots):.4f}, r={Phi_n(r_roots):.4f}")
        print()
    except:
        continue

# =============================================================
# Part 9: The n=2 Pythagorean identity in terms of moments
# =============================================================
print("\n--- Part 9: Connection to power sums ---\n")

print("""
For n=2 with roots a, b:
  S2 = (a-b)^2 = e_1^2 - 4e_2 = (a+b)^2 - 4ab  (where e_k are elem. sym.)
  S2 is a function of the first two power sums:
  S2 = 2(p_1^2/4 - e_2) = 2*sigma_2 (the variance * 2)

For MSS: the cumulants of r are the sums of cumulants of p and q.
  k_2(r) = k_2(p) + k_2(q)

For n=2: gap^2 = 2*k_2, so gap_r^2 = 2*(k_2(p)+k_2(q)) = gap_p^2 + gap_q^2.
This is WHY the n=2 case gives equality!

For n >= 3: Phi_n depends on ALL cumulants k_2, k_3, ..., k_n.
But the cumulants of r are still k_j(r) = k_j(p) + k_j(q).

INSIGHT: 1/Phi_n is a function of (k_2, ..., k_n) and we need superadditivity
under coordinate-wise addition of cumulant vectors.

This IS the cumulant approach that was declared failed. But maybe the
ROOT GEOMETRY gives insight into WHY the cumulant approach should work.

The point is: Phi_n as a function of roots is conceptually cleaner
than as a function of cumulants, because in roots the INTERLACING structure
of MSS convolution is visible.
""")

# =============================================================
# Part 10: Cauchy-Schwarz via n=3 explicit formula
# =============================================================
print("\n--- Part 10: n=3 Cauchy-Schwarz attempt ---\n")

# For n=3, Phi_3 = 2*(s^4 + 2s^3t + 3s^2t^2 + 2st^3 + t^4) / (s^2 t^2 (s+t)^2)
# where s = d1, t = d2 are the gaps.
# 1/Phi_3 = s^2 t^2 (s+t)^2 / (2*(s^4 + 2s^3t + 3s^2t^2 + 2st^3 + t^4))

# Note: s^4 + 2s^3t + 3s^2t^2 + 2st^3 + t^4
# = s^4 + t^4 + 2s^3t + 2st^3 + 3s^2t^2
# = (s^2+t^2)^2 + 2st(s^2+t^2) + s^2t^2
# = (s^2 + st + t^2)^2

s, t = sp.symbols('s t', positive=True)
num = s**4 + 2*s**3*t + 3*s**2*t**2 + 2*s*t**3 + t**4
factored = sp.factor(num)
print(f"  s^4 + 2s^3t + 3s^2t^2 + 2st^3 + t^4 = {factored}")

# So Phi_3 = 2(s^2+st+t^2)^2 / (s^2 t^2 (s+t)^2)
# 1/Phi_3 = s^2 t^2 (s+t)^2 / (2(s^2+st+t^2)^2)

print(f"\n  Phi_3 = 2(s^2+st+t^2)^2 / (s^2 t^2 (s+t)^2)")
print(f"  1/Phi_3 = s^2 t^2 (s+t)^2 / (2(s^2+st+t^2)^2)")

# Now for the MSS convolution of n=3 polynomials:
# p has gaps s_p, t_p; q has gaps s_q, t_q
# r = p boxplus_3 q has gaps s_r, t_r

# We need: 1/Phi_3(r) >= 1/Phi_3(p) + 1/Phi_3(q)
# i.e., s_r^2 t_r^2 (s_r+t_r)^2 / (2(s_r^2+s_rt_r+t_r^2)^2)
#     >= s_p^2 t_p^2 (s_p+t_p)^2 / (2(s_p^2+s_pt_p+t_p^2)^2)
#      + s_q^2 t_q^2 (s_q+t_q)^2 / (2(s_q^2+s_qt_q+t_q^2)^2)

# This is a relation between gaps of p, q, and r.
# The gaps of r depend on the FULL root configurations of p and q,
# not just their gaps (because MSS convolution depends on ALL coefficients).

# Let me study the relationship between gaps of r and those of p, q numerically.
print("\n  Studying gap transformation under MSS for n=3:")
np.random.seed(12345)
for trial in range(8):
    n = 3
    p_roots = np.sort(np.random.randn(n) * 2)
    q_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(p_roots)) < 0.2:
        p_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(q_roots)) < 0.2:
        q_roots = np.sort(np.random.randn(n) * 2)

    try:
        r_roots = mss_convolve_roots(p_roots, q_roots)
        sp_val, tp_val = np.diff(p_roots)
        sq_val, tq_val = np.diff(q_roots)
        sr_val, tr_val = np.diff(r_roots)

        # Spread = total gap = s+t
        Sp = sp_val + tp_val
        Sq = sq_val + tq_val
        Sr = sr_val + tr_val

        # n=2 analogy: gap_r^2 = gap_p^2 + gap_q^2
        # For n=3: is spread_r^2 >= spread_p^2 + spread_q^2?
        # spread = s + t = lambda_n - lambda_1

        print(f"  Trial {trial}: Sp={Sp:.3f}, Sq={Sq:.3f}, Sr={Sr:.3f}")
        print(f"    Sp^2+Sq^2 = {Sp**2+Sq**2:.3f}, Sr^2 = {Sr**2:.3f}")
        print(f"    Sr^2 >= Sp^2+Sq^2? {Sr**2 >= Sp**2+Sq**2}")
    except:
        continue

print("\n  NOTE: spread_r^2 >= spread_p^2 + spread_q^2 does NOT always hold.")
print("  So the n=2 Pythagorean identity doesn't extend to spreads.")

print("\nDone with Part 3.")
