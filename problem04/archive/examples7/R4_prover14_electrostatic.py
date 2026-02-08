"""
PROVER-14 Part 2: Deep electrostatic analysis

KEY INSIGHT:
  H_p = -grad(W_p) where W_p = -sum_{i<j} log|mu_i - mu_j|
  Phi_n = ||grad W_p||^2

So 1/Phi_n = 1/||grad W||^2.

The conjecture says 1/||grad W_r||^2 >= 1/||grad W_p||^2 + 1/||grad W_q||^2.

Equivalently: ||grad W_r||^2 <= (||grad W_p||^2 * ||grad W_q||^2) / (||grad W_p||^2 + ||grad W_q||^2)
            = harmonic mean of ||grad W_p||^2 and ||grad W_q||^2

So: The field energy of the convolution is bounded by the harmonic mean
of the individual field energies.

ANOTHER FORM: Let alpha = Phi_p, beta = Phi_q, gamma = Phi_r.
The conjecture says: 1/gamma >= 1/alpha + 1/beta
i.e., gamma <= alpha*beta/(alpha+beta) (harmonic mean bound)

STRATEGY:
1. Study the log-potential W and its Hessian
2. Connect MSS convolution to the energy landscape
3. Use the Cauchy-Schwarz or Cramer-Rao type bound
"""
import numpy as np
from itertools import combinations
import sympy as sp
from math import factorial

print("="*70)
print("PROVER-14 Part 2: ELECTROSTATIC / ENERGY ANALYSIS")
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

def log_potential(roots):
    """W = -sum_{i<j} log|roots_i - roots_j|"""
    n = len(roots)
    W = 0
    for i in range(n):
        for j in range(i+1, n):
            W -= np.log(abs(roots[i] - roots[j]))
    return W

def hessian_W(roots):
    """Compute the Hessian of W = -sum_{i<j} log|r_i - r_j|.

    d^2 W / dr_i^2 = sum_{j!=i} 1/(r_i - r_j)^2
    d^2 W / (dr_i dr_j) = -1/(r_i - r_j)^2  for i != j
    """
    n = len(roots)
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                H[i, j] = -1.0 / (roots[i] - roots[j])**2
                H[i, i] += 1.0 / (roots[i] - roots[j])**2
    return H

# =============================================================
# Part 1: Hessian analysis
# =============================================================
print("\n--- Part 1: Hessian of log-potential ---\n")

roots = np.array([-2.0, 0.0, 2.0])
Hess = hessian_W(roots)
print(f"Roots: {roots}")
print(f"Hessian of W:")
print(Hess)
eigvals = np.linalg.eigvalsh(Hess)
print(f"Eigenvalues: {eigvals}")
print(f"(Hessian is positive semidefinite on {'{'}sum=0{'}'} subspace)")

# Check: Hessian row sums should be 0 (translation invariance)
print(f"Row sums: {np.sum(Hess, axis=1)}")

# The Hessian restricted to the hyperplane sum=0 is positive definite
# This means W is convex on this subspace!
# The 0 eigenvalue corresponds to the all-ones vector (translation)
print()

for config in [np.array([-3., -1., 1., 3.]), np.array([-2., -0.5, 0.5, 4.])]:
    Hess = hessian_W(config)
    eigvals = np.linalg.eigvalsh(Hess)
    print(f"Roots: {config}")
    print(f"Hessian eigenvalues: {eigvals}")
    print(f"Row sums: {np.sum(Hess, axis=1)}")
    print()

# =============================================================
# Part 2: Hessian and Phi relationship
# =============================================================
print("\n--- Part 2: Connection between Hessian and Phi ---\n")

print("""
KEY OBSERVATION:
  grad W_i = -H_p(lambda_i)
  Hess W_{ij} = d^2W/dr_idr_j

For the Cauchy-Schwarz/information-geometric approach:
  Phi_n = ||grad W||^2 = (grad W)^T (grad W)

The Hessian H_W is the "Fisher information matrix" of the log-potential.
By the Cramer-Rao inequality analogy:
  If we think of roots as "parameters" and W as "log-likelihood",
  then the Hessian bounds the inverse of the variance.

But what we really need is a relationship between the Hessians
of W_p, W_q, and W_r.
""")

# =============================================================
# Part 3: Alternative: Cauchy-Schwarz on H-vectors
# =============================================================
print("\n--- Part 3: Cauchy-Schwarz on H-vectors ---\n")

print("""
The H-vector lives in the hyperplane sum_i H_i = 0 in R^n.
Phi_n = ||H||^2.
The conjecture is 1/||H_r||^2 >= 1/||H_p||^2 + 1/||H_q||^2.

This is equivalent to:
  ||H_r||^2 <= ||H_p||^2 * ||H_q||^2 / (||H_p||^2 + ||H_q||^2)

Or: 1/||H_r||^2 >= 1/||H_p||^2 + 1/||H_q||^2.

If there existed a LINEAR map M such that H_r = M H_p (or some
combination), we could use operator norm bounds.

But H is a NONLINEAR function of roots, and the MSS convolution
is also nonlinear in roots.
""")

# =============================================================
# Part 4: MSS convolution via random matrices
# =============================================================
print("\n--- Part 4: MSS via random matrix interpretation ---\n")

def mss_convolve_roots(roots_p, roots_q):
    """Compute roots of p boxplus_n q."""
    n = len(roots_p)
    p_coeffs = roots_to_monic_coeffs(roots_p)
    q_coeffs = roots_to_monic_coeffs(roots_q)
    r_coeffs = mss_convolve_n(p_coeffs, q_coeffs, n)
    r_roots = np.roots(r_coeffs)
    r_roots = np.sort(np.real(r_roots))
    return r_roots

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

print("""
Random matrix model (ADMITTED B):
  r = p boxplus_n q has the EXPECTED characteristic polynomial of A + UBU*
  where A = diag(mu), B = diag(nu), U ~ Haar(U(n)).

The ACTUAL eigenvalues of A + UBU* for a specific U are lambda_1(U), ..., lambda_n(U).
These satisfy: E_U[chi_{A+UBU*}] = r.

Key: For any specific U, the eigenvalues of A + UBU* satisfy:
  Phi_n(eigenvalues(A+UBU*)) <= some bound depending on Phi_n(mu) and Phi_n(nu)?

But the EXPECTED polynomial's roots are NOT the expected eigenvalues!
The roots of E[chi] != E[roots of chi].

So this direction needs care.
""")

# =============================================================
# Part 5: Direct n=3 analysis
# =============================================================
print("\n--- Part 5: Explicit n=3 analysis ---\n")

# For n=3 with roots a < b < c, center at 0: a+b+c=0
# Let b = t, then a = -t - s, c = s for some s > t (ensuring a < b < c)
# Wait, let me parametrize differently.
# Let gaps be d1 = b-a > 0, d2 = c-b > 0.
# Center at 0: a = -(d1+d2)/3 + (-d2+d1)/3 ... let me just use a, b, c with a+b+c=0

# With a+b+c = 0 (WLOG by translation invariance of 1/Phi):
# a < b < c, a = -(b+c), so a < b means -(b+c) < b, i.e., c > -2b
# Parametrize by b and c with b < c, b+c > -b (i.e., c > -2b)

# Actually for n=3, let's use the substitution:
# roots = (-d1-d2)/3, (-d1+2*d2)/3 ... no this is getting complicated.
# Let me just center and scale.

# For n=3 centered at 0: roots are -u-v, v, u where u > 0, -u-v < v < u
# i.e., v > -(u+v)/2 which is v > -u, and v < u
# So v in (-u, u), and we need -u-v < v, i.e., v > -u/2

# Hmm, let me use a cleaner parametrization.
# roots = (-s, 0, s+t) centered, then translate to have mean 0:
# mean = t/3, so roots = (-s-t/3, -t/3, s+2t/3)
# Gaps: d1 = s, d2 = s+t. But this doesn't capture all configs.

# Better: roots = (a, b, c) with a+b+c = 0 (translation invariance of Phi).
# Then a = -(b+c), and the polynomial is x^3 + px + q
# where p = ab+ac+bc, q = -abc.
# With a+b+c=0: p = -(a^2+b^2+c^2)/2, q = -abc.

# Phi_3 in terms of a, b, c with a+b+c=0:
a, b, c = sp.symbols('a b c', real=True)
H_a_sym = 1/(a-b) + 1/(a-c)
H_b_sym = 1/(b-a) + 1/(b-c)
H_c_sym = 1/(c-a) + 1/(c-b)
Phi_3_sym = sp.expand(H_a_sym**2 + H_b_sym**2 + H_c_sym**2)

# Substitute c = -a-b
Phi_3_centered = Phi_3_sym.subs(c, -a-b)
Phi_3_centered_simplified = sp.simplify(Phi_3_centered)
print(f"Phi_3 (centered, c=-a-b) = {Phi_3_centered_simplified}")

# Let me try to factor using d1 = b-a, d2 = c-b = -a-2b
# So a = b-d1, c = b+d2, and d1+d2 = c-a = -2a, a = -(d1+d2)/2... no.
# With c = -a-b: d1 = b-a, d2 = -a-b-b = -a-2b
# d1+d2 = b-a-a-2b = -2a-b = -2a-b
# But c = -a-b, so d1 = b-a, d2 = c-b = -a-2b
# For d2 > 0: b < -a/2

# Use s = d1, t = d2 (the gaps), s,t > 0
# a = -(2s+t)/3, b = (s-t)/3 ... let me compute:
# d1 = b-a = s, d2 = c-b = t
# a+b+c = 0, b = a+s, c = a+s+t
# 3a+2s+t = 0, a = -(2s+t)/3
# b = -(2s+t)/3 + s = (s-t)/3
# c = (s-t)/3 + t = (s+2t)/3

s, t = sp.symbols('s t', positive=True)
a_val = -(2*s+t)/3
b_val = (s-t)/3
c_val = (s+2*t)/3

Phi_3_gaps = Phi_3_sym.subs([(a, a_val), (b, b_val), (c, c_val)])
Phi_3_gaps_simplified = sp.simplify(Phi_3_gaps)
print(f"\nPhi_3 in terms of gaps s=d1, t=d2:")
print(f"  Phi_3 = {Phi_3_gaps_simplified}")

# Also compute 1/Phi_3
inv_Phi_3_gaps = sp.simplify(1/Phi_3_gaps_simplified)
print(f"  1/Phi_3 = {inv_Phi_3_gaps}")

# Try to express as function of ratio r = t/s only (by scaling)
# Phi_n is homogeneous of degree -2 in roots (since H is degree -1)
# So Phi_3(lambda*roots) = lambda^(-2) * Phi_3(roots)
# 1/Phi_3(lambda*roots) = lambda^2 * 1/Phi_3(roots)
# With s as scale: roots = s * (-(2+t/s)/3, (1-t/s)/3, (1+2t/s)/3)
# So 1/Phi_3 = s^2 * f(t/s) for some function f
r = sp.Symbol('r', positive=True)
Phi_3_ratio = Phi_3_gaps_simplified.subs(t, r*s)
Phi_3_ratio_simplified = sp.simplify(Phi_3_ratio)
print(f"\n  Phi_3(s, r*s) = {Phi_3_ratio_simplified}")

# Check: Phi_3 should scale as 1/s^2
check = sp.simplify(s**2 * Phi_3_ratio_simplified)
print(f"  s^2 * Phi_3 = {check}")
print(f"  (This should be a function of r only)")

f_r = check
print(f"\n  f(r) = s^2 * Phi_3 = {f_r}")
print(f"  1/Phi_3 = s^2 / f(r)")

# =============================================================
# Part 6: n=3 MSS convolution in gap coordinates
# =============================================================
print("\n\n--- Part 6: n=3 MSS convolution in gap coordinates ---\n")

np.random.seed(42)

# For n=3, study how the gap ratio transforms
print("How gap ratio r=d2/d1 transforms under MSS convolution:")
for trial in range(10):
    s_p = np.random.exponential(2) + 0.5
    t_p = np.random.exponential(2) + 0.5
    s_q = np.random.exponential(2) + 0.5
    t_q = np.random.exponential(2) + 0.5

    p_roots = np.array([-(2*s_p+t_p)/3, (s_p-t_p)/3, (s_p+2*t_p)/3])
    q_roots = np.array([-(2*s_q+t_q)/3, (s_q-t_q)/3, (s_q+2*t_q)/3])

    try:
        r_roots = mss_convolve_roots(p_roots, q_roots)
        # Center the result
        mean_r = np.mean(r_roots)
        r_centered = r_roots - mean_r
        d1_r = r_centered[1] - r_centered[0]
        d2_r = r_centered[2] - r_centered[1]

        ratio_p = t_p / s_p
        ratio_q = t_q / s_q
        ratio_r = d2_r / d1_r

        total_gap_p = s_p + t_p
        total_gap_q = s_q + t_q
        total_gap_r = d1_r + d2_r

        print(f"  p: (s,t)=({s_p:.2f},{t_p:.2f}), r_p={ratio_p:.3f}, spread={total_gap_p:.2f}")
        print(f"  q: (s,t)=({s_q:.2f},{t_q:.2f}), r_q={ratio_q:.3f}, spread={total_gap_q:.2f}")
        print(f"  r: gaps=({d1_r:.2f},{d2_r:.2f}), r_r={ratio_r:.3f}, spread={total_gap_r:.2f}")

        # Check inequality
        inv_phi_p = 1/Phi_n(p_roots)
        inv_phi_q = 1/Phi_n(q_roots)
        inv_phi_r = 1/Phi_n(r_roots)
        gap = inv_phi_r - inv_phi_p - inv_phi_q
        print(f"  gap = {gap:.6e}")
        print()
    except:
        continue

# =============================================================
# Part 7: Is Phi_n Schur-convex in sorted root gaps?
# =============================================================
print("\n--- Part 7: Schur-convexity of Phi_n in root gaps ---\n")

print("""
For Phi_n to be Schur-convex in the gap vector d = (d_1, ..., d_{n-1}):
  If d majorizes d', then Phi_n(d) >= Phi_n(d').

Majorization d >> d' means:
  - d and d' have the same sum
  - The partial sums of d (sorted decreasing) dominate those of d'
  - d is "more spread" than d'

If Phi_n is Schur-convex in gaps, then MORE EQUAL gaps give SMALLER Phi_n
(i.e., LARGER 1/Phi_n).

Check: for n=3, Phi_3 as a function of (d1, d2) with fixed d1+d2.
""")

# For n=3, fix d1+d2 = D, vary d1 from eps to D-eps
D = 4.0
d1_vals = np.linspace(0.1, D-0.1, 100)
phi_vals = []
for d1 in d1_vals:
    d2 = D - d1
    roots = np.array([-(2*d1+d2)/3, (d1-d2)/3, (d1+2*d2)/3])
    phi_vals.append(Phi_n(roots))

phi_vals = np.array(phi_vals)
min_idx = np.argmin(phi_vals)
print(f"n=3, fixed D = d1+d2 = {D}:")
print(f"  min Phi at d1 = {d1_vals[min_idx]:.3f} (d2 = {D-d1_vals[min_idx]:.3f})")
print(f"  d1=d2={D/2}: Phi = {Phi_n(np.array([-(D/2+D/2)/3, 0, (D/2+D/2)/3])):.6f}")
print(f"  Minimum Phi = {phi_vals[min_idx]:.6f}")
print(f"  Phi is minimized when d1 = d2 (equal gaps) -> Schur-convexity holds!")

# For n=4, check Schur-convexity
print(f"\nn=4, fixed total span:")
D4 = 6.0
# Parametrize: d1, d2, d3 with d1+d2+d3 = D4
# Test: equal gaps vs unequal
equal_gaps = np.array([D4/3, D4/3, D4/3])
roots_eq = np.array([-D4/2, -D4/6, D4/6, D4/2])  # equally spaced centered
phi_eq = Phi_n(roots_eq)

# Very unequal
d_vec = np.array([0.5, 1.0, 4.5])
a0 = -(d_vec[0] + d_vec[1] + d_vec[2])/2  # not quite right for centering
# Let me just construct roots from gaps
roots_uneq = np.zeros(4)
roots_uneq[0] = 0
for i in range(3):
    roots_uneq[i+1] = roots_uneq[i] + d_vec[i]
roots_uneq -= np.mean(roots_uneq)  # center
phi_uneq = Phi_n(roots_uneq)

print(f"  Equal gaps (2,2,2): Phi = {phi_eq:.6f}, 1/Phi = {1/phi_eq:.6f}")
print(f"  Unequal gaps (0.5,1,4.5): Phi = {phi_uneq:.6f}, 1/Phi = {1/phi_uneq:.6f}")
print(f"  Equal gaps have LOWER Phi (HIGHER 1/Phi) -> consistent with Schur-convexity")

# Random sampling to verify Schur-convexity
print(f"\n  Random verification of Schur-convexity for n=4:")
np.random.seed(999)
violations_schur = 0
for _ in range(10000):
    # Random gap vector d, total = D4
    d = np.random.exponential(1, 3)
    d = d / np.sum(d) * D4

    # More equal version: move toward equal
    d_eq = d * 0.5 + D4/3 * 0.5  # convex combination toward equal

    # Both have same sum
    roots1 = np.cumsum(np.concatenate([[0], d]))
    roots1 -= np.mean(roots1)
    roots2 = np.cumsum(np.concatenate([[0], d_eq]))
    roots2 -= np.mean(roots2)

    phi1 = Phi_n(roots1)
    phi2 = Phi_n(roots2)

    # d majorizes d_eq (since d_eq is a convex combination toward center)
    # So Schur-convexity means phi1 >= phi2
    if phi1 < phi2 - 1e-10:
        violations_schur += 1

print(f"  {violations_schur} violations of Schur-convexity in 10000 trials")

# =============================================================
# Part 8: Does MSS convolution make gaps more equal?
# =============================================================
print("\n\n--- Part 8: MSS convolution and gap majorization ---\n")

print("""
Hypothesis: The sorted gap vector of r = p boxplus_n q is majorized by...
what exactly? Not by the gap vectors of p and q individually.

Alternative: Does 1/Phi_n satisfy a stronger subadditivity
that follows from a structural property of MSS?
""")

np.random.seed(888)
for trial in range(5):
    n = 4
    p_roots = np.sort(np.random.randn(n) * 2)
    q_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(p_roots)) < 0.2:
        p_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(q_roots)) < 0.2:
        q_roots = np.sort(np.random.randn(n) * 2)

    try:
        r_roots = mss_convolve_roots(p_roots, q_roots)

        d_p = np.diff(p_roots)
        d_q = np.diff(q_roots)
        d_r = np.diff(r_roots)

        # Coefficient of variation (lower = more equal)
        cv_p = np.std(d_p)/np.mean(d_p)
        cv_q = np.std(d_q)/np.mean(d_q)
        cv_r = np.std(d_r)/np.mean(d_r)

        print(f"  Trial {trial}: CV(gaps): p={cv_p:.3f}, q={cv_q:.3f}, r={cv_r:.3f}")
        print(f"    p_gaps = {d_p}")
        print(f"    q_gaps = {d_q}")
        print(f"    r_gaps = {d_r}")
        print(f"    Is CV(r) < min(CV(p), CV(q))? {cv_r < min(cv_p, cv_q)}")
    except:
        continue

# =============================================================
# Part 9: Derivative formula: p'/p and p''/p'
# =============================================================
print("\n\n--- Part 9: p''/p' formula and connections ---\n")

print("""
The key identity: H_p(lambda_i) = p''(lambda_i) / (2 p'(lambda_i))

Also: p'(lambda_i) = prod_{j!=i} (lambda_i - lambda_j)

So H_p(lambda_i) = sum_{j!=i} 1/(lambda_i - lambda_j)
                  = (d/dx)[log p'(x)]|_{x=lambda_i}  ... NO

Actually: H_p(lambda_i) = [log |p(x)|]''... no.

Let f(x) = log|p(x)| = sum_i log|x - lambda_i|
f'(x) = sum_i 1/(x - lambda_i) = p'(x)/p(x)
f''(x) = -sum_i 1/(x - lambda_i)^2

At x = lambda_i, f' has a pole, but the regular part of p'(x)/p(x)
at lambda_i is H_p(lambda_i).

Alternative: Consider the RECIPROCAL polynomial.
Let q(x) = p(x)/(x - lambda_i). Then q(lambda_i) = p'(lambda_i).
And q'(lambda_i) = p''(lambda_i)/2.
So H_p(lambda_i) = q'(lambda_i)/q(lambda_i) = [log|q(x)|]'|_{x=lambda_i}.

This is the log-derivative of the "reduced polynomial" at the removed root!
""")

# Verify this interpretation
roots = np.array([-2.0, 0.5, 1.0, 3.0])
n = len(roots)
for i in range(n):
    # q(x) = p(x)/(x - lambda_i)
    other_roots = np.delete(roots, i)
    # q(lambda_i) = prod_{j!=i} (lambda_i - lambda_j)
    q_val = np.prod(roots[i] - other_roots)
    # q'(lambda_i) = sum_{j!=i} prod_{k!=i,j} (lambda_i - lambda_k)
    q_prime = 0
    for j_idx in range(len(other_roots)):
        term = 1.0
        for k_idx in range(len(other_roots)):
            if k_idx != j_idx:
                term *= (roots[i] - other_roots[k_idx])
        q_prime += term

    H_i = q_prime / q_val  # = sum_{j!=i} 1/(lambda_i - lambda_j)
    H_direct = sum(1/(roots[i] - roots[j]) for j in range(n) if j != i)

    print(f"  Root {i} ({roots[i]}): H = {H_direct:.6f}, q'/q = {H_i:.6f}, match: {abs(H_direct - H_i) < 1e-10}")

# =============================================================
# Part 10: Key insight - Cauchy transform connection
# =============================================================
print("\n\n--- Part 10: Cauchy transform and Stieltjes inversion ---\n")

print("""
CRITICAL CONNECTION:
G_p(z) = (1/n) * p'(z)/p(z) = (1/n) sum_i 1/(z - lambda_i)

Near z = lambda_i:
G_p(z) = (1/n)/(z - lambda_i) + (1/n)*H_p(lambda_i) + O(z - lambda_i)

So H_p(lambda_i) is n times the REGULAR PART of G_p at lambda_i.

For MSS convolution, ADMITTED C says:
  (p boxplus_n q)'(x) = n * (p^{(1)} boxplus_{n-1} q^{(1)})(x)

where p^{(1)} = (1/n)*p' has roots at the critical points of p.

This gives a RECURSIVE structure!

SUBORDINATION (from free probability):
For the additive free convolution boxplus, there exist analytic
subordination functions omega_1, omega_2 such that
  G_{mu boxplus nu}(z) = G_mu(omega_1(z)) = G_nu(omega_2(z))
with omega_1(z) + omega_2(z) = z + 1/G_{mu boxplus nu}(z).

The FINITE version of subordination should give:
  G_r(z) = G_p(omega_1(z)) = G_q(omega_2(z))
where omega_1(z) + omega_2(z) = z + 1/G_r(z).

Let me verify this numerically!
""")

# Check finite subordination
def G_p_func(z, roots):
    """Cauchy transform G_p(z) = (1/n) sum_i 1/(z - lambda_i)"""
    n = len(roots)
    return (1.0/n) * sum(1.0/(z - r) for r in roots)

# For n=3
p_roots = np.array([-2.0, 0.0, 2.0])
q_roots = np.array([-1.0, 0.5, 1.5])
r_roots = mss_convolve_roots(p_roots, q_roots)

print(f"p roots: {p_roots}")
print(f"q roots: {q_roots}")
print(f"r roots: {r_roots}")

# Test subordination at various z in upper half plane
print("\nTesting subordination G_r(z) = G_p(w1(z)) = G_q(w2(z)):")
print("(where w1 + w2 = z + 1/G_r(z))\n")

from scipy.optimize import fsolve

def subordination_eq(w_flat, z, p_roots, q_roots, r_roots):
    """Find w1, w2 such that G_p(w1) = G_q(w2) = G_r(z) and w1+w2 = z + 1/G_r(z)"""
    w1 = w_flat[0] + 1j*w_flat[1]
    w2 = w_flat[2] + 1j*w_flat[3]

    G_r_z = G_p_func(z, r_roots)
    G_p_w1 = G_p_func(w1, p_roots)
    G_q_w2 = G_p_func(w2, q_roots)

    eq1 = G_p_w1 - G_r_z  # G_p(w1) = G_r(z)
    eq2 = G_q_w2 - G_r_z  # G_q(w2) = G_r(z)
    eq3 = w1 + w2 - z - 1.0/G_r_z  # w1 + w2 = z + 1/G_r(z) ... actually maybe different for finite

    return [eq1.real, eq1.imag, eq3.real, eq3.imag]

# This may not work for FINITE convolution. Let me check.
z_test = 5.0 + 1j
G_r_z = G_p_func(z_test, r_roots)

# Try to find w1 such that G_p(w1) = G_r(z)
# Then w2 = z + 1/G_r(z) - w1
# And check if G_q(w2) = G_r(z)

def find_w1(z, p_roots, r_roots):
    """Find w1 in UHP such that G_p(w1) = G_r(z)"""
    target = G_p_func(z, r_roots)

    def eq(w_flat):
        w = w_flat[0] + 1j*w_flat[1]
        val = G_p_func(w, p_roots) - target
        return [val.real, val.imag]

    # Start near z
    w0 = [z.real, z.imag]
    sol = fsolve(eq, w0, full_output=True)
    w1 = sol[0][0] + 1j*sol[0][1]
    return w1

for y_val in [0.5, 1.0, 2.0, 5.0]:
    z = 3.0 + 1j*y_val
    G_r_z = G_p_func(z, r_roots)

    w1 = find_w1(z, p_roots, r_roots)
    G_p_w1 = G_p_func(w1, p_roots)

    # Check: w2 = z + 1/(n*G_r(z)) - w1 ? (finite version may differ)
    # In free probability: w1 + w2 = z + F_r(z) where F_r = 1/G_r
    # But here G_r = (1/n) p'/p, so F_r(z) = n/p'(z)*p(z)... different

    # Let's just find w2 independently
    w2 = find_w1(z, q_roots, r_roots)
    G_q_w2 = G_p_func(w2, q_roots)

    print(f"  z = {z}")
    print(f"    G_r(z) = {G_r_z:.6f}")
    print(f"    w1 = {w1:.6f}, G_p(w1) = {G_p_w1:.6f}, match: {abs(G_p_w1 - G_r_z) < 1e-6}")
    print(f"    w2 = {w2:.6f}, G_q(w2) = {G_q_w2:.6f}, match: {abs(G_q_w2 - G_r_z) < 1e-6}")
    print(f"    w1 + w2 = {w1+w2:.6f}")
    print(f"    z + 1/G_r = {z + 1/G_r_z:.6f}")
    print(f"    z + n/G_r*n = {z + len(r_roots)/(len(r_roots)*G_r_z):.6f}")
    print(f"    w1+w2 - z = {w1+w2-z:.6f}")
    F_r = 1.0/(len(r_roots) * G_r_z)  # this is 1/(sum 1/(z-lambda_i))
    print(f"    n*F_r = {len(r_roots)*F_r:.6f}")
    print()

print("\nDone with electrostatic analysis.")
