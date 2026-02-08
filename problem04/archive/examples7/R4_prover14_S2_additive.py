"""
PROVER-14 Part 4: Exploiting S2 additivity

CRITICAL DISCOVERY: S2 = sum_{i<j} (lambda_i - lambda_j)^2 is EXACTLY ADDITIVE
under MSS convolution: S2(r) = S2(p) + S2(q).

This is because S2 = n * k_2 (where k_2 is the second finite free cumulant),
and cumulants are additive under MSS.

Actually, S2 = n * (second centered moment) = n * k_2. Let me verify.

PLAN:
1. Verify the S2 = n*k2 relationship and understand precisely
2. Express Phi_n in terms of S2 and HIGHER order invariants
3. Use the fact that S2 is additive to reduce the problem
4. The conjecture then becomes: for fixed S2, how does Phi_n depend on
   the higher order structure?

The key insight may be:
  1/Phi_n is a function F(S2, S3, S4, ...) where S_k are generalized
  symmetric functions of roots. S2 is additive, and the inequality
  should follow from convexity/concavity of F in the other arguments.
"""
import numpy as np
from itertools import combinations
from math import factorial
import sympy as sp

print("="*70)
print("PROVER-14 Part 4: S2 ADDITIVITY AND CONSEQUENCES")
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

def S2(roots):
    """Sum of squared pairwise differences."""
    return sum((roots[i]-roots[j])**2 for i in range(len(roots)) for j in range(i+1, len(roots)))

def Sm2(roots):
    """Sum of reciprocal squared pairwise differences."""
    return sum(1/(roots[i]-roots[j])**2 for i in range(len(roots)) for j in range(i+1, len(roots)))

# =============================================================
# Part 1: Verify S2 = n * k2
# =============================================================
print("\n--- Part 1: S2 and cumulants ---\n")

# S2 = sum_{i<j} (lambda_i - lambda_j)^2
# = (1/2) sum_{i,j} (lambda_i - lambda_j)^2
# = (1/2) * [2n sum_i lambda_i^2 - 2 (sum_i lambda_i)^2]  -- by expansion
# = n * sum_i lambda_i^2 - (sum_i lambda_i)^2
# = n * [sum lambda_i^2 - (sum lambda_i)^2/n]
# = n * (p_2 - p_1^2/n)  where p_k = sum lambda_i^k

# The second cumulant k2 satisfies:
# For the empirical distribution (1/n) sum delta_{lambda_i}:
# variance = (1/n) sum lambda_i^2 - ((1/n) sum lambda_i)^2
# = p_2/n - p_1^2/n^2

# So S2 = n * (p_2 - p_1^2/n) = n * p_2 - p_1^2 = n * n * variance = n^2 * var
# Actually: S2 = n * sum lambda_i^2 - (sum lambda_i)^2

# The FINITE FREE cumulant k_2 (Marcus-Spielman-Srivastava) for a monic polynomial
# of degree n is defined differently from the classical variance.
# k_2 = -a_2/(n choose 2) + a_1^2/(n*(n-1))  ... need to check.

# For MSS: if p(x) = x^n + a_1 x^{n-1} + a_2 x^{n-2} + ..., then
# the finite free cumulants are defined so that kappa_j(p boxplus q) = kappa_j(p) + kappa_j(q).
# kappa_1 = -a_1/n = mean of roots
# kappa_2 relates to a_2 via e_2 = a_2 = sum_{i<j} lambda_i lambda_j

# Actually: e_1 = -a_1, e_2 = a_2.
# S2 = n*p_2 - p_1^2 = n*(e_1^2 - 2e_2) - e_1^2 = (n-1)*e_1^2 - 2n*e_2
# where p_1 = sum lambda_i = e_1, p_2 = sum lambda_i^2 = e_1^2 - 2e_2.

# Check additivity numerically
np.random.seed(42)
print("Verifying S2 additivity under MSS:")
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
        s2p = S2(p_roots)
        s2q = S2(q_roots)
        s2r = S2(r_roots)
        print(f"  Trial {trial}: S2(p)={s2p:.6f}, S2(q)={s2q:.6f}, S2(r)={s2r:.6f}, S2(p)+S2(q)={s2p+s2q:.6f}, diff={s2r-s2p-s2q:.2e}")
    except:
        continue

# Also check S3 = sum_{i<j} (lambda_i - lambda_j)^3 ... not symmetric
# Instead check S_k = sum_{i<j} |lambda_i - lambda_j|^k for k=4
print("\nChecking S4 = sum_{i<j} (lambda_i - lambda_j)^4:")
np.random.seed(42)
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
        def S4(roots):
            return sum((roots[i]-roots[j])**4 for i in range(len(roots)) for j in range(i+1, len(roots)))
        s4p = S4(p_roots)
        s4q = S4(q_roots)
        s4r = S4(r_roots)
        print(f"  Trial {trial}: S4(p)={s4p:.4f}, S4(q)={s4q:.4f}, S4(r)={s4r:.4f}, S4(p)+S4(q)={s4p+s4q:.4f}, ratio={s4r/(s4p+s4q):.4f}")
    except:
        continue

# =============================================================
# Part 2: Express Phi_n in terms of power sums of differences
# =============================================================
print("\n\n--- Part 2: Phi_n decomposition ---\n")

print("""
Phi_n = sum_i [sum_{j!=i} 1/(lambda_i - lambda_j)]^2

Let's expand:
Phi_n = sum_i sum_{j!=i} 1/(lambda_i - lambda_j)^2
      + sum_i sum_{j!=i, k!=i, j!=k} 1/((lambda_i-lambda_j)(lambda_i-lambda_k))

= A + B where:
  A = sum_{i!=j} 1/(lambda_i-lambda_j)^2 = 2 * Sm2
  B = sum_i sum_{j!=i,k!=i,j!=k} 1/((lambda_i-lambda_j)(lambda_i-lambda_k))

Now A = 2*Sm2. What about B?

By partial fractions:
  1/((lambda_i-lambda_j)(lambda_i-lambda_k))
  = [1/(lambda_k-lambda_j)] * [1/(lambda_i-lambda_j) - 1/(lambda_i-lambda_k)]

Summing over i (with i != j, i != k):
  sum_{i!=j,k} 1/((lambda_i-lambda_j)(lambda_i-lambda_k))
  = [1/(lambda_k-lambda_j)] * [sum_{i!=j,k} 1/(lambda_i-lambda_j) - sum_{i!=j,k} 1/(lambda_i-lambda_k)]
  = [1/(lambda_k-lambda_j)] * [(H(lambda_j) - 1/(lambda_j-lambda_k)) - (H(lambda_k) - 1/(lambda_k-lambda_j))]
  = [1/(lambda_k-lambda_j)] * [H(lambda_j) - H(lambda_k) + 2/(lambda_k-lambda_j)]

Hmm, this is getting circular. Let me try a different approach.
""")

# Actually, let me just use the identity:
# Phi_n = 2*Sm2 + B
# Verify the decomposition numerically

print("Phi_n vs 2*Sm2 for various configurations:")
for roots in [np.array([-2.,0.,2.]), np.array([-3.,-1.,1.,3.]),
              np.array([-4.,-1.,0.,2.,5.])]:
    phi = Phi_n(roots)
    sm2 = Sm2(roots)
    B = phi - 2*sm2
    print(f"  roots: {roots}")
    print(f"  Phi = {phi:.6f}, 2*Sm2 = {2*sm2:.6f}, B = {B:.6f}")
    print(f"  B/Phi = {B/phi:.4f}")
    print()

# For n=2: B = 0 since there's only one pair
# For n=3: B != 0

# =============================================================
# Part 3: Key relationship: Phi * S2
# =============================================================
print("\n--- Part 3: Product Phi_n * S2 ---\n")

print("Study Phi_n * S2 under MSS:")
np.random.seed(42)
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

        ps = Phi_n(p_roots) * S2(p_roots)
        qs = Phi_n(q_roots) * S2(q_roots)
        rs = Phi_n(r_roots) * S2(r_roots)

        print(f"  Trial {trial}: Phi*S2: p={ps:.4f}, q={qs:.4f}, r={rs:.4f}")
    except:
        continue

# =============================================================
# Part 4: Cauchy-Schwarz on S2 and Phi
# =============================================================
print("\n\n--- Part 4: Cauchy-Schwarz via S2 ---\n")

print("""
Since S2(r) = S2(p) + S2(q), we can think of MSS convolution
as ADDING the "spread" parameters.

Now, by Cauchy-Schwarz:
  S2(r) = S2(p) + S2(q)

And we want: 1/Phi(r) >= 1/Phi(p) + 1/Phi(q)

If Phi(r) and S2(r) were related by Phi(r) = c/S2(r) (as for n=2 where
Phi_2 = 2/gap^2 = 2*n/S2 ... wait for n=2: S2 = (a-b)^2, Phi = 2/(a-b)^2 = 2/S2.
So Phi * S2 = 2 for n=2.)

For n >= 3, Phi * S2 is NOT constant. But maybe there's a lower bound.

BY CAUCHY-SCHWARZ (applied to vectors in R^n):
  Phi_n * (something related to S2) >= (something)^2

Actually, the Cauchy-Schwarz inequality gives:
  (sum_i a_i b_i)^2 <= (sum_i a_i^2)(sum_i b_i^2)

With a_i = H_i and b_i = 1:
  (sum_i H_i)^2 <= (sum_i H_i^2)(sum_i 1) = n * Phi_n
But sum_i H_i = 0, so this is trivial.

With a_i = H_i * lambda_i and b_i = 1:
  (sum_i H_i * lambda_i)^2 <= Phi_n * (sum_i lambda_i^2)

What is sum_i H_i * lambda_i?
""")

# Compute sum_i H_i * lambda_i
np.random.seed(42)
print("Computing sum_i H_i * lambda_i:")
for roots in [np.array([-2.,0.,2.]), np.array([-3.,-1.,1.,3.]),
              np.array([-4.,-1.,0.,2.,5.])]:
    H = H_values(roots)
    S = np.sum(H * roots)
    # WLOG center at 0
    roots_c = roots - np.mean(roots)
    H_c = H_values(roots_c)  # H is translation-invariant
    S_c = np.sum(H_c * roots_c)
    n = len(roots)
    print(f"  roots: {roots}")
    print(f"  sum H_i lambda_i = {S:.6f}")
    print(f"  (centered) sum H_i lambda_i = {S_c:.6f}")
    # Note: H is translation invariant since d_{ij} = lambda_i - lambda_j is
    # But lambda_i changes, so sum H_i lambda_i is NOT translation invariant
    # Unless... sum H_i = 0, so sum H_i * (lambda_i + c) = sum H_i * lambda_i
    # AH YES! Since sum H_i = 0, the sum H_i lambda_i IS translation invariant.
    print(f"  (Should be same as uncentered: {abs(S - S_c) < 1e-10})")
    print()

# =============================================================
# Part 5: The key identity sum H_i lambda_i = -(n-1)
# =============================================================
print("\n--- Part 5: Identity for sum H_i * lambda_i ---\n")

# sum_i H_i lambda_i = sum_i lambda_i sum_{j!=i} 1/(lambda_i - lambda_j)
# = sum_{i!=j} lambda_i / (lambda_i - lambda_j)
# = sum_{i<j} [lambda_i/(lambda_i-lambda_j) + lambda_j/(lambda_j-lambda_i)]
# = sum_{i<j} [lambda_i/(lambda_i-lambda_j) - lambda_j/(lambda_i-lambda_j)]
# = sum_{i<j} (lambda_i - lambda_j)/(lambda_i-lambda_j)
# = sum_{i<j} 1
# = n(n-1)/2

# Wait that's wrong. Let me redo:
# sum_{i!=j} lambda_i / (lambda_i - lambda_j)
# For each ordered pair (i,j) with i != j:
# = sum_{i<j} [lambda_i/(lambda_i-lambda_j) + lambda_j/(lambda_j-lambda_i)]
# = sum_{i<j} [(lambda_i - lambda_j + lambda_j)/(lambda_i-lambda_j) + lambda_j/(lambda_j-lambda_i)]
# = sum_{i<j} [1 + lambda_j/(lambda_i-lambda_j) + lambda_j/(lambda_j-lambda_i)]
# = sum_{i<j} [1 + lambda_j/(lambda_i-lambda_j) - lambda_j/(lambda_i-lambda_j)]
# = sum_{i<j} 1
# = n(n-1)/2

print("Verifying: sum_i H_i * lambda_i = n(n-1)/2")
for roots in [np.array([-2.,0.,2.]), np.array([-3.,-1.,1.,3.]),
              np.array([-4.,-1.,0.,2.,5.]), np.array([0.1, 0.5, 1.2, 3.7, 8.1])]:
    H = H_values(roots)
    n = len(roots)
    val = np.sum(H * roots)
    expected = n*(n-1)/2
    print(f"  n={n}, roots={roots}: sum H_i*lambda_i = {val:.6f}, expected {expected:.1f}")

# WRONG! Let me recheck for n=3:
# roots = [-2, 0, 2]
# H = [-0.75, 0, 0.75]
# sum H_i * lambda_i = (-0.75)(-2) + 0*0 + 0.75*2 = 1.5 + 1.5 = 3.0
# n(n-1)/2 = 3. YES!

# Let me verify my algebra:
# Actually:
# sum_{i,j: i!=j} lambda_i / (lambda_i - lambda_j)
# = sum_{i<j} [lambda_i/(lambda_i-lambda_j) + lambda_j/(lambda_j-lambda_i)]
# The second term: lambda_j / (lambda_j - lambda_i) = lambda_j / (-(lambda_i - lambda_j))
# So sum = sum_{i<j} (lambda_i - lambda_j)/(lambda_i - lambda_j) = sum_{i<j} 1
# Wait: lambda_i/(lambda_i-lambda_j) + lambda_j/(lambda_j-lambda_i)
# = lambda_i/(lambda_i-lambda_j) - lambda_j/(lambda_i-lambda_j)
# = (lambda_i - lambda_j)/(lambda_i - lambda_j) = 1
# So sum = number of pairs = n(n-1)/2. Confirmed!

print("\nIDENTITY CONFIRMED: sum_i H_p(lambda_i) * lambda_i = n(n-1)/2")
print("This is a UNIVERSAL identity, independent of the specific roots!")

# =============================================================
# Part 6: Higher moment identities
# =============================================================
print("\n\n--- Part 6: Higher moment identities ---\n")

print("Computing sum_i H_i * lambda_i^k for various k:")
for roots in [np.array([-2.,0.,2.]), np.array([-3.,-1.,1.,3.]),
              np.array([-1., 0., 2., 5.])]:
    H = H_values(roots)
    n = len(roots)
    print(f"\n  roots = {roots} (n={n}):")
    for k in range(4):
        val = np.sum(H * roots**k)
        print(f"    k={k}: sum H_i * lambda_i^{k} = {val:.6f}")

# k=0: sum H_i = 0 (known)
# k=1: sum H_i lambda_i = n(n-1)/2 (just proved)
# k=2: ?

# sum_{i!=j} lambda_i^2 / (lambda_i - lambda_j)
# = sum_{i<j} [lambda_i^2/(lambda_i-lambda_j) + lambda_j^2/(lambda_j-lambda_i)]
# = sum_{i<j} [lambda_i^2/(lambda_i-lambda_j) - lambda_j^2/(lambda_i-lambda_j)]
# = sum_{i<j} (lambda_i^2 - lambda_j^2)/(lambda_i-lambda_j)
# = sum_{i<j} (lambda_i + lambda_j)
# = (n-1) sum_i lambda_i = (n-1) * p_1

print("\n\nVerifying: sum_i H_i * lambda_i^2 = (n-1) * p_1 where p_1 = sum lambda_i")
for roots in [np.array([-2.,0.,2.]), np.array([-3.,-1.,1.,3.]),
              np.array([-1., 0., 2., 5.])]:
    H = H_values(roots)
    n = len(roots)
    val = np.sum(H * roots**2)
    expected = (n-1) * np.sum(roots)
    print(f"  n={n}, roots={roots}: sum H_i*lambda_i^2 = {val:.6f}, (n-1)*p_1 = {expected:.6f}")

# For k=3: sum_{i!=j} lambda_i^3 / (lambda_i - lambda_j)
# = sum_{i<j} (lambda_i^3 - lambda_j^3)/(lambda_i-lambda_j)
# = sum_{i<j} (lambda_i^2 + lambda_i*lambda_j + lambda_j^2)
# = (n-1)*p_2 + (p_1^2 - p_2)/2  ... hmm

print("\nFor k=3:")
for roots in [np.array([-2.,0.,2.]), np.array([-3.,-1.,1.,3.]),
              np.array([-1., 0., 2., 5.])]:
    H = H_values(roots)
    n = len(roots)
    val = np.sum(H * roots**3)
    p_1 = np.sum(roots)
    p_2 = np.sum(roots**2)
    expected = (n-1)*p_2 + (p_1**2 - p_2)/2
    # Actually: sum_{i<j} (lambda_i^2 + lambda_i*lambda_j + lambda_j^2)
    manual = sum(roots[i]**2 + roots[i]*roots[j] + roots[j]**2
                 for i in range(n) for j in range(i+1, n))
    print(f"  n={n}, roots={roots}: sum H_i*lambda_i^3 = {val:.6f}, manual sum = {manual:.6f}")

# =============================================================
# Part 7: Phi in terms of Sm2 and Sm11
# =============================================================
print("\n\n--- Part 7: Phi_n decomposition into Sm2 and cross terms ---\n")

# Phi_n = sum_i (sum_{j!=i} 1/d_{ij})^2
# = sum_i sum_{j!=i} 1/d_{ij}^2 + sum_i sum_{j!=i,k!=i,j!=k} 1/(d_{ij} * d_{ik})

# First part: sum_i sum_{j!=i} 1/d_{ij}^2 = sum_{i!=j} 1/d_{ij}^2 = 2*Sm2

# Cross part: B = sum_i sum_{j!=i,k!=i,j!=k} 1/(d_{ij}*d_{ik})
# = sum_{j!=k} sum_{i!=j,i!=k} 1/(d_{ij}*d_{ik})

# For a given pair (j,k) with j != k:
# sum_{i!=j,i!=k} 1/((lambda_i-lambda_j)(lambda_i-lambda_k))

# Using partial fractions:
# 1/((x-a)(x-b)) = 1/(a-b) * (1/(x-a) - 1/(x-b))  for a != b

# So sum_{i!=j,i!=k} 1/((lambda_i-lambda_j)(lambda_i-lambda_k))
# = 1/(lambda_j-lambda_k) * sum_{i!=j,i!=k} [1/(lambda_i-lambda_j) - 1/(lambda_i-lambda_k)]
# = 1/(lambda_j-lambda_k) * [(H_j - 1/(lambda_j-lambda_k)) - (H_k - 1/(lambda_k-lambda_j))]
# = 1/(lambda_j-lambda_k) * [H_j - H_k - 1/(lambda_j-lambda_k) - 1/(lambda_k-lambda_j)]
# = 1/(lambda_j-lambda_k) * [H_j - H_k - 2/(lambda_j-lambda_k)]
# = (H_j - H_k)/(lambda_j - lambda_k) - 2/(lambda_j-lambda_k)^2

# So B = sum_{j!=k} [(H_j - H_k)/(lambda_j - lambda_k) - 2/(lambda_j-lambda_k)^2]
# = sum_{j!=k} (H_j - H_k)/(lambda_j - lambda_k) - 2 * sum_{j!=k} 1/(lambda_j-lambda_k)^2
# = 2 * sum_{j<k} (H_j - H_k)/(lambda_j - lambda_k) - 4*Sm2

# And Phi_n = 2*Sm2 + B = 2*Sm2 + 2 * sum_{j<k} (H_j-H_k)/(lambda_j-lambda_k) - 4*Sm2
# = 2 * sum_{j<k} (H_j - H_k)/(lambda_j - lambda_k) - 2*Sm2

print("Verifying: Phi_n = 2 * sum_{j<k} (H_j-H_k)/(lambda_j-lambda_k) - 2*Sm2")
for roots in [np.array([-2.,0.,2.]), np.array([-3.,-1.,1.,3.]),
              np.array([-1., 0., 2., 5.])]:
    H = H_values(roots)
    n = len(roots)
    phi = Phi_n(roots)
    sm2 = Sm2(roots)
    cross = 0
    for j in range(n):
        for k in range(j+1, n):
            cross += (H[j] - H[k]) / (roots[j] - roots[k])

    computed = 2*cross - 2*sm2
    print(f"  roots={roots}: Phi={phi:.6f}, 2*cross - 2*Sm2 = {computed:.6f}, match: {abs(phi-computed)<1e-10}")

# So Phi_n = 2 * T - 2*Sm2 where T = sum_{j<k} (H_j - H_k)/(lambda_j - lambda_k)

# Now T = sum_{j<k} (H_j - H_k)/(lambda_j - lambda_k)
# = sum_{j<k} sum_{i!=j} 1/((lambda_j-lambda_i)(lambda_j-lambda_k))
#   - sum_{j<k} sum_{i!=k} 1/((lambda_k-lambda_i)(lambda_j-lambda_k))  ... getting complex

# =============================================================
# Part 8: The DISCRETE LAPLACIAN connection
# =============================================================
print("\n\n--- Part 8: Discrete Laplacian ---\n")

print("""
OBSERVATION: T = sum_{j<k} (H_j - H_k) / (lambda_j - lambda_k)

This is a discrete analogue of the integral of (H(x) - H(y))/(x-y) dx dy,
which is related to the Dirichlet form of H.

In fact, let us define the "Cauchy bilinear form" or "Hilbert-Schmidt form":
  <f, g>_C = sum_{i<j} (f(lambda_i) - f(lambda_j))(g(lambda_i) - g(lambda_j)) / (lambda_i - lambda_j)^2

Then T = <H, id>_C where id(x) = x (identity function).

Hmm, that's not quite right either. Let me think about this differently.

The key object (H_j - H_k)/(lambda_j - lambda_k) is like a "discrete derivative"
of H evaluated at the pair (j,k). Since H is already a discrete derivative
of the log-potential W, T involves the SECOND derivative of W.
""")

# =============================================================
# Part 9: Scale-invariant form of the conjecture
# =============================================================
print("\n--- Part 9: Scale-invariant form ---\n")

# Since Phi_n scales as lambda^{-2} under root scaling lambda_i -> c*lambda_i,
# and S2 scales as lambda^2, the product Phi_n * S2 is SCALE-INVARIANT.

# The conjecture 1/Phi_r >= 1/Phi_p + 1/Phi_q can be rewritten:
# Phi_r <= Phi_p * Phi_q / (Phi_p + Phi_q)

# Multiply both sides by S2_r = S2_p + S2_q:
# Phi_r * S2_r <= (Phi_p * Phi_q / (Phi_p + Phi_q)) * (S2_p + S2_q)

# Using the scale-invariant quantity Q = Phi * S2:
# Q_r = Phi_r * S2_r

# The conjecture becomes:
# Q_r <= (Phi_p * Phi_q / (Phi_p + Phi_q)) * (S2_p + S2_q)
# This doesn't simplify cleanly.

# BUT: in terms of Q = Phi * S2:
# 1/Phi_r = S2_r / Q_r = (S2_p + S2_q) / Q_r
# 1/Phi_p = S2_p / Q_p
# 1/Phi_q = S2_q / Q_q

# So the conjecture is:
# (S2_p + S2_q) / Q_r >= S2_p / Q_p + S2_q / Q_q

print("Scale-invariant form using Q = Phi * S2:")
print("Conjecture: (S2_p + S2_q) / Q_r >= S2_p/Q_p + S2_q/Q_q\n")

np.random.seed(42)
for trial in range(8):
    n = 4
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
        Qp = Phi_n(p_roots) * S2(p_roots)
        Qq = Phi_n(q_roots) * S2(q_roots)
        Qr = Phi_n(r_roots) * S2(r_roots)
        s2p = S2(p_roots)
        s2q = S2(q_roots)

        lhs = (s2p + s2q) / Qr
        rhs = s2p/Qp + s2q/Qq
        print(f"  Trial {trial}: Q_p={Qp:.4f}, Q_q={Qq:.4f}, Q_r={Qr:.4f}")
        print(f"    LHS={(s2p+s2q)/Qr:.6f}, RHS={s2p/Qp + s2q/Qq:.6f}, gap={lhs-rhs:.6e}")
    except:
        continue

# =============================================================
# Part 10: Key test: Is Q_r <= max(Q_p, Q_q)?
# =============================================================
print("\n\n--- Part 10: Monotonicity of Q = Phi * S2 ---\n")

# If Q_r <= Q_p + Q_q, combined with S2 additivity...
# The conjecture (S2_p+S2_q)/Q_r >= S2_p/Q_p + S2_q/Q_q
# This is Titu's lemma / Engel form if we think of it as:
# sum_i a_i^2 / b_i >= (sum a_i)^2 / (sum b_i)
# with a_i = sqrt(S2_i), b_i = Q_i

# By Cauchy-Schwarz (Titu/Engel):
# S2_p/Q_p + S2_q/Q_q >= (sqrt(S2_p) + sqrt(S2_q))^2 / (Q_p + Q_q)

# And we need (S2_p + S2_q)/Q_r >= this.
# This would follow if Q_r <= Q_p + Q_q AND (S2_p+S2_q)(Q_p+Q_q) >= (sqrt(S2_p)+sqrt(S2_q))^2 * Q_r
# The second condition is: (S2_p+S2_q)(Q_p+Q_q) >= (S2_p+S2_q+2sqrt(S2_p*S2_q)) * Q_r

# Hmm, let me just check if Q_r <= Q_p + Q_q holds:
print("Is Q_r <= Q_p + Q_q?")
np.random.seed(42)
violations_Q = 0
for trial in range(5000):
    n = 4
    p_roots = np.sort(np.random.randn(n) * 2)
    q_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(p_roots)) < 0.1:
        p_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(q_roots)) < 0.1:
        q_roots = np.sort(np.random.randn(n) * 2)
    try:
        r_roots = mss_convolve_roots(p_roots, q_roots)
        if np.min(np.diff(r_roots)) < 1e-6:
            continue
        Qp = Phi_n(p_roots) * S2(p_roots)
        Qq = Phi_n(q_roots) * S2(q_roots)
        Qr = Phi_n(r_roots) * S2(r_roots)
        if Qr > Qp + Qq + 1e-8:
            violations_Q += 1
    except:
        continue
print(f"  {violations_Q} violations of Q_r <= Q_p + Q_q")

# Also check: Q_r <= alpha * max(Q_p, Q_q)?
print("\nIs Q_r <= max(Q_p, Q_q)?")
violations_Qmax = 0
for trial in range(5000):
    n = 4
    p_roots = np.sort(np.random.randn(n) * 2)
    q_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(p_roots)) < 0.1:
        p_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(q_roots)) < 0.1:
        q_roots = np.sort(np.random.randn(n) * 2)
    try:
        r_roots = mss_convolve_roots(p_roots, q_roots)
        if np.min(np.diff(r_roots)) < 1e-6:
            continue
        Qp = Phi_n(p_roots) * S2(p_roots)
        Qq = Phi_n(q_roots) * S2(q_roots)
        Qr = Phi_n(r_roots) * S2(r_roots)
        if Qr > max(Qp, Qq) + 1e-8:
            violations_Qmax += 1
    except:
        continue
print(f"  {violations_Qmax} violations of Q_r <= max(Q_p, Q_q)")

# What about Q_r vs harmonic mean?
print("\nQ_r vs 2*Qp*Qq/(Qp+Qq) (harmonic mean):")
min_ratio = float('inf')
for trial in range(5000):
    n = 4
    p_roots = np.sort(np.random.randn(n) * 2)
    q_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(p_roots)) < 0.1:
        p_roots = np.sort(np.random.randn(n) * 2)
    while np.min(np.diff(q_roots)) < 0.1:
        q_roots = np.sort(np.random.randn(n) * 2)
    try:
        r_roots = mss_convolve_roots(p_roots, q_roots)
        if np.min(np.diff(r_roots)) < 1e-6:
            continue
        Qp = Phi_n(p_roots) * S2(p_roots)
        Qq = Phi_n(q_roots) * S2(q_roots)
        Qr = Phi_n(r_roots) * S2(r_roots)
        harm = 2*Qp*Qq/(Qp+Qq)
        ratio = Qr / harm if harm > 0 else float('inf')
        min_ratio = min(min_ratio, ratio)
    except:
        continue
print(f"  min(Q_r / harmonic_mean(Q_p,Q_q)) = {min_ratio:.6f}")

print("\nDone.")
