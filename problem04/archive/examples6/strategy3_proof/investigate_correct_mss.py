"""
Correct MSS finite free additive convolution.

The CORRECT formula from Marcus-Spielman-Srivastava (2015), Interlacing Families II:

For monic polynomials p(x) = sum_{k=0}^n (-1)^k e_k(p) x^{n-k} and q similarly,

(p boxplus_n q)(x) = sum_{k=0}^n (-1)^k e_k(r) x^{n-k}

where (Definition 4.1, Theorem 4.5):

e_k(r) = sum_{j=0}^{k} binom(n-j, k-j) / binom(n, k) * (-1)^{k-j} * ... no

Actually from the MSS definition of "additive convolution" for REAL-ROOTED polynomials:

The correct formula uses the MIXED CHARACTERISTIC POLYNOMIAL.
For diagonal matrices A = diag(a_1,...,a_n) and B = diag(b_1,...,b_n):

mu[A + B](x) = expected char poly of A + U*B*U^T over random signing matrices U.

Wait no, that's a different convolution. Let me be precise.

MSS (2015) define the finite free ADDITIVE convolution as:
Given p,q of degree n, p boxplus_n q is computed via:

Method 1: Expected characteristic polynomial over Haar unitary:
  (p boxplus_n q)(x) = E_U[det(xI - (A + UBU*))]

Method 2: Via the "additive mixed characteristic polynomial":
  mu[A, B](x) = sum_{S subset [n]} det(xI - A_S - B_{S^c})
  (appropriately normalized)

For computational purposes, Method 1 with high-accuracy MC is reliable.
But we can also use an EXACT formula.

From Anderson (2014), the exact formula for the expected characteristic polynomial
of A + UBU* over Haar U on U(n) is:

E_U[det(xI - A - UBU*)] = sum_{sigma in S_n} (1/n!) * ... no, this isn't tractable.

Actually, the key identity (Theorem 1.1 of MSS 2015, also Proposition 2.3 of
Ravichandran survey):

For the SYMMETRIC ADDITIVE CONVOLUTION (over orthogonal group):
  This is NOT what we want.

For the UNITARY ADDITIVE CONVOLUTION:
  (p boxplus_n q)(x) = sum_{k=0}^n (-1)^k e_k x^{n-k}
  where
  e_k = sum_j (-1)^{k-j} * C(k,j) * C(n-j, k-j)^{-1} * ...

OK I realize the issue. The "hat convolution" IS correct for the specific definition
of finite free convolution in Marcus (2021), BUT it's a DIFFERENT convolution from
the expected characteristic polynomial over Haar unitary.

Let me check: for the MARCUS finite free convolution (which he calls D_n[p,q]):
  hat_e_k(D_n[p,q]) = sum_{j=0}^k hat_e_j(p) * hat_e_{k-j}(q)
  where hat_e_k = e_k / C(n,k)

This is the "finite free additive convolution" from Definition 2.1 in Marcus (2021),
which equals the expected characteristic polynomial of A + UBU* where:
- A,B are diagonal with eigenvalues = roots of p,q
- U is from the SYMMETRIC GROUP (random permutation matrix), not Haar unitary!

The HAAR UNITARY version is DIFFERENT and is what appears in the MSS interlacing
families paper. This is the one relevant for our problem.

For the UNITARY version, the formula is more complex. Let me implement the
correct one using integration or the Weingarten function approach.

Actually, upon further reflection, the relevant convolution for the Fisher
information problem is the one defined via subordination, which corresponds
to the FREE ADDITIVE CONVOLUTION in free probability (infinite n limit).

For FINITE n, the relevant convolution is the one from MSS that preserves
real-rootedness, which IS the expected char poly over Haar unitary.

Let me implement this correctly using numerical integration.
"""

import numpy as np
from math import comb, factorial
from itertools import combinations
from scipy import integrate

def elem_sym_poly(roots, k):
    n = len(roots)
    if k == 0: return 1.0
    if k > n: return 0.0
    return sum(np.prod([roots[i] for i in subset])
               for subset in combinations(range(n), k))


def boxplus_haar_mc(roots_p, roots_q, samples=500000):
    """Compute p boxplus_n q via Monte Carlo over Haar unitary group."""
    n = len(roots_p)
    A = np.diag(np.array(roots_p, dtype=float))
    B = np.diag(np.array(roots_q, dtype=float))

    sum_ek = np.zeros(n+1)

    for _ in range(samples):
        # Random Haar unitary
        Z = (np.random.randn(n, n) + 1j * np.random.randn(n, n)) / np.sqrt(2)
        Q, R_mat = np.linalg.qr(Z)
        d = np.diagonal(R_mat)
        ph = d / np.abs(d)
        U = Q * ph[np.newaxis, :]

        M = A + U @ B @ U.conj().T
        eigs = np.sort(np.linalg.eigvalsh(M))

        for k in range(n+1):
            sum_ek[k] += elem_sym_poly(eigs, k)

    avg_ek = sum_ek / samples

    # Build polynomial and find roots
    coeffs = np.zeros(n+1)
    for k in range(n+1):
        coeffs[k] = (-1)**k * avg_ek[k]
    roots_r = np.sort(np.real(np.roots(coeffs)))
    return roots_r, avg_ek


def boxplus_permutation(roots_p, roots_q):
    """Compute the PERMUTATION convolution (Marcus D_n) exactly.

    D_n[p,q](x) = (1/n!) sum_{sigma in S_n} prod_{i=1}^{n} (x - a_i - b_{sigma(i)})
    """
    n = len(roots_p)
    from itertools import permutations

    if n > 7:  # Too slow for large n
        return None, None

    sum_ek = np.zeros(n+1)
    count = 0

    for perm in permutations(range(n)):
        shifted_roots = np.array([roots_p[i] + roots_q[perm[i]] for i in range(n)])
        for k in range(n+1):
            sum_ek[k] += elem_sym_poly(shifted_roots, k)
        count += 1

    avg_ek = sum_ek / count

    coeffs = np.zeros(n+1)
    for k in range(n+1):
        coeffs[k] = (-1)**k * avg_ek[k]
    roots_r = np.sort(np.real(np.roots(coeffs)))
    return roots_r, avg_ek


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


# ================================================================
# First: compare HAAR vs PERMUTATION convolutions for small n
# ================================================================

np.random.seed(42)

print("="*70)
print("COMPARING HAAR UNITARY vs PERMUTATION CONVOLUTIONS")
print("="*70)

test_cases = [
    (np.array([-1.0, 1.0]), np.array([-2.0, 2.0]), "n=2 sym"),
    (np.array([0.0, 1.0]), np.array([0.0, 3.0]), "n=2 asym"),
    (np.array([-2.0, 0.0, 2.0]), np.array([-3.0, 0.0, 3.0]), "n=3 sym"),
    (np.array([-1.0, 0.5, 2.0]), np.array([-0.5, 1.0, 3.0]), "n=3 asym"),
    (np.array([-3.0, -1.0, 1.0, 3.0]), np.array([-2.0, -0.5, 0.5, 2.0]), "n=4 sym"),
]

for roots_p, roots_q, label in test_cases:
    print(f"\n--- {label} ---")

    # Permutation convolution (exact)
    roots_perm, ek_perm = boxplus_permutation(roots_p, roots_q)

    # Haar unitary MC
    roots_haar, ek_haar = boxplus_haar_mc(roots_p, roots_q, 200000)

    if roots_perm is not None:
        print(f"  Permutation roots: {roots_perm}")
        print(f"  Haar MC roots:     {roots_haar}")
        print(f"  Difference:        {np.max(np.abs(roots_perm - roots_haar)):.6f}")

        # Inequality check for permutation
        Phi_p = Phi_n(roots_p)
        Phi_q = Phi_n(roots_q)
        Phi_r_perm = Phi_n(roots_perm)
        Phi_r_haar = Phi_n(roots_haar)

        ineq_perm = 1/Phi_r_perm - 1/Phi_p - 1/Phi_q
        ineq_haar = 1/Phi_r_haar - 1/Phi_p - 1/Phi_q

        print(f"  Phi_p={Phi_p:.6f}, Phi_q={Phi_q:.6f}")
        print(f"  Phi_r(perm)={Phi_r_perm:.6f}, ineq(perm)={ineq_perm:.8f}")
        print(f"  Phi_r(haar)={Phi_r_haar:.6f}, ineq(haar)={ineq_haar:.8f}")

# For n=2: the permutation convolution is:
# D_2[p,q](x) = (1/2)[(x-a-c)(x-b-d) + (x-a-d)(x-b-c)]
# = (1/2)[(x^2 - (a+c+b+d)x + (a+c)(b+d)) + (x^2 - (a+d+b+c)x + (a+d)(b+c))]
# = x^2 - (a+b+c+d)x + (1/2)[(a+c)(b+d) + (a+d)(b+c)]
# = x^2 - (a+b+c+d)x + (1/2)[ab+ad+bc+cd + ab+ac+bd+cd]
# = x^2 - (a+b+c+d)x + ab + (1/2)(ac+ad+bc+bd) + cd
# = x^2 - (a+b+c+d)x + ab + cd + (1/2)(a+b)(c+d)

# For the Haar convolution: this IS the same for n=2!
# Because for n=2, the average over U(2) equals the average over permutations
# (the Weingarten function simplifies).

# Actually, for n=2, BOTH convolutions give the same result.
# The difference appears at n >= 3.

print("\n\n" + "="*70)
print("WHICH CONVOLUTION SATISFIES THE FISHER INEQUALITY?")
print("="*70)

# Test the permutation convolution
print("\nPermutation convolution tests:")
perm_pass = 0
perm_total = 0
for trial in range(100):
    n = np.random.choice([2, 3, 4, 5])
    if n > 5:
        continue
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5

    try:
        roots_r, _ = boxplus_permutation(roots_p, roots_q)
        if roots_r is None:
            continue
        if not np.all(np.isreal(roots_r)):
            continue

        gaps = np.diff(roots_r)
        if np.any(gaps < 0.01):
            continue

        ineq = 1/Phi_n(roots_r) - 1/Phi_n(roots_p) - 1/Phi_n(roots_q)
        perm_total += 1
        if ineq >= -1e-8:
            perm_pass += 1
        else:
            if perm_total <= 5 or ineq < -0.1:
                print(f"  FAIL n={n}: ineq={ineq:.8f}")
                print(f"    p={roots_p}, q={roots_q}, r={roots_r}")
    except:
        pass

print(f"\nPermutation: {perm_pass}/{perm_total} pass")

# Test the Haar convolution
print("\nHaar convolution tests (200k MC):")
haar_pass = 0
haar_total = 0
for trial in range(50):
    n = np.random.choice([2, 3, 4])
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5

    try:
        roots_r, _ = boxplus_haar_mc(roots_p, roots_q, 200000)
        if not np.all(np.isreal(roots_r)):
            continue

        gaps = np.diff(roots_r)
        if np.any(gaps < 0.01):
            continue

        ineq = 1/Phi_n(roots_r) - 1/Phi_n(roots_p) - 1/Phi_n(roots_q)
        haar_total += 1
        if ineq >= -1e-6:  # Larger tolerance for MC
            haar_pass += 1
        else:
            print(f"  FAIL n={n}: ineq={ineq:.8f}")
    except:
        pass

print(f"\nHaar: {haar_pass}/{haar_total} pass")


# ================================================================
# Check what the ACTUAL MSS boxplus_n is
# ================================================================
print("\n\n" + "="*70)
print("THE KEY QUESTION: WHAT IS boxplus_n?")
print("="*70)
print("""
In the MSS (2015) paper, the convolution IS the expected char poly
over Haar unitary. But in Marcus's 2021 survey on "Finite Free Probability",
he defines a DIFFERENT convolution using the normalized coefficients:
  hat_e_k(p boxplus q) = sum_j hat_e_j(p) * hat_e_{k-j}(q)

These are THE SAME for n=2 but DIFFERENT for n >= 3.

The Marcus (2021) "finite free additive convolution" is actually the
expected characteristic polynomial over RANDOM PERMUTATION MATRICES
(i.e., over the symmetric group S_n), not over the full unitary group.

For the Fisher information inequality, the relevant convolution is
the one that arises from the subordination structure. We need to check
which one that is.

Actually, looking at this more carefully: the "finite free convolution"
in the sense of Voiculescu/free probability is the N->infinity limit.
For FINITE n, the MSS convolution (over Haar unitary) is the relevant one.
""")

# Verify: for n=2, permutation = Haar (since S_2 = Z_2 and averages agree)
print("Verification: n=2, permutation == Haar:")
roots_p = np.array([-1.0, 3.0])
roots_q = np.array([-2.0, 5.0])

r_perm, _ = boxplus_permutation(roots_p, roots_q)
# Exact: D_2 has roots mean +/- sqrt((a-b)^2/4 + (c-d)^2/4)
mean = sum(roots_p) + sum(roots_q)
mean /= 2
s = roots_p[1] - roots_p[0]
t = roots_q[1] - roots_q[0]
spread = np.sqrt(s**2/4 + t**2/4)
r_exact = np.array([mean - spread, mean + spread])
print(f"  Permutation: {r_perm}")
print(f"  n=2 formula: {r_exact}")
print(f"  Match: {np.allclose(r_perm, r_exact)}")
