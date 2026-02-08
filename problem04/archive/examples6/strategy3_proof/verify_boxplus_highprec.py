"""
High-precision Monte Carlo for the borderline case p=q={-3,-1,1,3}.
The 100K MC showed a tiny negative gap (-2.69e-4) which could be MC noise.
Run with 1M samples to determine if this is real.
"""
import numpy as np
from math import factorial, comb
from itertools import combinations

np.random.seed(12345)

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

def haar_unitary(n):
    Z = (np.random.randn(n, n) + 1j * np.random.randn(n, n)) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    d = np.diag(R)
    ph = d / np.abs(d)
    return Q @ np.diag(ph)

def boxplus_correct(roots_p, roots_q):
    """Formula A/B: the correct MSS formula."""
    n = len(roots_p)
    poly_p = np.poly(roots_p)
    poly_q = np.poly(roots_q)
    a = poly_p
    b = poly_q
    c = np.zeros(n + 1)
    for k in range(n + 1):
        for i in range(k + 1):
            j = k - i
            if i <= n and j <= n:
                w = factorial(n - i) * factorial(n - j) / (factorial(n) * factorial(n - k))
                c[k] += w * a[i] * b[j]
    return np.sort(np.real(np.roots(c))), c

# Case: p = q = {-3, -1, 1, 3}
rp = np.array([-3.0, -1.0, 1.0, 3.0])
rq = np.array([-3.0, -1.0, 1.0, 3.0])

phi_p = Phi_n(rp)
phi_q = Phi_n(rq)
print(f"Phi_n(p) = {phi_p}, 1/Phi_n(p) = {1/phi_p}")
print(f"1/Phi(p) + 1/Phi(q) = {2/phi_p}")

# Formula result
roots_formula, poly_formula = boxplus_correct(rp, rq)
phi_formula = Phi_n(roots_formula)
gap_formula = 1/phi_formula - 2/phi_p
print(f"\nFormula A/B:")
print(f"  roots: {roots_formula}")
print(f"  Phi_n(r) = {phi_formula}")
print(f"  gap = {gap_formula:.10e}")

# High-precision Monte Carlo
N_samples = 1000000
print(f"\nRunning Monte Carlo with N={N_samples}...")
A = np.diag(rp.astype(complex))
B = np.diag(rq.astype(complex))
n = 4
poly_sum = np.zeros(n + 1)
for trial in range(N_samples):
    U = haar_unitary(n)
    M = A + U @ B @ U.conj().T
    poly_sum += np.real(np.poly(M))

poly_mc = poly_sum / N_samples
poly_mc[0] = 1.0
roots_mc = np.sort(np.real(np.roots(poly_mc)))
phi_mc = Phi_n(roots_mc)
gap_mc = 1/phi_mc - 2/phi_p

print(f"\nMonte Carlo (N={N_samples}):")
print(f"  poly: {poly_mc}")
print(f"  roots: {roots_mc}")
print(f"  Phi_n(r) = {phi_mc}")
print(f"  gap = {gap_mc:.10e}")

print(f"\nFormula poly: {poly_formula}")
print(f"MC poly:      {poly_mc}")
print(f"Diff:         {poly_formula - poly_mc}")
print(f"Max coeff err: {np.max(np.abs(poly_formula - poly_mc)):.6e}")

# Also compute standard error of the MC polynomial coefficients
# by batching
batch_size = 10000
n_batches = N_samples // batch_size
batch_polys = np.zeros((n_batches, n+1))
np.random.seed(12345)
for b in range(n_batches):
    poly_batch = np.zeros(n + 1)
    for _ in range(batch_size):
        U = haar_unitary(n)
        M = A + U @ B @ U.conj().T
        poly_batch += np.real(np.poly(M))
    batch_polys[b] = poly_batch / batch_size

se = np.std(batch_polys, axis=0) / np.sqrt(n_batches)
print(f"\nMC standard errors: {se}")
print(f"Formula - MC in units of SE: {(poly_formula - poly_mc) / (se + 1e-20)}")
