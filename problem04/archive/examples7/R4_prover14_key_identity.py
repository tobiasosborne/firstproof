"""
PROVER-14 Part 5: KEY IDENTITY and Proof Strategy

CRITICAL DISCOVERY: Phi_n = 2 * Sm2 = 2 * sum_{i<j} 1/(lambda_i - lambda_j)^2

This SIMPLIFIES Phi_n enormously! Instead of a sum of squares of sums,
it's just twice the sum of inverse squared differences.

PROOF of Phi_n = 2*Sm2:
  Phi_n = sum_i H_i^2 = sum_i (sum_{j!=i} 1/(lambda_i - lambda_j))^2
        = sum_i [sum_{j!=i} 1/d_ij^2 + sum_{j!=i,k!=i,j!=k} 1/(d_ij * d_ik)]

  The cross term: for each i, sum_{j!=i,k!=i,j!=k} 1/(d_ij * d_ik)
  = sum_{j!=i} (1/d_ij) * sum_{k!=i,k!=j} (1/d_ik)
  = sum_{j!=i} (1/d_ij) * (H_i - 1/d_ij)

  So cross_i = H_i * sum_{j!=i} 1/d_ij - sum_{j!=i} 1/d_ij^2
             = H_i^2 - sum_{j!=i} 1/d_ij^2

  Wait, that gives: H_i^2 = sum_{j!=i} 1/d_ij^2 + cross_i
  And cross_i = H_i^2 - sum_{j!=i} 1/d_ij^2
  This is circular.

  Let me try differently.

  B = sum_i sum_{j!=i,k!=i,j!=k} 1/(d_{ij} * d_{ik})

  Each triple (i,j,k) with j!=k, j!=i, k!=i appears once.
  By partial fractions:
    1/((x-a)(x-b)) = (1/(a-b)) * (1/(x-b) - 1/(x-a))

  For the triple (i,j,k):
    1/(d_{ij} * d_{ik}) = 1/((lambda_i-lambda_j)(lambda_i-lambda_k))
    = 1/(lambda_k-lambda_j) * (1/(lambda_i-lambda_j) - 1/(lambda_i-lambda_k))
    = 1/(lambda_k-lambda_j) * (1/d_{ij} - 1/d_{ik})

  Summing over all triples:
  B = sum over ordered (j,k) with j!=k, summing over i!=j,i!=k:
    = sum_{j!=k} 1/(lambda_k-lambda_j) * sum_{i!=j,k} (1/(lambda_i-lambda_j) - 1/(lambda_i-lambda_k))

  The inner sum = (sum_{i!=j} 1/(lambda_i-lambda_j) - 1/(lambda_k-lambda_j))
                - (sum_{i!=k} 1/(lambda_i-lambda_k) - 1/(lambda_j-lambda_k))
                = (H_j - 1/(lambda_k-lambda_j)) - (H_k + 1/(lambda_k-lambda_j))
                = H_j - H_k - 2/(lambda_k-lambda_j)

  Wait, 1/(lambda_j-lambda_k) = -1/(lambda_k-lambda_j), so:
  = H_j - 1/d_{jk} - H_k - 1/(lambda_j-lambda_k)
  Hmm, let me be more careful with signs.

  Let me just verify numerically and use the result.
"""
import numpy as np
from itertools import combinations
from math import factorial
import sympy as sp

print("="*70)
print("PROVER-14 Part 5: KEY IDENTITY Phi_n = 2*Sm2")
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

def Sm2(roots):
    """sum_{i<j} 1/(lambda_i - lambda_j)^2"""
    n = len(roots)
    return sum(1/(roots[i]-roots[j])**2 for i in range(n) for j in range(i+1, n))

def S2(roots):
    n = len(roots)
    return sum((roots[i]-roots[j])**2 for i in range(n) for j in range(i+1, n))

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
# Part 1: Verify Phi_n = 2*Sm2
# =============================================================
print("\n--- Part 1: Verify Phi_n = 2*Sm2 ---\n")

np.random.seed(42)
print("Verifying Phi_n = 2*Sm2 for various n and random roots:\n")
for n in [2, 3, 4, 5, 6, 7]:
    violations = 0
    for trial in range(100):
        roots = np.sort(np.random.randn(n) * 3)
        while np.min(np.diff(roots)) < 0.1:
            roots = np.sort(np.random.randn(n) * 3)
        phi = Phi_n(roots)
        sm2 = Sm2(roots)
        if abs(phi - 2*sm2) > 1e-8:
            violations += 1
            print(f"  VIOLATION at n={n}: Phi={phi:.10f}, 2*Sm2={2*sm2:.10f}")
    if violations == 0:
        # Print one example
        roots = np.sort(np.random.randn(n) * 3)
        while np.min(np.diff(roots)) < 0.1:
            roots = np.sort(np.random.randn(n) * 3)
        phi = Phi_n(roots)
        sm2 = Sm2(roots)
        print(f"  n={n}: Phi={phi:.10f}, 2*Sm2={2*sm2:.10f}, diff={phi-2*sm2:.2e} -> MATCH")

# =============================================================
# Part 2: PROVE Phi_n = 2*Sm2 symbolically
# =============================================================
print("\n\n--- Part 2: Symbolic proof of Phi_n = 2*Sm2 ---\n")

# For n=3
a, b, c = sp.symbols('a b c', real=True)
H_a = 1/(a-b) + 1/(a-c)
H_b = 1/(b-a) + 1/(b-c)
H_c = 1/(c-a) + 1/(c-b)

Phi_3 = sp.expand(H_a**2 + H_b**2 + H_c**2)
Sm2_3 = 1/(a-b)**2 + 1/(a-c)**2 + 1/(b-c)**2

diff_3 = sp.simplify(Phi_3 - 2*Sm2_3)
print(f"n=3: Phi_3 - 2*Sm2 = {diff_3}")

# For n=4
x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4', real=True)
roots_sym = [x1, x2, x3, x4]

H_sym = []
for i in range(4):
    hi = sum(1/(roots_sym[i] - roots_sym[j]) for j in range(4) if j != i)
    H_sym.append(hi)

Phi_4 = sum(h**2 for h in H_sym)
Sm2_4 = sum(1/(roots_sym[i] - roots_sym[j])**2 for i in range(4) for j in range(i+1, 4))

diff_4 = sp.simplify(sp.expand(Phi_4 - 2*Sm2_4))
print(f"n=4: Phi_4 - 2*Sm2 = {diff_4}")

print("""
PROOF of Phi_n = 2*Sm2 for all n:

Phi_n = sum_i H_i^2 = sum_i (sum_{j!=i} 1/(lambda_i - lambda_j))^2

Expand the square:
= sum_i [sum_{j!=i} 1/(lambda_i-lambda_j)^2 + 2*sum_{j<k, j!=i, k!=i} 1/((lambda_i-lambda_j)(lambda_i-lambda_k))]

= sum_i sum_{j!=i} 1/(lambda_i-lambda_j)^2 + 2*sum_i sum_{j<k, j!=i, k!=i} 1/((lambda_i-lambda_j)(lambda_i-lambda_k))

The first sum = sum_{i!=j} 1/(lambda_i-lambda_j)^2 = 2*Sm2.

So we need to show the second sum = 0:
  B = sum_i sum_{j<k, j!=i, k!=i} 1/((lambda_i-lambda_j)(lambda_i-lambda_k)) = 0

For each unordered triple {i,j,k}, the contribution to B comes from
three terms (one for each vertex i of the triple):
  1/((lambda_i-lambda_j)(lambda_i-lambda_k))
  + 1/((lambda_j-lambda_i)(lambda_j-lambda_k))
  + 1/((lambda_k-lambda_i)(lambda_k-lambda_j))

Let a = lambda_i, b = lambda_j, c = lambda_k.
The sum is:
  1/((a-b)(a-c)) + 1/((b-a)(b-c)) + 1/((c-a)(c-b))

Putting over common denominator (a-b)(b-c)(c-a):
= (b-c)/((a-b)(b-c)(a-c)) + (a-c)/((b-a)(b-c)(c-a)) + (a-b)/((c-a)(c-b)(a-b))

Wait, let me compute directly:
  (b-c)/((a-b)(a-c)(b-c)) ... no.

Common denominator = (a-b)(a-c)(b-c).
Numerator:
  (b-c) [from first term: (a-b)(a-c) in denom, missing factor is (b-c)]
  Actually, let D = (a-b)(b-c)(c-a).

  1/((a-b)(a-c)) = 1/((a-b)(-(c-a))) = -1/((a-b)(c-a)) = (b-c)/D
  1/((b-a)(b-c)) = 1/(-(a-b)(b-c)) = -(c-a)/D = (a-c)/D... hmm.

  Let me use D = (a-b)(b-c)(a-c).
  1/((a-b)(a-c)) = (b-c)/D
  1/((b-a)(b-c)) = 1/(-(a-b)(b-c)) = -(a-c)/D
  1/((c-a)(c-b)) = 1/(-(a-c)(-(b-c))) = 1/((a-c)(b-c)) = (a-b)/D

  Sum = [(b-c) - (a-c) + (a-b)] / D
      = [b-c-a+c+a-b] / D
      = 0 / D = 0.  QED!
""")

# Verify
print("Symbolic verification for a triple:")
a, b, c = sp.symbols('a b c')
expr = 1/((a-b)*(a-c)) + 1/((b-a)*(b-c)) + 1/((c-a)*(c-b))
print(f"  1/((a-b)(a-c)) + 1/((b-a)(b-c)) + 1/((c-a)(c-b)) = {sp.simplify(expr)}")

print("""
THEOREM (Proved): Phi_n = 2 * Sm2 = 2 * sum_{i<j} 1/(lambda_i - lambda_j)^2

PROOF: Expand Phi_n = sum_i H_i^2 = sum_i (sum_{j!=i} 1/d_{ij})^2.
The diagonal terms give 2*Sm2. The cross terms, summed over all
unordered triples {i,j,k}, give 0 because for each triple:
  1/((a-b)(a-c)) + 1/((b-a)(b-c)) + 1/((c-a)(c-b)) = 0.
QED.
""")

# =============================================================
# Part 3: Reformulated conjecture
# =============================================================
print("\n--- Part 3: Reformulated conjecture ---\n")

print("""
REFORMULATED CONJECTURE:
  1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q)

Since Phi_n = 2*Sm2, this becomes:
  1/(2*Sm2(r)) >= 1/(2*Sm2(p)) + 1/(2*Sm2(q))

i.e., 1/Sm2(r) >= 1/Sm2(p) + 1/Sm2(q)

where Sm2(p) = sum_{i<j} 1/(mu_i - mu_j)^2 for roots mu of p.

This is EQUIVALENT to:
  Sm2(r) <= Sm2(p) * Sm2(q) / (Sm2(p) + Sm2(q))

i.e., Sm2(r) <= harmonic_mean(Sm2(p), Sm2(q)) / 2

Or equivalently: Sm2 is SUBHARMONIC under MSS convolution.

The quantity Sm2 = sum_{i<j} 1/(lambda_i - lambda_j)^2 is the
TOTAL INVERSE-SQUARE INTERACTION ENERGY between the roots.

In electrostatics: if charges are at lambda_1,...,lambda_n and the
pair potential is V(r) = 1/r^2 (inverse-square, not Coulomb),
then Sm2 is the total potential energy.
""")

# =============================================================
# Part 4: Sm2 in terms of cumulants
# =============================================================
print("\n--- Part 4: Sm2 in terms of power sums ---\n")

# For n=3, gaps s=d1, t=d2:
# Sm2 = 1/s^2 + 1/t^2 + 1/(s+t)^2
# We showed Phi_3 = 2(s^2+st+t^2)^2 / (s^2*t^2*(s+t)^2)
# And 2*Sm2 = 2/s^2 + 2/t^2 + 2/(s+t)^2 = 2[(s+t)^2*t^2 + (s+t)^2*s^2 + s^2*t^2] / (s^2*t^2*(s+t)^2)
# = 2[t^2(s+t)^2 + s^2(s+t)^2 + s^2*t^2] / denom
# = 2[(s^2+t^2)(s+t)^2 + s^2*t^2] / denom
# = 2[s^4+2s^3t+s^2t^2+s^2t^2+2st^3+t^4+s^2t^2] / denom
# = 2[s^4+2s^3t+3s^2t^2+2st^3+t^4] / denom
# = 2(s^2+st+t^2)^2 / denom  -> matches!

s, t = sp.symbols('s t', positive=True)
sm2_3 = 1/s**2 + 1/t**2 + 1/(s+t)**2
sm2_3_expanded = sp.simplify(sm2_3)
phi_3_expected = (s**2+s*t+t**2)**2 / (s**2*t**2*(s+t)**2)
print(f"n=3: Sm2 = {sm2_3_expanded}")
print(f"      = (s^2+st+t^2)^2 / (s^2*t^2*(s+t)^2)? {sp.simplify(sm2_3 - phi_3_expected) == 0}")

# 1/Sm2 for n=3:
inv_sm2_3 = sp.simplify(1/sm2_3)
print(f"1/Sm2 = {inv_sm2_3}")

# =============================================================
# Part 5: Direct proof attempt for n=3
# =============================================================
print("\n\n--- Part 5: n=3 proof attempt ---\n")

print("""
For n=3, the conjecture becomes:
  1/Sm2(r) >= 1/Sm2(p) + 1/Sm2(q)

where Sm2(p) = 1/s_p^2 + 1/t_p^2 + 1/(s_p+t_p)^2
      Sm2(q) = 1/s_q^2 + 1/t_q^2 + 1/(s_q+t_q)^2
and s_r, t_r are the gaps of r = p boxplus_3 q.

The key challenge: how do s_r, t_r relate to (s_p, t_p, s_q, t_q)?

For n=3: p has 3 roots, so p is determined by (mean, s_p, t_p).
By translation invariance of 1/Phi, we can center at 0.

p boxplus_3 q: need the full MSS formula.
The result r depends not just on (s_p, t_p, s_q, t_q) but also on
the relative positions of the roots (the means).

Actually, by translation invariance of boxplus_n:
(p(x-a)) boxplus_n (q(x-b)) = (p boxplus_n q)(x - (a+b))
[assuming mean shifts additively, which follows from kappa_1 additivity]

So WLOG mean(p) = mean(q) = 0, and then mean(r) = 0.

With mean 0: p is determined by (s_p, t_p) and q by (s_q, t_q).
""")

# For n=3, centered polynomials:
# p(x) = (x - a_p)(x - b_p)(x - c_p) with a+b+c = 0
# Roots: a = -(2s+t)/3, b = (s-t)/3, c = (s+2t)/3
# p(x) = x^3 + a_2 x + a_3

# The MSS formula for n=3:
# c_0 = 1
# c_1 = a_1 + b_1 (both 0 since centered)
# c_2 = (2!*1!)/(3!*1!) * (a_1*b_1) + a_2 + b_2 = (1/3)*a_1*b_1 + a_2 + b_2
# With a_1 = b_1 = 0: c_2 = a_2 + b_2
# c_3 = ... complex

# Actually for centered polynomials:
# p(x) = x^3 + e_2_p x + e_3_p where
# e_2_p = lambda_1*lambda_2 + lambda_1*lambda_3 + lambda_2*lambda_3 = -(s_p^2 + s_p*t_p + t_p^2)/3
# (using the parametrization above)
# e_3_p = -lambda_1*lambda_2*lambda_3

# Let me compute numerically
print("Computing MSS for n=3 centered polynomials:")
def centered_roots(s, t):
    """Roots of centered (mean 0) polynomial with gaps s, t."""
    return np.array([-(2*s+t)/3, (s-t)/3, (s+2*t)/3])

for sp_val, tp_val, sq_val, tq_val in [(2, 2, 1, 1), (3, 1, 2, 2), (1, 4, 2, 1)]:
    p_roots = centered_roots(sp_val, tp_val)
    q_roots = centered_roots(sq_val, tq_val)
    r_roots = mss_convolve_roots(p_roots, q_roots)
    r_gaps = np.diff(r_roots)

    sm2_p = Sm2(p_roots)
    sm2_q = Sm2(q_roots)
    sm2_r = Sm2(r_roots)

    inv_gap = 1/sm2_r - 1/sm2_p - 1/sm2_q

    print(f"  p: gaps=({sp_val},{tp_val}), q: gaps=({sq_val},{tq_val})")
    print(f"  r: gaps=({r_gaps[0]:.4f},{r_gaps[1]:.4f})")
    print(f"  Sm2: p={sm2_p:.6f}, q={sm2_q:.6f}, r={sm2_r:.6f}")
    print(f"  1/Sm2: p={1/sm2_p:.6f}, q={1/sm2_q:.6f}, r={1/sm2_r:.6f}")
    print(f"  gap = 1/Sm2(r) - 1/Sm2(p) - 1/Sm2(q) = {inv_gap:.6e}")
    print()

# =============================================================
# Part 6: Connection to Cauchy-Schwarz via Sm2 * S2
# =============================================================
print("\n--- Part 6: Cauchy-Schwarz bound Sm2 * S2 >= n(n-1)/2 ---\n")

print("""
By Cauchy-Schwarz on the pairwise distances:
  (sum_{i<j} 1/(d_ij)^2) * (sum_{i<j} d_ij^2) >= (sum_{i<j} 1)^2 = [n(n-1)/2]^2

So Sm2 * S2 >= n^2(n-1)^2/4.

Equality iff all d_ij are equal, which requires equally spaced roots.
(But for n >= 3, equally spaced roots have d_ij = k*(lambda_{n}-lambda_1)/(n-1)
for various k, so equality doesn't hold for n >= 3.)

Actually, Cauchy-Schwarz gives:
  (sum_{i<j} a_{ij}^2)(sum_{i<j} b_{ij}^2) >= (sum_{i<j} a_{ij} b_{ij})^2

With a_{ij} = 1/|d_{ij}|, b_{ij} = |d_{ij}|:
  Sm2 * S2 >= (sum_{i<j} 1)^2 = [n(n-1)/2]^2

So Sm2 >= [n(n-1)/2]^2 / S2.

Since S2 is additive: S2(r) = S2(p) + S2(q), this gives:
  Sm2(r) >= [n(n-1)/2]^2 / (S2(p) + S2(q))

And: 1/Sm2(r) <= (S2(p) + S2(q)) / [n(n-1)/2]^2

Similarly: 1/Sm2(p) <= S2(p) / [n(n-1)/2]^2
          1/Sm2(q) <= S2(q) / [n(n-1)/2]^2

But adding: 1/Sm2(p) + 1/Sm2(q) <= (S2(p)+S2(q)) / [n(n-1)/2]^2

So BOTH 1/Sm2(r) and 1/Sm2(p) + 1/Sm2(q) are bounded above by the
same quantity! This means the inequality could go either way.
The Cauchy-Schwarz bound is too crude.
""")

np.random.seed(42)
print("Verifying Sm2 * S2 >= [n(n-1)/2]^2:")
for n in [3, 4, 5]:
    min_product = float('inf')
    bound = (n*(n-1)/2)**2
    for trial in range(1000):
        roots = np.sort(np.random.randn(n) * 3)
        while np.min(np.diff(roots)) < 0.1:
            roots = np.sort(np.random.randn(n) * 3)
        product = Sm2(roots) * S2(roots)
        min_product = min(min_product, product)
    print(f"  n={n}: min(Sm2*S2)={min_product:.4f}, bound={bound:.4f}, ratio={min_product/bound:.4f}")

# =============================================================
# Part 7: The CONVEXITY of 1/Sm2
# =============================================================
print("\n\n--- Part 7: Convexity properties ---\n")

print("""
KEY REFORMULATION (using Phi = 2*Sm2):

1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q)
<=> 1/Sm2(r) >= 1/Sm2(p) + 1/Sm2(q)

Sm2 = sum_{i<j} 1/(lambda_i - lambda_j)^2 = sum_{i<j} f(lambda_i - lambda_j)
where f(x) = 1/x^2.

This is a function of PAIRWISE differences. Under MSS convolution,
the pairwise differences transform in a specific way.

The function 1/Sm2 is CONCAVE in what?

Note: For n=2, Sm2 = 1/(lambda_1-lambda_2)^2, and
1/Sm2 = (lambda_1-lambda_2)^2 = S2.
Since S2 is additive, 1/Sm2 = S2 is additive: 1/Sm2(r) = 1/Sm2(p) + 1/Sm2(q).

For n >= 3, 1/Sm2 is NOT linear in S2 (it depends on the distribution
of gaps, not just their total square sum). But it is SUPER-additive.

Insight: 1/Sm2 as a function of the gap configuration.
Let d = (d_1, ..., d_{n-1}) be the sorted gaps with d_i > 0.
Then lambda_i = lambda_1 + sum_{k=1}^{i-1} d_k, and
  d_{ij} = lambda_i - lambda_j = sum_{k=j}^{i-1} d_k for i > j.

Sm2 = sum_{i>j} 1/(sum_{k=j}^{i-1} d_k)^2

This is a function of the gap vector d. We've shown it's Schur-convex
(more equal gaps = smaller Sm2).
""")

# =============================================================
# Part 8: Titu's lemma / Engel form
# =============================================================
print("\n--- Part 8: Titu's lemma approach ---\n")

print("""
TITU'S LEMMA: sum_i x_i^2/y_i >= (sum x_i)^2 / (sum y_i)
for y_i > 0.

Applied to our problem:
Let x_i = sqrt(S2_i), y_i = Sm2_i * S2_i = Q_i.

Then: sum_i S2_i / Q_i >= (sum sqrt(S2_i))^2 / (sum Q_i)

But we need (sum S2_i) / Q_r >= sum S2_i/Q_i.
This would follow if Q_r <= sum Q_i ... but Q_r <= Q_p + Q_q is too weak.

Actually: we showed Q_r <= max(Q_p, Q_q) numerically!
And the conjecture in Q-form is:
  (S2_p + S2_q) / Q_r >= S2_p/Q_p + S2_q/Q_q

If Q_r <= min(Q_p, Q_q), then LHS >= (S2_p+S2_q)/min(Q) >= S2_p/Q_p + S2_q/Q_q?
No, that's not right either.

Let me try yet another angle.
""")

# Actually, let me verify if Q_r <= min(Q_p, Q_q) (STRONGER than max)
print("Is Q_r <= min(Q_p, Q_q)?")
np.random.seed(42)
violations_min = 0
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
        Qp = Sm2(p_roots) * S2(p_roots)
        Qq = Sm2(q_roots) * S2(q_roots)
        Qr = Sm2(r_roots) * S2(r_roots)
        if Qr > min(Qp, Qq) + 1e-8:
            violations_min += 1
    except:
        continue
print(f"  {violations_min} violations of Q_r <= min(Q_p, Q_q)")

# =============================================================
# Part 9: The CORRECT sufficient condition
# =============================================================
print("\n\n--- Part 9: Sufficient condition ---\n")

print("""
The conjecture: 1/Sm2(r) >= 1/Sm2(p) + 1/Sm2(q)
where S2(r) = S2(p) + S2(q).

Let alpha = S2(p), beta = S2(q), so S2(r) = alpha + beta.

Write Sm2 = Sm2(S2, shape) where "shape" captures the distribution
of gaps relative to their total S2.

For a fixed S2, define the SHAPE FACTOR:
  Q(d) = Sm2(d) * S2(d)  (dimensionless, scale-invariant)

Then Sm2 = Q / S2, and 1/Sm2 = S2/Q.

The conjecture becomes:
  (alpha+beta) / Q_r >= alpha/Q_p + beta/Q_q

i.e., (alpha+beta) * Q_p * Q_q >= alpha * Q_q * Q_r + beta * Q_p * Q_r
i.e., (alpha+beta) * Q_p * Q_q >= Q_r * (alpha * Q_q + beta * Q_p)
i.e., Q_r <= Q_p * Q_q * (alpha+beta) / (alpha * Q_q + beta * Q_p)

By AM-GM: alpha*Q_q + beta*Q_p >= 2*sqrt(alpha*beta*Q_p*Q_q)
So: RHS <= Q_p*Q_q*(alpha+beta) / (2*sqrt(alpha*beta*Q_p*Q_q))
         = sqrt(Q_p*Q_q) * (alpha+beta) / (2*sqrt(alpha*beta))
         = sqrt(Q_p*Q_q) * (sqrt(alpha/beta) + sqrt(beta/alpha)) / 2

This is >= sqrt(Q_p*Q_q).

So a SUFFICIENT condition is: Q_r <= sqrt(Q_p * Q_q)
(geometric mean bound on Q).
""")

# Check: Q_r <= sqrt(Q_p * Q_q)?
print("Is Q_r <= sqrt(Q_p * Q_q)?")
np.random.seed(42)
violations_geom = 0
min_ratio_geom = float('inf')
for trial in range(10000):
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
        Qp = Sm2(p_roots) * S2(p_roots)
        Qq = Sm2(q_roots) * S2(q_roots)
        Qr = Sm2(r_roots) * S2(r_roots)
        geom = np.sqrt(Qp * Qq)
        ratio = Qr / geom
        min_ratio_geom = min(min_ratio_geom, ratio)
        if Qr > geom + 1e-8:
            violations_geom += 1
    except:
        continue
print(f"  {violations_geom} violations, min(Q_r/sqrt(Q_p*Q_q)) = {min_ratio_geom:.6f}")

# Check: exact sufficient condition Q_r <= Q_p*Q_q*(alpha+beta)/(alpha*Q_q+beta*Q_p)
print("\nVerifying EXACT sufficient condition:")
violations_exact = 0
for trial in range(10000):
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
        Qp = Sm2(p_roots) * S2(p_roots)
        Qq = Sm2(q_roots) * S2(q_roots)
        Qr = Sm2(r_roots) * S2(r_roots)
        alpha = S2(p_roots)
        beta = S2(q_roots)
        bound = Qp * Qq * (alpha+beta) / (alpha * Qq + beta * Qp)
        if Qr > bound + 1e-8:
            violations_exact += 1
    except:
        continue
print(f"  {violations_exact} violations of Q_r <= Q_p*Q_q*(S2_p+S2_q)/(S2_p*Q_q+S2_q*Q_p)")
print("  (This is EQUIVALENT to the original conjecture.)")

print("\nDone.")
