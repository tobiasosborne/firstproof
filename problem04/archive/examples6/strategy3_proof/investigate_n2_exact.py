"""
Exact n=2 investigation. For n=2 we can compute everything analytically.

For degree 2: p(x) = (x - a)(x - b), q(x) = (x - c)(x - d).
The finite free additive convolution r = p boxplus_2 q has the property:
  e_1(r) = e_1(p) + e_1(q) = (a+b) + (c+d)
  kappa_2(p) = e_1(p)^2/2 - e_2(p) = (a+b)^2/2 - ab = (a-b)^2/2
  (Actually for n=2, the finite free cumulant kappa_2^{(2)} = e_2/C(2,2) + ... let me compute directly)

Actually for n=2 the MSS convolution is:
  p boxplus_2 q: roots are the eigenvalues of [[a, t], [t, d]] + [[c, s], [s, b]]
  Wait, that's not right. Let me use the coefficient formula.

For n=2: p(x) = x^2 - (a+b)x + ab, q(x) = x^2 - (c+d)x + cd.
The MSS convolution is defined by:
  (p boxplus_n q)(x) = sum_{k=0}^{n} (-1)^k * sum_{j=0}^{k} C(k,j)*C(n-k,n-k)/C(n,k) ...

Let me just use the direct definition: for n=2,
  p boxplus_2 q has coefficients (in x^2 + a_1 x + a_0 form):
  a_1(r) = a_1(p) + a_1(q)  [sum of means is additive]

  For the constant term, the finite free cumulant satisfies additivity:
  kappa_2^{(n)}(p boxplus_n q) = kappa_2^{(n)}(p) + kappa_2^{(n)}(q)

  For n=2: kappa_2^{(2)}(p) = -a_0(p) + a_1(p)^2/2 = -e_2 + e_1^2/2 = (a-b)^2/2
  Wait, let me check: for p(x) = x^2 - e_1 x + e_2 where e_1 = a+b, e_2 = ab,
  the variance-like quantity is e_1^2/(n*(n-1)) - e_2/(n-1) for the empirical measure...
  Actually let me just look at this concretely.

For n=2, p boxplus_2 q is:
  r(x) = x^2 - (a+b+c+d)x + [ab + cd + (a+b)(c+d)/1 - ... ]

  Hmm, let me just compute it numerically very accurately by using the matrix model.
  For n=2, r = p boxplus_2 q can be computed EXACTLY as:
  The expected char poly of [[a,0],[0,b]] + U[[c,0],[0,d]]U* averaged over U in SU(2).

  Actually, for 2x2 this has a nice closed form. The eigenvalues of
  A + UBU* where A = diag(a,b), B = diag(c,d), U in SU(2) are:
  lambda_pm = (a+b+c+d)/2 +/- sqrt(((a-b)/2)^2 + ((c-d)/2)^2 + (a-b)(c-d)/2 * cos(2*theta))
  where theta parameterizes U.

  The expected char poly averages over theta in [0, 2pi]:
  E[lambda_1 + lambda_2] = a+b+c+d (always)
  E[lambda_1 * lambda_2] = E[det(A + UBU*)] = ...

  For 2x2: det(A + UBU*) = det(A) + det(B) + Tr(A)*Tr(B) - Tr(AB') where B' = UBU*
  Wait: det(M) = (Tr M)^2/2 - Tr(M^2)/2 for 2x2.

  Actually: det(A + UBU*) = det(A) + det(B) + tr(adj(A)*UBU*) where adj is the adjugate...

  Let me just use: for 2x2 Hermitian matrices,
  det(A + UBU*) = (a+b+c+d)^2/2 - ((a+b+c+d)^2 - ... this is getting messy.

  SIMPLEST: p boxplus_2 q for p(x) = (x-a)(x-b), q(x) = (x-c)(x-d) is:
  r(x) = x^2 - (a+b+c+d)x + [ab + cd + (a+b)(c+d)/2 - (a^2+b^2+c^2+d^2)/2 + (a+b+c+d)^2/4 - ...]

  OK I'll just compute directly.
"""

import numpy as np

def boxplus_2_exact(a, b, c, d):
    """Compute roots of (x-a)(x-b) boxplus_2 (x-c)(x-d) exactly.

    For n=2, the finite free convolution is:
    - Mean is additive: mean(r) = mean(p) + mean(q)
    - Variance is additive: var(r) = var(p) + var(q) where var = (root_1 - root_2)^2 / 4

    Actually for the finite free convolution with n=2:
    The finite free cumulants kappa_k^{(n)} satisfy:
    kappa_1^{(2)} = e_1/2 (half the sum of roots = mean)
    kappa_2^{(2)} relates to the variance.

    For the empirical spectral distribution of a 2x2 matrix with eigenvalues a,b:
    mu = (delta_a + delta_b)/2, so mean = (a+b)/2, variance = (a-b)^2/4.

    The finite free R-transform for n=2:
    R^{(2)}_p(z) = kappa_1^{(2)} + kappa_2^{(2)} * z
    where kappa_1^{(2)} = (a+b)/2, kappa_2^{(2)} = (a-b)^2/4 (I think).

    Then for r = p boxplus_2 q:
    kappa_k^{(2)}(r) = kappa_k^{(2)}(p) + kappa_k^{(2)}(q)

    So mean(r) = (a+b)/2 + (c+d)/2 = (a+b+c+d)/2
    var(r) = (a-b)^2/4 + (c-d)^2/4

    Roots of r: mean +/- sqrt(var) = (a+b+c+d)/2 +/- sqrt((a-b)^2/4 + (c-d)^2/4)

    Let me verify with Monte Carlo.
    """
    mean_r = (a + b + c + d) / 2
    var_r = (a - b)**2 / 4 + (c - d)**2 / 4

    nu_1 = mean_r - np.sqrt(var_r)
    nu_2 = mean_r + np.sqrt(var_r)
    return nu_1, nu_2


def boxplus_2_montecarlo(a, b, c, d, samples=50000):
    """Monte Carlo verification."""
    A = np.diag([a, b])
    B = np.diag([c, d])

    sum_e1 = 0.0
    sum_e2 = 0.0

    for _ in range(samples):
        Z = (np.random.randn(2, 2) + 1j * np.random.randn(2, 2)) / np.sqrt(2)
        Q, R = np.linalg.qr(Z)
        dd = np.diagonal(R)
        ph = dd / np.abs(dd)
        U = Q * ph[np.newaxis, :]

        M = A + U @ B @ U.conj().T
        eigs = np.sort(np.real(np.linalg.eigvalsh(M)))
        sum_e1 += eigs[0] + eigs[1]
        sum_e2 += eigs[0] * eigs[1]

    e1 = sum_e1 / samples
    e2 = sum_e2 / samples

    # roots of x^2 - e1*x + e2
    disc = e1**2 - 4*e2
    if disc < 0:
        disc = 0
    nu_1 = (e1 - np.sqrt(disc)) / 2
    nu_2 = (e1 + np.sqrt(disc)) / 2
    return nu_1, nu_2


# Test the exact formula
print("Verifying exact n=2 formula against Monte Carlo:")
np.random.seed(42)
for _ in range(5):
    a, b = sorted(np.random.randn(2) * 2)
    c, d = sorted(np.random.randn(2) * 2)

    exact = boxplus_2_exact(a, b, c, d)
    mc = boxplus_2_montecarlo(a, b, c, d)

    print(f"  p: ({a:.4f}, {b:.4f}), q: ({c:.4f}, {d:.4f})")
    print(f"  Exact: ({exact[0]:.6f}, {exact[1]:.6f})")
    print(f"  MC:    ({mc[0]:.6f}, {mc[1]:.6f})")
    print(f"  Error: {abs(exact[0]-mc[0]):.6f}, {abs(exact[1]-mc[1]):.6f}")
    print()


# Now do exact n=2 analysis
print("\n" + "="*60)
print("EXACT n=2 ANALYSIS")
print("="*60)

def full_n2_analysis(a, b, c, d, label=""):
    """Complete exact analysis for n=2."""
    assert a < b and c < d

    nu_1, nu_2 = boxplus_2_exact(a, b, c, d)

    # H values
    H_p_a = 1/(a - b)
    H_p_b = 1/(b - a)
    H_q_c = 1/(c - d)
    H_q_d = 1/(d - c)
    H_r_1 = 1/(nu_1 - nu_2)
    H_r_2 = 1/(nu_2 - nu_1)

    Phi_p = H_p_a**2 + H_p_b**2  # = 2/(b-a)^2
    Phi_q = H_q_c**2 + H_q_d**2  # = 2/(d-c)^2
    Phi_r = H_r_1**2 + H_r_2**2  # = 2/(nu_2-nu_1)^2

    print(f"\n--- {label} ---")
    print(f"p roots: a={a}, b={b}, gap={b-a}")
    print(f"q roots: c={c}, d={d}, gap={d-c}")
    print(f"r roots: nu_1={nu_1:.8f}, nu_2={nu_2:.8f}, gap={nu_2-nu_1:.8f}")

    # For n=2: Phi_n = 2/gap^2
    gap_p = b - a
    gap_q = d - c
    gap_r = nu_2 - nu_1

    print(f"Phi_p = 2/{gap_p}^2 = {Phi_p:.8f}")
    print(f"Phi_q = 2/{gap_q}^2 = {Phi_q:.8f}")
    print(f"Phi_r = 2/{gap_r}^2 = {Phi_r:.8f}")

    # The inequality: 1/Phi_r >= 1/Phi_p + 1/Phi_q
    # = gap_r^2/2 >= gap_p^2/2 + gap_q^2/2
    # = gap_r^2 >= gap_p^2 + gap_q^2

    # gap_r = 2*sqrt(var_r) = 2*sqrt((b-a)^2/4 + (d-c)^2/4) = sqrt((b-a)^2 + (d-c)^2)
    # So gap_r^2 = (b-a)^2 + (d-c)^2 = gap_p^2 + gap_q^2
    # EQUALITY!

    print(f"\ngap_r^2 = {gap_r**2:.10f}")
    print(f"gap_p^2 + gap_q^2 = {gap_p**2 + gap_q**2:.10f}")
    print(f"Difference: {gap_r**2 - (gap_p**2 + gap_q**2):.2e}")

    inv_Phi_r = 1.0/Phi_r
    inv_sum = 1.0/Phi_p + 1.0/Phi_q
    print(f"\n1/Phi_r = {inv_Phi_r:.10f}")
    print(f"1/Phi_p + 1/Phi_q = {inv_sum:.10f}")
    print(f"Difference: {inv_Phi_r - inv_sum:.2e}")
    print(f"EQUALITY for n=2!")

    # Now check the AB formulation
    A = Phi_p - Phi_r
    B = Phi_q - Phi_r
    AB = A * B
    h4 = Phi_r**2

    print(f"\nA = Phi_p - Phi_r = {A:.10f}")
    print(f"B = Phi_q - Phi_r = {B:.10f}")
    print(f"AB = {AB:.10f}")
    print(f"h^4 = Phi_r^2 = {h4:.10f}")
    print(f"AB - h^4 = {AB - h4:.2e}")

    # For n=2: Phi_p = 2/s^2, Phi_q = 2/t^2, Phi_r = 2/(s^2+t^2)
    # where s = gap_p, t = gap_q
    # A = 2/s^2 - 2/(s^2+t^2) = 2*t^2/(s^2*(s^2+t^2))
    # B = 2/t^2 - 2/(s^2+t^2) = 2*s^2/(t^2*(s^2+t^2))
    # AB = 4/(s^2+t^2)^2 = Phi_r^2/1 = h4
    # So AB = h4 EXACTLY for n=2. Equality!

    s, t = gap_p, gap_q
    A_exact = 2*t**2/(s**2*(s**2+t**2))
    B_exact = 2*s**2/(t**2*(s**2+t**2))
    AB_exact = 4/(s**2+t**2)**2
    h4_exact = (2/(s**2+t**2))**2

    print(f"\nExact: A = {A_exact:.10f}, B = {B_exact:.10f}")
    print(f"Exact: AB = {AB_exact:.10f}, h4 = {h4_exact:.10f}")
    print(f"AB = h4: {np.isclose(AB_exact, h4_exact)}")

    return {'Phi_p': Phi_p, 'Phi_q': Phi_q, 'Phi_r': Phi_r}


# Test cases
full_n2_analysis(-1, 1, -2, 2, "n=2 symmetric")
full_n2_analysis(0, 1, 0, 3, "n=2 asymmetric")
full_n2_analysis(-5, 5, -1, 1, "n=2 extreme ratio")

print("\n\n" + "="*60)
print("KEY CONCLUSION FOR n=2:")
print("The inequality 1/Phi_r >= 1/Phi_p + 1/Phi_q is EQUALITY at n=2.")
print("This explains why Monte Carlo got it wrong -- the MC was inaccurate.")
print("="*60)

# Now verify that the earlier Monte Carlo test was wrong
print("\n\nVerifying Monte Carlo was inaccurate for n=2:")
nu1_exact, nu2_exact = boxplus_2_exact(-1, 1, -2, 2)
print(f"Exact roots: {nu1_exact:.8f}, {nu2_exact:.8f}")
print(f"Exact gap: {nu2_exact - nu1_exact:.8f}")
print(f"sqrt(4+16) = sqrt(20) = {np.sqrt(20):.8f}")

# The MC gave roots -2.2334659, 2.2334659 => gap = 4.4669318
# Exact gives gap = sqrt((2)^2 + (4)^2) = sqrt(20) = 4.472136...
print(f"MC gap was: 4.4669318, exact gap: {np.sqrt(20):.8f}")
print(f"MC error in gap: {abs(4.4669318 - np.sqrt(20)):.6f}")
print(f"This MC error in gap propagates to ~{2*abs(4.4669318 - np.sqrt(20))/np.sqrt(20)*100:.1f}% error in Phi_r")
