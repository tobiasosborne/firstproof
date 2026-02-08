"""
Debug the boxplus implementation by checking n=2 case analytically.
"""
import numpy as np
from math import factorial, comb

def boxplus_n2(a, b, c, d):
    """
    For n=2, p(x) = (x-a)(x-b), q(x) = (x-c)(x-d).

    Standard coefficients: p(x) = x^2 - (a+b)x + ab
    Normalized: p(x) = x^2 - C(2,1)*a_1*x + C(2,2)*a_2 = x^2 - 2*a_1*x + a_2
    So a_1 = (a+b)/2, a_2 = ab.
    Similarly b_1 = (c+d)/2, b_2 = cd.

    MSS formula: c_k = sum_{i+j=k} [(2-i)!(2-j)!/(2!(2-k)!)] * a_i * b_j

    c_0 = a_0*b_0 * [2!*2!/(2!*2!)] = 1*1*1 = 1
    c_1 = a_1*b_0 * [(2-1)!(2-0)!/(2!(2-1)!)] + a_0*b_1 * [(2-0)!(2-1)!/(2!(2-1)!)]
         = a_1*1 * [1!*2!/(2!*1!)] + 1*b_1 * [2!*1!/(2!*1!)]
         = a_1 * [2/(2)] + b_1 * [2/(2)]
         = a_1 + b_1
    c_2 = a_2*b_0 * [(2-2)!(2-0)!/(2!(2-2)!)] + a_1*b_1 * [(2-1)!(2-1)!/(2!(2-2)!)]
          + a_0*b_2 * [(2-0)!(2-2)!/(2!(2-2)!)]
        = a_2*1 * [0!*2!/(2!*0!)] + a_1*b_1 * [1!*1!/(2!*0!)] + 1*b_2 * [2!*0!/(2!*0!)]
        = a_2 * [2/2] + a_1*b_1 * [1/2] + b_2 * [2/2]
        = a_2 + (1/2)*a_1*b_1 + b_2

    So r(x) = x^2 - 2*c_1*x + c_2
    = x^2 - 2*(a_1+b_1)*x + (a_2 + (1/2)*a_1*b_1 + b_2)
    = x^2 - (a+b+c+d)*x + (ab + (1/2)*(a+b)/2*(c+d)/2 + cd)
    = x^2 - (a+b+c+d)*x + (ab + (a+b)(c+d)/8 + cd)

    Roots via quadratic formula:
    r_1, r_2 = ((a+b+c+d) +/- sqrt((a+b+c+d)^2 - 4*(ab + (a+b)(c+d)/8 + cd))) / 2

    Discriminant = (a+b+c+d)^2 - 4*(ab + (a+b)(c+d)/8 + cd)
    = (a+b)^2 + 2(a+b)(c+d) + (c+d)^2 - 4ab - (a+b)(c+d)/2 - 4cd
    = (a-b)^2 + (3/2)(a+b)(c+d) + (c-d)^2
    Hmm, that doesn't simplify to (a-b)^2 + (c-d)^2 which we'd need for Pythagorean.

    Wait, let me recheck. 1/Phi_2(p) = (b-a)^2/2 for n=2.
    The conjecture says 1/Phi_2(r) >= 1/Phi_2(p) + 1/Phi_2(q)
    i.e., (r_2-r_1)^2/2 >= (b-a)^2/2 + (d-c)^2/2
    i.e., discriminant of r >= (b-a)^2 + (d-c)^2.

    Discriminant = (a+b+c+d)^2 - 4*(ab + (a+b)(c+d)/8 + cd)

    Let me expand more carefully with s = a+b, t = c+d:
    = s^2 + 2st + t^2 - 4ab - st/2 - 4cd
    = (a-b)^2 + 2ab + 2st + (c-d)^2 + 2cd - 4ab - st/2 - 4cd
    = (a-b)^2 + (c-d)^2 + (3/2)st - 2ab - 2cd

    For (a-b)^2 + (c-d)^2 case we'd need 3st/2 - 2ab - 2cd = 0 for equality.
    With a=0,b=1,c=0,d=1: s=1,t=1, 3/2 - 0 - 0 = 3/2 != 0.

    So the n=2 case does NOT give equality? That contradicts the claim.
    Let me recheck what boxplus should be for n=2...
    """

    n = 2
    # a_k coefficients
    a1 = (a + b) / 2
    a2 = a * b
    b1 = (c + d) / 2
    b2 = c * d

    print(f"a_0=1, a_1={a1}, a_2={a2}")
    print(f"b_0=1, b_1={b1}, b_2={b2}")

    # c_k
    c0 = 1.0
    c1_direct = a1 * (factorial(1)*factorial(2)/(factorial(2)*factorial(1))) \
              + b1 * (factorial(2)*factorial(1)/(factorial(2)*factorial(1)))
    c1 = a1 + b1
    print(f"c_1 = {c1}")

    # c_2 more carefully:
    # i=2,j=0: (2-2)!(2-0)!/(2!(2-2)!) = 0!*2!/(2!*0!) = 1*2/(2*1) = 1
    # i=1,j=1: (2-1)!(2-1)!/(2!(2-2)!) = 1!*1!/(2!*0!) = 1/(2) = 0.5
    # i=0,j=2: (2-0)!(2-2)!/(2!(2-2)!) = 2!*0!/(2!*0!) = 1
    c2 = a2 * 1 + a1*b1 * 0.5 + b2 * 1
    print(f"c_2 = {c2}")
    print(f"  = a2 + a1*b1/2 + b2 = {a2} + {a1*b1/2} + {b2}")

    # r(x) = x^2 - 2*c1*x + c2
    poly_r = [1, -2*c1, c2]
    roots_r = np.roots(poly_r)
    roots_r = np.sort(np.real(roots_r))
    print(f"r(x) = x^2 - {2*c1}x + {c2}")
    print(f"roots of r: {roots_r}")

    gap_r = roots_r[1] - roots_r[0]
    gap_p = abs(b - a)
    gap_q = abs(d - c)

    print(f"gap_p = {gap_p}, gap_q = {gap_q}, gap_r = {gap_r}")
    print(f"gap_r^2 = {gap_r**2}")
    print(f"gap_p^2 + gap_q^2 = {gap_p**2 + gap_q**2}")

    disc = 4*c1**2 - 4*c2
    print(f"discriminant = {disc}")
    print(f"(a-b)^2 + (c-d)^2 = {(a-b)**2 + (c-d)**2}")

    # So gap_r^2 = discriminant = 4*c1^2 - 4*c2
    # = 4*(a1+b1)^2 - 4*(a2 + a1*b1/2 + b2)
    # = 4*a1^2 + 8*a1*b1 + 4*b1^2 - 4*a2 - 2*a1*b1 - 4*b2
    # = 4*a1^2 - 4*a2 + 6*a1*b1 + 4*b1^2 - 4*b2
    # = (a+b)^2 - 4*ab + 6*(a+b)(c+d)/4 + (c+d)^2 - 4*cd
    # = (a-b)^2 + 3(a+b)(c+d)/2 + (c-d)^2
    # This is > (a-b)^2 + (c-d)^2 when (a+b)(c+d) > 0
    #         < (a-b)^2 + (c-d)^2 when (a+b)(c+d) < 0
    # So the inequality can go EITHER way for n=2!

    # THIS MEANS THE CONJECTURE IS WRONG, or the formula is wrong!

    # Let me check: maybe the MSS formula should have DIFFERENT normalization.
    # Marcus-Spielman-Srivastava 2015 define the finite free convolution differently.
    # Perhaps the correct formula is the "rectangular" or "symmetric" one.

    return roots_r

print("=" * 70)
print("n=2 ANALYTICAL TEST")
print("=" * 70)

# Case 1: centered roots
print("\nCase 1: a=-1, b=1, c=-1, d=1")
boxplus_n2(-1, 1, -1, 1)

print("\nCase 2: a=0, b=1, c=0, d=1")
boxplus_n2(0, 1, 0, 1)

print("\nCase 3: a=-2, b=2, c=-1, d=1")
boxplus_n2(-2, 2, -1, 1)

print("\nCase 4: a=-1, b=1, c=0, d=2")
boxplus_n2(-1, 1, 0, 2)

print("\n" + "=" * 70)
print("CHECKING: What is the CORRECT MSS convolution formula?")
print("=" * 70)

# The MSS 2015 paper: for CHARACTERISTIC polynomials of random matrices
# p_A(x) = det(xI - A), p_B(x) = det(xI - B)
# E[p_{A+UBU*}(x)] = p_A(x) boxplus p_B(x)
# For rank-1 projections this gives the expected characteristic polynomial.

# For n=2, A = diag(a, b), B = diag(c, d), U = random unitary
# A + UBU* has eigenvalues that are the roots of the expected characteristic poly.

# Let me check: what is E[det(xI - A - UBU*)] for 2x2 matrices?
# Using Weingarten calculus or direct integration.

# For 2x2 unitary U = [[cos t, -sin t],[sin t, cos t]] parametrized by t in [0, 2pi):
# A + UBU* = diag(a,b) + U*diag(c,d)*U^{-1}
# Actually UBU* means U @ B @ U^{dagger}
# For U = [[cos t, -sin t],[sin t, cos t]]:
# UBU* = U @ [[c,0],[0,d]] @ U^T (since U is real orthogonal)
# = [[c*cos^2(t)+d*sin^2(t), (c-d)*cos(t)*sin(t)],
#    [(c-d)*cos(t)*sin(t), c*sin^2(t)+d*cos^2(t)]]

# A + UBU* = [[a + c*cos^2+d*sin^2, (c-d)*cos*sin],
#              [(c-d)*cos*sin, b + c*sin^2+d*cos^2]]

# char poly = (x - a - c*cos^2-d*sin^2)(x - b - c*sin^2-d*cos^2) - (c-d)^2*cos^2*sin^2
# = x^2 - (a+b+c+d)*x + [stuff]

# For the E[] w.r.t. Haar measure on O(2) (i.e., uniform t in [0,2pi)):
# E[cos^2(t)] = 1/2, E[sin^2(t)] = 1/2, E[cos^2*sin^2] = 1/8, E[cos^4]=3/8
# Actually for U(2) Haar measure, we need to be more careful.

# For U(2) Haar, parametrize by U = [[e^{ia}*cos t, -e^{ib}*sin t],
#                                      [e^{ic}*sin t, e^{id}*cos t]]
# Actually, the key result is:
# E[p_{A+UBU*}] uses the MIXED discriminant formula.

# Let me just use the known result from Marcus-Spielman-Srivastava (Theorem 4.4):
# For the finite free convolution, the normalized coefficients satisfy:
# c_k = sum_{S in C(n,k)} prod_{i in S} alpha_i * prod_{j in S} beta_j
# Wait, that's the multiplicative version.

# From MSS: the additive version uses the "finite R-transform":
# If p(x) = sum_{k=0}^n (-1)^k e_k(lambda) x^{n-k} and similarly for q,
# then r(x) = sum_{k=0}^n (-1)^k c_k x^{n-k} where
# c_k / C(n,k) = sum_{i+j=k} [C(n,i) * C(n,j)]^{-1} * C(n-i,j)^{-1} * C(n,k) * e_i * f_j / C(n,k)
# Hmm, this is confusing. Let me look up the exact formula.

# From Theorem 4.4 of MSS 2015 "Interlacing families II":
# If p(x) = sum_{k=0}^n a_k x^{n-k} and q(x) = sum_{k=0}^n b_k x^{n-k}
# with a_0 = b_0 = 1, then
# (p boxplus q)(x) = sum_{k=0}^n c_k x^{n-k} where
# c_k = sum_{j=0}^k C(k,j) / C(n,j) * a_j * b_{k-j}
# (with the convention that a_k, b_k are the coefficients of x^{n-k})

# Wait, there are MULTIPLE conventions. Let me use the simplest.

# From Definition 4.1 of MSS 2015:
# p boxplus q (x) = sum_{k=0}^n C(n,k)^{-1} * D^k [p] * D^{n-k} [q] evaluated at x
# where D = d/dx.
# NO, that's not right either.

# Let me try: use the formula from Marcus (2021), "Polynomial convolutions and (finite) free probability"
# or the standard definition via the finite R-transform.

# Standard definition:
# Write p(x) = sum_k (-1)^k e_k x^{n-k} where e_k = e_k(roots_p)
# Write q(x) = sum_k (-1)^k f_k x^{n-k}
# The "finite free cumulants" are r_k(p) = e_k / C(n,k)
# Then r_k(p boxplus q) = r_k(p) + r_k(q)
# i.e., (p boxplus q) has e_k(r) = C(n,k) * (e_k(p)/C(n,k) + e_k(q)/C(n,k))
#      = e_k(p) + e_k(q)
# So (p boxplus q)(x) = x^n + sum_{k=1}^n (-1)^k (e_k(p) + e_k(q)) x^{n-k}

# WAIT — that's WAY too simple. That would mean the convolution is just
# adding the elementary symmetric polynomials. For n=2:
# p(x) = x^2 - (a+b)x + ab
# q(x) = x^2 - (c+d)x + cd
# r(x) = x^2 - (a+b+c+d)x + (ab+cd)

# roots of r: ((a+b+c+d) +/- sqrt((a+b+c+d)^2 - 4(ab+cd))) / 2
# disc = (a+b+c+d)^2 - 4(ab+cd) = (a+b)^2 + 2(a+b)(c+d) + (c+d)^2 - 4ab - 4cd
#       = (a-b)^2 + 2(a+b)(c+d) + (c-d)^2

# This is STILL not Pythagorean. So the conjecture 1/Phi >= 1/Phi + 1/Phi is
# NOT expected to hold with this definition either.

# Actually, wait. The CORRECT finite free convolution from MSS is NOT the
# simple additive one above. It uses a more subtle formula.

# From MSS "Interlacing families II" (Theorem 4.4):
# For p,q of degree n:
# (p boxplus_n q)(x) = D_n^{-1} [ D_n[p] * D_n[q] ]
# where D_n is the "finite difference" operator.
# This is quite different from just adding e_k's.

# Actually, the CORRECT formula from Marcus-Spielman-Srivastava uses:
# "Expected characteristic polynomial" of A + UBU* where U is Haar unitary.

# Let me use the DIRECT matrix approach for n=2:
print("\nDirect matrix computation for n=2:")
a, b = -1.0, 1.0  # eigenvalues of A
c, d = -1.0, 1.0  # eigenvalues of B

A = np.diag([a, b])
B_diag = np.diag([c, d])

# Average over many random unitaries
N_samples = 100000
poly_sum = np.zeros(3)  # coefficients of x^2, x^1, x^0
for _ in range(N_samples):
    # Random 2x2 unitary (Haar measure)
    Z = (np.random.randn(2, 2) + 1j * np.random.randn(2, 2)) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    # Fix phases to get Haar measure
    d_phases = np.diag(R)
    d_phases = d_phases / np.abs(d_phases)
    U = Q @ np.diag(d_phases)

    M = A + U @ B_diag @ U.conj().T
    coeffs = np.poly(M)  # characteristic polynomial coefficients
    poly_sum += np.real(coeffs)

poly_avg = poly_sum / N_samples
print(f"  E[char poly] coefficients: {poly_avg}")
print(f"  E[char poly] roots: {np.sort(np.real(np.roots(poly_avg)))}")

# Compare with MSS formula
print(f"\n  MSS formula (my implementation):")
rr_mss = np.zeros(3)
e_p = [1, a+b, a*b]  # e_0, e_1, e_2 of p
e_q = [1, c+d, c*d]  # e_0, e_1, e_2 of q

# Using the formula from my boxplus:
# a_k = e_k / C(n,k)
a_k = [e_p[k] / comb(2, k) for k in range(3)]
b_k = [e_q[k] / comb(2, k) for k in range(3)]
print(f"  a_k = {a_k}")
print(f"  b_k = {b_k}")

# c_k = sum_{i+j=k} [(n-i)!(n-j)!/(n!(n-k)!)] * a_i * b_j
c_k = [0, 0, 0]
for k in range(3):
    for i in range(k+1):
        j = k - i
        w = factorial(2-i) * factorial(2-j) / (factorial(2) * factorial(2-k))
        c_k[k] += w * a_k[i] * b_k[j]

print(f"  c_k = {c_k}")

# r(x) = x^2 - C(2,1)*c_1*x + C(2,2)*c_2 = x^2 - 2*c_1*x + c_2
poly_mss = [1, -2*c_k[1], c_k[2]]
print(f"  MSS poly: {poly_mss}")
print(f"  MSS roots: {np.sort(np.real(np.roots(poly_mss)))}")

# Simple additive formula:
# r(x) = x^2 - (e_1(p)+e_1(q))x + (e_2(p)+e_2(q))
poly_add = [1, -(e_p[1]+e_q[1]), e_p[2]+e_q[2]]
print(f"\n  Simple additive poly: {poly_add}")
print(f"  Simple additive roots: {np.sort(np.real(np.roots(poly_add)))}")

# The finite free cumulant formula (Marcus 2021):
# r_k(p) = e_k(p) / C(n,k)
# r_k(p boxplus q) = r_k(p) + r_k(q)
# e_k(r) = C(n,k) * (r_k(p) + r_k(q))
r_p = [e_p[k] / comb(2, k) for k in range(3)]
r_q = [e_q[k] / comb(2, k) for k in range(3)]
e_r = [comb(2, k) * (r_p[k] + r_q[k]) for k in range(3)]
poly_cumulant = [1, -e_r[1], e_r[2]]
print(f"\n  Finite free cumulant poly: {poly_cumulant}")
print(f"  Cumulant roots: {np.sort(np.real(np.roots(poly_cumulant)))}")
