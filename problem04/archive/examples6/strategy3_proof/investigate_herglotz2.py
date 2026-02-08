#!/usr/bin/env python3
"""
Part 2: Deeper investigation of omega_1 structure.

Key finding from Part 1: omega_1 maps C^+ to C^+ but NOT with
Im(omega_1(z)) >= Im(z). The function phi = omega_1 - id does NOT
have all-positive residues, resolving the Nevanlinna contradiction.

Wait -- but the Nevanlinna representation theorem says that EVERY
function mapping C^+ to C^+ has the form alpha + beta*z + integral...
with beta >= 0 and positive measure. This should apply to omega_1.

RESOLUTION: omega_1 maps C^+ to C^+ with Im(omega_1(z)) > 0 but
possibly Im(omega_1(z)) < Im(z). In the Nevanlinna form:
  omega_1(z) = a + b*z + sum c_j/(d_j - z), c_j > 0
Then phi = omega_1 - z = a + (b-1)*z + sum c_j/(d_j-z).

If b = 1 (matching the z behavior at infinity), then phi = a + sum c_j/(d_j-z)
and phi' = sum c_j/(d_j-z)^2 > 0 everywhere. Contradicts phi'(nu_k) = 0.

If b < 1, then phi(z) = a + (b-1)*z + sum c_j/(d_j-z), which has
phi'(z) = (b-1) + sum c_j/(d_j-z)^2.
For phi'(nu_k) = 0: sum c_j/(d_j-nu_k)^2 = 1 - b.
This is possible if b < 1 since the left side is positive!

And for large z: omega_1(z) = b*z + ..., so omega_1(z)/z -> b.
If b < 1, then omega_1(z) grows slower than z.

Let's check: what is b = lim_{z->inf} omega_1(z)/z?
From G_r(z) = G_p(omega_1(z)):
For large z: G_r(z) ~ 1/z, G_p(w) ~ 1/w.
So 1/z ~ 1/omega_1(z), hence omega_1(z) ~ z. This gives b = 1!

WAIT. Let me be more careful.
G_r(z) = (1/n) * r'(z)/r(z).
For large z with r(z) = z^n + ...:
r'(z)/r(z) = n/z + ...
So G_r(z) = 1/z + O(1/z^2).

Similarly G_p(w) = 1/w + O(1/w^2).
G_r(z) = G_p(omega_1(z)) means 1/z + O(1/z^2) = 1/omega_1(z) + O(1/omega_1(z)^2).
So omega_1(z) ~ z + O(1) as z -> infinity. Hence b = 1.

So the Nevanlinna representation has b = 1, and phi = a + sum c_j/(d_j - z) with c_j > 0.
And phi' = sum c_j/(d_j-z)^2 > 0 everywhere.
But phi'(nu_k) = 0. Contradiction.

The ONLY resolution: omega_1 does NOT have the standard Nevanlinna representation
because omega_1 is NOT Herglotz in the standard sense.

Let me verify numerically that omega_1 maps C^+ to C^+, but may also map some
points in C^+ to points ON the real line (boundary of C^+). If omega_1 maps some
points of C^+ to R, it is NOT strictly Herglotz but maps C^+ to closure(C^+).

The Nevanlinna representation still holds for maps to closure(C^+), but the
measure can have atoms and the representation can degenerate.

Actually, a FINITE BLASCHKE PRODUCT is an example: B(z) maps C^+ to C^+ but
B maps the unit circle to itself (boundary). A finite Blaschke product of
degree n has exactly n zeros.

For our setting: omega_1 is a DEGREE-n rational function mapping C^+ to C^+_closure.
If omega_1 maps some points of R to R (which it does: omega_1(nu_k) = lambda_k),
then it maps the REAL LINE to itself (generically).

Actually, EVERY rational function with real coefficients maps R to R (or infinity).
So omega_1 maps R to R. The question is whether it maps C^+ to C^+ or C^+_closure.

A real rational function maps C^+ to C^+ iff all its critical points in C^+ have
positive imaginary part. Actually, a real rational function maps C^+ to
either C^+ or C^-, or maps C^+ onto C (depending on the arrangement of roots
and poles on R).

Let me think about this differently.

ACTUAL RESOLUTION: omega_1 has degree n but it is NOT of the form
a + bz + sum c_j/(d_j - z) (which has degree max(1, m+1) with m simple poles).
That form is the general Nevanlinna representation for Herglotz functions of
degree at most m+1.

But omega_1 can be a general degree-n rational function P(z)/Q(z) with
deg P = n, deg Q = n-1 (say), which can map C^+ to C^+_closure without being
of the Nevanlinna form!

WAIT: The Nevanlinna representation theorem says that EVERY holomorphic function
on C^+ with Im(f(z)) >= 0 has that integral representation. For a rational
function, the integral reduces to a finite sum. So omega_1 MUST have that form
IF it maps C^+ to closure(C^+).

The form is: f(z) = alpha + beta*z + sum_{j=1}^m c_j/(d_j - z) with c_j > 0.
This is a function of degree max(1, m+1) (if beta != 0) or m (if beta = 0).
For beta = 1 and m poles: degree is m+1 or m (depending on how you count).

Actually, the Nevanlinna form a + z + sum c_j/(d_j-z) has degree m+1 (where m
is the number of poles). For omega_1 of degree n, we need m = n-1.

Let me reconsider. Write omega_1(z) = (numerator of degree n)/(denominator of degree n-1).
In partial fraction form:
  omega_1(z) = z + b + sum_{j=1}^{n-1} c_j/(z - d_j)

NOTE: I'm using the convention c_j/(z - d_j), not c_j/(d_j - z). Let me use:
  omega_1(z) = z + b + sum_{j=1}^{n-1} c_j/(z - d_j)

For omega_1 to be Herglotz: Im(omega_1(z)) >= 0 when Im(z) > 0.
Im(omega_1(z)) = Im(z) + sum_j c_j * Im(1/(z-d_j))
= Im(z) + sum_j c_j * (-Im(z))/|z-d_j|^2
= Im(z) * (1 - sum_j c_j/|z-d_j|^2)

For this to be >= 0 for all z in C^+:
1 - sum_j c_j/|z-d_j|^2 >= 0 for all z with Im(z) > 0.

This is NOT the same as requiring all c_j > 0!
If some c_j < 0, the sum could still be <= 1.

Let me re-derive the Nevanlinna form. For f(z) Herglotz and rational:
f(z) = a + bz + sum c_j/(z - d_j) with:
- b >= 0
- d_j real
- If b > 0: c_j can have any sign as long as the total imaginary part stays positive
- WAIT no. The standard Nevanlinna representation has:
  f(z) = a + bz + integral [1/(t-z) + t/(1+t^2)] dmu(t) with mu POSITIVE.
  For rational: integral becomes sum c_j * [1/(d_j-z) + d_j/(1+d_j^2)]
  = sum c_j/(d_j-z) + constant
  So f(z) = (new alpha) + bz + sum c_j/(d_j - z) with all c_j >= 0.

Hmm, with the convention 1/(d_j - z):
Im(1/(d_j - z)) = Im(z)/|d_j - z|^2 (positive for z in C^+).
So Im(f(z)) = Im(z)*b + sum c_j * Im(z)/|d_j-z|^2
= Im(z) * (b + sum c_j/|d_j-z|^2)
Since b >= 0 and c_j >= 0, this is >= 0.

With the other convention 1/(z - d_j):
Im(1/(z - d_j)) = -Im(z)/|z - d_j|^2 (NEGATIVE for z in C^+).
So Im(f(z)) = Im(z)*b + sum c_j * (-Im(z))/|z-d_j|^2
= Im(z) * (b - sum c_j/|z-d_j|^2)
For this to be >= 0: b >= sum c_j/|z-d_j|^2 for all z in C^+.

But |z-d_j|^2 = (Re(z)-d_j)^2 + Im(z)^2. As Im(z) -> 0 and Re(z) -> d_j,
|z-d_j|^2 -> 0, so c_j/|z-d_j|^2 -> infinity.
Therefore b >= sum c_j/|z-d_j|^2 is IMPOSSIBLE for any c_j > 0 as z -> d_j.

So the convention matters. With 1/(d_j - z) and c_j > 0, we get Herglotz.
With 1/(z - d_j), the residues must be NEGATIVE (i.e., c_j < 0 in that convention)
to keep Im positive.

Let me use the correct convention:
  omega_1(z) = alpha + z + sum_{j=1}^{n-1} c_j / (d_j - z)  with c_j > 0.

Then:
  phi(z) = omega_1(z) - z = alpha + sum_j c_j/(d_j - z).
  phi'(z) = sum_j c_j/(d_j - z)^2 * (-(-1)) = sum_j c_j/(d_j-z)^2.

Wait, d/dz [c_j/(d_j - z)] = c_j/(d_j-z)^2. Yes, positive for real z != d_j.

So phi'(z) = sum_j c_j/(d_j-z)^2 > 0 for all real z not at a pole.
And phi'(nu_k) = 0 requires sum_j c_j/(d_j-nu_k)^2 = 0 with all terms > 0.
Impossible.

DEFINITIVE CONCLUSION: omega_1 CANNOT be Herglotz (even to closure(C^+))
if omega_1'(nu_k) = 1 and omega_1 has the z + partial_fractions form.

Let me verify numerically. Maybe omega_1 does NOT map C^+ to C^+.
"""

import numpy as np
from math import factorial
from itertools import combinations

def elem_sym_poly(roots, k):
    n = len(roots)
    if k == 0: return 1.0
    if k > n: return 0.0
    return sum(np.prod([roots[i] for i in subset])
               for subset in combinations(range(n), k))

def poly_coeffs_from_roots(roots):
    n = len(roots)
    ek = [elem_sym_poly(roots, k) for k in range(n+1)]
    return [(-1)**k * ek[k] for k in range(n+1)]

def boxplus_mss(roots_p, roots_q):
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
    r_roots = np.sort(np.real(np.roots(c)))
    return r_roots, c

def H_values(roots):
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H


def compute_all_branches(z, roots_p, roots_r):
    """Compute ALL n branches of omega_1 at a complex point z."""
    n = len(roots_p)
    p_poly = np.poly(roots_p)
    r_poly = np.poly(roots_r)
    p_prime = np.polyder(p_poly)
    r_prime = np.polyder(r_poly)

    rz = np.polyval(r_poly, z)
    rpz = np.polyval(r_prime, z)

    # Solve r(z)*p'(w) - r'(z)*p(w) = 0
    p_prime_padded = np.zeros(n+1, dtype=complex)
    p_prime_padded[1:] = p_prime
    combined = rz * p_prime_padded - rpz * np.array(p_poly, dtype=complex)
    roots_w = np.roots(combined)

    return roots_w

np.random.seed(42)

print("="*70)
print("DEFINITIVE TEST: Is omega_1 Herglotz?")
print("="*70)

# For each test case, compute omega_1(z) for z in C^+ and check Im(omega_1(z))
for trial in range(5):
    n = 3
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 1.0:
            roots_p[i] = roots_p[i-1] + 1.0
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 1.0:
            roots_q[i] = roots_q[i-1] + 1.0

    roots_r, c = boxplus_mss(roots_p, roots_q)
    raw = np.roots(c)
    if np.any(np.abs(np.imag(raw)) > 0.01):
        continue
    roots_r = np.sort(np.real(raw))
    if np.any(np.diff(roots_r) < 0.1):
        continue

    print(f"\nTrial {trial}: p={np.round(roots_p,3)}, r={np.round(roots_r,3)}")

    # Test on a grid in C^+
    min_im_ratio = np.inf
    for x in np.linspace(min(roots_r)-3, max(roots_r)+3, 30):
        for y in [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0]:
            z = complex(x, y)
            branches = compute_all_branches(z, roots_p, roots_r)

            # Select the "correct" branch (the one with Im > 0 that connects
            # to the physical sheet). Use the criterion: omega_1(z) should
            # be closest to the identity map for large |z|.
            upper = branches[np.imag(branches) > -y*0.1]
            if len(upper) == 0:
                continue

            # For the physical branch, take the one closest to z
            # (since omega_1(z) ~ z for large z)
            physical = upper[np.argmin(np.abs(upper - z))]

            im_ratio = np.imag(physical) / y
            if im_ratio < min_im_ratio:
                min_im_ratio = im_ratio

    print(f"  min(Im(omega_1)/Im(z)) over grid = {min_im_ratio:.8f}")
    if min_im_ratio < 0:
        print(f"  *** omega_1 MAPS SOME C^+ POINTS TO C^-! ***")
    elif min_im_ratio < 1:
        print(f"  omega_1 maps C^+ to C^+ but Im(omega_1) < Im(z) sometimes")
        print(f"  (phi = omega_1 - id maps some C^+ points below real line)")
    else:
        print(f"  omega_1 is STRICTLY Herglotz: Im(omega_1) >= Im(z)")


# ================================================================
# KEY INSIGHT: omega_1 might be "degree n" in a DIFFERENT sense
# ================================================================
print("\n\n" + "="*70)
print("EXAMINING omega_1 AS A MULTI-VALUED FUNCTION")
print("="*70)

# omega_1 is defined by G_r(z) = G_p(omega_1(z)).
# G_p(w) = (1/n) * p'(w)/p(w) = (1/n) * sum_k 1/(w - lambda_k)
# Setting G_p(w) = g (given):
# (1/n) * sum_k 1/(w - lambda_k) = g
# sum_k 1/(w - lambda_k) = n*g
# This is a degree-n equation in w (after clearing denominators).
# It has n solutions.

# omega_1 is the BRANCH that connects lambda_k to nu_k.
# On the real line between consecutive poles of G_p (which are at lambda_k),
# G_p is monotone decreasing from +inf to -inf.
# Similarly for G_r on the intervals between nu_k.

# As z crosses a pole of G_r (at nu_k), G_r(z) jumps from -inf to +inf.
# The equation G_p(w) = G_r(z) has n solutions in w.
# The "physical" branch omega_1(z) selects the solution that varies continuously.

# Near z = nu_k: G_r(z) has a pole, and G_r(z) sweeps from -inf to +inf.
# For each value of G_r, there's one w in each interval (lambda_k, lambda_{k+1}).
# As z -> nu_k from below: G_r(z) -> -inf, and omega_1(z) -> the w in some interval
# that approaches lambda_k from below.
# As z -> nu_k from above: G_r(z) -> +inf, and omega_1(z) -> lambda_k from above.

# This means omega_1 is CONTINUOUS at nu_k with omega_1(nu_k) = lambda_k.
# It's NOT multi-valued on the real line.

# BUT in the complex plane, the equation G_p(w) = G_r(z) has n branches.
# omega_1 is one specific branch, and the others are omega_1^{(2)}, ..., omega_1^{(n)}.

# For the branch omega_1 mapping C^+ to C^+:
# We showed numerically that Im(omega_1(z)) > 0 when Im(z) > 0.
# But Im(omega_1(z)) can be LESS than Im(z).

# A real rational function omega_1(z) = P(z)/Q(z) with real coefficients
# maps C^+ to C^+ iff it has all its critical points in C^+_closure (?)
# Actually, a real rational function maps C^+ to C^+ iff between any two
# consecutive real poles, there is exactly one real zero (interlacing property).

n = 3
roots_p = np.array([-2.0, 0.0, 3.0])
roots_q = np.array([-1.0, 1.0, 2.0])
roots_r, c = boxplus_mss(roots_p, roots_q)
roots_r = np.sort(np.real(roots_r))

print(f"\nExample: n={n}")
print(f"p roots (lambda): {roots_p}")
print(f"q roots (mu):     {roots_q}")
print(f"r roots (nu):     {np.round(roots_r, 6)}")

# Compute omega_1 on the real line at many points
x_grid = np.linspace(min(roots_r)-5, max(roots_r)+5, 2000)
p_poly = np.poly(roots_p)
r_poly = np.poly(roots_r)
p_prime = np.polyder(p_poly)
r_prime = np.polyder(r_poly)

omega1_real = np.full(len(x_grid), np.nan)
for i, x in enumerate(x_grid):
    branches = compute_all_branches(x + 0j, roots_p, roots_r)
    # Select real branch closest to identity
    real_branches = branches[np.abs(np.imag(branches)) < 0.01]
    if len(real_branches) == 0:
        continue
    real_branches = np.real(real_branches)
    # The physical branch is the one that connects lambda_k at nu_k
    best = min(real_branches, key=lambda w: abs(w - x))
    omega1_real[i] = best

phi_real = omega1_real - x_grid

# Find poles of omega_1 by looking for discontinuities
dphi = np.diff(phi_real)
dx = np.diff(x_grid)
phi_deriv = dphi / dx

# The poles of omega_1 are where omega_1(x) -> +/- infinity
# These are the real solutions of Q(z) = 0 where omega_1 = P/Q.

# For our formula: omega_1 satisfies r(z)*p'(omega_1) - r'(z)*p(omega_1) = 0
# This defines omega_1 implicitly. The poles of omega_1 occur when the denominator
# of the rational function vanishes.

# Actually, omega_1 is NOT a rational function of degree n in general!
# omega_1 is defined by an ALGEBRAIC equation of degree n.
# It might NOT be globally rational -- it could be multi-valued.

# In the INFINITE (free probability) case, omega_1 IS rational (or at least holomorphic).
# In the FINITE MSS case, omega_1 is a branch of an algebraic function.

# Let me check: is omega_1 actually RATIONAL for the finite case?
# A rational function of degree n would be determined by 2n+1 parameters.
# omega_1(nu_k) = lambda_k (n conditions)
# omega_1'(nu_k) = 1 (n conditions)
# Plus behavior at infinity: omega_1(z) ~ z + c_0 + c_1/z + ...
# That's 2n + 2 conditions for a degree-n rational function (2n+1 parameters).
# OVERDETERMINED! So omega_1 might NOT be rational.

# Let me test: compute omega_1 at many points and try to fit a rational function.

# For n=3: a degree-3 rational function P(z)/Q(z) has deg(P)=3, deg(Q)=2 (if leading
# behavior is z). So P = z^3 + p2*z^2 + p1*z + p0, Q = z^2 + q1*z + q0.
# 5 free parameters. Conditions:
# omega_1(nu_k) = lambda_k (3 conditions)
# omega_1'(nu_k) = 1 (3 conditions) ... total 6, but 5 parameters. Overdetermined.

# Unless the MSS structure makes these 6 conditions consistent with 5 params.

# Actually wait: I said deg(P) = 3, deg(Q) = 2. But omega_1 ~ z + const as z -> inf.
# So P/Q = z + c_0 + c_1/z + ... We need P(z)/Q(z) - z -> c_0 as z -> inf.
# P(z) = z*Q(z) + R(z) where deg(R) <= 1. So P = z*Q + a*z + b.
# = z*(z^2 + q1*z + q0) + a*z + b = z^3 + q1*z^2 + (q0+a)*z + b.
# omega_1(z) = z + (a*z + b)/Q(z) = z + (az+b)/(z^2+q1*z+q0).
# For large z: omega_1 ~ z + a/z + ... So c_0 = 0? That contradicts our earlier finding
# that c_0 = -mean(mu).

# Hmm, let me redo: P(z) = z * Q(z) + R(z) where deg R < deg Q = 2, so deg R <= 1.
# omega_1(z) = P(z)/Q(z) = z + R(z)/Q(z).
# R(z)/Q(z) -> 0 as z -> inf (since deg R < deg Q). But we found c_0 != 0!

# This means deg Q <= 1, or the leading coefficient of Q is not 1... let me re-think.

# omega_1(z) = z + c_0 + O(1/z). For rational: P(z)/Q(z) = z + c_0 + ...
# P(z) = (z + c_0)*Q(z) + R(z) where deg R < deg Q.
# So omega_1(z) = (z+c_0) + R(z)/Q(z).
# deg(omega_1) = max(deg P, deg Q).
# For omega_1 of degree n (as defined by the algebraic equation of degree n):
# We need to determine deg P and deg Q.

# For n = 3, the algebraic equation is degree 3 in w. A branch may or may not
# be rational. Let's just test numerically.

print("\nFitting rational function to omega_1...")

# Collect (z, omega_1(z)) pairs for real z values
valid = ~np.isnan(omega1_real) & ~np.isinf(omega1_real)
z_data = x_grid[valid]
w_data = omega1_real[valid]

# Remove points near discontinuities
keep = np.ones(len(z_data), dtype=bool)
for i in range(1, len(z_data)):
    if abs(w_data[i] - w_data[i-1]) > 10:
        keep[max(0,i-3):min(len(keep),i+4)] = False
z_data = z_data[keep]
w_data = w_data[keep]

# Try to fit omega_1(z) = (z^3 + a2*z^2 + a1*z + a0)/(z^2 + b1*z + b0)
# This is a degree-3/degree-2 rational function (5 parameters).
from scipy.optimize import least_squares

def rational_residual(params, z, w):
    a2, a1, a0, b1, b0 = params
    num = z**3 + a2*z**2 + a1*z + a0
    den = z**2 + b1*z + b0
    return (num/den - w) * np.abs(den)  # weight by |den| to avoid pole issues

res = least_squares(rational_residual, [0, 0, 0, 0, 1],
                    args=(z_data, w_data), method='lm')
a2, a1, a0, b1, b0 = res.x
print(f"Fit: omega_1(z) = (z^3 + {a2:.4f}*z^2 + {a1:.4f}*z + {a0:.4f}) / "
      f"(z^2 + {b1:.4f}*z + {b0:.4f})")
print(f"Residual norm: {np.sqrt(np.mean(res.fun**2)):.6e}")

# Check the fit at known points
for k in range(n):
    num = roots_r[k]**3 + a2*roots_r[k]**2 + a1*roots_r[k] + a0
    den = roots_r[k]**2 + b1*roots_r[k] + b0
    fit_val = num/den
    print(f"  omega_1(nu_{k}) = {fit_val:.6f} (should be {roots_p[k]:.6f}, "
          f"error = {abs(fit_val - roots_p[k]):.2e})")

# Check derivative
for k in range(n):
    z = roots_r[k]
    num = z**3 + a2*z**2 + a1*z + a0
    num_d = 3*z**2 + 2*a2*z + a1
    den = z**2 + b1*z + b0
    den_d = 2*z + b1
    omega1_d = (num_d*den - num*den_d) / den**2
    print(f"  omega_1'(nu_{k}) = {omega1_d:.6f} (should be 1.0)")

# Poles of the rational function
poles = np.roots([1, b1, b0])
print(f"\nPoles of omega_1: {poles}")
print(f"Poles of phi = omega_1 - z: same as omega_1")

# phi(z) = omega_1(z) - z = (z^3 + a2*z^2 + a1*z + a0 - z*(z^2+b1*z+b0)) / (z^2+b1*z+b0)
# = ((a2-b1)*z^2 + (a1-b0)*z + a0) / (z^2 + b1*z + b0)
# This is a degree-2/degree-2 rational function (if a2 != b1).
# If a2 = b1, it reduces to degree 1/degree 2.

phi_num_coeffs = [a2-b1, a1-b0, a0]
print(f"\nphi(z) = ({a2-b1:.4f}*z^2 + {a1-b0:.4f}*z + {a0:.4f}) / (z^2 + {b1:.4f}*z + {b0:.4f})")
print(f"Leading coeff of phi numerator: {a2-b1:.4f}")
print(f"(This should be -mean(mu) = {-np.mean(roots_q):.4f} if phi(z) -> c_0)")

# phi(z) -> (a2-b1) as z -> infinity.
# We expect c_0 = -mean(mu).
print(f"\nphis(inf) = a2-b1 = {a2-b1:.6f}")
print(f"-mean(mu) = {-np.mean(roots_q):.6f}")


# ================================================================
# Now let's derive the PARTIAL FRACTION form of phi
# ================================================================
print("\n\n" + "="*70)
print("PARTIAL FRACTION DECOMPOSITION OF phi")
print("="*70)

# phi(z) = (c2*z^2 + c1*z + c0) / (z^2 + b1*z + b0)
# = c2 + (c1 - c2*b1)*z + (c0 - c2*b0)) / (z^2 + b1*z + b0)
# = c2 + ((c1-c2*b1)*z + (c0-c2*b0)) / ((z-p1)(z-p2))
# where p1, p2 are the poles.

c2 = a2 - b1
c1 = a1 - b0
c0_phi = a0

print(f"phi(z) = {c2:.4f} + ({c1-c2*b1:.4f}*z + {c0_phi-c2*b0:.4f}) / (z^2 + {b1:.4f}*z + {b0:.4f})")

if np.all(np.abs(np.imag(poles)) < 0.01):
    poles_real = np.real(np.sort(poles))
    p1, p2 = poles_real

    # Partial fractions: R(z)/((z-p1)(z-p2)) = A1/(z-p1) + A2/(z-p2)
    # R(z) = (c1-c2*b1)*z + (c0_phi - c2*b0)
    R_p1 = (c1-c2*b1)*p1 + (c0_phi - c2*b0)
    R_p2 = (c1-c2*b1)*p2 + (c0_phi - c2*b0)

    A1 = R_p1 / (p1 - p2)
    A2 = R_p2 / (p2 - p1)

    print(f"\nphi(z) = {c2:.4f} + {A1:.4f}/(z - {p1:.4f}) + {A2:.4f}/(z - {p2:.4f})")
    print(f"Residue at z={p1:.4f}: {A1:.4f}")
    print(f"Residue at z={p2:.4f}: {A2:.4f}")

    # For Herglotz (with convention 1/(d-z)):
    # phi(z) = c2 + (-A1)/(p1-z) + (-A2)/(p2-z)
    # Herglotz requires -A1 > 0 and -A2 > 0, i.e., A1 < 0 and A2 < 0.
    print(f"\nIn Herglotz form phi(z) = {c2:.4f} + {-A1:.4f}/({p1:.4f}-z) + {-A2:.4f}/({p2:.4f}-z)")
    print(f"  For Herglotz: need c_1 = {-A1:.4f} > 0: {'YES' if -A1 > 0 else 'NO'}")
    print(f"  For Herglotz: need c_2 = {-A2:.4f} > 0: {'YES' if -A2 > 0 else 'NO'}")

    if A1 < 0 and A2 < 0:
        print("\n  phi IS in Herglotz form! But then phi' > 0 everywhere.")
        print("  This contradicts phi'(nu_k) = 0.")
        print("  Let's check phi' at nu_k directly...")
        for k in range(n):
            z = roots_r[k]
            phi_d = (-A1)/(p1-z)**2 + (-A2)/(p2-z)**2
            print(f"    phi'(nu_{k}) = {phi_d:.8f} (should be 0)")
    else:
        print("\n  phi is NOT in Herglotz form. Some residues have wrong sign.")
        print("  This is consistent with phi'(nu_k) = 0.")

# ================================================================
# The ACTUAL form of phi and alpha
# ================================================================
print("\n\n" + "="*70)
print("DIRECT COMPUTATION OF alpha_k FROM PARTIAL FRACTIONS")
print("="*70)

# phi(z) = c_0 + A1/(z-p1) + A2/(z-p2)  (where some A_j may be positive)
# phi'(z) = -A1/(z-p1)^2 - A2/(z-p2)^2
# phi''(z) = 2*A1/(z-p1)^3 + 2*A2/(z-p2)^3

# alpha_k = phi''(nu_k)/2 = A1/(nu_k-p1)^3 + A2/(nu_k-p2)^3

# h_k = H_r(nu_k) = sum_{j!=k} 1/(nu_k - nu_j)

# <h, alpha> = sum_k h_k * alpha_k
# = sum_k [sum_{j!=k} 1/(nu_k-nu_j)] * [A1/(nu_k-p1)^3 + A2/(nu_k-p2)^3]
# = A1 * sum_k h_k/(nu_k-p1)^3 + A2 * sum_k h_k/(nu_k-p2)^3

# So <h,alpha> = A1 * S(p1) + A2 * S(p2)
# where S(p) = sum_k H_r(nu_k) / (nu_k - p)^3

if np.all(np.abs(np.imag(poles)) < 0.01):
    poles_real = np.real(np.sort(poles))
    p1, p2 = poles_real

    h = H_values(roots_r)
    alpha_direct = h * 0
    for k in range(n):
        alpha_direct[k] = A1/(roots_r[k]-p1)**3 + A2/(roots_r[k]-p2)**3

    alpha_from_H = H_values(roots_p) - h

    print(f"alpha from partial fractions: {alpha_direct}")
    print(f"alpha from H_p - H_r:        {alpha_from_H}")
    print(f"Match: {np.allclose(alpha_direct, alpha_from_H, atol=0.01)}")

    h_alpha = np.dot(h, alpha_from_H)
    print(f"\n<h,alpha> = {h_alpha:.8f}")

    S_p1 = sum(h[k]/(roots_r[k]-p1)**3 for k in range(n))
    S_p2 = sum(h[k]/(roots_r[k]-p2)**3 for k in range(n))
    h_alpha_decomp = A1*S_p1 + A2*S_p2
    print(f"<h,alpha> = A1*S(p1) + A2*S(p2) = {A1:.4f}*{S_p1:.4f} + {A2:.4f}*{S_p2:.4f} = {h_alpha_decomp:.8f}")

    print(f"\nS(p1) = sum h_k/(nu_k-p1)^3 = {S_p1:.6f}")
    print(f"S(p2) = sum h_k/(nu_k-p2)^3 = {S_p2:.6f}")
    print(f"A1 = {A1:.6f}, A2 = {A2:.6f}")


# ================================================================
# ALTERNATIVE APPROACH: Use the identity for <h, alpha> directly
# ================================================================
print("\n\n" + "="*70)
print("ALTERNATIVE: CONTOUR INTEGRAL APPROACH")
print("="*70)

# <h, alpha> = sum_k H_r(nu_k) * (H_p(lambda_k) - H_r(nu_k))
#            = sum_k H_r(nu_k) * H_p(lambda_k) - Phi_n(r)
#
# = sum_k [sum_{j!=k} 1/(nu_k-nu_j)] * [sum_{l!=k} 1/(lambda_k-lambda_l)] - Phi_n(r)
#
# This double sum is complicated. Can we relate it to a contour integral?
#
# Consider the function f(z) = nG_r(z) * [nG_p(omega_1(z)) - nG_r(z)]
# = nG_r(z) * [nG_r(z) - nG_r(z)] = 0 ... no, G_r(z) = G_p(omega_1(z)), so this is 0.
#
# Try: f(z) = (1/2pi i) * oint [nG_r(z)]^2 * phi'(z) dz
# where the contour encloses all nu_k.
#
# Near nu_k: G_r(z) = 1/(n(z-nu_k)) + H_r(nu_k)/n + O(z-nu_k)
# phi(z) = phi(nu_k) + phi''(nu_k)/2 * (z-nu_k)^2 + O((z-nu_k)^3)
# phi'(z) = phi''(nu_k) * (z-nu_k) + O((z-nu_k)^2) = 2*alpha_k*(z-nu_k) + O((z-nu_k)^2)
#
# [nG_r(z)]^2 = 1/(z-nu_k)^2 + 2*H_r(nu_k)/(z-nu_k) + O(1)
#
# [nG_r(z)]^2 * phi'(z) = [1/(z-nu_k)^2 + 2*H_r(nu_k)/(z-nu_k) + O(1)]
#                         * [2*alpha_k*(z-nu_k) + O((z-nu_k)^2)]
# = 2*alpha_k/(z-nu_k) + O(1) + 4*H_r(nu_k)*alpha_k + O(z-nu_k)
# Residue at nu_k = 2*alpha_k
#
# Hmm, the contour integral gives sum of residues = 2*sum_k alpha_k.
# This doesn't directly give <h,alpha>.
#
# Try: f(z) = [nG_r(z)]^2 * phi(z)
# = [1/(z-nu_k)^2 + 2H_r(nu_k)/(z-nu_k) + ...] * [phi(nu_k) + phi''(nu_k)/2*(z-nu_k)^2 + ...]
# = phi(nu_k)/(z-nu_k)^2 + 2*H_r(nu_k)*phi(nu_k)/(z-nu_k) + ... + phi''(nu_k)/2 + ...
# Residue at nu_k = 2*H_r(nu_k)*phi(nu_k) = 2*h_k*(lambda_k - nu_k)
#
# This gives sum_k h_k*(lambda_k - nu_k), not <h,alpha>.
#
# Try: f(z) = [nG_r'(z)] * phi(z)
# G_r'(z) = -(1/n)*sum_k 1/(z-nu_k)^2
# nG_r'(z) = -sum_k 1/(z-nu_k)^2
# f(z) = [-sum_k 1/(z-nu_k)^2] * phi(z)
# Near nu_l: f has pole of order 2 from the 1/(z-nu_l)^2 term.
# Residue of 1/(z-nu_l)^2 * phi(z) at z = nu_l:
# = d/dz[phi(z)] at z = nu_l = phi'(nu_l) = 0! (since omega_1'(nu_k)=1)
#
# So f(z) actually has REMOVABLE singularity at each nu_k from the 1/(z-nu_k)^2 * phi(z) term!
# The other terms (cross terms from j != k) give simple poles:
# f(z) near nu_l: -sum_{k!=l} 1/(nu_l-nu_k)^2 * phi(nu_l) + (analytic)
# Wait, this isn't quite right. Let me be more careful.

# f(z) = -[sum_k 1/(z-nu_k)^2] * phi(z)
# This is a product of two functions. The first has poles of order 2 at each nu_k.
# phi(z) is analytic at nu_k (since omega_1(nu_k) = lambda_k is finite).
#
# The ONLY poles of f are at nu_k (from G_r') and at the poles of phi (= poles of omega_1).
#
# Residue of f at nu_k (from the pole of G_r'):
# We need the residue of [1/(z-nu_k)^2] * phi(z) at z = nu_k, PLUS the residues
# from the other terms [1/(z-nu_j)^2] * phi(z) at z = nu_k (which are 0 since
# 1/(z-nu_j)^2 is analytic at nu_k for j != k).
#
# Res_{z=nu_k} [1/(z-nu_k)^2 * phi(z)] = phi'(nu_k) = 0.
#
# So f(z) has NO residue at any nu_k! All residues come from the poles of phi.

# This means: for a contour enclosing all nu_k but no poles of phi:
# oint f(z) dz = 2*pi*i * sum of residues at nu_k = 0.
# And for a contour enclosing everything:
# oint f(z) dz = 0 (by residue theorem for a rational function decaying at infinity)
# So sum of residues at poles of phi = 0 as well.

# This doesn't help. Let me try a DIFFERENT integrand.

# THE KEY INTEGRAND:
# Consider I = (1/2pi i) oint [nG_r(z)]^3 * phi(z) dz (enclosing all nu_k)
#
# Near nu_k: nG_r(z) = 1/(z-nu_k) + H_r(nu_k) + c_2*(z-nu_k) + ...
# where c_2 = sum_{j!=k} 1/(nu_k-nu_j)^3 / ... actually let me just compute:
# nG_r(z) = sum_j 1/(z - nu_j)
# Near z = nu_k: nG_r(z) = 1/(z-nu_k) + sum_{j!=k} 1/(z-nu_j)
# = 1/(z-nu_k) + sum_{j!=k} [1/(nu_k-nu_j) + (z-nu_k)/(nu_k-nu_j)^2 + ...]
# = 1/(z-nu_k) + H_r(nu_k) + (z-nu_k)*sum_{j!=k} 1/(nu_k-nu_j)^2 + O((z-nu_k)^2)
#
# Let H'_k = sum_{j!=k} 1/(nu_k-nu_j)^2 (derivative of H at the root)
# Then nG_r(z) = 1/epsilon + h_k + H'_k*epsilon + ... where epsilon = z - nu_k.
#
# [nG_r(z)]^3 = 1/epsilon^3 + 3*h_k/epsilon^2 + 3*(h_k^2 + H'_k)/epsilon + ...
# Wait: (a+b+c*e+...)^3 = a^3 + 3a^2*b + 3a^2*c*e + 3a*b^2 + ...
# With a = 1/e, b = h_k, c = H'_k*e:
# = 1/e^3 + 3h_k/e^2 + 3(h_k^2 + H'_k)/e + ...

# [nG_r(z)]^3 * phi(z) where phi(z) = phi_k + alpha_k*e^2 + ...
# = (phi_k/e^3 + 3h_k*phi_k/e^2 + [3(h_k^2+H'_k)*phi_k + alpha_k]/e + ...)
#
# Residue at nu_k = 3(h_k^2 + H'_k)*phi_k + alpha_k
# where phi_k = lambda_k - nu_k.
#
# This is too complicated and doesn't directly give <h,alpha>.

# Let me try the simplest possible approach:
# I = (1/2pi i) oint n^2 * G_r(z) * G_r'(z) * phi(z) dz  (enclosing all nu_k)
#
# n*G_r(z) = 1/(z-nu_k) + h_k + H'_k*(z-nu_k) + ...
# n*G_r'(z) = -1/(z-nu_k)^2 - H'_k + ...
#
# n^2 * G_r(z) * G_r'(z) = [1/(z-nu_k) + h_k + ...] * [-1/(z-nu_k)^2 - H'_k + ...]
# = -1/(z-nu_k)^3 - H'_k/(z-nu_k) - h_k/(z-nu_k)^2 + ...
# = -1/(z-nu_k)^3 - h_k/(z-nu_k)^2 - H'_k/(z-nu_k) + ...
#
# Times phi(z) = phi_k + alpha_k*(z-nu_k)^2 + ...:
# Residue at nu_k from -1/(z-nu_k)^3 * phi(z): coefficient of (z-nu_k)^2 in phi = alpha_k
# Residue at nu_k from -h_k/(z-nu_k)^2 * phi(z): coefficient of (z-nu_k) in phi = phi'(nu_k) = 0
# Residue at nu_k from -H'_k/(z-nu_k) * phi(z): phi(nu_k) = phi_k = lambda_k - nu_k
#
# Total residue at nu_k = alpha_k - H'_k * phi_k
#
# So I = sum_k (alpha_k - H'_k * phi_k).
# This gives sum alpha_k, NOT <h,alpha>.

# I'm not finding a clean contour integral for <h,alpha>.
# Let me try one more thing:

# I = sum_k Res_{z=nu_k} [nG_r(z)]^2 * phi'(z) * g(z) dz
# where g(z) is chosen to extract h_k from the pole.
#
# Near nu_k: [nG_r(z)]^2 = 1/e^2 + 2h_k/e + (h_k^2 + 2H'_k) + ...
# phi'(z) = 2*alpha_k*e + ...
#
# [nG_r]^2 * phi' = 2*alpha_k/e + 2*alpha_k*(h_k^2+2H'_k)*e + ... + 4*h_k*alpha_k + ...
# Actually: (1/e^2 + 2h_k/e + ...)(2*alpha_k*e + ...) = 2*alpha_k/e + 4*h_k*alpha_k + ...
# Residue at nu_k = 2*alpha_k.
#
# If I multiply by nG_r(z) to get [nG_r]^3 * phi':
# (1/e^3 + ...)(2*alpha_k*e + ...) = 2*alpha_k/e^2 + ...
# Residue is the COEFFICIENT of 1/e, which involves phi''' and h_k.

# This is getting nowhere fast. Let me try a completely different approach.

print("\nContour integral approach: appears to NOT directly yield <h,alpha>.")
print("The problem is that <h,alpha> mixes the pole structure of G_r with phi''.")
print("No simple contour integral seems to isolate this combination.")

# ================================================================
# TRY: MATRIX/RESOLVENT APPROACH
# ================================================================
print("\n\n" + "="*70)
print("MATRIX RESOLVENT APPROACH")
print("="*70)

# Key idea: Think of h_k and alpha_k in terms of the Cauchy matrix.
#
# Let C_{ij} = 1/(nu_i - nu_j) for i != j and C_{ii} = 0.
# Then h_k = sum_j C_{kj}.
# And alpha_k = H_p(lambda_k) - H_r(nu_k) = sum_{j!=k} [1/(lambda_k-lambda_j) - 1/(nu_k-nu_j)]
#
# <h, alpha> = sum_k (sum_j C_{kj}) * alpha_k
#
# Can we write this as a trace?
# h = C * 1  (where 1 is the all-ones vector)
# <h, alpha> = (C*1)^T * alpha = 1^T * C^T * alpha = 1^T * C * alpha (C is antisymmetric? No.)
#
# C_{kj} = 1/(nu_k - nu_j) is ANTI-symmetric: C_{kj} = -C_{jk}.
# So h_k = sum_j C_{kj} and h = C*1 where C is anti-symmetric.
# <h, alpha> = (C*1)^T * alpha = 1^T * C^T * alpha = -1^T * C * alpha = -sum_{k,j} C_{kj}*alpha_j / no wait.
# = sum_k (sum_j C_{kj}) * alpha_k = sum_{k,j} C_{kj} * alpha_k
# = sum_{k < j} C_{kj}*(alpha_k - alpha_j)  (using antisymmetry of C)
# = sum_{k < j} (alpha_k - alpha_j) / (nu_k - nu_j)

# AH HA! This is a key reformulation:
# <h,alpha> = sum_{k < j} (alpha_k - alpha_j) / (nu_k - nu_j)

# Since nu_1 < nu_2 < ... < nu_n, the denominators are negative for k < j.
# nu_k - nu_j < 0 for k < j.

# So <h,alpha> >= 0 iff sum_{k<j} (alpha_k - alpha_j)/(nu_k - nu_j) >= 0.
# Equivalently: sum_{k<j} (alpha_j - alpha_k)/(nu_j - nu_k) >= 0.
# I.e., sum_{k<j} (alpha_j - alpha_k)/(nu_j - nu_k) >= 0.

# This is a sum of terms (alpha_j - alpha_k)/(nu_j - nu_k).
# Each term is like a "difference quotient" of alpha as a function of nu.

# If alpha were a CONVEX function of nu (alpha(nu_k) convex in k),
# then (alpha_j - alpha_k)/(nu_j - nu_k) would be increasing in j for fixed k,
# and... hmm, that doesn't directly give the sum being positive.

# Actually, the sum sum_{k<j} (alpha_j - alpha_k)/(nu_j - nu_k)
# = sum_{k<j} Delta_{kj}(alpha) / Delta_{kj}(nu)
# where Delta_{kj}(f) = f_j - f_k.

# This is related to the DIVIDED DIFFERENCES of alpha with respect to nu.

# For n=3, verify:
h_direct = H_values(roots_r)
alpha_direct = H_values(roots_p) - h_direct

ha_direct = np.dot(h_direct, alpha_direct)
ha_antisym = 0
for k in range(n):
    for j in range(k+1, n):
        ha_antisym += (alpha_direct[j] - alpha_direct[k]) / (roots_r[j] - roots_r[k])

print(f"\n<h,alpha> = {ha_direct:.8f}")
print(f"sum_{{k<j}} (alpha_j - alpha_k)/(nu_j - nu_k) = {ha_antisym:.8f}")
print(f"Match: {abs(ha_direct - ha_antisym) < 1e-8}")

# Let's investigate the terms
print(f"\nIndividual terms (alpha_j - alpha_k)/(nu_j - nu_k):")
for k in range(n):
    for j in range(k+1, n):
        term = (alpha_direct[j] - alpha_direct[k]) / (roots_r[j] - roots_r[k])
        print(f"  k={k}, j={j}: ({alpha_direct[j]:.4f} - {alpha_direct[k]:.4f})/({roots_r[j]:.4f} - {roots_r[k]:.4f}) = {term:.6f}")

# ================================================================
# LARGE-SCALE TEST: Is (alpha_j - alpha_k)/(nu_j - nu_k) always positive?
# ================================================================
print("\n\n" + "="*70)
print("TEST: ARE INDIVIDUAL TERMS (alpha_j - alpha_k)/(nu_j - nu_k) POSITIVE?")
print("="*70)

np.random.seed(42)
all_terms_positive = 0
some_negative = 0
total = 0

for trial in range(500):
    n = np.random.choice([3, 4, 5])
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        raw = np.roots(np.poly(roots_r))
        if np.any(np.abs(np.imag(raw)) > 0.01):
            continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01):
            continue

        h = H_values(roots_r)
        u = H_values(roots_p)
        alpha = u - h

        total += 1
        all_pos = True
        for k in range(n):
            for j in range(k+1, n):
                term = (alpha[j] - alpha[k]) / (roots_r[j] - roots_r[k])
                if term < -1e-8:
                    all_pos = False
                    break
            if not all_pos:
                break

        if all_pos:
            all_terms_positive += 1
        else:
            some_negative += 1
    except:
        pass

print(f"All terms positive: {all_terms_positive}/{total}")
print(f"Some term negative: {some_negative}/{total}")

# If individual terms can be negative, the sum being positive is non-trivial.
# Let's check what sign pattern alpha has.

print("\n\n--- Alpha sign patterns ---")
np.random.seed(42)
for trial in range(10):
    n = 4
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.5:
            roots_p[i] = roots_p[i-1] + 0.5
    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.5:
            roots_q[i] = roots_q[i-1] + 0.5

    try:
        roots_r, _ = boxplus_mss(roots_p, roots_q)
        raw = np.roots(np.poly(roots_r))
        if np.any(np.abs(np.imag(raw)) > 0.01):
            continue
        roots_r = np.sort(np.real(raw))
        if np.any(np.diff(roots_r) < 0.01):
            continue

        h = H_values(roots_r)
        u = H_values(roots_p)
        alpha = u - h

        signs = ['+'if a > 0 else '-' for a in alpha]
        print(f"  alpha signs: {signs}, alpha = {np.round(alpha,3)}")
        print(f"  h signs:     {['+'if x>0 else '-' for x in h]}, h = {np.round(h,3)}")
        print(f"  <h,alpha> = {np.dot(h,alpha):.6f}")
        print()
    except:
        pass


# ================================================================
# THE DIVIDED DIFFERENCE REFORMULATION
# ================================================================
print("\n" + "="*70)
print("THE DIVIDED DIFFERENCE REFORMULATION")
print("="*70)
print("""
KEY IDENTITY:
<h,alpha> = sum_{k<j} (alpha_j - alpha_k) / (nu_j - nu_k)

where alpha_k = H_p(lambda_k) - H_r(nu_k) = u_k - h_k.

Substituting: alpha_j - alpha_k = (u_j - h_j) - (u_k - h_k) = (u_j - u_k) - (h_j - h_k)

So: <h,alpha> = sum_{k<j} [(u_j-u_k) - (h_j-h_k)] / (nu_j-nu_k)
             = sum_{k<j} (u_j-u_k)/(nu_j-nu_k) - sum_{k<j} (h_j-h_k)/(nu_j-nu_k)

The second sum is sum_{k<j} (h_j-h_k)/(nu_j-nu_k) = sum_k h_k * sum_{j!=k} 1/(nu_k-nu_j)
Wait, let me redo this carefully.

sum_{k<j} (h_j - h_k)/(nu_j - nu_k)
= sum_{k<j} h_j/(nu_j-nu_k) - sum_{k<j} h_k/(nu_j-nu_k)
= sum_j h_j * sum_{k<j} 1/(nu_j-nu_k) - sum_k h_k * sum_{j>k} 1/(nu_j-nu_k)
= sum_j h_j * sum_{k<j} 1/(nu_j-nu_k) + sum_k h_k * sum_{j>k} 1/(nu_k-nu_j)

Hmm, with the change k<->j in the second sum:
= sum_j h_j * [sum_{k<j} 1/(nu_j-nu_k) + sum_{l>j} 1/(nu_j-nu_l)]
= sum_j h_j * sum_{m!=j} 1/(nu_j-nu_m)
= sum_j h_j * h_j
= ||h||^2

So: sum_{k<j} (h_j-h_k)/(nu_j-nu_k) = ||h||^2 = Phi_r

And: <h,alpha> = sum_{k<j} (u_j-u_k)/(nu_j-nu_k) - Phi_r

Equivalently: <h,alpha> + Phi_r = sum_{k<j} (u_j-u_k)/(nu_j-nu_k) = <h,u>

So <h,alpha> = <h,u> - ||h||^2. This is just the definition. Circular.

BUT: the reformulation <h,alpha> = sum_{k<j} (alpha_j-alpha_k)/(nu_j-nu_k)
is still useful because it expresses <h,alpha> in terms of "divided differences"
of alpha with respect to nu.
""")

# Verify the identity
h = H_values(roots_r)
u = H_values(roots_p)
alpha = u - h

lhs = np.dot(h, alpha)
rhs = 0
for k in range(n):
    for j in range(k+1, n):
        rhs += (u[j] - u[k]) / (roots_r[j] - roots_r[k])
rhs -= sum(h**2)

print(f"<h,alpha> = {lhs:.8f}")
print(f"sum (u_j-u_k)/(nu_j-nu_k) - Phi_r = {rhs:.8f}")
print(f"Match: {abs(lhs - rhs) < 1e-8}")
