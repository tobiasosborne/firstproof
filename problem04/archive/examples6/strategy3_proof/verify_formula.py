"""
Verify which MSS formula is correct by checking multiple n=2 cases exactly.

For n=2: p(x)=(x-a)(x-b), q(x)=(x-c)(x-d).
We KNOW r = p boxplus_2 q has roots:
  mean = (a+b+c+d)/2 +/- sqrt((a-b)^2/4 + (c-d)^2/4)

And the inequality is EXACT EQUALITY.

Check that the formula reproduces this.
"""

import numpy as np
from math import comb
from itertools import combinations

def elem_sym_poly(roots, k):
    n = len(roots)
    if k == 0: return 1.0
    if k > n: return 0.0
    return sum(np.prod([roots[i] for i in subset])
               for subset in combinations(range(n), k))

def boxplus_hat(roots_p, roots_q):
    """Formula 1: hat-e convolution."""
    n = len(roots_p)
    ep = [elem_sym_poly(roots_p, k) for k in range(n+1)]
    eq = [elem_sym_poly(roots_q, k) for k in range(n+1)]

    hat_ep = [ep[k] / comb(n, k) for k in range(n+1)]
    hat_eq = [eq[k] / comb(n, k) for k in range(n+1)]

    hat_er = np.zeros(n+1)
    for k in range(n+1):
        for j in range(k+1):
            hat_er[k] += hat_ep[j] * hat_eq[k-j]

    er = [hat_er[k] * comb(n, k) for k in range(n+1)]

    coeffs = np.zeros(n+1)
    for k in range(n+1):
        coeffs[k] = (-1)**k * er[k]

    roots = np.sort(np.real(np.roots(coeffs)))
    return roots, er

def boxplus_2_exact(a, b, c, d):
    """Exact n=2 formula."""
    mean_r = (a + b + c + d) / 2
    var_r = (a - b)**2 / 4 + (c - d)**2 / 4
    return np.array([mean_r - np.sqrt(var_r), mean_r + np.sqrt(var_r)])


# Test many n=2 cases
np.random.seed(42)
print("Checking hat-convolution formula against exact n=2:")
max_error = 0
for trial in range(100):
    roots_p = np.sort(np.random.randn(2) * 3)
    if abs(roots_p[1] - roots_p[0]) < 0.1:
        roots_p[1] = roots_p[0] + 0.5
    roots_q = np.sort(np.random.randn(2) * 3)
    if abs(roots_q[1] - roots_q[0]) < 0.1:
        roots_q[1] = roots_q[0] + 0.5

    r_formula, _ = boxplus_hat(roots_p, roots_q)
    r_exact = boxplus_2_exact(roots_p[0], roots_p[1], roots_q[0], roots_q[1])

    err = np.max(np.abs(r_formula - r_exact))
    max_error = max(max_error, err)

    if err > 1e-10:
        print(f"  Trial {trial}: ERROR={err:.2e}")
        print(f"    p={roots_p}, q={roots_q}")
        print(f"    Formula: {r_formula}")
        print(f"    Exact:   {r_exact}")

print(f"Max error over 100 trials: {max_error:.2e}")

# If formula is correct for n=2, then the failures were from the subordination computation.
# Let's check the subordination more carefully.

print("\n" + "="*60)
print("CHECKING SUBORDINATION COMPUTATION")
print("="*60)

def H_values(roots):
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if j != i:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H


# For n=2, let's compute everything analytically
def n2_full_analysis(a, b, c, d):
    """Full analytical computation for n=2."""
    assert a < b and c < d

    s = b - a  # gap of p
    t = d - c  # gap of q

    nu1, nu2 = boxplus_2_exact(a, b, c, d)
    R = nu2 - nu1  # gap of r = sqrt(s^2 + t^2)

    Phi_p = 2/s**2
    Phi_q = 2/t**2
    Phi_r = 2/R**2

    # Subordination: omega_1(nu_k) = lambda_{sigma(k)} with sigma = identity
    # omega_1(nu_1) = a, omega_1(nu_2) = b

    # omega_1'(nu_k) = 1 (proven)

    # omega_1''(nu_k):
    # From implicit differentiation of G_r(z) = G_p(omega_1(z)):
    # G_r'(z) = G_p'(omega_1(z)) * omega_1'(z)
    # G_r''(z) = G_p''(omega_1(z)) * (omega_1'(z))^2 + G_p'(omega_1(z)) * omega_1''(z)

    # At z = nu_k: omega_1'(nu_k) = 1
    # omega_1''(nu_k) = [G_r''(nu_k) - G_p''(lambda_k)] / G_p'(lambda_k)

    # For degree n polynomial p(x) = (x-a)(x-b):
    # G_p(z) = (1/2) * p'(z)/p(z) = (1/2) * [(z-b)+(z-a)] / [(z-a)(z-b)]
    #        = (1/2) * [2z - (a+b)] / [(z-a)(z-b)]

    # But wait, the definition is G_p(z) = (1/n)*p'(z)/p(z).
    # For n=2: G_p(z) = (1/2)*[1/(z-a) + 1/(z-b)]

    # G_p'(z) = -(1/2)*[1/(z-a)^2 + 1/(z-b)^2]
    # G_p''(z) = (1/2)*[2/(z-a)^3 + 2/(z-b)^3] = [1/(z-a)^3 + 1/(z-b)^3]

    # At z = nu_k = a (i.e., omega_1(nu_1) = a):
    # G_p has a POLE at w=a. So we need to use the F(z,w) approach.

    # Actually we already computed omega_1'' via the F approach. Let me just
    # verify numerically for n=2.

    # For n=2 with p=(x-a)(x-b), r=(x-nu1)(x-nu2):
    # F(z,w) = r'(z)*p(w) - p'(w)*r(z)
    #        = (2z - nu1 - nu2)*(w-a)*(w-b) - (2w - a - b)*(z-nu1)*(z-nu2)

    # At (nu1, a): r(nu1)=0, p(a)=0
    # F_w = r'(nu1)*p'(a) = (2*nu1 - nu1 - nu2)*(a-b) ... wait
    # r'(nu1) = nu1 - nu2 = -R
    # p'(a) = a - b = -s

    # F_w = r'(nu1)*p'(a) = (-R)*(-s) = R*s
    # F_z = -p'(a)*r'(nu1) = -(-s)*(-R) = -R*s
    # omega_1'(nu1) = -F_z/F_w = R*s/(R*s) = 1 ✓

    # F_zz = -p'(a)*r''(nu1) = s * 2 = 2s  (since r''(x) = 2 for degree 2)
    # F_zw = r''(nu1)*p'(a) - p''(a)*r'(nu1) = 2*(-s) - 2*(-R) = -2s + 2R
    # F_ww = r'(nu1)*p''(a) = (-R)*2 = -2R

    # omega_1''(nu1) = -(F_zz + 2*F_zw + F_ww) / F_w
    #               = -(2s + 2(-2s + 2R) + (-2R)) / (Rs)
    #               = -(2s - 4s + 4R - 2R) / (Rs)
    #               = -(-2s + 2R) / (Rs)
    #               = (2s - 2R) / (Rs)
    #               = 2(s - R) / (Rs)

    omega1_pp_1 = 2*(s - R) / (R*s)

    # Similarly at (nu2, b):
    # r'(nu2) = nu2 - nu1 = R
    # p'(b) = b - a = s
    # F_w = R*s
    # F_zz = -s*2 = -2s
    # F_zw = 2*s - 2*R
    # F_ww = R*2 = 2R
    # omega_1''(nu2) = -(-2s + 2(2s-2R) + 2R) / (Rs)
    #               = -(-2s + 4s - 4R + 2R) / (Rs)
    #               = -(2s - 2R) / (Rs)
    #               = -2(s-R) / (Rs)
    #               = 2(R-s) / (Rs)

    omega1_pp_2 = 2*(R - s) / (R*s)

    alpha_1 = omega1_pp_1 / 2  # = (s - R) / (Rs)
    alpha_2 = omega1_pp_2 / 2  # = (R - s) / (Rs)

    # Similarly for omega_2:
    omega2_pp_1 = 2*(t - R) / (R*t)
    omega2_pp_2 = 2*(R - t) / (R*t)

    beta_1 = omega2_pp_1 / 2  # = (t - R) / (Rt)
    beta_2 = omega2_pp_2 / 2  # = (R - t) / (Rt)

    # h values: H_r(nu_k)
    h_1 = 1/(nu1 - nu2)  # = -1/R
    h_2 = 1/(nu2 - nu1)  # = 1/R
    h = np.array([-1/R, 1/R])

    alpha = np.array([alpha_1, alpha_2])
    beta = np.array([beta_1, beta_2])

    # u = h + alpha
    u_1 = h_1 + alpha_1  # = -1/R + (s-R)/(Rs) = (-s + s - R)/(Rs) = -R/(Rs) = -1/s
    u_2 = h_2 + alpha_2  # = 1/R + (R-s)/(Rs) = (s + R - s)/(Rs) = R/(Rs) = 1/s
    u = np.array([u_1, u_2])

    # Check: H_p(a) = 1/(a-b) = -1/s, H_p(b) = 1/(b-a) = 1/s ✓

    # <h, alpha> = h_1*alpha_1 + h_2*alpha_2
    #            = (-1/R)*(s-R)/(Rs) + (1/R)*(R-s)/(Rs)
    #            = (-(s-R) + (R-s)) / (R^2 * s)
    #            = (-s+R + R-s) / (R^2 * s)
    #            = 2(R-s) / (R^2 * s)

    h_alpha = 2*(R - s) / (R**2 * s)

    # Since R = sqrt(s^2 + t^2) > s (when t > 0), we have R > s, so <h,alpha> > 0 ✓
    # But wait, R > s iff t > 0, which is always true. So <h,alpha> > 0 always at n=2.

    # <h, beta>:
    h_beta = 2*(R - t) / (R**2 * t)
    # Similarly R > t, so <h,beta> > 0 always.

    # A = ||u||^2 - ||h||^2 = 2/s^2 - 2/R^2 = 2(R^2 - s^2)/(s^2*R^2) = 2*t^2/(s^2*R^2)
    A_val = 2*t**2/(s**2*R**2)

    # B = 2*s^2/(t^2*R^2)
    B_val = 2*s**2/(t**2*R**2)

    # AB = 4/(R^4) = (2/R^2)^2 = Phi_r^2 = ||h||^4
    AB_val = A_val * B_val  # = 4/(R^4)
    h4_val = (2/R**2)**2  # = 4/R^4

    print(f"  s={s:.4f}, t={t:.4f}, R={R:.6f}")
    print(f"  Phi_p={Phi_p:.6f}, Phi_q={Phi_q:.6f}, Phi_r={Phi_r:.6f}")
    print(f"  alpha = [{alpha_1:.6f}, {alpha_2:.6f}]")
    print(f"  beta  = [{beta_1:.6f}, {beta_2:.6f}]")
    print(f"  <h,alpha> = {h_alpha:.6f} (>0: {h_alpha > 0})")
    print(f"  <h,beta>  = {h_beta:.6f} (>0: {h_beta > 0})")
    print(f"  A = {A_val:.6f}, B = {B_val:.6f}")
    print(f"  AB = {AB_val:.10f}, h4 = {h4_val:.10f}, diff = {AB_val - h4_val:.2e}")

    return {
        'h_alpha': h_alpha, 'h_beta': h_beta,
        'A': A_val, 'B': B_val, 'AB': AB_val, 'h4': h4_val,
        's': s, 't': t, 'R': R,
    }

    # (h_1, h_2 are local variables in the function above)

print("\nn=2 analytical cases:")
print("\nCase 1: s=2, t=4")
n2_full_analysis(-1, 1, -2, 2)
print("\nCase 2: s=1, t=3")
n2_full_analysis(0, 1, 0, 3)

# Now check what goes wrong in the subordination numerical computation
print("\n\n" + "="*60)
print("DIAGNOSING SUBORDINATION NUMERICAL ISSUES")
print("="*60)

# The issue is that the "hat convolution" formula might NOT be the MSS formula.
# Let me check: for n=3, does the hat convolution give real roots?

from itertools import combinations

np.random.seed(123)
complex_root_count = 0
for trial in range(200):
    n = 3
    roots_p = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_p[i] - roots_p[i-1] < 0.3:
            roots_p[i] = roots_p[i-1] + 0.3

    roots_q = np.sort(np.random.randn(n) * 2)
    for i in range(1, n):
        if roots_q[i] - roots_q[i-1] < 0.3:
            roots_q[i] = roots_q[i-1] + 0.3

    roots_r, _ = boxplus_hat(roots_p, roots_q)
    if np.any(np.abs(np.imag(np.roots(np.poly(roots_r)))) > 0.1):
        complex_root_count += 1

    # More directly:
    ep = [elem_sym_poly(roots_p, k) for k in range(n+1)]
    eq = [elem_sym_poly(roots_q, k) for k in range(n+1)]
    hat_ep = [ep[k] / comb(n, k) for k in range(n+1)]
    hat_eq = [eq[k] / comb(n, k) for k in range(n+1)]
    hat_er = np.zeros(n+1)
    for k in range(n+1):
        for j in range(k+1):
            hat_er[k] += hat_ep[j] * hat_eq[k-j]
    er = [hat_er[k] * comb(n, k) for k in range(n+1)]
    coeffs = np.zeros(n+1)
    for k in range(n+1):
        coeffs[k] = (-1)**k * er[k]
    raw_roots = np.roots(coeffs)
    if np.any(np.abs(np.imag(raw_roots)) > 0.01):
        complex_root_count += 1
        if complex_root_count <= 3:
            print(f"  Trial {trial}: complex roots! {raw_roots}")
            print(f"    p={roots_p}, q={roots_q}")
            print(f"    e_k(p)={ep}, e_k(q)={eq}, e_k(r)={er}")

print(f"\nComplex root cases: {complex_root_count}/200")
print("If > 0, the formula is WRONG for n>=3.")

# Let's check against MC for a specific n=3 case
print("\n\nDirect n=3 comparison: formula vs MC")
roots_p = np.array([-2.0, 0.0, 2.0])
roots_q = np.array([-3.0, 0.0, 3.0])

r_formula, er_formula = boxplus_hat(roots_p, roots_q)
print(f"Formula roots: {r_formula}")
print(f"Formula e_k:   {er_formula}")

# MC with high accuracy
from investigate_exact import boxplus_mc
mc_ek = boxplus_mc(roots_p, roots_q, 500000)
mc_roots = np.sort(np.real(np.roots(
    [(-1)**k * mc_ek[k] for k in range(4)]
)))
print(f"MC roots (500k): {mc_roots}")
print(f"MC e_k:          {list(mc_ek)}")

# Also try a different n=3 case
roots_p = np.array([-1.0, 0.5, 2.0])
roots_q = np.array([-0.5, 1.0, 3.0])

r_formula, er_formula = boxplus_hat(roots_p, roots_q)
mc_ek = boxplus_mc(roots_p, roots_q, 500000)
mc_roots = np.sort(np.real(np.roots(
    [(-1)**k * mc_ek[k] for k in range(4)]
)))
print(f"\nAnother n=3 case:")
print(f"Formula roots: {r_formula}")
print(f"MC roots:      {mc_roots}")
print(f"Error:         {np.max(np.abs(r_formula - mc_roots)):.6f}")
