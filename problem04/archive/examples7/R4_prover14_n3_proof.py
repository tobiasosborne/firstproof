"""
PROVER-14 Part 7: Direct n=3 proof

For n=3 centered, MSS adds (sigma, tau) coordinate-wise.
We need: 1/Sm2(sigma_p+sigma_q, tau_p+tau_q) >= 1/Sm2(sigma_p,tau_p) + 1/Sm2(sigma_q,tau_q)

1/Sm2 is NOT concave, so we need a different approach.

KEY IDEA: Write 1/Sm2 = S2/Q = -6*sigma/Q(rho) where rho = 27*tau^2/(-4*sigma^3).

Then: 1/Sm2 = -6*sigma / Q(rho)

where Q depends only on rho, and sigma < 0.

The conjecture becomes:
  -6*(sigma_p+sigma_q) / Q(rho_r) >= -6*sigma_p/Q(rho_p) + (-6*sigma_q)/Q(rho_q)

i.e., (sigma_p+sigma_q)/Q(rho_r) <= sigma_p/Q(rho_p) + sigma_q/Q(rho_q)

(Note: sigma < 0 so dividing flips the inequality.)

Let a = -sigma_p > 0, b = -sigma_q > 0. Then -sigma_r = a + b.
Let Qp = Q(rho_p), Qq = Q(rho_q), Qr = Q(rho_r).

The conjecture becomes: (a+b)/Qr >= a/Qp + b/Qq.
This is Titu/Engel with numerators a, b and denominators Qp, Qq, Qr.

Actually more precisely: (a+b)/Qr >= a/Qp + b/Qq
<=> (a+b)*Qp*Qq >= Qr*(a*Qq + b*Qp)
<=> Qr <= Qp*Qq*(a+b) / (a*Qq + b*Qp)
This is the weighted harmonic mean bound we already found.

Alternative: ALONG THE LINE (sigma, tau) = t*(sigma_0, tau_0),
1/Sm2 is linear in t (since Sm2 = Q/(−6*sigma) = Q(rho)/(-6*t*sigma_0) and rho is independent of t).
So 1/Sm2 = t*(-6*sigma_0/Q(rho_0)).

This means: 1/Sm2 is homogeneous of degree 1 in (sigma, tau)!
(Since roots scale as lambda -> c*lambda when (sigma,tau) -> (c^2*sigma, c^3*tau),
and 1/Sm2 scales as c^2.)

Wait, that's not quite homogeneity of degree 1 in (sigma, tau) because
sigma has weight 2 and tau has weight 3.

Actually: if we scale (sigma, tau) -> (c*sigma, c^{3/2}*tau), then
rho = 27*tau^2/(-4*sigma^3) -> 27*c^3*tau^2/(-4*c^3*sigma^3) = rho.
And 1/Sm2 = -6*sigma/Q(rho) -> -6*c*sigma/Q(rho) = c * (1/Sm2).

So 1/Sm2 is homogeneous of degree 1 under the WEIGHTED scaling
(sigma, tau) -> (c*sigma, c^{3/2}*tau).

This is NOT the standard notion of concavity but rather a
QUASI-CONCAVITY on weighted cones.

NEW APPROACH: Use the parametrization (a, rho) where a = -sigma, rho = 27*tau^2/(-4*sigma^3).
Then 1/Sm2 = 6*a / Q(rho).
And under MSS: a_r = a_p + a_q, but rho_r depends on both a and rho values.
"""
import numpy as np
from itertools import combinations
from math import factorial
import sympy as sp

print("="*70)
print("PROVER-14 Part 7: DIRECT n=3 PROOF")
print("="*70)

def Sm2(roots):
    n = len(roots)
    return sum(1/(roots[i]-roots[j])**2 for i in range(n) for j in range(i+1, n))

def S2(roots):
    n = len(roots)
    return sum((roots[i]-roots[j])**2 for i in range(n) for j in range(i+1, n))

def Q(roots):
    return Sm2(roots) * S2(roots)

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
# Part 1: Express Sm2 in terms of sigma, tau for n=3
# =============================================================
print("\n--- Part 1: Sm2(sigma, tau) for n=3 ---\n")

# The cubic x^3 + sigma*x + tau = 0 with roots a < b < c (a+b+c=0)
# Discriminant: Delta = -4*sigma^3 - 27*tau^2 > 0

# Using the trigonometric solution:
# roots = 2*sqrt(-sigma/3) * cos(theta + 2k*pi/3) for k=0,1,2
# where cos(3*theta) = 3*sqrt(3)*tau / (2*(-sigma)^{3/2})

# Let R = sqrt(-sigma/3), so roots = 2R*cos(theta+2kpi/3)

# Pairwise differences:
# root_i - root_j = 2R*(cos(theta+2i*pi/3) - cos(theta+2j*pi/3))
# = -4R*sin((theta+(i+j)*pi/3))*sin((i-j)*pi/3)

# For (i,j) = (0,1), (0,2), (1,2):
# d_{01} = -4R*sin(theta+pi/3)*sin(-pi/3) = 4R*sin(theta+pi/3)*sin(pi/3) = 2R*sqrt(3)*sin(theta+pi/3)
# d_{02} = -4R*sin(theta+2pi/3)*sin(-2pi/3) = 4R*sin(theta+2pi/3)*sin(2pi/3) = 2R*sqrt(3)*sin(theta+2pi/3)
# d_{12} = -4R*sin(theta+pi)*sin(-pi/3) = 4R*sin(theta+pi)*sin(pi/3) = -2R*sqrt(3)*sin(theta)

# Wait, I need to be more careful. Let me use specific indices.
# Let the three roots be:
# c = 2R*cos(theta)          (largest, for theta in (0, pi/3))
# b = 2R*cos(theta + 2pi/3)  (middle)
# a = 2R*cos(theta + 4pi/3)  (smallest)

# where cos(3*theta) = 3*sqrt(3)*tau / (2*(-sigma)^{3/2}) = -tau / (2*R^3/sqrt(3))
# Actually cos(3*theta) = -tau*sqrt(3)/(2*R^3) ... let me just derive.

# For x^3 + sigma*x + tau = 0, substituting x = 2R*cos(phi):
# 8R^3*cos^3(phi) + 2R*sigma*cos(phi) + tau = 0
# Using cos^3 = (3cos+cos3)/4:
# 8R^3*(3cos(phi)+cos(3phi))/4 + 2R*sigma*cos(phi) + tau = 0
# 2R^3*(3cos(phi)+cos(3phi)) + 2R*sigma*cos(phi) + tau = 0
# (6R^3 + 2R*sigma)*cos(phi) + 2R^3*cos(3phi) + tau = 0

# Set 6R^3 + 2R*sigma = 0 -> R^2 = -sigma/3 -> R = sqrt(-sigma/3). Good.
# Then: 2R^3*cos(3phi) + tau = 0
# cos(3phi) = -tau/(2R^3) = -tau / (2*(-sigma/3)^{3/2})
# = -tau * 3^{3/2} / (2*(-sigma)^{3/2})
# = -3*sqrt(3)*tau / (2*(-sigma)^{3/2})

# And rho = 27*tau^2 / (-4*sigma^3) = cos^2(3*phi) * 4*(-sigma)^3 / (9*3) * 27 / (-4*sigma^3)
# Hmm, let me just compute:
# cos(3*phi) = -3*sqrt(3)*tau / (2*(-sigma)^{3/2})
# cos^2(3*phi) = 27*tau^2 / (4*(-sigma)^3) = rho
# So cos^2(3*phi) = rho, giving cos(3*phi) = +/- sqrt(rho)

# The roots are 2R*cos(phi), 2R*cos(phi+2pi/3), 2R*cos(phi+4pi/3)
# with phi in (0, pi/3) for distinct ordered roots.

print("Roots in trigonometric form:")
print("  x_k = 2*sqrt(-sigma/3) * cos(phi + 2k*pi/3), k=0,1,2")
print("  where cos(3*phi) = -3*sqrt(3)*tau / (2*(-sigma)^{3/2})")
print("  and phi in (0, pi/3)")
print()

# Now compute Sm2 in terms of R and phi:
# Pairwise differences squared:
# (x_i - x_j)^2 = 4R^2*(cos(phi+2i*pi/3) - cos(phi+2j*pi/3))^2
# = 4R^2 * 4*sin^2((phi+(i+j)*pi/3)) * sin^2((i-j)*pi/3)

# Sm2 = sum_{i<j} 1/(x_i-x_j)^2 = sum_{i<j} 1 / [4R^2 * 4*sin^2(...)*sin^2(...)]
# = 1/(16*R^4) * sum_{i<j} 1/(sin^2((phi+(i+j)*pi/3))*sin^2((i-j)*pi/3))

# This is getting complicated. Let me compute numerically and fit.

# Actually, I realize we can express Sm2 more cleanly.
# Sm2 = 1/d01^2 + 1/d02^2 + 1/d12^2
# where d_{ij} are the 3 pairwise gaps.

# Using trigonometric parametrization with R and phi:
# gaps: s = x_0 - x_1 = 2R*(cos(phi) - cos(phi+2pi/3))
# = 2R*2*sin(phi+pi/3)*sin(pi/3) = 4R*sin(phi+pi/3)*sqrt(3)/2 = 2R*sqrt(3)*sin(phi+pi/3)

# t = x_0 - x_2 = 2R*(cos(phi) - cos(phi+4pi/3))
# = 2R*2*sin(phi+2pi/3)*sin(2pi/3) = 2R*sqrt(3)*sin(phi+2pi/3)

# And x_1 - x_2 = 2R*sqrt(3)*sin(phi) [by similar computation]

# So d01 = 2R*sqrt(3)*sin(phi+pi/3), d02 = 2R*sqrt(3)*sin(phi+2pi/3), d12 = 2R*sqrt(3)*sin(phi)

# Sm2 = 1/(12R^4) * [1/sin^2(phi+pi/3) + 1/sin^2(phi+2pi/3) + 1/sin^2(phi)]

# S2 = d01^2 + d02^2 + d12^2 = 12R^4 * [sin^2(phi+pi/3) + sin^2(phi+2pi/3) + sin^2(phi)]

# Note: sin^2(x) + sin^2(x+2pi/3) + sin^2(x+4pi/3) = 3/2 (identity!)
# So S2 = 12R^4 * 3/2 = 18R^4 ... wait, S2 has units of length^2, not length^4.
# S2 = sum d_{ij}^2 = 12R^2 * [sin^2(phi+pi/3) + sin^2(phi+2pi/3) + sin^2(phi)]
# = 12R^2 * 3/2 = 18R^2

# Wait: d_{ij}^2 = (2R*sqrt(3)*sin(...))^2 = 12R^2 * sin^2(...)
# So S2 = 12R^2 * sum sin^2 = 12R^2 * 3/2 = 18R^2

# But also S2 = -6*sigma = -6*(-3R^2) = 18R^2. Checks out!

# And Sm2 = 1/(12R^2) * [1/sin^2(phi+pi/3) + 1/sin^2(phi+2pi/3) + 1/sin^2(phi)]

# Let F(phi) = 1/sin^2(phi) + 1/sin^2(phi+pi/3) + 1/sin^2(phi+2pi/3)

# Then Sm2 = F(phi) / (12R^2) = F(phi) / (2*S2/3) = 3*F(phi) / (2*S2)

# So Q = Sm2 * S2 = 3*F(phi)/2

# And 1/Sm2 = 2*S2 / (3*F(phi)) = 12R^2 / F(phi) = -4*sigma / F(phi)

# So Q = 3*F(phi)/2 depends ONLY on phi.
# And 1/Sm2 = -4*sigma / F(phi) = 4*|sigma| / F(phi).

print("\n--- KEY FORMULA ---")
print("F(phi) = 1/sin^2(phi) + 1/sin^2(phi+pi/3) + 1/sin^2(phi+2pi/3)")
print("Q = 3*F(phi)/2")
print("1/Sm2 = 4*|sigma| / F(phi)")
print("phi in (0, pi/3), cos(3*phi) = -3*sqrt(3)*tau / (2*|sigma|^{3/2})")
print()

# Verify
print("Numerical verification:")
for sigma_val in [-1.0, -3.0, -10.0]:
    for tau_val in [0.0, 0.5, 1.0]:
        disc = -4*sigma_val**3 - 27*tau_val**2
        if disc <= 0:
            continue
        roots = np.sort(np.real(np.roots([1, 0, sigma_val, tau_val])))
        R = np.sqrt(-sigma_val/3)
        cos3phi = -3*np.sqrt(3)*tau_val / (2*(-sigma_val)**1.5)
        if abs(cos3phi) > 1:
            continue
        phi = np.arccos(cos3phi) / 3

        F_phi = 1/np.sin(phi)**2 + 1/np.sin(phi+np.pi/3)**2 + 1/np.sin(phi+2*np.pi/3)**2
        Q_computed = 3*F_phi/2
        Q_actual = Q(roots)
        inv_sm2_computed = 4*abs(sigma_val) / F_phi
        inv_sm2_actual = 1/Sm2(roots)

        print(f"  sigma={sigma_val}, tau={tau_val}: Q={Q_actual:.6f} vs {Q_computed:.6f}, "
              f"1/Sm2={inv_sm2_actual:.6f} vs {inv_sm2_computed:.6f}")

# =============================================================
# Part 2: The conjecture in (sigma, phi) coordinates
# =============================================================
print("\n\n--- Part 2: Conjecture in (sigma, phi) coordinates ---\n")

print("""
sigma_r = sigma_p + sigma_q  (i.e., a_r = a_p + a_q where a = -sigma)
tau_r = tau_p + tau_q

1/Sm2 = 4*a / F(phi)

Conjecture: 4*a_r / F(phi_r) >= 4*a_p/F(phi_p) + 4*a_q/F(phi_q)
i.e., a_r/F(phi_r) >= a_p/F(phi_p) + a_q/F(phi_q)
i.e., (a_p+a_q)/F(phi_r) >= a_p/F(phi_p) + a_q/F(phi_q)

Now phi_r is determined by cos(3*phi_r) = -3*sqrt(3)*tau_r / (2*a_r^{3/2}).

With tau_r = tau_p + tau_q:
  cos(3*phi_r) = -3*sqrt(3)*(tau_p+tau_q) / (2*(a_p+a_q)^{3/2})

This depends on all four parameters (a_p, a_q, tau_p, tau_q).

The condition (a_p+a_q)/F(phi_r) >= a_p/F(phi_p) + a_q/F(phi_q)
can be rearranged as:
  F(phi_r) <= F_W where F_W = (a_p+a_q) / (a_p/F(phi_p) + a_q/F(phi_q))
            = weighted harmonic mean of F(phi_p) and F(phi_q) with weights a_p, a_q

So the conjecture reduces to:
  F(phi_r) <= WH(F(phi_p), F(phi_q); a_p, a_q)

where WH is the weighted harmonic mean and phi_r is determined by the addition formulas.
""")

# Study F(phi)
phi_vals = np.linspace(0.01, np.pi/3 - 0.01, 100)
F_vals = [1/np.sin(phi)**2 + 1/np.sin(phi+np.pi/3)**2 + 1/np.sin(phi+2*np.pi/3)**2 for phi in phi_vals]

# F(phi) has a minimum at phi = pi/6 (by symmetry, since the three angles are
# phi, phi+pi/3, phi+2pi/3, which are equally spaced when phi = pi/6 - pi/6 ... hmm)
# Actually at phi = pi/6: sin(pi/6) = 1/2, sin(pi/6+pi/3) = sin(pi/2) = 1,
# sin(pi/6+2pi/3) = sin(5pi/6) = 1/2
# F(pi/6) = 4 + 1 + 4 = 9

min_idx = np.argmin(F_vals)
print(f"F(phi) minimum: phi = {phi_vals[min_idx]:.6f} = {phi_vals[min_idx]*3/np.pi:.4f}*pi/3")
print(f"F_min = {F_vals[min_idx]:.6f}")
print(f"F(pi/6) = {1/np.sin(np.pi/6)**2 + 1/np.sin(np.pi/6+np.pi/3)**2 + 1/np.sin(np.pi/6+2*np.pi/3)**2:.6f}")

# So Q_min = 3*9/2 = 13.5. Matches our earlier finding!

# F is CONVEX on (0, pi/3)?
print(f"\nIs F(phi) convex on (0, pi/3)?")
dF = np.diff(F_vals) / np.diff(phi_vals)
ddF = np.diff(dF) / np.diff(phi_vals[:-1])
print(f"  min(F'') = {np.min(ddF):.6f}, max(F'') = {np.max(ddF):.6f}")
print(f"  F is {'CONVEX' if np.min(ddF) >= -1e-6 else 'NOT CONVEX'} on (0, pi/3)")

# =============================================================
# Part 3: Study phi_r as a function of phi_p, phi_q
# =============================================================
print("\n\n--- Part 3: phi_r as function of phi_p, phi_q ---\n")

print("""
phi is related to (a, tau) by:
  cos(3*phi) = -3*sqrt(3)*tau / (2*a^{3/2})

And tau = -2*a^{3/2}*cos(3*phi) / (3*sqrt(3))

Under addition:
  a_r = a_p + a_q
  tau_r = tau_p + tau_q = -2/(3*sqrt(3)) * [a_p^{3/2}*cos(3*phi_p) + a_q^{3/2}*cos(3*phi_q)]

So: cos(3*phi_r) = -3*sqrt(3)*tau_r / (2*a_r^{3/2})
   = [a_p^{3/2}*cos(3*phi_p) + a_q^{3/2}*cos(3*phi_q)] / (a_p+a_q)^{3/2}

Let alpha = a_p^{3/2}, beta = a_q^{3/2}, gamma = (a_p+a_q)^{3/2}.
Then cos(3*phi_r) = [alpha*cos(3*phi_p) + beta*cos(3*phi_q)] / gamma.

Note: alpha + beta < gamma (since (a+b)^{3/2} > a^{3/2}+b^{3/2} for a,b > 0).
So cos(3*phi_r) is a CONTRACTED weighted average of cos(3*phi_p) and cos(3*phi_q).
""")

# Verify
np.random.seed(42)
print("Verifying the phi_r formula:")
for trial in range(5):
    a_p = np.random.exponential(2) + 0.5
    a_q = np.random.exponential(2) + 0.5
    phi_p = np.random.uniform(0.1, np.pi/3 - 0.1)
    phi_q = np.random.uniform(0.1, np.pi/3 - 0.1)

    sigma_p = -a_p
    sigma_q = -a_q
    tau_p = -2*a_p**1.5*np.cos(3*phi_p) / (3*np.sqrt(3))
    tau_q = -2*a_q**1.5*np.cos(3*phi_q) / (3*np.sqrt(3))

    sigma_r = sigma_p + sigma_q
    tau_r = tau_p + tau_q
    a_r = -sigma_r

    disc = -4*sigma_r**3 - 27*tau_r**2
    if disc <= 0:
        continue

    cos3phi_r = -3*np.sqrt(3)*tau_r / (2*a_r**1.5)
    if abs(cos3phi_r) > 1:
        continue
    phi_r = np.arccos(cos3phi_r) / 3

    # Check formula
    cos3phi_r_formula = (a_p**1.5*np.cos(3*phi_p) + a_q**1.5*np.cos(3*phi_q)) / a_r**1.5
    print(f"  Trial {trial}: cos(3*phi_r) = {cos3phi_r:.8f}, formula = {cos3phi_r_formula:.8f}, match: {abs(cos3phi_r - cos3phi_r_formula) < 1e-8}")

# =============================================================
# Part 4: Reduced conjecture for n=3
# =============================================================
print("\n\n--- Part 4: REDUCED CONJECTURE FOR n=3 ---\n")

print("""
REDUCED CONJECTURE (n=3):
  F(phi_r) <= (a_p+a_q) / (a_p/F(phi_p) + a_q/F(phi_q))

where:
  F(phi) = csc^2(phi) + csc^2(phi+pi/3) + csc^2(phi+2*pi/3)
  cos(3*phi_r) = [a_p^{3/2}*cos(3*phi_p) + a_q^{3/2}*cos(3*phi_q)] / (a_p+a_q)^{3/2}
  a_p, a_q > 0, phi_p, phi_q, phi_r in (0, pi/3)

SIMPLIFICATION: Set psi = 3*phi so psi in (0, pi).
  cos(psi_r) = [a_p^{3/2}*cos(psi_p) + a_q^{3/2}*cos(psi_q)] / (a_p+a_q)^{3/2}

Note the CONTRACTION: a_p^{3/2} + a_q^{3/2} < (a_p+a_q)^{3/2} (strict concavity of x^{3/2}).
So the weighted average that gives cos(psi_r) uses weights summing to
LESS than 1. The remainder (1 - (a_p^{3/2}+a_q^{3/2})/(a_p+a_q)^{3/2}) is >0.

This means: cos(psi_r) is CLOSER TO ZERO than the weighted avg of cos(psi_p), cos(psi_q).
Equivalently: psi_r is CLOSER TO pi/2 (symmetric case) than the weighted average.
This is why Q_r decreases!
""")

# Verify: check that cos(psi_r) is indeed contracted toward 0
np.random.seed(42)
print("Verifying contraction of cos(psi_r) toward 0:")
for trial in range(10):
    a_p = np.random.exponential(2) + 0.5
    a_q = np.random.exponential(2) + 0.5
    psi_p = np.random.uniform(0.1, np.pi - 0.1)
    psi_q = np.random.uniform(0.1, np.pi - 0.1)

    gamma = (a_p + a_q)**1.5
    cos_psi_r = (a_p**1.5 * np.cos(psi_p) + a_q**1.5 * np.cos(psi_q)) / gamma
    w_avg = (a_p**1.5 * np.cos(psi_p) + a_q**1.5 * np.cos(psi_q)) / (a_p**1.5 + a_q**1.5)

    contraction = (a_p**1.5 + a_q**1.5) / gamma  # < 1
    print(f"  contraction={contraction:.4f}, |cos(psi_r)|={abs(cos_psi_r):.4f}, |weighted_avg|={abs(w_avg):.4f}, contracted: {abs(cos_psi_r) <= abs(w_avg) + 1e-10}")

# =============================================================
# Part 5: F as function of psi = 3*phi
# =============================================================
print("\n\n--- Part 5: F as function of psi ---\n")

# F(phi) = csc^2(phi) + csc^2(phi+pi/3) + csc^2(phi+2*pi/3)
# = csc^2(psi/3) + csc^2(psi/3+pi/3) + csc^2(psi/3+2*pi/3)

# G(psi) = F(psi/3) = csc^2(psi/3) + csc^2(psi/3+pi/3) + csc^2(psi/3+2*pi/3)

# G(psi) depends on psi = 3*phi, and cos(psi) = cos(3*phi) = ...
# cos(psi) is what determines the root structure.

# Key: G is symmetric about psi = pi/2 (since phi = pi/6 gives equal gaps)
# G(pi - psi) = G(psi) by the symmetry phi -> pi/3 - phi.

# Let's plot G(psi) vs cos(psi)
psi_vals = np.linspace(0.01, np.pi - 0.01, 200)
G_vals = []
for psi in psi_vals:
    phi = psi/3
    G = 1/np.sin(phi)**2 + 1/np.sin(phi+np.pi/3)**2 + 1/np.sin(phi+2*np.pi/3)**2
    G_vals.append(G)
G_vals = np.array(G_vals)

print("G(psi) vs cos(psi):")
for i in range(0, len(psi_vals), 20):
    print(f"  cos(psi)={np.cos(psi_vals[i]):.4f}, G={G_vals[i]:.4f}")

# Is G convex/concave as a function of cos(psi)?
cos_vals = np.cos(psi_vals)
# Interpolate G as function of u = cos(psi)
# Since G(psi) = G(pi - psi), and cos(psi) is decreasing:
# G is an even function of cos(psi) centered at cos(psi) = 0

u_vals = np.cos(psi_vals)  # from 1 to -1
# G(u) for u = cos(psi)
# Check: is G convex in u?
# d^2G/du^2
du = np.diff(u_vals)
dG = np.diff(G_vals)
dG_du = dG / du
d2G = np.diff(dG_du) / du[:-1]

print(f"\nIs G convex in u = cos(psi)?")
print(f"  min(d^2G/du^2) = {np.min(d2G):.4f}")
print(f"  G is {'CONVEX' if np.min(d2G) >= -1e-2 else 'NOT CONVEX'} in cos(psi)")

# Is G convex in cos^2(psi)?
rho_vals = np.cos(psi_vals)**2  # this is the "rho" parameter
# But rho values are not monotone... use the first half
half = len(psi_vals)//2
psi_half = psi_vals[:half]
rho_half = np.cos(psi_half)**2  # decreasing from 1 to ~0
G_half = G_vals[:half]

# Sort by rho
sorted_idx = np.argsort(rho_half)
rho_sorted = rho_half[sorted_idx]
G_sorted = G_half[sorted_idx]

drho = np.diff(rho_sorted)
dG_drho = np.diff(G_sorted) / drho
d2G_rho = np.diff(dG_drho) / drho[:-1]

print(f"\nIs G convex in rho = cos^2(psi)?")
print(f"  min(d^2G/drho^2) = {np.min(d2G_rho):.4f}")
print(f"  G is {'CONVEX' if np.min(d2G_rho) >= -1e-2 else 'NOT CONVEX'} in rho")

# =============================================================
# Part 6: Explicit computation: is G convex in u = cos(psi)?
# =============================================================
print("\n\n--- Part 6: Symbolic G(u) ---\n")

psi_sym = sp.Symbol('psi', positive=True)
phi_sym = psi_sym / 3
G_sym = 1/sp.sin(phi_sym)**2 + 1/sp.sin(phi_sym + sp.pi/3)**2 + 1/sp.sin(phi_sym + 2*sp.pi/3)**2

# Substitute u = cos(psi) = cos(3*phi)
# This requires expressing sin(phi), sin(phi+pi/3), etc in terms of u = cos(3*phi)
# Using multiple angle formula: cos(3x) = 4cos^3(x) - 3cos(x)

# Let c = cos(phi). Then cos(3*phi) = 4c^3 - 3c = u.
# sin^2(phi) = 1 - c^2
# sin^2(phi+pi/3) = sin^2(phi)*cos^2(pi/3) + 2*sin(phi)*cos(phi)*sin(pi/3)*cos(pi/3) + cos^2(phi)*sin^2(pi/3)
# = (1-c^2)/4 + 2*sqrt(1-c^2)*c*sqrt(3)/4 + 3c^2/4
# = (1-c^2+3c^2)/4 + sqrt(3)*c*sqrt(1-c^2)/2
# = (1+2c^2)/4 + sqrt(3)*c*sqrt(1-c^2)/2

# This is getting messy. Let me use a symbolic computation.

c = sp.Symbol('c', real=True)  # cos(phi)
s = sp.sqrt(1 - c**2)  # sin(phi)

# sin(phi+pi/3) = sin(phi)cos(pi/3) + cos(phi)sin(pi/3) = s/2 + c*sqrt(3)/2
sp1 = s/2 + c*sp.sqrt(3)/2

# sin(phi+2*pi/3) = sin(phi)cos(2pi/3) + cos(phi)sin(2pi/3) = -s/2 + c*sqrt(3)/2
sp2 = -s/2 + c*sp.sqrt(3)/2

G_c = 1/s**2 + 1/sp1**2 + 1/sp2**2
G_c_simplified = sp.simplify(G_c)
print(f"G(c) = {G_c_simplified}")

# Now u = 4c^3 - 3c. We want G as a function of u.
# But c is a function of u... this is the Chebyshev relation.

# For numerical purposes, just check convexity of G as function of u
# via the formula G(c(u)) where c satisfies 4c^3 - 3c = u

# Actually, for the PROOF, what we need is:
# Given cos(psi_r) = (alpha*cos(psi_p) + beta*cos(psi_q))/gamma
# where gamma >= alpha + beta (contraction),
# show F(phi_r) <= weighted_harmonic_mean(F(phi_p), F(phi_q); a_p, a_q)

# =============================================================
# Part 7: Direct numerical proof approach for n=3
# =============================================================
print("\n\n--- Part 7: Exhaustive numerical verification for n=3 ---\n")

# For n=3, the conjecture has 4 free parameters (a_p, a_q, phi_p, phi_q)
# after centering and using the (a, phi) parametrization.
# But by scale invariance, we can set a_p + a_q = 1.
# So 3 free parameters: a_p in (0,1), phi_p in (0,pi/3), phi_q in (0,pi/3).

np.random.seed(42)
violations = 0
valid = 0
min_gap = float('inf')

for trial in range(100000):
    a_p = np.random.uniform(0.01, 0.99)
    a_q = 1 - a_p
    phi_p = np.random.uniform(0.01, np.pi/3 - 0.01)
    phi_q = np.random.uniform(0.01, np.pi/3 - 0.01)

    # Compute tau
    tau_p = -2*a_p**1.5*np.cos(3*phi_p) / (3*np.sqrt(3))
    tau_q = -2*a_q**1.5*np.cos(3*phi_q) / (3*np.sqrt(3))

    sigma_p = -a_p
    sigma_q = -a_q
    sigma_r = sigma_p + sigma_q
    tau_r = tau_p + tau_q
    a_r = a_p + a_q  # = 1

    # Check domain
    disc_r = -4*sigma_r**3 - 27*tau_r**2
    if disc_r <= 0:
        continue

    # Compute phi_r
    cos3phi_r = -3*np.sqrt(3)*tau_r / (2*a_r**1.5)
    if abs(cos3phi_r) > 1 - 1e-10:
        continue
    phi_r = np.arccos(np.clip(cos3phi_r, -1, 1)) / 3

    # Compute F values
    def F_func(phi):
        return 1/np.sin(phi)**2 + 1/np.sin(phi+np.pi/3)**2 + 1/np.sin(phi+2*np.pi/3)**2

    F_p = F_func(phi_p)
    F_q = F_func(phi_q)
    F_r = F_func(phi_r)

    # 1/Sm2 = 4a/F
    inv_sm2_p = 4*a_p / F_p
    inv_sm2_q = 4*a_q / F_q
    inv_sm2_r = 4*a_r / F_r

    gap = inv_sm2_r - inv_sm2_p - inv_sm2_q
    valid += 1
    min_gap = min(min_gap, gap)

    if gap < -1e-8:
        violations += 1

print(f"n=3 exhaustive test: {violations} violations in {valid} valid trials")
print(f"min gap = {min_gap:.6e}")

# =============================================================
# Part 8: When is the gap smallest?
# =============================================================
print("\n\n--- Part 8: Near-equality cases ---\n")

np.random.seed(42)
near_equality = []

for trial in range(200000):
    a_p = np.random.uniform(0.01, 0.99)
    a_q = 1 - a_p
    phi_p = np.random.uniform(0.01, np.pi/3 - 0.01)
    phi_q = np.random.uniform(0.01, np.pi/3 - 0.01)

    tau_p = -2*a_p**1.5*np.cos(3*phi_p) / (3*np.sqrt(3))
    tau_q = -2*a_q**1.5*np.cos(3*phi_q) / (3*np.sqrt(3))

    sigma_r = -a_p - a_q
    tau_r = tau_p + tau_q
    a_r = a_p + a_q

    disc_r = -4*sigma_r**3 - 27*tau_r**2
    if disc_r <= 0:
        continue

    cos3phi_r = -3*np.sqrt(3)*tau_r / (2*a_r**1.5)
    if abs(cos3phi_r) > 1 - 1e-10:
        continue
    phi_r = np.arccos(np.clip(cos3phi_r, -1, 1)) / 3

    def F_func(phi):
        return 1/np.sin(phi)**2 + 1/np.sin(phi+np.pi/3)**2 + 1/np.sin(phi+2*np.pi/3)**2

    F_p = F_func(phi_p)
    F_q = F_func(phi_q)
    F_r = F_func(phi_r)

    inv_sm2_p = 4*a_p / F_p
    inv_sm2_q = 4*a_q / F_q
    inv_sm2_r = 4*a_r / F_r

    gap = inv_sm2_r - inv_sm2_p - inv_sm2_q

    if gap < 0.001:
        near_equality.append((gap, a_p, phi_p, phi_q))

near_equality.sort()
print("Smallest gaps:")
for gap, a_p, phi_p, phi_q in near_equality[:10]:
    print(f"  gap={gap:.6e}, a_p={a_p:.4f}, phi_p={phi_p:.4f}={phi_p*3/np.pi:.4f}*pi/3, phi_q={phi_q:.4f}={phi_q*3/np.pi:.4f}*pi/3")

print("""
Near-equality happens when:
- a_p and a_q are very unequal (one dominates)
- OR when phi_p, phi_q are both close to pi/6 (equal gaps, symmetric)
""")

print("\nDone.")
