#!/usr/bin/env python3
"""
Find the RIGHT contour integral for <h,alpha>.

KEY INSIGHT: We proved sum_k h_k*(nu_k - lambda_k) >= 0 via the contour
integral of [nG_r(z)]^2 * [z - omega_1(z)].

Now we need sum_k h_k * (H_p(lambda_k) - H_r(nu_k)) >= 0.

IDEA: Replace the "position" displacement (nu_k - lambda_k) with the
"Hilbert transform" displacement (H_p(lambda_k) - H_r(nu_k)).

What function, evaluated at nu_k, gives H_p(lambda_k) - H_r(nu_k)?

Answer: Define Gamma(z) = nG_p(omega_1(z)) - nG_r(z) + 1/(z-omega_1(z))
Wait, this isn't right because nG_p(omega_1(z)) = nG_r(z).

Let's think more carefully. We need a function g(z) such that:
g(nu_k) = H_p(lambda_k) - H_r(nu_k) = alpha_k = omega_1''(nu_k)/2

So g(z) = omega_1''(z)/2 would do it, and we want
<h,alpha> = sum_k h_k * g(nu_k) = (1/2) * sum_k h_k * omega_1''(nu_k)

The question is: can we express this as a contour integral whose sign
can be determined?

Using [nG_r]^2 as the kernel (with Res = 2*h_k at nu_k):
sum_k h_k * g(nu_k) = (1/2) * sum_k Res_{nu_k}[[nG_r(z)]^2 * g(z)]
BUT the residue of [nG_r]^2 * g at nu_k is NOT 2*h_k*g(nu_k) in general!
It's: Res = 2*h_k*g(nu_k) + g'(nu_k) (from the expansion of 1/(z-nu_k)^2*g(z))

So: sum_k [2*h_k*g(nu_k) + g'(nu_k)] = (1/2pi i) oint [nG_r]^2 * g dz + ...

Actually, Res_{nu_k} [f(z)*g(z)] where f has a double pole:
f(z) = a_{-2}/(z-nu_k)^2 + a_{-1}/(z-nu_k) + ...
g(z) = g_0 + g_1*(z-nu_k) + ...
f*g = a_{-2}*g_0/(z-nu_k)^2 + [a_{-2}*g_1 + a_{-1}*g_0]/(z-nu_k) + ...
Res = a_{-2}*g_1 + a_{-1}*g_0

For [nG_r]^2: a_{-2} = 1, a_{-1} = 2*h_k
For g(z) = omega_1''(z)/2: g_0 = alpha_k, g_1 = omega_1'''(nu_k)/4

So Res = omega_1'''(nu_k)/4 + 2*h_k*alpha_k

Hmm, the omega_1''' term complicates things.

BETTER IDEA: Use [nG_r(z)]^2 * Phi(z) where Phi has a specific structure.

OR: use the IDENTITY omega_1'(nu_k) = 1 to do integration by parts.
"""

import numpy as np
from math import factorial

np.random.seed(42)

def boxplus_n(p_coeffs, q_coeffs, n):
    a, b = p_coeffs, q_coeffs
    r = np.zeros(n + 1)
    r[0] = 1.0
    for k in range(1, n + 1):
        c_k = 0.0
        for i in range(0, k + 1):
            j = k - i
            if i <= n and j <= n:
                coeff = (factorial(n - i) * factorial(n - j)) / (factorial(n) * factorial(n - k))
                c_k += coeff * a[i] * b[j]
        r[k] = c_k
    return r

def H_values(roots):
    n = len(roots)
    H = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i != j:
                H[i] += 1.0 / (roots[i] - roots[j])
    return H

def generate_random_poly(n, spread=5.0):
    roots = np.sort(np.random.uniform(-spread, spread, n))
    for i in range(1, n):
        if roots[i] - roots[i-1] < 0.5:
            roots[i] = roots[i-1] + 0.5 + np.random.uniform(0, 0.5)
    return roots


print("=" * 70)
print("KEY APPROACH: Contour integral of [nG_r]^2 * [nG_r' / nG_r - nG_p' / nG_p at omega_1]")
print("=" * 70)

# Actually, let me try a much more direct approach.
#
# IDENTITY: <h,alpha> = <h,u> - ||h||^2
#           = sum_k H_r(nu_k) * H_p(lambda_k) - sum_k H_r(nu_k)^2
#
# We know sum_k H_r(nu_k)^2 = Phi_n(r).
# We need to evaluate sum_k H_r(nu_k) * H_p(lambda_k).
#
# KEY: H_r(nu_k) = Res_{z=nu_k}[(nG_r(z))^2 / 2] ... no.
# Res_{z=nu_k}[(nG_r(z))^2] = 2*H_r(nu_k).
# So sum_k H_r(nu_k)^2 = (1/2)*sum_k Res[nG_r^2]*H_r(nu_k)
#
# But H_r(nu_k) = value of (nG_r(z) - 1/(z-nu_k)) at z = nu_k
# This is not a simple function of z.
#
# HOWEVER: there's a slick identity.
# sum_k H_r(nu_k)^2 = -(1/2pi i) oint (r''/r)^2 / 4 ... etc
#
# Let's take a step back and think about what FUNCTION evaluates to
# H_p(lambda_k) when you plug in z = nu_k.
#
# At z = nu_k, omega_1(z) = lambda_k.
# H_p(lambda_k) is the regular part of nG_p at w = lambda_k.
# Since nG_p(w) = 1/(w-lambda_k) + H_p(lambda_k) + O(w-lambda_k),
# we have H_p(lambda_k) = lim_{w->lambda_k} [nG_p(w) - 1/(w-lambda_k)].
#
# Now consider the function:
# F(z) = nG_p(omega_1(z)) - 1/(omega_1(z) - omega_1(nu_k))   restricted to z != nu_k
# At z = nu_k: omega_1(z) -> lambda_k, so
# F(nu_k) = lim [nG_p(omega_1(z)) - 1/(omega_1(z)-lambda_k)]
# = H_p(lambda_k) (the regular part!)
#
# But nG_p(omega_1(z)) = nG_r(z), so:
# F(z) = nG_r(z) - 1/(omega_1(z) - lambda_k)
# At z = nu_k: 1/(omega_1(z)-lambda_k) = 1/((z-nu_k)*omega_1'(nu_k)+...) = 1/(z-nu_k) + ...
# (using omega_1'(nu_k) = 1)
# So F(z) = nG_r(z) - 1/(z-nu_k) + O(1) near z = nu_k
# F(nu_k) = [regular part of nG_r at nu_k] + [correction from omega_1 expansion]
# = H_r(nu_k) + [correction]
#
# More carefully:
# omega_1(z) - lambda_k = (z-nu_k) + alpha_k*(z-nu_k)^2 + ...
# 1/(omega_1(z)-lambda_k) = 1/((z-nu_k)(1+alpha_k*(z-nu_k)+...))
# = (1/(z-nu_k)) * (1 - alpha_k*(z-nu_k) + ...)
# = 1/(z-nu_k) - alpha_k + ...
#
# nG_r(z) = 1/(z-nu_k) + H_r(nu_k) + ...
#
# F(z) = [1/(z-nu_k) + H_r(nu_k) + ...] - [1/(z-nu_k) - alpha_k + ...]
# = H_r(nu_k) + alpha_k + ... = H_p(lambda_k) + ...
# So F(nu_k) = H_p(lambda_k) = u_k. Correct!
#
# This gives us:
# u_k = H_p(lambda_k) = [nG_r(z) - 1/(omega_1(z)-lambda_k)]|_{z=nu_k}
#
# Now: sum_k h_k * u_k = sum_k H_r(nu_k) * F(nu_k)
# where F is a function that's regular at all nu_k.
#
# Using the residue identity:
# sum_k h_k * F(nu_k) = sum_k Res_{nu_k}[(nG_r)^2/2 * ???]
# No, this picks up F'(nu_k) terms too.
#
# SIMPLER: Use
# sum_k F(nu_k) = sum_k Res_{nu_k}[nG_r(z)*F(z)]
# = (1/2pi i) oint nG_r(z)*F(z) dz - (residues at other poles)
#
# This gives sum u_k = sum H_p(lambda_k) = 0, which is trivially true.
#
# For sum_k h_k * F(nu_k), we need a kernel with residue h_k.
# We already discussed that Res_{nu_k}[(nG_r)^2] = 2*h_k but then
# Res[(nG_r)^2*F] picks up both h_k*F(nu_k) and F'(nu_k).
#
# ACTUALLY, there's a cleaner way:
# sum_k f(nu_k)*g(nu_k) for functions f,g regular at nu_k is hard to express
# as a single contour integral. BUT:
#
# sum_k h_k*u_k = sum_k H_r(nu_k)*H_p(lambda_k)
#
# Consider the function:
# nG_r(z) * nG_p(omega_1(z)) = [nG_r(z)]^2  (since nG_p(omega_1(z)) = nG_r(z))
#
# This is just [nG_r]^2, which we already know.
# sum Res_{nu_k}[(nG_r)^2] = 2*sum h_k = 0
# (contour integral of 1/z^2 at infinity vanishes)
#
# What about nG_r(z) * [nG_p(omega_1(z)) - 1/(omega_1(z) - w)] for some parameter w?

# Let me try yet another approach entirely.

print("\n" + "=" * 70)
print("APPROACH: Direct proof via the implicit function")
print("=" * 70)

# The subordination equation is: nG_r(z) = nG_p(omega_1(z))
# i.e., r'(z)/r(z) = p'(omega_1(z))/p(omega_1(z))
# i.e., r'(z)*p(omega_1(z)) = p'(omega_1(z))*r(z)
#
# Differentiate with respect to z:
# r''(z)*p(omega_1) + r'(z)*p'(omega_1)*omega_1' = p''(omega_1)*omega_1'*r(z) + p'(omega_1)*r'(z)
#
# Wait, this is getting messy. Let me use a cleaner notation.
# Let G = nG_r, and recall G(z) = sum 1/(z-nu_k).
# Let omega = omega_1.
#
# G(z) = nG_p(omega(z))   ... (*)
#
# The key quantities:
# h_k = H_r(nu_k) = [regular part of G at nu_k]
# u_k = H_p(lambda_k) = [regular part of nG_p at lambda_k]
# alpha_k = u_k - h_k = omega''(nu_k)/2
#
# <h,u> = sum_k h_k*u_k, ||h||^2 = sum_k h_k^2
# <h,alpha> = <h,u> - ||h||^2

# Can we write <h,u> as a contour integral?
# <h,u> = sum_k h_k * u_k
# = sum_k [regular part of G at nu_k] * [regular part of nG_p at omega(nu_k)]
#
# Near z = nu_k: G(z) = 1/(z-nu_k) + h_k + O(z-nu_k)
# Near z = nu_k: nG_p(omega(z)) = 1/(z-nu_k) + u_k + O(z-nu_k)
#   (because omega(z) ~ lambda_k + (z-nu_k) near nu_k, so 1/(omega-lambda_k) ~ 1/(z-nu_k))
#   Wait, we showed: nG_p(omega(z)) = G(z) = 1/(z-nu_k) + h_k + ...
#   So u_k must equal h_k + alpha_k, and the expansion of nG_p(omega(z)) = G(z) near nu_k
#   is 1/(z-nu_k) + h_k + ... NOT 1/(z-nu_k) + u_k + ...
#   Let me reconsider.
#
# nG_p(w) near w = lambda_k: nG_p(w) = 1/(w-lambda_k) + u_k + O(w-lambda_k)
# omega(z) near z = nu_k: omega(z) = lambda_k + (z-nu_k) + alpha_k*(z-nu_k)^2 + ...
# omega(z) - lambda_k = (z-nu_k)(1 + alpha_k*(z-nu_k) + ...)
# nG_p(omega(z)) = 1/(omega-lambda_k) + u_k + O(omega-lambda_k)
# = 1/((z-nu_k)(1+alpha_k(z-nu_k)+...)) + u_k + O(z-nu_k)
# = [1/(z-nu_k)](1 - alpha_k(z-nu_k) + ...) + u_k + O(z-nu_k)
# = 1/(z-nu_k) - alpha_k + u_k + O(z-nu_k)
# = 1/(z-nu_k) + h_k + O(z-nu_k)   (since u_k - alpha_k = h_k)
# Consistent! So G(z) = nG_r(z) = 1/(z-nu_k) + h_k + ...

# The PRODUCT G(z) * G(z) near nu_k:
# G(z)^2 = 1/(z-nu_k)^2 + 2h_k/(z-nu_k) + (h_k^2 + c) + ...
# Res at nu_k is 2h_k.
# sum Res = 2*sum h_k = 0 (by large-z behavior)

# Consider instead: G(z) * [d/dz(G(z))]
# = (1/2) d/dz [G(z)^2]
# This is an exact derivative, so its contour integral is 0.

# TRY: Look at the integral from a DIFFERENT angle.
# <h,alpha> = (1/2) sum_k h_k * omega''(nu_k)
#
# omega''(z) = -2*sum_j m_j/(z-p_j)^3  where p_j are poles of omega_1, m_j > 0
#
# h_k * omega''(nu_k) = -2*h_k*sum_j m_j/(nu_k-p_j)^3
#
# <h,alpha> = -sum_j m_j * sum_k h_k/(nu_k-p_j)^3
#           = -sum_j m_j * [-(1/4)*F''(p_j)]   (from our earlier identity)
#           = (1/4)*sum_j m_j * F''(p_j)
# where F(p) = r''(p)/r(p)  (WAIT this has poles at roots of r!)
#
# Hmm, this conflicts with our formula from the earlier script.
# Let me recheck the signs.
# alpha_k = omega''(nu_k)/2
# omega(z) = z + c + sum_j m_j/(z-p_j), so omega''(z) = sum_j 2*m_j/(z-p_j)^3
# alpha_k = sum_j m_j/(nu_k-p_j)^3
# <h,alpha> = sum_k h_k * sum_j m_j/(nu_k-p_j)^3
#           = sum_j m_j * sum_k h_k/(nu_k-p_j)^3
#           = sum_j m_j * A_j
# where A_j = sum_k h_k/(nu_k-p_j)^3

# Earlier we showed sum_k h_k/(nu_k-p) = -(1/2)*r''(p)/r(p) = -(1/2)*F(p)
# d/dp: sum_k h_k/(nu_k-p)^2 = (1/2)*F'(p)
# d/dp again: 2*sum_k h_k/(nu_k-p)^3 = (1/2)*F''(p)
# Wait: d/dp[1/(nu_k-p)^2] = 2/(nu_k-p)^3
# So d/dp[sum h_k/(nu_k-p)^2] = sum 2*h_k/(nu_k-p)^3 = 2*A_j
# And (1/2)*F''(p) = 2*A_j
# So A_j = (1/4)*F''(p)

# So <h,alpha> = sum_j m_j * (1/4)*F''(p_j) = (1/4)*sum_j m_j * F''(p_j)

# Now: F(p) = r''(p)/r(p) = 2*sum_k h_k/(p-nu_k) (as shown)
# F'(p) = -2*sum_k h_k/(p-nu_k)^2
# F''(p) = 4*sum_k h_k/(p-nu_k)^3

# So A_j = (1/4)*4*sum_k h_k/(p_j-nu_k)^3 = sum_k h_k/(p_j-nu_k)^3
# And since 1/(p_j-nu_k)^3 = -1/(nu_k-p_j)^3:
# A_j = -sum_k h_k/(nu_k-p_j)^3

# SIGN ERROR! Let me recheck.
# omega''(z) = sum_j 2*m_j/(z-p_j)^3
# At z = nu_k: omega''(nu_k) = sum_j 2*m_j/(nu_k-p_j)^3
# alpha_k = omega''(nu_k)/2 = sum_j m_j/(nu_k-p_j)^3

# <h,alpha> = sum_k h_k * sum_j m_j/(nu_k-p_j)^3
#           = sum_j m_j * [sum_k h_k/(nu_k-p_j)^3]

# Now use: sum_k h_k/(nu_k-p) = -(1/2)*F(p)
# d/dp: -sum_k h_k/(nu_k-p)^2 = -(1/2)*F'(p), so sum_k h_k/(nu_k-p)^2 = (1/2)*F'(p)
# d/dp: -2*sum_k h_k/(nu_k-p)^3 = (1/2)*F''(p), so sum_k h_k/(nu_k-p)^3 = -(1/4)*F''(p)

# YES: A_j = sum_k h_k/(nu_k-p_j)^3 = -(1/4)*F''(p_j)
# So <h,alpha> = sum_j m_j * [-(1/4)*F''(p_j)] = -(1/4)*sum_j m_j*F''(p_j)

# For <h,alpha> >= 0, we need sum_j m_j*F''(p_j) <= 0.

# Now F(p) = r''(p)/r(p) has poles at the roots nu_k of r.
# Between consecutive roots, r doesn't vanish, and F is a smooth function.
# The poles p_j of omega_1 interlace with the roots of r (they lie in the gaps
# between roots, or outside).

# Q: Is F concave (F'' <= 0) between the roots of r?
# F = r''/r, and between roots r has constant sign.
# This is a deep question about the convexity of r''/r.

# Let's compute F'' between roots numerically and check its sign.

print("\n" + "=" * 70)
print("Testing concavity of F(p) = r''/r between roots of r")
print("=" * 70)

for n in [3, 4, 5]:
    np.random.seed(42 + n)
    count_concave = 0
    count_total = 0
    for trial in range(100):
        r_roots = generate_random_poly(n)
        r_poly = np.poly(r_roots)
        r_d2 = np.polyder(np.polyder(r_poly))
        r_d3 = np.polyder(r_d2)
        r_d4 = np.polyder(r_d3)

        # Compute F'' at many points between roots
        for gap in range(n - 1):
            p_vals = np.linspace(r_roots[gap] + 0.01, r_roots[gap+1] - 0.01, 20)
            for p in p_vals:
                # F(p) = r''(p)/r(p)
                rp = np.polyval(r_poly, p)
                rp2 = np.polyval(r_d2, p)
                rp3 = np.polyval(r_d3, p)
                rp4 = np.polyval(r_d4, p)

                # F' = (r''' * r - r'' * r') / r^2
                rp1 = np.polyval(np.polyder(r_poly), p)
                Fp = (rp3 * rp - rp2 * rp1) / rp**2

                # F'' = d/dp [(r''' r - r'' r') / r^2]
                # numerically:
                h = 1e-6
                pp = p + h
                pm = p - h
                Fp_plus = (np.polyval(r_d3, pp)*np.polyval(r_poly, pp) -
                           np.polyval(r_d2, pp)*np.polyval(np.polyder(r_poly), pp)) / np.polyval(r_poly, pp)**2
                Fp_minus = (np.polyval(r_d3, pm)*np.polyval(r_poly, pm) -
                            np.polyval(r_d2, pm)*np.polyval(np.polyder(r_poly), pm)) / np.polyval(r_poly, pm)**2
                Fpp = (Fp_plus - Fp_minus) / (2*h)

                count_total += 1
                if Fpp <= 1e-4:
                    count_concave += 1

    print(f"  n={n}: F'' <= 0 in {count_concave}/{count_total} = {count_concave/count_total:.3f} "
          f"of points between consecutive roots")

# Also check outside
print("\nF'' outside roots:")
for n in [3, 4, 5]:
    np.random.seed(42 + n)
    count_concave = 0
    count_total = 0
    for trial in range(100):
        r_roots = generate_random_poly(n)
        r_poly = np.poly(r_roots)
        r_d2 = np.polyder(np.polyder(r_poly))
        r_d3 = np.polyder(r_d2)

        for p in np.linspace(r_roots[0] - 5, r_roots[0] - 0.1, 10):
            rp = np.polyval(r_poly, p)
            h = 1e-6
            pp = p + h
            pm = p - h
            Fp_plus = (np.polyval(r_d3, pp)*np.polyval(r_poly, pp) -
                       np.polyval(r_d2, pp)*np.polyval(np.polyder(r_poly), pp)) / np.polyval(r_poly, pp)**2
            Fp_minus = (np.polyval(r_d3, pm)*np.polyval(r_poly, pm) -
                        np.polyval(r_d2, pm)*np.polyval(np.polyder(r_poly), pm)) / np.polyval(r_poly, pm)**2
            Fpp = (Fp_plus - Fp_minus) / (2*h)
            count_total += 1
            if Fpp <= 1e-4:
                count_concave += 1

        for p in np.linspace(r_roots[-1] + 0.1, r_roots[-1] + 5, 10):
            rp = np.polyval(r_poly, p)
            h = 1e-6
            pp = p + h
            pm = p - h
            Fp_plus = (np.polyval(r_d3, pp)*np.polyval(r_poly, pp) -
                       np.polyval(r_d2, pp)*np.polyval(np.polyder(r_poly), pp)) / np.polyval(r_poly, pp)**2
            Fp_minus = (np.polyval(r_d3, pm)*np.polyval(r_poly, pm) -
                        np.polyval(r_d2, pm)*np.polyval(np.polyder(r_poly), pm)) / np.polyval(r_poly, pm)**2
            Fpp = (Fp_plus - Fp_minus) / (2*h)
            count_total += 1
            if Fpp <= 1e-4:
                count_concave += 1

    print(f"  n={n}: F'' <= 0 outside roots: {count_concave}/{count_total} = {count_concave/count_total:.3f}")


print("\n" + "=" * 70)
print("CONTOUR INTEGRAL for <h,alpha> via residues at poles of omega_1")
print("=" * 70)

# We established: <h,alpha> = -(1/4)*sum_j m_j*F''(p_j)
# where F(p) = r''(p)/r(p) and p_j are poles of omega_1 with residues m_j > 0.
#
# The integral of [nG_r(z)]^2 * [z - omega_1(z)] gave us <h,delta>.
# The integral of [nG_r(z)]^2 * [nG_r(z) - something related to omega]
# should give us <h,alpha>...
#
# KEY NEW IDEA: instead of looking at omega_1(z) - z (position shift),
# look at the LOG DERIVATIVE shift.
#
# Define: L(z) = d/dz log(r(z)/p(omega_1(z)))
# = r'(z)/r(z) - [p'(omega_1(z))*omega_1'(z)]/p(omega_1(z))
# = nG_r(z) - nG_p(omega_1(z))*omega_1'(z)
# = nG_r(z) - nG_r(z)*omega_1'(z)   (using nG_p(omega_1) = nG_r)
# = nG_r(z)*(1 - omega_1'(z))
# = -nG_r(z)*Psi'(z)  where Psi = omega_1 - id, Psi' = omega_1' - 1
#
# Near z = nu_k:
# nG_r(z) = 1/(z-nu_k) + h_k + ...
# Psi'(z) = 0 + Psi''(nu_k)*(z-nu_k) + ... = -2*alpha_k*(z-nu_k) + ...
# L(z) = -(1/(z-nu_k)+h_k+...)*(−2*alpha_k*(z-nu_k)+...)
# = 2*alpha_k + 2*h_k*alpha_k*(z-nu_k) + ...
# So L is REGULAR at nu_k with L(nu_k) = 2*alpha_k.
#
# Now: sum_k alpha_k = sum_k Res_{nu_k}[nG_r(z)] * alpha_k
# = sum_k Res_{nu_k}[nG_r(z) * L(z)/(2)] ... no, L is regular at nu_k.
# Actually: Res_{nu_k}[nG_r * L/2] = L(nu_k)/2 = alpha_k.
# So sum alpha_k = (1/2pi i) oint nG_r(z)*L(z)/2 dz + (poles at p_j)
# = (1/2)*sum Res_{nu_k} + (1/2)*sum Res_{p_j}
# Since sum alpha_k = 0, this is consistent.
#
# For <h,alpha>: we need sum_k h_k*alpha_k.
# Consider: (1/2pi i) oint [nG_r(z)]^2 * L(z)/2 dz
# = sum_k Res_{nu_k}[(nG_r)^2 * L/2] + sum_j Res_{p_j}[...]
#
# Res_{nu_k}[(nG_r)^2*L/2]:
# (nG_r)^2 = 1/(z-nu_k)^2 + 2h_k/(z-nu_k) + ...
# L/2 = alpha_k + h_k*alpha_k*(z-nu_k) + ... (I need the L'(nu_k)/2 coefficient)
# Actually L(z) = -nG_r(z)*Psi'(z)
# L(nu_k) = 2*alpha_k (computed above)
# L'(z) = -nG_r'*Psi' - nG_r*Psi''
# L'(nu_k) = -nG_r'(nu_k)*0 - nG_r(nu_k)*Psi''(nu_k)
# But nG_r has a pole at nu_k, so nG_r(nu_k) is not defined.
# The product nG_r*Psi'' near nu_k:
# = [1/(z-nu_k)+h_k+...]*(Psi''(nu_k) + ...) = Psi''(nu_k)/(z-nu_k) + h_k*Psi''(nu_k) + ...
# So L'(z) near nu_k has a pole! L is NOT regular at nu_k in the sense I thought.
# Wait, I computed L = -nG_r*Psi'. Near nu_k:
# nG_r ~ 1/(z-nu_k), Psi' ~ -2*alpha_k*(z-nu_k)
# Product ~ -2*alpha_k (regular). YES, L is regular.
# L'(z) = d/dz[-nG_r(z)*Psi'(z)]:
# = [-nG_r'(z)*Psi'(z)] + [-nG_r(z)*Psi''(z)]
# Near nu_k: nG_r'(z) ~ -1/(z-nu_k)^2 + ..., Psi'(z) ~ -2*alpha_k*(z-nu_k) + ...
# First term: (1/(z-nu_k)^2)*(-2*alpha_k*(z-nu_k)) = -2*alpha_k/(z-nu_k)  -- has a pole!
# Second term: -(1/(z-nu_k)+h_k)*Psi''(nu_k) + ... = -Psi''(nu_k)/(z-nu_k) + ...
#            = 2*alpha_k/(z-nu_k) + ...   (since Psi''(nu_k) = -omega''(nu_k) = -2*alpha_k)
# Wait: Psi = z - omega_1, so Psi'' = -omega_1''. omega_1'' = 2*alpha, so Psi'' = -2*alpha.
# Second term = -(1/(z-nu_k))*(-2*alpha_k) + ... = 2*alpha_k/(z-nu_k) + ...
# Sum: L'(z) = (-2*alpha_k + 2*alpha_k)/(z-nu_k) + ... = regular!
# So L'(nu_k) = ... some regular value involving higher-order terms.
#
# This is getting complicated. Let me just compute numerically.

# KEY RESULT: <h,alpha> = -(1/4)*sum_j m_j*F''(p_j)
# F'' can be positive or negative between roots.
# BUT: the ACTUAL poles p_j of omega_1 are NOT arbitrary points.
# They depend on the specific subordination function.

# Let me verify numerically that <h,alpha> = -(1/4)*sum_j m_j*F''(p_j)
# by actually finding the poles and residues of omega_1.

# For the Herglotz function omega_1(z) = z + c + sum_j m_j/(z-p_j):
# At z = infinity: omega_1(z) = z + c + O(1/z)
# The constant c and the poles p_j, residues m_j characterize omega_1.
#
# omega_1 has n poles (it's a degree-n rational function, since it comes from
# solving the degree-n equation r'(z)p(w) = p'(w)r(z) for w).
# Wait: omega_1 is degree (n-1) as a rational function of z.
# The equation nG_r(z) = nG_p(w) gives a relation where the RHS has n-1 roots in w.
# For z large, omega_1(z) ~ z, so omega_1 is a degree-1/degree-0 rational map at
# leading order, but with n-2 additional simple poles from the other branches.
# Actually omega_1 is a degree-n rational function with numerator and denominator
# both degree n, but with leading coefficient ratio 1.
# No: if omega(z) = z + c + sum_{j=1}^{n-1} m_j/(z-p_j), then n-1 poles.

# For n=3: omega_1 has 2 poles. Let me find them.

print("\nFinding poles and residues of omega_1 for n=3:")
n = 3
np.random.seed(100)
p_roots = generate_random_poly(n)
q_roots = generate_random_poly(n)
p_coeffs = np.poly(p_roots)
q_coeffs = np.poly(q_roots)
r_coeffs = boxplus_n(p_coeffs, q_coeffs, n)
r_roots_raw = np.roots(r_coeffs)
r_roots = np.sort(np.real(r_roots_raw))

print(f"p_roots = {p_roots}")
print(f"r_roots = {r_roots}")
print(f"H_p = {H_values(p_roots)}")
print(f"H_r = {H_values(r_roots)}")

# omega_1(z) solves r'(z)*p(w) = p'(w)*r(z)
# This is a polynomial in w of degree n, with leading coefficient r'(z) (for p(w) = w^n + ...)
# minus p'(w) * r(z) = n*w^{n-1}*r(z) + ...
# Combined: r'(z)*w^n - n*r(z)*w^{n-1} + ... = 0
# So it's degree n in w, with one root at w = omega_1(z).
# The OTHER roots correspond to the other branches of the subordination function
# (the ones that map C+ to C-).

# For the Herglotz function omega_1(z) with omega_1(z) ~ z at infinity:
# Write omega_1(z) = (a(z) * z + b(z)) / (c(z) * z + d(z)) ... no, simpler.
# omega_1(z) = z + c_0 + sum m_j/(z - p_j)
# = [z*(z-p_1)(z-p_2)...(z-p_{n-1}) + c_0*prod(z-p_j) + sum_j m_j*prod_{l!=j}(z-p_l)] / prod(z-p_j)
# The denominator is prod_{j=1}^{n-1} (z-p_j), degree n-1.
# The numerator has degree n (from the z*prod term).
# So omega_1 = (degree n) / (degree n-1), which is a degree-n rational function.

# To find p_j: look at where omega_1 blows up.
# From the equation r'(z)p(w) = p'(w)r(z), at a pole z = p_j of omega_1,
# w = omega_1(z) -> infinity.
# If w -> infinity: p(w) ~ w^n, p'(w) ~ n*w^{n-1}
# So r'(p_j)*w^n ~ n*r(p_j)*w^{n-1}, i.e., r'(p_j)*w ~ n*r(p_j)
# So w ~ n*r(p_j)/r'(p_j) which diverges only if r'(p_j) = 0.
# Wait: that means p_j are the CRITICAL POINTS of r, i.e., roots of r'.
# r has degree n, so r' has degree n-1 and has n-1 roots.
# These are exactly the n-1 poles of omega_1!

r_deriv = np.polyder(np.poly(r_roots))
critical_pts = np.sort(np.real(np.roots(r_deriv)))
print(f"\nCritical points of r (roots of r'): {critical_pts}")
print(f"These should be the poles of omega_1!")

# Verify: at z = critical point c_j of r, omega_1(c_j) should diverge.
# More precisely, near z = c_j:
# r'(z) = r''(c_j)*(z-c_j) + ..., r(z) = r(c_j) + r'(c_j)*(z-c_j) + ... = r(c_j) + O(z-c_j)^2
# Wait: r'(c_j) = 0, so r(z) = r(c_j) + (1/2)*r''(c_j)*(z-c_j)^2 + ...
# From the equation: r'(z)*p(omega) = p'(omega)*r(z)
# r''(c_j)*(z-c_j)*p(omega) ~ p'(omega)*r(c_j)
# If omega -> infinity: p(omega) ~ omega^n, p'(omega) ~ n*omega^{n-1}
# r''(c_j)*(z-c_j)*omega^n ~ n*r(c_j)*omega^{n-1}
# omega ~ n*r(c_j) / (r''(c_j)*(z-c_j))
# So omega_1(z) ~ n*r(c_j)/(r''(c_j)*(z-c_j)) near z = c_j

# This means the residue m_j of omega_1 at p_j = c_j is:
# m_j = -n*r(c_j)/r''(c_j)  (with the convention omega = m/(z-p) + ... = -m*(z-p)^{-1}...
# Hmm careful: omega_1(z) ~ n*r(c_j)/(r''(c_j)*(z-c_j)) means
# the leading term is [n*r(c_j)/r''(c_j)] * 1/(z-c_j)
# So m_j = n*r(c_j)/r''(c_j)
# For this to be positive (Herglotz property), we need r(c_j)/r''(c_j) > 0.

# Between consecutive roots nu_k and nu_{k+1}:
# r(c_j) has sign (-1)^{n-1-k} (since r changes sign at each root)
# Wait: r is monic of degree n. For z just above nu_1 (smallest root): r > 0 if n is odd
# Actually: r(z) = prod(z-nu_i), so for z > nu_n: r > 0.
# Between nu_{k} and nu_{k+1} (where k is 0-indexed from 1 to n):
# r has the sign (-1)^{n-k-1} if we count from the right.
# The critical point c_j lies between nu_j and nu_{j+1} (0-indexed: between j-th and (j+1)-th root)
# r''(c_j): at a local extremum, r'' has the opposite sign of r (for a local max) or same (local min)
# Actually: if r(c_j) > 0 (local max), then r''(c_j) < 0.
# If r(c_j) < 0 (local min), then r''(c_j) > 0.
# So r(c_j)*r''(c_j) < 0 at every critical point c_j!
# Therefore m_j = n*r(c_j)/r''(c_j) < 0 ... which contradicts m_j > 0!

# Wait, this can't be right. Let me reconsider.
# omega_1(z) = z + const + sum m_j/(z - p_j) where m_j > 0
# omega_1(z) near z = p_j: omega_1(z) ~ m_j/(z-p_j) for the singular part
# But we computed omega_1(z) ~ n*r(c_j)/(r''(c_j)*(z-c_j))
# If r(c_j)/r''(c_j) < 0, then omega_1 ~ [negative]/(z-c_j)
# This means m_j > 0 requires a sign flip... let me recheck.
#
# omega_1(z) = z + const + sum m_j/(z-p_j), and near p_j:
# omega_1(z) ~ m_j/(z-p_j) (the z + const + other terms are regular at p_j)
# But we computed the FULL omega_1(z) ~ n*r(c_j)/(r''(c_j)) * 1/(z-c_j)
# This includes the z term! So z + m_j/(z-c_j) ~ c_j + m_j/(z-c_j) + (z-c_j)
# For large m_j/(z-c_j), the dominant term is m_j/(z-c_j).
# So m_j = n*r(c_j)/r''(c_j).
# For the Herglotz property m_j > 0, we need r(c_j)/r''(c_j) > 0.
# But r(c_j)/r''(c_j) < 0 at every critical point!
# This is a CONTRADICTION!

# Unless omega_1 is NOT of the form z + const + sum m_j/(z-p_j).
# For a Herglotz function f: C+ -> C+, the representation is:
# f(z) = az + b + integral [1/(t-z) + t/(1+t^2)] dmu(t) where a >= 0, b real, mu positive measure
# For a RATIONAL Herglotz function with poles on the real line:
# f(z) = az + b + sum m_j/(p_j - z)  with a >= 0, m_j > 0
# Note: 1/(p_j - z) = -1/(z - p_j)
# So f(z) = az + b - sum m_j/(z - p_j)

# AH! I had the sign wrong!
# omega_1(z) = z + c - sum_j m_j/(z - p_j) with m_j > 0
# So omega_1(z) near p_j: omega_1(z) ~ -m_j/(z-p_j)
# And we need -m_j = n*r(c_j)/r''(c_j)
# So m_j = -n*r(c_j)/r''(c_j) > 0 since r(c_j)/r''(c_j) < 0.

# Now: omega_1''(z) = d^2/dz^2 [z + c - sum m_j/(z-p_j)]
# = -sum m_j * d^2/dz^2 [1/(z-p_j)]
# = -sum m_j * 2/(z-p_j)^3

# So alpha_k = omega_1''(nu_k)/2 = -sum_j m_j/(nu_k - p_j)^3

# And <h,alpha> = sum_k h_k * alpha_k = -sum_j m_j * sum_k h_k/(nu_k - p_j)^3
#              = -sum_j m_j * A_j
# where A_j = sum_k h_k/(nu_k - p_j)^3

# Using sum_k h_k/(nu_k-p)^3 = -(1/4)*F''(p):
# <h,alpha> = -sum_j m_j * [-(1/4)*F''(p_j)] = (1/4)*sum_j m_j*F''(p_j)

# For <h,alpha> >= 0: need sum_j m_j*F''(p_j) >= 0
# i.e., F'' >= 0 (CONVEXITY of F) at the poles p_j = critical points of r.

# Let's check this!
print("\n*** CORRECTED SIGN: omega_1(z) = z + c - sum m_j/(z-p_j) ***")
print("*** Poles are at CRITICAL POINTS of r ***")
print("*** m_j = -n*r(c_j)/r''(c_j) > 0 ***")
print()

for n in [3, 4, 5]:
    np.random.seed(42 + n)
    count_convex_all = 0
    count_total = 0
    for trial in range(200):
        r_roots = generate_random_poly(n)
        r_poly = np.poly(r_roots)
        r_d1 = np.polyder(r_poly)
        r_d2 = np.polyder(r_d1)
        r_d3 = np.polyder(r_d2)

        crit_pts = np.sort(np.real(np.roots(r_d1)))

        all_convex = True
        for c_j in crit_pts:
            # Compute F''(c_j) numerically
            h = 1e-6
            F_plus = np.polyval(r_d2, c_j+h) / np.polyval(r_poly, c_j+h)
            F_0 = np.polyval(r_d2, c_j) / np.polyval(r_poly, c_j)
            F_minus = np.polyval(r_d2, c_j-h) / np.polyval(r_poly, c_j-h)
            Fpp = (F_plus - 2*F_0 + F_minus) / h**2

            if Fpp < -1e-2:
                all_convex = False

        count_total += 1
        if all_convex:
            count_convex_all += 1

    print(f"  n={n}: F''(c_j) >= 0 at ALL critical points in {count_convex_all}/{count_total} cases")


print("\n\nMore detailed check:")
for n in [3, 4, 5, 6]:
    np.random.seed(42 + n)
    min_Fpp = float('inf')
    for trial in range(500):
        r_roots = generate_random_poly(n)
        r_poly = np.poly(r_roots)
        r_d1 = np.polyder(r_poly)
        r_d2 = np.polyder(r_d1)

        crit_pts = np.sort(np.real(np.roots(r_d1)))

        for c_j in crit_pts:
            h = 1e-6
            F_plus = np.polyval(r_d2, c_j+h) / np.polyval(r_poly, c_j+h)
            F_0 = np.polyval(r_d2, c_j) / np.polyval(r_poly, c_j)
            F_minus = np.polyval(r_d2, c_j-h) / np.polyval(r_poly, c_j-h)
            Fpp = (F_plus - 2*F_0 + F_minus) / h**2

            if Fpp < min_Fpp:
                min_Fpp = Fpp

    print(f"  n={n}: min F''(c_j) across 500 trials = {min_Fpp:.6f}")

# Also verify m_j > 0
print("\nVerifying m_j = -n*r(c_j)/r''(c_j) > 0:")
for n in [3, 4, 5]:
    np.random.seed(42 + n)
    count_positive = 0
    count_total = 0
    for trial in range(200):
        r_roots = generate_random_poly(n)
        r_poly = np.poly(r_roots)
        r_d1 = np.polyder(r_poly)
        r_d2 = np.polyder(r_d1)

        crit_pts = np.sort(np.real(np.roots(r_d1)))

        for c_j in crit_pts:
            m_j = -n * np.polyval(r_poly, c_j) / np.polyval(r_d2, c_j)
            count_total += 1
            if m_j > 0:
                count_positive += 1

    print(f"  n={n}: m_j > 0 in {count_positive}/{count_total} cases")
