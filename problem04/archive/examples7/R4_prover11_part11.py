"""
PROVER-11 Part 11: Symmetric case (s=t, k3p=-k3q) numerator analysis.

For s = t = K, k3p = -k3q = m, k4p = k4q = w:
The gap numerator is:
  N = 256*K^6*m^2 + 48*K^5*w^2 - 352*K^4*m^2*w + 128*K^3*m^4
    + 8*K^3*w^3 + 48*K^2*m^2*w^2 - 16*K*m^4*w - K*w^4 - 2*m^2*w^3

Domain: K > 0, w < 4K^2, 2m^2 < 4K^3 + K*w
Denominator: 24*K*(w-8K^2)*(w-4K^2)*(4K^3+Kw-2m^2) > 0 when all factors correct sign.

Let me scale: m = K^{3/2}*mu, w = K^2*nu. Then N = K^6 * N_reduced.
"""
from sympy import (symbols, expand, factor, collect, cancel,
                   Rational, simplify, S, Poly)
import numpy as np

K = symbols('K', positive=True)
m, w = symbols('m w', real=True)
mu, nu = symbols('mu nu', real=True)

N = (256*K**6*m**2 + 48*K**5*w**2 - 352*K**4*m**2*w + 128*K**3*m**4
     + 8*K**3*w**3 + 48*K**2*m**2*w**2 - 16*K*m**4*w - K*w**4 - 2*m**2*w**3)

# Substitute m = K^{3/2}*mu, w = K^2*nu
N_scaled = expand(N.subs({m: K**S(3,2)*mu, w: K**2*nu}))
# Factor out K^9 (since N has cumulant weight ... let me check)
# m^2 contributes K^3, K^6*m^2 = K^9.
# Actually let me just divide by appropriate power.
# All terms: K^6*m^2 = K^6*(K^3*mu^2) = K^9*mu^2
# K^5*w^2 = K^5*K^4*nu^2 = K^9*nu^2. ✓
# So N = K^9 * Nred
Nred = simplify(N_scaled / K**9)
print(f"N_reduced(mu, nu) = {expand(Nred)}")

# Domain in scaled vars: nu < 4, 2*mu^2 < 4 + nu
# i.e., mu^2 < (4+nu)/2

# Factor
print(f"\nFactored: {factor(Nred)}")

# Collect by mu
Nred_exp = expand(Nred)
Nred_poly = Poly(Nred_exp, mu)
print(f"\nAs polynomial in mu:")
for monom, coeff in sorted(Nred_poly.as_dict().items()):
    print(f"  mu^{monom[0]}: {factor(coeff)}")

# Collect by nu
Nred_poly_nu = Poly(Nred_exp, nu)
print(f"\nAs polynomial in nu:")
for monom, coeff in sorted(Nred_poly_nu.as_dict().items()):
    print(f"  nu^{monom[0]}: {factor(coeff)}")

# Try SOS decomposition
# N = 256*mu^2 + 48*nu^2 - 352*mu^2*nu + 128*mu^4 + 8*nu^3 + 48*mu^2*nu^2
#   - 16*mu^4*nu - nu^4 - 2*mu^2*nu^3
# Group: terms with mu^4: 128*mu^4 - 16*mu^4*nu = 16*mu^4*(8-nu) [positive since nu<4]
# Terms with mu^2: 256*mu^2 - 352*mu^2*nu + 48*mu^2*nu^2 - 2*mu^2*nu^3
#                = mu^2*(256 - 352*nu + 48*nu^2 - 2*nu^3)
# Terms without mu: 48*nu^2 + 8*nu^3 - nu^4

print("\n\nGrouping:")
print(f"  mu^4 terms: 16*mu^4*(8 - nu) [positive for nu < 4 < 8]")

coeff_mu2 = 256 - 352*nu + 48*nu**2 - 2*nu**3
print(f"  mu^2 coefficient: {coeff_mu2}")
print(f"    = {factor(coeff_mu2)}")

coeff_mu0 = 48*nu**2 + 8*nu**3 - nu**4
print(f"  mu^0 terms: {coeff_mu0}")
print(f"    = {factor(coeff_mu0)}")

# Check sign of mu^2 coefficient on nu < 4:
# 256 - 352*nu + 48*nu^2 - 2*nu^3
# = -2*(nu^3 - 24*nu^2 + 176*nu - 128)
# Roots of nu^3 - 24*nu^2 + 176*nu - 128 = 0?
from sympy import solve, Rational as R
nu_sym = symbols('nu_sym')
roots = solve(nu_sym**3 - 24*nu_sym**2 + 176*nu_sym - 128, nu_sym)
print(f"\n  Roots of cubic: {[float(r) for r in roots]}")
# If all roots > 4, then on nu < 4 the cubic is negative, so coeff is positive.

# Check mu^0 terms: nu^2*(48 + 8*nu - nu^2) = nu^2*(8+nu)(6-... wait
coeff_mu0_fac = factor(coeff_mu0)
print(f"\n  mu^0 factored: {coeff_mu0_fac}")

# Numerical check: is mu^2 coefficient positive on (-4, 4)?
import numpy as np
nu_test = np.linspace(-3.99, 3.99, 10000)
coeff_vals = 256 - 352*nu_test + 48*nu_test**2 - 2*nu_test**3
print(f"\n  mu^2 coeff: min={coeff_vals.min():.4f} at nu={nu_test[np.argmin(coeff_vals)]:.4f}")
print(f"              max={coeff_vals.max():.4f}")

# mu^0 coefficient
mu0_vals = 48*nu_test**2 + 8*nu_test**3 - nu_test**4
print(f"\n  mu^0 coeff: min={mu0_vals.min():.4f} at nu={nu_test[np.argmin(mu0_vals)]:.4f}")
