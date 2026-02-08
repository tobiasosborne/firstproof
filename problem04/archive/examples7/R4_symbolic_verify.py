"""
PROVER-9: Symbolic verification and proof attempt for R_4 superadditivity.

Step 1: Verify the R_4 formula from first principles.
Step 2: Check the decomposition 1/Phi_4 = (1/12)*k2 + R_4(k2,k3,k4).
Step 3: Attempt analytical proof of R_4 superadditivity.
"""
from sympy import (symbols, Rational, simplify, factor, expand, cancel,
                   collect, sqrt, pprint, latex, Matrix, det, S, Poly,
                   together, numer, denom, degree, Symbol, oo, limit,
                   series, diff, solve, groebner, resultant, discriminant,
                   Expr)
from sympy import Function
import sympy

print("=" * 70)
print("STEP 1: Verify R_4 formula from elementary symmetric polynomials")
print("=" * 70)

# For a centered degree-4 polynomial: x^4 + a2*x^2 + a3*x + a4
# where a1 = 0 (centered), a2 = e2, a3 = -e3, a4 = e4
# (Vieta: p(x) = x^4 - e1*x^3 + e2*x^2 - e3*x + e4)

e2, e3, e4 = symbols('e2 e3 e4', real=True)

# The discriminant of x^4 + e2*x^2 - e3*x + e4 (centered monic quartic):
# disc_4 = 256*e4^3 - 128*e2^2*e4^2 + 144*e2*e3^2*e4 - 27*e3^4 + 16*e2^4*e4 - 4*e2^3*e3^2
#
# But actually for p(x) = x^4 + a*x^2 + b*x + c, the discriminant is:
# D = 256c^3 - 128a^2*c^2 + 144a*b^2*c - 27b^4 + 16a^4*c - 4a^3*b^2
# with a = e2, b = -e3, c = e4

a, b, c = e2, -e3, e4
disc_quartic = 256*c**3 - 128*a**2*c**2 + 144*a*b**2*c - 27*b**4 + 16*a**4*c - 4*a**3*b**2
disc_quartic = expand(disc_quartic.subs([(a, e2), (b, -e3), (c, e4)]))
print("\nDisc(x^4 + e2*x^2 - e3*x + e4) =")
print(f"  {disc_quartic}")

# N_4 = Phi_4 * disc (the numerator)
# From the problem statement:
N4_given = -8*e2**5 - 36*e2**2*e3**2 - 64*e2**3*e4 - 432*e3**2*e4 + 384*e2*e4**2

print(f"\nN_4 (given) = {N4_given}")

# So Phi_4 = N_4 / disc_4
# And 1/Phi_4 = disc_4 / N_4

# Now convert to cumulants:
# k2 = -2*e2/3, k3 = 2*e3, k4 = -32*e4/3 + 8*e2^2/9
# So: e2 = -3*k2/2, e3 = k3/2, e4 = -3*k4/32 + 3*k2^2/16

k2, k3, k4 = symbols('k2 k3 k4', real=True)

e2_expr = Rational(-3, 2) * k2
e3_expr = k3 / 2
e4_expr = Rational(-3, 32) * k4 + Rational(3, 16) * k2**2

print("\n\nCumulant -> elementary symmetric substitutions:")
print(f"  e2 = {e2_expr}")
print(f"  e3 = {e3_expr}")
print(f"  e4 = {e4_expr}")

subs_dict = {e2: e2_expr, e3: e3_expr, e4: e4_expr}

# Compute disc_4 in cumulants
disc_in_k = expand(disc_quartic.subs(subs_dict))
print(f"\ndisc_4 in cumulants = {disc_in_k}")
disc_in_k_factored = factor(disc_in_k)
print(f"disc_4 factored = {disc_in_k_factored}")

# Compute N_4 in cumulants
N4_in_k = expand(N4_given.subs(subs_dict))
print(f"\nN_4 in cumulants = {N4_in_k}")
N4_in_k_factored = factor(N4_in_k)
print(f"N_4 factored = {N4_in_k_factored}")

# Compute 1/Phi_4 = disc_4/N_4 in cumulants
print("\n\n" + "=" * 70)
print("STEP 2: Verify 1/Phi_4 = k2/12 + R_4")
print("=" * 70)

# 1/Phi_4 = disc_4/N_4
# We want to check: 1/Phi_4 = k2/12 + R_4
# So R_4 = disc_4/N_4 - k2/12

inv_phi4 = cancel(disc_in_k / N4_in_k)
print(f"\n1/Phi_4 = {inv_phi4}")

# Extract R_4 = 1/Phi_4 - k2/12
R4 = simplify(inv_phi4 - k2/12)
R4 = cancel(R4)
print(f"\nR_4 = 1/Phi_4 - k2/12 = {R4}")

R4_num = numer(R4)
R4_den = denom(R4)
print(f"\nR_4 numerator = {expand(R4_num)}")
print(f"R_4 denominator = {expand(R4_den)}")

# Factor them
R4_num_f = factor(R4_num)
R4_den_f = factor(R4_den)
print(f"\nR_4 numerator factored = {R4_num_f}")
print(f"R_4 denominator factored = {R4_den_f}")

# Compare with the given formula
R4_num_given = -16*k2**3*k3**2 - 4*k2**2*k4**2 + 20*k2*k3**2*k4 - 8*k3**4 - k4**3
R4_den_given = 24*(16*k2**5 - 8*k2**2*k3**2 - k2*k4**2 + 2*k3**2*k4)

# Check if our computed formula matches the given one
diff_num = simplify(R4_num * R4_den_given - R4_den * R4_num_given)
print(f"\nNumerator cross-check (should be 0): {simplify(diff_num)}")

if simplify(diff_num) == 0:
    print("  VERIFIED: R_4 formula matches!")
else:
    print("  MISMATCH: R_4 formula does not match, investigating...")
    # Let's check what our formula actually is
    print(f"  Our R4_num = {expand(R4_num)}")
    print(f"  Given R4_num = {expand(R4_num_given)}")
    print(f"  Our R4_den = {expand(R4_den)}")
    print(f"  Given R4_den = {expand(R4_den_given)}")
    # Try to find the ratio
    ratio = cancel(R4_num * R4_den_given / (R4_den * R4_num_given))
    print(f"  Ratio = {ratio}")
