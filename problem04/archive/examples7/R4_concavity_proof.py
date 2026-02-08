"""
PROVER-9: Attempt proof via perspective function / joint concavity.

KEY INSIGHT: For k3=0, R_4 = -k4^2 / (24*k2*(4*k2^2-k4))

This is of the form h(k2, k4) = -k4^2 / (24*(4*k2^3 - k2*k4))

Superadditivity means: h(s+t, u+v) >= h(s,u) + h(t,v)
i.e., -h is subadditive.

-h(k2,k4) = k4^2 / (24*(4*k2^3 - k2*k4))
           = k4^2 / (24*k2*(4*k2^2 - k4))

For f subadditive: f(x+y) <= f(x)+f(y).

A function is subadditive if it is concave and f(0) >= 0.
But -h(0,0) = 0, so IF -h is concave, then -h is subadditive. Done.

Let's check: is -h(k2, k4) = k4^2 / (24*k2*(4*k2^2 - k4)) CONCAVE?
On the domain k2 > 0, k4 < 4*k2^2.

Actually wait: for subadditivity from concavity, we need:
f(x+y) <= f(x) + f(y) where f(0)=0.
This follows because f is concave => f(x+y) = f(x+y) + f(0) <= 2f((x+y)/2 + 0/2)...
No, that doesn't work. Let me be more careful.

Concavity + f(0) = 0 gives:
f(lambda*x) >= lambda*f(x) for lambda in [0,1]
i.e., f is superadditive in the radial direction.

But f(x+y) >= f(x) + f(y) (superadditivity) does NOT follow from concavity alone.
In fact, concavity + f(0)=0 gives SUBadditivity: f(x+y) <= f(x)+f(y).

Hmm, but we want h superadditive, i.e., -h subadditive.
-h concave + -h(0)=0 => -h subadditive => h superadditive!

So the question is: is -h(k2, k4) = k4^2 / (24*k2*(4*k2^2-k4)) concave on the domain?
"""
import numpy as np
from sympy import symbols, Rational, simplify, factor, expand, cancel, Matrix, det, diff, Poly, hessian

k2, k4 = symbols('k2 k4', real=True)

# f(k2, k4) = k4^2 / (24*k2*(4*k2^2 - k4))
f = k4**2 / (24*k2*(4*k2**2 - k4))

print("f(k2, k4) = k4^2 / (24*k2*(4*k2^2 - k4))")
print("= -R_4|_{k3=0}")
print()
print("Checking if f is concave on {k2 > 0, k4 < 4*k2^2}...")
print()

# Compute Hessian
f_22 = diff(f, k2, k2)
f_24 = diff(f, k2, k4)
f_44 = diff(f, k4, k4)

H = Matrix([[f_22, f_24], [f_24, f_44]])

print("Hessian of f:")
print(f"  f_k2k2 = {simplify(f_22)}")
print(f"  f_k2k4 = {simplify(f_24)}")
print(f"  f_k4k4 = {simplify(f_44)}")

det_H = simplify(det(H))
print(f"\n  det(H) = {simplify(det_H)}")

# For concavity, we need the Hessian to be negative semidefinite:
# f_k2k2 <= 0 AND det(H) >= 0

# Let's check f_k4k4 first (simpler)
f_44_simplified = simplify(f_44)
print(f"\n  f_k4k4 = {f_44_simplified}")

# Substitute k4 = w*k2^2 for simplification
w = symbols('w', real=True)
f_44_sub = simplify(f_44.subs(k4, w*k2**2))
print(f"  f_k4k4 at k4=w*k2^2 = {simplify(f_44_sub)}")

# Check: for concavity, f_k4k4 should be <= 0
# f(k2,k4) = k4^2 / (24*k2*(4*k2^2-k4))
# f_k4 = [2*k4*(4*k2^2-k4) + k4^2] / (24*k2*(4*k2^2-k4)^2)
#       = k4*(8*k2^2-k4) / (24*k2*(4*k2^2-k4)^2)

f_k4 = simplify(diff(f, k4))
print(f"\n  f_k4 = {f_k4}")

# f_k4k4: differentiate again
# This is getting complicated. Let's just evaluate numerically.

print("\n\nNumerical check of Hessian eigenvalues:")
import numpy as np

# Grid of (k2, w) values where k4 = w*k2^2
for k2_val in [0.5, 1.0, 2.0, 5.0]:
    for w_val in [-3.0, -1.0, 0.0, 1.0, 2.0, 3.0, 3.9]:
        k4_val = w_val * k2_val**2
        if 4*k2_val**2 - k4_val <= 0:
            continue

        # Evaluate Hessian numerically
        h11 = float(f_22.subs({k2: k2_val, k4: k4_val}))
        h12 = float(f_24.subs({k2: k2_val, k4: k4_val}))
        h22 = float(f_44.subs({k2: k2_val, k4: k4_val}))

        H_num = np.array([[h11, h12], [h12, h22]])
        eigs = np.linalg.eigvalsh(H_num)

        if eigs[1] > 1e-10:  # Positive eigenvalue => NOT concave
            print(f"  k2={k2_val}, w={w_val}: eigs = {eigs[0]:.6f}, {eigs[1]:.6f}  NOT CONCAVE!")
        else:
            pass  # concave at this point

# Check for a single point in detail
print("\nDetailed check at k2=1, k4=0:")
k2_v, k4_v = 1.0, 0.0
h11 = float(f_22.subs({k2: k2_v, k4: k4_v}))
h12 = float(f_24.subs({k2: k2_v, k4: k4_v}))
h22 = float(f_44.subs({k2: k2_v, k4: k4_v}))
print(f"  H = [[{h11:.6f}, {h12:.6f}], [{h12:.6f}, {h22:.6f}]]")
eigs = np.linalg.eigvalsh(np.array([[h11, h12], [h12, h22]]))
print(f"  eigenvalues: {eigs}")

print("\nDetailed check at k2=1, k4=2:")
k2_v, k4_v = 1.0, 2.0
h11 = float(f_22.subs({k2: k2_v, k4: k4_v}))
h12 = float(f_24.subs({k2: k2_v, k4: k4_v}))
h22 = float(f_44.subs({k2: k2_v, k4: k4_v}))
print(f"  H = [[{h11:.6f}, {h12:.6f}], [{h12:.6f}, {h22:.6f}]]")
eigs = np.linalg.eigvalsh(np.array([[h11, h12], [h12, h22]]))
print(f"  eigenvalues: {eigs}")

print("\nDetailed check at k2=1, k4=-2:")
k2_v, k4_v = 1.0, -2.0
h11 = float(f_22.subs({k2: k2_v, k4: k4_v}))
h12 = float(f_24.subs({k2: k2_v, k4: k4_v}))
h22 = float(f_44.subs({k2: k2_v, k4: k4_v}))
print(f"  H = [[{h11:.6f}, {h12:.6f}], [{h12:.6f}, {h22:.6f}]]")
eigs = np.linalg.eigvalsh(np.array([[h11, h12], [h12, h22]]))
print(f"  eigenvalues: {eigs}")

# Check across many random points
print("\n\nSystematic concavity check (random sampling):")
np.random.seed(42)
not_concave_count = 0
total_checks = 0
for _ in range(10000):
    k2_v = np.random.uniform(0.1, 10)
    w_v = np.random.uniform(-3.9, 3.9)
    k4_v = w_v * k2_v**2
    if 4*k2_v**2 - k4_v <= 0:
        continue
    total_checks += 1

    h11 = float(f_22.subs({k2: k2_v, k4: k4_v}))
    h12 = float(f_24.subs({k2: k2_v, k4: k4_v}))
    h22 = float(f_44.subs({k2: k2_v, k4: k4_v}))

    H_num = np.array([[h11, h12], [h12, h22]])
    eigs = np.linalg.eigvalsh(H_num)

    if eigs[1] > 1e-10:
        not_concave_count += 1
        if not_concave_count <= 5:
            print(f"  NOT CONCAVE at k2={k2_v:.3f}, k4={k4_v:.3f}: eigs={eigs}")

if not_concave_count == 0:
    print(f"  CONCAVE at all {total_checks} test points!")
else:
    print(f"\n  NOT concave at {not_concave_count}/{total_checks} points")
    print("  Concavity approach FAILS for the k3=0 case.")
    print("  Need a different proof strategy.")
