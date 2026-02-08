"""
PROVER-11 Part 5: Proper domain verification for k3=0 case.

Key finding from Part 4: the "violations" were caused by sampling outside
the domain (u > 4s^2). Let me verify with proper domain constraints.

Domain for R4(K2, 0, K4):
  K2 > 0
  K4 < 4*K2^2  (factor 1 positive)
  K4 > -4*K2^2  (factor 2 positive, with K3=0)

So the proper check: u ∈ (-4s^2, 4s^2), v ∈ (-4t^2, 4t^2).
AND for the sum: u+v ∈ (-4(s+t)^2, 4(s+t)^2).

Also: R4 has a sign issue. When K4 < 0 and |K4| > 4*K2^2:
  4*K2^2 - K4 > 0  (since K4 < 0)
  But R4 = -K4^2/(24*K2*(4K2^2-K4))
  K4^2 > 0, K2 > 0, 4K2^2-K4 > 0
  So R4 < 0 always. ✓

Let me also check: is R4 actually positive sometimes?
R4 = -u^2/(24*s*(4s^2-u)).
If u > 0: 4s^2-u > 0, so R4 = -u^2/(positive) < 0.
If u < 0: 4s^2-u > 4s^2 > 0, so R4 = -u^2/(positive) < 0.
If u = 0: R4 = 0.
So R4 <= 0 always when k3=0. ✓
"""
import numpy as np

np.random.seed(42)
violations = 0
n_trials = 500000
min_gap = 1e10
valid_trials = 0

for _ in range(n_trials):
    s_val = np.random.exponential(1) + 0.01
    t_val = np.random.exponential(1) + 0.01

    # Strict domain: u in (-4s^2 + eps, 4s^2 - eps)
    u_max = 4*s_val**2 - 0.0001
    u_min = -4*s_val**2 + 0.0001
    v_max = 4*t_val**2 - 0.0001
    v_min = -4*t_val**2 + 0.0001

    u_val = np.random.uniform(u_min, u_max)
    v_val = np.random.uniform(v_min, v_max)

    sr = s_val + t_val
    ur = u_val + v_val

    # Check sum domain
    if ur >= 4*sr**2 - 0.0001 or ur <= -4*sr**2 + 0.0001:
        continue

    # Verify each point is in domain
    assert u_val < 4*s_val**2 and u_val > -4*s_val**2
    assert v_val < 4*t_val**2 and v_val > -4*t_val**2

    valid_trials += 1

    def R4_k30(K2, K4):
        return -K4**2 / (24 * K2 * (4*K2**2 - K4))

    r_sum = R4_k30(sr, ur)
    r_p = R4_k30(s_val, u_val)
    r_q = R4_k30(t_val, v_val)
    gap = r_sum - r_p - r_q
    min_gap = min(min_gap, gap)

    if gap < -1e-10:
        violations += 1
        if violations <= 3:
            print(f"VIOLATION: s={s_val:.6f}, t={t_val:.6f}, u={u_val:.6f}, v={v_val:.6f}")
            print(f"  u/4s^2={u_val/(4*s_val**2):.4f}, v/4t^2={v_val/(4*t_val**2):.4f}")
            print(f"  R_p={r_p:.10f}, R_q={r_q:.10f}, R_sum={r_sum:.10f}, gap={gap:.6e}")
            # Check Cauchy-Schwarz
            A = 4*s_val**2 - u_val
            B = 4*t_val**2 - v_val
            C = 4*sr**2 - ur
            print(f"  sA+tB={s_val*A+t_val*B:.6f}, (s+t)C={sr*C:.6f}")
            lhs = ur**2 / (sr*C)
            rhs = u_val**2/(s_val*A) + v_val**2/(t_val*B)
            print(f"  (u+v)^2/((s+t)C)={lhs:.10f}, u^2/(sA)+v^2/(tB)={rhs:.10f}")

print(f"\n{violations} violations in {valid_trials} valid trials (out of {n_trials} total)")
print(f"Minimum gap: {min_gap:.6e}")

if violations == 0:
    print("\n*** k3=0 CASE: VERIFIED (0 violations in 500K trials) ***")
    print("\nThe Cauchy-Schwarz proof is correct:")
    print("  Step 1: s*A + t*B <= (s+t)*C  (proved: difference = 12st(s+t)-tu-sv >= 8st(s+t) > 0)")
    print("  Step 2: CS inequality => (u+v)^2/(s*A+t*B) <= u^2/(s*A) + v^2/(t*B)")
    print("  Combined: (u+v)^2/((s+t)*C) <= u^2/(s*A) + v^2/(t*B)")
