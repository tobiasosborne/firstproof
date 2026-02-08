# Proof Strategy Reports for Hard Lemma 2

Three independent strategies for proving the main inequality in the Fisher superadditivity conjecture.

---

## Strategy 1: Contour Integral & Electrostatic Energy (Agent 1)

**Core idea:** Express Phi_n(p) as a contour integral via the residue formula
`Phi_n(p) = (1/2pi i) oint p''(z)^2 / (4 p'(z) p(z)) dz`, then use subordination
functions as changes of variables to relate the three integrals for p, q, r.

**Key insight:** The subordination functions provide exact changes of variables between
integration contours. The condition omega_i'(nu_k) = 1 makes the change of variables
isometric at the poles.

**Novel contribution:** The contour integral representation of Phi_n appears to be new
in this context. Also discovered the explicit formula for alpha_k via implicit
differentiation of F(z,w) = r'(z)p(w) - p'(w)r(z) = 0.

**Status:** Promising but incomplete. The Schwarz-Pick argument doesn't directly apply
(omega maps C^+ to C^+, not to a smaller domain), and the integrand may not have
definite sign on complex contours.

---

## Strategy 2: Finite Cramér-Rao / Random Matrix (Agent 2)

**Core idea:** Identify 1/Phi_n(p) as the Cramér-Rao bound (minimum variance) for
estimating the location parameter of the squared-Vandermonde density Delta(lambda)^2
on root configurations. The MSS convolution corresponds to an HCIZ-mediated convolution
of these densities, and the Stam inequality (information-theoretic convolution inequality)
gives superadditivity directly.

**Key insight:** The HCIZ integral satisfies translation covariance:
I_n(A + tI, B) = exp(t tr(B)) I_n(A, B), which is exactly the structure needed for
the Stam inequality proof.

**Novel contribution:** Reframes the entire problem in information-theoretic language,
bypassing subordination entirely.

**Status:** Elegant but has significant gaps. The squared Vandermonde is not normalizable
without a confining potential (needs conditioning on tr = const). The connection between
the MSS expected characteristic polynomial and the mode of the HCIZ-convolved density
needs Jensen's inequality for 1/Phi_n, requiring convexity that is unverified.

---

## Strategy 3: Lorentzian Polynomials / Subordination Partition of Unity (Agent 3)

**Core idea:** Use the partition of unity omega_1'(gamma_k) + omega_2'(gamma_k) = 1
at the roots of r to decompose Phi_n(r) = Phi_n^(1) + Phi_n^(2), then prove a
"weighted subordination contraction" Key Lemma and apply Cauchy-Schwarz.

**KEY RESULT: CONDITIONAL PROOF FOUND!**

The inequality 1/Phi_n(r) >= 1/Phi_n(p) + 1/Phi_n(q) follows from:

**Key Lemma:** sum_k H_p(omega_1(gamma_k))^2 omega_1'(gamma_k) <= Phi_n(p)

**Proof of main theorem from Key Lemma:**
1. From chain rule: H_r(gamma_k) = H_p(omega_1(gamma_k)) omega_1'(gamma_k)
   So: sum_k H_r(gamma_k)^2 / omega_1'(gamma_k) = sum_k H_p(omega_1(gamma_k))^2 omega_1'(gamma_k) <= Phi_n(p)
2. By Cauchy-Schwarz (Titu/Engel form):
   sum_k H_r^2/omega_1' >= (sum H_r^2)^2 / (sum H_r^2 omega_1')
   So: 1/Phi_n(p) <= 1/(sum H_r^2/omega_1') <= (sum H_r^2 omega_1')/Phi_n(r)^2
3. Similarly: 1/Phi_n(q) <= (sum H_r^2 omega_2')/Phi_n(r)^2
4. Adding: 1/Phi_n(p) + 1/Phi_n(q) <= (Phi_n^(1) + Phi_n^(2))/Phi_n(r)^2 = 1/Phi_n(r)

**CRITICAL DEPENDENCIES:**
- omega_1'(gamma_k) + omega_2'(gamma_k) = 1 (needs verification at finite n)
- The Key Lemma (weighted subordination contraction)
- H_r(gamma_k) = H_p(omega_1(gamma_k)) * omega_1'(gamma_k) (NOT omega_1'=1 as in the
  chain rule from node 1.6! This uses the GENERAL chain rule before specializing.)

**WARNING:** The chain rule in node 1.6 established omega_1'(nu_k) = 1, which would make
omega_1' + omega_2' = 2, NOT 1. This contradicts the partition of unity. The resolution
may be that the partition of unity holds for a DIFFERENT normalization of the subordination
functions, or that this strategy requires rethinking the relationship between the chain
rule and the partition of unity.

---

## Comparison and Recommended Next Steps

| Strategy | Novelty | Feasibility | Key Gap |
|----------|---------|-------------|---------|
| 1 (Contour) | High | Medium | Integrand sign control |
| 2 (Cramér-Rao) | Very High | Low-Medium | Normalization + Jensen |
| 3 (Partition) | Medium-High | **HIGH** | Partition of unity + Key Lemma |

**Strategy 3 is the most promising.** It gives an explicit conditional proof modulo
two lemmas. However, there is a critical tension: the established chain rule says
omega_1'(nu_k) = 1, while the partition of unity needs omega_1' + omega_2' = 1.
These cannot both hold unless omega_1' = omega_2' = 1/2, which seems too special.

**Immediate next step:** Verify numerically at n=2,3 whether:
(a) omega_1'(nu_k) + omega_2'(nu_k) = 1 or = 2 at the roots
(b) The Key Lemma holds
(c) Whether there is a different decomposition that reconciles both facts
