# HANDOFF: Hybrid Fisher Superadditivity Proof (examples8)

**Created:** 2026-02-08
**Origin:** Synthesis of examples6 (subordination) and examples7 (cumulant/heat flow)

---

## 1. The Conjecture

For all monic real-rooted polynomials p, q of degree n >= 2 with simple roots:

```
1/Phi_n(p ⊞_n q) >= 1/Phi_n(p) + 1/Phi_n(q)
```

This is the finite-n analog of Voiculescu's free Stam inequality (1998).

---

## 2. Why This Proof Tree Exists

Two independent projects (examples6, examples7) attacked this conjecture:
- **ex6**: Subordination approach — deep infrastructure (chain rule, corrections), but Hard Lemma 2 open
- **ex7**: Cumulant decomposition + heat flow — rich structural identities, but polynomial manipulation exhausted

Both proved n=2,3 and n=4-symmetric independently. Both hit the same wall at general n>=4.

**This tree combines their strengths:**
- ex6's subordination provides the "why" (L² contraction framework)
- ex7's de Bruijn/heat flow provides the "how" (entropy + concavity tools)
- The classical Voiculescu/Shlyakhtenko-Tao/Costa template provides the roadmap

---

## 3. Three Proof Paths (any one suffices)

### Path A: Subordination + L² Contraction (PRIMARY — nodes 1.3.*)
- Construct finite subordination ω_p, ω_q (Hard Lemma 1 from ex6)
- Use chain rule H_r(λ_i) = H_p(μ_σ(i)) - α_i (proved in ex6)
- L² Pythagoras: Φ_p = Φ_r + J_p where J_p = 2⟨h,α⟩ + ‖α‖²
- **KEY STEP (node 1.3.4)**: Prove J_p·J_q ≥ ‖h‖⁴ from Herglotz coupling
- This is the finite analog of Shlyakhtenko-Tao's analytic proof

### Path B: De Bruijn + 1/Φ Concavity (SECONDARY — nodes 1.4.*)
- De Bruijn identity dS/dt = Φ_n (validated in ex7)
- Gaussian splitting (proved in ex7)
- **KEY STEP (node 1.4.2)**: Prove 1/Φ_n(p⊞G_t) concave in t
- Derive Stam from concavity (following Costa/Villani)
- This is the finite analog of Costa's entropy power concavity

### Path C: Entropy Power Inequality (TERTIARY — nodes 1.5.*)
- Finite EPI: N(p⊞q) ≥ N(p)+N(q), numerically validated (13K+ trials)
- **KEY STEP (node 1.5.1-1.5.2)**: Prove finite EPI
- Derive Stam via differentiation (classical EPI→Stam reduction)

---

## 4. What Is Already Proved (import from ex6/ex7)

| Result | Source | Node |
|--------|--------|------|
| n=2 equality | Both | 1.2 |
| n=3 full proof | Both (independent) | 1.2 |
| n=4 symmetric | Both (independent) | 1.2 |
| Φ_n = 2·Sm2 | ex7 PROVER-14 | 1.1 |
| S2 additivity | ex7 PROVER-14 | 1.1 |
| De Bruijn identity | ex7 PROVER-13 | 1.4.1 |
| Gaussian splitting | ex7 PROVER-13 | 1.4.4 |
| ω'(λ_i) = 1 | ex6 chain rule | 1.3.2 |
| C_n = 4/(n²(n-1)) | ex7 PROVER-10b | 1.1 |

## 5. What Is Open (the hard steps)

| Hard Step | Node | Description | Numerics |
|-----------|------|-------------|----------|
| Subordination existence | 1.3.1 | Finite Herglotz ω_p, ω_q | Verified n=2-5 |
| Herglotz coupling | 1.3.4 | J_p·J_q ≥ ‖h‖⁴ | Equivalent to conjecture |
| 1/Φ concavity | 1.4.2 | d²/dt²(1/Φ) ≤ 0 | 0/590 violations |
| Finite EPI | 1.5.1 | N(p⊞q) ≥ N(p)+N(q) | 0/13,770 violations |

## 6. Exhausted Approaches (DO NOT RETRY)

From ex6: ⟨h,α⟩ ≥ 0 is FALSE, partition of unity is FALSE, shape factor bound KILLED
From ex7: Joint concavity of -R_4 FAILS, SOS FAILS, perspective function BLOCKED, monotone gap KILLED, all polynomial manipulation exhausted for n≥4 general

## 7. Key External References

- MSS: Marcus-Spielman-Srivastava, Interlacing Families II (2015)
- Cumulant additivity: Arizmendi-Perales (2018)
- Free Stam: Voiculescu, Invent. Math. 132 (1998)
- Analytic proof: Shlyakhtenko-Tao, arXiv:2009.01882 (2020)
- Entropy monotonicity: Gribinski, arXiv:1907.12826 (2019)

## 8. Recommended Agent Strategy

1. **First wave**: One prover on 1.3.1 (subordination existence), one on 1.4.2 (1/Φ concavity computation). These are independent.
2. **Second wave**: If subordination exists, attack 1.3.4 (coupling lemma). If concavity computed, attack 1.4.3 (Stam derivation).
3. **Verification**: Each claimed result gets an adversarial verifier before acceptance.
4. **Fallback**: If Path A and B both stall, Path C (EPI via random matrix/Brascamp-Lieb) is independent.

## 9. Source Documents

- **Conjecture source**: `examples6/fisher_subordination_proof.md`
- **ex6 full status**: `examples6/HANDOFF.md`
- **ex7 full status**: `examples7/HANDOFF.md`
- **Numerical infrastructure**: `examples7/*.py` (28K+ lines of verification scripts)
