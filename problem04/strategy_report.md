# Strategy Report: Hybrid Fisher Superadditivity Proof

## Priority-Ordered Work Items

### TIER 1: Highest Impact (attack first)

**1. Node 1.4.2 — 1/Φ Concavity Along Heat Flow** (Path B key step)
- **Why first**: Most concrete and computationally tractable. Reduces to showing Φ·Φ'' ≥ 2(Φ')² using root dynamics. All ingredients are available (Dyson dynamics, Φ_n=2·Sm2).
- **Agent type**: Prover with strong calculus/ODE skills
- **Approach**: Compute Φ', Φ'' explicitly via dr_i/dt = H(r_i), then verify the Cauchy-Schwarz structure
- **Fallback**: If exact proof fails, try perturbative proof near Gaussian (large t) + monotonicity argument

**2. Node 1.3.1 — Finite Subordination Existence** (Path A infrastructure)
- **Why**: Unlocks entire Path A. Without this, nodes 1.3.2-1.3.5 are conditional.
- **Agent type**: Prover with complex analysis background
- **Approach**: For G_r(z) = Σ 1/(z-λ_i), solve G_p(w) = G_r(z) for w = ω_p(z). Show solution is a degree-(n-1) rational Herglotz function via Rouché's theorem + argument principle.
- **Key reference**: Biane (1998) finite free subordination, Belinschi-Bercovici (2007)

### TIER 2: High Impact (attack after Tier 1 progress)

**3. Node 1.3.4 — Herglotz Coupling Lemma** (Path A key step)
- **Depends on**: 1.3.1 (subordination existence)
- **Why**: This is THE hard step of Path A. Equivalent to the full conjecture given subordination.
- **Approach**: Express J_p, J_q as quadratic forms in Herglotz residues. Use coupling ω_p+ω_q = z+F_r(z) to constrain. Apply matrix AM-GM or Schur complement.
- **Risk**: May be as hard as the original conjecture in a different language

**4. Node 1.4.3 — Stam from Concavity + Splitting** (Path B synthesis)
- **Depends on**: 1.4.2 (concavity proved)
- **Why**: If concavity is proved, this step should be routine (classical template).
- **Approach**: Follow Costa/Villani derivation. May need careful analysis of t→0 limit.

### TIER 3: Independent / Exploratory

**5. Node 1.5.1 — Finite Free EPI** (Path C)
- **Independent** of Paths A and B
- **Approach 1**: Random matrix representation + Jensen's inequality
- **Approach 2**: HCIZ formula + Prékopa-Leindler
- **Approach 3**: Direct computation for small n, then pattern

**6. Node 1.3.4.2 — Coupling Constraint** (Path A detail)
- **Key sub-question**: Is ω_p(z)+ω_q(z) = z + F_r(z) correct at finite n?
- **Approach**: Verify numerically for n=2,3,4, then prove algebraically from MSS structure

## Recommended Agent Deployment

### Wave 1 (Parallel, independent)
| Agent | Node | Task |
|-------|------|------|
| PROVER-1 | 1.4.2 | Compute Φ', Φ'' via root dynamics, prove concavity |
| PROVER-2 | 1.3.1 | Construct finite subordination, prove Herglotz property |
| PROVER-3 | 1.5.1-1.5.2 | Attack finite EPI via random matrix / HCIZ |

### Wave 2 (After Wave 1 results)
| Agent | Node | Task |
|-------|------|------|
| VERIFIER-1 | Wave 1 outputs | Adversarial verification of all claims |
| PROVER-4 | 1.3.4 or 1.4.3 | Attack whichever key step is unblocked |

## Critical Warnings (from ex6/ex7 experience)

1. **DO NOT** assume ⟨h,α⟩ ≥ 0 — it is FALSE (ex6 Session 129)
2. **DO NOT** assume monotone gap along heat flow — KILLED (ex7 Session 132, 44% violations)
3. **DO NOT** try polynomial manipulation for n≥4 general — exhausted by 4 prover agents in ex7
4. **DO NOT** use the wrong boxplus formula — use ONLY c_k = Σ_{i+j=k} [(n-i)!(n-j)!/(n!(n-k)!)] a_i b_j
5. **DO NOT** assume partition of unity ω₁'+ω₂'=1 — FALSE (correct: each ω'=1 independently)

## Success Metrics

- **Full proof**: Any one of the three KEY HARD STEPS (1.3.4, 1.4.2, 1.5.1) proved → conjecture proved
- **Partial success**: Subordination existence (1.3.1) proved → major infrastructure for future work
- **Minimum viable**: 1/Φ concavity computation completed (even if Cauchy-Schwarz closure fails) → valuable mathematical content
