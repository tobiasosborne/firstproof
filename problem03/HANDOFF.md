# Problem 03 Handoff — Markov Chain with ASEP Polynomial Stationary Distribution

## What This Problem Asks

Let λ = (λ₁ > ⋯ > λₙ ≥ 0) be a restricted partition with distinct parts (unique part 0, no part 1). Let S_n(λ) be the set of compositions of λ. Does there exist a **nontrivial** Markov chain on S_n(λ) whose stationary distribution is:

> π(μ) = F\*_μ(x₁,…,xₙ; q=1, t) / P\*_λ(x₁,…,xₙ; q=1, t)

where "nontrivial" means the transition probabilities are NOT described using the F\*_μ polynomials?

## Answer

**YES.** The inhomogeneous multispecies t-PushTASEP provides such a chain.

## Current State (Session 1)

- **Proof tree:** 9 nodes (root + 8 children), all pending
- **0 VALIDATED**, 0 refuted, 0 archived
- **7 external references** registered (AMW24, CMW22, BDW25, FM07, KS96, ABW23, W22)
- **6 definitions** registered (partition, restricted, S_n(λ), F\*_μ, P\*_λ, Markov chain)
- **Run `af status`** from this directory to see the full tree.

### Proof Strategy (Hybrid from 3 independent analyses)

Three independent proof strategies were generated:
- **Strategy A** (algebraic combinatorics): t-PushTASEP + multiline queue weights
- **Strategy B** (stochastic processes): t-PushTASEP + Ferrari-Martin lumping
- **Strategy C** (Hecke algebra / Yang-Baxter): R-matrix swap chain + transfer matrix eigenvector

All three agreed on YES and converged on the t-PushTASEP as the chain. The registered proof tree is a hybrid taking the strongest elements of each.

### Proof Tree Structure

| Node | Description | Status | Risk |
|------|-------------|--------|------|
| 1 | Root conjecture | pending | — |
| 1.1 | State space setup: S_n(λ) ↔ particle configs on ring | pending | Low |
| 1.2 | Polynomial identification: F\*_μ via Hecke operators, f_μ as t-PushTASEP weights | pending | Low |
| 1.3 | Positivity and normalization: π(μ) ≥ 0, sums to 1 | pending | Medium |
| 1.4 | Chain construction: inhomogeneous multispecies t-PushTASEP | pending | Low |
| 1.5 | Stationarity: multiline process (Ferrari-Martin/AMW24) + YBE transfer matrix | pending | Medium |
| 1.6 | **Interpolation = homogeneous ratio at q=1** | pending | **HIGH** |
| 1.7 | Nontriviality: local rates vs global polynomials | pending | Low |
| 1.8 | Conclusion (QED) | pending | Low |

### Critical Vulnerabilities

**Node 1.6 is the hardest step.** All three strategies flagged it. The claim is:

> F\*_μ(x; q=1, t) / P\*_λ(x; q=1, t) = f_μ(x; q=1, t) / P_λ(x; q=1, t)

The current argument: both ratios sum to 1, and the t-PushTASEP has a unique stationary distribution, so they must agree. But this requires proving that F\*_μ/P\*_λ **also satisfies the stationarity equations** — not just that it sums to 1. Two backup approaches:

1. **Hecke algebra route (from Strategy C):** The Demazure-Lusztig relations T_i F\*_μ = F\*_{s_i μ} hold for the full interpolation polynomials (not just top degree). If the transfer matrix eigenvector argument works for F\*_μ, stationarity follows directly and the ratio identity becomes a corollary.

2. **Direct verification route:** At q=1, explicitly check that the interpolation corrections g_μ = F\*_μ − f_μ satisfy g_μ/f_μ = G/P_λ (same proportional correction for each μ), forcing the ratios to agree.

**Node 1.3 positivity** is also nontrivial: the signed multiline queue formula for interpolation ASEP polynomials does not have manifestly positive terms. At q=1, positivity follows from the probabilistic interpretation (stationary probabilities are nonneg), but this is circular if we haven't yet established that F\*_μ/P\*_λ = f_μ/P_λ.

### Next Steps (Session 2)

1. **Priority 1:** Sub-refine node 1.6 — break the interpolation=homogeneous ratio claim into verifiable sub-steps. Consider the Hecke algebra approach.
2. **Priority 2:** Verify nodes 1.1, 1.4, 1.7 (low-risk definitional nodes) to build up validated foundation.
3. **Priority 3:** Refine node 1.5 with explicit details of the multiline process construction.
4. **Priority 4:** Address positivity (node 1.3) — may depend on resolution of 1.6.

### Key References

| Tag | Paper | Role in proof |
|-----|-------|---------------|
| AMW24 | Ayyer-Martin-Williams, arXiv:2310.09740 | t-PushTASEP stationary distribution = f_μ at q=1 |
| CMW22 | Corteel-Mandelshtam-Williams, arXiv:1811.01024 | ASEP polynomials via multiline queues |
| BDW25 | Ben Dali-Williams, arXiv:2510.02587 | Interpolation ASEP polynomials, Hecke relations (Prop 2.10) |
| FM07 | Ferrari-Martin, Ann. Probab. 35(3) | Multiline construction for multispecies TASEP |
| KS96 | Knop (1997), Sahi (1996) | Interpolation Macdonald polynomials |
| ABW23 | Aggarwal-Borodin-Wheeler, arXiv:2309.11865 | Yang-Baxter → stationary measures |
| W22 | Williams survey, arXiv:2202.00214 | Hopping particles and positivity |

### Lessons Learned

1. **All three strategies converged** on the same answer and same chain — strong signal that YES is correct.
2. **The interpolation vs. homogeneous distinction** is the mathematical crux. The problem deliberately asks about F\*_μ (interpolation), not f_μ (homogeneous). This is likely the "trap" that makes the problem research-level.
3. **The restricted partition condition** (no part 1, unique part 0) is needed for species labels to be well-separated and for the multiline queue / positivity arguments to work cleanly.
4. **q=1 is the magic specialization** — it collapses the double affine Hecke algebra to the affine Hecke algebra, making the stochastic R-matrix interpretation possible.
