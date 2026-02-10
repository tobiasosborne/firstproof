# Problem 03 Handoff — Markov Chain with ASEP Polynomial Stationary Distribution

## What This Problem Asks

Let λ = (λ₁ > ⋯ > λₙ ≥ 0) be a restricted partition with distinct parts (unique part 0, no part 1). Let S_n(λ) be the set of compositions of λ. Does there exist a **nontrivial** Markov chain on S_n(λ) whose stationary distribution is:

> π(μ) = F\*_μ(x₁,…,xₙ; q=1, t) / P\*_λ(x₁,…,xₙ; q=1, t)

where "nontrivial" means the transition probabilities are NOT described using the F\*_μ polynomials?

## Answer

**YES.** The inhomogeneous multispecies t-PushTASEP provides such a chain.

## Current State (Session 2 — Verification Wave 1 Complete)

- **Proof tree:** 9 nodes (root + 8 children), all pending
- **0 VALIDATED**, 0 refuted, 0 archived
- **97 open challenges** across all 8 leaf nodes
- **ALL 8 leaf nodes CHALLENGED** — none accepted
- Run `af status` from this directory to see the full tree

### Verification Wave 1 Results (Session 2)

All 8 leaf nodes (1.1–1.8) were verified breadth-first by independent verifier subagents. Every node was challenged. The proof as currently stated has **fundamental structural problems** beyond individual node defects.

| Node | Description | Challenges | Severity | Key Issue |
|------|-------------|-----------|----------|-----------|
| 1.1 | State space setup | 2 major | CHALLENGED | "Well-separated" undefined; missing site-occupancy spec |
| 1.2 | Polynomial identification | 1 critical, 4 major | CHALLENGED | **Wrong arXiv for AMW24**; F\*/f\* conflation; notation chaos |
| 1.3 | Positivity & normalization | 2 critical, 2 major | CHALLENGED | **Circular dependency** on 1.6; F\*_μ nonnegativity is open |
| 1.4 | Chain construction | 2 critical, 4 major | CHALLENGED | **Cascade mechanism omitted**; ambiguous dynamics |
| 1.5 | Stationarity | 3 critical, 3 major | CHALLENGED | **Proof method is fictional** — no FM multiline process for t-PushTASEP at general t |
| 1.6 | Interpolation = homogeneous ratio | 4 critical | CHALLENGED | **Fatal logical fallacy** ("both sum to 1 ∴ equal"); open research problem |
| 1.7 | Nontriviality | 3 critical, 1 major | CHALLENGED | **False locality claim**; no rigorous exclusion argument |
| 1.8 | Conclusion | 1 critical, 1 major | CHALLENGED | Depends on all siblings, none validated |

### Systemic Problems Found

1. **Wrong reference throughout**: AMW24 is cited as arXiv:2310.09740 (Ayyer-Martin, t=0 only). The correct paper for t-PushTASEP at general t is likely arXiv:2403.10485 (Ayyer-Martin-Williams). This must be verified and corrected.

2. **Node 1.5 proof method doesn't exist**: The Ferrari-Martin multiline Markov process construction does NOT generalize to t-PushTASEP at general t. AMW24 (2403.10485) proves stationarity via algebraic Macdonald polynomial methods, not a multiline Markov process. Node 1.5 must be completely rewritten.

3. **Node 1.6 contains a textbook fallacy**: "Both ratios sum to 1, therefore they're equal" is logically invalid (1 constraint on n! unknowns). The handoff from Session 1 already flagged this. BDW25 Remark 1.17 defers the interpolation probabilistic interpretation to a forthcoming paper [BDW] that has NOT appeared.

4. **Circular dependency**: Node 1.3 (positivity of F\*_μ) depends on Node 1.6 (ratio identity), which in turn would need 1.3 for its uniqueness argument. The proof tree must be restructured.

5. **Notation inconsistency across papers**: AMW24 uses F_η (uppercase, no star), BDW25 distinguishes f\*_μ (algebraic) from F\*_μ (combinatorial). The proof conflates these freely without establishing mappings.

6. **No node declares dependencies**: Every node is listed as independent, but logically they form a chain with complex interdependencies.

### Proof Tree Structure (Updated Risk Assessment)

| Node | Description | Status | Risk |
|------|-------------|--------|------|
| 1 | Root conjecture | pending | — |
| 1.1 | State space setup | CHALLENGED (2 major) | **Medium** — fixable |
| 1.2 | Polynomial identification | CHALLENGED (1 crit, 4 major) | **High** — wrong reference, notation mess |
| 1.3 | Positivity and normalization | CHALLENGED (2 crit, 2 major) | **Critical** — circular, open problem |
| 1.4 | Chain construction | CHALLENGED (2 crit, 4 major) | **High** — cascade omitted, CT/DT mismatch |
| 1.5 | Stationarity | CHALLENGED (3 crit, 3 major) | **Critical** — fictional proof method |
| 1.6 | Interpolation = homogeneous ratio | CHALLENGED (4 crit) | **Critical** — logical fallacy, open problem |
| 1.7 | Nontriviality | CHALLENGED (3 crit, 1 major) | **High** — false premise, no proof |
| 1.8 | Conclusion | CHALLENGED (1 crit, 1 major) | Blocked by all above |

### Next Steps (Session 3)

**The proof needs a substantial rewrite, not just patching.** Recommended approach:

1. **Fix the reference**: Verify arXiv:2403.10485 is the correct AMW24 paper. Update the external reference.

2. **Rewrite Node 1.4** (chain construction): Include the full cascade mechanism, specify vacancy behavior, address CT vs DT, establish irreducibility. This is fixable.

3. **Rewrite Node 1.5** (stationarity): Cite AMW24 Theorem 1.1 directly. Describe the actual algebraic proof method. Do NOT invoke Ferrari-Martin for general t.

4. **Restructure the proof tree to break the circularity**:
   - Node 1.5 proves: stationary dist of t-PushTASEP ∝ f_μ(x;q=1,t) (homogeneous) [cite AMW24]
   - Node 1.6 proves: F\*_μ/P\*_λ = f_μ/P_λ at q=1 [this is the hard step]
   - Node 1.3 becomes a COROLLARY of 1.5 + 1.6 (positivity follows from ratio identity + f_μ ≥ 0)

5. **Attack Node 1.6 seriously**: This is the crux. Three possible approaches:
   - (a) **Hecke algebra**: Show T_i-relations for F\*_μ at q=1 imply F\*_μ/P\*_λ satisfies balance equations → uniqueness gives ratio identity
   - (b) **BDW25 Theorem 7.1**: Check if the factorization of interpolation Macdonald polynomials at q=1 implies a universal factor F\*_μ = f_μ · C(x,t) with C independent of μ
   - (c) **Direct n=2 verification**: Start with λ=(2,0), use BDW25 Example 1.16 to check if the ratio identity holds. If it fails, the answer might be NO (or the proof strategy is wrong).

6. **Clarify nontriviality (Node 1.7)**: Define "described using" formally, then provide a rigorous argument (degree counting, or simply argue the t-PushTASEP is not the trivial chain P(μ→ν) = π(ν)).

### Key References (Corrected)

| Tag | Paper | Role in proof | Note |
|-----|-------|---------------|------|
| AMW24 | Ayyer-Martin-Williams, **arXiv:2403.10485** (?) | t-PushTASEP stationary distribution | **VERIFY THIS arXiv NUMBER** — 2310.09740 is t=0 only |
| CMW22 | Corteel-Mandelshtam-Williams, arXiv:1811.01024 | ASEP polynomials via multiline queues | |
| BDW25 | Ben Dali-Williams, arXiv:2510.02587 | Interpolation ASEP polys, Hecke relations | Key: Prop 2.10, Prop 2.15, Thm 1.15, Thm 7.1, Remark 1.17 |
| FM07 | Ferrari-Martin, Ann. Probab. 35(3) | Multiline construction — **ordinary TASEP only** | Does NOT apply to t-PushTASEP |
| KS96 | Knop (1997), Sahi (1996) | Interpolation Macdonald polynomials | |
| ABW23 | Aggarwal-Borodin-Wheeler, arXiv:2309.11865 | Yang-Baxter → stationary measures | |
| W22 | Williams survey, arXiv:2202.00214 | Hopping particles and positivity | |

### Lessons Learned (Session 2)

1. **Every "low-risk" node had significant defects.** The initial risk assessment was wildly optimistic. Only 1.1 has purely fixable issues; all others have critical gaps.
2. **The proof method in Node 1.5 was fabricated** — no Ferrari-Martin multiline process exists for general t. The AI prover hallucinated a proof technique.
3. **Node 1.6's argument is a textbook logical fallacy.** This was predicted by Session 1's handoff but was never addressed.
4. **The interpolation ASEP polynomial positivity at q=1 is an open research problem** (BDW25 Remark 1.17 defers it to a forthcoming paper).
5. **Notation hygiene matters**: Three papers use incompatible conventions (F vs f, * vs no-*, η vs μ). A clean notation table is essential before any rewrite.
6. **The overall proof strategy (YES via t-PushTASEP) is likely correct** — the answer is almost certainly YES — but the current proof tree is too defective to salvage by patching. A structured rewrite is needed.
