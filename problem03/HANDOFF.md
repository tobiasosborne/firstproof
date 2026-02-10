# Problem 03 Handoff — Markov Chain with ASEP Polynomial Stationary Distribution

## What This Problem Asks

Let λ = (λ₁ > ⋯ > λₙ ≥ 0) be a restricted partition with distinct parts (unique part 0, no part 1). Let S_n(λ) be the set of compositions of λ. Does there exist a **nontrivial** Markov chain on S_n(λ) whose stationary distribution is:

> π(μ) = F\*_μ(x₁,…,xₙ; q=1, t) / P\*_λ(x₁,…,xₙ; q=1, t)

where "nontrivial" means the transition probabilities are NOT described using the F\*_μ polynomials?

## Answer

**YES.** The inhomogeneous multispecies t-PushTASEP provides such a chain.

## Current State (Session 3 — Prover Wave Complete)

- **Proof tree:** 9 nodes (root + 8 children), all pending
- **0 VALIDATED**, 0 refuted, 0 archived
- **97 challenges filed in Session 2, ALL 97 RESOLVED in Session 3**
- **0 open challenges** — all nodes rewritten
- Run `af status` from this directory to see the full tree

### Session 3: Prover Wave Results

All 8 leaf nodes (1.1–1.8) were rewritten by independent prover subagents to address the 97 challenges from the Session 2 verification wave. Every systemic problem identified in Session 2 has been addressed:

| Node | Description | Challenges Resolved | Key Changes |
|------|-------------|-------------------|-------------|
| 1.1 | State space setup | 7 (4M, 3m) | Precise definitions: ring Z_n with periodic boundary, species labels {0}∪{k≥2}, one occupant per site, \|S_n(λ)\|=n! |
| 1.2 | Polynomial identification | 13 (2C, 8M, 3m) | Full notation table (BDW25/AMW24/CMW22), correct arXiv:2403.10485, six numbered claims with citations |
| 1.3 | Positivity & normalization | 11 (4C, 6M, 1m) | **Restructured as corollary** of 1.5+1.6. Circularity broken. Derives ratio positivity, not standalone F\*_μ≥0 |
| 1.4 | Chain construction | 17 (4C, 10M, 3m) | Full cascade mechanism, geometric displacement formula, vacancy behavior, irreducibility, CT→DT bridge |
| 1.5 | Stationarity | 18 (7C, 10M, 1m) | **Complete rewrite**: cites AMW24 Thm 1.1 correctly, algebraic proof method, no Ferrari-Martin fiction |
| 1.6 | Ratio identity (CRUX) | 16 (11C, 5M) | **Complete rewrite**: 3-step Hecke stationarity argument (BDW25 Prop 2.10 + AMW24 Thm 1.1 + Perron-Frobenius) |
| 1.7 | Nontriviality | 13 (7C, 5M, 1m) | Formal definition of "nontrivial" (≠ trivial i.i.d. sampler), rigorous sparsity proof with n=3 computation |
| 1.8 | Conclusion | 2 (1C, 1M) | Full dependency chain declared, logical synthesis of all nodes, parameter domains |

### Systemic Problems from Session 2 — Resolution Status

1. **Wrong AMW24 reference** → FIXED. All nodes now cite arXiv:2403.10485 (Ayyer-Martin-Williams). The external reference entry still has the old arXiv (af tool lacks update-external command), but all node statements cite the correct paper inline.

2. **Node 1.5 fictional proof method** → FIXED. Ferrari-Martin claim completely removed. Now correctly cites AMW24 Theorem 1.1 with algebraic/combinatorial proof method.

3. **Node 1.6 logical fallacy** → FIXED. "Both sum to 1 ∴ equal" replaced by 3-step argument: (1) Hecke relations for f\*_μ (BDW25 Prop 2.10), (2) balance equation transfer from f_μ to f\*_μ, (3) Perron-Frobenius uniqueness forces proportionality.

4. **Circular dependency** → FIXED. Node 1.3 restructured as corollary of 1.5+1.6. Dependency flow: 1.1→1.2→1.4→1.5→1.6→1.3→1.7→1.8.

5. **Notation inconsistency** → FIXED. Node 1.2 now contains explicit cross-paper notation table. All nodes use consistent BDW25 notation (f\*_μ, F\*_μ, f_μ).

6. **No dependencies declared** → FIXED. Every node now explicitly declares its dependencies.

### Proof Tree Structure (Post-Rewrite)

| Node | Description | Status | Dependencies | Risk |
|------|-------------|--------|-------------|------|
| 1 | Root conjecture | pending | 1.1-1.8 | — |
| 1.1 | State space setup | REWRITTEN (0 open) | — | **Low** |
| 1.2 | Polynomial identification | REWRITTEN (0 open) | 1.1 | **Low** |
| 1.3 | Positivity & normalization | REWRITTEN (0 open) | 1.5, 1.6 | **Low** (now a corollary) |
| 1.4 | Chain construction | REWRITTEN (0 open) | 1.1 | **Medium** — cascade details may need verification |
| 1.5 | Stationarity (AMW24 Thm 1.1) | REWRITTEN (0 open) | 1.2, 1.4 | **Medium** — external theorem citation |
| 1.6 | Ratio identity (CRUX) | REWRITTEN (0 open) | 1.2, 1.4, 1.5 | **High** — Step 2 (Hecke transfer) is nontrivial; see caveat below |
| 1.7 | Nontriviality | REWRITTEN (0 open) | 1.1, 1.2, 1.4, 1.5 | **Low** — sparsity argument is concrete |
| 1.8 | Conclusion | REWRITTEN (0 open) | 1.1-1.7 | **Low** — synthesis node |

### Critical Caveat: Node 1.6 Step 2

The ratio identity proof (Node 1.6) hinges on **Step 2**: the claim that the AMW24 balance equation proof for f_μ transfers to f\*_μ because both satisfy the same Hecke relations (BDW25 Prop 2.10). This is the deepest mathematical claim in the proof. BDW25 Remark 1.17 defers the interpolation probabilistic interpretation to a forthcoming paper [BDW] that has not appeared, suggesting the experts consider this transfer nontrivial. The node acknowledges this caveat and identifies BDW25 Theorem 7.1 (factorization) as an alternative route. A verifier should scrutinize Step 2 with particular care.

### Next Steps (Session 4)

**Verification wave 2**: All 8 leaf nodes are ready for re-verification. Priority:

1. **Node 1.6** (highest risk) — Scrutinize the Hecke stationarity transfer argument (Step 2)
2. **Node 1.4** — Verify cascade mechanism matches AMW24 precisely
3. **Node 1.5** — Verify AMW24 Theorem 1.1 is cited faithfully
4. **Nodes 1.1, 1.2, 1.3, 1.7, 1.8** — Should be relatively quick to verify

After verification wave 2, address any new challenges, then attempt to validate leaf nodes and work upward to root.

### Key References (Verified)

| Tag | Paper | arXiv | Role in proof |
|-----|-------|-------|---------------|
| AMW24 | Ayyer-Martin-Williams | **2403.10485** | t-PushTASEP stationarity (Thm 1.1), irreducibility (Prop 2.4) |
| BDW25 | Ben Dali-Williams | 2510.02587 | Interpolation ASEP polys: Def 1.2, Thm 1.15, Prop 2.10, Prop 2.15, Thm 2.3, Thm 7.1, Remark 1.17 |
| CMW22 | Corteel-Mandelshtam-Williams | 1811.01024 | ASEP polynomials via multiline queues |
| FM07 | Ferrari-Martin | Ann. Probab. 35(3) | **Ordinary TASEP only** — does NOT apply to t-PushTASEP |
| KS96 | Knop (1997), Sahi (1996) | — | Interpolation Macdonald polynomials |

### Lessons Learned (Session 3)

1. **The Hecke stationarity argument (Node 1.6) is a plausible original mathematical argument** — not just a citation of existing work. It uses BDW25 Prop 2.10 + AMW24 Thm 1.1 + Perron-Frobenius in a novel combination. This makes it the most vulnerable node to further challenges.
2. **Restructuring 1.3 as a corollary of 1.5+1.6 elegantly breaks the circularity** without needing to prove the open problem of F\*_μ nonnegativity directly.
3. **The nontriviality argument via sparsity (Node 1.7) is clean and concrete** — exhibiting a zero in the transition matrix is much stronger than informal "local vs global" heuristics.
4. **Notation discipline pays off** — the cross-paper notation table in Node 1.2 makes all subsequent nodes readable.
5. **Provers should always resolve ALL challenges**, not just critical ones. Leaving minor challenges open accumulates technical debt.
