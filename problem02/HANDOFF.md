# Problem 02 Handoff — Whittaker Functions for Rankin–Selberg Integrals

## What This Problem Asks

Let F be a non-archimedean local field, Π a generic irreducible admissible representation of GL_{n+1}(F) in its ψ⁻¹-Whittaker model. For each generic irreducible admissible π of GL_n(F) with conductor ideal 𝔮, set u_Q = I_{n+1} + Q·E_{n,n+1} where Q generates 𝔮⁻¹. **Must there exist a single W ∈ W(Π, ψ⁻¹) such that for every such π, some V ∈ W(π, ψ) makes the local Rankin–Selberg integral finite and nonzero for all s ∈ ℂ?**

## Answer

**YES.** The essential Whittaker function (new vector) W° of Π works.

## Current State (Session 2)

- **Proof tree:** 17 nodes (expanded from 7 in session 1)
- **Validated:** 2 (Nodes 1.1, 1.2)
- **Pending verification:** 15 new/amended nodes
- **Challenges:** 14 raised in session 2, ALL 14 RESOLVED by provers
- **Critical path:** `1 → 1.4 → 1.4.4 → 1.4.4.2` (depth 4)
- **Run `af status`** from `problem02/` to see the full tree.

## Session 2 Summary

### Verification Pass (breadth-first on all 7 original nodes)
- **Node 1** (root): Structurally sound, no challenges
- **Node 1.1** (Commutation Identity): VALIDATED ✓
- **Node 1.2** (Algebraic Characterization): VALIDATED ✓
- **Node 1.3** (Test Vectors): 1 critical challenge — unramified case not handled
- **Node 1.4** (Iwasawa Unfolding): 4 challenges — torus collapse fails for n≥2, K-integral not a Gauss sum, missing deps, incomplete justification
- **Node 1.5** (Gauss Sum): 5 challenges — undefined objects, type errors, wrong mechanism for n≥2
- **Node 1.6** (QED): 4 challenges — false universality, missing deps, incomplete factorization, k never given

### Prover Wave (addressed all 14 challenges)
All four challenged nodes were substantially rewritten:

**Node 1.3** → Case split on conductor:
- Ramified (c(π) ≥ 1): V = V° (new vector)
- Unramified (c(π) = 0): V = V_0 (compact Kirillov support)
- New children: 1.3.1 (V_0 construction), 1.3.2 (monomial proof for unramified case)

**Node 1.4** → "Iwasawa Unfolding and Conductor Analysis (Ramified Case)":
- Dropped "Gauss sum" terminology; K-integral now "matrix-coefficient integral with partial additive twist"
- W-factorization proved explicitly
- New children: 1.4.1 (W-factorization), 1.4.2 (case a vanishing), 1.4.3 (conductor analysis), 1.4.4 (torus sum reduction), 1.4.4.1 (supercuspidal: single torus point), 1.4.4.2 (non-supercuspidal ramified: epsilon factor approach)

**Node 1.5** → "Nonvanishing of Surviving Terms (Ramified Case)":
- Removed undefined τ(π,ψ_Q) and "multiplicative character of π"
- New children: 1.5.1 (K-integral nonvanishing, supercuspidal), 1.5.2 (epsilon factor nonvanishing, non-supercuspidal)

**Node 1.6** → Three-case QED assembly with explicit dependencies, complete factorization, and explicit k(Π,π) values

## Proof Tree Structure

```
1 [pending] Root conjecture
├── 1.1 [VALIDATED] Commutation identity
├── 1.2 [VALIDATED] Algebraic characterization (monomial iff)
├── 1.3 [pending] Test vector choice (case split)
│   ├── 1.3.1 [pending] Unramified test vector V_0 construction
│   └── 1.3.2 [pending] Monomial proof for unramified case
├── 1.4 [pending] Iwasawa unfolding (ramified case)
│   ├── 1.4.1 [pending] W-factorization
│   ├── 1.4.2 [pending] Case (a) vanishing
│   ├── 1.4.3 [pending] Conductor analysis + fiber decomposition
│   └── 1.4.4 [pending] Torus sum reduction
│       ├── 1.4.4.1 [pending] Supercuspidal: single torus point
│       └── 1.4.4.2 [pending] Non-supercuspidal ramified (epsilon factor) ← CRITICAL PATH
├── 1.5 [pending] Nonvanishing of surviving terms (ramified)
│   ├── 1.5.1 [pending] K-integral nonvanishing (supercuspidal)
│   └── 1.5.2 [pending] Epsilon factor nonvanishing (non-supercuspidal)
└── 1.6 [pending/qed] Three-case conclusion
```

## What the Next Agent Should Do

### Priority 1: Verify the new child nodes (breadth-first)
15 pending nodes need adversarial verification. Run verifiers breadth-first:
- Depth 2: 1.3.1, 1.3.2, 1.4.1, 1.4.2, 1.4.3, 1.4.4, 1.5.1, 1.5.2
- Depth 3: 1.4.4.1, 1.4.4.2
- Then re-verify amended parents: 1.3, 1.4, 1.5, 1.6
- Finally: root node 1

### Priority 2: Scrutinize Node 1.4.4.2 (critical path)
The non-supercuspidal ramified case is the weakest link. The prover acknowledged:
- The connection between the u_Q-twisted integral and ε(s, Π×π, ψ) is "expected but not fully verified"
- An alternative modified test vector approach is proposed but deferred
- This node needs the most rigorous verification

### Priority 3: Scrutinize Node 1.3.2 (unramified monomial proof)
The argument that the twist vanishes and the integral collapses to W°(I)·∫_K V_0(k) dk needs checking:
- Property (P1): torus support collapse via intersection of Casselman-Shalika and Kirillov supports
- Property (P3): positivity of the K-integral at a = I_n

### Known Remaining Weaknesses
1. **Node 1.4.4.2:** Non-supercuspidal ramified case — epsilon factor identification not fully proved
2. **Node 1.5.2:** Depends on 1.4.4.2's output being well-defined
3. **Dependency declarations:** Nodes 1.4 and 1.5 have textual (not metadata) dependency declarations due to af tool limitations

## Definitions in Scope

- non_archimedean_local_field, generic_representation, Whittaker_model
- conductor_ideal, Rankin_Selberg_integral, upper_triangular_unipotent

## External References

| Name | Source |
|------|--------|
| Whittaker models | Cogdell, Piatetski-Shapiro (2004) |
| Rankin-Selberg theory | Jacquet, Piatetskii-Shapiro, Shalika (1983) |
| JPSS Rankin-Selberg | JPSS (1983), Amer. J. Math. 105(2), 367-464 |
| Matringe essential Whittaker | Matringe (2013), Doc. Math. 18, 1191-1214 |
| Miyauchi newforms | Miyauchi (2014), J. Math. Soc. Japan 66, 17-24 |
| Godement-Jacquet epsilon factors | Godement-Jacquet (1972), LNM 260 |
| Tate thesis | Tate (1950), Princeton PhD thesis |

## af Tool Quick Reference

```bash
af status                          # see tree
af get <id>                        # node details
af claim <id> --owner X --role Y   # claim a node
af refine <id> --owner X -s "..."  # add proof content
af release <id> --owner X          # release claim
af defs                            # list definitions
af externals                       # list references
af challenges                      # list all challenges
af jobs                            # see available work
```

Note: Run all commands from the `problem02/` directory.
