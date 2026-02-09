# Problem 02 Handoff — Whittaker Functions for Rankin–Selberg Integrals

## What This Problem Asks

Let F be a non-archimedean local field, Π a generic irreducible admissible representation of GL_{n+1}(F) in its ψ⁻¹-Whittaker model. For each generic irreducible admissible π of GL_n(F) with conductor ideal 𝔮, set u_Q = I_{n+1} + Q·E_{n,n+1} where Q generates 𝔮⁻¹. **Must there exist a single W ∈ W(Π, ψ⁻¹) such that for every such π, some V ∈ W(π, ψ) makes the local Rankin–Selberg integral finite and nonzero for all s ∈ ℂ?**

## Answer

**YES.** The essential Whittaker function (new vector) W° of Π works.

## Current State (Session 3 — INTERRUPTED)

- **Proof tree:** 17 nodes
- **Validated:** 4 (Nodes 1.1, 1.2, 1.4.1, 1.4.2)
- **Challenged (wave 2):** 6 nodes with 20 new open challenges
- **Not yet verified (wave 2):** 7 nodes remaining
- **Run `af status`** from `problem02/` to see the full tree.
- **Run `af challenges`** to see all 20 open challenges.

## Session 3 Summary (Verification Wave 2 — PARTIAL)

### Completed Verifications (8 of 15 nodes)

| Node | Depth | Verdict | Challenges | Key Issues |
|------|-------|---------|------------|------------|
| 1.3.1 | 2 | CHALLENGED | 2 critical, 1 major | Kirillov coords false off mirabolic; Casselman-Shalika scope; K-projection dilemma |
| 1.3.2 | 2 | CHALLENGED | 2 critical, 1 major, 1 minor | Same Kirillov/CS issues as 1.3.1; incoherent φ₀ definition |
| 1.4.1 | 2 | **ACCEPTED** | — | W° factorization is clean and correct |
| 1.4.2 | 2 | **ACCEPTED** | — | Case (a) vanishing is rigorous |
| 1.4.3 | 2 | CHALLENGED | 1 critical, 2 major, 1 minor | GL_n(o/q) factorization invalid for m_n<0; missing step for 1≤m_n≤c-1; missing case m_n=c |
| 1.4.4 | 2 | CHALLENGED | 1 critical, 3 major, 1 minor | Zero analysis invalid (Laurent poly HAS zeros); pole analysis hand-wavy; Schur poly wrong for ramified |
| 1.5.1 | 2 | CHALLENGED | 1 critical, 1 major, 1 minor | Nonvanishing circular for n≥2; fiber sum domain wrong |
| 1.5.2 | 2 | CHALLENGED | 2 critical, 1 major, 2 minor | Test vector lit doesn't cover general case; "nonzero at one point ≠ monomial" error |

### Not Yet Verified (7 nodes remaining)
- **Depth 3:** 1.4.4.1 (supercuspidal torus collapse), 1.4.4.2 (non-supercuspidal epsilon)
- **Depth 1 parents:** 1.3, 1.4, 1.5, 1.6
- **Root:** 1

## What the Next Agent Should Do

### Priority 1: FINISH Verification Wave 2 (7 remaining nodes)
Continue the breadth-first verification sweep. Order:
1. **1.4.4.1** (depth 3) — supercuspidal single torus point. Key question: is φ° really supported on a single torus coset for GL_n, n≥3?
2. **1.4.4.2** (depth 3) — non-supercuspidal ramified epsilon factor. CRITICAL PATH.
3. **1.3** (depth 1) — re-verify amended parent (test vector case split)
4. **1.4** (depth 1) — re-verify amended parent (Iwasawa unfolding)
5. **1.5** (depth 1) — re-verify amended parent (nonvanishing)
6. **1.6** (depth 1) — re-verify QED node
7. **1** (depth 0) — root node

Use the same verifier protocol: `af claim <id> --owner verifier-w2 --role verifier`, then accept/challenge, then `af release <id> --owner verifier-w2`.

### Priority 2: Launch Prover Wave for ALL Challenged Nodes
After verification wave 2 completes, launch provers for all nodes with open challenges. The 6 challenged nodes so far are: **1.3.1, 1.3.2, 1.4.3, 1.4.4, 1.5.1, 1.5.2** (plus any from the remaining 7 verifications).

### Systematic Issues to Address (Recurring across multiple nodes)

1. **Kirillov coordinate identification** (1.3.1, 1.3.2): The Kirillov model realizes π|_{P_n} on functions on F^{n-1}\{0}, where P_n is the mirabolic. For elements k ∈ K_n NOT in P_n (positive measure subset), evaluating V_0(ak) requires the full representation action π(k), not pointwise evaluation of φ₀. This invalidates the entire unramified torus collapse for n≥2.

2. **Casselman-Shalika scope** (1.3.1, 1.3.2): The Casselman-Shalika formula applies only to the spherical Whittaker function of UNRAMIFIED Π. For ramified Π, W° has different torus support (Matringe formulas). The dominance condition m_i ≥ 0 is not established for all generic Π.

3. **K-projection dilemma** (1.3.2): For unramified π, the K-average of V_0 is proportional to the spherical vector V°. If nonzero, the integral reproduces L(s, Π×π) which has poles. If zero, nonvanishing fails. Fundamental obstruction.

4. **Test vector literature** (1.5.2): Humphries (2021) is GL_2-specific. Assing-Blomer (2024) doesn't establish the universal identity for GL_{n+1}×GL_n. The identification I(s) = C·ε(s, Π×π, ψ) is unproven.

5. **J_K(0) nonvanishing** (1.5.1): The "maximal Fourier content at conductor level" claim has no proof for n≥2. Need either: Bushnell-Kutzko type theory, finite group rep theory on GL_n(F_q), or test vector theorems.

## Proof Tree Structure

```
1 [pending] Root conjecture
├── 1.1 [VALIDATED] Commutation identity
├── 1.2 [VALIDATED] Algebraic characterization (monomial iff)
├── 1.3 [pending, NOT YET RE-VERIFIED] Test vector choice (case split)
│   ├── 1.3.1 [CHALLENGED, 3 open] Unramified test vector V_0 construction
│   └── 1.3.2 [CHALLENGED, 4 open] Monomial proof for unramified case
├── 1.4 [pending, NOT YET RE-VERIFIED] Iwasawa unfolding (ramified case)
│   ├── 1.4.1 [VALIDATED] W-factorization ✓
│   ├── 1.4.2 [VALIDATED] Case (a) vanishing ✓
│   ├── 1.4.3 [CHALLENGED, 4 open] Conductor analysis + fiber decomposition
│   └── 1.4.4 [CHALLENGED, 5 open] Torus sum reduction
│       ├── 1.4.4.1 [pending, NOT YET VERIFIED] Supercuspidal: single torus point
│       └── 1.4.4.2 [pending, NOT YET VERIFIED] Non-supercuspidal ramified ← CRITICAL PATH
├── 1.5 [pending, NOT YET RE-VERIFIED] Nonvanishing of surviving terms (ramified)
│   ├── 1.5.1 [CHALLENGED, 3 open] K-integral nonvanishing (supercuspidal)
│   └── 1.5.2 [CHALLENGED, 5 open] Epsilon factor nonvanishing (non-supercuspidal)
└── 1.6 [pending, NOT YET RE-VERIFIED] Three-case conclusion
```

## Session History

### Session 1
- Initial 7 nodes created
- Nodes 1.1, 1.2 validated

### Session 2
- Verification wave 1: 14 challenges raised on nodes 1.3, 1.4, 1.5, 1.6
- Prover wave: All 14 challenges resolved; 10 new child nodes created
- Tree expanded to 17 nodes

### Session 3 (current, interrupted)
- Verification wave 2: 8 of 15 nodes verified
- 2 newly validated (1.4.1, 1.4.2), 6 newly challenged (20 open challenges)
- 7 nodes remain to be verified before prover wave can begin

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
