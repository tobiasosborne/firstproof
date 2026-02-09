# Problem 02 Handoff — Whittaker Functions for Rankin–Selberg Integrals

## What This Problem Asks

Let F be a non-archimedean local field, Π a generic irreducible admissible representation of GL_{n+1}(F) in its ψ⁻¹-Whittaker model. For each generic irreducible admissible π of GL_n(F) with conductor ideal 𝔮, set u_Q = I_{n+1} + Q·E_{n,n+1} where Q generates 𝔮⁻¹. **Must there exist a single W ∈ W(Π, ψ⁻¹) such that for every such π, some V ∈ W(π, ψ) makes the local Rankin–Selberg integral finite and nonzero for all s ∈ ℂ?**

## Answer

**YES.** The essential Whittaker function (new vector) W° of Π works.

## Current State (Session 1)

- **Proof tree:** 7 nodes (1 root + 6 children)
- **All PENDING** — no adversarial verification yet
- **5 external references** added
- **6 definitions** inherited from initialization
- **Run `af status`** from `problem02/` to see the full tree.

## Proof Strategy

The proof has 6 steps, registered as nodes 1.1–1.6:

### Node 1.1 — Commutation Identity (algebraic reduction)
For g ∈ GL_n(F): W(diag(g,1)·u_Q) = ψ⁻¹(Q·g_{nn})·W(diag(g,1)).

This is the foundational reduction. Conjugating u_Q past diag(g,1) produces n'(g,Q) ∈ N_{n+1} whose only superdiagonal entry contributing to ψ⁻¹ is Q·g_{nn} at position (n,n+1). The twist ψ⁻¹(Q·g_{nn}) is left-N_n-invariant, so the integrand remains well-defined.

**Difficulty:** Low. Pure matrix algebra + definition of Whittaker model.

### Node 1.2 — Algebraic Characterization of "Finite and Nonzero ∀s"
A rational function R(q⁻ˢ) is finite and nonzero for all s ∈ ℂ iff R is a nonzero monomial c·q⁻ᵏˢ.

Combined with JPSS theory (the integral is rational in q⁻ˢ), this sets the concrete target: show the integral is a nonzero monomial.

**Difficulty:** Low. Standard algebra of Laurent polynomials.

### Node 1.3 — Choice of Test Vectors
Set W = W° (essential Whittaker function of Π, depends only on Π). For each π, set V = V° (essential Whittaker function of π). Key property: W°(I_{n+1}) ≠ 0 (Matringe-Miyauchi).

**Difficulty:** Low. Definitions + citation of Matringe (2013) and Miyauchi (2014).

### Node 1.4 — Iwasawa Unfolding and Gauss Sum Extraction (CORE COMPUTATION)
Using Iwasawa decomposition g = nak, the integral unfolds to a torus sum. The additive twist from Node 1.1 combined with K₁(𝔮)-invariance of V° produces a K-integral that is a generalized Gauss sum. Conductor-level analysis:
- val(Q·ϖ^{m_n}) > 0 → twist trivial → K-integral vanishes (orthogonality)
- val(Q·ϖ^{m_n}) < −c(π) → oscillation too rapid → K-integral vanishes (cancellation)
- val(Q·ϖ^{m_n}) = −c(π) → conductors match → K-integral = nonzero Gauss sum

Only the matched level survives → torus sum collapses to a single term → monomial.

**Difficulty:** HIGH. This is the technical heart. Requires careful Iwasawa computation, tracking of support conditions for W° and V° on the torus, and the three-case conductor analysis.

### Node 1.5 — Gauss Sum Nonvanishing
The Gauss sum τ(π, ψ_Q) ≠ 0 because conductor of ψ_Q(x) = ψ(Qx) exactly matches conductor of π. Classical result (Tate thesis for GL₁; Godement-Jacquet for GL_n). The torus value W°(diag(a,1)) is nonzero by Matringe-Miyauchi.

**Difficulty:** Medium. Requires precise statement and citation of epsilon factor theory.

### Node 1.6 — QED: Conclusion and Uniformity
Assembly: integral = c(Π,π)·q⁻ᵏ⁽ᐩ'ᵖⁱ⁾ˢ with c ≠ 0. This is a nonzero monomial (Node 1.2), hence finite and nonzero for all s. W° depends only on Π (Node 1.3); conductor matching is automatic from the definition of Q.

**Difficulty:** Low (given Nodes 1.1–1.5).

## What the Next Agent Should Do

### Priority 1: Refine Node 1.1 (Commutation Identity)
This is the easiest node and provides the foundation for everything else. Write out the explicit matrix computation. Verify left-N_n-invariance of the twist.

### Priority 2: Refine Node 1.4 (Iwasawa Unfolding) — CRITICAL
This is the hardest step. The agent should:
1. Write the full Iwasawa decomposition of the integral
2. Substitute the commutation identity from 1.1
3. Perform the K-integral explicitly, separating into the three conductor cases
4. Show the torus sum collapses

**Known vulnerability:** The torus-sum collapse may not be as clean for general n as it is for n = 1. The GL(2)×GL(1) case is straightforward (see backup below). For general n, the support properties of W° and V° on the torus (Matringe-Miyauchi formulas) must be invoked carefully.

### Priority 3: Refine Node 1.5 (Gauss Sum Nonvanishing)
State the precise epsilon factor identity. Cite Godement-Jacquet or Tate.

### Backup Strategy (if new vector W° fails for general n)
If the new vector W° does not produce a monomial for all π, use the **Kirillov model approach**: choose W with Kirillov-model restriction φ = 1_{𝔬ⁿ} (characteristic function of 𝔬ⁿ). This lies in K(Π) for any generic Π by the Bernstein-Zelevinsky embedding S(Fⁿ \ {0}) ↪ K(Π). The compact support forces a finite torus sum, and V can be chosen (depending on π) to isolate a single level.

### Unramified Case (c(π) = 0) — Special Handling
When π is unramified, Q is a unit and the twist ψ⁻¹(Q·g_{nn}) is trivial on GL_n(𝔬). The integral reduces to the standard RS integral. If W = W° produces L(s, Π×π) (which may have poles), choose V ≠ V° from W(π, ψ) with compact Kirillov support to get a monomial instead.

## Key Insights from Strategy Development

Three independent strategies were evaluated:

1. **Strategy A (L-functions + essential vectors):** Most concrete. Proposes W = W° and explicit Matringe-Miyauchi computation. Risk: torus collapse hand-waved for general n.

2. **Strategy B (Kirillov model + Gauss sums):** Best foundational computation (the commutation identity). Uses Kirillov model rather than new vector, giving more flexibility. Risk: Claim 6.2 for n ≥ 2 unfinished.

3. **Strategy C (Bernstein center + Baire category):** Most abstract. Existential proof via distributional argument. Risk: multiple gaps, lowest confidence (65%).

The synthesized strategy takes the **commutation identity from B**, the **explicit new-vector choice from A**, and keeps the **Kirillov model fallback from B/C** in reserve.

## Key Pitfall: Why u_Q Is Not Just a Right-Translation

It is tempting to treat the u_Q-modified integral as simply a standard RS integral with R(u_Q)W in place of W. While formally correct, this obscures the mechanism: R(u_Q)W depends on Q, which depends on π. The commutation identity (Node 1.1) is more useful because it separates the pi-dependent part (the additive twist ψ⁻¹(Q·g_{nn})) from the pi-independent part (W(diag(g,1))).

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
```

Note: Run all commands from the `problem02/` directory.
