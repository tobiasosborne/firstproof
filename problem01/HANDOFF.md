# Problem 01 Handoff — Φ⁴₃ Measure Equivalence Under Smooth Shifts

## What This Problem Asks

Let μ be the Φ⁴₃ measure on D'(T³) (three-dimensional torus), ψ : T³ → ℝ smooth and nonzero, T_ψ(u) = u + ψ the shift map. **Are μ and T_ψ*μ equivalent (mutually absolutely continuous)?**

## Answer

**YES.** The measures are equivalent.

## Current State

- **Proof tree:** 21 nodes, 3 levels deep, all `pending` verification. Registered in the `af` workspace in this directory.
- **No proofs written yet** — only the claim structure and strategy are in place.
- **Run `af status -d .`** from this directory to see the full tree.

## Proof Strategy Summary

The proof constructs an explicit Radon-Nikodym derivative d(T_ψ*μ)/dμ via UV regularization and a limiting argument. It has four stages:

### Stage A: Regularized RN Derivative (Nodes 1.2, 1.3)

For the UV-regularized measure μ_ε = Z_ε⁻¹ exp(-V_ε(φ)) dμ₀ (where μ₀ is the GFF), we have μ_ε ~ μ₀ trivially (explicit positive density) and T_ψ*μ₀ ~ μ₀ by Cameron-Martin (ψ ∈ C^∞ ⊂ H¹ = CM space). So T_ψ*μ_ε ~ μ_ε with explicit RN derivative:

```
R_ε(φ) = exp(-ΔV_ε(φ) + ⟨(-Δ+m²)ψ, φ⟩ - ½‖ψ‖²_{H¹})
```

where ΔV_ε = V_ε(φ−ψ) − V_ε(φ).

### Stage B: Wick Expansion (Nodes 1.4, 1.4.1–1.4.4)

The interaction difference decomposes by Wick calculus:

```
ΔV_ε(φ) = −4λ ∫ψ_ε :φ_ε³: dx          (cubic, converges)
          + 6λ ∫ψ_ε² :φ_ε²: dx          (quadratic, converges)
          + (−4λ⟨ψ_ε³,φ_ε⟩ − 2δm_ε²⟨ψ_ε,φ_ε⟩)  (linear, DIVERGES via mass counterterm)
          + constants
```

**Critical point:** The Wick constant C_ε(0) cancels exactly in the quartic shift — no new UV divergences arise in the Wick-power terms. The ONLY divergent field-dependent term is `2δm_ε²⟨ψ,φ⟩` (linear in φ), coming from the mass counterterm δm_ε² ~ λ² log(1/ε).

### Stage C: Renormalization and Convergence (Nodes 1.5, 1.6)

Split the exponent: Ψ_ε = Ψ_ε^{ren} + L_ε + K_ε where:
- **Ψ_ε^{ren}:** convergent Wick-power terms → limit Ψ^{ren}(φ) under μ
- **L_ε = 2δm_ε²⟨ψ,φ⟩:** divergent linear term
- **K_ε:** divergent constant, fixed by normalization E_{μ_ε}[R_ε] = 1

The normalization constraint absorbs the divergence. Key estimates needed:
1. Smeared Wick powers ∫ψ:φ³:dx, ∫ψ²:φ²:dx are well-defined under μ (regularity theory)
2. They have **finite exponential moments** under μ (sub-Gaussian tails, Barashkov-Gubinelli)
3. {R_ε} is **uniformly integrable** w.r.t. {μ_ε} (de la Vallée-Poussin + Hölder)
4. R_ε → R in L¹(μ) with E_μ[R] = 1

### Stage D: Conclusion (Nodes 1.7, 1.8)

- R = exp(real-valued r.v.) > 0 μ-a.s. → T_ψ*μ ≪ μ
- Replace ψ → −ψ → μ ≪ T_ψ*μ
- Therefore T_ψ*μ ~ μ. QED.

## What the Next Agent Should Do

The proof tree is **claims only** — every node is `pending`. The next step is to **fill in the actual proofs** for each node, working bottom-up from the leaves. Priority order:

1. **Node 1.4.1–1.4.2** (Wick shift algebra) — These are straightforward algebraic computations. Verify and accept.
2. **Node 1.4.3** (Full interaction difference) — Combine 1.4.1 + 1.4.2. Straightforward.
3. **Node 1.4.4** (UV divergence analysis) — Needs regularity estimates for Wick powers on T³.
4. **Node 1.5.1–1.5.3** (Exponent decomposition and convergence) — Mostly bookkeeping once 1.4 is done.
5. **Node 1.5.4** (Absorption of divergent linear term) — **Key technical node.** Must show the normalized R_ε stays bounded in L^p despite individual terms diverging.
6. **Node 1.6.1–1.6.4** (Convergence) — **Hardest part.** Needs:
   - Regularity of :φ³:, :φ²: under μ (cite Hairer/GIP)
   - Exponential integrability (cite Barashkov-Gubinelli sub-Gaussian tails)
   - Uniform integrability argument
   - Weak convergence + UI → L¹ convergence
7. **Nodes 1.7, 1.8** (Positivity + symmetry) — Easy once 1.6 is done.

## Key Pitfall to Watch

**μ ⊥ μ₀** (the Φ⁴₃ measure is SINGULAR w.r.t. the GFF). This was proved by Barashkov-Gubinelli (2021) and Hairer. The naive argument "μ ~ μ₀ ~ T_ψ*μ₀, so μ ~ T_ψ*μ" is **wrong**. The proof must work directly with the regularized measures and take limits, never assuming μ ~ μ₀.

## Definitions and References in the af Workspace

**11 definitions:** `Phi43_measure`, `shift_map`, `pushforward_measure`, `measure_equivalence`, `three_torus`, `Wick_ordering`, `mass_counterterm`, `Radon_Nikodym_derivative`, `sub_Gaussian_tails`, `UV_regularization`, `Cameron_Martin_space`

**6 external references:**
- Hairer (2014) — Regularity structures / Φ⁴₃ construction
- Gubinelli-Imkeller-Perkowski (2015) — Paracontrolled distributions
- Barashkov-Gubinelli (2020) — Variational method for Φ⁴₃ (Duke Math. J.)
- Barashkov-Gubinelli (2021) — Girsanov construction + singularity proof (EJP)
- Barashkov-Gubinelli (2021) — Sub-Gaussian tails
- Bogachev (1998) — Cameron-Martin theorem

## af Tool Quick Reference

```bash
af status -d /home/tobiasosborne/Projects/firstproof/problem01   # see tree
af get <id> -d ...                                                # node details
af claim <id> --owner prover-1 --role prover -d ...               # claim a node
af refine <id> --owner prover-1 -s "proof text" -d ...            # add proof content
af release <id> --owner prover-1 -d ...                           # release claim
af defs -d ...                                                    # list definitions
af externals -d ...                                               # list references
```

Note: `af types` and `af inferences` do NOT accept `-d`; you must `cd` into the directory first.
