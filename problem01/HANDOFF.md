# Problem 01 Handoff — Φ⁴₃ Measure Equivalence Under Smooth Shifts

## What This Problem Asks

Let μ be the Φ⁴₃ measure on D'(T³) (three-dimensional torus), ψ : T³ → ℝ smooth and nonzero, T_ψ(u) = u + ψ the shift map. **Are μ and T_ψ*μ equivalent (mutually absolutely continuous)?**

## Answer

**YES.** The measures are equivalent.

## Current State (Session 2)

- **Proof tree:** 59 nodes total (expanded from original 21 by prover refinements)
- **14 VALIDATED** by adversarial verification
- **4 REFUTED** (all repaired or in repair)
- **5 ARCHIVED** (superseded by repairs)
- **36 PENDING** (including child proof nodes, parent wrappers, and remaining work)
- **Run `af status -d .`** from this directory to see the full tree.

### Validated Nodes (Adversarially Verified)

| Node | Description | Status |
|------|-------------|--------|
| 1.2 | Setup: μ_ε ~ μ₀, Cameron-Martin, T_ψ*μ_ε ~ μ_ε | VALIDATED |
| 1.3 | Regularized RN derivative formula | VALIDATED |
| 1.4.1 | Quartic Wick shift (algebraic identity) | VALIDATED |
| 1.4.2 | Quadratic Wick shift (algebraic identity) | VALIDATED |
| 1.4.3 | Full interaction difference (combining 1.4.1 + 1.4.2) | VALIDATED |
| 1.4.4 | UV divergence analysis (regularity estimates) | VALIDATED |
| 1.5.1 | Decomposition of exponent: Ψ_ε = Ψ^ren + L_ε + K_ε | VALIDATED |
| 1.5.2 | Normalization constraint: K_ε absorbs divergence | VALIDATED |
| 1.5.3 | Convergence of renormalized part Ψ^ren in L^p(μ) | VALIDATED |
| 1.5.4 | Absorption of divergent linear term (corrected) | VALIDATED |
| 1.6.1 | Well-definedness of smeared Wick powers under μ (corrected) | VALIDATED |
| 1.6.2.5 | k=1 exponential integrability (sub-Gaussian, BG 2020) | VALIDATED |

### Refutation History (HIGH-VALUE FINDINGS)

#### Refutation 1: Node 1.5.4 (Absorption of divergent linear term)
- **Original claim:** R_ε remains bounded in L^p(μ) for suitable p
- **Error found:** For any fixed p > 1, E_{μ_ε}[R_ε^p] ~ exp(C(log 1/ε)²) → ∞. Uniform integrability was claimed but never proven.
- **Repair:** Scoped to per-ε L^p only, deferred UI to node 1.6. Re-verified and VALIDATED.

#### Refutation 2: Node 1.6.1 (Well-definedness of Wick powers)
- **Errors found:** (a) Wrong duality condition (α+β > 0 vs α > β), (b) Borderline Besov pairing (regularity sum = 0), (c) Reversed Besov inclusion direction
- **Repair:** Fixed all three with correct indices and relabeling. Re-verified and VALIDATED.

#### Refutation 3: Node 1.6.2 (Exponential integrability, attempt 1)
- **Error found:** Decomposed Wick powers into raw powers (introducing divergent C_ε(0) ~ 1/ε corrections), applied Fatou's lemma across changing measures μ_ε without coupling.
- **Root cause:** Wick-to-raw decomposition destroys cancellations; Fatou requires fixed measure.

#### Refutation 4: Node 1.6.2 (Exponential integrability, attempt 2)
- **Error found:** Claimed BG (2020) sub-Gaussian concentration extends to enhanced data (φ, :φ²:, :φ³:). BG only proves concentration for φ itself.
- **Root cause:** Enhanced data lift φ → (φ, :φ²:, :φ³:) is not Lipschitz; Malliavin gradient of ⟨f, :φ^k:⟩ is φ-dependent and unbounded.

## What the Next Agent Should Do

### IMMEDIATE: Complete Node 1.6.2 (Exponential Integrability) — CRITICAL BLOCKER

A **third repair attempt** (node 1.6.2.10) is submitted and AWAITING VERIFICATION. Strategy:
1. **A priori estimate:** ‖:φ_ε^k:‖_{C^{-k/2-δ}} ≤ C(1 + V_ε)^{k/4} (treats Wick powers as single distributional objects)
2. **Young's inequality** with exponents (4/k, 4/(4-k)): quartic potential V_ε dominates degree-k perturbation since k ≤ 3 < 4
3. **Skorokhod coupling** on enhanced data space X (Polish) for μ_ε → μ
4. **Vitali convergence** on the common probability space

**Key vulnerability:** The a priori estimate (AP) in node 1.6.2.10.2 needs careful verification. It claims the Besov norm of :φ_ε^k: is controlled by V_ε^{k/4}. A verifier should check whether this is actually established in Mourrat-Weber (2017) or Gubinelli-Hofmanova (2019), or whether it conflates global and local norms.

**Alternative approaches if attempt 3 fails:**
- **Nelson hypercontractivity:** L^p moment growth rate under the Φ⁴₃ measure (with LSI from Bauerschmidt-Bodineau 2019) → exponential integrability for all α when k < 4
- **Direct variational formula:** Boue-Dupuis representation from BG (2020) applied to E_μ[exp(α⟨f,:φ^k:⟩)]
- **Lattice approximation:** Prove on finite lattice (where everything is finite-dimensional), then transfer to continuum with uniform bounds

### THEN: Complete Remaining 1.6.x Nodes

After 1.6.2 is validated:

1. **Node 1.6.3 (Uniform integrability of R_ε)** — This is the hardest remaining node. The original plan (L^p bound for fixed p > 1) was shown to FAIL by the 1.5.4 refutation. Need alternative:
   - **Scheffé's lemma approach:** If R_ε → R μ-a.s. and E_μ[R] = 1, then L^1 convergence follows without UI
   - **Relative entropy bounds:** H(T_ψ*μ_ε | μ_ε) bounded → UI
   - **Variational / Girsanov framework:** BG (2021) Girsanov construction may give R_ε as a stochastic exponential with built-in convergence

2. **Node 1.6.4 (Passage to limit)** — Depends on 1.6.3.

### THEN: Easy Remaining Nodes

3. **Nodes 1.7, 1.8** (Strict positivity + symmetry ψ → −ψ) — Straightforward once 1.6 is done.
4. **Parent wrapper nodes** (1.4, 1.5, 1.6, 1.1, 1) — Just summarize validated children.

## Proof Strategy Summary

The proof constructs an explicit Radon-Nikodym derivative d(T_ψ*μ)/dμ via UV regularization and a limiting argument. It has four stages:

### Stage A: Regularized RN Derivative (Nodes 1.2, 1.3) — VALIDATED

R_ε(φ) = exp(Ψ_ε(φ)) where Ψ_ε = −ΔV_ε + ⟨(−Δ+m²)ψ, φ⟩ − ½‖ψ‖²_{H¹}

### Stage B: Wick Expansion (Nodes 1.4.1–1.4.4) — ALL VALIDATED

ΔV_ε decomposes into cubic (converges), quadratic (converges), linear (DIVERGES via δm_ε²), and constant terms. The Wick constant C_ε(0) cancels exactly.

### Stage C: Renormalization (Nodes 1.5.1–1.5.4) — ALL VALIDATED

Ψ_ε = Ψ^{ren}_ε + L_ε + K_ε, where Ψ^{ren} converges in L^p(μ), L_ε diverges linearly, K_ε is a divergent constant absorbed by normalization. Sub-Gaussian tails uniform in ε ensure per-ε finiteness.

**Key finding (from refutation):** For any fixed p > 1, E_{μ_ε}[R_ε^p] diverges as exp(C(p²−1)(log 1/ε)²). Standard L^p UI approach is impossible.

### Stage C': Convergence (Nodes 1.6.1–1.6.4) — IN PROGRESS

- 1.6.1 (Wick power regularity): VALIDATED
- 1.6.2 (Exponential integrability): Third attempt pending verification
- 1.6.3 (Uniform integrability / L^1 convergence): NOT YET ATTEMPTED — needs new approach
- 1.6.4 (Passage to limit): NOT YET ATTEMPTED

### Stage D: Conclusion (Nodes 1.7, 1.8) — NOT YET ATTEMPTED

## Key Pitfalls Discovered

### 1. μ ⊥ μ₀ (Singular measures)
The Φ⁴₃ measure is SINGULAR w.r.t. the GFF. Never assume μ ~ μ₀.

### 2. Do NOT decompose Wick powers into raw powers
Converting :φ^k: = φ^k − C(0)·(lower terms) introduces divergent C_ε(0) ~ 1/ε corrections. Always treat Wick-ordered quantities as single distributional objects in their Besov spaces.

### 3. No Fatou across changing measures
Fatou's lemma requires a FIXED measure. To pass from μ_ε to μ, use Skorokhod coupling on a common probability space.

### 4. The L^p route to UI is blocked
E_{μ_ε}[R_ε^p] diverges for any fixed p > 1. Uniform integrability of {R_ε} requires an alternative approach (Scheffé, relative entropy, or variational/Girsanov).

### 5. BG concentration is for φ only
Barashkov-Gubinelli (2020/2021) sub-Gaussian concentration applies to the field φ and Lipschitz functionals thereof. It does NOT extend to the enhanced data (:φ²:, :φ³:) without additional proof.

### 6. Besov duality requires strict regularity gap
The pairing C^α × C^{−β} → ℝ requires α > β (or s₁ + s₂ > 0 in standard notation), NOT α + β > 0. Always use ‖f‖_{C^{k/2+2δ}} (not k/2+δ) to pair with :φ^k: in C^{−k/2−δ}.

## Challenges Filed in af Workspace

Multiple challenges have been filed by verifiers. Key open challenges:
- `ch-9084a5aa987e9a2a` (minor): Unproven variance lower bound in 1.5.4
- `ch-f5e33a864eca31af` (minor): Incorrect Fourier formula for mass counterterm in 1.4.4
- `ch-3b149fd41c002559` (minor): Incomplete justification for ⟨ψ,φ⟩ ≠ 0 μ-a.s. in 1.4.4

## Definitions and References in the af Workspace

**11 definitions:** `Phi43_measure`, `shift_map`, `pushforward_measure`, `measure_equivalence`, `three_torus`, `Wick_ordering`, `mass_counterterm`, `Radon_Nikodym_derivative`, `sub_Gaussian_tails`, `UV_regularization`, `Cameron_Martin_space`

**9+ external references** (original 6 + additions):
- Hairer (2014) — Regularity structures / Φ⁴₃ construction
- Gubinelli-Imkeller-Perkowski (2015) — Paracontrolled distributions
- Barashkov-Gubinelli (2020) — Variational method for Φ⁴₃ (Duke Math. J.)
- Barashkov-Gubinelli (2021) — Girsanov construction + singularity proof (EJP)
- Barashkov-Gubinelli (2021) — Sub-Gaussian tails
- Bogachev (1998) — Cameron-Martin theorem
- Mourrat-Weber (2017) — Global solutions for Φ⁴₃
- Gubinelli-Hofmanova (2019) — Global solutions Phi^4_3
- Bauerschmidt-Bodineau (2019) — Log-Sobolev for Φ⁴ (lattice)
- Nelson (1973) — Hypercontractivity
- Skorokhod representation / Vitali convergence — standard references

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
