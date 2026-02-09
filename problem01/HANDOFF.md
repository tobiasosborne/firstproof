# Problem 01 Handoff — Φ⁴₃ Measure Equivalence Under Smooth Shifts

## What This Problem Asks

Let μ be the Φ⁴₃ measure on D'(T³) (three-dimensional torus), ψ : T³ → ℝ smooth and nonzero, T_ψ(u) = u + ψ the shift map. **Are μ and T_ψ*μ equivalent (mutually absolutely continuous)?**

## Answer

**YES.** The measures are equivalent.

## Current State (Session 4)

- **Proof tree:** ~82 nodes total (new child 1.6.4.3.3 from repair)
- **~20 VALIDATED** by adversarial verification (including all Stage A, B, C nodes + key Stage C' nodes)
- **~9 REFUTED** (all repaired or archived — 1.6.4.3 refuted this session)
- **~7 ARCHIVED** (superseded by repairs)
- **Remaining:** ~46 pending (many are auto-generated children or exploratory)
- **Run `af status -d .`** from this directory to see the full tree.

### Key Progress This Session (Session 4)

#### 1.6.4.3 Tilted Measure Convergence — REFUTED AND REPAIRED

**Refutation 10: Node 1.6.4.3 (tilted measure convergence μ_ε^L → μ)**

The original claim that the L_eps-tilted measures μ_ε^L converge weakly to μ was **REFUTED** with two critical flaws:

**Critical flaw 1 (Tightness is FALSE):** L_ε = 2δm_ε² ⟨ψ_ε, φ_ε⟩ has coefficient δm_ε² ~ C log(1/ε) → ∞. Under μ_ε^L, the expectation E_{μ_ε^L}[⟨ψ, φ⟩] diverges (at least as (log 1/ε)^{1/3} from competition between quartic confinement λφ⁴ and linear tilt βφ). Since ‖φ‖_{C^{-1/2-δ}} ≥ c|⟨ψ, φ⟩|/‖ψ‖, the Besov norm diverges in expectation, contradicting tightness.

**Critical flaw 2 (Wrong limit identification):** Even if tight, V_ε^L differs from V_ε by a DIVERGENT perturbation. Under χ = φ − ψ, V_ε^L becomes λ:(χ+ψ)⁴: + δm_ε²:χ²:, which expands to the standard Φ⁴₃ potential for χ PLUS cubic/quadratic/linear terms (4λψ:χ³:, 6λψ²:χ²:, etc.) that break Z₂ symmetry and change the limiting measure. Φ⁴₃ uniqueness says the limit is independent of the UV regularization scheme — it does NOT say different potentials produce the same measure.

**4 challenges filed:** ch-fd54aa98b705108d (tightness failure), ch-24f408a6def8934e (wrong limit), ch-5d9cbd5cfd3117eb (one-point insertion unjustified), ch-f8e2b15b70f0b268 (missing dependencies).

**Repair (node 1.6.4.3, 7th amendment):** The tilted-measure-convergence approach was ABANDONED entirely. Node 1.6.4.3 now contains:
1. What IS proved: T_ψ*μ_ε → T_ψ*μ weakly (trivial), ρ_ε → ρ weakly (Vitali), and the covariance decomposition d_ε = Cov_{Q_ε}(G, W_ε)/E_{Q_ε}[W_ε]
2. What NEEDS proving: lim d_ε = 0, deferred to child node 1.6.4.3.3

**New child node 1.6.4.3.3:** Proposes a Boue-Dupuis variational approach:
- The shift T_ψ corresponds to adding a deterministic drift h in BG's stochastic control formulation
- Girsanov's theorem gives the density M_h between the shifted and unshifted optimal drifts
- The identification M_h = exp(Ψ^{ren})/Z (conditional on terminal value) follows from the first variation of the BG functional
- L_ε disappears because it's a regularization artifact — in the continuum BG theory, the drift change h is fixed

**Status of 1.6.4.3.3: PENDING VERIFICATION — this is the #1 priority for Session 5.**

### Validated Nodes (All Sessions — unchanged from Session 3)

| Node | Description | Status |
|------|-------------|--------|
| 1.2 | Setup: μ_ε ~ μ₀, Cameron-Martin, T_ψ*μ_ε ~ μ_ε | VALIDATED |
| 1.3 | Regularized RN derivative formula | VALIDATED |
| 1.4.1 | Quartic Wick shift (algebraic identity) | VALIDATED |
| 1.4.2 | Quadratic Wick shift (algebraic identity) | VALIDATED |
| 1.4.3 | Full interaction difference | VALIDATED |
| 1.4.4 | UV divergence analysis | VALIDATED |
| 1.5.1 | Decomposition of exponent: Ψ_ε = Ψ^ren + L_ε + K_ε | VALIDATED |
| 1.5.2 | Normalization constraint | VALIDATED |
| 1.5.3 | Convergence of Ψ^ren in L^p(μ) | VALIDATED |
| 1.5.4 | Absorption of divergent linear term (scoped) | VALIDATED |
| 1.6.1 | Well-definedness of smeared Wick powers | VALIDATED |
| 1.6.2.5 | k=1 exponential integrability (sub-Gaussian) | VALIDATED |
| 1.6.2.10.2.1 | Partition function ratio | VALIDATED |
| 1.6.2.10.2.2.1 | BG structural requirements | VALIDATED |
| 1.6.2.10.2.3.3 | Uniform Z_eps^{(alpha)} bound | VALIDATED |
| 1.6.2.10.2.4 | Exponential integrability QED (regularized) | VALIDATED |
| 1.6.2.10.3 | Skorokhod passage from mu_eps to mu | VALIDATED |
| 1.6.2.10.4 | Full exponential integrability QED | VALIDATED |
| 1.6.3.4 | UI of exp(Psi^{ren}) on Skorokhod space | VALIDATED |

### Refutation History (HIGH-VALUE FINDINGS)

#### Sessions 1-2 Refutations
- 1.5.4: Incorrect L^p uniform integrability claim → scoped to per-eps
- 1.6.1: Wrong duality conditions → fixed indices
- 1.6.2 attempts 1-2: Wick-to-raw decomposition, BG concentration overreach

#### Session 3 Refutations

**Refutation 5: Node 1.6.2.10.1 (L^p growth rate)**
- Circular reasoning: invoked parent node's conclusion
- Unjustified mu_0-to-mu_eps comparison after Brascamp-Lieb fails
- Node may be vestigial (sibling 1.6.2.10.2 doesn't use its output)

**Refutation 6: Node 1.6.2.10.2 (exponential bounds, attempt 3)**
- A priori estimate ||:phi_eps^k:||_{C^{-k/2-δ}} ≤ C(1+V_eps)^{k/4} is UNJUSTIFIED
- Cited references give dynamical PDE estimates, not equilibrium pathwise bounds
- Norm-type mismatch: Besov supremum vs spatial integral
- Deterministic lower bound V_eps >= -C_0 (uniform in eps) is FALSE

**Refutation 7: Node 1.6.2.10.2.2 (structural requirements)**
- Cubic coupling is RELEVANT (engineering dim +3/2), not irrelevant
- Spatial inhomogeneity of coefficients not addressed
- Both repaired in 1.6.2.10.2.2.1

**Refutation 8: Node 1.6.2.10.2.3 (main Boue-Dupuis step)**
- Adapted drift treated as deterministic (cross-terms don't vanish)
- Cubic again called irrelevant (not propagated from fix)
- Hand-waving at BG mechanism

**Refutation 9: Node 1.6.2.10.2.3.1 (1st repair of main step)**
- Fabricated BG Propositions 4.1 and 4.3 (don't exist as described)
- Fictitious running-coupling Polchinski flow (BG actually use Boue-Dupuis + paracontrolled + Gamma-convergence)
- Repaired in 1.6.2.10.2.3.3

#### Session 4 Refutation

**Refutation 10: Node 1.6.4.3 (tilted measure convergence)**
- Tightness of {μ_ε^L} is FALSE — divergent linear tilt overwhelms quartic confinement
- Limit identification is WRONG — Φ⁴₃ uniqueness ≠ different potentials giving same measure
- The L_ε-tilted measures shift the mean by ~(log 1/ε)^{1/3}, breaking tightness
- Under χ = φ−ψ, the tilted potential generates Z₂-breaking cubic terms → different limit
- Repaired by abandoning tilted-measure approach; new BG variational approach in 1.6.4.3.3

## What the Next Agent Should Do

### IMMEDIATE: Verify node 1.6.4.3.3 (Boue-Dupuis identification) — CRITICAL

This is the repaired core argument. It claims T_ψ*μ = ρμ via the BG variational framework. **Known vulnerabilities to check:**

1. **Part 4 (conditional expectation identity):** E_P[M_h | X_1^{u*} = φ] = exp(Ψ^{ren}(φ))/Z. This is the mathematical heart. Is it justified or just sketched? The stochastic integral ∫⟨h_s, dW_s⟩ depends on the full Brownian path, not just the terminal value X_1. The conditional expectation over paths terminating at φ needs rigorous justification.

2. **Part 5 (Gamma-convergence extension):** The claim that BG's Gamma-convergence (Section 6) absorbs L_ε + K_ε into the limit for the SHIFTED variational problem V(· + ψ) is stated as "standard perturbation theory" without proof. This needs checking — does the BG Gamma-convergence machinery actually extend to shifted potentials?

3. **Part 2 (drift representation):** ψ = ∫₀¹ e^{-(1-s)A} h_s ds for deterministic h. Is the surjectivity claim correct? Is h ∈ L²([0,1]; L²(T³)) sufficient?

4. **Part 3 (Girsanov):** Is the Novikov condition satisfied for M_h? Since h is deterministic with ‖h‖ < ∞, this should be fine, but verify.

5. **Sibling nodes 1.6.4.4 and 1.6.4.5** still reference the OLD approach (Boltzmann ratio with μ_ε^L → μ). If 1.6.4.3.3 is accepted, these may need restructuring since the BG approach provides the identification directly.

### THEN: Verify remaining 1.6.3 nodes
- 1.6.3.2, 1.6.3.3, 1.6.3.8 (straightforward, likely correct)

### THEN: Prove and verify nodes 1.7, 1.8
- **1.7 (Strict positivity):** R = exp(Ψ^{ren})/Z > 0 μ-a.s. (should be easy — exponential is always positive)
- **1.8 (Symmetry ψ → −ψ):** Apply the same argument with −ψ to get μ << T_ψ*μ

### FINALLY: Close parent wrapper nodes
- Nodes 1.4, 1.5, 1.6, 1.1, 1 — summarize validated children

## Proof Strategy Summary

### Stage A: Regularized RN Derivative (Nodes 1.2, 1.3) — VALIDATED
### Stage B: Wick Expansion (Nodes 1.4.1–1.4.4) — ALL VALIDATED
### Stage C: Renormalization (Nodes 1.5.1–1.5.4) — ALL VALIDATED
### Stage C': Convergence (Nodes 1.6.1–1.6.4) — IN PROGRESS
- 1.6.1 (Wick power regularity): VALIDATED
- 1.6.2 (Exponential integrability): EFFECTIVELY RESOLVED (4th repair validated)
- 1.6.3 (UI of exp(Ψ^{ren})): Key node 1.6.3.4 VALIDATED, others pending verification
- 1.6.4 (Passage to limit): Original approach REFUTED; BG variational repair in 1.6.4.3.3 pending verification
### Stage D: Conclusion (Nodes 1.7, 1.8) — NOT YET ATTEMPTED

## Key Pitfalls Discovered (Updated)

### 1–10: (Sessions 1-3)
1. μ ⊥ μ₀ — never assume μ ~ μ₀
2. Do NOT decompose Wick powers into raw powers
3. No Fatou across changing measures — use Skorokhod coupling
4. The L^p route to UI is blocked — E_{μ_ε}[R_ε^p] diverges for any p > 1
5. BG concentration is for φ only — not enhanced data
6. Besov duality requires strict regularity gap — α > β
7. Cubic coupling is RELEVANT in d=3 (engineering dim +3/2)
8. BG (2020) mechanism is NOT a Polchinski running-coupling flow
9. Adapted drifts in Boue-Dupuis are NOT deterministic
10. Do NOT fabricate citations

### 11. Divergent linear tilts BREAK tightness (NEW — Session 4)
The L_ε-tilted measures μ_ε^L := exp(L_ε)dμ_ε / Z have DIVERGENT means. Since L_ε = 2δm_ε²⟨ψ_ε, φ_ε⟩ with δm_ε² ~ log(1/ε) → ∞, the tilt pushes E_{μ_ε^L}[⟨ψ, φ⟩] → ∞ (competition: quartic confinement λφ⁴ vs linear tilt βφ gives mean ~β^{1/3}). This kills tightness in C^{-1/2-δ}. Do NOT assume exponential tilts with divergent coefficients preserve tightness.

### 12. Φ⁴₃ uniqueness ≠ different potentials giving the same measure (NEW — Session 4)
Uniqueness of Φ⁴₃ (Hairer 2014, GH 2019, BG 2020) says the limit is independent of the UV regularization scheme. It does NOT say different interaction potentials produce the same limit. A Φ⁴₃ measure with additional smooth lower-order perturbations (cubic, linear terms) is a DIFFERENT measure. In particular, Z₂-symmetry-breaking perturbations (like cubic terms 4λψ:φ³:) produce a non-Z₂-symmetric limit.

## Definitions and References in the af Workspace

**11 definitions** (unchanged from Session 2)

**References** (expanded):
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
