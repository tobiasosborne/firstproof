# Problem 01 Handoff — Φ⁴₃ Measure Equivalence Under Smooth Shifts

## What This Problem Asks

Let μ be the Φ⁴₃ measure on D'(T³) (three-dimensional torus), ψ : T³ → ℝ smooth and nonzero, T_ψ(u) = u + ψ the shift map. **Are μ and T_ψ*μ equivalent (mutually absolutely continuous)?**

## Answer

**YES.** The measures are equivalent.

## Current State (Session 3)

- **Proof tree:** ~80 nodes total (expanded by 4th repair of exponential integrability + new 1.6.3/1.6.4 content)
- **~20 VALIDATED** by adversarial verification (including all Stage A, B, C nodes + key Stage C' nodes)
- **~8 REFUTED** (all repaired or archived)
- **~7 ARCHIVED** (superseded by repairs)
- **Remaining:** ~45 pending (many are auto-generated children or exploratory)
- **Run `af status -d .`** from this directory to see the full tree.

### Key Progress This Session (Session 3)

#### 1.6.2 Exponential Integrability — EFFECTIVELY RESOLVED
The critical blocker from Session 2 (node 1.6.2, exponential integrability of Wick powers) has been resolved through a 4th repair attempt using the **tilted potential / Boue-Dupuis approach**:

| Node | Description | Status |
|------|-------------|--------|
| 1.6.2.10.2.1 | Partition function ratio reformulation | VALIDATED |
| 1.6.2.10.2.2.1 | Structural requirements for BG (cubic subordination + spatial inhomogeneity) | VALIDATED |
| 1.6.2.10.2.3.3 | Uniform upper bound on Z_eps^{(alpha)} via BG (2020) Theorem 1.1 | VALIDATED |
| 1.6.2.10.2.4 | QED: M(alpha,f,k) = exp(C_+) | VALIDATED |
| 1.6.2.10.3 | Skorokhod passage from mu_eps to mu | VALIDATED |
| 1.6.2.10.4 | QED combining k=1,2,3 | VALIDATED |

**4th repair approach:** Express E_{mu_eps}[exp(alpha <f, :phi_eps^k:>)] = Z_eps^{(alpha)}/Z_eps. The tilted partition function Z_eps^{(alpha)} is bounded uniformly in eps because V_eps^{(alpha)} is a valid Phi^4_3-type potential (same quartic coupling lambda > 0, lower-order perturbation absorbed by Young's inequality since k ≤ 3 < 4). Apply BG (2020) Theorem 1.1 to the tilted potential.

#### 1.6.3 Uniform Integrability — PARTIALLY PROVED
Key results established:
| Node | Description | Status |
|------|-------------|--------|
| 1.6.3.4 | UI of exp(Psi_eps^{ren}) on Skorokhod space | VALIDATED |
| 1.6.3.2 | Skorokhod coupling + a.s. convergence of Psi^{ren} | pending (needs verification) |
| 1.6.3.3 | Z = E_mu[exp(Psi^{ren})] finite and positive | pending (needs verification) |
| 1.6.3.8 | QED: exp(Psi^{ren})/Z is a well-defined density rho | pending (needs verification) |

**Key finding:** The divergent L_eps tilt prevents direct L^1 convergence of R_eps. The identification T_psi*mu = rho is deferred to 1.6.4.

#### 1.6.4 Passage to Limit — PROOF WRITTEN, NEEDS VERIFICATION
A proof has been written using the **Boltzmann ratio + tilted measure convergence** approach:
- Rewrite regularized identity in Boltzmann ratio form (1.6.4.2)
- Show L_eps-tilted measures mu_eps^L converge weakly to mu (1.6.4.3) — **KEY VULNERABILITY**
- Apply Vitali convergence under tilted measures (1.6.4.4)
- Assemble via weak convergence (1.6.4.5)

### Validated Nodes (All Sessions)

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
| 1.6.2.10.3 | Skorokhod passage to mu | VALIDATED |
| 1.6.2.10.4 | Full exponential integrability QED | VALIDATED |
| 1.6.3.4 | UI of exp(Psi^{ren}) on Skorokhod space | VALIDATED |

### Refutation History (HIGH-VALUE FINDINGS)

#### Sessions 1-2 Refutations (see previous HANDOFF for details)
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

## What the Next Agent Should Do

### IMMEDIATE: Verify node 1.6.4.3 (tilted measure convergence) — CRITICAL

This is the most vulnerable new node. It claims mu_eps^L → mu weakly, where mu_eps^L is the L_eps-tilted regularized Phi^4_3 measure. The argument:
1. Complete the square: mu_eps^L is a recentered Phi^4_3 measure
2. Mass counterterm is universal (doesn't depend on smooth centering)
3. Tightness from quartic coercivity
4. Uniqueness of Phi^4_3 limit

**Key vulnerability:** Does Phi^4_3 uniqueness (Hairer 2014 Theorem 1.1) handle recentered potentials? This needs checking.

### THEN: Verify remaining 1.6.3 and 1.6.4 nodes
- 1.6.3.2, 1.6.3.3, 1.6.3.8 (straightforward, likely correct)
- 1.6.4.1, 1.6.4.2, 1.6.4.4, 1.6.4.4.1, 1.6.4.5

### THEN: Prove and verify nodes 1.7, 1.8
- **1.7 (Strict positivity):** R = exp(Psi^{ren})/Z > 0 μ-a.s. (should be easy — exponential is always positive)
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
- 1.6.3 (UI of exp(Psi^{ren})): Key node 1.6.3.4 VALIDATED, others pending verification
- 1.6.4 (Passage to limit): Proof written, all nodes pending verification
### Stage D: Conclusion (Nodes 1.7, 1.8) — NOT YET ATTEMPTED

## Key Pitfalls Discovered (Updated)

### 1–6: (Same as Session 2 — see below)
1. μ ⊥ μ₀ — never assume μ ~ μ₀
2. Do NOT decompose Wick powers into raw powers
3. No Fatou across changing measures — use Skorokhod coupling
4. The L^p route to UI is blocked — E_{μ_ε}[R_ε^p] diverges for any p > 1
5. BG concentration is for φ only — not enhanced data
6. Besov duality requires strict regularity gap — α > β

### 7. Cubic coupling is RELEVANT in d=3 (NEW)
Engineering dimension of ∫φ³dx is d − 3(d−2)/2 = 3/2 > 0. Do NOT call it "irrelevant." Use quartic subordination via Young's inequality instead.

### 8. BG (2020) mechanism is NOT a Polchinski running-coupling flow (NEW)
BG use: Boue-Dupuis variational formula + paracontrolled a priori estimates (Section 3/8) + Gamma-convergence (Section 6). They do NOT track running couplings through a scale-by-scale Polchinski flow. The Polchinski equation appears only as the HJB equation for the value function (a duality, not the proof mechanism).

### 9. Adapted drifts in Boue-Dupuis are NOT deterministic (NEW)
In the Boue-Dupuis formula, u is an adapted process, so I_u = ∫₀¹ u_s ds is RANDOM. Cross-terms E[:X^j: h^{n-j}] do NOT vanish for j ≥ 1. The BG approach handles adapted drifts through a priori estimates, not by naively taking expectations.

### 10. Do NOT fabricate citations (NEW)
BG (2020) Section 4 is "Three dimensions" not "Gamma-convergence." There is no "Proposition 4.1" or "Proposition 4.3" tracking running couplings. Actual theorem numbers: Theorem 1 (main), Theorem 2 (Boue-Dupuis). Always verify citations.

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
