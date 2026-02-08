# HANDOFF: Fisher Superadditivity Proof Tree

**Session date:** 2026-02-08 (Sessions 130-132)
**Tool:** `af` (Adversarial Proof Framework) in `examples7/`
**Agents used:** 25 total (12 provers, 13 verifiers) across 6 waves

---

## 1. The Conjecture

For all monic real-rooted polynomials `p, q` of degree `n` with simple roots:

```
1/Phi_n(p ⊞_n q) >= 1/Phi_n(p) + 1/Phi_n(q)
```

where `⊞_n` is the MSS finite free additive convolution and
`Phi_n(p) = sum_i H_p(lambda_i)^2` is the finite free Fisher information.

**Source document:** `examples7/fisher_subordination_proof2.md`

---

## 2. Breakthrough: Cumulant Decomposition

```
1/Phi_n = C_n * kappa_2 + R_n(kappa_2, ..., kappa_n)
```

- **C_n = 4/(n^2*(n-1))** [VALIDATED by verifier-10, verifier-10b]
  - Original claim C_n = 2/(n(n-1)) was WRONG (off by factor 2/n)
  - C_2=1, C_3=2/9, C_4=1/12, C_5=1/25, C_6=1/45
- **kappa_k** = finite free cumulants (additive under MSS convolution)
- **R_n** = nonlinear correction (always <= 0, weight-2 homogeneous)
- Conjecture REDUCES TO: R_n is superadditive

---

## 3. Current Proof State

### VALIDATED NODES (adversarially verified)
- **Node 1.5.2.5**: n=3 proof COMPLETE (Cauchy-Schwarz on k3^2/k2^2)
- **Node 1.5.2.6.1**: R_4 k3=0 case PROVED analytically [VALIDATED by verifier-11]
  - Proof: Key Lemma (sA+tB <= (s+t)C) + Cauchy-Schwarz
  - 15 independent checks, ~1.5M numerical trials, 0 violations
- **Node 1.5.5.1**: General C_n formula + R_n pattern [VALIDATED by verifier-10b]
  - C_n = 4/(n^2(n-1)) for n=2,...,10
  - R_n superadditivity: 0 violations for n=3,4,5,6
  - R_4 exact rational formula verified

### OPEN (main remaining challenge)
- **Full R_4 analytical proof (k3 nonzero)**: NOT CLOSED
  - Joint concavity FAILS (Hessian indefinite)
  - Term-by-term partial fraction subadditivity FAILS
  - PROVER-12 tried decomposition into Term_A + Term_B — neither individually subadditive
  - Numerically verified: 0/447K violations
  - Promising: SOS certificate, perturbation from k3=0, matrix Cauchy-Schwarz

### NUMERICAL ONLY
- R_5, R_6 superadditivity: 0 violations (no analytical proof)
- General R_n pattern: numerically confirmed for n=2,...,10

---

## 4. Agent Summary (this session pair 130-131)

### Wave 1 (from crashed orchestrator, completed)
- PROVER-8: Rewrote node 1.5.2, resolved 4 challenges
- PROVER-9: R_4 structural analysis (node 1.5.2.6)
- PROVER-10: General R_n pattern, Hermite proof of C_n
- VERIFIER-8: n=3 proof correct
- VERIFIER-9: R_4 formula confirmed, C_n error found

### Wave 2 (Session 130)
- **PROVER-10b** (a9e9996): C_n for n=2,...,10, R_n superadditivity → node 1.5.5.1
- **PROVER-11** (abecf90): k3=0 PROVED, full case numerically verified → node 1.5.2.6.1
- **VERIFIER-10** (aae363f): C_n formula audit, node 1.5.2 amended

### Wave 3 (Session 131)
- **VERIFIER-11** (aa10331): k3=0 proof VALID → node 1.5.2.6.1 VALIDATED
- **VERIFIER-10b** (adb15a3): C_n + R_n claims VALID → node 1.5.5.1 VALIDATED
- **PROVER-12** (a07f251): Full R_4 attempt — did NOT close proof, found Term_A/B decomp fails

---

## 5. Key Files

### Reports
- `R4_proof_result.md` — PROVER-11's k3=0 proof writeup
- `R4_proof_attempt.md` — PROVER-9's structural analysis
- `R4_verification.md` — VERIFIER-9's R_4 audit
- `Rn_general_report.md` — PROVER-10b's general pattern report
- `n3_proof_verification.md` — VERIFIER-8's n=3 audit
- `verifier11_report.md` — VERIFIER-11's k3=0 verification
- `verifier10b_report.md` — VERIFIER-10b's C_n verification

### Verification Scripts
- `R4_prover11_final.py` — Consolidated R_4 verification (all tests pass)
- `Rn_prover10b.py` — Comprehensive R_n computation
- `verifier11_check.py` — 15-check adversarial k3=0 verification
- `verifier10b_check.py` — Independent C_n/R_n verification

---

## 6. Orchestrator Assessment: The Polynomial Approach Has Hit a Wall

**Node requiring radical revision: 1.5.2.6** (R_4 superadditivity)

The current approach — decompose R_n into polynomials, then prove non-negativity
of the gap polynomial via algebraic manipulation — has been exhaustively tried
by 4 prover agents and **every standard technique has failed**:

| Technique | Result | Why it fails |
|-----------|--------|--------------|
| Joint concavity of -R_4 | FAILED | Hessian is indefinite |
| SOS on gap numerator | FAILED | Mixed-sign cross terms |
| Perspective function | BLOCKED | Denominator not linear in k4 |
| Titu/Engel Cauchy-Schwarz | INAPPLICABLE | Denominator nonlinear (degree 10) |
| Partial fraction subadditivity | FAILED | -K4/(24K2) term not subadditive |
| Term_A + Term_B decomposition | FAILED | Neither piece individually subadditive |

This pattern — where the inequality is numerically rock-solid (0/2M+ violations)
but resists all polynomial-level proofs — strongly suggests **the inequality is
true for structural reasons that pure algebraic manipulation cannot access.**

### The analogy with classical information theory

The classical analog, Stam's inequality (1/I(X+Y) >= 1/I(X) + 1/I(Y)), is NOT
proved by expanding in cumulants and doing polynomial arithmetic. It's proved via:
- **de Bruijn's identity** (Fisher info = derivative of entropy along heat flow)
- **Data processing inequality** (Fisher info decreases under channels)
- **Optimal transport** arguments

The finite free version likely needs a similarly structural argument.

### Directions for novel prover agents

The next orchestrator should spawn provers that explore **non-polynomial** approaches:

1. **Finite free heat flow / Ornstein-Uhlenbeck**: Is there a discrete analog of
   de Bruijn's identity for Phi_n? The MSS convolution with a Gaussian polynomial
   should monotonically decrease Phi_n — if this can be established, it may yield
   the inequality via a semigroup argument.

2. **Monotonicity under natural operations**: Does Phi_n satisfy a data processing
   inequality? If there's a "channel" operation on polynomials that contracts Phi_n,
   and MSS convolution factors through it, the inequality follows.

3. **Representation theory of S_n**: MSS convolution comes from expected
   characteristic polynomials of random matrices. The cumulants kappa_k have
   algebraic meaning in terms of S_n representations. Perhaps R_n superadditivity
   follows from a representation-theoretic identity.

4. **Log-convexity / Schur-convexity of root gaps**: The inequality may follow
   from properties of the root interlacing that MSS convolution guarantees,
   rather than from the cumulant expansion.

5. **Direct n=4 via different parametrization**: Instead of (k2,k3,k4), use
   the roots directly. Phi_4 has a clean expression in terms of root gaps
   delta_ij = lambda_i - lambda_j. The MSS convolution has known root interlacing
   properties. Perhaps the inequality follows from these geometric constraints.

6. **McCrimmon-style operator identity**: In Jordan algebra theory, many
   identities that resist direct polynomial proof follow from operator identities
   (U, T, L operators). Is there an operator framework for Phi_n?

### What to preserve from current work

The cumulant decomposition framework IS correct and useful:
- C_n = 4/(n^2(n-1)) is validated
- The reduction to R_n superadditivity is clean
- The k3=0 case and n=3 case are proved and can serve as base cases
- The numerical infrastructure is solid for testing new ideas

The framework should be KEPT but the proof strategy for the core inequality
needs to come from outside the "polynomial gap non-negativity" paradigm.

---

## 7. DO NOT RETRY (exhaustively failed approaches)

- Joint concavity of -R_4 or f(u,v)
- Titu/Engel or standard Cauchy-Schwarz on the fraction
- Perspective function / convex conjugate
- Term-by-term partial fraction subadditivity
- Direct SOS decomposition of gap numerator in (s,t,u,v,m,w) variables
- Any approach that starts with "expand the gap polynomial and show it's non-negative"

---

## 8. Wave 4 Results: Novel Non-Polynomial Approaches (Session 132)

Three provers explored structural approaches beyond polynomial manipulation:

### PROVER-13 (aed99d8): Heat Flow / De Bruijn
- **De Bruijn identity VALID**: dS/dt = Phi_n(p) where S = sum_{i<j} log|r_i - r_j|
- **Root dynamics**: dr_i/dt = H_p(r_i) (Dyson Brownian motion drift)
- **EPI analog VALID**: N(p ⊞ q) >= N(p) + N(q) with N = exp(2S/m), 0 violations / 13,770 trials
- **n=2 EPI proved analytically**: N(p ⊞ q) = N(p) + N(q) with EQUALITY
- **Gaussian splitting VALID**: (p ⊞ G_s) ⊞ (q ⊞ G_t) = (p ⊞ q) ⊞ G_{s+t} (exact)
- **1/Phi concavity**: 1/Phi(p ⊞ G_t) concave in t (0 violations / 446 trials)
- **Monotone gap PROPOSED but INVALID** (see Wave 5 below)

### PROVER-14 (a2d8116): Root Geometry / Electrostatics
- **PROVED**: Phi_n = 2 * Sm2 where Sm2 = sum_{i<j} 1/(lambda_i - lambda_j)^2
  - Proof via triple identity: 1/((a-b)(a-c)) + cyclic = 0
- **S2 additivity PROVED**: S2(p ⊞ q) = S2(p) + S2(q) (exact algebraic identity)
- **Universal identity PROVED**: sum_i H_i * lambda_i = n(n-1)/2
- **Q = Sm2*S2 scale-invariant**: Q_r <= max(Q_p, Q_q) (0 violations / 50K trials)
- **n=3 trigonometric parametrization**: G(c) = 9/sin^2(3*arccos(c)), convex

### PROVER-15 (a27cea6): Literature Survey
- Conjecture = exact finite-n analog of Voiculescu's free Stam inequality (1998)
- Most promising framework: finite subordination + L^2 projection
- Key refs: Shlyakhtenko-Tao 2020, Gribinski 2019, Marcus-Spielman-Srivastava 2022

---

## 9. Wave 5 Results: Adversarial Verification (Session 132)

### VERIFIER-12 (ac01689): Audit PROVER-13 — HIGH VALUE FINDINGS
| Claim | Verdict |
|-------|---------|
| De Bruijn identity | **VALID** (algebraic backbone confirmed) |
| EPI N(p⊞q)>=N(p)+N(q) | **VALID** (0/13,770 violations; n=2 proved analytically) |
| Monotone gap | **INVALID** — gap NOT monotonically decreasing! |
| Gaussian splitting | **VALID** (exact identity; but cumulant formula is WRONG) |

**Critical finding**: The monotone gap proof mechanism is BROKEN.
- gap_H(t) = same-noise gap: monotone decreasing but GOES NEGATIVE
- gap_K(t) = double-noise gap: non-negative but NOT monotone (15/137 pairs increase)
- PROVER-13's cumulant formula kappa_k = (-1)^k*a_k*n^{k-1}/C(n,k) is WRONG
- Despite this, the CONJECTURE appears TRUE (gap always non-negative)

### VERIFIER-13 (a302546): Audit PROVER-14
| Claim | Verdict |
|-------|---------|
| Phi_n = 2*Sm2 | **VALID** (triple identity + combinatorial counting) |
| S2 additivity | **VALID** (exact algebraic proof from MSS formula) |
| sum H_i*lambda_i = n(n-1)/2 | **VALID** (elementary pairing argument) |
| Q_r <= max(Q_p,Q_q) | **PARTIALLY VALID** (numerically supported but NOT proved; does NOT imply conjecture) |
| n=3 trigonometric param | **MOSTLY VALID** (G(c) formula correct; convexity confirmed; n=3 remains unproved) |

### VERIFIER-14 (a20154f): Audit monotone gap proof strategy — STILL RUNNING

---

## 10. Validated Structural Results (new this session)

These are now PROVED (algebraic identities, not just numerical):
1. **Phi_n = 2*Sm2** (via triple identity cancellation)
2. **S2(p ⊞ q) = S2(p) + S2(q)** (from MSS coefficient formula)
3. **sum_i H_i*lambda_i = n(n-1)/2** (pairing argument)
4. **Gaussian splitting exact**: (p ⊞ G_s) ⊞ (q ⊞ G_t) = (p ⊞ q) ⊞ G_{s+t}
5. **n=2 equality**: 1/Phi(p⊞q) = 1/Phi(p) + 1/Phi(q) for ALL degree-2 polynomials

---

## 11. Revised Orchestrator Assessment (Session 132)

The monotone gap proof strategy is DEAD (killed by VERIFIER-12). However, Wave 4 discovered extraordinary structural results. The most promising next directions:

1. **EPI-to-Stam reduction**: The EPI N(p⊞q) >= N(p)+N(q) is numerically validated to 13K+ trials. If proved, does it IMPLY Stam? In classical theory, EPI implies Stam. Check if this holds here.

2. **De Bruijn + concavity**: 1/Phi(p⊞G_t) is concave in t. Combined with de Bruijn identity dS/dt = Phi_n, this gives d²S/dt² <= 0 (entropy concavity). This is the finite analog of Costa's strengthening.

3. **S2 additivity + Phi=2*Sm2**: Since S2 is additive, the conjecture 1/Phi_r >= 1/Phi_p + 1/Phi_q becomes S2_r/(2*Q_r) >= S2_p/(2*Q_p) + S2_q/(2*Q_q), i.e., (S2_p+S2_q)/Q_r >= S2_p/Q_p + S2_q/Q_q, i.e., Q_r <= harmonic_mean(Q_p, Q_q) weighted by S2.

4. **n=3 complete proof**: With trigonometric parametrization + MSS contraction in cos(3phi), n=3 is tantalizingly close. The contraction factor (a_p^{3/2}+a_q^{3/2})/(a_p+a_q)^{3/2} < 1 pushes Q_r toward minimum.

---

## 12. Previous Sessions

### Session 132 (this): Waves 4-5. Novel approaches (heat flow, root geometry, literature). Monotone gap KILLED. EPI+structural results validated.
### Session 131: Wave 3 — 2 verifiers validated, 1 prover didn't close. Orchestrator assessment added.
### Session 130: Wave 2 — k3=0 proved, C_n verified, R_n pattern confirmed
### Session 129: Crashed orchestrator recovery, Wave 1 results collected
