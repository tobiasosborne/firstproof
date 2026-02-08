# HANDOFF — Fisher Subordination Proof (examples6)

## What This Is

An adversarial proof formalization of the **finite free Fisher information superadditivity conjecture**:

> For monic real-rooted polynomials p, q of degree n with simple roots,
> 1/Phi_n(p boxplus_n q) >= 1/Phi_n(p) + 1/Phi_n(q)

Built with the `af` CLI tool (Adversarial Proof Framework). There are **two** af workspaces:

1. **`fisher_proof/`** — Original proof tree (44 nodes, 61% complete). Established definitions, subordination, chain rule. Hard Lemma 2 remains open.
2. **`strategy3_proof/`** — Modified Strategy 3 attack (25+ nodes, 6 validated, 5 challenges resolved). Extensively tested in Sessions 127-129.

The source conjecture document is `fisher_subordination_proof.md`.

---

## Session 130 Results (Latest)

### MAJOR RESULT: n=3 PROVED ALGEBRAICALLY

The Fisher superadditivity conjecture is **proved for n=3**. Independently verified.

**Proof sketch:** For centered cubics (e_1=0), MSS boxplus is additive in coefficients: E_r=E_p+E_q, F_r=F_p+F_q where E=-e_2>0, F=e_3. Formula: 1/Phi_3=(4E^3-27F^2)/(18E^2). The excess is a positive definite quadratic form in (Fp,Fq) with det=2Ep^3*Eq^3*(Ep+Eq)^2>0. Equality iff both polynomials have equally-spaced roots.

**Files:** `prove_n3_symbolic.py` (proof), `verify_n3_proof.py` (independent verification), `verification_n3_proof.md` (audit report, 56/56 checks pass).

### n=4: Partial Progress

- **1/Phi_4 = -disc(f)/(4·I·J)** where I=e2^2+12e4>0, J=2e2^3-8e2e4+9e3^2<0
- **Symmetric subcase (e3=0) PROVED** via strict concavity of phi(t)=t(1-4t)/(1+12t)
- General n=4 BLOCKED: MSS cross term (1/6)e2p·e2q in g_4 breaks coefficient additivity; 659-term excess
- 500K random trials: zero violations
- **Files:** `prove_n4_coefficient.py`, `prove_n4_symmetric_proof.py`, `findings_n4_coefficient.md`

### Approaches Investigated and Killed (Session 130)

| Approach | Status | Why |
|---|---|---|
| Shape factor SF(r) <= min(SF(p),SF(q)) | **KILLED** | FALSE: 42.7% violations at n=3 |
| Coefficient additivity for n>=4 | **BLOCKED** | Cross terms in g_k for k>=4 |
| Schur convexity of Phi in gaps | **BLOCKED for n>=4** | Phi NOT Schur-convex in gaps for n>=4 |

### Structural Facts Discovered (Session 130)

1. **Var(r) = Var(p) + Var(q) exactly** (proved for all n via e_2 additivity of centered boxplus)
2. **Centered root majorization**: roots of r ALWAYS majorize those of p and q (tested n=2-8)
3. **SF(r) <= max(SF(p), SF(q))** appears true (proved for n=3, 60K+ trials n=3-8)
4. **Reformulation chain audited**: 56/56 checks pass, 5 minor gaps (none fatal)

### Most Promising Remaining Directions (Updated)

1. **n=4 general case via SOS/DSOS**: The 659-term excess might be certifiable as sum-of-squares via semidefinite programming
2. **Free probability / random matrix approach**: Since r = E_U[chi_{A+UBU*}], use trace inequality methods with the CORRECT formula
3. **Induction via derivative identity with simultaneous-level argument**: F_k >= 0 at ALL derivative levels, suggesting a proof across all levels at once
4. **Direct subordination-based argument**: Use omega_1, omega_2 properties more deeply

---

## Quick Start for Next Orchestrator

```bash
cd examples6/strategy3_proof
af status          # 25+ node proof tree
af challenges      # open challenges
af jobs            # available work
```

**READ THIS ENTIRE HANDOFF BEFORE SPAWNING AGENTS.**

---

## Session 129 Results (Current)

### CRITICAL DISCOVERY 1: `<h,alpha> >= 0` is FALSE

The sub-lemma `<h,alpha> >= 0` (node 1.7.2) was **disproved** by independent verification:
- 34 counterexamples in 10,000 trials at n=4 (~0.34% failure rate)
- Confirmed by corrector agent (~1% at n=3,4,5,6) and verifier agent independently
- Concrete counterexample (n=4): roots_p=[-2.957,-2.657,-0.975,0.015], roots_q=[-3.461,1.388,1.688,2.668] gives <h,alpha>=-0.192
- TRUE only at n=2 (proved analytically)
- **Node 1.7.2 marked as FALSE in the proof tree**

### CRITICAL DISCOVERY 2: Boxplus Formula Conflict Resolved

Two different boxplus formulas exist in the literature:
- **CORRECT (MSS):** `c_k = sum_{i+j=k} [(n-i)!(n-j)! / (n!(n-k)!)] * a_i * b_j` — verified against Monte Carlo (100K+ Haar unitaries)
- **EQUIVALENT:** `g_k = sum_{i+j=k} C(n-j,i)/C(n,i) * e_i(p) * e_j(q)` — same formula, different notation
- **WRONG:** `hat_e_k(r) = sum_j hat_e_j(p) * hat_e_{k-j}(q)` where `hat_e_k = e_k/C(n,k)` — differs for non-centered polynomials

An agent claimed counterexamples at n=4 using the WRONG formula. Definitive Monte Carlo verification (100K Haar unitaries, p=q={-5,-1,1,5}) showed the correct formula gives NO violation. **The conjecture is NOT refuted.** Challenge raised on node 1.10.3.

### Challenges Resolved (Session 129)
| Challenge | Node | Resolution |
|---|---|---|
| ch-640cc7d38a4 (critical) | 1.6.1 | Sigma corrected to order-preserving. <h,alpha> mostly >= 0 but not always. |
| ch-6eba2806382 (minor) | 1.5 | Positivity node 1.5.2 added. |
| ch-94f57d6b50c (critical) | 1.7.1 | <h,alpha> >= 0 is FALSE. Node amended. |
| ch-9007eb39a0c (major) | 1.7.3 | "CRITICAL DISCOVERY" → "OBSERVATION". Pythagorean only at n=2,3. |

### Proof Approaches Attempted and Their Status

| Approach | Node | Status | Why |
|---|---|---|---|
| <h,alpha> >= 0 sub-lemma | 1.7.2 | **DEAD** | Disproved by counterexamples |
| Herglotz convexity | 1.7.2 | **DEAD** | omega_1 is algebraic, not rational Nevanlinna |
| Contour integral / residues | 1.7.2 | **DEAD** | Proves <h,delta>>=0 but not <h,alpha>>=0 |
| Induction via r'=n*(p^{(1)} boxplus q^{(1)}) | 1.10.1 | **STUCK** | Correction term has indeterminate sign (~40% negative) |
| Direct algebraic identity for AB-h^4 | 1.10.2 | **STUCK** | No manifestly non-negative expression found |
| AM-GM via A+B >= 2*Phi_r | 1.10.2 | **DEAD** | AM-GM goes wrong direction: (A+B)/2 >= sqrt(AB) is UPPER bound |
| Random matrix trace inequality | 1.10.3 | **INVALID** | Used wrong boxplus formula; "refutation" withdrawn |

### Useful Mathematical Facts Discovered

1. **Cauchy matrix decomposition (VERIFIED):**
   `<h,alpha> = sum_{k<l} (h_k - h_l)(phi_k - phi_l) / [(nu_k - nu_l)(lambda_k - lambda_l)]`
   where phi_k = nu_k - lambda_k.

2. **Divided difference identity (VERIFIED):**
   `<h,alpha> = sum_{k<l} H_r[k,l] * (1/D_{kl} - 1)`
   where H_r[k,l] is the divided difference of H_r, D_{kl} = (lambda_k-lambda_l)/(nu_k-nu_l).

3. **Position displacement positivity (PROVED via contour integral):**
   `<h, delta> >= 0` where delta_k = nu_k - lambda_k (NOT the same as <h,alpha>).

4. **Generating function identity (VERIFIED):**
   `sum_k H_r(nu_k)/(nu_k - p) = -(1/2)*r''(p)/r(p)` for p not a root of r.

5. **A + B >= 2*Phi_r (NUMERICALLY TRUE, proved for n=2,3):**
   Equivalently Phi_p + Phi_q >= 4*Phi_r. But does NOT imply AB >= Phi_r^2 via AM-GM (wrong direction).

6. **Induction decomposition (EXACT):**
   `F_n = F_{n-1} + (delta_r - delta_p - delta_q)` where delta_p = 1/Phi_n(p) - 1/Phi_{n-1}(p^{(1)}).

---

## The Mathematical Situation (Updated)

### What Is Proved (validated, no gaps)

1. **Partition of unity is FALSE.** omega_1'(nu_k) = omega_2'(nu_k) = 1, sum = 2 not 1.
2. **Vector reformulation.** u = h + alpha, v = h + beta from chain rule.
3. **Clean reformulation.** Target: AB >= ||h||^4 where A = Phi_p - Phi_r, B = Phi_q - Phi_r.
4. **n=2 equality.** Proved algebraically. AB = ||h||^4 exactly.
5. **Numerical confirmation.** 0 violations in 5000+ trials with CORRECT MSS formula (verified by Monte Carlo).
6. **<h,alpha> >= 0 is FALSE.** Counterexamples at n >= 3 (~0.3-1%).
7. **A > 0 and B > 0 always.** Fisher info decreases under convolution.

### What Is Open (the hard part)

8. **Prove AB >= ||h||^4.** No proof exists. All attempted approaches have failed or are stuck.

---

## Most Promising Remaining Directions

### Direction 1: Schur convexity / majorization (NEW)

The Herglotz prover discovered that for random (non-MSS) sorted point pairs, <h,alpha> >= 0 fails 71% of the time. This means the MSS structure is **essential**. The MSS convolution has the specific property that centered r-roots majorize centered p-roots. Phi is Schur-concave. Can this structure be exploited more directly?

### Direction 2: n=3 direct computation (TRACTABLE)

For n=3, the conjecture might be provable by direct algebraic computation. The Gram matrix of (h, alpha, beta) is singular at n=3 (3 vectors in R^3 with MSS constraint force linear dependence). This constraint might close the inequality.

### Direction 3: Two-level simultaneous argument

The induction decomposition shows F_k >= 0 at ALL derivative levels simultaneously. Rather than proving each level from the next, look for an argument that works across ALL levels at once — perhaps via a generating function or spectral identity.

### Direction 4: Herglotz representation of omega_1

omega_1(z) = z + c_0 + sum c_j/(d_j - z) with c_j of specific signs. The positions d_j are between the roots of r. Even though individual c_j may be negative, the COMBINED effect forces AB >= ||h||^4. This is structurally similar to a positive-definite kernel argument.

### Direction 5: Prove 1/Phi_r >= 1/Phi_p + 1/Phi_q directly

Skip the AB >= ||h||^4 reformulation entirely. Work directly with:
  sum 1/(sum_{j!=i} 1/(nu_i - nu_j))^2 >= sum 1/(sum_{j!=i} 1/(lambda_i - lambda_j))^2 + sum 1/(sum_{j!=i} 1/(mu_i - mu_j))^2
This is unwieldy but might yield to Schur-convexity or convexity of 1/Phi arguments.

---

## What NOT to Do (Updated Traps)

1. **Do NOT assume <h,alpha> >= 0.** It is FALSE.
2. **Do NOT use AM-GM to go from A+B >= 2c to AB >= c^2.** AM-GM goes the wrong way.
3. **Do NOT use the hat-e convolution formula `hat_e_k(r) = sum hat_e_j(p)*hat_e_{k-j}(q)`.** It is WRONG for non-centered polynomials. Use the MSS formula or verify against Monte Carlo.
4. **Do NOT assume partition of unity.** omega_1' + omega_2' = 2, NOT 1.
5. **Do NOT expect naive Cauchy-Schwarz to close.** Gives GM bound, need HM bound.
6. **Do NOT try induction without addressing the indeterminate-sign correction term.**

---

## Orchestration Notes

- **Subagent rule:** Every job must be a fresh subagent. Max 5 parallel.
- **Boxplus formula:** ALWAYS use `c_k = sum_{i+j=k} [(n-i)!(n-j)! / (n!(n-k)!)] * a_i * b_j` or equivalently `g_k = sum_{i+j=k} C(n-j,i)/C(n,i) * e_i(p) * e_j(q)`. VERIFY against Monte Carlo if in doubt.
- **Verifier standard:** Accept ONLY airtight steps. Challenge anything with gaps.
- **Prover standard:** Be honest about gaps. A sorry with notes beats a false claim.

---

## Files (Updated)

```
examples6/strategy3_proof/
  # Core
  meta.json, ledger/

  # Verification scripts (CORRECT formula)
  verify_sigma_corrected.py        # Order-preserving sigma verification
  verify_boxplus_definitive.py     # DEFINITIVE: all 3 formulas vs Monte Carlo
  verify_boxplus_highprec.py       # 1M-sample Monte Carlo for edge case
  verify_claims_171_173.py         # Verifier: nodes 1.7.1, 1.7.3
  verify_claims_171_173_final.py   # Final version with correct MSS formula

  # Investigation scripts
  investigate_induction_v4.py      # Induction with correct boxplus (FINAL)
  investigate_algebraic_identity.py # Direct AB - h^4 identity search
  investigate_matrix_approach.py   # Matrix trace approach (USES WRONG FORMULA)
  investigate_AplusB.py            # A+B >= 2*Phi_r investigation
  investigate_herglotz*.py         # Herglotz convexity attempts

  # Proof writeups
  proof_h_alpha_herglotz.md        # Herglotz approach (BLOCKED)
  proof_h_alpha_residues.md        # Residue approach (useful identities found)
  verification_wave1_proofs.md     # Verifier report on Wave 1 proofs
  findings_induction.md            # Induction findings (correct formula)
  findings_algebraic.md            # Algebraic identity findings
  findings_matrix_info.md          # Matrix approach (WRONG FORMULA - see challenge)
  findings_AplusB.md               # A+B >= 2*Phi_r findings
```
