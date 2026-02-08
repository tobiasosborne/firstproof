# Archive Manifest — Problem 04: Fisher Superadditivity

**Date:** 2026-02-08
**Canonical proof tree:** `../` (problem04 root, 24-node af tree)

---

## Archive Structure

### `problem04_original/`
The bare af proof tree that existed in problem04 before consolidation.
- 1 node, 2 externals, 6 ledger entries

### `examples6/` — Subordination Approach
**Origin:** `~/Projects/af-tests/examples6/`
**Strategy:** Finite subordination functions ω₁, ω₂ → chain rule → vector inequality
**434 files** including:

| Path | Description |
|------|-------------|
| `HANDOFF.md` | Session 129-130 status, approaches killed, directions |
| `fisher_subordination_proof.md` | Source mathematics (320 lines): definitions, subordination, chain rule, Hard Lemma 2 |
| `fisher_proof/` | Original af proof tree (44 nodes, 61% complete) |
| `fisher_proof/ledger/` | 44 ledger entries |
| `strategy3_proof/` | Modified Strategy 3 attack (25+ nodes, 6 validated) |
| `strategy3_proof/prove_n3_symbolic.py` | **n=3 COMPLETE PROOF** |
| `strategy3_proof/prove_n4_symmetric_proof.py` | **n=4 symmetric PROOF** |
| `strategy3_proof/prove_n4_coefficient.py` | n=4 formula derivation |
| `strategy3_proof/verification_n3_proof.md` | Independent audit: 56/56 checks pass |
| `strategy3_proof/audit_reformulation_chain.md` | Full chain audit |
| `strategy3_proof/findings_*.md` | 9 findings documents (approaches tried) |
| `strategy3_proof/investigate_*.py` | Exploration scripts |
| `strategy3_proof/verify_*.py` | Verification scripts (correct MSS formula) |

**Key results:** n=3 proved, n=4 symmetric proved, ⟨h,α⟩≥0 DISPROVED, partition of unity corrected, boxplus formula conflict resolved

### `examples7/` — Cumulant / Heat Flow Approach
**Origin:** `~/Projects/af-tests/examples7/`
**Strategy:** Cumulant decomposition 1/Φ = C_n·κ₂ + R_n → prove R_n superadditive; then heat flow / de Bruijn / EPI
**238 files** including:

| Path | Description |
|------|-------------|
| `HANDOFF.md` | Sessions 130-132 status, 25 agents deployed |
| `fisher_subordination_proof2.md` | Conjecture statement + MSS definitions |
| `R4_proof_result.md` | κ₃=0 analytical proof (PROVER-11) |
| `Rn_general_report.md` | C_n formula + R_n pattern (PROVER-10b) |
| `prover14_root_geometry_report.md` | Φ_n=2·Sm2 identity (MAJOR) |
| `prover15_literature_report.md` | Literature survey: Voiculescu, Shlyakhtenko-Tao |
| `verifier11_report.md` | κ₃=0 proof VALIDATED |
| `verifier12_report.md` | Heat flow audit — monotone gap KILLED |
| `verifier13_report.md` | Root geometry audit |
| `verifier14_report.md` | Monotone gap DEFINITIVELY REFUTED |
| `*.py` (63 files) | ~28K lines verification/proof scripts |
| `R4_prover11_final.py` | κ₃=0 verification (all tests pass) |
| `Rn_prover10b.py` | Comprehensive R_n computation |

**Key results:** n=3 proved (independent), n=4 κ₃=0 proved, C_n=4/(n²(n-1)) proved, Φ_n=2·Sm2 proved, S2 additivity proved, de Bruijn validated, Gaussian splitting proved, EPI numerically validated (13K+ trials), monotone gap KILLED

### `examples8_canonical/`
The canonical hybrid proof tree (snapshot before copy to problem04 root).
- 24 nodes, 14 externals, strategy report
- This is a backup; the live version is in the problem04 root

---

## Reconstruction Guide

To reconstruct any previous proof state:
- `examples6/fisher_proof/`: run `af status` in that directory
- `examples6/strategy3_proof/`: run `af status` in that directory
- `examples7/`: run `af status` in that directory (155 ledger entries)
- `problem04_original/`: the bare tree before consolidation

All Python scripts are self-contained and can be run with `python3 <script>.py` (require numpy, sympy).

## Exhausted Approaches (DO NOT RETRY)

| Approach | Source | Finding |
|----------|--------|---------|
| ⟨h,α⟩ ≥ 0 | ex6 | FALSE (~0.3-1% violations n≥3) |
| Monotone gap along heat flow | ex7 | FALSE (44% violation rate) |
| Shape factor SF(r) ≤ min | ex6 | FALSE (42.7% violations) |
| Partition of unity ω₁'+ω₂'=1 | ex6 | FALSE (correct: each ω'=1) |
| Joint concavity of -R₄ | ex7 | Hessian indefinite |
| SOS on gap numerator | ex7 | Mixed-sign cross terms |
| Coefficient additivity n≥4 | ex6 | Cross terms in g_k for k≥4 |
| AM-GM from A+B≥2Φ_r | ex6 | Wrong direction |
| Polynomial manipulation for n≥4 general | both | Exhausted by 4+ prover agents |
