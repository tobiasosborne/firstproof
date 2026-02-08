# Problem 10: Preconditioned CG for CP Decomposition with RKHS Constraints

**Field:** Numerical Linear Algebra / Tensor Decomposition

**Authors:** Tamara G. Kolda (MathSci.ai) and Rachel Ward (UT Austin)

## Setup

Given a **d-way tensor T ∈ ℝ^{n₁ × n₂ × ⋯ × n_d}** with missing entries (unaligned data), we consider computing a **CP decomposition of rank r** where some modes are infinite-dimensional and constrained to lie in a **Reproducing Kernel Hilbert Space (RKHS)**.

We use an **alternating optimization** approach, and focus on the **mode-k subproblem** for an infinite-dimensional mode.

### Notation

- **N = ∏ᵢ nᵢ**: product of all mode sizes
- **n ≡ n_k**: size of mode k
- **M = ∏_{i≠k} nᵢ**: product of all dimensions except mode k
- **n ≪ M** (mode k is small relative to the other modes)
- **q ≪ N**: number of observed entries (sparse observations)
- **T ∈ ℝ^{n × M}**: mode-k unfolding of the tensor with missing entries set to zero
- **S ∈ ℝ^{N × q}**: selection matrix (subset of columns of I_N) such that S^T vec(T) selects the q known entries
- **Z = A_d ⊙ ⋯ ⊙ A_{k+1} ⊙ A_{k-1} ⊙ ⋯ ⊙ A₁ ∈ ℝ^{M × r}**: Khatri-Rao product of all factor matrices except mode k
- **B = TZ ∈ ℝ^{n × r}**: the MTTKRP (Matricized Tensor Times Khatri-Rao Product)
- **K ∈ ℝ^{n × n}**: positive semidefinite RKHS kernel matrix for mode k
- **A_k = KW** where **W ∈ ℝ^{n × r}** is the unknown

### The Linear System

The system to be solved for W is:

> [(Z ⊗ K)^T S S^T (Z ⊗ K) + λ(I_r ⊗ K)] vec(W) = (I_r ⊗ K) vec(B)

This is a system of size **nr × nr**. A standard direct solver costs O(n³r³), and explicitly forming the matrix is an additional expense.

## Question

Explain how an **iterative preconditioned conjugate gradient (PCG)** linear solver can be used to solve this problem more efficiently. Specifically:

1. **Describe the method** and choice of preconditioner.
2. **Explain in detail** how the matrix-vector products are computed and why this works.
3. **Provide complexity analysis**.

### Assumptions
- n, r < q ≪ N
- Avoid any computation of order N

## Key Concepts

- **CP decomposition**: Expressing a tensor as a sum of r rank-1 tensors: T ≈ Σ_{j=1}^{r} a₁ⱼ ⊗ a₂ⱼ ⊗ ⋯ ⊗ a_dj.
- **RKHS (Reproducing Kernel Hilbert Space)**: A Hilbert space of functions where evaluation functionals are continuous. The kernel matrix K encodes inner products.
- **Khatri-Rao product (⊙)**: Column-wise Kronecker product of matrices.
- **Kronecker product (⊗)**: The standard matrix Kronecker product.
- **MTTKRP**: Matricized Tensor Times Khatri-Rao Product, a key operation in tensor decomposition.
- **Conjugate gradient (CG)**: An iterative method for solving symmetric positive definite linear systems.
- **Preconditioner**: A matrix that approximates the inverse of the system matrix to accelerate CG convergence.
- **Selection matrix S**: Encodes the observation pattern for missing data.

## Constraints

The proof/explanation should be roughly five pages or fewer. All computations must avoid O(N) cost.
