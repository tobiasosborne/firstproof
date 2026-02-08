# Problem 9: Algebraic Relations on Quadrilinear Determinantal Tensors

**Field:** Tensor Analysis / Algebraic Geometry

**Author:** Joe Kileel (University of Texas at Austin)

## Setup

Let **n ≥ 5**.

Let **A⁽¹⁾, …, A⁽ⁿ⁾ ∈ ℝ^{3×4}** be Zariski-generic matrices.

For α, β, γ, δ ∈ [n], construct the tensor **Q^{(αβγδ)} ∈ ℝ^{3×3×3×3}** so that its (i, j, k, ℓ) entry for 1 ≤ i, j, k, ℓ ≤ 3 is:

> Q^{(αβγδ)}_{ijkℓ} = det [A⁽α⁾(i,:); A⁽β⁾(j,:); A⁽γ⁾(k,:); A⁽δ⁾(ℓ,:)]

where A(i,:) denotes the i-th row of matrix A, and semicolon denotes vertical concatenation (forming a 4×4 matrix whose determinant is taken).

## Question

Does there exist a **polynomial map F : ℝ^{81n⁴} → ℝ^N** satisfying all three of the following properties?

### Property 1: Independence from A
The map **F** does not depend on A⁽¹⁾, …, A⁽ⁿ⁾.

### Property 2: Bounded degree
The degrees of the coordinate functions of **F** do not depend on n.

### Property 3: Characterization of scaling
Let λ ∈ ℝ^{n×n×n×n} satisfy λ_{αβγδ} ≠ 0 precisely for α, β, γ, δ ∈ [n] that are not all identical. Then:

> F(λ_{αβγδ} · Q^{(αβγδ)} : α, β, γ, δ ∈ [n]) = 0

holds **if and only if** there exist u, v, w, x ∈ (ℝ*)ⁿ such that:

> λ_{αβγδ} = u_α · v_β · w_γ · x_δ

for all α, β, γ, δ ∈ [n] that are not all identical.

## Key Concepts

- **Zariski-generic**: The matrices lie outside a proper algebraic subvariety; i.e., they satisfy no special polynomial relations. Equivalently, the result holds for "generic" matrices.
- **Determinantal tensor Q^{(αβγδ)}**: A 3×3×3×3 tensor whose entries are 4×4 determinants formed by stacking rows from four of the given matrices.
- **Polynomial map F**: A map whose coordinate functions are polynomials in the entries of the input tensors.
- **Rank-1 scaling**: The condition λ_{αβγδ} = u_α v_β w_γ x_δ means that the 4-way array λ has tensor rank 1 (restricted to non-identical indices).
- **(ℝ*)ⁿ**: Vectors in ℝⁿ with all nonzero entries.

## Constraints

The proof should be roughly five pages or fewer.
