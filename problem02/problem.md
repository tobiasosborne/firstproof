# Problem 2: Existence of Whittaker Functions for Rankin–Selberg Integrals

**Field:** Representation Theory / Automorphic Forms

**Author:** Paul D. Nelson (Aarhus University)

## Setup

Let **F** be a non-archimedean local field with ring of integers **𝔬**.

Let **N_r** denote the subgroup of GL_r(F) consisting of upper-triangular unipotent elements.

Let **ψ : F → ℂ×** be a nontrivial additive character of conductor 𝔬, identified in the standard way with a generic character of N_r.

Let **Π** be a generic irreducible admissible representation of GL_{n+1}(F), realized in its ψ⁻¹-Whittaker model W(Π, ψ⁻¹).

## Question

Must there exist **W ∈ W(Π, ψ⁻¹)** with the following property?

### The Property

Let **π** be a generic irreducible admissible representation of GL_n(F), realized in its ψ-Whittaker model W(π, ψ).

Let **𝔮** denote the conductor ideal of π, let **Q ∈ F×** be a generator of 𝔮⁻¹, and set:

> u_Q := I_{n+1} + Q · E_{n,n+1} ∈ GL_{n+1}(F)

where E_{i,j} is the matrix with a 1 in the (i,j)-entry and 0 elsewhere.

For some **V ∈ W(π, ψ)**, the **local Rankin–Selberg integral**:

> ∫_{N_n \ GL_n(F)} W(diag(g, 1) · u_Q) · V(g) · |det g|^{s - 1/2} dg

is **finite and nonzero** for all s ∈ ℂ.

## Key Concepts

- **Non-archimedean local field**: A locally compact totally disconnected field (e.g., ℚ_p or 𝔽_q((t))).
- **Generic irreducible admissible representation**: An irreducible smooth representation of GL_r(F) that admits a Whittaker model.
- **Whittaker model W(Π, ψ⁻¹)**: The space of functions on GL_{n+1}(F) obtained by the ψ⁻¹-Whittaker functional applied to Π.
- **Conductor ideal 𝔮**: The ideal measuring the ramification of π.
- **Local Rankin–Selberg integral**: A zeta integral pairing Whittaker functions of GL_{n+1} and GL_n, central to the theory of automorphic L-functions.
- **E_{i,j}**: The elementary matrix with 1 in position (i,j) and 0 elsewhere.

## Constraints

The proof should be roughly five pages or fewer.
