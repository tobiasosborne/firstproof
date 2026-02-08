# Problem 6: Existence of Large ε-Light Subsets in Graphs

**Field:** Spectral Graph Theory

**Author:** Daniel Spielman (Yale University)

## Setup

For a graph **G = (V, E)**:

- Let **G_S = (V, E(S, S))** denote the graph with the same vertex set, but only the edges between vertices in S (i.e., both endpoints in S).
- Let **L** be the Laplacian matrix of G.
- Let **L_S** be the Laplacian matrix of G_S.

### Definition: ε-light set

A set of vertices S is **ε-light** if the matrix **εL − L_S** is positive semidefinite.

## Question

Does there exist a **constant c > 0** such that for every graph G and every ε ∈ (0, 1):

> V contains an ε-light subset S of size at least **c · ε · |V|**?

## Key Concepts

- **Graph Laplacian L**: For a graph G = (V, E), the matrix L = D − A where D is the diagonal degree matrix and A is the adjacency matrix. L is positive semidefinite.
- **Induced subgraph Laplacian L_S**: The Laplacian of the subgraph induced by S, embedded in ℝ^{|V| × |V|} (with zeros for vertices not in S).
- **ε-light**: S is ε-light means εL − L_S ≽ 0 (positive semidefinite). Intuitively, the edges within S contribute at most an ε fraction of the total Laplacian's quadratic form.
- **Positive semidefinite**: A symmetric matrix M is PSD if x^T M x ≥ 0 for all vectors x.

## Constraints

The proof should be roughly five pages or fewer.
