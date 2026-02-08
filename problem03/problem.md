# Problem 3: Markov Chain with ASEP Polynomial Stationary Distribution

**Field:** Algebraic Combinatorics

**Author:** Lauren Williams (Harvard University)

## Setup

Let **λ = (λ₁ > ⋯ > λₙ ≥ 0)** be a partition with distinct parts.

Assume moreover that λ is **restricted**, meaning:
- It has a unique part of size 0
- No part of size 1

Let **S_n(λ)** denote the set of compositions associated with λ.

For μ ∈ S_n(λ), define:

- **F*_μ(x₁, …, xₙ; q, t)**: the interpolation ASEP polynomial
- **P*_λ(x₁, …, xₙ; q, t)**: the interpolation Macdonald polynomial

## Question

Does there exist a **nontrivial** Markov chain on S_n(λ) whose stationary distribution is given by:

> F*_μ(x₁, …, xₙ; q=1, t) / P*_λ(x₁, …, xₙ; q=1, t)   for μ ∈ S_n(λ)

If so, prove that the Markov chain you construct has the desired stationary distribution.

By **"nontrivial"** we mean that the transition probabilities of the Markov chain should **not** be described using the polynomials F*_μ(x₁, …, xₙ; q, t).

## Key Concepts

- **Partition with distinct parts**: A sequence λ₁ > λ₂ > ⋯ > λₙ ≥ 0 of non-negative integers.
- **Restricted partition**: A partition with a unique part of size 0 and no part of size 1.
- **S_n(λ)**: The set of compositions that are permutations of the parts of λ.
- **Interpolation ASEP polynomial F*_μ**: A family of polynomials arising from the asymmetric simple exclusion process (ASEP), parametrized by compositions μ.
- **Interpolation Macdonald polynomial P*_λ**: A specialization of the Macdonald polynomial theory related to interpolation.
- **Markov chain**: A stochastic process on S_n(λ) with transition matrix whose rows sum to 1.
- **Stationary distribution**: A probability distribution π such that πP = π where P is the transition matrix.

## Constraints

The proof should be roughly five pages or fewer.
