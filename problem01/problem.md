# Problem 1: Equivalence of the Φ⁴₃ Measure Under Smooth Shifts

**Field:** Stochastic Analysis / Quantum Field Theory

**Author:** Martin Hairer (EPFL and Imperial)

## Setup

Let **T³** be the three-dimensional unit-size torus.

Let **μ** be the **Φ⁴₃ measure** on the space of distributions D'(T³).

Let **ψ : T³ → ℝ** be a smooth function that is **not identically zero**.

Define the **shift map** T_ψ : D'(T³) → D'(T³) by:

> T_ψ(u) = u + ψ

where we use the usual identification of smooth functions as distributions.

Let **T_ψ* μ** denote the **pushforward** of μ under T_ψ.

## Question

Are the measures **μ** and **T_ψ* μ** equivalent?

Here, **equivalence of measures** means they have the same null sets (i.e., they are mutually absolutely continuous).

## Key Concepts

- **Φ⁴₃ measure**: The probability measure on D'(T³) associated with the Φ⁴₃ Euclidean quantum field theory in three dimensions. This is a singular measure constructed via renormalization techniques.
- **Pushforward measure**: (T_ψ* μ)(A) = μ(T_ψ⁻¹(A)) for measurable sets A.
- **Mutual absolute continuity**: μ ≪ T_ψ* μ and T_ψ* μ ≪ μ (equivalently, they have the same null sets).

## Constraints

The proof should be roughly five pages or fewer.
