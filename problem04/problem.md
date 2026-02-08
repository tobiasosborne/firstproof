# Problem 4: Superadditivity of Inverse Φ Under Free Additive Convolution

**Field:** Free Probability / Spectral Theory

**Author:** Nikhil Srivastava (UC Berkeley)

## Setup

Let p(x) and q(x) be two **monic polynomials of degree n**:

> p(x) = Σ_{k=0}^{n} a_k · x^{n-k}    where a₀ = 1

> q(x) = Σ_{k=0}^{n} b_k · x^{n-k}    where b₀ = 1

### Finite Free Additive Convolution

Define **p ⊞ₙ q(x)** to be the polynomial:

> (p ⊞ₙ q)(x) = Σ_{k=0}^{n} c_k · x^{n-k}

where the coefficients c_k are given by:

> c_k = Σ_{i+j=k} [(n-i)! · (n-j)!] / [n! · (n-k)!] · aᵢ · bⱼ

for k = 0, 1, …, n.

### The Functional Φₙ

For a monic polynomial p(x) = ∏_{i≤n} (x - λᵢ), define:

> Φₙ(p) := Σ_{i≤n} (Σ_{j≠i} 1/(λᵢ - λⱼ))²

and Φₙ(p) := ∞ if p has a multiple root.

## Question

Is it true that if p(x) and q(x) are **monic real-rooted polynomials** of degree n, then:

> 1/Φₙ(p ⊞ₙ q) ≥ 1/Φₙ(p) + 1/Φₙ(q)?

## Key Concepts

- **Monic polynomial**: A polynomial with leading coefficient 1.
- **Real-rooted polynomial**: A polynomial all of whose roots are real.
- **Finite free additive convolution (⊞ₙ)**: An operation on degree-n monic polynomials that is the finite analog of free additive convolution in free probability. It preserves real-rootedness.
- **Φₙ(p)**: A functional measuring the "spread" of roots via the sum of squared Cauchy transforms evaluated at the roots. Related to the logarithmic energy of the roots.
- **Superadditivity of 1/Φₙ**: The claimed inequality says that 1/Φₙ is superadditive under ⊞ₙ, analogous to how variance is additive under free convolution.

## Constraints

The proof should be roughly five pages or fewer.
