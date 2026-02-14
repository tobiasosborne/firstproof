# Official Solution: Problem 7: Lattices with 2-Torsion and Rationally Acyclic Universal Covers (Shmuel Weinberger)

> Source: *First Proof Solutions and Comments*, Abouzaid, Blumberg, Hairer, Kileel, Kolda, Nelson, Spielman, Srivastava, Ward, Weinberger, Williams. February 14, 2026.

---

## Authors' Commentary on AI-Generated Solutions


In the no internet version, Theorem 4 and in the internet version it is Lemma 5, are false (they are the same statement). The counterexample is  $\mathbb{R}^1$  and f is a translation. It has no fixed points, but its Lefschetz number in their sense is -1.

All proofs by AI's I've seen only use finite complex and Poincaré duality. However, Fowler's paper shows that if  $\Gamma$  is a lattice in a linear semisimple group G, then taking a homomorphism from  $\Gamma$  to a finite group  $\Delta$ , with kernel  $\Gamma_0$  torsion free, the product  $M^3 \times (K \setminus G/\Gamma_0 \times E\Delta)/\Delta$ , where  $E\Delta$  is a contractible space with free  $\Delta$  action, and  $M^3$  is any closed hyperbolic 3-manifold, has the rational type of a finite complex, and satisfies Rational Poincaré duality. It has fundamental group  $\pi_1(M^3) \times \Gamma$  which is a lattice in  $SO(3,1) \times G$ . This shows that all such proofs must fail.

Some proofs try to use "multiplicativity of Euler characteristic in finite covers". This is false for infinite complexes with finitely generated homology over  $\mathbb{Q}$ . The simplest example I know is the following: Consider the universal cover of  $\mathbb{R}P^2$  wedge an infinite number of  $S^2$ 's. It has an involution, and  $\pi_2$  is  $\mathbb{Z}[-1] + \mathbb{Z}[\mathbb{Z}/2]^{\infty}$ . ( $\mathbb{Z}[-1]$  is  $\mathbb{Z}$  acted on by the involution by multiplication by -1.) This module is, after tensoring with  $\mathbb{Z}[1/2]$  a free  $\mathbb{Z}[1/2][\mathbb{Z}/2]$  module, so one can use a free basis to equivariantly attach  $D^3 \times \mathbb{Z}/2$ 's to kill the homology (=homotopy). The new space will be rationally acyclic, and both it and its quotient under  $\mathbb{Z}/2$  will be, and will have rational Euler characteristic = 1.

---

## Official Solution

## Fowler's theorem for involutions.

## Sylvain Cappell, S.Weinberger, and M. Yan

Fowler, in his Ph.D. thesis, proved that if Γ is a uniform lattice in a real semisimple group with odd torsion in Γ then there is no compact closed manifold M whose universal cover is rationally acyclic. A proof can be found in [W2]. We show that the same is true for Γ with 2 torsion.

Without loss of generality (by considering a normal subgroup of finite index), it suLices to prove this for the special case where Γ = π ⋊ Z2 for a torsion free group π, a lattice in G, for which there is an involution on M = K\G/π (by isometries with the locally symmetric metric) whose fixed set F is not empty. (F might be disconnected; for simplicity we will write what follows just for the connected case – there are no diLerences in the general case.)

Now suppose that Xm is a manifold with fundamental group Γ, Y its 2-fold cover, and suppose that the universal cover of X (and therefore Y) are rationally acyclic. We will consider the symmetric signatures of Y in the (symmetric = quadratic L-group) L(**R**π), where **R** is the real numbers. There is an equivalence f:Y -> M which (while not degree one) gives an equivalence of symmetric signatures (because over **R**, all degrees have square roots, so the symmetric signature is only sensitive to the sign of the degree of the map). Since the Novikov conjecture is true for π, the assembly map from Hm(Bπ; L(**R**)) -> Lm(**R**π) is injective, and this detects in the degree m piece Hm(Bπ; **Z**) the class that these manifolds represent in group homology. It follows that this map is degree one. f\*[Y] = [M].

Now we use a cobordism argument from [W1]. We now consider the image of the fundamental class of any manifold Z with fundamental group π involution inducing this automorphism of π and the image of [Z] in Hm(BΓ; **Ζ**<sup>2</sup> ). It follows from standard equivariant homotopy theory that Z has an equivariant map, g, to M, and thus there is a map from its fixed set ZZ2 -> F. We claim that g\*[Z] = g\* [ZZ2] where we make use of the map from Z2 x π1F -> Γ (and the periodicity on the group homology of Z2 to raise the dimension from that of F to dim M).

This cobordism is between Z and a projective space bundle over ZZ2 - namely the projectivized normal bundle to ZZ2. (The fundamental class of the latter is the desired element by the Leray-Hirsch theorem.) It is explicitly Z x [0,1] and on Z x {1} mod out in the complement of the equivariant regular neighborhood of ZZ2 the Z/2 action.

Thus for Y, this image is 0, since the action is free. For M however, this is always nonzero. The action by Z2 by isometries has fixed set which is aspherical and indeed the Borel

construction for the action on M shows that Z2xF -> Γ induces an injection on homology in dimension dim(M/Z2) (and an isomorphism in higher dimensions, see [B]). Since the fundamental class of an aspherical manifold is always nontrivial in its group homology, we have a contradiction.

## References

[B] A.Borel, A seminar on transformation groups, Princeton University Press 1960

[W1] S.Weinberger, Group actions and higher signatures II, CPAM 1987

[W2] S.Weinberger, Variations on a theorem of Borel, Cambridge University Press 2022
