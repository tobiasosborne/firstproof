# Official Solution: Problem 8: Lagrangian Smoothing of Polyhedral Lagrangian Surfaces (Mohammed Abouzaid)

> Source: *First Proof Solutions and Comments*, Abouzaid, Blumberg, Hairer, Kileel, Kolda, Nelson, Spielman, Srivastava, Ward, Weinberger, Williams. February 14, 2026.

---

## Authors' Commentary on AI-Generated Solutions


The best two solutions produced during testing both correctly identified the existence of a local smoothing near every vertex; the proof uses essentially the same basic linear algebra argument that appears in the human solution. The proof then proceeds to perform a local-to-global gluing argument. It was a priori clear that there must be a gap in this argument because the LLM solution refers to the existence of a linear symplectic transformation that brings a neighbourhood of each vertex and each edge into a standard position, but fails to discuss the compatibility between these choices. In the case of the solution produced by the model which was not discouraged to use the internet, the error was finally identified, after a careful reading, in Step 3 of the Proof of Theorem 1: the LLM system asserted that one can choose disjoint neighbourhoods of the edges and of the vertices. In the other case, the error is in Step 2: the model performs a local move near vertices, which changes the local geometry near the edges, invalidating the application of the edge move.

The errors in these solutions can be repaired at the cost of significant computations of changes of coordinates, which would become extremely burdensome in any generalisation. The point of the solution we provide is to obtain a proof which avoids (most of) the hard work, and which experts can readily generalise to other symplectic manifolds (in any dimension).

---

## Official Solution

## Mohammed Abouzaid

## February 3, 2026

Remark 1. This note is expanded from a short motivating discussion in a research paper that is supposed to develop a theory of polyhedral Lagrangian submanifolds for the purpose of being able to use computers to explore conjectures in symplectic topology. It includes some details that would normally be omitted (e.g. the proof of Lemma 1, which is a linear algebra exercise, and much of the explanation about closed 1-forms). The paper does not cite any references as the reader is assumed to be able to deduce all asserted results from standard references, e.g. [1, 2].

I would like to thank Kyler Siegel and Umut Varolgunes for helpful discussions around this circle of ideas.

For the purpose of this note, we equip R <sup>4</sup> with coordinates (q1, q2, p1, p2), and with the standard symplectic form ω = dp<sup>1</sup> ∧ dq<sup>1</sup> + dp<sup>2</sup> ∧ dq2.

Definition 1. A polyhedral Lagrangian surface in R 4 is a finite polyhedral complex all of whose faces are Lagrangians, and which is a topological submanifold of R 4 .

Proposition 1. If K is a polyhedral Lagrangian surface with the property that exactly 4 faces meet at every vertex, then there is a Hamiltonian isotopy K<sup>t</sup> of smooth Lagrangian submanifolds, parameterised by (0, 1], extending to a topological isotopy, parametrised by [0, 1], with endpoint K<sup>0</sup> = K.

In order to prove this result, we need two preliminary results: a local statement asserting triviality near each vertex, and a global statement implying the compatibility of these local trivialisations.

Lemma 1. For each embedding R <sup>2</sup> → R <sup>4</sup> which is linear on the four quadrants with Lagrangian image, and whose image Σ is not contained in a plane, there is a linear symplectic transformation of R <sup>4</sup> which maps Σ to the product of the union of the positive coordinate axes in R 2 p1q<sup>1</sup> and R 2 p2q<sup>2</sup> .

Proof. Let (v1, v2, u1, u2) denote tangent vectors at the origin to the edges of Σ, ordered so that cyclically adjacent vectors span the faces of Σ. The pairings ω(v<sup>i</sup> , ui) cannot vanish, for otherwise ω would identically vanish on a 3-dimensional linear subspace. By swapping the pair of coordinates (v<sup>i</sup> , ui) if necessary, we may assume that both pairings are strictly positive, and by rescaling we may assume that they are 1. We conclude that the vectors (v1, v2, u1, u2) form a standard symplectic basis for R 4 , and that the mapping ∂p<sup>i</sup> → v<sup>i</sup> and ∂q<sup>i</sup> → u<sup>i</sup> is the desired linear transformation.

In the plane R 2 pq, the symplectic pairing projects the union of the positive axes homeomorphically to the dual of the line p = q. Taking the product, and applying the previous Lemma, we conclude:

Corollary 1. There exists a linear Lagrangian plane L ⊂ R 4 so that the symplectic pairing R <sup>4</sup> → L <sup>∨</sup> defines a homeomorphism Σ → L ∨.

The previous corollary in particular equips Σ with a smooth structure arising from its projection to L <sup>∨</sup>. This smooth structure will be fixed for the remainder of the discussion.

Given a choice of plane L, we say that a Lagrangian Λ ⊂ R 4 is graphical if the symplectic pairing defines a diffeomorphism Λ ∼= L <sup>∨</sup>. If Σ were smooth, the standard description of Lagrangians in cotangent bundles would imply that such Lagrangians bijectively correspond to smooth closed 1 forms, which, because Σ is contractible and hence every closed form on it is exact, can be identified with smooth functions modulo addition of constants. We shall formulate a replacement for this correspondence that accounts for the singularities of Σ.

To this end, let us choose further a Lagrangian splitting of the projection R <sup>4</sup> → L <sup>∨</sup>; we shall later see that our constructions are independent of this choice. The splitting gives a direct sum decomposition R <sup>4</sup> ∼= L ⊕ L <sup>∨</sup> (polarization), with respect to which the image of each quadrant is graphical over L <sup>∨</sup>. Graphical (linear) Lagrangians bijectively correspond to quadratic forms, so we obtain quadratic forms {qij}i,j∈± on L <sup>∨</sup> whose graphs contain the corresponding faces of Σ. The restriction of the quadratic forms associated to any two faces agree to first order along the images in L <sup>∨</sup> of the edges of Σ. Via the identification Σ ∼= L <sup>∨</sup> from the previous corollary, we write q<sup>Σ</sup> for the C 1 -function on Σ whose restriction to each face is given by the composition of qij with the projection to L <sup>∨</sup>. We use this to obtain an explicit description of the desired local smoothings, which will be essential in establishing the required global smoothability:

Definition 2. The space S (Σ) of smoothing functions for Σ is the space of C 1 functions f : Σ → R satisfying the property that the function on f + q<sup>Σ</sup> is infinitely differentiable.

It follows immediately from the definition that S (Σ) is invariant under addition of smooth functions, which will be used in the next result:

Lemma 2. The space of smoothing functions S (Σ) depends only on L (and not on the splitting of the projection R <sup>4</sup> → L <sup>∨</sup>).

Proof. A different choice of complementary subspaces correspond to adding a quadratic form q ′ to qij , and the corresponding smooth function on Σ to qΣ.

We shall now associate a graphical Lagrangian to each smoothing function: the construction relies on the fact that the union of all translates of L passing through a face of Σ is canonically symplectomorphic to the cotangent bundle of Σ, with the cotangent fibre at z ∈ Σ corresponding to the translate of L passing through z. In this way, a smoothing function f determines a Lagrangian Λdf ⊂ R 4 , piecewise as the graph of the restriction of the differential df to each face.

Lemma 3. The assignment f 7→ Λdf determines a bijective correspondence between graphical Lagrangians and smoothing functions on Σ up to addition of constants.

Proof. In terms of the polarization from the discussion preceding Definition 2, the Lagrangian Λdf corresponds to the graph of the differential of the function f + q<sup>Σ</sup> considered as a function on L ∨ via the projection map, because each face of Σ is the graph of dqij . The result now follows from the fact that graphical Lagrangians over L <sup>∨</sup> are graphs of differentials of smooth functions.

Note that while the proof uses the polarization, the construction does not. As in Lemma 2, we conclude that this bijection depends only on the choice of Lagrangian L.

The above completes our local analysis near vertices. Near edges, the analysis is much simpler:

Lemma 4. If Σ consists of a pair of linear Lagrangian half-planes in R <sup>4</sup> meeting along a line ℓ, then the space of Lagrangian subspaces L, satisfying the property that the symplectic pairing Σ → L ∨ is a homeomorphism, is contractible.

Proof. The submanifold Σ is equivalent by (affine) linear symplectic transformations to the symplectic product of the real axis in an R 2 factor with the piecewise Lagrangian consisting of the positive axes in another. If the projection Σ → L <sup>∨</sup> is a homeomorphism, then L must be transverse to both Lagrangian half-planes comprising Σ. This implies that the symplectic reduction of L along ℓ (i.e. the image under the quotient by ℓ of the intersection of L with the symplectic annihilator ℓ ⊥) is a line transverse to two coordinate lines in ℓ <sup>⊥</sup>/ℓ ∼= R 2 , and Σ projects homeomorphically to L ∨ if and only if this reduction intersects the interior of the positive quadrant, which is a contractible condition. The argument is completed by noting that the space of Lagrangian lifts of a line ℓ ′ in R 2 is contractible: any two lifts to ℓ <sup>⊥</sup> differ by the graph of a map from ℓ ′ to ℓ, and L is determined up to contractible choice by L ∩ ℓ <sup>⊥</sup>, since it must lie in the symplectic orthogonal of this line, and the space of planes in R 3 containing a given line (in this case L ∩ ℓ <sup>⊥</sup>) and avoiding another line (in this case ℓ) is contractible.

Extending Definitions 2 and 3 verbatim to the case of a pair of edges, we obtain the analogue of Lemma 3, using a splitting into factors as in the above proof.

In the global setting, we cannot work with translates with a single Lagrangian, so we need to consider a family L<sup>z</sup> of Lagrangian planes, passing through each point z ∈ Σ, which are not necessarily translates of each other. We shall require four properties of such a family, the first three of which are easy to state:

- 1. L<sup>z</sup> consists of translates of a single Lagrangian near the origin.
- 2. L<sup>z</sup> varies smoothly along the edges.
- 3. L<sup>z</sup> varies smoothly along the faces.

To formulate the last property, say that σ and σ ′ are faces meeting along an edge τ , and let z be a point on τ . The choice of L<sup>z</sup> determines an identification

$$T_z \sigma \cong L_z^{\vee} \cong T_z \sigma'$$

which is compatible with the inclusion of Tzτ on both sides. A matched normal field along τ is a choice of sections of T σ|<sup>τ</sup> and T σ′ |<sup>τ</sup> which are inward pointing, and are opposite vectors under the above identification. For simplicity, we require this normal field, at the origin τ , to point along the direction of the edge of σ (or σ ′ ) which meets τ . Because the faces of Σ are flat, this choice therefore determines an embedding τ × [0, ϵ) → σ, which is a collar neighbourhood (and similarly for σ ′ ).

Definition 3. A conormal fibration dual to Σ is a family L<sup>z</sup> of (affine)-linear Lagrangian planes in R 4 , parametrised by z ∈ Σ, satisfying the above three properties and so that, in a collar of each edge, the Lagrangians in the normal direction are translates of the Lagrangians along the edge.

The choice of collars in the above construction determines a smooth structure on Σ by using negative coordinates on one of the collars as well as the identification (−ϵ, 0]∪[0, ϵ) ∼= (−ϵ, ϵ). This is an a priori different way of constructing a smooth structure than our earlier formulation, and the next result asserts the compatibility of these contructions; in this setting, we choose an affine-linear Lagrangian Λ<sup>z</sup> passing through z, which is transverse to Lz, and consider the (locally defined) map from Σ to Λ<sup>z</sup> which assigns to z ′ ∈ Σ near z the intersection points L<sup>z</sup> ′ ∩Λ<sup>z</sup> which is unique because L<sup>z</sup> is close to L<sup>z</sup> ′.

Lemma 5. The projection map to Λ<sup>z</sup> is a local diffeomorphism.

*Proof.* The only case that needs to be discussed is when z lies on an edge  $\tau$ . The condition that  $L_{z'}$  be given by translates along the collar direction implies that this map may be written along the collar of  $\tau$  in a face  $\sigma$  as  $(t,s) \mapsto \gamma(t) + s \cdot \nu_{\sigma}(t)$ , where t is the coordinate along  $\tau$  and  $s \in [0,\epsilon)$  is the coordinate in the normal direction. The requirement that the normal fields are matched is equivalent to the condition that  $\nu_{\sigma} = -\nu_{\sigma'}$  if  $\sigma$  and  $\sigma'$  are the two faces meeting along  $\tau$ . The smoothness of the map is immediate from this description.

Whenever the family  $L_z$  does not consist of translates, the Lagrangians  $L_z$  will have non-empty intersections. However, such intersections always take place outside some open neighbourhood  $\nu\Sigma$  of  $\Sigma$ , which we now fix. As before, the fibration  $L_z$  determines a projection  $\nu\Sigma \to \Sigma$ . We say that a Lagrangian is *graphical* with respect to  $L_z$  if it is contained in this neighbourhood, and its projection to  $\Sigma$  is a diffeomorphism.

**Lemma 6.** Every graphical Lagrangian with respect to  $L_z$  arises as the graph of a smoothing function. Moreover, any smoothing function whose differential is sufficiently small defines a graphical Lagrangian.

*Proof.* The correspondence between graphical Lagrangians and smoothing functions is local on  $\Sigma$ . It thus suffices to consider a point  $z \in \Sigma$ , and observe that a Lagrangian plane  $L_z^{\vee}$  which is transverse to  $L_z$  at z will also be transverse to nearby fibres, so that a neighbourhood of z in  $\nu\Sigma$  is modelled after the conormal bundle of  $L_z^{\vee}$ , by Weinstein's tubular neighbourhood theorem. The result then follows by the standard construction of Lagrangians as graphs of closed 1-forms.

In order for the previous result to be helpful, we need to be able to produce the desired functions; this is not completely obvious because the space of smoothing functions is not invariant under rescaling:

**Lemma 7.** There exist smoothing functions of arbitrarily small  $C^1$ -norm.

*Proof.* As a preliminary step, choose a partition of unity  $\sum_{\sigma} \chi_{\sigma} = 1$  on  $\Sigma$ , of bounded  $C^k$ -norms for all k, indexed by the strata of  $\Sigma$ , so that  $\chi_{\sigma}$  vanishes outside a small neighbourhood of  $\sigma$  and its restriction to  $\sigma$  is identically 1 in the complement of a small neighbourhood of the boundary of  $\sigma$ . If  $\chi_{\sigma}^{\epsilon}$  is the composition of  $\chi_{\sigma}$  with the dilation of the plane by  $1/\epsilon$ , we obtain a family of partitions of unity which are uniformly bounded, and whose  $C^1$ -norms are bounded by a constant multiple of  $1/\epsilon$ .

We now choose a Lagrangian plane  $\Lambda_{\sigma}$  which contains each stratum  $\sigma \subset \Sigma$ , and which is transverse to L, and let  $f_{\sigma}$  denote the corresponding smoothing function. Note that the tangency conditions imply that the functions  $f_{\sigma}$  and  $f_{\sigma'}$  agree to first order along  $\sigma \cap \sigma'$ . Let  $f^{\epsilon}$  denote the function  $\sum \chi_{\sigma}^{\epsilon} f_{\sigma}$ . The fact that  $f_{\sigma}^{\epsilon}$  is a family of smoothing functions follows from the partition of unity, and the fact that the  $C^1$ -norm is bounded follows from the product rule and the observation that, while the norm of the gradient of  $\chi_{\sigma}^{\epsilon}$  grows like  $1/\epsilon$ , it is supported in a region where the difference between  $f_{\sigma}$  and  $f_{\sigma'}$  is bounded by a constant multiple of  $\epsilon^2$ .

We now proceed with the global part of the argument, and thus return to the setting where K is a polyhedral Lagrangian surface in  $\mathbb{R}^4$ . The first step is to globalise the choice of L:

**Definition 4.** A conormal fibration dual to K is a smoothly varying family  $L_z$  of (affine)-linear Lagrangian planes in  $\mathbb{R}^4$ , parametrised by  $z \in K$ , which locally satisfies the properties from Definition 3.

**Lemma 8.** The surface K admits a dual conormal fibration which, near vertex, agrees with the choice given by Corollary 1.

Proof. Lemma 4 implies that the choices near the vertices may be extended to the edges. Choosing a normal vector field to one of the faces that meets along an edge determines matched normals, and the extension to the interior of the faces is then standard, as the space of Lagrangian planes transverse to a given one is contractible.

The conormal fibration determines a subset S (K) of the space of C 1 -functions consisting of those functions which are smooth in the interior of each face, and which are smoothing functions in the sense of Definition 2 near each edge and vertex.

Lemma 9. There exist smoothing functions for K of arbitrarily small C 1 -norm.

Proof. Choose a partition of unity P α ρ<sup>α</sup> = 1 on K, indexed by the strata of K, so that ρ<sup>α</sup> is supported in the open star of α (the union of all strata adjacent to it). Lemma 7 asserts the existence of smoothing functions f<sup>α</sup> of arbitrarily small C 1 -norm defined on the open star of α. The function P α ραf<sup>α</sup> satisfies the desired property.

We now arrive at the proof of the main result, which mostly consists of assembling together all the previous steps:

Proof of Proposition 1. We have a neighbourhood νK of K in R 4 in which the conormal fibres L<sup>z</sup> are disjoint. The statement of Lemma 6 and its proof apply verbatim to this space, replacing K by Σ. The existence of sufficiently many global smoothing functions is guaranteed by Lemma 9.

As a consequence, we obtain a sequence K<sup>i</sup> of smooth embedded Lagrangians, which are all isotopic to K by a piecewise smooth isotopy and converge to it, that are moreover graphs of differentials of smooth functions (over each other) with respect to the fibration {Lz}. This graphical description yields a smooth Hamiltonian path of graphical Lagrangians connecting K<sup>i</sup> to Ki+1, and smoothing the concatenation of these paths yields the desired result.

# References

- [1] M. W. Hirsch. Differential Topology. Graduate Texts in Mathematics. Springer, 1976.
- [2] D. McDuff and D. Salamon. Introduction to Symplectic Topology. Oxford Mathematical Monographs. Oxford University Press, third edition, 2017.
