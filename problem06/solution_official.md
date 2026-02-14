# Official Solution: Problem 6: Existence of Large Epsilon-Light Subsets in Graphs (Daniel Spielman)

> Source: *First Proof Solutions and Comments*, Abouzaid, Blumberg, Hairer, Kileel, Kolda, Nelson, Spielman, Srivastava, Ward, Weinberger, Williams. February 14, 2026.

---

## Authors' Commentary on AI-Generated Solutions


Gemini asserted that it presented a proof of the existence of a constant that satisfied Question 6. But, after some correct statements, it presented a very vague explanation of how the proof could be finished. To me, it seems unlikely that the approach can be turned into a correct proof.

ChatGPT 5.2 Pro asserted that it could not answer the question. So, it instead offered a correct upper bound of 1/2 on the constant, if it exists.

---

## Official Solution

# Light Sets of Vertices Daniel Spielman, Jan. 30, 2026

Throughout this note, G = (V, E, w) will be a weighted graph with n vertices. For an edge  $(s,t) \in E$ , we let w(s,t) be its weight. For two vertex sets, S and T, the subgraph  $G_{S,T}$  of G has vertex set V, but only the edges going between vertices in S and T. We write  $G_S$  for the graph that only contains the edges between vertices in S.

The matrix L is the Laplacian of G, which we recall may be defined by

$$L = \sum_{(s,t)\in E} w(s,t)(\boldsymbol{\delta}_s - \boldsymbol{\delta}_t)(\boldsymbol{\delta}_s - \boldsymbol{\delta}_t)^T,$$

where  $\delta_s$  is the elementary unit vector with a 1 in position s. We let  $L_S$  denote the Laplacian of  $G_S$ . As  $G_S$  and G have been defined to have the same vertex set,  $L_S$  has the same dimension as L.

**Lemma 0.1.** For every weighted graph G = (V, E, w) with n vertices, and for every  $0 < \epsilon < 1$ , there is an  $S \subseteq V$  of size at least  $\epsilon n/42$  so that

$$\epsilon L \succcurlyeq L_S$$
.

We call such a set of vertices S an  $\epsilon$ -light set. A set S is 0-light if and only if it is independent, and we could view lightness as a qualitative measure of independence. We might have called it "spectral independence," if that term were not already in use.

This lemma was proved by Daniel Spielman while working on the paper "Sparsified Cholesky Solvers for SDD linear systems", written with Richard Peng and Yin-Tat Lee [LPS15]. We decided not to include the lemma in that paper because, while it could be used to obtain interesting variants of some results, it was not necessary for the main results in that paper. That paper evolved into the paper "Sparsified Cholesky and Multigrid Solvers for Connection Laplacians," written with Rasmus Kyng, Yin Tat Lee, Richard Peng and Sushant Sachdeva [KLP+16].

# 1 Proof Strategy

We define  $L_{S,T}$  to be the Laplacian of  $G_{S,T}$ . For a vertex t and a subset of vertices S, we define  $L_{S,t}$  to be the Laplacian of  $G_{S,\{t\}}$ .

For a matrix L, we write its pseudo-inverse as  $L^{\dagger}$ . We write  $L^{\dagger/2}$  for the square root of the pseudo-inverse. We will prove the following statement that is equivalent to Lemma 0.1

$$\left\| L^{\dagger/2} L_S L^{\dagger/2} \right\| \le \epsilon.$$

We will find it convenient to multiply all Laplacian matrices on the left and right by  $L^{\dagger/2}$ . So, we define

$$\widetilde{L}_S = L^{\dagger/2} L_S L^{\dagger/2}, \quad \widetilde{L}_{S,T} = L^{\dagger/2} L_{S,T} L^{\dagger/2}, \quad \widetilde{L}_{S,t} = L^{\dagger/2} L_{S,t} L^{\dagger/2},$$

and recall that  $L^{\dagger/2}LL^{\dagger/2}\stackrel{\mathrm{def}}{=}\Pi$  is a symmetric projection matrix.

We are going to build up S in a greedy fashion. We will begin with a singleton set, and then add one vertex at a time. As we add vertices to S, we will need to maintain bounds on two quantities: a modification of the upper barrier function from [BSS12] and the sum of the leverage scores of edges between S and V \ S.

The leverage score of an edge (s, t) is defined to be w(s, t) times the effective resistance between s and t:

$$\ell(s,t) = w(s,t)(\boldsymbol{\delta}_s - \boldsymbol{\delta}_t)^T L^{\dagger}(\boldsymbol{\delta}_s - \boldsymbol{\delta}_t) = \operatorname{Tr}\left(w(s,t)(\boldsymbol{\delta}_s - \boldsymbol{\delta}_t)(\boldsymbol{\delta}_s - \boldsymbol{\delta}_t)^T L^{\dagger}\right) = \operatorname{Tr}\left(L_{\{s\},\{t\}}L^{\dagger}\right)$$

.

.

For vertices s and t for which (s, t) is not an edge, we define ℓ(s, t) = 0. For subsets of vertices S and T, we define

$$\ell(S,T) \stackrel{\text{def}}{=} \sum_{s \in S} \sum_{t \in T} \ell(s,t) = \sum_{s \in S} \sum_{t \in T: (s,t) \in E} \ell(s,t),$$

and

$$\ell(S) \stackrel{\text{def}}{=} \ell(S, V - S).$$

Claim 1.1. For S and T subsets of vertices, ℓ(S, T) = Tr <sup>L</sup>eS,T .

Proof. From the definition of the Laplacian of a graph, we have LS,T = P s∈S P <sup>t</sup>∈<sup>T</sup> L{s},{t} . So,

$$\operatorname{Tr}\left(\widetilde{L}_{S,T}\right) = \operatorname{Tr}\left(L^{\dagger/2}L_{S,T}L^{\dagger/2}\right) = \operatorname{Tr}\left(L_{S,T}L^{\dagger}\right)$$
$$= \sum_{s \in S} \sum_{t \in T} \operatorname{Tr}\left(L_{\{s\},\{t\}}L^{\dagger}\right) = \sum_{s \in S} \sum_{t \in T} \ell(s,t) = \ell(S,T).$$

We modify the BSS barrier function to make it better suited to matrices of rank at most σ by only incorporating the largest σ eigenvalues of the matrix. For a matrix A with eigenvalues λ<sup>1</sup> ≥ λ<sup>2</sup> ≥ · · · ≥ λn, and a u > λ1, we define

$$\Phi_{\sigma}^{u}(A) \stackrel{\text{def}}{=} \sum_{i=1}^{\sigma} \frac{1}{u - \lambda_{i}}.$$

If u ≤ λ1, we define Φ<sup>u</sup> σ (A) = ∞. We overload the definition of Φ by setting

$$\Phi_{\sigma}^{u}(S) \stackrel{\text{def}}{=} \Phi_{\sigma}^{u}(\widetilde{L}_{S}).$$

Our objective is to find a set S of size σ so that Φ<sup>ϵ</sup> σ (S) < ∞.

We deal with this barrier function by considering a modified trace of a matrix that only sums the largest σ eigenvalues of its argument:

$$\operatorname{Tr}_{\sigma}(A) \stackrel{\operatorname{def}}{=} \sum_{i=1}^{\sigma} \lambda_{i},$$

where the eigenvalues of A are λ<sup>1</sup> ≥ λ<sup>2</sup> ≥ · · · ≥ λn. We then have Φ<sup>u</sup> σ (A) = Tr<sup>σ</sup> (uI − A) −1 In all cases we consider, the argument of Tr<sup>σ</sup> is a diagonalizable matrix with real eigenvalues.

For the rest of this note, define

$$\delta \stackrel{\mathrm{def}}{=} \frac{21}{n}, \quad \phi \stackrel{\mathrm{def}}{=} \frac{n}{21}, \quad \mathrm{and} \quad \sigma \stackrel{\mathrm{def}}{=} \lfloor \epsilon n/42 \rfloor \,.$$

We will prove Lemma 0.1 by iteratively applying the following lemma.

Lemma 1.2. If |S| ≤ σ, ℓ(S) ≤ 4 |S|, and Φ u σ (S) ≤ ϕ, then there is a t ̸∈ S so that

$$\Phi_{\sigma}^{u+\delta}(S \cup \{t\}) \le \phi$$
 and  $\ell(S \cup \{t\}) \le \ell(S) + 4$ .

Proof. Lemma 2.1 says that for more than half the t ̸∈ S, ℓ(S ∪ {t}) ≤ ℓ(S) + 4. And, under the conditions of the lemma, Lemma 2.5 says that for at least half the t ̸∈ S, Φ<sup>u</sup> σ (S ∪ {t}) ≤ ϕ. So, there is a t ̸∈ S that satisfies both conditions.

Proof of Lemma 0.1. Set u<sup>0</sup> = ϵ/2 and let S<sup>0</sup> = {v0} an arbitrary v<sup>0</sup> ∈ V . As GS<sup>0</sup> has no edges,

$$\Phi_{\sigma}^{u_0}(S_0) = \sigma/u_0 \le \frac{n}{21} = \phi.$$

By applying Lemma 1.2 σ times, we inductively construct a set S of σ + 1 vertices so that ℓ(S) ≤ 4σ and Φu0+σδ σ (S) <sup>≤</sup> <sup>ϕ</sup>. This implies that all of the eigenvalues of <sup>L</sup>e<sup>S</sup> are at most

$$u_0 + \sigma \delta = \frac{\epsilon}{2} + \sigma \frac{21}{n} \le \epsilon.$$

# 2 Proofs

Lemma 2.1. Let S ⊂ V . Then, for more than half the t not in S,

$$\ell(S \cup \{t\}) \le \ell(S) + 4.$$

Proof. Recall ℓ(S ∪ {t}) = ℓ(S ∪ {t} , V − (S ∪ {t})). For t ̸∈ S, we use the inequality

$$\ell(S \cup \{t\}, V - (S \cup \{t\})) \le \ell(S \cup \{t\}, V - S) = \ell(S) + \ell(t, V - S).$$

So, it suffices to show that for more than half the t ̸∈ S, ℓ(t, V − S) ≤ 4. This follows from the non-negativity of ℓ and Claim 2.2 which shows that

$$\sum_{t \in V - S} \ell(t, V - S) < 2 |V - S|.$$

Claim 2.2. For every T ⊂ V ,

$$\sum_{t \in T} \ell(t, T) \le 2(|T| - 1).$$

Proof.

$$\sum_{t \in T} \ell(t,T) = \sum_{t \in T} \operatorname{Tr} \left( L_{\{t\},T} L^{\dagger} \right) = 2 \operatorname{Tr} \left( L_T L^{\dagger} \right).$$

To show that  $\operatorname{Tr}(L_T L^{\dagger}) < |T|$ , observe that  $L_T \leq L$ , so all the eigenvalues of  $L_T L^{\dagger}$  are between 0 and 1. Because  $L_T$  has rank at most |T| - 1, at most |T| - 1 eigenvalues of  $L_T L^{\dagger}$  are non-zero.

For convenience, we now state a few key properties of the function  $\text{Tr}_{\sigma}$  of a matrix. We begin with its defect: it is not additive. But, Ky Fan's eigenvalue inequality (see Theorem 4.3.47a of [HJ12]) tells us that it is subadditive:

$$\operatorname{Tr}_{\sigma}(A+B) \le \operatorname{Tr}_{\sigma}(A) + \operatorname{Tr}_{\sigma}(B)$$
. (1)

Most of the properties of  $\text{Tr}_{\sigma}$  that we find helpful follow from the fact that, for matrices A and B, AB has the same non-zero eigenvalues as BA, counted with multiplicity.

**Proposition 2.3.** For symmetric matrices A and B,

- a.  $\operatorname{Tr}_{\sigma}(A) = \max_{U} \operatorname{Tr}(UAU^{T})$ , where the maximum is taken over all orthogonal matrices of rank  $\sigma$ .
- b. If A is positive semidefinite, then  $\operatorname{Tr}_{\sigma}(AB) = \operatorname{Tr}_{\sigma}(BA)$ .
- c. If A and B are positive semidefinite, then  $\operatorname{Tr}_{\sigma}(AB) \geq 0$ .
- d. If  $A \leq B$ , then  $\operatorname{Tr}_{\sigma}(A) \leq \operatorname{Tr}_{\sigma}(B)$ .
- e. If C is positive semidefinite and  $A \leq B$ , then  $\operatorname{Tr}_{\sigma}(AC) \leq \operatorname{Tr}_{\sigma}(BC)$ .

Proof. Part a is Ky Fan's maximum principle, proved in [Fan49]. Part b is a direct consequence of the facts that AB has n real eigenvalues if A is positive semidefinite, and AB and BA have the same non-zero eigenvalues. Part c follows from the fact that all eigenvalues of the product of positive semidefinite matrices are non-negative. Part d follows from using (1) to show  $\operatorname{Tr}_{\sigma}(A) \leq \operatorname{Tr}_{\sigma}(B) + \operatorname{Tr}_{\sigma}(A-B) \leq \operatorname{Tr}_{\sigma}(B)$ , using the fact that A-B is negative semidefinite and so  $\operatorname{Tr}_{\sigma}(A-B) \leq 0$ . To derive part e from part d, let V be a matrix so that  $V^{T}V = C$ , and apply b to show the conclusion is equivalent to  $\operatorname{Tr}_{\sigma}(VAV^{T}) \leq \operatorname{Tr}_{\sigma}(VBV^{T})$ , which follows from  $VAV^{T} \leq VBV^{T}$ .

Note that  $\widetilde{L}_{S \cup \{t\}} = \widetilde{L}_S + \widetilde{L}_{S,t}$ . To show that we can choose a  $t \notin S$  that does not increase the barrier function, we employ the following adaptation of Lemma 19 of [SHS15], which in turn is an adaptation of Lemma 3.3 from [BSS12]. We include a proof for completeness.

**Lemma 2.4.** Let A and B be positive semidefinite matrices,  $\delta > 0$ , and let  $M = (u + \delta)I - A$ . If  $\Phi_{\sigma}^{u}(A) < \infty$  and

$$\frac{\operatorname{Tr}_{\sigma}\left(M^{-2}B\right)}{\Phi_{\sigma}^{u}(A) - \Phi_{\sigma}^{u+\delta}(A)} + \operatorname{Tr}_{\sigma}\left(M^{-1}B\right) < 1,\tag{2}$$

then  $\Phi_{\sigma}^{u+\delta}(A+B) \leq \Phi_{\sigma}^{u}(A)$ .

*Proof.* Our assumption that  $\Phi_{\sigma}^{u}(A) < \infty$  implies that  $M, M^{-1}$ , and  $M^{-2}$  are all positive definite. Thus, Proposition 2.3c implies that both terms in (2) are non-negative. Let C be a matrix for which  $B = CC^T$ , and so by Proposition 2.3b  $\operatorname{Tr}_{\sigma}(M^{-1}B) = \operatorname{Tr}_{\sigma}(C^TM^{-1}C) < 1$ . Recall  $\Phi_{\sigma}^{u+\delta}(A+B) = \operatorname{Tr}_{\sigma}((M-CC^T)^{-1})$ . By the Sherman-Morrison-Woodbury formula,

$$(M - CC^{T})^{-1} = M^{-1} + M^{-1}C(I - C^{T}M^{-1}C)^{-1}C^{T}M^{-1}.$$

As  $||C^T M^{-1}C|| \leq \operatorname{Tr}_{\sigma}(C^T M^{-1}C) < 1$ , we know that right-hand term is positive definite, and thus all eigenvalues of A + B are less than  $u + \delta$ . Now, (1) implies

$$\Phi_{\sigma}^{u+\delta}(A+B) \le \operatorname{Tr}_{\sigma}(M^{-1}) + \operatorname{Tr}_{\sigma}(M^{-1}C(I-C^{T}M^{-1}C)^{-1}C^{T}M^{-1}).$$

By Proposition 2.3b,

$$\operatorname{Tr}_{\sigma} \left( M^{-1}C(I - C^{T}M^{-1}C)^{-1}C^{T}M^{-1} \right) = \operatorname{Tr}_{\sigma} \left( (I - C^{T}M^{-1}C)^{-1}C^{T}M^{-2}C \right)$$

As  $||C^T M^{-1} C|| \le \operatorname{Tr}_{\sigma} (C^T M^{-1} C) < 1$ ,  $(I - C^T M^{-1} C)^{-1} \le (1 - \operatorname{Tr}_{\sigma} (C^T M^{-1} C))^{-1} I$ , and by Proposition 2.3d,

$$\operatorname{Tr}_{\sigma}\left((I - C^{T}M^{-1}C)^{-1}C^{T}M^{-2}C\right) \le \frac{\operatorname{Tr}_{\sigma}\left(C^{T}M^{-2}C\right)}{1 - \operatorname{Tr}_{\sigma}\left(C^{T}M^{-1}C\right)}.$$

Writing  $\operatorname{Tr}_{\sigma}(M^{-1}) = \Phi_{\sigma}^{u}(A) - (\Phi_{\sigma}^{u}(A) - \Phi_{\sigma}^{u+\delta}(A))$ , we obtain

$$\Phi_{\sigma}^{u+\delta}(A+B) \le \Phi_{\sigma}^{u}(A) - (\Phi_{\sigma}^{u}(A) - \Phi_{\sigma}^{u+\delta}(A)) + \frac{\operatorname{Tr}_{\sigma}\left(C^{T}M^{-2}C\right)}{1 - \operatorname{Tr}_{\sigma}\left(C^{T}M^{-1}C\right)},$$

which (2) and Proposition 2.3b imply is at most  $\Phi_{\sigma}^{u}(A)$ .

We will apply this result with  $A = \widetilde{L}_S$  and  $B = \widetilde{L}_{S,t}$ . When these terms, along with u and  $\delta$  are given, it will be convenient to write

$$U(S,t) \stackrel{\text{def}}{=} \frac{\operatorname{Tr}_{\sigma} \left( M^{-2} \widetilde{L}_{S,t} \right)}{\Phi_{\sigma}^{u}(S) - \Phi_{\sigma}^{u+\delta}(S)} + \operatorname{Tr}_{\sigma} \left( M^{-1} \widetilde{L}_{S,t} \right).$$

**Lemma 2.5.** If  $|S| \leq \sigma$ ,  $\Phi_{\sigma}^{u}(S) \leq \phi$ , and  $\ell(S) \leq 4|S|$ , then for at least half the  $t \notin S$ ,

*Proof.* We will prove that

$$\sum_{t \notin S} U(S, t) \le \frac{5}{\delta} + 5\phi.$$

As U(S,t) is non-negative, this implies that for at least half the  $t \notin S$ ,

$$U(S,t) \le \frac{2}{n-|S|} \left(\frac{5}{\delta} + 5\phi\right) \le \frac{2}{n} \frac{42}{41} \left(\frac{5n}{21} + \frac{5n}{21}\right) < 1.$$

We need to upper bound the terms  $\operatorname{Tr}_{\sigma}\left(M^{p}\widetilde{L}_{S,t}\right)$  for  $p \in \{-1,-2\}$ . We do this by breaking each term into two parts. Let  $\Pi_{S}$  be the symmetric projection onto the span of  $\widetilde{L}_{S}$  and let  $\Pi_{T} = I - \Pi_{S}$ . As  $M = (u + \delta)(\Pi_{S} + \Pi_{T}) - \widetilde{L}_{S}$ ,  $\Pi_{T}\Pi_{S} = \Pi_{T}\widetilde{L}_{S} = 0$ , and  $\Pi_{S}^{p} = \Pi_{S}$ ,

$$M^{p} = (u+\delta)^{p}\Pi_{T} + \left((u+\delta)\Pi_{S} - \widetilde{L}_{S}\right)^{p}.$$

By the subadditivity of  $\text{Tr}_{\sigma}$  we conclude

$$\operatorname{Tr}_{\sigma}\left(M^{p}\widetilde{L}_{S,t}\right) \leq \operatorname{Tr}_{\sigma}\left((u+\delta)^{p}\Pi_{T}\widetilde{L}_{S,t}\right) + \operatorname{Tr}_{\sigma}\left(\left((u+\delta)\Pi_{S} - \widetilde{L}_{S}\right)^{p}\widetilde{L}_{S,t}\right).$$

The term invovling  $\Pi_S$  is addressed by Claim 2.6, which says

$$\sum_{t \notin S} \operatorname{Tr}_{\sigma} \left( \left( (u + \delta) \Pi_{S} - \widetilde{L}_{S} \right)^{p} \widetilde{L}_{S, t} \right) \leq \operatorname{Tr}_{\sigma} \left( M^{p} \right).$$

For the other term, we recall that  $\Pi_T$  and  $\widetilde{L}_{S,t}$  are positive semidefinite and so their product has only non-negative eigenvalues to show

$$\operatorname{Tr}_{\sigma}\left((u+\delta)^{p}\Pi_{T}\widetilde{L}_{S,t}\right) \leq \operatorname{Tr}\left((u+\delta)^{p}\Pi_{T}\widetilde{L}_{S,t}\right) = (u+\delta)^{p}\operatorname{Tr}\left(\Pi_{T}\widetilde{L}_{S,t}\right) \leq (u+\delta)^{p}\operatorname{Tr}\left(\widetilde{L}_{S,t}\right).$$

Claim 1.1 tells us that this equals  $(u + \delta)^p \ell(S, t)$ , giving

$$\sum_{t \notin S} \operatorname{Tr}_{\sigma} \left( (u+\delta)^{p} \Pi_{T} \widetilde{L}_{S,t} \right) \leq (u+\delta)^{p} \sum_{t \notin S} \ell(S,t) = (u+\delta)^{p} \ell(S) \leq (u+\delta)^{p} 4 |S|.$$

To combine these terms, note that all the eigenvalues of M are at most  $(u + \delta)$ , and thus for p < 0 all the eigenvalues of  $M^p$  are at least  $(u + \delta)^p$ . This tells us that  $\operatorname{Tr}_{\sigma}(M^p) \geq \sigma(u + \delta)^p \geq |S|(u + \delta)^p$ . We conclude that

$$\sum_{t \notin S} \operatorname{Tr}_{\sigma} \left( M^{p} \widetilde{L}_{S,t} \right) \leq 5 \operatorname{Tr}_{\sigma} \left( M^{p} \right).$$

To finish, we return to

$$\sum_{t \notin S} U(S, t) = \sum_{t \notin S} \frac{\operatorname{Tr}_{\sigma} \left( M^{-2} \widetilde{L}_{S, t} \right)}{\Phi_{\sigma}^{u}(S) - \Phi_{\sigma}^{u + \delta}(S)} + \sum_{t \notin S} \operatorname{Tr}_{\sigma} \left( M^{-1} \widetilde{L}_{S, t} \right) \leq \frac{5 \operatorname{Tr}_{\sigma} \left( M^{-2} \right)}{\Phi_{\sigma}^{u}(S) - \Phi_{\sigma}^{u + \delta}(S)} + 5 \operatorname{Tr}_{\sigma} \left( M^{-1} \right).$$

The right-hand term is at most  $5\Phi_{\sigma}^{u+\delta}(S)$ , and Claim 2.7 shows that the left-hand term is at most  $\frac{5}{\delta}$ . Summing these together gives the result.

Claim 2.6. Assume that  $|S| \leq \sigma$ . For  $M = (u + \delta)I - \widetilde{L}_S$ , and nonzero real p,

$$\sum_{t \notin S} \operatorname{Tr}_{\sigma} \left( \left( (u + \delta) \Pi_{S} - \widetilde{L}_{S} \right)^{p} \widetilde{L}_{S, t} \right) \leq \operatorname{Tr}_{\sigma} \left( M^{p} \right).$$

*Proof.* Because both  $\widetilde{L}_{S,t}$  and  $\left((u+\delta)\Pi_S - \widetilde{L}_S\right)^p$  are positive semidefinite, the eigenvalues of their product are nonnegative, and so

$$\operatorname{Tr}_{\sigma}\left(\left((u+\delta)\Pi_{S}-\widetilde{L}_{S}\right)^{p}\widetilde{L}_{S,t}\right) \leq \operatorname{Tr}\left(\left((u+\delta)\Pi_{S}-\widetilde{L}_{S}\right)^{p}\widetilde{L}_{S,t}\right).$$

As  $\sum_{t \notin S} \widetilde{L}_{S,t} = \widetilde{L}_{S,T} \preccurlyeq I$ , Proposition 2.3d implies

$$\begin{split} \sum_{t \notin S} \operatorname{Tr} \left( \left( (u + \delta) \Pi_S - \widetilde{L}_S \right)^p \widetilde{L}_{S,t} \right) &= \operatorname{Tr} \left( \left( (u + \delta) \Pi_S - \widetilde{L}_S \right)^p \widetilde{L}_{S,T} \right) \\ &\leq \operatorname{Tr} \left( \left( (u + \delta) \Pi_S - \widetilde{L}_S \right)^p \right) &= \operatorname{Tr} \left( \Pi_S \left( (u + \delta) I - \widetilde{L}_S \right)^p \Pi_S \right) &= \operatorname{Tr} \left( \Pi_S M^p \Pi_S \right). \end{split}$$

By Ky Fan's maximum principle (Proposition 2.3a) this latter term is at most  $\operatorname{Tr}_{\sigma}(M^p)$ .

#### Claim 2.7.

$$\Phi_{\sigma}^{u}(S) - \Phi_{\sigma}^{u+\delta}(S) \ge \delta \operatorname{Tr}_{\sigma} (M^{-2}).$$

*Proof.* Let  $\lambda_1, \ldots, \lambda_{\sigma}$  be the largest  $\sigma$  eigenvalues of  $\widetilde{L}_S$ . Then,

$$\Phi_{\sigma}^{u}(S) - \Phi_{\sigma}^{u+\delta}(S) = \sum_{i=1}^{\sigma} \frac{1}{u - \lambda_{i}} - \sum_{i=1}^{\sigma} \frac{1}{u + \delta - \lambda_{i}}$$

$$= \sum_{i=1}^{\sigma} \frac{\delta}{(u - \lambda_{i})(u + \delta - \lambda_{i})}$$

$$\geq \sum_{i=1}^{\sigma} \frac{\delta}{(u + \delta - \lambda_{i})^{2}}.$$

$$= \delta \operatorname{Tr}_{\sigma} (M^{-2}).$$

## References

[BSS12] Joshua Batson, Daniel A Spielman, and Nikhil Srivastava. Twice-Ramanujan sparsifiers. SIAM Journal on Computing, 41(6):1704–1721, 2012.

- [Fan49] Ky Fan. On a theorem of Weyl concerning eigenvalues of linear transformations I. Proceedings of the National Academy of Sciences of the United States of America, 35(11):652, 1949.
- [HJ12] Roger A Horn and Charles R Johnson. *Matrix analysis*. Cambridge university press, 2012.
- [KLP+16] Rasmus Kyng, Yin Tat Lee, Richard Peng, Sushant Sachdeva, and Daniel A Spielman. Sparsified Cholesky and multigrid solvers for connection Laplacians. In Proceedings of the forty-eighth annual ACM symposium on Theory of Computing, pages 842–850. ACM, 2016.

- [LPS15] Yin Tat Lee, Richard Peng, and Daniel A. Spielman. Sparsified Cholesky solvers for SDD linear systems. CoRR, abs/1506.08204, 2015.
- [SHS15] Marcel K De Carli Silva, Nicholas JA Harvey, and Cristiane M Sato. Sparse sums of positive semidefinite matrices. ACM Transactions on Algorithms (TALG), 12(1):1– 17, 2015.

