# Official Solution: Problem 1: Equivalence of the Phi^4_3 Measure Under Smooth Shifts (Martin Hairer)

> Source: *First Proof Solutions and Comments*, Abouzaid, Blumberg, Hairer, Kileel, Kolda, Nelson, Spielman, Srivastava, Ward, Weinberger, Williams. February 14, 2026.

---

## Authors' Commentary on AI-Generated Solutions


<span id="page-5-1"></span>In this case, a note with a very short sketch of proof (far short of the level of detail one would expect for a published article) was posted on the author's homepage some time ago. The answer given by GPT-Pro simply quotes that note, claiming that it contains a detailed proof of the result. This is incorrect and it is despite the LLM being specifically instructed to comply with "mathematics publication" levels of scholarship. (Taking for granted a result that is merely stated in an unpublished note with a very rough sketch of proof is not considered acceptable in the mathematics literature.)

Another behaviour we observed was that the LLM would take as a premise the (wrong!) statement that the Φ 4 <sup>3</sup> measure is equivalent to the free field measure, from which it then

correctly deduces the (incorrect) claim that the Φ 4 <sup>3</sup> measure is quasi-invariant under smooth shifts.

---

## Official Solution

# (Lack of) quasi-shift invariance of the $\Phi^4_3$ measure

Martin Hairer

## 1 Question

Let  $\mathbb{T}^3$  be the three dimensional unit size torus and let  $\mu$  be the  $\Phi_3^4$  measure on the space of distributions  $\mathcal{D}'(\mathbb{T}^3)$ . Let  $\psi:\mathbb{T}^3\to\mathbb{R}$  be a smooth function that is not identically zero and let  $T_\psi:\mathcal{D}'(\mathbb{T}^3)\to\mathcal{D}'(\mathbb{T}^3)$  be the shift map given by  $T_\psi(u)=u+\psi$  (with the usual identification of smooth functions as distributions). Prove or disprove the statement "the measures  $\mu$  and  $T_\psi^*\mu$  are equivalent". Here, equivalence of measures is in the sense of having the same null sets and  $T_\psi^*$  denotes the pushforward under  $T_\psi$ .

#### 2 Some Context

One of the very few interacting quantum field theories that can be rigorously constructed is the so-called (bosonic)  $\Phi^4$  theory in (space-time) dimensions 2 and 3. It has long been known that in dimension 2 and finite volume there is a natural identification between the Hilbert space of the interacting theory and that of the corresponding free theory. On the other hand, Glimm [Gli68] observed that this is no longer the case in dimension 3. At the level of the corresponding Euclidean theories (which are represented by probability measures on the space of Schwartz distributions on the corresponding space-time), this translates into the fact that the  $\Phi^4$  measure  $\mu$  and the corresponding free field measure  $\nu$  are equivalent in dimension 2 but mutually singular in dimension 3. In fact, there is a sense in which the dimension that delimits between the two behaviors is 8/3. It is then natural to ask in which dimensions  $\mu$  has the weaker property that  $\mu$  and  $T^*_{\psi}\mu$  are equivalent for smooth  $\psi$ . Here it turns out that the borderline dimension is 3, and the question probes on which side it falls.

#### 2.1 An incorrect heuristic

Regarding the proof, a tempting heuristic is to use the fact that one should think of  $\mu$  as having the density with respect to "Lebesgue measure on  $\mathcal{D}'$ " (which of course doesn't exist) proportional to

$$\exp\left(-\int_{\mathbb{T}^3} \left(\frac{1}{2} |\nabla \Phi(x)|^2 + \frac{1}{4} |\Phi(x)|^4 - \frac{C}{2} |\Phi(x)|^2\right) dx\right),\,$$

Notations

where C is a (diverging) constant of the form  $C=3c_1-9c_2$ , where  $c_1$  is the expectation of  $|\Phi(x)|^2$  under the free field measure  $\nu$  (which is of course infinite) and  $c_2$  is an additional logarithmically divergent constant. The density of  $T_{\psi}^*\mu$  with respect to  $\mu$  is then formally given by

$$\exp\left(-\int_{\mathbb{T}^3} \left(\frac{1}{2} |\nabla \psi(x)|^2 + \frac{1}{4} |\psi(x)|^4 + \Phi(x) \Delta \psi(x) - \Phi(x) \psi^3(x) - \psi(x) (\Phi^3(x) - C\Phi(x)) + \frac{\psi^2(x)}{2} (3\Phi^2(x) - C)\right) dx\right),$$

Since the terms on the first line are well-defined for smooth  $\psi$  and one expects  $\Phi^3 - 3c_1\Phi$  and  $\Phi^2 - c_1$  to be quite well-behaved, the additional logarithmically divergent terms proportional to  $c_2$  cause this "density" to diverge, suggesting (correctly) that  $\mu$  and  $T_{\psi}^*\mu$  are mutually singular.

There are at least two problems with such an approach. First,  $\Phi^3 - 3c_1\Phi$  does actually *not* define a random distribution, whether  $\Psi$  is distributed according to  $\mu$  or to the free field  $\nu$  (which guides the heuristic). This is because if it were, it would have a covariance behaving like  $|x-y|^{-3}$  around the diagonal, which is not integrable in dimension 3. The second problem is that such an argument suggests that, if  $\mu_n = \exp(-f_n)\nu$  for some "nice" probability measure  $\nu$  and functions  $f_n$  that fail to converge to a "nice" limit, then  $\mu_n$  fails to converge to a limit  $\mu$ . This of course is not true: for a suitable (diverging) sequence of constants  $c_n$ , the sequence  $f_n(x) = c_n + n\cos(nx)$  is such that if  $\nu$  is Lebesgue measure on [0,1], then  $\mu_n$  converges weakly to Lebesgue measure even though the log-densities  $f_n$  fail to converge. Any proof needs to be based on a different approach or to satisfactorily address these problems.

## Acknowledgement

The proof presented below is a simplified version of an argument that will appear as part of a joint publication with Jacopo Peroni. This in turn is based on the article [HKN24].

## **3** Notations

We fix a space-time white noise  $\xi$  on  $\mathbb{R} \times \mathbb{T}^3$ . We define † as the stationary solution to the linear equation

$$(\partial_t + 1 - \Delta)! = \mathcal{E}, \quad \text{on } \mathbb{R} \times \mathbb{T}^3.$$

(We use the convention that symbols represent random space-time distributions rather than elements of a regularity structure.) Starting from this process, we define  $\lor$  and  $\checkmark$  as its Wick square and cube respectively, which are given by

$$V = \lim_{N \to \infty} H_2(\uparrow_N, c_N) , \qquad \psi = \lim_{N \to \infty} H_3(\uparrow_N, c_N) ,$$

where  $| \cdot \rangle_N = P_N |$  and  $c_N = \mathbb{E} |_N^2$  (which is constant in space and time). Here,  $P_N$  denotes the projection onto Fourier modes with  $|k| \leq N$  and  $H_n$  denotes the nth

Hermite polynomial normalised such that  $H_0 \equiv 1$ ,  $H'_n = nH_{n-1}$ , and  $\mathbb{E}H_n(Z,1) = 0$  for a normal random variable Z. The first convergence takes place in the space of continuous functions of time with values in  $C^{-1-2\kappa}$ , while the second convergence takes place in the space-time parabolic space  $C^{-\frac{3}{2}-3\kappa}$ .

With these notations in place, we define Y as the stationary solution to

$$(\partial_t + 1 - \Delta)Y = V$$

and similarly for  $\psi$ . For a more comprehensive and pedagogical introduction to the general tree-like notation, we refer the reader to [MWX17]. We also write " $\prec$ " for Bony's paraproduct (in space) as defined for example in [GIP15, Sec. 2.1] and, given a random N-dependent process w, we will sometimes use the physicists shorthand notation : $w^k$ : instead of  $H_k(w, c_N)$ .

## 4 Solution and Proof

The statement is **false**. In particular, for any smooth function  $\psi \not\equiv 0$  and any choice of the parameters involved in the definition of  $\mu$  (mass and coupling constant, provided that the latter is non-zero), the measures  $\mu$  and  $T_{\psi}^*\mu$  are mutually singular.

For notational simplicity, we fix the mass and the coupling constant to 1, but this has no incidence on the proof. Our main starting point is the following statement, a proof of which can be found for example in [HM18] and [EW24, Lemma 4.19] for (4.1), combined with [CC18] for (4.2) (see Ansatz 2.11 there). Throughout this proof,  $\kappa > 0$  is chosen small enough ( $\kappa = 1/100$  is certainly sufficient).

**Proposition 4.1.** There exists a stationary process v that is almost surely continuous in time with values in  $C^{1-2\kappa}(\mathbb{T}^d)$  and such that the process

$$u = \uparrow - \dot{\Psi} + v , \qquad (4.1)$$

is stationary with fixed time distribution equal to  $\mu$ . Furthermore, the process v is such that

$$v = -3((v - \dot{\Upsilon}) \prec \Upsilon) + v^{\sharp} , \qquad (4.2)$$

where  $v^{\sharp}$  is continuous with values in  $C^{1+4\kappa}(\mathbb{T}^d)$ .

It was furthermore shown in [HKN24, Lemmas 3.1 & 3.4] (but see [Hai14] for a similar result using a slightly different regularisation) that the processes † and  $\psi$  are almost surely continuous in time with values in  $C^{-\frac{1}{2}-\kappa}$  and  $C^{\frac{1}{2}-3\kappa}$  respectively.

Before we proceed, we remind some notations and preliminary results. First of all, we define the additional diverging constant

$$c_{N,2} := \mathbb{E}\left[\bigvee_{N}\bigvee_{N}\right],\tag{4.3}$$

where  $\bigvee_N := P_N \bigvee$  and  $\bigvee_N := P_N \bigvee$ . The main ingredient of our proof is the event

$$B^{\gamma} := \left\{ u \in \mathcal{D}' : \lim_{N \to \infty} \left\langle (\log N)^{-\gamma} \left( H_3 \left( P_N u; c_N \right) + 9 c_{N,2} P_N u \right), \psi \right\rangle_{\mathbb{T}^3} = 0 \right\},\,$$

which will be used to distinguish between the shifted and the non-shifted measures.

We will also use the following two technical lemmas whose proofs can be found in Section 4.1. These are very similar to [HKN24, Lemma 3.11 and Lemma 3.12].

**Lemma 4.2.** Let  $\gamma > \frac{1}{2}$ . Then, for any fixed t > 0,

$$\lim_{N \to \infty} (\log N)^{-\gamma} : (\uparrow_N)^3 : (t) = 0$$

almost surely in  $C^{-\frac{3}{2}}(\mathbb{T}^3)$  and in  $L^p(\Omega; C^{-\frac{3}{2}}(\mathbb{T}^3))$  for any p > 0.

**Lemma 4.3.** For N large, one has  $c_{N,2} \gtrsim \log N$ .

The following results are essentially standard, but we recall their statements for later reference.

**Lemma 4.4.** For any polynomial P, the expression  $\uparrow_N P(\Lsh_N)$  converges almost surely to some finite limit in  $C^{-\frac{1}{2}-\kappa}$ .

*Proof.* By paralinearisation and standard commutator estimates (see [GIP15, Lems 2.4 & 2.6]) it suffices to consider the case P(x) = x. This is by now standard, see for example [CC18, Sec. 4.4].

**Lemma 4.5.** Let v be a process satisfying the decomposition (4.2). Then, the expressions  $|\cdot|_N^2: \dot{\vee}_N - 3c_{N,2}|_N$  and  $|\cdot|_N^2: v_N + 3c_{N,2}(v_N - \dot{\vee}_N)$  both converge almost surely to finite limits in  $C^{-1-2\kappa}$  as  $N \to \infty$ .

*Proof.* Regarding the first expression, its convergence was essentially for example in [CC18, Sec. 4.6]. (The approximation used there is slightly different, but the differences are unimportant.) Regarding the second expression, the claim follows from [CC18, Sec. 4.5] (modulo again unimportant changes in the approximation scheme), combined with the commutator estimate [GIP15, Lem. 2.4]. □

We now turn to the proof of the main claim. For this, we first claim that if u is as in (4.1), then, for any fixed t, one has  $u(t) \in B^{\gamma}$ . Indeed, writing  $u_N$  as a shorthand for

 $P_N u$  and expanding the Wick power, we have

$$(\log N)^{-\gamma} H_3(u_N; c_N) = (\log N)^{-\gamma} H_3((\uparrow_N - \psi_N + v_N); c_N)$$

$$= (\log N)^{-\gamma} \sum_{i=0}^3 \binom{3}{i} : (\uparrow_N)^i : (-\psi_N + v_N)^{3-i}$$

$$= (\log N)^{-\gamma} \sum_{i=0}^3 \sum_{j=0}^{3-i} \binom{3}{i} \binom{3-i}{j} : (\uparrow_N)^i : (-\psi_N)^j (v_N)^{3-i-j}$$

$$= (\log N)^{-\gamma} : \uparrow_N^3 : -3(\log N)^{-\gamma} : \uparrow_N^2 : \psi_N + 3(\log N)^{-\gamma} : \uparrow_N^2 : v_N$$

$$+ (\log N)^{-\gamma} \sum_{\substack{0 \le i+j \le 3 \\ (i,j) \ne (3,0), (2,1), (2,0)}} \binom{3}{i} \binom{3-i}{j} : (\uparrow_N)^i : (-\psi_N)^j (v_N)^{3-i-j}.$$

The first term  $(\log N)^{-\gamma}$ :  ${}_{N}^{3}$ : and the terms present in the last sum all converge to 0 by Lemma 4.2 (given that  $\gamma > \frac{1}{2}$ ), standard product estimates (e.g. [Bon81, Theorem 2.5] or [MWX17, Proposition 2.3]) and Lemma 4.4.

It therefore remains to show that  $-:\uparrow_N^2: \psi_N + :\uparrow_N^2: v_N + 3c_{N,2}u_N$  also converges to zero almost surely in the sense of distributions. We rewrite this term as

$$:\uparrow_N^2:v_N+3c_{N,2}(v_N-\psi_N)-(:\uparrow_N^2:\psi_N-3c_{N,2}\uparrow_N)$$
.

By Lemma 4.5 we know that this expression converges to an element of  $C^{-1-2\kappa}(\mathbb{T}^d)$ , whence we conclude that

$$\langle (\log N)^{-\gamma} \left( -: !_N^2 : \psi_N + : !_N^2 : v_N + 3c_{N,2}u_N \right), \psi \rangle \xrightarrow{N \to \infty} 0$$

almost surely, thus proving that  $\mu(B^{\gamma}) = 1$ .

In order to conclude the proof, it suffices to show that  $u + \psi \notin B^{\gamma}$ . For this, we expand similarly to before the expression appearing in this event as

$$(\log N)^{-\gamma} H_3 ((u_N + \psi_N); c_N) + (\log N)^{-\gamma} 9c_{N,2} (u_N + \psi_N)$$

$$= (\log N)^{-\gamma} \sum_{i=0}^{3} {3 \choose i} : (u_N)^i : (\psi_N)^{3-i} + (\log N)^{-\gamma} 9c_{N,2} (u_N + \psi_N)$$

$$= (\log N)^{-\gamma} : (u_N)^3 : + (\log N)^{-\gamma} 9c_{N,2} u_N$$

$$+ 3(\log N)^{-\gamma} : (u_N)^2 : (\psi_N) + 3(\log N)^{-\gamma} (u_N) (\psi_N)^2$$

$$+ (\log N)^{-\gamma} (\psi_N)^3 + (\log N)^{-\gamma} 9c_{N,2} \psi_N.$$

The sum of the first two terms was just shown to converge to 0 almost surely in  $C^{-\frac{3}{2}-3\kappa}(\mathbb{T}^d)$  for  $N\to\infty$ .

Since  $:u_N^2:$  and  $u_N$  both converge to finite distributional limits almost surely by Lemma 4.4, the next three terms also converge to 0 almost surely.

Concerning the last element however, we know from Lemma 4.3 that

$$(\log N)^{-\gamma} c_{N,2} \gtrsim (\log N)^{-\gamma+1} .$$

Since the contribution of this term to the expression in the event  $B^{\gamma}$  is given by

$$9(\log N)^{-\gamma} c_{N,2} \langle \psi_N, \psi \rangle$$
,

and since  $\langle \psi_N, \psi \rangle \to \|\psi\|^2 > 0$ , this diverges, whence we conclude that  $u + \psi \notin B^{\gamma}$  and therefore  $(T_{\psi}^* \mu)(B^{\gamma}) = 0$ , so that  $T_{\psi}^* \mu$  and  $\mu$  are mutually singular.

## 4.1 Proof of the lemmas

*Proof of Lemma 4.2.* We use the embedding  $W^{\beta,2p}\hookrightarrow W^{\beta-\frac{d}{2p},\infty}\hookrightarrow \mathcal{C}^{\beta-\frac{d}{2p}}$ , with  $\beta=-\frac{d}{2}$ . Using the definition of  $W^{\beta,2p}$  norm and the equivalence of moments for Gaussian polynomials, one has

$$\mathbb{E}\left[\left\|(\log N)^{-\gamma}:(\uparrow_N)^J:\right\|_{W^{-\frac{3}{2},2p}}^{2p}\right] \lesssim \int_{\mathbb{T}^d} \mathbb{E}\left[(\log N)^{-2\gamma}\left|\langle\nabla\rangle^{-\frac{3}{2}}:(\uparrow_N)^3:\right|^2\right]^p dx \right].$$

Since one has

$$\mathbb{E}\left[\left(\log N\right)^{-2\gamma}\left|\left\langle\nabla\right\rangle^{-\frac{3}{2}}:\left(\uparrow_{N}\right)^{3}:\right|^{2}\right] \lesssim \left(\log N\right)^{-2\gamma} \sum_{|\omega_{i}| \leq N} \left\langle\omega_{1} + \dots + \omega_{3}\right\rangle^{-3} \prod_{i=1}^{3} \left\langle\omega_{i}\right\rangle^{-2}$$
$$\lesssim \left(\log N\right)^{-2\gamma} \sum_{r_{1}=0}^{N} \frac{r_{1}^{2}}{\left(1 + r_{1}^{2}\right)^{\frac{5}{2}}} r_{1}^{2} \lesssim \left(\log N\right)^{-2\gamma+1},$$

the desired result follows from a standard Borel–Cantelli argument.

Next, we prove Lemma 4.3, which provides a lower bound on the parameter  $\gamma$ . This bound ensures that the event  $A^{\gamma}$  (or  $B^{\gamma}$ ) is distinguishable under the shifted measure when compared to the non-shifted one.

*Proof of Lemma 4.3.* Expanding the definition of  $c_{N,2} := \mathbb{E} [\bigvee_{N} \bigvee_{N}]$ , we get

$$c_{N,2} = 2 \sum_{\substack{\omega_1 + \omega_2 = \omega_3 \\ |\omega_i| \le N}} \int_{\mathbb{R}} \hat{P}_{t-u}(\omega_3) \int_{\mathbb{R}} \hat{P}_{t-u_1}(\omega_1) \hat{P}_{u-u_1}(-\omega_1) du_1$$

$$\times \int_{\mathbb{R}} \hat{P}_{t-u_2}(\omega_2) \hat{P}_{u-u_2}(-\omega_2) du_2 du$$

$$\simeq \sum_{\substack{\omega_1 + \omega_2 = \omega_3 \\ |\omega_i| \le N}} \int_{\mathbb{R}} e^{-|t-u|\langle\omega_3\rangle^2} \frac{e^{-|t-u|\langle\omega_1\rangle^2}}{\langle\omega_1\rangle^2} \frac{e^{-|t-u|\langle\omega_2\rangle^2}}{\langle\omega_2\rangle^2} du$$

$$\gtrsim \sum_{|\omega_i| < N} \frac{1}{\langle\omega_1\rangle^2} \frac{1}{\langle\omega_1\rangle^2} \frac{1}{\langle\omega_1\rangle^2 + \langle\omega_2\rangle^2 + \langle\omega_1 + \omega_2\rangle^2}$$

$$\gtrsim \sum_{|\omega_i| \le N} \frac{1}{1 + |\omega_1|^2} \frac{1}{1 + |\omega_2|^2} \frac{1}{1 + |\omega_1|^2 \vee |\omega_2|^2}$$

$$\gtrsim \sum_{|\omega_1| \le |\omega_2| \le N} \frac{1}{1 + |\omega_1|^2} \frac{1}{1 + |\omega_2|^4} .$$

Bounding the sum by an integral, we finally conclude that this expression is bounded from below by a multiple of

$$\int_0^N \frac{r^2}{1+r^2} \int_r^\infty \frac{s^2}{1+s^4} \, ds \, dr \gtrsim \int_0^N \frac{1}{1+r^2} \, dr \simeq \log N \;,$$

as claimed.  $\Box$ 

#### References

- [Bon81] Jean-Michel Bony. Calcul symbolique et propagation des singularités pour les équations aux dérivées partielles non linéaires. *Ann. Sci. École Norm. Sup.* (4), 14(2):209–246, 1981.
- [CC18] Rémi Catellier and Khalil Chouk. Paracontrolled distributions and the 3-dimensional stochastic quantization equation. *Ann. Probab.*, 46(5):2621–2679, 2018.
- [EW24] Salvador Esquivel and Hendrik Weber. A priori bounds for the dynamic fractional  $\phi^4$  model on  $\mathbb{T}^3$  in the full subcritical regime. *arXiv preprint*, 2024.
- [GIP15] Massimiliano Gubinelli, Peter Imkeller, and Nicolas Perkowski. Paracontrolled distributions and singular PDEs. *Forum Math. Pi*, 3:e6, 75pp, 2015.
- [Gli68] James Glimm. Boson fields with the : $\Phi^4$ : interaction in three dimensions. *Comm. Math. Phys.*, 10(1):1–47, 1968.
- [Hai14] Martin Hairer. A theory of regularity structures. *Invent. Math.*, 198(2):269–504, 2014.
- [HKN24] Martin Hairer, Seiichiro Kusuoka, and Hirotatsu Nagoji. Singularity of solutions to singular SPDEs. *arXiv preprint*, 2024.
- [HM18] Martin Hairer and Konstantin Matetski. Discretisations of rough stochastic PDEs. *Ann. Probab.*, 46(3):1651–1709, 2018.
- [MWX17] Jean-Christophe Mourrat, Hendrik Weber, and Weijun Xu. Construction of  $\Phi^4_3$  diagrams for pedestrians. In From particle systems to partial differential equations, volume 209 of Springer Proc. Math. Stat., pages 1–46. Springer, Cham, 2017.
