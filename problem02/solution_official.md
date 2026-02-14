# Official Solution: Problem 2: Existence of Whittaker Functions for Rankin-Selberg Integrals (Paul Nelson)

> Source: *First Proof Solutions and Comments*, Abouzaid, Blumberg, Hairer, Kileel, Kolda, Nelson, Spielman, Srivastava, Ward, Weinberger, Williams. February 14, 2026.

---

## Authors' Commentary on AI-Generated Solutions


In some attempts, the LLM constructed *W* depending on *π*, but the problem asks for a single *W* that works for all *π*. This is a critical condition; without it, the problem is much easier and the solution is well-known. In some (but not all) cases, the LLM noted that it had solved a weaker problem.

In the best attempt in our trial runs, ChatGPT 5.2 Pro identified a suitable choice of *W* and reduced (as in our solution) to exhibiting *V* for which the integral R GL*n*(o) *V* (*g*)*ψ*(−*Qgnn*) *dg* does not vanish. This nonvanishing is the key point.

ChatGPT then attempted to choose *V* so that the integrand is constant on its support, which, if possible, would make the nonvanishing clear. This strategy is unviable. For instance, when *n* = 1, *V* must be (a nonzero multiple of) a character of *F* <sup>×</sup> and the integral is a normalized Gauss sum; in particular, the integrand is typically non-constant. For larger *n*, the unviability follows similarly by considering the action of the center.

To identify the specific error in the attempted solution, we look for the first place asserting stronger support properties of *V* than are generally true. The culprit is the support condition claimed in the "standard Howe-vector existence result," which never holds: it contradicts the fact that *V* has a central character.

---

## Official Solution

#### 1. Statement

Let F be a non-archimedean local field with ring of integers  $\mathfrak{o}$ . Let  $\psi: F \to \mathbb{C}^{\times}$  be a nontrivial additive character of conductor  $\mathfrak{o}$ . We write

$$G_r := \operatorname{GL}_r(F),$$

and let  $N_r < G_r$  denote the subgroup of upper-triangular unipotent elements. We embed  $G_n \hookrightarrow G_{n+1}$  as the upper-left block. We write  $E_{ij}$  for the matrix with a 1 in the (i,j)-entry and 0 elsewhere.

A more precise form of the following "lemma" will appear in forthcoming joint work with Subhajit Jana. It says informally that pure unipotent translates of fixed vectors in the Whittaker model of a representation of  $G_{n+1}$  may serve as test vectors for Rankin–Selberg integrals against all representations of  $G_n$  with a given conductor.

**Theorem 1.** Let  $\Pi$  be a generic irreducible admissible representation of  $G_{n+1}$ , realized in its  $\psi^{-1}$ -Whittaker model  $W(\Pi, \psi^{-1})$ . Then there exists  $W \in W(\Pi, \psi^{-1})$  with the following property. Let  $\pi$  be a generic irreducible admissible representation of  $G_n$ , realized in its  $\psi$ -Whittaker model  $W(\pi, \psi)$ . Let  $\mathfrak{q}$  denote the conductor ideal of  $\pi$ , let  $Q \in F^{\times}$  be a generator of  $\mathfrak{q}^{-1}$ , and set

$$u_Q := I_{n+1} + Q E_{n,n+1} \in G_{n+1}.$$

There exists  $V \in \mathcal{W}(\pi, \psi)$  so that the local Rankin–Selberg integral

$$\int_{N_n \setminus G_n} W(\operatorname{diag}(g, 1)u_Q) V(g) |\det g|^{s - \frac{1}{2}} dg$$

is finite and nonzero for all  $s \in \mathbb{C}$ .

#### 2. Context

Rankin–Selberg local zeta integrals arise as proportionality factors relating global Rankin–Selberg integrals and L-functions. The above result provides test vectors, obtained via pure translates of fixed vectors, that work simultaneously for all representations of the smaller group having some given conductor. Such results are sometimes useful in global applications because they relate problems concerning L-functions (subconvexity, moment asymptotics, ...) to problems concerning automorphic forms (quantitative equidistribution, ...). The n=1 case follows from standard properties of Gauss sums and stationary phase analysis in one variable; it has been applied in, e.g., [7, 6]. For general n, [2] contains a similar result, but with an average over many unipotent translates rather than just one.

#### 3. Proof

We first sketch the argument. The basic idea is to apply the Godement–Jacquet functional to the Whittaker function on the smaller group. This is readily seen to relate the unipotent-shifted Rankin–Selberg integral to an integral involving a translate of the standard congruence subgroup  $K_1(\mathfrak{q}) \leq \operatorname{GL}_n(\mathfrak{o})$ , consisting of matrices whose last row is congruent to  $(0,\ldots,0,1)$  modulo  $\mathfrak{q}$ . We then conclude via newvector theory.

Turning to details, we recall that F is a non-archimedean local field, with ring of integers  $\mathfrak{o}$ . We denote by  $\mathfrak{p}$  the maximal ideal and q the residue field cardinality. We set  $K_r := \mathrm{GL}_r(\mathfrak{o})$  and equip  $G_r$  and  $N_r$  with the Haar measures assigning volume

1

one to K<sup>r</sup> and N<sup>r</sup> ∩ Kr, respectively. As in the theorem statement, we write Π (resp. π) for a generic irreducible representation of Gn+1 (resp. Gn).

We continue to denote by q the conductor ideal of π, defined to be the smallest ideal for which π has a nonzero vector fixed by K1(q). We choose a generator Q for q −1 , so that |Q| = [o : q]. We recall (see [4, 5]) that |Q| (and hence q) may also be characterized in terms of the local ε-factor of π:

$$\varepsilon(\frac{1}{2} + s, \pi, \psi) = |Q|^{-s} \varepsilon(\frac{1}{2}, \pi, \psi). \tag{1}$$

We recall the functional equation of Godement–Jacquet [3, Theorem 3.3].

Lemma 2. Let f be a matrix coefficient of π, and let ϕ ∈ S(Mn(F)). For s ∈ C, the local zeta integral

$$Z(\phi, f, s) := \int_{G_n} \phi(g) f(g) |\det g|^{\frac{n-1}{2} + s} dg,$$
 (2)

converges absolutely for ℜ(s) sufficiently large. It extends to a meromorphic function on the complex plane for which the ratio

$$\frac{Z(\phi, f, s)}{L(s, \pi)}$$

is holomorphic. It satisfies the local functional equation

$$\gamma(s, \pi, \psi) Z(\phi, f, s) = Z(\phi^{\wedge}, f^{\vee}, 1 - s), \tag{3}$$

where

$$\gamma(s, \pi, \psi) = \varepsilon(s, \pi, \psi) \frac{L(1 - s, \tilde{\pi})}{L(s, \pi)},$$

with π˜ the contragredient of π, and where the Fourier transform is defined by

$$f^{\vee}(g) := f(g^{-1}),$$

$$\phi^{\wedge}(x) := \int_{M_n(F)} \phi(y) \psi(\operatorname{trace}(xy)) \, dy,$$

with M<sup>n</sup> the space of n×n matrices and the Haar measure normalized to be self-dual with respect to ψ. Moreover, both of the zeta integrals in (3) converge absolutely provided that, e.g., π is unitary and generic and ℜ(s) = 1/2.

We recall that a matrix coefficient of π is a linear combination of functions of the form f(g) = ℓ(gv), where v ∈ π and ℓ lies in the contragredient of π (i.e., the admissible dual). The conclusions of Lemma 2 remain valid for more general coefficients of π. For instance, suppose more generally that f is of the same form, but with ℓ allowed to be any linear functional on π (not necessarily in the admissible dual). Given ϕ as above, we may choose a compact open subgroup U of G<sup>n</sup> under which ϕ is bi-invariant. The integrals in question do not change if we then replace f by its two-sided average with respect to U, which has the effect of replacing v by its average v <sup>U</sup> ∈ π <sup>U</sup> and ℓ with its projection ℓ <sup>U</sup> to the dual of π <sup>U</sup> , extended by zero on the kernel of the averaging operator π → π <sup>U</sup> . In particular, by specializing to the case that ℓ is a Whittaker functional on π, we see that such identities remain valid when f is a Whittaker function for π.

We denote by S e (F <sup>×</sup>) the space of all Schwartz–Bruhat functions β ∈ S(F ×) such that β(xy) = β(x) whenever |y| = 1, or equivalently, for which β(x) depends

only upon |x|. We note that each  $\beta \in \mathcal{S}^e(F^\times)$  satisfies the Mellin inversion formula

$$\beta(y) = \int_{(\sigma)} \tilde{\beta}(s)|y|^s ds, \qquad \tilde{\beta}(s) := \int_{F^{\times}} \beta(y)|y|^{-s} d^{\times}y. \tag{4}$$

For  $\beta \in \mathcal{S}^e(F^\times)$ , we define the transform  $\beta^{\sharp} := \beta^{\sharp,\pi}$  of  $\beta$  by

$$\beta^{\sharp}(y) := \int_{(\sigma)} \frac{\tilde{\beta}(s)|y|^{-s} ds}{\gamma(\frac{1}{2} + s, \pi, \psi)},$$

initially for  $\sigma$  large enough

**Lemma 3.** Define  $\beta$  via Mellin inversion (4) by

$$\tilde{\beta}(s) := \frac{\varepsilon(\frac{1}{2} + s, \pi, \psi)}{L(\frac{1}{2} + s, \pi)}.$$

Then:

- (1)  $\beta$  is supported on  $\{y: |Q| \leq |y| \leq |Q|q^n\}$  and takes the value  $\varepsilon(\frac{1}{2}, \pi, \psi)$  on  $\{y: |y| = |Q|\}.$  (2)  $\beta^{\sharp}$  is supported on  $\{y: 1 \le |y| \le |q|\}$  and takes the value  $\varepsilon(\frac{1}{2}, \pi, \psi)$  on  $\{y: |y| = |Q|\}.$   $\{y: |y| = 0\}$  and takes the value 1 on  $\{y: |y| = 1\} = \mathfrak{o}^{\times}.$

*Proof.* We appeal to the characterization (1) of |Q|. We note first that  $\beta^{\sharp}$  has Mellin transform

$$\widetilde{\beta}^{\sharp}(s) = \frac{1}{L(\frac{1}{2} + s, \tilde{\pi})}.$$

Since the inverse L-values appearing above are monic polynomials in  $q^{-s}$  of degree at most n, we see by Mellin inversion that  $\beta$  and  $\beta^{\sharp}$  have the claimed properties.  $\square$ 

**Lemma 4.** Assume that  $\pi$  is unitary and generic. We then have the identity of absolutely convergent integrals

$$\int_{G_n} \phi(g) f(g) \beta(\det g) |\det g|^{\frac{n}{2}} dg = \int_{G_n} \phi^{\wedge}(g) f^{\vee}(g) \beta^{\sharp}(\det g) |\det g|^{\frac{n}{2}} dg.$$
 (5)

*Proof.* Starting with the left hand side, we insert the Mellin expansion of  $\beta$ , with  $\sigma = 0$ . The resulting double integral over g and s converges absolutely, so we may swap the order. We recognize the result as the integral  $\int_{(0)} \beta(s) Z(\phi, f, \frac{1}{2} +$ s) ds involving the Godement–Jacquet zeta integral (2). We now apply the local functional equation and expand the result as

$$\int_{(0)} \frac{\tilde{\beta}(s)}{\gamma(\frac{1}{2}+s,\pi,\psi)} \left( \int_{G_n} \phi^{\wedge}(g) f^{\vee}(g) |\det g|^{\frac{n}{2}-s} dg \right) ds.$$

This double integral again converges absolutely, so we may rearrange it to obtain the stated identity. 

For the same reasons as indicated following the statement of Lemma 2, such identities persist for more general coefficients than matrix coefficients, and in particular, when f is a Whittaker function.

Recall that we embed  $G_n \hookrightarrow G_{n+1}$  as the upper-left block. We set

$$W_0(g) := \int_{N_n} 1_{K_n}(xg) \, \psi(x) \, dx, \tag{6}$$

which defines a Whittaker function on  $G_n$  and extends, by the theory of the Kirillov model [1], to an element of  $\mathcal{W}(\Pi, \psi^{-1})$  on  $G_{n+1}$ .

For  $x \in F$  and  $y \in F^{\times}$ , we set

$$d_y := \text{diag}(1, \dots, 1, y) \in G_n \hookrightarrow G_{n+1}, \qquad u_x := I_{n+1} + x E_{n, n+1} \in N_{n+1}.$$

We then define

$$t_Q := d_Q^{-1} u_Q = u_1 d_Q^{-1}.$$

**Lemma 5.** There exist  $\beta \in \mathcal{S}^e(F^{\times})$  and  $\phi \in \mathcal{S}(M_n(F))$  so that for all  $g \in G_n$ , we have

$$\int_{N_n} \beta(\det xg)\phi(xg)\psi(x) dx = \varepsilon(\frac{1}{2}, \pi, \psi)W_0(gt_Q)$$
 (7)

and

$$\beta^{\sharp}(\det g)\phi^{\wedge}(g) = |Q|^n 1_{K_1(\mathfrak{g})}(g). \tag{8}$$

Proof. We set

$$\phi_0 := 1_{M_n(\mathfrak{o})},$$

$$\phi(x) := \psi(-x_{nn})\phi_0(xd_Q^{-1}). \tag{9}$$

and take  $\beta$  as in Lemma 3, so that in particular,

$$\beta|_{Q_0} = \varepsilon(\frac{1}{2}, \pi, \psi) 1_{Q_0 \times} \tag{10}$$

and

$$\beta^{\sharp}|_{\mathfrak{g}} = 1_{\mathfrak{g}^{\times}}.\tag{11}$$

We must verify the relations (7) and (8).

We start with (7). Recall from (6) that  $W_0$  is the  $\psi^{-1}$ -Whittaker function  $W_0(g) = \int_{N_n} 1_{K_n}(xg)\psi(x) dx$ . In particular,

$$W_0(gt_Q) = W_0(gu_1d_Q^{-1}) = \psi(-g_{nn})W_0(gd_Q^{-1}). \tag{12}$$

Using this identity, we may rewrite the desired relation (7) as

$$\int_{N_n} \beta(\det(xg))\phi(xg)\psi(x) dx = \varepsilon(\frac{1}{2}, \pi, \psi)\psi(-g_{nn})W_0(gd_Q^{-1}).$$
 (13)

We verify this as follows. First, we see from the definition (9) and the identity  $(xg)_{nn} = g_{nn}$  that for  $x \in N_n$  and  $g \in G_n$ , we have

$$\phi(xg) = \psi(-g_{nn})\phi_0(xgd_Q^{-1}). \tag{14}$$

Next, we have

$$\begin{split} \beta(\det g)\phi_0(gd_Q^{-1}) &= \varepsilon(\tfrac{1}{2},\pi,\psi) \mathbf{1}_{Q\mathfrak{o}^\times}(\det g)\phi_0(gd_Q^{-1}) \\ &= \varepsilon(\tfrac{1}{2},\pi,\psi) \mathbf{1}_{K_n}(gd_Q^{-1}). \end{split}$$

(In the first step, we use that  $\phi_0(gd_Q^{-1})$  is nonzero only if  $\det(g) \in Q\mathfrak{o}$  and apply (10). In the second step, we use that  $1_{K_n}(g) = 1_{\mathfrak{o}^{\times}}(\det g)\phi_0(g)$  and  $\det(d_Q) = Q$ , which gives  $1_{Q\mathfrak{o}^{\times}}(\det g)\phi_0(gd_Q^{-1}) = 1_{K_n}(gd_Q^{-1})$ .) Combining the above identities, we obtain

$$\beta(\det(xg))\phi(xg) = \varepsilon(\frac{1}{2}, \pi, \psi)\psi(-g_{nn})1_{K_n}(xgd_Q^{-1}).$$

Integrating both sides against  $\psi(x) dx$  gives (13), as required.

We verify (8) as follows (here  $E_{ij}$  denotes the elementary matrix):

$$\beta^{\sharp}(\det g)\phi^{\wedge}(g) = 1_{\mathfrak{o}^{\times}}(\det g)\phi^{\wedge}(g)$$

$$= 1_{\mathfrak{o}^{\times}}(\det g)|Q|^{n}\phi_{0}^{\wedge}(d_{Q}(g - E_{nn}))$$

$$= |Q|^{n}1_{\mathfrak{o}^{\times}}(\det g)1_{M_{n}(\mathfrak{o})}(d_{Q}(g - E_{nn}))$$

$$= |Q|^{n}1_{K_{1}(\mathfrak{g})}(g).$$

Here, for the first step, we observed that  $\phi^{\wedge}(x)$  is nonzero only if  $x \in E_{nn} + d_Q^{-1}M_n(\mathfrak{o}) \subseteq M_n(\mathfrak{o})$ , so that, in particular, det  $x \in \mathfrak{o}$ ; we then applied (11). For the second step, we applied the general Fourier analytic calculation

$$\phi^{\wedge}(x) = |Q|^n \phi_0^{\wedge}(d_O(x - E_{nn})). \tag{15}$$

For the third, we applied the Fourier self-duality  $\phi_0^{\wedge} = \phi_0 = 1_{M_n(\mathfrak{o})}$ . For the final step, we use that  $K_1(\mathfrak{q})$  consists of all  $x \in M_n(F)$  for which  $d_Q(x - E_{nn}) \in M_n(\mathfrak{o})$  and  $\det x \in \mathfrak{o}^{\times}$ .

For  $W \in \mathcal{W}(\Pi, \psi^{-1})$ ,  $V \in \mathcal{W}(\pi, \psi)$ , and  $s \in \mathbb{C}$ , we define the Rankin–Selberg integral

$$\ell_{RS}(s, W, V) := \int_{N_n \setminus G_n} W(\operatorname{diag}(g, 1)) V(g) |\det g|^{s - \frac{1}{2}} dg.$$
 (16)

The following result verifies Theorem 1 in a more precise form.

**Proposition 6.** Let  $W_0 \in \mathcal{W}(\Pi, \psi^{-1})$  be such that for all  $g \in G_n$ , we have

$$W_0(g) = \int_{N_n} 1_{K_n}(xg)\psi(x) dx.$$

Let  $V \in \mathcal{W}(\pi, \psi)$  denote the normalized newvector (i.e., the unique  $K_1(\mathfrak{q})$ -invariant vector for which V(1) = 1, see [4, 5]). Then for all  $s \in \mathbb{C}$ , we have

$$\ell_{\rm RS}(s, u_O W_0, d_O V) = c|Q|^{-\frac{n}{2}},\tag{17}$$

where

$$c := \varepsilon(\frac{1}{2}, \pi, \psi)^{-1} |Q|^n \operatorname{vol}(K_1(\mathfrak{q})) \approx 1.$$
(18)

*Proof.* We note first that, by a change of variables, we have the homogeneity property

$$\ell_{RS}(s, u_Q W_0, d_Q V) = |Q|^{-\left(s - \frac{1}{2}\right)} \ell_{RS}(s, t_Q W_0, V). \tag{19}$$

In view of this, the desired identity (17) is equivalent to

$$\ell_{RS}(s, t_Q W_0, V) = c|Q|^{s - \frac{n+1}{2}}.$$
(20)

Next, since  $W_0$  is supported on  $\det^{-1}(\mathfrak{o}^{\times})$ , we see that the translate  $t_QW_0$  is supported on  $\det^{-1}(Q\mathfrak{o}^{\times})$ , so the left hand side of (20) is a constant multiple of  $|Q|^s$ . For this reason, it suffices to verify (20) for (say)  $s = \frac{n+1}{2}$ , where our task is to check that  $\ell_{\text{RS}}(\frac{n+1}{2}, t_QW_0, V) = c$ . Inserting definitions and unfolding, we obtain,

with 
$$f(g) := V(g)$$
,

$$\begin{split} \varepsilon(\frac{1}{2},\pi,\psi)\ell_{\mathrm{RS}}(\frac{n+1}{2},t_{Q}W_{0},V) &\stackrel{(16)}{=} \varepsilon(\frac{1}{2},\pi,\psi) \int_{N_{n}\backslash G_{n}} W_{0}(gt_{Q})V(g)|\mathrm{det}(g)|^{\frac{n}{2}}\,dg \\ &\stackrel{(7)}{=} \int_{G_{n}} \phi(g)f(g)\beta(\det g)|\det g|^{n/2}\,dg \\ &\stackrel{(5)}{=} \int_{G_{n}} \phi^{\wedge}(g)f^{\vee}(g)\beta^{\sharp}(\det g)|\det g|^{n/2}\,dg \\ &\stackrel{(8)}{=} |Q|^{n} \int_{K_{1}(\mathfrak{q})} V(g^{-1})|\det g|^{n/2}\,dg \\ &= |Q|^{n}\operatorname{vol}(K_{1}(\mathfrak{q})), \end{split}$$

where in the final step, we use the K1(q)-invariance of V , the normalization V (1) = 1, and the fact that |det g| = 1 on K1(q). Thus (20) holds. □

#### References

- [1] Joseph N. Bernstein. P-invariant distributions on GL(N) and the classification of unitary representations of GL(N) (non-Archimedean case). In Lie group representations, II (College Park, Md., 1982/1983), volume 1041 of Lecture Notes in Math., pages 50–102. Springer, Berlin, 1984.
- [2] Andrew R. Booker, M. Krishnamurthy, and Min Lee. Test vectors for Rankin-Selberg Lfunctions. J. Number Theory, 209:37–48, 2020.
- [3] Roger Godement and Herv´e Jacquet. Zeta functions of simple algebras. Lecture Notes in Mathematics, Vol. 260. Springer-Verlag, Berlin, 1972.
- [4] H. Jacquet, I. I. Piatetski-Shapiro, and J. Shalika. Conducteur des repr´esentations du groupe lin´eaire. Math. Ann., 256(2):199–214, 1981.
- [5] Nadir Matringe. Essential Whittaker functions for GL(n). Doc. Math., 18:1191–1214, 2013.
- [6] Philippe Michel and Akshay Venkatesh. The subconvexity problem for GL2. Publ. Math. Inst. Hautes Etudes Sci. ´ , (111):171–271, 2010.
- [7] Peter Sarnak. Fourth moments of Gr¨ossencharakteren zeta functions. Comm. Pure Appl. Math., 38(2):167–178, 1985.

