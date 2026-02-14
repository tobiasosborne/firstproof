# Official Solution: Problem 9: Algebraic Relations on Quadrilinear Determinantal Tensors (Joe Kileel)

> Source: *First Proof Solutions and Comments*, Abouzaid, Blumberg, Hairer, Kileel, Kolda, Nelson, Spielman, Srivastava, Ward, Weinberger, Williams. February 14, 2026.

---

## Authors' Commentary on AI-Generated Solutions


The best LLM answer found during testing was NoInternet-040226. This is an essentially correct answer. It constructs the same algebraic relations as in my own answer, namely the various 5 × 5 minors of the four 3*n* × 27*n* <sup>3</sup> flattenings of the block tensor assembling together the *Q*(*αβγδ*) . The proof by the LLM that the algebraic relations satisfy the desired properties differs from my own argument. The LLM considers a torus action on an appropriate Grassmannian, argues the stabilizer of a generic point is 1-dimensional, and uses this to show separability of *λ* in a somewhat fidgety way. By contrast, I directly constrain *λ* by considering certain selected algebraic relations. Some other LLM answers produced during testing were incorrect, and claimed that no algebraic relations exist that satisfy the desired properties. Those answers seemed to get confused about the question setup midway through. My question is closely related to a work I published with Miao and Lerman in 2024 ([https://proceedings.neurips.cc/paper\\_files/paper/2024/hash/](https://proceedings.neurips.cc/paper_files/paper/2024/hash/80cddcdd52c84d19b8b4a27a8e8c17d8-Abstract-Conference.html) [80cddcdd52c84d19b8b4a27a8e8c17d8-Abstract-Conference.html](https://proceedings.neurips.cc/paper_files/paper/2024/hash/80cddcdd52c84d19b8b4a27a8e8c17d8-Abstract-Conference.html)). Indeed, it is a fourthorder variant of Theorem 2 in that paper which concerns the third-order case. Therefore, if LLMs locate and understand that paper they would have a warm-start for this question.

---

## Official Solution

# Question

Let  $n \geq 5$ . Let  $A^{(1)}, \ldots, A^{(n)} \in \mathbb{R}^{3 \times 4}$  be Zariski-generic. For  $\alpha, \beta, \gamma, \delta \in [n]$ , construct  $Q^{(\alpha\beta\gamma\delta)} \in \mathbb{R}^{3 \times 3 \times 3 \times 3}$  so that its  $(i,j,k,\ell)$  entry for  $1 \leq i,j,k,\ell \leq 3$  is given by  $Q^{(\alpha\beta\gamma\delta)}_{ijk\ell} = \det[A^{(\alpha)}(i,:);A^{(\beta)}(j,:);A^{(\gamma)}(k,:);A^{(\delta)}(\ell,:)]$ . Here A(i,:) denotes the ith row of a matrix A, and semicolon denotes vertical concatenation. We are interested in algebraic relations on the set of tensors  $\{Q^{(\alpha\beta\gamma\delta)}: \alpha, \beta, \gamma, \delta \in [n]\}$ .

More precisely, does there exist a polynomial map  $\mathbf{F}: \mathbb{R}^{81n^4} \to \mathbb{R}^N$  that satisfies the following three properties?

- The map **F** does not depend on  $A^{(1)}, \dots A^{(n)}$ .
- The degrees of the coordinate functions of F do not depend on n.
- Let  $\lambda \in \mathbb{R}^{n \times n \times n \times n}$  satisfy  $\lambda_{\alpha\beta\gamma\delta} \neq 0$  for precisely  $\alpha, \beta, \gamma, \delta \in [n]$  that are not identical. Then  $\mathbf{F}(\lambda_{\alpha\beta\gamma\delta}Q^{(\alpha\beta\gamma\delta)}: \alpha, \beta, \gamma, \delta \in [n]) = 0$  holds if and only if there exist  $u, v, w, x \in (\mathbb{R}^*)^n$  such that  $\lambda_{\alpha\beta\gamma\delta} = u_{\alpha}v_{\beta}w_{\gamma}x_{\delta}$  for all  $\alpha, \beta, \gamma, \delta \in [n]$  that are not identical.

JK acknowledges support from NSF DMS 2309782, NSF CISE-IIS 2312746, DE SC0025312, and the Sloan Foundation.

# Answer (from work by Daniel Miao, Gilad Lerman, Joe Kileel)

Yes, such algebraic relations do exist. Assemble the various tensors  $\{Q^{(\alpha\beta\gamma\delta)}:\alpha,\beta,\gamma,\delta\in[n]\}$  into one tensor  $\mathbf{Q}\in\mathbb{R}^{3n\times3n\times3n}$ , thought of as an  $n\times n\times n\times n$  block tensor where the  $(\alpha,\beta,\gamma,\delta)$ -block is  $Q^{(\alpha\beta\gamma\delta)}\in\mathbb{R}^{3\times3\times3\times3}$ . Let  $\mathbf{F}$  be the polynomial map sending  $\{Q^{(\alpha\beta\gamma\delta)}:\alpha,\beta,\gamma,\delta\in[n]\}$  to the  $5\times 5$  minors of the four  $3n\times27n^3$  matrix flattenings of  $\mathbf{Q}$ . We will prove that  $\mathbf{F}$  satisfies the desired properties.

A key point is to discover the following algebraic identity.

**Lemma 1.** Consider  $\mathbf{Q} \in \mathbb{R}^{3n \times 3n \times 3n}$  as above. It admits a Tucker tensor decomposition

$$\mathbf{Q} = \mathcal{C} \times_1 \mathbf{A} \times_2 \mathbf{A} \times_3 \mathbf{A} \times_4 \mathbf{A}, \tag{1}$$

for  $C \in \mathbb{R}^{4 \times 4 \times 4 \times 4}$  and  $\mathbf{A} \in \mathbb{R}^{3n \times 4}$ . Explicitly, we can take

$$C_{abcd} = \begin{cases} sgn(abcd) & if \ a, b, c, d \in [4] \ are \ distinct \\ 0 & otherwise, \end{cases}$$

where sgn is parity of a permutation, and **A** to be the vertical concatenation  $[A^{(1)}; \ldots; A^{(n)}]$ .

*Proof.* Let  $[n] \times [3]$  stand for the indices of **Q** in each mode and for the row indices of **A**. By definition of Tucker product, for all  $(\alpha, i), (\beta, j), (\gamma, k), (\delta, \ell) \in [n] \times [3]$  we have

$$(\mathcal{C} \times_1 \mathbf{A} \times_2 \mathbf{A} \times_3 \mathbf{A} \times_4 \mathbf{A})_{(\alpha,i),(\beta,j),(\gamma,k),(\delta,\ell)} = \sum_{a,b,c,d \in [4]} \mathcal{C}_{abcd} \mathbf{A}_{(\alpha,i),a} \mathbf{A}_{(\beta,j),b} \mathbf{A}_{(\gamma,k),c} \mathbf{A}_{(\delta,\ell),d}$$

$$=\sum_{a,b,c,d\in[4]\text{ distinct}}\operatorname{sgn}(abcd)A_{ia}^{(\alpha)}A_{jb}^{(\beta)}A_{kc}^{(\alpha)}A_{\ell d}^{(\alpha)}=\det\Big[A^{(\alpha)}(i,:);A^{(\beta)}(j,:);A^{(\gamma)}(k,:);A^{(\delta)}(\ell,:)\Big]$$

$$=Q_{ijk\ell}^{(\alpha\beta\gamma\delta)} = \mathbf{Q}_{(\alpha,i),(\beta,j),(\gamma,k),(\delta,\ell)}.$$

The lemma explains why **F** captures algebraic relations between the tensors  $\{Q^{(\alpha\beta\gamma\delta)}: \alpha, \beta, \gamma, \delta \in [n]\}$ . Indeed, the block tensor **Q** has multilinear rank bounded by (4,4,4,4) due to the Tucker decomposition in (1). Therefore, all  $5 \times 5$  minors in **F** vanish.

Below we break up the proof of the third property into two directions. The other properties are clear. Throughout the proof, for  $\lambda \in \mathbb{R}^{n \times n \times n \times n}$  we let  $\lambda \odot_b \mathbf{Q} \in \mathbb{R}^{3n \times 3n \times 3n \times 3n}$  denote blockwise scalar multiplication, i.e., the  $(\alpha, \beta, \gamma, \delta)$ -block of  $\lambda \odot_b \mathbf{Q}$  is  $\lambda_{\alpha\beta\gamma\delta}Q^{(\alpha\beta\gamma\delta)} \in \mathbb{R}^{3\times3\times3\times3}$ . Roughly speaking, we need to show that a blockwise scaling of  $\mathbf{Q}$  preserves multilinear rank if and only if the scaling is a rank-1 tensor off the diagonal.

#### "If" Direction

This follows easily from Lemma 1. Assume  $\lambda \in \mathbb{R}^{n \times n \times n \times n}$  agrees off-diagonal with  $u \otimes v \otimes w \otimes x$  for  $u, v, w, x \in (\mathbb{R}^*)^n$  and is 0 on the diagonal. Then

$$\lambda \odot_b \mathbf{Q} = (u \otimes v \otimes w \otimes x) \odot_b \mathbf{Q},$$

because the diagonal blocks of  $\mathbf{Q}$  vanish. That is,  $Q^{(\alpha\alpha\alpha\alpha)}=0$  since each entry of  $Q^{(\alpha\alpha\alpha\alpha)}$  is the determinant of a matrix with a repeated row. Note that blockwise scalar product with a rank-1 tensor with nonzero entries is equivalent to Tucker product with invertible matrices:

$$(u \otimes v \otimes w \otimes w) \odot_b \mathbf{Q} = \mathbf{Q} \times_1 D_u \times_2 D_v \times_3 D_w \times_4 D_x.$$

Here  $D_u \in \mathbb{R}^{3n \times 3n}$  is the diagonal matrix triplicating the entries of u and likewise for  $D_v, D_w, D_x$ . Thus  $\lambda \odot_b \mathbf{Q}$  and  $\mathbf{Q}$  have the same multilinear rank, and from the lemma  $\mathbf{F}(\lambda_{\alpha\beta\gamma\delta}Q^{(\alpha\beta\gamma\delta)}: \alpha, \beta, \gamma, \delta \in [n]) = 0$ .

## "Only If" Direction

The converse takes more work. Let  $\lambda \in \mathbb{R}^{n \times n \times n}$  have nonzero entries precisely off the diagonal and assume  $\mathbf{F}(\lambda_{\alpha\beta\gamma\delta}Q^{(\alpha\beta\gamma\delta)}:\alpha,\beta,\gamma,\delta\in[n])=0$ . We further assume  $\lambda_{\alpha 111}=\lambda_{1\beta 11}=\lambda_{11\gamma 1}=\lambda_{111\delta}=1$  for all  $\alpha,\beta,\gamma,\delta\in\{2,\ldots,n\}$ . We reduce to this case by replacing  $\lambda$  by its entrywise product with  $\bar{u}\otimes\bar{v}\otimes\bar{w}\otimes\bar{x}$ , where

$$\bar{u}_{\alpha} = \begin{cases} 1 & \text{for } \alpha = 1\\ \lambda_{\alpha 1 1 1}^{-1} & \text{for } \alpha \in \{2, \dots, n\}, \end{cases}$$

and  $\bar{v}, \bar{w}, \bar{x}$  are defined similarly using the second, third and fourth modes respectively. The replacement preserves the multilinear rank of  $\lambda \odot_b \mathbf{Q}$  and whether or not  $\lambda$  agrees off-diagonal with a rank-1 tensor. Hence it is without loss of generality.

Through some explicit calculations, we will prove there exists  $c \in \mathbb{R}^*$  such that

- $\lambda_{\alpha\beta\gamma\delta} = c$  if exactly two of  $\alpha, \beta, \gamma, \delta$  equal 1
- $\lambda_{\alpha\beta\gamma\delta} = c^2$  if exactly one of  $\alpha, \beta, \gamma, \delta$  equals 1
- $\lambda_{\alpha\beta\gamma\delta} = c^3$  if none of  $\alpha, \beta, \gamma, \delta$  equal 1 and  $\alpha, \beta, \gamma, \delta$  are not identical.

This will establish the "only if" direction, as setting u = v = w = (1, c, ..., c) and  $x = (\frac{1}{c}, 1, ..., 1)$  gives  $\lambda_{\alpha\beta\gamma\delta} = u_{\alpha}v_{\beta}w_{\gamma}x_{\delta}$  whenever  $\alpha, \beta, \gamma, \delta$  are not identical. Our proof strategy is to examine appropriate coordinates of  $\mathbf{F}(\lambda_{\alpha\beta\gamma\delta}Q^{(\alpha\beta\gamma\delta)}:\alpha,\beta,\gamma,\delta\in[n])=0$  in order to constrain  $\lambda$ . Equivalently, we will consider the vanishing of the determinants of certain well-chosen  $5\times 5$  submatrices of the flattenings of  $\lambda \odot_b \mathbf{Q}$ . Write  $\mathbf{Q}_{(1)}$  and  $(\lambda \odot_b \mathbf{Q})_{(1)}$  for mode-1 flattenings in  $\mathbb{R}^{3n\times 27n^3}$ . Rows correspond to the first tensor mode and are indexed by  $(\alpha,i)\in[n]\times[3]$ , while columns correspond to the other modes and are indexed by  $((\beta,j),(\gamma,k),(\delta,\ell))\in([n]\times[3])^3$ .

Step 1: The first submatrix of  $(\lambda \odot_b \mathbf{Q})_{(1)}$  we consider has column indices  $((\alpha, 1), (1, 3), (1, 2)), ((1, 2), (\beta, 2), (1, 1)), ((1, 2), (\beta, 3), (1, 1)), ((1, 3), (\beta, 3), (1, 2)), ((1, 1), (\beta, 1), (1, 3))$  and row indices  $(1, 1), (1, 2), (1, 3), (\alpha, 1), (\alpha, 2)$ , where  $\alpha, \beta \in \{2, \ldots, n\}$ . Explicitly, the submatrix is

$$\begin{bmatrix} Q_{1132}^{(1\alpha 11)} & Q_{1221}^{(11\beta 1)} & Q_{1231}^{(11\beta 1)} & Q_{1332}^{(11\beta 1)} & Q_{1113}^{(11\beta 1)} \\ Q_{2132}^{(1\alpha 11)} & Q_{2221}^{(11\beta 1)} & Q_{2231}^{(11\beta 1)} & Q_{2332}^{(11\beta 1)} & Q_{2113}^{(11\beta 1)} \\ Q_{3132}^{(1\alpha 11)} & Q_{3221}^{(11\beta 1)} & Q_{3231}^{(11\beta 1)} & Q_{3332}^{(11\beta 1)} & Q_{3113}^{(11\beta 1)} \\ \lambda_{\alpha\alpha 11}Q_{1132}^{(\alpha\alpha 11)} & \lambda_{\alpha 1\beta 1}Q_{1221}^{(\alpha 1\beta 1)} & \lambda_{\alpha 1\beta 1}Q_{1231}^{(\alpha 1\beta 1)} & \lambda_{\alpha 1\beta 1}Q_{1332}^{(\alpha 1\beta 1)} & \lambda_{\alpha 1\beta 1}Q_{1113}^{(\alpha 1\beta 1)} \\ \lambda_{\alpha\alpha 11}Q_{2132}^{(\alpha\alpha 11)} & \lambda_{\alpha 1\beta 1}Q_{2221}^{(\alpha 1\beta 1)} & \lambda_{\alpha 1\beta 1}Q_{2231}^{(\alpha 1\beta 1)} & \lambda_{\alpha 1\beta 1}Q_{2332}^{(\alpha 1\beta 1)} & \lambda_{\alpha 1\beta 1}Q_{2113}^{(\alpha 1\beta 1)} \end{bmatrix}$$

which we abbreviate as

$$\begin{bmatrix} * & * & * & * & * & * \\ * & * & * & * &$$

with asterisk denoting the corresponding entry in  $\mathbf{Q}_{(1)}$ . We view the determinant of (2) as a polynomial with respect to  $\lambda$ . It has degree  $\leq 2$  in the variables  $\lambda_{\alpha\alpha 11}$ ,  $\lambda_{\alpha 1\beta 1}$ . Observe that if  $\lambda_{\alpha 1\beta 1} = 0$ , the bottom two rows of the matrix are linearly independent. Also if  $\lambda_{\alpha 1\beta 1} - \lambda_{\alpha\alpha 11} = 0$ , then (2) equals a  $5 \times 5$  submatrix of  $\mathbf{Q}_{(1)}$  with rows operations performed; therefore (2) is rank-deficient. It follows that the determinant of (2) takes the form

$$s\lambda_{\alpha 1\beta 1}(\lambda_{\alpha 1\beta 1} - \lambda_{\alpha \alpha 11}).$$

Here the scale  $s = s(A^{(1)}, A^{(\alpha)}, A^{(\beta)})$  is a polynomial in the A-matrices. Due to polynomiality, s is nonzero Zariski-generically if we can exhibit a *single* instance of matrices  $A^{(1)}, A^{(\alpha)}, A^{(\beta)}$  where the determinant of (2) does not vanish identically for all  $\lambda_{\alpha 1\beta 1}, \lambda_{\alpha \alpha 11}$ . Furthermore, we just need an instance with  $\alpha = \beta$ , as this corresponds to a specialization of the case  $\alpha \neq \beta$ . Computational verification with a random numerical instance of  $A^{(1)}, A^{(\alpha)}$  proves the non-vanishing (see attached code). Recalling the standing assumptions, we deduce  $\lambda_{\alpha 1\beta 1} = \lambda_{\alpha \alpha 11}$ .

We apply the same argument to modewise permutations of  $\lambda \odot_b \mathbf{Q}$  and  $\mathbf{Q}$ , and obtain

$$\lambda_{\pi(\alpha_1\beta_1)} = \lambda_{\pi(\alpha_1\beta_1)}$$
 for all  $\alpha, \beta \in \{2, \dots, n\}$  and permutations  $\pi$ .

The argument goes through as  $\pi \cdot \mathbf{Q}$  and  $\pi \cdot (\lambda \odot_b \mathbf{Q})$  have multilinear ranks bounded by (4, 4, 4, 4) and  $\pi \cdot \mathbf{Q} = \operatorname{sgn}(\pi)\mathbf{Q}$ . So (2) looks the same but with indices permuted and possibly a sign flip.

We now see that  $\lambda$ -entries with two 1-indices agree. Indeed, taking  $\alpha = \beta$  above gives  $\lambda_{\pi_1(\alpha 1 \alpha 1)} = \lambda_{\pi_2(\alpha \alpha 1 1)}$  for all  $\pi_1$  and  $\pi_2$  that fix  $(\alpha \alpha 1 1)$  and  $(\alpha 1 \alpha 1)$  respectively. So,  $\lambda_{\alpha \alpha 1 1} = \lambda_{\pi(\alpha \alpha 1 1)}$  for all  $\pi$ . Taking  $\alpha \neq \beta$  gives  $\lambda_{\alpha \alpha 1 1} = \lambda_{\pi(\alpha 1 \beta 1)} = \lambda_{\beta \beta 1 1}$  for all  $\pi$ . Together, there exists  $c \in \mathbb{R}^*$  such that  $c = \lambda_{\pi(\alpha \beta 1 1)}$  for all  $\alpha, \beta \in \{2, \ldots, n\}$  and permutations  $\pi$ .

Step 2: Next we consider the submatrix of  $(\lambda \odot_b \mathbf{Q})_{(1)}$  with column indices  $((\beta, 1), (\gamma, 3), (1, 2)), ((1, 2), (\beta, 2), (1, 1)), ((1, 2), (\beta, 3), (1, 1)), ((1, 3), (\beta, 3), (1, 2)), ((1, 1), (\beta, 1), (1, 3))$  and row indices  $(1, 1), (1, 2), (1, 3), (\alpha, 1), (\alpha, 2)$ , where  $\alpha, \beta, \gamma \in \{2, \ldots, n\}$ . It looks like

$$\begin{bmatrix} c* & * & * & * & * \\ c* & * & * & * & * \\ c* & * & * & * & * \\ \lambda_{\alpha\beta\gamma1}* & c* & c* & c* & c* \\ \lambda_{\alpha\beta\gamma1}* & c* & c* & c* & c* \end{bmatrix},$$
(3)

where asterisks denote corresponding entries in  $\mathbf{Q}_{(1)}$ . As a polynomial in c and  $\lambda_{\alpha\beta\gamma1}$ , the determinant of (3) is a scalar multiple of  $c(c^2 - \lambda_{\alpha\beta\gamma1})$ . This is because the polynomial has degree  $\leq 3$ , if c = 0 then the bottom two rows of (3) are linearly dependent, and if  $c^2 = \lambda_{\alpha\beta\gamma1}$  then (3) is a  $5 \times 5$  submatrix of  $\mathbf{Q}_{(1)}$  with row and column operations performed. The scale is a polynomial in  $A^{(1)}, A^{(\alpha)}, A^{(\beta)}, A^{(\gamma)}$ . It is Zariski-generically nonzero if we exhibit one instance of A-matrices such that the determinant of (2) does not vanish for all  $c, \lambda_{\alpha\beta\gamma1}$ . Further, it suffices to find an instance where  $\alpha = \beta = \gamma$ , as all other cases specialize to this. Computational verification with a random numerical instance of  $A^{(1)}, A^{(\alpha)}$  proves the non-vanishing. It follows that  $c^2 = \lambda_{\alpha\beta\gamma1}$ . Appealing to symmetry like before,  $c^2 = \lambda_{\pi(\alpha\beta\gamma1)}$  for all  $\alpha, \beta, \gamma \in \{2, \ldots, n\}$  and permutations  $\pi$ . Summarizing, all  $\lambda$ -entries with a single 1-index equal  $c^2$ .

Step 3: Consider the submatrix of  $(\lambda \odot \mathbf{Q})_{(1)}$  with columns  $((\beta, 1), (\gamma, 3), (\delta, 2)), ((1, 2), (\alpha, 2), (1, 1)), ((1, 2), (\alpha, 3), (1, 1)), ((1, 3), (\alpha, 3), (1, 2)), ((1, 1), (\alpha, 1), (1, 3))$  and rows  $(1, 1), (1, 2), (1, 3), (\alpha, 1), (\alpha, 2)$ , where  $\alpha, \beta, \gamma, \delta \in \{2, \ldots, n\}$  and  $\alpha, \delta$  are distinct. The submatrix looks like

$$\begin{bmatrix} c^{2} * & * & * & * & * \\ c^{2} * & * & * & * & * \\ c^{2} * & * & * & * & * \\ \lambda_{\alpha\beta\gamma\delta} * & c * & c * & c * & c * \\ \lambda_{\alpha\beta\gamma\delta} * & c * & c * & c * & c * \end{bmatrix}. \tag{4}$$

The determinant of (4) is  $c(c^3 - \lambda_{\alpha\beta\gamma\delta})$  multiplied by a polynomial in  $A^{(1)}, A^{(\alpha)}, A^{(\beta)}, A^{(\gamma)}, A^{(\delta)}$ . The most specialized case is  $\alpha = \beta = \gamma$ . Computer verification with a random numerical instance proves the polynomial is not identically zero. We deduce that  $c^3 = \lambda_{\alpha\beta\gamma\delta}$ . By symmetry,  $c^3 = \lambda_{\pi(\alpha\beta\gamma\delta)}$  for all  $\alpha, \beta, \gamma, \delta \in \{2, ..., n\}$  with  $\alpha, \delta$  distinct and all permutations  $\pi$ . In other words,  $\lambda$ -entries with no 1-indices and non-identical indices equal  $c^3$ .

Steps 1, 2 and 3 show that  $\lambda$  takes the announced form. So,  $\lambda$  is rank-1 off the diagonal. This finishes the "only if" direction. Overall, we have proven that the  $5 \times 5$  minors of the  $3n \times 27n^3$  flattenings of  $\mathbf{Q}$  give algebraic relations on  $\{Q^{(\alpha\beta\gamma\delta)}: \alpha, \beta, \gamma, \delta \in [n]\}$  with the desired properties.
