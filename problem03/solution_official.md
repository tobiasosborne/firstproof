# Official Solution: Problem 3: Markov Chain with ASEP Polynomial Stationary Distribution (Lauren Williams)

> Source: *First Proof Solutions and Comments*, Abouzaid, Blumberg, Hairer, Kileel, Kolda, Nelson, Spielman, Srivastava, Ward, Weinberger, Williams. February 14, 2026.

---

## Authors' Commentary on AI-Generated Solutions


The best solution that LLM's produced for Question 3 in our internal experiments was to use the Metropolis-Hastings algorithm to produce a Markov chain whose stationary distribution had the desired formula. However, by design, the Metropolis-Hastings algorithm uses the desired formula to define its transition rates. This algorithm can be used to cook up a Markov chain with any desired distribution. Hence this is considered a "trivial" solution to the problem (which specifically asked that the transition probabilities not be described in terms of the interpolation polynomials). Sometimes the LLM's would give a slight variant of the above trivial solution where they would replace the interpolation polynomials by an equivalent formula for them (the signed multiline queue formula of Ben Dali–Williams).

Another common response given by LLM's was to change the problem to a related but different, and already-solved problem, namely, to replace interpolation ASEP and interpolation Macdonald polynomials by ASEP and Macdonald polynomials. In this case the solution to this problem is the *t*-Push TASEP and was given in a paper by Ayyer, Martin, and Williams.

---

## Official Solution

# A PROBABILISTIC INTERPRETATION FOR INTERPOLATION MACDONALD POLYNOMIALS

#### HOUCINE BEN DALI AND LAUREN KIYOMI WILLIAMS

#### 1. The problem

Let  $\lambda = (\lambda_1 > \cdots > \lambda_n \ge 0)$  be a partition with distinct parts. Assume moreover that  $\lambda$  is restricted, in the sense that it has a unique part of size 0 and no part of size 1. Does there exist a nontrivial Markov chain on  $S_n(\lambda)$  whose stationary distribution is given by

$$\frac{F_{\mu}^*(x_1,\ldots,x_n;q=1,t)}{P_{\lambda}^*(x_1,\ldots,x_n;q=1,t)} \text{ for } \mu \in S_n(\lambda)$$

where  $F_{\mu}^{*}(x_{1}, \ldots, x_{n}; q, t)$  and  $P_{\lambda}^{*}(x_{1}, \ldots, x_{n}; q, t)$  are the interpolation ASEP polynomial and interpolation Macdonald polynomial, respectively? If so, prove that the Markov chain you construct has the desired stationary distribution. By "nontrivial" we mean that the transition probabilities of the Markov chain should not be described using the polynomials  $F_{\mu}^{*}(x_{1}, \ldots, x_{n}; q, t)$ .

Date: January 30, 2026.

## 2. The solution

The answer to the question is yes, as we explain below. For  $1 \le k \le n$ , we define

(1) 
$$\mathfrak{p}_k := \frac{t^{-n+1}(1-t)}{x_k - t^{-n+2}} \in \mathbb{Q}(t, x_1, \dots, x_n) \quad \text{and} \quad \mathfrak{q}_k := \frac{(1-t)x_k}{x_k - t^{-n+2}} \in \mathbb{Q}(t, x_1, \dots, x_n).$$

If 0 < t < 1 and  $x_i > t^{-n+1}$  for  $1 \le i \le n$ , then  $\mathfrak{p}_k$  and  $\mathfrak{q}_k$  are probabilities.

**Definition 2.1.** Fix a partition  $\lambda = (\lambda_1 \geqslant \cdots \geqslant \lambda_n)$  with  $\lambda_n = 0$ . The *interpolation* t-Push TASEP with content  $\lambda$  is a Markov chain on  $S_n(\lambda)$ ; we think of its states as configurations of particles on a ring labeled by  $\lambda_1, \ldots, \lambda_n$ , where state  $\eta$  corresponds to having a particle labeled  $\eta_j$  at position j. Moreover, there is a bell attached to each particle. The transitions from  $\eta \in S_n(\lambda)$  are as follows.

(Step 0) The bell at position j rings with probability

$$P_{j} = \frac{\prod_{k < j} \left( x_{k} - \frac{1}{t^{n-2}} \right) \prod_{k > j} \left( x_{k} - \frac{1}{t^{n-1}} \right)}{e_{n-1}^{*}(\boldsymbol{x}; t)},$$

where  $e_{n-1}^*(\boldsymbol{x};t) = \sum_{j=1}^n \prod_{k < j} \left( x_k - \frac{1}{t^{n-2}} \right) \prod_{k > j} \left( x_k - \frac{1}{t^{n-1}} \right)$ .

- (Step 1) The particle at position j, say with label a, is activated, and starts traveling clockwise according to the rules of the t-Push TASEP. That is, suppose there are m "weaker" particles in the system, i.e. particles whose labels are less than a, including vacancies (label 0). Then with probability  $\frac{t^{k-1}}{[m]_t}$  the activated particle will move to the location of the kth of these weaker particles. If this location contains a particle with positive label, then that particle becomes active, and chooses a weaker particle to displace in the same way. The procedure continues until the active particle arrives at a vacancy. At the end of this step, position j is vacant, and we regard this vacancy as a particle labeled a := 0.
- (Step 2) The particle labeled a := 0 now goes to position 1 and starts traveling clockwise. When it gets to site k for  $1 \le k \le j-1$  containing a particle with label  $b \ge 0$ , it skips over that site with probability  $1 \mathfrak{p}_k$  if  $b \ge a$ , and  $1 \mathfrak{q}_k$  if b < a; otherwise it settles at that site, activating/ displacing the site's particle. Once it activates a new particle, the old particle settles at site k and the new active particle continues to travel clockwise towards position j, activating a new particle according to the rule above. The active particle stops once it displaces/activates another particle or arrives at position j, in which case it settles in position j.

We denote the resulting configuration by  $\nu$  and the transition probability by  $\mathbb{P}(\eta, \nu)$ .

Moreover, we let  $\mathbb{P}_{\lambda,j}^{(1)} = \mathbb{P}_{j}^{(1)}$  and  $\mathbb{P}_{\lambda,j}^{(2)} = \mathbb{P}_{j}^{(2)}$  denote the transition probabilities associated with (Step 1) and (Step 2), respectively. We then have, for  $\mu, \nu \in S_n(\lambda)$ ,

$$\mathbb{P}(\mu,\nu) = \sum_{1\leqslant j\leqslant n} P_j \sum_{\rho\in S_n(\lambda): \rho_j=0} \mathbb{P}_j^{(1)}(\mu,\rho) \mathbb{P}_j^{(2)}(\rho,\nu).$$

**Theorem 2.2.** In the interpolation t-Push TASEP with content  $\lambda = (\lambda_1, ..., \lambda_n)$  and parameters  $\mathbf{x} = (x_1, ..., x_n)$  and t, the stationary probability of  $\mu \in S_n(\lambda)$  is given by

$$\pi_{\lambda}^*(\mu) = \frac{F_{\mu}^*(\boldsymbol{x}; 1, t)}{P_{\lambda}^*(\boldsymbol{x}; 1, t)}.$$

#### 3. The proof

Recall the notion of classical two-line queues from [CMW22] and signed two-line queues from [BDW25] together with their weight functions. (Here we specialize q = 1.)

Let  $\mathcal{Q}^{\eta}_{\kappa}$  denote the set of classical two-line queues with top row  $\eta = (\eta_1, \dots, \eta_n)$  and bottom row  $\kappa = (\kappa_1, \dots, \kappa_n)$ , and let  $a^{\eta}_{\kappa}$  denote the weight generating function of  $\mathcal{Q}^{\eta}_{\kappa}$ .

(2) 
$$a_{\kappa}^{\eta} = a_{\kappa}^{\eta}(t) := \sum_{Q \in \mathcal{Q}_{\kappa}^{\eta}} \operatorname{wt}_{\operatorname{pair}}(Q).$$

Let  $\mathcal{G}^{\alpha}_{\mu}$  denote the set of signed two-line queues with top row  $\alpha = (\alpha_1, \ldots, \alpha_n)$  and bottom row  $\mu = (\mu_1, \ldots, \mu_n)$ , and let  $b^{\alpha}_{\mu}$  denote the weight generating function of  $\mathcal{G}^{\alpha}_{\mu}$ .

(3) 
$$b_{\mu}^{\alpha} = b_{\mu}^{\alpha}(t) := \sum_{Q \in \mathcal{G}_{\mu}^{\alpha}} \operatorname{wt}_{\operatorname{pair}}(Q).$$

Let  $\operatorname{wt}(Q) := \operatorname{wt}_{\operatorname{pair}}(Q) \operatorname{wt}_{\operatorname{ball}}(Q)$  be the product of the pair weight and the ball weight. We obtain

(4) 
$$\operatorname{wt}_{\alpha} b_{\mu}^{\alpha} = \sum_{Q \in \mathcal{G}_{\mu}^{\alpha}} \operatorname{wt}(Q), \quad \text{where} \quad \operatorname{wt}_{\alpha} := \prod_{k: \alpha_{k} > 0} x_{k} \prod_{k: \alpha_{k} < 0} \frac{-1}{t^{n-1}}.$$

**Definition 3.1.** Given a signed two-line queue  $Q \in \mathcal{G}^{\alpha}_{\mu}$ , we associate to it an *unsigned* version  $\bar{Q}$  obtained by forgetting the signs of the balls in the top row. The composition we read in the bottom row (respectively the top row) of  $\bar{Q}$  is  $\mu$  (respectively  $\|\alpha\|$ ), where

$$\|\alpha\| = (|\alpha_1|, \ldots, |\alpha_n|).$$

We then define  $\bar{\mathcal{G}}^{\kappa}_{\mu}$  as the set of paired ball systems obtained by applying this operation on  $Q \in \mathcal{G}^{\alpha}_{\mu}$ , where  $\alpha \in \mathbb{Z}^n$  satisfying  $\|\alpha\| = \kappa$ .

This leads us to define the following weights. Fix  $\bar{Q} \in \bar{\mathcal{G}}_{\mu}^{\kappa}$ :

• A nontrivial pairing p in  $\bar{Q}$  has the weight

(5) 
$$\operatorname{wt}(p) = (1 - t)t^{\operatorname{skip}(p)}.$$

• Let B be a ball labeled a > 0 in column k and such that the ball below is labeled b (If B has a vacancy below it, we take b = 0.) We define the weight of B by:

(6) 
$$\operatorname{wt}(B) := \begin{cases} x_k - \frac{1}{t^{n-1}} & \text{if } b = a, \\ x_k & \text{if } b > a, \\ \frac{1}{t^{n-1}} & \text{if } b < a. \end{cases}$$

The weight of  $\bar{Q}$  is defined by

$$\operatorname{wt}(\bar{Q}) := \prod_{B \text{ in the top row}} \operatorname{wt}(B) \prod_{p \text{ nontrivial pairing}} \operatorname{wt}(p).$$

We then have the following lemma.

**Lemma 3.2.** Fix a partition  $\lambda$  with distinct parts and two compositions  $\kappa, \mu \in S_n(\lambda)$ . Let  $\overline{Q} \in \overline{\mathcal{G}}^{\kappa}_{\mu}$ . Then

$$\operatorname{wt}(\bar{Q}) = \sum_{Q} \operatorname{wt}(Q),$$

where the sum is taken over all signed two-line queues Q from which  $\bar{Q}$  is obtained by forgetting signs.

*Proof.* We consider all the possible ways of "adding signs" to the balls in the top row of  $\overline{Q}$  to obtain a signed two-line queue. Fix such a ball B labeled a > 0:

- if B has below it a vacancy or a ball labeled b < a, then we must assign a sign to B.
- if B has a ball labeled b > a below it, then we must assign a + sign to B.
- if B has a ball labeled b = a below it, then we can give B a + or sign.

We then check that the possible signs for each ball B is consistent with the choice of weights in Equation (6). In particular, one notices that when a ball B is given a - sign, the ball weight should be multiplied by -1 when we go from  $\bar{Q}$  to Q, but the weight of the pairing connected to B is also multiplied by -1.

Given  $\kappa \in S_n(\nu)$ , we define  $c_{\nu}^{\kappa}$  by

(7) 
$$c_{\nu}^{\kappa} := \sum_{\alpha: \|\alpha\| = \kappa} \operatorname{wt}_{\alpha} b_{\nu}^{\alpha}.$$

We get the following corollary obtained by combining Equation (4) and Lemma 3.2.

**Lemma 3.3.** Fix  $\lambda$  a partition with distinct parts, and  $\kappa, \mu \in S_n(\lambda)$ . Then

$$c^{\kappa}_{\mu} = \sum_{\bar{Q} \in \bar{\mathcal{G}}^{\kappa}_{\kappa}} \operatorname{wt}(\bar{Q}).$$

Since  $\lambda$  has distinct parts,  $\bar{\mathcal{G}}^{\kappa}_{\nu}$  is either empty or contains exactly one element.

Fix a weakly order-preserving function  $\phi : \mathbb{N} \to \mathbb{N}$ . Fix two partitions  $\lambda$  and  $\kappa$  such that  $\phi(\lambda) = \kappa$ . For  $\eta \in S_n(\kappa)$ , define

$$G_{\eta}^*(\boldsymbol{x};t) := \sum_{\rho \in S_n(\lambda): \phi(\rho) = \eta} F_{\rho}^*(\boldsymbol{x};1,t).$$

Let  $G_{\eta}$  be the top homogeneous part of  $G_{\eta}^*$ .

The following is an analogue of [AMW25, Theorem 4.18], and can be proved in essentially the same way, using interpolation analogues of results from [AS19].

**Theorem 3.4.** Fix  $\lambda$  and  $\kappa$  as above. For all  $\eta \in S_n(\kappa)$ , we have at q = 1 that

$$\frac{G_{\eta}^*(\boldsymbol{x};t)}{P_{\lambda}^*(\boldsymbol{x};1,t)} = \frac{F_{\eta}^*(\boldsymbol{x};1,t)}{P_{\kappa}^*(\boldsymbol{x};1,t)}.$$

Given a composition  $\rho$ , let  $\rho^- := (\rho_1^-, \dots, \rho_n^-)$ , where  $\rho_i^- = \max(\rho_i - 1, 0)$ .

**Corollary 3.5.** Consider a composition  $\rho$  with  $\rho_i \neq 1$  for any  $1 \leq i \leq n$ . Let k be the number of non-zero parts of  $\rho$ . Set  $\eta = \rho^-$ . We then have at q = 1,

$$F_{\rho}^{*}(\boldsymbol{x};1,t) = F_{\eta}^{*}(\boldsymbol{x};1,t) \cdot e_{k}^{*}(\boldsymbol{x};t).$$

Proof. Let  $\lambda$  and  $\kappa$  be the two partitions obtained by reordering  $\rho$  and  $\eta$ , respectively. Consider the weakly order-preserving function  $\phi: i \mapsto \max(i-1,0)$ . We then have  $\phi(\rho) = \eta$ . Since  $\lambda$  does not have parts of size 1, and  $\phi$  is bijective from  $\{0,2,3,\ldots\}$  to  $\{0,1,2,\ldots\}$ , then  $\rho$  is the unique composition in  $S_n(\lambda)$  such that  $\phi(\rho) = \eta$  and we have  $G_{\eta}^* = F_{\rho}^*$ . It follows then from Theorem 3.4 that

$$\frac{F_{\rho}^{*}(\boldsymbol{x};1,t)}{P_{\lambda}^{*}(\boldsymbol{x};1,t)} = \frac{F_{\eta}^{*}(\boldsymbol{x};1,t)}{P_{\kappa}^{*}(\boldsymbol{x};1,t)}.$$

We now recall that at q = 1, we have from [Dol17, BDW25] that

(8) 
$$P_{\lambda}^{*}(x_{1},\ldots,x_{n};1,t) = \prod_{1 \leq i \leq \lambda_{1}} P_{\lambda'_{i}}^{*}(x_{1},\ldots,x_{n};1,t) = \prod_{1 \leq i \leq \lambda_{1}} e_{\lambda'_{i}}^{*}(x_{1},\ldots,x_{n};t),$$

where  $\lambda'$  is the partition conjugate to  $\lambda$ . Using this plus the fact that  $\kappa$  is obtained from  $\lambda$  by removing the largest column (of size k), we get that

$$\frac{P_{\lambda}^*(\boldsymbol{x}; 1, t)}{P_{\kappa}^*(\boldsymbol{x}; 1, t)} = e_k^*(\boldsymbol{x}; t),$$

which implies that  $F_{\rho}^*(\boldsymbol{x}; 1, t) = F_{\eta}^*(\boldsymbol{x}; 1, t) \cdot e_k^*(\boldsymbol{x}; t)$ .

**Proposition 3.6.** Fix  $\rho, \nu \in S_n(\lambda)$ , and let j be the index such that  $\rho_j = 0$ . We have

$$\mathbb{P}_{j}^{(2)}(\rho,\nu) = \frac{c_{\nu}^{\rho}}{\prod_{k< j} \left(x_{k} - \frac{1}{t^{n-2}}\right) \prod_{k> j} \left(x_{k} - \frac{1}{t^{n-1}}\right)},$$

or equivalently,

$$P_{j} \cdot \mathbb{P}_{j}^{(2)}(\rho, \nu) = \frac{c_{\nu}^{\rho}}{e_{n-1}^{*}},$$

where  $c_{\nu}^{\rho}$  is the coefficient from Equation (7), i.e. the generating function for the set  $\bar{\mathcal{G}}_{\nu}^{\rho}$ .

The idea of the proof below is that a signed two-line queue encodes Step 2 of the interpolation t-Push TASEP.

Proof. Note that (Step 2) of Definition 2.1 is encoded by an element of a set  $\bar{\mathcal{G}}^{\rho}_{\nu}$  (see Definition 3.1). Indeed, the transition in (Step 2) from the configuration  $\rho$  to the configuration  $\nu$  is possible if and only there is an element  $\bar{Q}$  in  $\bar{\mathcal{G}}^{\rho}_{\nu}$  (recall that this set contains at most one element). More precisely, a particle labeled a > 0 which moved from position  $k \in [n]$  to a position k', corresponds to a non trivial pairing in  $\bar{Q}$  connecting a ball labeled a in column k of the top row to a ball labeled a in column k' of the bottom row. Particles which do not move correspond to trivial pairings.

We now claim that  $\operatorname{wt}(\bar{Q})$  divided by  $D := \prod_{k < j} \left( x_k - \frac{1}{t^{n-2}} \right) \prod_{k > j} \left( x_k - \frac{1}{t^{n-1}} \right)$  gives  $\mathbb{P}_j^{(2)}(\rho, \nu)$ . We will prove the claim below by showing that each ball or pairing weight in  $\operatorname{wt}(\bar{Q})$ , divided by one of the factors in D, equals one of the skipping/ displacement probabilities from Item (Step 2) (whose product is  $\mathbb{P}_j^{(2)}(\rho, \nu)$ ). Note that in what follows, instead of associating the weight  $(1-t)t^{\operatorname{skip}(p)}$  to each nontrivial pairing, we will associate (1-t) to the top ball in each nontrivial pairing, and a factor of t to each skipped ball.

- Each ball in column k > j of  $\bar{Q}$  is necessarily trivially paired, since no ball in position k > j get skipped or displaced in (Step 2). In  $\bar{Q}$  this ball gets weight  $x_k \frac{1}{t^{n-1}}$ ; when we divide this weight by the kth factor of D, we get 1, which corresponds to the fact that balls in position k > j do not contribute to  $\mathbb{P}_j^{(2)}(\rho, \nu)$ .
- A ball in  $\overline{Q}$  labeled b in column k < j which is trivially paired, and which is not skipped by a ball a > b, also has weight  $x_k \frac{1}{t^{n-1}}$ . When we divide this weight by the kth factor of D, we get  $1 \mathfrak{p}_k$  (see (1)). This is what we desired, because such a trivial pairing in  $\overline{Q}$  corresponds to a particle labeled b which is skipped over by a particle with a smaller label, and hence contributes  $1 \mathfrak{p}_k$  to  $\mathbb{P}_j^{(2)}(\rho, \nu)$ .
- A ball in  $\bar{Q}$  labeled b in column k < j which is trivially paired, and which is skipped by a ball a > b, gets a weight  $t(x_k \frac{1}{t^{n-1}})$ . When we divide this weight by the kth factor of D, we get  $1 \mathfrak{q}_k$  (see (1)). This is what we desired, because such a trivial pairing corresponds to a particle labeled b skipped over by a particle with a larger label, and hence contributes  $1 \mathfrak{q}_k$  to  $\mathbb{P}_j^{(2)}(\rho, \nu)$ .
- A ball labeled b in the top row of  $\bar{Q}$  in column k < j which has a ball labeled a < b below it gets a weight  $(1-t)\frac{1}{t^{n-1}}$  (the factor (1-t) is the nontrivial pairing weight). When we divide this weight by the kth factor of D, we get  $\mathfrak{p}_k$ . This is what we desired, because this pairing corresponds to a particle labeled b being displaced by a particle with a smaller label, and hence contributing  $\mathfrak{p}_k$  to  $\mathbb{P}_j^{(2)}(\rho,\nu)$ .
- A ball labeled b in the top row of  $\overline{Q}$  in column k < j which has a ball labeled a > b below it gets a weight  $(1-t)x_k$  (the factor (1-t) is the nontrivial pairing weight). When we divide this weight by the kth factor of D, we get  $\mathfrak{q}_k$ . This is what we desired, because this pairing corresponds to a particle labeled b being displaced by a particle with a larger label, and hence contributing  $\mathfrak{q}_k$  to  $\mathbb{P}_j^{(2)}(\rho,\nu)$ .

**Proposition 3.7.** If  $\lambda$  is restricted, and  $\mu, \nu \in S_n(\lambda)$ , then

$$\mathbb{P}(\mu,\nu) = \sum_{\rho \in S_n(\lambda)} \frac{a_\rho^\mu c_\nu^\rho}{e_{n-1}^*}.$$

*Proof.* Combining [AMW25, Lemma 5.4] and Proposition 3.6, we get

$$\mathbb{P}(\mu,\nu) = \sum_{1 \leq j \leq n} P_j \sum_{\rho \in S_n(\lambda): \rho_j = 0} \mathbb{P}_j^{(1)}(\mu,\rho) \mathbb{P}_j^{(2)}(\rho,\nu)$$

$$= \sum_{1 \leq j \leq n} \sum_{\rho \in S_n(\lambda): \rho_j = 0} \frac{a_\rho^\mu c_\nu^\rho}{e_{n-1}^*}$$

$$= \sum_{\rho \in S_n(\lambda)} \frac{a_\rho^\mu c_\nu^\rho}{e_{n-1}^*}.$$

Proof of Theorem 2.2. Fix a restricted partition  $\lambda$ . Let  $\nu \in S_n(\lambda)$ . From [BDW25, Theorem 1.15 and Lemma 5.6], we have

$$F_{\nu}^{*}(\boldsymbol{x};1,t) = \sum_{\eta \in \mathbb{N}^{n}} F_{\nu}^{*\eta}(\boldsymbol{x};t) F_{\eta^{-}}^{*}(\boldsymbol{x};1,t),$$

where

$$F_{\nu}^{*\eta}(\boldsymbol{x};t) := \sum_{\alpha \in \mathbb{Z}^n} b_{\nu}^{\alpha} \operatorname{wt}_{\alpha} a_{\|\alpha\|}^{\eta} = \sum_{\kappa \in \mathbb{N}^n} a_{\kappa}^{\eta} c_{\nu}^{\kappa}.$$

But we know from Corollary 3.5 that

$$F_{\eta^{-}}^{*}(\boldsymbol{x};1,t) = \frac{F_{\eta}^{*}(\boldsymbol{x};1,t)}{e_{\eta-1}^{*}(\boldsymbol{x};t)},$$

we use here the fact that  $\eta$  has a unique part of size 0. Hence

$$F_{\nu}^{*}(\boldsymbol{x};1,t) = \sum_{\eta \in \mathbb{N}^{n}} F_{\eta}^{*}(\boldsymbol{x};1,t) \sum_{\kappa \in \mathbb{N}^{n}} \frac{a_{\kappa}^{\eta} c_{\nu}^{\kappa}}{e_{n-1}^{*}(\boldsymbol{x};t)},$$

which can be rewritten using the transition probabilities of the interpolation t-Push TASEP (Proposition 3.7) we get

$$F_{\nu}^{*}(\boldsymbol{x};1,t) = \sum_{\eta \in \mathbb{N}^{n}} F_{\eta}^{*}(\boldsymbol{x};1,t) \mathbb{P}(\eta,\nu).$$

This proves that  $F_{\mu}^{*}(\boldsymbol{x}; 1, t)$  are proportional to the stationary distribution of the interpolation t-Push TASEP  $\pi_{\lambda}^{*}(\mu)$ . Finally, we use the fact that  $P_{\lambda}^{*} = \sum_{\mu \in S_{n}(\lambda)} F_{\mu}^{*}$  to deduce that  $\frac{F_{\mu}^{*}(\boldsymbol{x}; 1, t)}{P_{\lambda}^{*}(\boldsymbol{x}; 1, t)} = \pi_{\lambda}^{*}(\mu)$ .

#### REFERENCES

[AMW25] Arvind Ayyer, James Martin, and Lauren Williams, The inhomogeneous t-PushTASEP and Macdonald polynomials at q=1, Annales de l'Institut Henri Poincaré D (2025).

[AS19] Per Alexandersson and Mehtaab Sawhney, Properties of non-symmetric Macdonald polynomials at q=1 and q=0, Ann. Comb. 23 (2019), no. 2, 219–239. MR 3962853

[BDW25] Houcine Ben Dali and Lauren Williams, A combinatorial formula for Interpolation Macdonald polynomials , Preprint arXiv:2510.02587, 2025.

[CMW22] Sylvie Corteel, Olya Mandelshtam, and Lauren Williams, From multiline queues to Macdonald polynomials via the exclusion process, Amer. J. Math. 144 (2022), no. 2, 395–436. MR 4401508

[Doł17] Maciej Dołęga, Strong factorization property of Macdonald polynomials and higher-order Macdonald's positivity conjecture, J. Algebraic Combin. 46 (2017), no. 1, 135–163. MR 3666415

DEPARTMENT OF MATHEMATICS, HARVARD UNIVERSITY, CAMBRIDGE, MA, AND CENTER FOR MATHEMATICAL SCIENCES AND APPLICATIONS, HARVARD UNIVERSITY, CAMBRIDGE, MA *Email address*: bendali@math.harvard.edu

DEPARTMENT OF MATHEMATICS, HARVARD UNIVERSITY, CAMBRIDGE, MA *Email address*: williams@math.harvard.edu
