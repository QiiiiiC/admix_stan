# Relaxed graph with parameters as constants: the mathematics

Equation numbers are the ones cited in `stan/*.stan`.  Time is in generations
backwards from the present; genetic length in cM; $N$ is a haploid size (two
lineages in one deme coalesce at rate $1/N$).

## 1. Likelihood of the data given the slot matrices

**Data.**  $L$ leaves with $n_i$ haploid samples each.  For every leaf pair $i\le j$
and length bin $b=[u_b,v_b)$: the IBD segment count $k_{b,ij}$ over a genome of
$C$ cM.  For the SNP part: the centred allele-frequency covariance $\hat w_{ij}$
with standard error $s_{ij}$.

**Parameters** (unknown constants, not random variables):
$$\theta=(\Delta_1,\dots,\Delta_K,\ f,\ N),\qquad \Delta_k\ge 1,\qquad t_k=\sum_{m\le k}\Delta_m,\quad t_0=0,\quad t_{K+1}=t_K+T_{\max}. \tag{1}$$

**Demes and events.**  $D$ demes, $K$ events at $t_1<\dots<t_K$, epoch $k=[t_{k-1},t_k)$,
$k=1,\dots,K+1$.  Event $k$ moves every lineage between demes with a row-stochastic
matrix $M_k$: $(M_k)_{ab}=\Pr(\text{deme } b \text{ just older than } t_k \mid \text{deme } a \text{ just younger})$.
Within an epoch nothing moves.

**Pair state.**  For one lineage from leaf $i$ and one from leaf $j$,
$W_{ab}(t)=\Pr(\text{lineage 1 in } a,\ \text{lineage 2 in } b,\ \text{not coalesced by } t)$:
$$W(0)=e_ie_j^{\top};\qquad W_{aa}(t)=W_{aa}(t_{k-1})\,e^{-(t-t_{k-1})/N},\ \ W_{ab}(t)=W_{ab}(t_{k-1})\ (a\ne b),\quad t\in\text{epoch }k; \tag{2}$$
$$W(t_k^{+}) = M_k^{\top}\,W(t_k^{-})\,M_k \qquad(\text{the two lineages move independently}). \tag{3}$$

**Coalescence density** of the pair, for $t$ in epoch $k$:
$$\varphi_{ij}(t)=\frac1N\sum_{a=1}^{D}W_{aa}(t_{k-1})\,e^{-(t-t_{k-1})/N}. \tag{4}$$

**IBD segments.**  Given a pair's TMRCA $t$, under the SMC the segment ends arrive at
rate $2t$ per Morgan, i.e. $t/50$ per cM, so segments start at rate $t/50$ per cM
and their lengths are Exponential with rate $t/50$.  The expected number of
segments per pair per cM with length in bin $b$ is
$$\nu_{b,ij}=\int_0^{\infty}\varphi_{ij}(t)\,\frac{t}{50}\Big(e^{-u_bt/50}-e^{-v_bt/50}\Big)\,dt. \tag{5}$$
Within an epoch with $W_{aa}(t_{k-1})=\omega$ the integral is closed-form; with
$c_u=-(1/N+u/50)$,
$$\frac{\omega}{50N}\int_{t_{k-1}}^{t_k} t\,e^{-(t-t_{k-1})/N-ut/50}\,dt=\frac{\omega}{50N}\Big[e^{-(t-t_{k-1})/N-ut/50}\,\frac{c_ut-1}{c_u^2}\Big]_{t_{k-1}}^{t_k}, \tag{6}$$
minus the same with $v$ (`segment_density`).  The observed counts are Poisson:
$$k_{b,ij}\sim\text{Poisson}\big(\lambda_{b,ij}\big),\qquad \lambda_{b,ij}=\nu_{b,ij}\,C\,P_{ij},\qquad P_{ii}=\tbinom{n_i}{2},\ P_{ij}=n_in_j .$$

**SNP drift** accumulates while the two lineages share a deme, up to the root and
without removing coalesced mass (Brownian drift):
$$\Omega_{ij}=\sum_{k=1}^{K}\sum_{a}\tilde W_{aa}^{(k)}\,\frac{\Delta_k}{N}, \tag{7}$$
where $\tilde W$ follows (3) with no decay.  Double centring
$\bar\Omega_{ij}=\Omega_{ij}-\bar\Omega_{i\cdot}-\bar\Omega_{\cdot j}+\bar\Omega_{\cdot\cdot}$
and $\hat w_{ij}\sim\mathcal N(\bar\Omega_{ij},s_{ij}^2)$ (composite, independent of the IBD term).

**Log-likelihood** of the data given $\theta$ and the matrices $M=(M_1,\dots,M_K)$:
$$\ell(\theta,M)=\sum_{i\le j}\sum_b\log\text{Pois}\big(k_{b,ij}\mid\lambda_{b,ij}\big)+\sum_{i\le j}\log\mathcal N\big(\hat w_{ij}\mid\bar\Omega_{ij},s_{ij}^2\big). \tag{8}$$
This is exactly the likelihood of `methods/stan/mixed_Nfixed_poisson.stan`, written
over demes instead of graph nodes (tested: equal at every binary vertex).

## 2. The skeleton and its discrete graphs

Skeleton $[\text{admixture},\ \text{merge}\times(L-1),\ \text{root}]$, so $K=L+1$
slots and $D=L+1$ demes: the leaves $1..L$ and a ghost $g=L+1$ (the admixture's
second source).  $e_a$ is the $a$-th unit vector.

*Admixture slot* ($k=1$); choice $x\in\{1..L\}$ = the admixed leaf:
$$A^{(x)}(f)=I-e_xe_x^{\top}+e_x\big(f\,e_x+(1-f)\,e_g\big)^{\top}. \tag{9}$$
Leaf $x$'s lineage stays in deme $x$ with probability $f$ (source 1) and moves to the
ghost with $1-f$ (source 2).

*Merge slot* ($k=2..K-1$); choice $p=\{a<b\}$ among the $\binom D2$ pairs:
$$B^{(a,b)}=I-e_be_b^{\top}+e_be_a^{\top}\qquad(\text{deme } b \text{ moves into deme } a). \tag{10}$$
Always merging into the lower index removes the relabelling duplicate $a\leftrightarrow b$.

*Root slot* ($k=K$): **compulsory, not a choice** — every deme moves into deme 1,
$$R=\mathbf 1\,e_1^{\top}, \tag{11}$$
so every pair coalesces in the final epoch whatever the earlier slots did.

A **vertex** $v=(x,p_2,\dots,p_{K-1})$ is one choice per slot, with
$M(v)=(A^{(x)}(f),B^{(p_2)},\dots,B^{(p_{K-1})},R)$.  For $L=3$: $3\times6\times6=108$
vertices.  Two things make different vertices the same graph:
* the admixture's two sources are exchangeable: swapping which one is "source 1"
  (and $f\leftrightarrow1-f$) is the same graph;
* a merge that names an already-empty deme moves nothing, so the root then joins
  three lineages at once.

The 108 vertices are **45 distinct graphs** (15 per admixed leaf): 33 binary (21
with two mirror vertices, 12 with one — the bubbles, whose mirror is themselves)
and 12 with a three-way root (3 or 6 vertices each).
Note that with $\theta$ shared (section 3) the two mirror vertices of a graph are
*not* equivalent at a given $f$; they are the same graph at $f$ and $1-f$.  The prior
is uniform on slot choices,
$$p(v)=\prod_k|C_k|^{-1}=1/108,\qquad p(G)=|\{v: v\mapsto G\}|/108 .$$

## 3. EM with $\theta$ as constants (the discrete model, exact)

The graph is the only random variable.  The **marginal likelihood of $\theta$**
sums it out:
$$p(x\mid\theta)=\sum_{v}p(v)\,e^{\ell(\theta,M(v))},\qquad \hat\theta=\arg\max_\theta\log p(x\mid\theta). \tag{12}$$
The graph posterior at $\hat\theta$ is a softmax over vertices:
$$p(v\mid x,\hat\theta)=\frac{p(v)\,e^{\ell(\hat\theta,v)}}{\sum_{v'}p(v')\,e^{\ell(\hat\theta,v')}},\qquad p(G\mid x,\hat\theta)=\sum_{v\mapsto G}p(v\mid x,\hat\theta). \tag{13}$$

**EM.**  With $r^{(m)}_v=p(v\mid x,\theta^{(m)})$ (E-step), the M-step is
$$\theta^{(m+1)}=\arg\max_\theta\ Q(\theta\mid\theta^{(m)}),\qquad Q(\theta\mid\theta^{(m)})=\sum_v r^{(m)}_v\,\ell(\theta,v), \tag{14}$$
and $\log p(x\mid\theta^{(m)})$ never decreases.  By Fisher's identity
$\nabla_\theta\log p(x\mid\theta)=\sum_v r_v(\theta)\nabla_\theta\ell(\theta,v)$, so plain gradient ascent on
(12) has the same fixed points; `discrete.stan` maximises (12) directly with L-BFGS
(the E-step is inside the gradient).

*Profile, for comparison* (every graph with its own parameters):
$\hat\ell_v=\max_\theta\ell(\theta,v)$, one optimisation per vertex (`discrete.stan` with $V=1$).

What (12)–(13) are and are not:
* $\theta$ is **shared** by all vertices, and every vertex has the same parameters
  (all $K$ slot times, $f$, $N$), so no graph gains by having more parameters; the
  slot times mean the same thing in every graph ("the second event happens at $t_2$").
* It is a plug-in: $p(G\mid x,\hat\theta)$, not the full Bayes
  $\int p(G\mid x,\theta)\,p(\theta\mid x)\,d\theta$.  There is no prior on $\theta$ and no Occam factor
  for it.  A graph whose own best $\theta$ is far from $\hat\theta$ is scored at $\hat\theta$,
  not at its own optimum (that is what the profile shows).
* The E-step costs $|V|$ likelihood evaluations: fine for 108 vertices, not for large
  skeletons.  That is what the relaxation removes.

## 4. The relaxed graph

**Gumbel-max.**  Give each slot a vector of logits with iid standard Gumbel priors,
$$y_{k,c}\sim\text{Gumbel}(0,1)\ \text{iid}\quad\Longrightarrow\quad H_k=\arg\max_c y_{k,c}\sim\text{Uniform}(C_k), \tag{15}$$
so the hard choice $H=(H_1,\dots,H_{K-1})$ has exactly the prior $p(v)$ of section 2.

**Soft slot matrices.**  At temperature $\lambda>0$, $w_k=\text{softmax}(y_k/\lambda)$ and
$$M_1(w)=\sum_x w_{1,x}A^{(x)}(f),\qquad M_k(w)=\sum_p w_{k,p}B^{(p)},\qquad M_K=R. \tag{16}$$
Entry-wise: row $x$ of $M_1$ keeps $1-w_{1,x}(1-f)$ and sends $w_{1,x}(1-f)$ to the ghost;
row $b$ of a merge slot keeps $1-\sum_{a<b}w_{k,ab}$ and sends $w_{k,ab}$ to each lower $a$.
Convex combinations of row-stochastic matrices are row-stochastic, so every $w$ is a
valid pulse-migration history; the vertices are the one-hot $w$.  (In the interior,
only the products $w_{1,x}(1-f)$ enter, so $f$ and $w_1$ trade off along a ridge.)

**Relaxed model and its limit.**
$$p_\lambda(x\mid\theta)=\mathbb E_{y}\big[e^{\ell(\theta,M(w(y)))}\big]\ \xrightarrow{\ \lambda\to0\ }\ \sum_vp(v)\,e^{\ell(\theta,v)}=p(x\mid\theta), \tag{17}$$
because $w(y)\to\text{onehot}(H)$ almost surely (ties have probability 0), $\ell$ is
continuous in $M$, and $e^{\ell}$ is bounded (a Poisson pmf is at most 1, the SNP
density at most $\prod(2\pi s_{ij}^2)^{-1/2}$), so dominated convergence applies.  Likewise
$$p_\lambda(H=v\mid x,\theta)=\frac{\mathbb E_y\big[\mathbf 1\{H(y)=v\}\,e^{\ell(\theta,M(w(y)))}\big]}{p_\lambda(x\mid\theta)}\ \xrightarrow{\ \lambda\to0\ }\ p(v\mid x,\theta). \tag{18}$$
So the relaxation targets exactly the discrete EM of section 3 as $\lambda\to0$; at
$\lambda>0$ it also weighs interior (non-graph) histories.

**Monte Carlo EM.**  $y$ is the latent variable and its prior does not involve
$\theta$, so $Q(\theta\mid\theta^{(m)})=\mathbb E_{y\mid x,\theta^{(m)}}[\ell(\theta,M(w(y)))]+\text{const}$:
* E-step: draw $y^{(1..S)}\sim p_\lambda(y\mid x,\theta^{(m)})\propto\prod\text{Gumbel}(y)\,e^{\ell(\theta^{(m)},M(w(y)))}$
  with NUTS (`relaxed_estep.stan`, $\theta$ as data);
* M-step:
$$\theta^{(m+1)}=\arg\max_\theta\ \frac1S\sum_{s=1}^{S}\ell\big(\theta,M(w(y^{(s)}))\big)\qquad(\texttt{relaxed\_mstep.stan}); \tag{19}$$
* stop when $\max|\Delta\log\theta|$ is below a tolerance; then read off
  $p_\lambda(H=v\mid x,\hat\theta)$ as the share of E-step draws whose argmax pattern is $v$.

Each likelihood evaluation costs one $\ell$ whatever the number of vertices.  The
price is the E-step's geometry: as $\lambda\to0$ the posterior of $y$ becomes flat
inside each vertex's cell with steep walls between cells, so chains stop crossing
between graphs; at larger $\lambda$ they cross, but interior histories carry weight.
Comparing the relaxed read-off with the exact (13) across $\lambda$ measures that
trade-off.
