# Model

scHINT (single-cell gene expression **H**eritability with
gene-by-cell-context **INT**eractions) estimates, **for one gene at a
time**, how the variance of expression divides into

- a **main genetic** component that is the same in every cell,
- a **genotype-by-context (G×C)** component whose size changes with the
  cell’s state or type,
- an **individual** component shared by all cells of a donor but not
  genetic,
- **context** and **covariate** components, and
- cell-level residual noise.

The same Haseman-Elston (HE) machinery is used in three settings, which
differ only in what “context” is and in what links observations:

| setting | function | unit of observation | genetic kernel | context |
|----|----|----|----|----|
| cell-level | [`schint_cell()`](https://LidaWangPSU.github.io/scHINT/reference/schint_cell.md) | cell | cis-GRM between donors | continuous cell state (PCs, pseudotime, …) and/or categorical |
| sample-level | [`schint_sample()`](https://LidaWangPSU.github.io/scHINT/reference/schint_sample.md) | donor × cell type (pseudobulk) | cis-GRM between donors | cell type |
| perturbation | [`schint_pert()`](https://LidaWangPSU.github.io/scHINT/reference/schint_pert.md) | cell | “same perturbation” | cell state *or* cell type (one at a time) |

## 1. Generative model (cell-level, G×Cell-State)

Let \\y\_{im}\\ be the normalized expression of the gene in cell \\m\\
of individual \\i\\ (\\i=1,\dots,n\\; \\m=1,\dots,M_i\\),
\\g_i\in\mathbb{R}^p\\ the standardized cis-genotypes (by default all
SNPs within ±500 kb of the gene), \\c\_{im}\in\mathbb{R}^K\\ the
cell-state features (for example the leading within-cell-type principal
components) and \\o\_{im}\in\mathbb{R}^L\\ further covariates. Then

\\ y\_{im}= g_i^\top\beta + (g_i\otimes c\_{im})^\top\gamma +
c\_{im}^\top\alpha + o\_{im}^\top\eta + e_i + \varepsilon\_{im}, \\

with independent random effects

\\ \beta\sim N\\\left(0,\tfrac{\sigma_G^2}{p}I_p\right),\quad
\gamma_k\sim N\\\left(0,\tfrac{\sigma\_{G\times
C,k}^2}{p}I_p\right),\quad e_i\sim N(0,\sigma_I^2),\quad
\varepsilon\_{im}\sim N(0,\sigma\_\varepsilon^2). \\

\\\beta\\ is the main cis-genetic effect, \\\gamma_k\\ is how that
effect is modulated along cell-state axis \\k\\, \\\alpha\\ and \\\eta\\
are the non-genetic effects of state and covariates, \\e_i\\ is the
donor effect shared across cells and \\\varepsilon\_{im}\\ is cell-level
noise. In the package \\\sigma_C^2\\ and \\\sigma_o^2\\ denote the
variance contributed by \\c^\top\alpha\\ and \\o^\top\eta\\ (treated as
random).

## 2. Haseman-Elston regression on pairs of cells

Define the cis-genetic relatedness (GRM)

\\K^G\_{ij}=\frac{g_i^\top g_j}{p}.\\

For two different cells \\(i,m)\neq(j,n)\\,

\\ \mathbb{E}(y\_{im}y\_{jn})=\sigma_G^2K^G\_{ij}
+\sum\_{k=1}^{K\_\text{int}}\sigma\_{G\times
C,k}^2K^G\_{ij}c\_{imk}c\_{jnk} +\sigma_I^2\mathbf 1(i=j)
+\sum\_{k=1}^{K}\sigma\_{C,k}^2c\_{imk}c\_{jnk}
+\sum\_{l=1}^{L}\sigma\_{O,l}^2o\_{iml}o\_{jnl}. \\

The residual term \\\sigma\_\varepsilon^2\mathbf 1(i=j,m=n)\\ vanishes
for distinct cells. Because every term is a known *kernel* times an
unknown variance, the variances are estimated by ordinary least squares
of the cross-product \\y\_{im}y\_{jn}\\ on the kernels over all pairs:

\\ y\_{im}y\_{jn}\sim 1 + K^G\_{ij} + \mathbf 1(i=j) + \sum_k
c\_{imk}c\_{jnk} + \sum_l o\_{iml}o\_{jnl} + \sum\_{k=1}^{K\_\text{int}}
K^G\_{ij}\\c\_{imk}c\_{jnk}. \\

The coefficients are the variance components \\\hat\sigma^2\\ (the
leading `1` absorbs the mean). The estimator is moment based: no
likelihood, no iteration, and it remains valid without Gaussian errors
(the variances may come out slightly negative; see
[Troubleshooting](https://LidaWangPSU.github.io/scHINT/articles/troubleshooting.md)).

**Mapping to package output** (`fit$coefficients$component`):

| symbol | component name | `group` |
|----|----|----|
| \\\sigma_G^2\\ | `G` (or `G1`, `G2`, … with several GRMs) | `G` |
| \\\sigma\_{G\times C,k}^2\\ | `G:<context variable>` | `GxC` |
| \\\sigma_I^2\\ | `I` | `I` |
| \\\sigma\_{C,k}^2\\ | `<context variable>` | `context` |
| \\\sigma\_{O,l}^2\\ | `<covariate>` | `covariate` |
| \\\sigma\_{G\times C}^2=\sum_k\sigma\_{G\times C,k}^2\\ | summary row `GxC` |  |

Variables are standardized (mean 0, variance 1) by default, so each
estimate is directly a variance contribution on the scale of
standardized expression.

### Several GRMs

Passing a list of GRMs in `grm` (e.g. built from SNPs split by MAF bin
or annotation with
[`split_snps()`](https://LidaWangPSU.github.io/scHINT/reference/split_snps.md)
and
[`make_grm()`](https://LidaWangPSU.github.io/scHINT/reference/make_grm.md))
replaces \\K^G\\ by \\K^{G,1},\dots,K^{G,M}\\, each with its own main
component \\\sigma^2\_{G,m}\\ and its own interaction
\\\sigma^2\_{G\times C,m,k}\\. Components of different GRMs enter the
normal equations through the element-wise products \\K^{G,m}\circ
K^{G,m'}\\.

### Categorical contexts

For a categorical context with levels \\t\\, the pair product
\\c\_{imk}c\_{jnk}\\ is replaced by \\\mathbf 1(t\_{im}=t\_{jn})\\: the
G×C kernel becomes \\K^G\_{ij}\mathbf 1(t\_{im}=t\_{jn})\\ (genetic
variance specific to pairs from the same level) and the main-effect
kernel becomes \\\mathbf 1(t\_{im}=t\_{jn})\\. With
`cat_mode = "per_level"` each level gets its own component, giving a
variance per level.

## 3. Sample-level G×Cell-Type model

For pseudobulk expression \\y\_{it}\\ of individual \\i\\ in cell type
\\t=1,\dots,T\\,

\\
y\_{it}=g_i^\top\beta+g_i^\top\gamma_t+\alpha_t+o\_{it}^\top\eta+e_i+\varepsilon\_{it},
\\

where \\\beta\\ is the genetic effect shared across cell types,
\\\gamma_t\\ the cell-type-specific deviation (\\\gamma_t\sim
N(0,\sigma^2\_{G\times CT}p^{-1}I_p)\\), \\\alpha_t\sim
N(0,\sigma_T^2)\\ the cell-type mean and \\e_i\\ a donor effect shared
across cell types. For two observations \\(i,s)\ne(j,t)\\

\\ \mathbb{E}(y\_{is}y\_{jt})=\sigma_G^2K^G\_{ij}+\sigma^2\_{G\times
CT}K^G\_{ij}\mathbf 1(s=t) +\sigma\_{CT}^2\mathbf
1(s=t)+\sigma_I^2\mathbf 1(i=j). \\

[`schint_sample()`](https://LidaWangPSU.github.io/scHINT/reference/schint_sample.md)
regresses \\y\_{is}y\_{jt}\\ on these four kernels (plus covariate
kernels). \\\sigma_G^2\\ is the cross-cell-type shared genetic variance
and \\\sigma^2\_{G\times CT}\\ the cell-type-specific genetic variance.
With `gxc = FALSE` only the shared-genetics model is fitted.

## 4. Perturbation model

For a perturbation screen the unit that links cells is the
**perturbation** rather than the donor. With \\p_m\\ the perturbation of
cell \\m\\ and \\c_m\in\mathbb{R}^K\\ the (standardized) cell state,

\\ y_m=b\_{p_m}+\sum_k c\_{mk}\alpha_k+\sum_k
c\_{mk}u\_{p_mk}+\varepsilon_m, \\

\\ b_p\sim N(0,\sigma_P^2),\quad \alpha_k\sim N(0,\sigma\_{C,k}^2),\quad
u\_{pk}\sim N(0,\sigma\_{P\times C,k}^2), \\

so \\\sigma^2_P\\ is the perturbation-associated variance shared across
states and \\\sigma^2\_{P\times C,k}\\ the variance of the perturbation
response along axis \\k\\. For two cells \\m\neq n\\

\\ y_my_n\sim 1+\mathbf 1(p_m=p_n)+\sum_k c\_{mk}c\_{nk}+\mathbf
1(p_m=p_n)\sum_k c\_{mk}c\_{nk}, \\

i.e. the genetic kernel is replaced by “same perturbation”. With a
categorical context (cell type) the last two terms use \\\mathbf
1(t_m=t_n)\\. Perturbations can be split into groups, each with its own
\\\sigma_P^2\\ and \\\sigma\_{P\times C}^2\\ (`perturb_group`). All
genes measured in the same cells share the same pair kernels, so
[`schint_pert()`](https://LidaWangPSU.github.io/scHINT/reference/schint_pert.md)
builds \\X^\top X\\ once and re-uses it for every gene.

## 5. Scalable estimation from sufficient statistics

Forming all cell pairs costs \\O(M^2K^2)\\ time and \\O(M^2)\\ memory
for \\M\\ cells, which is impossible for millions of cells. Every
regressor in the models above has the form

\\ X_r\[(i,m),(j,n)\] = \kappa_r(i,j)\\\sum_q
u\_{rq}(i,m)\\u\_{rq}(j,n), \\

a donor kernel \\\kappa_r\\ (all-ones, \\\mathbf 1(i=j)\\, a GRM, or
“same perturbation”) times a pair product of cell-level features \\u_r\\
(a column of ones, a standardized covariate, or a one-hot encoding). The
normal equations need only

\\ (X^\top X)\_{rs}=\sum\_{a\<b}X_r\[a,b\]X_s\[a,b\],\qquad (X^\top
y)\_r=\sum\_{a\<b}X_r\[a,b\]\\y_ay_b . \\

Write \\z_a=u_r(a)\odot u_s(a)\\ and let \\W_i=\sum\_{a\in i}z_a\\ be
its sum over the cells of donor \\i\\. Then

\\ (X^\top
X)\_{rs}=\tfrac12\Big\[\sum\_{i,j}(\kappa_r\circ\kappa_s)\_{ij}\\W_i^\top
W_j -\sum_a(\kappa_r\circ\kappa_s)\_{i(a)i(a)}\\\lVert
z_a\rVert^2\Big\], \\

and \\(X^\top y)\_r\\ is the same expression with \\z_a=u_r(a)\\y_a\\
and kernel \\\kappa_r\\. Only the donor-level sums \\W_i\\ are needed,
so the cost is \\O(N^2K^2)\\ time and \\O(N^2)\\ memory for \\N\\
donors, independent of the number of cells except for one linear pass.
Solving the small normal equations gives the variance components. This
is exactly equivalent to regressing over all pairs (the package tests
this against an explicit pairwise OLS).

## 6. Standard errors: donor-block jackknife

Because \\X^\top X\\ and \\X^\top y\\ are sums over donor pairs,
dropping a block \\B\\ of donors only removes terms with \\i\in B\\ or
\\j\in B\\ from the same quadratic forms. All leave-block-out versions
are therefore obtained from the already computed \\W_i\\ at almost no
extra cost; no re-fitting or bootstrapping is required. With \\B\\
blocks and leave-one-block-out estimates \\\hat\theta\_{(b)}\\,

\\
\widehat{\mathrm{SE}}(\hat\theta)=\sqrt{\tfrac{B-1}{B}\sum_b\big(\hat\theta\_{(b)}-\bar\theta\_{(\cdot)}\big)^2}.
\\

By default every donor is a block (`n_blocks = NULL`); in
[`schint_pert()`](https://LidaWangPSU.github.io/scHINT/reference/schint_pert.md)
blocks are perturbations. Standard errors are also returned for derived
quantities (sums of components, heritabilities), because they are
computed from the leave-block-out component vectors.

## 7. Heritability

`fit$summary` reports each component’s variance (`variance`, on the
standardized-expression scale, so it is also the fraction of the
variance of expression). Single-cell variance contains a large
cell-specific part that averaging over cells removes. To make the
estimates comparable to population-level heritability, `h2` removes the
cell-level residual \\\sigma\_\varepsilon^2\\ from the denominator:

\\ h^2_G=\frac{\sigma_G^2}{\sigma_G^2+\sigma^2\_{G\times
C}+\sigma_C^2+\sigma_o^2+\sigma_I^2},\qquad h^2\_{G\times
C}=\frac{\sigma\_{G\times C}^2}{\sigma_G^2+\sigma^2\_{G\times
C}+\sigma_C^2+\sigma_o^2+\sigma_I^2}. \\

This is the variance that would be explained if every cell were measured
without stochastic noise. The aggregate context-dependent component is
`GxC`; `G_total` is `G + GxC`. For the perturbation model the summary
gives the variance of each component (`P`, `PxC`, `P_total`, per-context
and per-group splits).

## 8. Practical notes

- Genotypes are standardized by allele frequency,
  \\\operatorname{GRM}=XX^\top/p\\.
- The response is centred and scaled; context and covariates are
  standardized (`scale_x` for the latter).
- A covariate also listed in `context` is not duplicated: context
  variables always carry a main-effect and an interaction term.
- The model is Gaussian; it does not model count sparsity. Use
  normalized, log-transformed (and residualized if desired) expression.

*Reference.* Wang L., Fadil C., Wang L., Van Dyke K., Gusev A.
*Single-cell gene expression heritability informed by gene ×
cell-context interactions* (manuscript).
