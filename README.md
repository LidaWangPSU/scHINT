# scHINT

[![R-CMD-check](https://github.com/LidaWangPSU/scHINT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LidaWangPSU/scHINT/actions)

**s**ingle-**c**ell gene expression **H**eritability with gene-by-cell-context **INT**eractions.

scHINT partitions the variance of a gene's expression into main genetic effects, genetic
effects that depend on the **cell state** or **cell type** (G×C), donor, context and
covariate components, and estimates them with a scalable Haseman-Elston regression
with jackknife standard errors. The same framework estimates **perturbation × context**
variance in perturb-seq screens.

Documentation: <https://LidaWangPSU.github.io/scHINT/> (model, users manual, tutorials)

## Install

```r
# install.packages("remotes")
remotes::install_github("LidaWangPSU/scHINT")
```

## Functions

| function | data | G×C term |
|---|---|---|
| `schint_cell()` | one row per cell | genotype × any continuous or categorical context (cell state, etc.) |
| `schint_sample()` | one row per donor × cell type | genotype × cell type (pooled, or one variance per cell type) |
| `schint_pert()` | one row per cell, perturbation screens | perturbation × cell state and/or cell type; many genes and perturbation groups |

Common options: `geno` (cis-SNP dosages) **or** `grm` (pre-computed GRM or list of GRMs);
`n_grm`/`grm_by` for several GRMs (default: one GRM from all cis-SNPs); `context`;
`covariates`; `jackknife = TRUE`.

## Quick start

```r
library(scHINT)
data(schint_example)        # simulated; see ?schint_example
ex <- schint_example

fit <- schint_cell(ex$cell, y = "GENE_A", id = "donor", geno = ex$geno,
                   context = "state",
                   covariates = c("pc1", "pc2", "age", "sex", "batch"),
                   jackknife = TRUE)
fit            # G, GxC, G_total, I, context, covariate (+ normalized h2, SEs)

fs <- schint_sample(ex$sample, y = "GENE_A", id = "donor", celltype = "celltype",
                    geno = ex$geno, covariates = c("age", "sex"))

data(schint_pert_example)
schint_pert(schint_pert_example, y = c("GENE_A", "GENE_B"), perturb = "perturb",
            context = c("state", "celltype"), control = "NT")
```

## Model

For cell $m$ of donor $i$, with standardized cis-genotypes $g_i$, cell-state features
$c_{im}$ and covariates $o_{im}$,

$$y_{im}=g_i^\top\beta+(g_i\otimes c_{im})^\top\gamma+c_{im}^\top\alpha+o_{im}^\top\eta+e_i+\varepsilon_{im}.$$

With $K^G_{ij}=g_i^\top g_j/p$, the expectation of the product of two different cells is

$$\mathbb E(y_{im}y_{jn})=\sigma_G^2K^G_{ij}+\sum_k\sigma^2_{G\times C,k}K^G_{ij}c_{imk}c_{jnk}+\sigma_I^2\mathbf 1(i=j)+\sum_k\sigma^2_{C,k}c_{imk}c_{jnk}+\sum_l\sigma^2_{O,l}o_{iml}o_{jnl},$$

so the variances are the coefficients of a regression of $y_{im}y_{jn}$ on known kernels.
Cell pairs are never formed: the normal equations are assembled from donor-level
sufficient statistics (cost $O(N^2K^2)$ for $N$ donors, independent of the number of
cells) and the same quantities give leave-donor-block-out estimates for the jackknife.
`h2` reports the normalized heritability $\sigma_G^2/(\sigma_G^2+\sigma_{G\times C}^2+\sigma_C^2+\sigma_o^2+\sigma_I^2)$
that is comparable to pseudobulk estimates. Cell-type and perturbation models replace the context
kernel by "same cell type" and the genetic kernel by "same perturbation". Full derivation:
[Model](https://LidaWangPSU.github.io/scHINT/articles/model.html).

## Helpers

`make_grm()`, `split_snps()`, `read_plink()`, `select_cis_snps()`, `read_grm_bin()`,
`write_grm_bin()`, `simulate_geno()`, `simulate_schint()`, `simulate_pert()`.

Example data are fully simulated; no individual-level data are distributed.

## Citation

Wang L, Fadil C, Wang L, Van Dyke K, Gusev A. *Single-cell gene expression heritability
informed by gene × cell-context interactions.* (manuscript)
