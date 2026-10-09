# scHINT

[![R-CMD-check](https://github.com/LidaWangPSU/scHINT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LidaWangPSU/scHINT/actions)

Documentation: <https://LidaWangPSU.github.io/scHINT/>

Haseman-Elston regression estimators of the heritability of gene expression
across cell states and cell types, with genotype-by-context interactions,
multiple GRMs, extra covariates and jackknife standard errors.

## Install

```r
# install.packages("remotes")
remotes::install_github("LidaWangPSU/scHINT")
```

## Main functions

| function | data | genotype-by-context term |
|---|---|---|
| `schint_cell()` | one row per cell | G x any continuous or categorical context (Cell state, etc.) |
| `schint_sample()` | one row per donor x cell type | G x cell type (pooled or per cell type) |
| `schint_pert()` | one row per cell, perturb-seq | perturbation x cell state and/or cell type (several genes at once, optional perturbation groups) |

Common options: `geno` (dosage matrix of the cis SNPs) **or** `grm` (pre-computed
GRM(s)); `n_grm` / `grm_by` to split SNPs into several GRMs (default: one GRM from
all cis-SNPs); `covariates`; `jackknife = TRUE` for delete-donor(-block) standard errors.

```r
library(scHINT)
data(schint_example)
ex <- schint_example

# cell level: G x continuous cell state, 3 GRMs split by MAF, jackknife SEs
fit <- schint_cell(ex$cell, y = "GENE_A", id = "donor", geno = ex$geno,
                   n_grm = 3, grm_by = "maf",
                   context = "state", covariates = c("pc1", "pc2", "age", "sex", "batch"),
                   jackknife = TRUE)
fit
fit$coefficients

# sample level: G x cell type, pre-computed GRM
fit2 <- schint_sample(ex$sample, y = "GENE_A", id = "donor", celltype = "celltype",
                      grm = make_grm(ex$geno), covariates = c("age", "sex"))
```

```r
data(schint_pert_example)
schint_pert(schint_pert_example, y = c("GENE_A", "GENE_B"), perturb = "perturb",
            context = c("state", "celltype"), covariates = c("pc1", "pc2"),
            control = "NT", jackknife = TRUE)
```

## Model

For two observations *a*, *b* from donors *i*, *j*, the regression is

    y_a y_b = b0 + sum_m s_m GRM_m[i,j] + t 1(i=j)
         + sum_v c_v u_va u_vb + sum_m sum_v g_mv GRM_m[i,j] u_va u_vb + e

where *u* are context variables (a factor contributes one component for pairs in
the same level) and covariates enter through main-effect terms. Estimates are
variance components on the scale of standardized `y`. Normal equations are built
from donor-aggregated sums, so the cost depends on the number of donors, not cells.

## Helpers
`make_grm()`, `split_snps()`, `read_plink()`, `select_cis_snps()`, `read_grm_bin()`,
`write_grm_bin()`, `simulate_schint()`, `simulate_pert()`.

For `schint_pert()` the donor kernel is replaced by "same perturbation".

Example data are fully simulated (`simulate_geno()`, `simulate_schint()`, `simulate_pert()`); no individual-level data are distributed.
