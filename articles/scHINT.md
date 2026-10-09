# Get started with scHINT

## 1. Introduction

**scHINT** (single-cell gene expression **H**eritability with
gene-by-cell-context **INT**eractions) quantifies how much of the
variation in a gene’s expression is explained by genetics, and how much
of that genetic effect **depends on the cellular context**. Most
expression-QTL and heritability analyses use pseudobulk profiles, which
average over cells and so cannot see genetic effects that change along a
continuous cell state. scHINT works directly on single-cell data.

**The main model is the cell-level G×Cell-State heritability model**
([`schint_cell()`](https://LidaWangPSU.github.io/scHINT/reference/schint_cell.md)):
for each gene, the variance of single-cell expression is partitioned
into

- a **main genetic** component (G), shared by all cells,
- a **G×Cell-State** component: genetic effects that vary continuously
  along cell-state axes such as principal components or pseudotime,
- a donor (individual) component, cell-state and covariate components,
  and cell-level noise.

Estimation is a Haseman-Elston (HE) regression on pairs of cells that is
assembled from donor-level sufficient statistics, so it scales to
millions of cells, and standard errors come from a donor-block jackknife
at almost no extra cost.

The same framework extends to other contexts:

| model | function | what it estimates |
|----|----|----|
| **G×Cell-State** (cell level) | [`schint_cell()`](https://LidaWangPSU.github.io/scHINT/reference/schint_cell.md) | main genetic and genotype × cell-state heritability from single-cell data |
| **G×Cell-Type** (pseudobulk) | [`schint_sample()`](https://LidaWangPSU.github.io/scHINT/reference/schint_sample.md) | cell-type-shared and cell-type-specific genetic variance from donor × cell-type expression |
| **Perturbation × context** | [`schint_pert()`](https://LidaWangPSU.github.io/scHINT/reference/schint_pert.md) | perturbation and perturbation × cell-state / cell-type variance in perturb-seq screens |

## 2. Models

Notation follows the manuscript (Wang *et al.*). Full derivations,
including the sufficient-statistics algorithm and the jackknife, are on
the [Model
page](https://LidaWangPSU.github.io/scHINT/articles/model.html).

### 2.1 G×Cell-State heritability (cell level)

For a gene, let \\y\_{im}\\ be the normalized expression in cell \\m\\
(\\m=1,\dots,M_i\\) of individual \\i\\ (\\i=1,\dots,n\\),
\\g_i\in\mathbb{R}^p\\ the standardized cis-genotypes (variants within
±500 kb of the gene), \\c\_{im}\in\mathbb{R}^K\\ the cell-state features
(e.g. the leading within-cell-type principal components) and \\o\_{im}\\
further covariates:

\\ y\_{im}=g_i^\top\beta+(g_i\otimes
c\_{im})^\top\gamma+c\_{im}^\top\alpha+o\_{im}^\top\eta+e_i+\varepsilon\_{im},
\\

\\ \beta\sim N\\\left(0,\sigma_G^2p^{-1}I_p\right),\quad \gamma_k\sim
N\\\left(0,\sigma\_{G\times C,k}^2p^{-1}I_p\right),\quad e_i\sim
N(0,\sigma_I^2),\quad \varepsilon\_{im}\sim N(0,\sigma\_\varepsilon^2).
\\

\\\sigma_G^2\\ is the contribution of main cis-genetic effects and
\\\sigma\_{G\times C,k}^2\\ that of genetic effects varying along
cell-state axis \\k\\. With \\K^G\_{ij}=g_i^\top g_j/p\\ the cis-GRM,
the expected product of two different cells is

\\ \mathbb{E}(y\_{im}y\_{jn})=\sigma_G^2K^G\_{ij}+\sum\_{k=1}^{K\_{\rm
int}}\sigma\_{G\times C,k}^2K^G\_{ij}c\_{imk}c\_{jnk} +\sigma_I^2\mathbf
1(i=j)+\sum\_{k=1}^{K}\sigma\_{C,k}^2c\_{imk}c\_{jnk}+\sum\_{l=1}^{L}\sigma\_{O,l}^2o\_{iml}o\_{jnl}.
\\

scHINT regresses \\y\_{im}y\_{jn}\\ on these kernels. The baseline model

\\ y\_{im}y\_{jn}\sim 1+K^G\_{ij}+\mathbf
1(i=j)+\sum\_{k=1}^{K}c\_{imk}c\_{jnk}+\sum\_{l=1}^{L}o\_{iml}o\_{jnl}
\\

estimates the main genetic component, and the interaction model adds
\\\sum\_{k=1}^{K\_{\rm int}}K^G\_{ij}c\_{imk}c\_{jnk}\\; the
coefficients of \\K^G\_{ij}\\ and \\K^G\_{ij}c\_{imk}c\_{jnk}\\ are
\\\hat\sigma_G^2\\ and \\\hat\sigma\_{G\times C,k}^2\\, and
\\\hat\sigma\_{G\times C}^2=\sum_k\hat\sigma\_{G\times C,k}^2\\.
Normalized heritabilities, comparable to pseudobulk analyses, exclude
the cell-level residual from the denominator:

\\ h_G^2=\frac{\sigma_G^2}{\sigma_G^2+\sigma\_{G\times
C}^2+\sigma_C^2+\sigma_o^2+\sigma_I^2},\qquad h\_{G\times
C}^2=\frac{\sigma\_{G\times C}^2}{\sigma_G^2+\sigma\_{G\times
C}^2+\sigma_C^2+\sigma_o^2+\sigma_I^2}. \\

Several GRMs (e.g. MAF bins) give one \\\sigma_G^2\\ and one
\\\sigma\_{G\times C}^2\\ per GRM.

### 2.2 G×Cell-Type heritability (pseudobulk)

For the pseudobulk expression \\y\_{it}\\ of individual \\i\\ in cell
type \\t=1,\dots,T\\ (Supplementary Note 1):

\\
y\_{it}=g_i^\top\beta+g_i^\top\gamma_t+\alpha_t+o\_{it}^\top\eta+e_i+\varepsilon\_{it},
\\

\\ \beta\sim N\\\left(0,\sigma_G^2p^{-1}I_p\right),\quad \gamma_t\sim
N\\\left(0,\sigma\_{G\times CT}^2p^{-1}I_p\right),\quad \alpha_t\sim
N(0,\sigma_T^2),\quad e_i\sim N(0,\sigma_I^2), \\

where \\\beta\\ is the genetic effect shared across cell types,
\\\gamma_t\\ the cell-type-specific deviation, \\\alpha_t\\ the
cell-type mean and \\e_i\\ a donor effect shared across cell types. For
two observations \\(i,s)\neq(j,t)\\

\\ \mathbb{E}(y\_{is}y\_{jt})=\sigma_G^2K^G\_{ij}+\sigma\_{G\times
CT}^2K^G\_{ij}\mathbf 1(s=t) +\sigma\_{CT}^2\mathbf
1(s=t)+\sigma_I^2\mathbf 1(i=j), \\

and the HE regression uses the kernels \\K^G\_{ij}\\, \\K^G\_{ij}\mathbf
1(s=t)\\, \\\mathbf 1(s=t)\\ and \\\mathbf 1(i=j)\\ (plus covariate
kernels). \\\sigma_G^2\\ is genetic variance shared across cell types
and \\\sigma\_{G\times CT}^2\\ is cell-type-specific genetic variance.

### 2.3 Cell-level and pseudobulk heritability

Pseudobulk averages out cell-specific variation, so cell-level and
pseudobulk heritability have different denominators (Supplementary Note
2). With \\m\\ cells per individual and
\\q\_{ik}=\alpha_1X_i+\alpha_2Z\_{ik}+\beta G_i+I_i+\varepsilon\_{ik}\\,

\\ h^2\_{\rm
cell}=\frac{\beta^2}{\alpha_1^2+\alpha_2^2+\beta^2+\sigma_I^2+\sigma\_\varepsilon^2},\qquad
h^2\_{\rm
PB}=\frac{\beta^2}{\alpha_1^2+\beta^2+\sigma_I^2+(\alpha_2^2+\sigma\_\varepsilon^2)/m}.
\\

scHINT reports `h2` (residual \\\sigma\_\varepsilon^2\\ removed from the
denominator, a lower bound of the pseudobulk value) and `h2_upper` (also
removing the cell-level covariate variance \\\alpha_2^2\\).

### 2.4 Perturbation × context heritability

For a gene, let \\p_m\\ be the perturbation of cell \\m\\ and
\\c_m\in\mathbb{R}^K\\ its standardized cell-state features:

\\
y_m=b\_{p_m}+\sum\_{k=1}^{K}c\_{mk}\alpha_k+\sum\_{k=1}^{K}c\_{mk}u\_{p_mk}+\varepsilon_m,
\\

\\ b\_{p}\sim N(0,\sigma_P^2),\quad \alpha_k\sim N(0,\sigma\_{{\rm
PC},k}^2),\quad u\_{pk}\sim N(0,\sigma\_{P\times{\rm PC},k}^2),\quad
\varepsilon_m\sim N(0,\sigma\_\varepsilon^2). \\

\\\sigma_P^2\\ is perturbation-associated variation shared across cell
states and \\\sigma\_{P\times{\rm PC},k}^2\\ the variance of the
perturbation response along axis \\k\\. For two cells \\m\neq n\\ the HE
regression is

\\ y_my_n\sim 1+\mathbf 1(p_m=p_n)+\sum\_{k=1}^{K}c\_{mk}c\_{nk}+\mathbf
1(p_m=p_n)\sum\_{k=1}^{K}c\_{mk}c\_{nk}. \\

With a categorical context (cell type) \\c\_{mk}c\_{nk}\\ is replaced by
\\\mathbf 1(t_m=t_n)\\. Because genes measured in the same cells share
their pair predictors, the design cross-products are computed once for
all genes.

## 3. Install

``` r

install.packages("remotes")
remotes::install_github("LidaWangPSU/scHINT")
```

scHINT needs only base R. (`BEDMatrix` is optional, to read PLINK files
with
[`read_plink()`](https://LidaWangPSU.github.io/scHINT/reference/read_plink.md).)
The examples below use the simulated data shipped with the package:

``` r

data(schint_example)        # genotypes + cell-level and sample-level expression
ex <- schint_example
```

## 4. Usage

### 4.1 G×Cell-State heritability: `schint_cell()`

**Input**

- `data` – one row per cell (of one cell type) with expression, donor
  ID, cell-state variables and covariates.
- `y`, `id` – names of the expression and donor-ID columns.
- `geno` – donors × cis-SNPs dosage matrix (row names = donor IDs),
  **or** `grm` – a pre-computed GRM (or list of GRMs). `n_grm`/`grm_by`
  split the SNPs into several GRMs (default: one GRM from all cis-SNPs).
- `context` – cell-state column(s) that interact with genetics (numeric
  = continuous, factor = categorical).
- `covariates` – other covariates (main effects only).
  `jackknife = TRUE` for standard errors.

``` r

head(ex$cell[, c("donor", "state", "pc1", "pc2", "age", "sex", "batch", "GENE_A")], 3)
#>   donor       state        pc1        pc2 age sex batch     GENE_A
#> 1  D001 -0.34413687 -0.6334822 0.02127604  63   F    B1 -1.1423329
#> 2  D001  0.04947898 -0.7615181 0.84821092  63   F    B1 -0.3653666
#> 3  D001 -0.84762422  1.2939331 0.24429519  63   F    B1  0.7018750
dim(ex$geno)
#> [1] 400 250
```

**Example**

``` r

fit <- schint_cell(ex$cell, y = "GENE_A", id = "donor", geno = ex$geno,
                   context = "state",
                   covariates = c("pc1", "pc2", "age", "sex", "batch"),
                   jackknife = TRUE, n_blocks = 50)
fit
#> scHINT model: 12000 observations, 400 donors, 1 GRM
#> 
#> Variance components (variance on the standardized-expression scale):
#>             term variance     h2 h2_upper      se   se_h2 se_h2_upper
#>                G  0.20500 0.4140    0.438 0.04700 0.08190      0.0860
#>              GxC  0.16300 0.3290    0.348 0.03890 0.04980      0.0526
#>          G_total  0.36900 0.7430    0.787 0.06680 0.07080      0.0732
#>                I  0.07670 0.1550    0.164 0.03580 0.07430      0.0789
#>          context  0.00704 0.0142       NA 0.00331 0.00617          NA
#>        covariate  0.04360 0.0879    0.093 0.01130 0.02600      0.0278
#>  total_explained  0.49600     NA       NA 0.06100      NA          NA
#> 
#> Standard errors: jackknife over 50 donor blocks
```

**Output** – `fit$summary` (grouped components) and `fit$coefficients`
(one row per regression term): `G` main genetic variance, `GxC`
G×Cell-State variance, `G_total = G + GxC`, `I` donor, `context` and
`covariate` components; `h2`/`h2_upper` normalized heritabilities;
jackknife `se*` columns. Several GRMs and several context variables add
per-GRM and per-variable rows.

``` r

# three MAF-binned GRMs, two cell-state axes
schint_cell(ex$cell, "GENE_A", "donor", geno = ex$geno, n_grm = 3, grm_by = "maf",
            context = c("state", "subtype"), covariates = c("pc1", "pc2"))$summary[, 1:3]
#>               term    variance    fraction
#> 1                G  0.21304077  0.21304077
#> 2                1  0.06190295  0.06190295
#> 3                2  0.05429398  0.05429398
#> 4                3  0.09684385  0.09684385
#> 5              GxC  0.15004511  0.15004511
#> 6              1:C  0.05324709  0.05324709
#> 7              2:C  0.04827433  0.04827433
#> 8              3:C  0.04852369  0.04852369
#> 9          G_total  0.36308588  0.36308588
#> 10       GxC:state  0.16876151  0.16876151
#> 11     GxC:subtype -0.01871639 -0.01871639
#> 12               I  0.09303695  0.09303695
#> 13         context  0.02630493  0.02630493
#> 14       covariate  0.02027540  0.02027540
#> 15 total_explained  0.50270316  0.50270316
```

### 4.2 G×Cell-Type heritability (pseudobulk): `schint_sample()`

**Input** – `data` with one row per donor × cell type (expression, donor
ID, cell-type column, optional covariates); the same genetic input as
above; `celltype` names the cell-type column.

``` r

head(ex$sample[, c("donor", "celltype", "age", "sex", "GENE_A")], 3)
#>   donor celltype age sex     GENE_A
#> 1  D001      CT1  63   F -1.2743896
#> 2  D001      CT2  63   F -0.4401739
#> 3  D001      CT3  63   F -2.2811747
```

**Example**

``` r

fs <- schint_sample(ex$sample, y = "GENE_A", id = "donor", celltype = "celltype",
                    geno = ex$geno, covariates = c("age", "sex"), jackknife = TRUE)
fs
#> scHINT model: 2000 observations, 400 donors, 1 GRM
#> 
#> Variance components (variance on the standardized-expression scale):
#>             term variance     h2     se  se_h2
#>                G   0.1020 0.1080 0.0237 0.0248
#>              GxC   0.0566 0.0601 0.0206 0.0204
#>          G_total   0.1580 0.1680 0.0295 0.0279
#>                I   0.0228 0.0243 0.0199 0.0211
#>          context   0.7320 0.7780 0.0216 0.0215
#>        covariate   0.0278 0.0295 0.0103 0.0106
#>  total_explained   0.9410     NA 0.0362     NA
#> 
#> Standard errors: jackknife over 400 donor blocks
```

**Output** – same structure as above: `G` is genetic variance shared
across cell types, `GxC` the cell-type-specific genetic variance, `I`
the donor component shared across cell types, `context` the cell-type
mean variance. `gxc = FALSE` fits shared genetics only;
`cat_mode = "per_level"` estimates a genetic variance for each cell
type.

### 4.3 Perturbation × context heritability: `schint_pert()`

**Input** – `data` with one row per cell, one or several expression
columns (`y`), the `perturb` column, `context` (cell state and/or cell
type), `covariates`, `control` labels to drop (e.g. `"NT"`) and
optionally `perturb_group`.

``` r

data(schint_pert_example)
head(schint_pert_example[, c("perturb", "perturb_group", "celltype", "state", "GENE_A")], 3)
#>   perturb perturb_group celltype     state     GENE_A
#> 1  PERT01        group1      CT2 1.1188078 0.07107567
#> 2  PERT01        group1      CT3 1.1119961 0.14960499
#> 3  PERT01        group1      CT1 0.4430754 0.38410583
```

**Example**

``` r

fp <- schint_pert(schint_pert_example, y = c("GENE_A", "GENE_B"), perturb = "perturb",
                  context = c("state", "celltype"), covariates = c("pc1", "pc2"),
                  control = "NT", jackknife = TRUE)
fp
#> scHINT perturbation model: 3470 cells, 60 perturbations, 2 genes
#> 
#>    gene            term variance      se
#>  GENE_A               P  0.02720 0.00866
#>  GENE_A             PxC  0.06180 0.01930
#>  GENE_A         P_total  0.08900 0.02080
#>  GENE_A       PxC:state  0.05070 0.01460
#>  GENE_A    PxC:celltype  0.01110 0.01130
#>  GENE_A         context  0.09800 0.01260
#>  GENE_A       covariate  0.04050 0.00617
#>  GENE_A total_explained  0.22700 0.02510
#>  GENE_B               P  0.04350 0.01510
#>  GENE_B             PxC  0.02450 0.01220
#>  GENE_B         P_total  0.06800 0.01880
#>  GENE_B       PxC:state -0.00277 0.00293
#>  GENE_B    PxC:celltype  0.02730 0.01260
#>  GENE_B         context  0.06640 0.01250
#>  GENE_B       covariate  0.05610 0.00832
#>  GENE_B total_explained  0.19100 0.02470
#> 
#> Standard errors: jackknife over 60 perturbation blocks
```

**Output** – long-format `fp$summary` / `fp$coefficients` with one block
per gene: `P` perturbation variance shared across states, `PxC`
perturbation × context variance (split into `PxC:state` and
`PxC:celltype`), `P_total`, and per-group rows when `perturb_group` is
used.

### More

[Users
manual](https://LidaWangPSU.github.io/scHINT/articles/users-manual.html)
· [Input formats and
outputs](https://LidaWangPSU.github.io/scHINT/articles/inputs-outputs.html)
·
[Tutorials](https://LidaWangPSU.github.io/scHINT/articles/tutorial-cell-state.html)
·
[Troubleshooting](https://LidaWangPSU.github.io/scHINT/articles/troubleshooting.html)

**Citation.** Wang L, Fadil C, Wang L, Van Dyke K, Gusev A. *Single-cell
gene expression heritability informed by gene × cell-context
interactions.* (manuscript)
