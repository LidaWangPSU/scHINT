# Getting started with scHINT

scHINT estimates how much of the variance of a gene’s expression is
explained by genetics (or by a perturbation), and how much of that
depends on the **cell state** or **cell type**, using Haseman-Elston
regression on pairs of observations.

The example data are fully simulated: `schint_example` (genotypes, a
cell-level table and a sample-level table) and `schint_pert_example`
(perturb-seq-like cells).

``` r

data(schint_example)
ex <- schint_example
str(ex, max.level = 1)
#> List of 5
#>  $ geno  : int [1:400, 1:250] 2 0 0 2 1 0 2 2 0 2 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ snp   :'data.frame':  250 obs. of  3 variables:
#>  $ cell  :'data.frame':  12000 obs. of  12 variables:
#>  $ sample:'data.frame':  2000 obs. of  10 variables:
#>  $ truth :'data.frame':  3 obs. of  6 variables:
head(ex$cell[, 1:9], 3)
#>     cell_id donor       state subtype        pc1        pc2 age sex batch
#> 1 cell00001  D001 -0.34413687       A -0.6334822 0.02127604  63   F    B1
#> 2 cell00002  D001  0.04947898       A -0.7615181 0.84821092  63   F    B1
#> 3 cell00003  D001 -0.84762422       A  1.2939331 0.24429519  63   F    B1
ex$truth
#>     gene    g  gxc    i  cov noise
#> 1 GENE_A 0.20 0.15 0.10 0.05   0.5
#> 2 GENE_B 0.25 0.00 0.10 0.05   0.6
#> 3 GENE_C 0.00 0.00 0.15 0.05   0.8
```

## Cell-level data: G x cell state

`geno` is the donor x SNP dosage matrix of the cis-SNPs of the gene. By
default one GRM is built from all SNPs. `context` is the cell-state
variable the genetic effect may depend on (numeric = continuous,
factor/character = categorical); `covariates` only enter as main
effects.

``` r

fit <- schint_cell(ex$cell, y = "GENE_A", id = "donor", geno = ex$geno,
                   context = "state",
                   covariates = c("pc1", "pc2", "age", "sex", "batch", "subtype"),
                   jackknife = TRUE)
fit
#> scHINT model: 12000 observations, 400 donors, 1 GRM
#> 
#> Variance components (fraction of Var(y)):
#>             term variance fraction      se se_fraction
#>                G 0.205000 0.205000 0.04490     0.04490
#>              GxC 0.163000 0.163000 0.04200     0.04200
#>          G_total 0.369000 0.369000 0.06240     0.06240
#>                I 0.076700 0.076700 0.03480     0.03480
#>          context 0.000348 0.000348 0.00312     0.00312
#>        covariate 0.069600 0.069600 0.01440     0.01440
#>  total_explained 0.515000 0.515000 0.06190     0.06190
#> 
#> Standard errors: jackknife over 400 donor blocks
```

- `G` - genetic variance independent of the cell state.
- `GxC` - genetic variance that changes along the context (here
  `state`).
- `G_total` - `G + GxC`; `I` - donor effect that is not genetic;
  `context`, `covariate` - main effects.

The per-component table is `fit$coefficients`. Variance components are
on the scale of standardized expression, so they are fractions of the
variance. A gene without genotype-by-state effects:

``` r

schint_cell(ex$cell, "GENE_B", "donor", geno = ex$geno, context = "state")$summary
#>              term      variance      fraction
#> 1               G  0.2431935652  0.2431935652
#> 2             GxC -0.0008234294 -0.0008234294
#> 3         G_total  0.2423701358  0.2423701358
#> 4               I  0.1077121143  0.1077121143
#> 5         context  0.0064094800  0.0064094800
#> 6 total_explained  0.3564917301  0.3564917301
```

### Several GRMs, or pre-computed GRMs

``` r

# split cis-SNPs into three MAF bins -> three GRMs
fit3 <- schint_cell(ex$cell, "GENE_A", "donor", geno = ex$geno,
                    n_grm = 3, grm_by = "maf", context = "state")
fit3$summary
#>               term    variance    fraction
#> 1                G 0.203760119 0.203760119
#> 2                1 0.049718321 0.049718321
#> 3                2 0.056917926 0.056917926
#> 4                3 0.097123872 0.097123872
#> 5              GxC 0.163916325 0.163916325
#> 6              1:C 0.071363745 0.071363745
#> 7              2:C 0.044393299 0.044393299
#> 8              3:C 0.048159282 0.048159282
#> 9          G_total 0.367676445 0.367676445
#> 10               I 0.092974593 0.092974593
#> 11         context 0.007036457 0.007036457
#> 12 total_explained 0.467687494 0.467687494

# or supply GRMs yourself (e.g. from GCTA via read_grm_bin())
K <- make_grm(ex$geno)
fitK <- schint_cell(ex$cell, "GENE_A", "donor", grm = K, context = "state")
```

`geno` can also be a list of matrices (one GRM each, e.g. SNPs by
annotation).

## Sample-level data: G x cell type

``` r

fs <- schint_sample(ex$sample, y = "GENE_A", id = "donor", celltype = "celltype",
                    geno = ex$geno, covariates = c("age", "sex", "pc1", "pc2"),
                    jackknife = TRUE)
fs
#> scHINT model: 2000 observations, 400 donors, 1 GRM
#> 
#> Variance components (fraction of Var(y)):
#>             term variance fraction     se se_fraction
#>                G   0.1020   0.1020 0.0237      0.0237
#>              GxC   0.0566   0.0566 0.0206      0.0206
#>          G_total   0.1580   0.1580 0.0295      0.0295
#>                I   0.0227   0.0227 0.0199      0.0199
#>          context   0.7320   0.7320 0.0216      0.0216
#>        covariate   0.0378   0.0378 0.0115      0.0115
#>  total_explained   0.9510   0.9510 0.0369      0.0369
#> 
#> Standard errors: jackknife over 400 donor blocks
```

`gxc = FALSE` fits shared genetics only; `cat_mode = "per_level"`
estimates a separate genetic variance for every cell type.

## Perturbation heritability

[`schint_pert()`](https://LidaWangPSU.github.io/scHINT/reference/schint_pert.md)
replaces the genetic kernel by “same perturbation”. Several genes are
fitted at once and the interaction can be a cell state, a cell type, or
both.

``` r

data(schint_pert_example)
fp <- schint_pert(schint_pert_example, y = c("GENE_A", "GENE_B", "GENE_C"),
                  perturb = "perturb", context = c("state", "celltype"),
                  covariates = c("pc1", "pc2"), control = "NT")
fp$summary[fp$summary$term %in% c("P", "PxC:state", "PxC:celltype"), ]
#>            gene         term      variance      fraction
#> GENE_A.1 GENE_A            P  0.0272190898  0.0272190898
#> GENE_A.4 GENE_A    PxC:state  0.0507430857  0.0507430857
#> GENE_A.5 GENE_A PxC:celltype  0.0110774848  0.0110774848
#> GENE_B.1 GENE_B            P  0.0435317075  0.0435317075
#> GENE_B.4 GENE_B    PxC:state -0.0027662529 -0.0027662529
#> GENE_B.5 GENE_B PxC:celltype  0.0272631582  0.0272631582
#> GENE_C.1 GENE_C            P  0.0665561848  0.0665561848
#> GENE_C.4 GENE_C    PxC:state  0.0009753755  0.0009753755
#> GENE_C.5 GENE_C PxC:celltype  0.0069372804  0.0069372804
```

Use `perturb_group = "perturb_group"` to estimate separate components
for groups of perturbations.

## Using your own data

1.  Genotypes of the cis window of one gene
    ([`read_plink()`](https://LidaWangPSU.github.io/scHINT/reference/read_plink.md) +
    [`select_cis_snps()`](https://LidaWangPSU.github.io/scHINT/reference/select_cis_snps.md)),
    as a donors x SNPs matrix whose row names match the donor column of
    your expression table.
2.  One row per cell (or donor x cell type) with expression of that
    gene, the donor ID, the context variable(s) and covariates.
3.  Call
    [`schint_cell()`](https://LidaWangPSU.github.io/scHINT/reference/schint_cell.md)
    /
    [`schint_sample()`](https://LidaWangPSU.github.io/scHINT/reference/schint_sample.md)
    once per gene.

Standard errors use a delete-donor (or `n_blocks` blocks of donors)
jackknife.
