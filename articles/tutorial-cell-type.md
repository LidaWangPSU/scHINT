# Tutorial: G×Cell-Type heritability (pseudobulk)

When expression is available as pseudobulk per donor and cell type,
[`schint_sample()`](https://LidaWangPSU.github.io/scHINT/reference/schint_sample.md)
asks how much genetic variance is **shared across cell types** (`G`) and
how much is **cell-type specific** (`GxC`).

``` r

head(ex$sample, 3)
#>   donor celltype        pc1        pc2 age sex batch     GENE_A    GENE_B
#> 1  D001      CT1 -0.1113267 -0.6001958  63   F    B1 -1.2743896 2.6097890
#> 2  D001      CT2 -0.2561393  0.6741177  63   F    B1 -0.4401739 0.8726951
#> 3  D001      CT3 -0.8792763  1.0525682  63   F    B1 -2.2811747 1.6666779
#>       GENE_C
#> 1  2.1917547
#> 2 -1.1131183
#> 3  0.1120468
table(ex$sample$celltype)
#> 
#> CT1 CT2 CT3 CT4 CT5 
#> 400 400 400 400 400
```

## Fit

``` r

fs <- schint_sample(ex$sample, y = "GENE_A", id = "donor", celltype = "celltype",
                    geno = ex$geno, covariates = c("age", "sex", "pc1", "pc2"),
                    jackknife = TRUE)
fs
#> scHINT model: 2000 observations, 400 donors, 1 GRM
#> 
#> Variance components (variance on the standardized-expression scale):
#>             term variance     h2     se  se_h2
#>                G   0.1020 0.1070 0.0237 0.0246
#>              GxC   0.0566 0.0595 0.0206 0.0202
#>          G_total   0.1580 0.1660 0.0295 0.0276
#>                I   0.0227 0.0238 0.0199 0.0209
#>          context   0.7320 0.7700 0.0216 0.0214
#>        covariate   0.0378 0.0398 0.0115 0.0116
#>  total_explained   0.9510     NA 0.0369     NA
#> 
#> Standard errors: jackknife over 400 donor blocks
```

- `G` - genetic variance shared by all cell types.
- `GxC` - genetic variance specific to a cell type (pairs from the same
  type).
- `I` - donor effect shared across cell types, not genetic.
- `context` - variance of cell-type mean differences.

## Shared-genetics model only

``` r

schint_sample(ex$sample, "GENE_A", "donor", "celltype", geno = ex$geno, gxc = FALSE)$summary
#>              term   variance   fraction        h2
#> 1               G 0.11261510 0.11261510 0.1293135
#> 2               I 0.02610727 0.02610727 0.0299784
#> 3         context 0.73214679 0.73214679 0.8407081
#> 4 total_explained 0.87086915 0.87086915        NA
```

## A separate genetic variance per cell type

``` r

fl <- schint_sample(ex$sample, "GENE_A", "donor", "celltype", geno = ex$geno,
                    cat_mode = "per_level", jackknife = TRUE)
subset(fl$coefficients, group == "GxC", select = c(component, variance, se_variance))
#>         component     variance se_variance
#> 5  G:celltype=CT1  0.015164751 0.009833371
#> 7  G:celltype=CT2  0.013769729 0.008339680
#> 9  G:celltype=CT3  0.017481713 0.008585974
#> 11 G:celltype=CT4  0.016568279 0.013131599
#> 13 G:celltype=CT5 -0.006398207 0.015131072
```

## Genes without cell-type specificity

``` r

schint_sample(ex$sample, "GENE_B", "donor", "celltype", geno = ex$geno,
              jackknife = TRUE)$summary[1:3, c("term", "variance", "se")]
#>      term     variance         se
#> 1       G 0.1282184619 0.03466213
#> 2     GxC 0.0005023073 0.01910134
#> 3 G_total 0.1287207692 0.03934649
```

## Notes

- The data should contain one row per donor × cell type. Several rows
  per pair trigger a warning (use
  [`schint_cell()`](https://LidaWangPSU.github.io/scHINT/reference/schint_cell.md)
  for cell-level data).
- Donors need not have all cell types; missing pairs are simply absent.
- Pseudobulk values should be normalized and, if desired, residualized
  for donor-level covariates before calling, or the covariates passed in
  `covariates`.
