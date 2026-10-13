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
#>             term variance     se
#>                G   0.1020 0.0237
#>              GxC   0.0566 0.0206
#>          G_total   0.1580 0.0295
#>                I   0.0227 0.0199
#>          context   0.7320 0.0216
#>        covariate   0.0378 0.0115
#>  total_explained   0.9510 0.0369
#> 
#> Standard errors: jackknife over 400 donor blocks
```

- `G` - genetic variance shared by all cell types.
- `GxC` - genetic variance specific to a cell type (pairs from the same
  type).
- `I` - donor effect shared across cell types, not genetic.
- `context` - variance of cell-type mean differences.

The variance components of pseudobulk expression are already at the
population level (the cell-level noise is averaged out when cells are
pooled), so no further normalization is applied and the output has no
`h2` column.

## Notes

- The data should contain one row per donor × cell type. Several rows
  per pair trigger a warning (use
  [`schint_cell()`](https://LidaWangPSU.github.io/scHINT/reference/schint_cell.md)
  for cell-level data).
- Donors need not have all cell types; missing pairs are simply absent.
- Pseudobulk values should be normalized and, if desired, residualized
  for donor-level covariates before calling, or the covariates passed in
  `covariates`.
