# Input formats and outputs

## Inputs

scHINT fits **one gene at a time**: you supply that gene’s expression,
cell/sample metadata and cis-genetic information.

### Expression and metadata (`data`)

A data frame (or tibble) with one row per observation.

| model | one row is | required columns |
|----|----|----|
| [`schint_cell()`](https://LidaWangPSU.github.io/scHINT/reference/schint_cell.md) | a cell (of one cell type) | expression `y`, donor `id`, optional `context`, `covariates` |
| [`schint_sample()`](https://LidaWangPSU.github.io/scHINT/reference/schint_sample.md) | a donor × cell-type pseudobulk | expression `y`, donor `id`, `celltype`, optional `covariates` |
| [`schint_pert()`](https://LidaWangPSU.github.io/scHINT/reference/schint_pert.md) | a cell | expression `y` (one or several genes), `perturb`, optional `context` (one cell-state or cell-type column), `covariates` |

``` r

head(ex$cell, 3)
#>     cell_id donor       state subtype        pc1        pc2 age sex batch
#> 1 cell00001  D001 -0.34413687       A -0.6334822 0.02127604  63   F    B1
#> 2 cell00002  D001  0.04947898       A -0.7615181 0.84821092  63   F    B1
#> 3 cell00003  D001 -0.84762422       A  1.2939331 0.24429519  63   F    B1
#>       GENE_A     GENE_B     GENE_C
#> 1 -1.1423329 -1.4138059 -0.6008276
#> 2 -0.3653666 -1.3328426 -1.2526841
#> 3  0.7018750 -0.5391143 -0.5859788
```

- **Expression** should be a normalized (e.g. log-normalized,
  depth-corrected) numeric value; it is centred and scaled internally.
  Non-finite values drop the row.
- **Context** columns: numeric = continuous (e.g. the leading principal
  components computed *within* a cell type, pseudotime, an activation
  score); factor/character = categorical (cell type, subtype,
  treatment). Numeric context is standardized.
- **Covariates** enter as main effects only. Use cell-level covariates
  (library size, PCs of expression, batch) and donor-level covariates
  (age, sex, genetic PCs). Categorical covariates are one-hot encoded as
  a “same level” component.
- Donor IDs must match the genetic input exactly (character comparison);
  donors without genotypes and rows with missing values are dropped with
  a message.

### Genetic input (`geno` **or** `grm`)

| argument | meaning |
|----|----|
| `geno` | donors × SNPs dosage matrix (0/1/2, `NA` allowed), **row names = donor IDs** – the cis SNPs of the gene; one GRM is built from all SNPs |
| `grm` | a pre-computed donor × donor GRM with donor IDs as dimnames, or a (named) list of GRMs for a multi-GRM model |

With `geno` one GRM is built from all supplied cis SNPs; for several
GRMs build them yourself (for example
`lapply(split_snps(geno, 3), make_grm)`) and pass the list to `grm`.
Select SNPs within ±500 kb of the gene
([`select_cis_snps()`](https://LidaWangPSU.github.io/scHINT/reference/select_cis_snps.md));
genotypes are standardized by allele frequency and the GRM is
\\XX^\top/p\\. GRMs computed by GCTA can be loaded with
[`read_grm_bin()`](https://LidaWangPSU.github.io/scHINT/reference/read_grm_bin.md).
PLINK files are read with
[`read_plink()`](https://LidaWangPSU.github.io/scHINT/reference/read_plink.md)
(needs the `BEDMatrix` package).

### Perturbation input

`perturb` is a column with the perturbation (target) of each cell.
Remove non-targeting controls with `control = "NT"` (they would
otherwise form a large “perturbation”); `min_cells` drops perturbations
with too few cells.

## Outputs

### `schint_cell()` and `schint_sample()`

Both return an object of class `"schint"`:

``` r

fit <- schint_cell(ex$cell, "GENE_A", "donor", geno = ex$geno, context = "state",
                   covariates = c("pc1", "pc2", "age", "sex", "batch"),
                   jackknife = TRUE, n_blocks = 50)
names(fit)
#>  [1] "coefficients" "summary"      "n_obs"        "n_donors"     "n_grm"       
#>  [6] "grm_names"    "context"      "covariates"   "cat_mode"     "XtX"         
#> [11] "Xty"          "jackknife"    "call"
```

**`fit$summary`** – grouped variance components (what you usually
report):

``` r

fit$summary
#>              term    variance         h2          se       se_h2
#> 1               G 0.205483423 0.41416751 0.047000781 0.081921231
#> 2             GxC 0.163326375 0.32919676 0.038891575 0.049752009
#> 3         G_total 0.368809799 0.74336427 0.066811617 0.070786682
#> 4               I 0.076689304 0.15457314 0.035847480 0.074297698
#> 5         context 0.007037157 0.01418393 0.003309055 0.006174124
#> 6       covariate 0.043599773 0.08787867 0.011296490 0.026030513
#> 7 total_explained 0.496136034         NA 0.060999075          NA
```

| column | meaning |
|----|----|
| `term` | `G`, `GxC` (sum over context variables), `G_total` = G + GxC, `I`, `context`, `covariate`, `total_explained`; with several GRMs also one row per GRM (`G1`, `G1:C`, …); with several context variables also `GxC:<variable>` |
| `variance` | variance contribution on the standardized-expression scale |
| `h2` | population-level heritability (cell-level residual removed from the denominator; see [Model](https://LidaWangPSU.github.io/scHINT/articles/model.md)) |
| `se`, `se_h2` | jackknife standard errors (only with `jackknife = TRUE`) |

**`fit$coefficients`** – one row per HE-regression term:

``` r

fit$coefficients
#>    component     group  grm     estimate    variance          se se_variance
#> 1  Intercept intercept <NA> -0.008875315          NA 0.005004164          NA
#> 2          G         G    G  0.206576226 0.205483423 0.047250741 0.047000781
#> 3          I         I <NA>  0.076689304 0.076689304 0.035847480 0.035847480
#> 4      state   context <NA>  0.007037744 0.007037157 0.003309331 0.003309055
#> 5    G:state       GxC    G  0.163913313 0.163326375 0.039031338 0.038891575
#> 6        pc1 covariate <NA>  0.010251003 0.010250149 0.002065749 0.002065576
#> 7        pc2 covariate <NA>  0.010026416 0.010025581 0.001797230 0.001797081
#> 8        age covariate <NA>  0.002508032 0.002507823 0.002897547 0.002897306
#> 9        sex covariate <NA>  0.013092156 0.013092156 0.010342335 0.010342335
#> 10     batch covariate <NA>  0.007724064 0.007724064 0.006181338 0.006181338
#>             z            p
#> 1  -1.7735860 7.613165e-02
#> 2   4.3719151 1.231614e-05
#> 3   2.1393220 3.240960e-02
#> 4   2.1266363 3.345032e-02
#> 5   4.1995310 2.674684e-05
#> 6   4.9623676 6.963904e-07
#> 7   5.5788154 2.421621e-08
#> 8   0.8655706 3.867257e-01
#> 9   1.2658801 2.055560e-01
#> 10  1.2495781 2.114537e-01
```

`estimate` is the regression coefficient (variance component of the
standardized variable), `variance` = `estimate` × the average diagonal
of its kernel (what the term contributes to the variance of expression,
which is 1). With jackknife: `se`, `z`, `p` (normal approximation) and
`se_variance`. `Intercept` is the pair-mean and carries no variance.

Other elements: `n_obs`, `n_donors`, `XtX`/`Xty` (sufficient
statistics), `jackknife$estimates` (leave-block-out estimates, blocks ×
components), `call`.

### `schint_pert()`

A `"schint_pert"` object with long-format `coefficients` and `summary`
tables with an extra `gene` column (terms: `P`, `PxC`, `P_total`,
`context`, `covariate`, `total_explained`).

### Extracting results

``` r

coef(fit)                      # named vector of regression coefficients
#>    Intercept            G            I        state      G:state          pc1 
#> -0.008875315  0.206576226  0.076689304  0.007037744  0.163913313  0.010251003 
#>          pc2          age          sex        batch 
#>  0.010026416  0.002508032  0.013092156  0.007724064
subset(fit$summary, term %in% c("G", "GxC"), select = c(term, h2, se_h2))
#>   term        h2      se_h2
#> 1    G 0.4141675 0.08192123
#> 2  GxC 0.3291968 0.04975201
```
