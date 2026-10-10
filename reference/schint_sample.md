# Sample-level G×Cell-Type heritability model

The sample-level (pseudobulk) counterpart of \[schint_cell()\]: each row
is one donor x cell-type expression value. Variance is partitioned into
a genetic component shared across cell types (\`G\`), cell-type-specific
genetics (\`G×Cell-Type\`, pairs from the same cell type), a donor
component shared across cell types (\`I\`), a cell-type component
(\`celltype\`) and covariates.

## Usage

``` r
schint_sample(
  data,
  y,
  id,
  celltype,
  geno = NULL,
  grm = NULL,
  covariates = NULL,
  gxc = TRUE,
  cat_mode = c("pooled", "per_level"),
  ind_effect = TRUE,
  celltype_main = TRUE,
  jackknife = FALSE,
  n_blocks = NULL,
  scale_x = TRUE,
  seed = 1
)
```

## Arguments

- data:

  Data frame with one row per donor x cell type.

- y:

  Name of the (numeric) expression column.

- id:

  Name of the donor ID column. IDs must match the row names of
  \`geno\`/\`grm\`.

- celltype:

  Name of the cell-type column.

- geno:

  Donors x SNPs dosage matrix of the cis SNPs, with donor IDs as row
  names; a single GRM is built from all SNPs. Give either \`geno\` or
  \`grm\`.

- grm:

  A pre-computed donor x donor GRM with donor IDs as dimnames, or a
  (named) list of GRMs for a multi-GRM model (one main genetic and one
  genotype-by-context component per GRM).

- covariates:

  Other covariates (donor-level or observation-level).

- gxc:

  Fit cell-type-specific genetics (default \`TRUE\`); \`FALSE\` gives
  the shared-genetics-only model.

- cat_mode:

  \`"pooled"\` (one \`G×Cell-Type\` component, default) or
  \`"per_level"\` (a separate genetic variance for each cell type).

- ind_effect:

  Include the shared-donor component \`I\` (default \`TRUE\`).

- celltype_main:

  Include the cell-type main-effect component (default \`TRUE\`).

- jackknife:

  If \`TRUE\`, delete-donor-block jackknife standard errors.

- n_blocks:

  Number of jackknife blocks over donors (default: one per donor).

- scale_x:

  Standardize continuous context and covariates (default \`TRUE\`).
  Expression \`y\` is always standardized to mean 0 and variance 1, so
  variance components are on the scale of standardized expression.

- seed:

  Seed for assigning donors to jackknife blocks.

## Value

An object of class \`"schint"\` with \`coefficients\`, \`summary\`
(grouped variance components), the jackknife replicates and the
sufficient statistics \`XtX\`, \`Xty\`.

## Examples

``` r
data(schint_example)
fit <- schint_sample(schint_example$sample, y = "GENE_A", id = "donor",
                     celltype = "celltype", geno = schint_example$geno,
                     covariates = c("age", "sex", "pc1", "pc2"))
fit
#> scHINT model: 2000 observations, 400 donors, 1 GRM
#> 
#> Variance components (variance on the standardized-expression scale):
#>             term variance     h2
#>                G   0.1020 0.1070
#>              GxC   0.0566 0.0595
#>          G_total   0.1580 0.1660
#>                I   0.0227 0.0238
#>          context   0.7320 0.7700
#>        covariate   0.0378 0.0398
#>  total_explained   0.9510     NA
```
