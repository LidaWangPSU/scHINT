# Cell-level GxCell-State heritability model

Haseman-Elston regression on pairs of cells that partitions the variance
of one gene's expression into a main genetic component (\`G\`, one per
GRM), genetic effects that change along a cell-state/context variable
(\`G x context\`), a shared-donor component (\`I\`), context main
effects and other covariates.

## Usage

``` r
schint_cell(
  data,
  y,
  id,
  geno = NULL,
  grm = NULL,
  n_grm = 1L,
  grm_by = c("maf", "position", "random"),
  context = NULL,
  covariates = NULL,
  cat_mode = c("pooled", "per_level"),
  ind_effect = TRUE,
  context_main = TRUE,
  jackknife = FALSE,
  n_blocks = NULL,
  scale_y = TRUE,
  scale_x = TRUE,
  seed = 1
)
```

## Arguments

- data:

  Data frame with one row per cell.

- y:

  Name of the (numeric) expression column.

- id:

  Name of the donor ID column. IDs must match the row names of
  \`geno\`/\`grm\`.

- geno:

  Genotypes of the (cis) SNPs: a donors x SNPs dosage matrix with row
  names, or a list of such matrices (one GRM each). Give either \`geno\`
  or \`grm\`.

- grm:

  A pre-computed donor x donor GRM with donor IDs as dimnames, or a
  (named) list of GRMs for a multi-GRM model.

- n_grm:

  Number of GRMs to build from a single \`geno\` matrix (default \`1\` =
  one GRM from all SNPs). Values \> 1 split the SNPs with
  \[split_snps()\].

- grm_by:

  How to split SNPs when \`n_grm \> 1\`: \`"maf"\`, \`"position"\` or
  \`"random"\`.

- context:

  Character vector of column names used as the cell context that
  interacts with genetics. Numeric columns are continuous contexts (e.g.
  a pseudotime or PC); factor/character columns are categorical.
  \`NULL\` fits a model without genotype-by-context terms.

- covariates:

  Character vector of other covariates (main effects only). Variables
  also listed in \`context\` are not duplicated.

- cat_mode:

  For categorical variables in \`context\`: \`"pooled"\` fits one shared
  component (pairs from the same level) and \`"per_level"\` fits one per
  level.

- ind_effect:

  Include the shared-donor component \`I\` (default \`TRUE\`).

- context_main:

  Include main-effect components for each context variable (default
  \`TRUE\`).

- jackknife:

  If \`TRUE\`, delete-donor-block jackknife standard errors.

- n_blocks:

  Number of jackknife blocks over donors (default: one per donor).

- scale_y, scale_x:

  Standardize \`y\` / continuous context and covariates (default
  \`TRUE\`).

- seed:

  Seed for assigning donors to jackknife blocks.

## Value

An object of class \`"schint"\` with \`coefficients\`, \`summary\`
(grouped variance components), the jackknife replicates and the
sufficient statistics \`XtX\`, \`Xty\`.

## Details

For a pair of cells \`a\`, \`b\` from donors \`i\`, \`j\` the model
regresses \`y_a \* y_b\` on \`1\`, \`GRM\[i,j\]\`, \`1(i == j)\`, \`u_a
\* u_b\` (context and covariate main effects) and \`GRM\[i,j\] \* u_a \*
u_b\` (genotype-by-context). The estimate for each component is its
variance contribution when \`y\` and the covariates are standardized.

## Examples

``` r
data(schint_example)
d <- schint_example$cell
fit <- schint_cell(d, y = "GENE_A", id = "donor",
                   geno = schint_example$geno,
                   context = "state", covariates = c("pc1", "pc2", "age", "sex", "batch"))
fit
#> scHINT model: 12000 observations, 400 donors, 1 GRM
#> 
#> Variance components (fraction of Var(y)):
#>             term variance fraction
#>                G  0.20500  0.20500
#>              GxC  0.16300  0.16300
#>          G_total  0.36900  0.36900
#>                I  0.07670  0.07670
#>          context  0.00704  0.00704
#>        covariate  0.04360  0.04360
#>  total_explained  0.49600  0.49600
```
