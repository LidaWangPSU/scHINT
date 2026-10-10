# Perturbation heritability with cell-state or cell-type interaction

The perturb-seq analogue of \[schint_cell()\]: the donor/genotype kernel
is replaced by "same perturbation" (cells carrying the same perturbation
are correlated, cells with different perturbations are not).
Haseman-Elston regression on pairs of cells partitions the variance of
each gene into a perturbation component \`P\`, a perturbation-by-context
component \`P x context\` (the context is one variable at a time: a
continuous cell state or a categorical cell type), context main effects
and covariates.

## Usage

``` r
schint_pert(
  data,
  y,
  perturb,
  context = NULL,
  covariates = NULL,
  perturb_group = NULL,
  control = NULL,
  min_cells = 2,
  cat_mode = c("pooled", "per_level"),
  context_main = TRUE,
  jackknife = FALSE,
  n_blocks = NULL,
  scale_x = TRUE,
  seed = 1
)
```

## Arguments

- data:

  Data frame with one row per cell.

- y:

  Character vector of expression column(s) (one per gene).

- perturb:

  Name of the column holding the perturbation (target) of each cell.

- context:

  Name of \*\*one\*\* column giving the cell context that interacts with
  the perturbation: numeric = cell state, factor/character = cell type.
  To study both, fit the model twice. \`NULL\` for no interaction.

- covariates:

  Other covariates (main effects only).

- perturb_group:

  Optional: a column name giving a perturbation group for each
  perturbation, or a named list of perturbation vectors. Each group then
  gets its own \`P\` and \`P x context\` components.

- control:

  Perturbation labels to drop (e.g. non-targeting controls). \`NA\` and
  \`""\` are always dropped.

- min_cells:

  Minimum cells per perturbation (perturbations with fewer are dropped).

- cat_mode:

  For categorical variables in \`context\`: \`"pooled"\` fits one shared
  component (pairs from the same level) and \`"per_level"\` fits one per
  level.

- context_main:

  Include main-effect components for each context variable (default
  \`TRUE\`).

- jackknife:

  Delete-perturbation-block jackknife standard errors.

- n_blocks:

  Number of jackknife blocks over perturbations (default one per
  perturbation).

- scale_x:

  Standardize continuous context and covariates (default \`TRUE\`).
  Expression \`y\` is always standardized to mean 0 and variance 1, so
  variance components are on the scale of standardized expression.

- seed:

  Seed for assigning donors to jackknife blocks.

## Value

An object of class \`"schint_pert"\` with long-format \`coefficients\`
and \`summary\` tables (one block of rows per gene).

## Details

For cells \`a\`, \`b\` the regressors are \`1\`, \`1(pert_a ==
pert_b)\`, \`u_a u_b\` and \`1(pert_a == pert_b) u_a u_b\`; for a
categorical context \`u_a u_b\` is the indicator that both cells are in
the same level (so \`P x celltype\` is "same perturbation and same cell
type"). Several genes can be fitted at once; the design-dependent
\`X'X\` is computed only once.

## Examples

``` r
data(schint_pert_example)
fit <- schint_pert(schint_pert_example, y = c("GENE_A", "GENE_B"),
                   perturb = "perturb", context = "state",
                   control = "NT", covariates = c("pc1", "pc2"))
fit
#> scHINT perturbation model: 3470 cells, 60 perturbations, 2 genes
#> 
#>    gene            term  variance
#>  GENE_A               P  0.031500
#>  GENE_A             PxC  0.051100
#>  GENE_A         P_total  0.082600
#>  GENE_A         context -0.000834
#>  GENE_A       covariate  0.040500
#>  GENE_A total_explained  0.122000
#>  GENE_B               P  0.053900
#>  GENE_B             PxC -0.001590
#>  GENE_B         P_total  0.052300
#>  GENE_B         context  0.000292
#>  GENE_B       covariate  0.056100
#>  GENE_B total_explained  0.109000
```
